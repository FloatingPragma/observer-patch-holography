"""Acquire independent source apertures, spectra and downstream material receiver records."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..common_receiver import measure as guide
from . import instrument

HERE = Path(__file__).parent
SHOTS = 2**27
SOURCE_WIDTH = float(np.sqrt(np.log(2 * 64 / 0.005) / (2 * SHOTS))) + 1e-12


def mark_counts(rng, psi, projectors, incomplete):
    draws = rng.multinomial(SHOTS, [0.25] * 4)
    probs = np.einsum("i,kij,j->k", psi.conj(), projectors, psi).real
    if incomplete:
        probs[0] = 0
    counts = rng.binomial(draws, np.rint(np.clip(probs, 0, 1) * 2**40) / 2**40)
    return {"edge_trials": draws.tolist(), "edge_clicks": counts.tolist()}


def run():
    rng = np.random.default_rng(261126)
    spectra, records = [], []
    for load in (1, 2):
        _, projectors, sorter = instrument.instruments(load)
        for flips in (1, 2):
            p = instrument.beat_probability(load, flips, 0.6)
            even = int(rng.binomial(SHOTS, round(p * 2**40) / 2**40))

            def gap(prob):
                return float(np.arccos(np.clip(2 * prob - 1, -1, 1)) / 0.6)

            rate = even / SHOTS
            spectra.append(
                dict(
                    load=load,
                    flips=flips,
                    even=even,
                    gap=gap(rate),
                    interval=[gap(rate + SOURCE_WIDTH), gap(rate - SOURCE_WIDTH)],
                )
            )
        for name, psi in instrument.preparations(load).items():
            weights = np.rint(abs(sorter @ psi) ** 2 * 2**40) / 2**40
            if weights.sum() != 1:
                raise ArithmeticError(
                    "prepared four-bin spectrum missed the dyadic grid"
                )
            spectral = rng.multinomial(SHOTS, weights)
            row = dict(
                load=load,
                stock=name,
                energy_ports=spectral.tolist(),
                excitation=float(spectral @ instrument.BITS / SHOTS),
                full=mark_counts(rng, psi, projectors, False),
            )
            if name in ("first dimer", "opposite dimer"):
                row["missing_first_edge"] = mark_counts(rng, psi, projectors, True)
            records.append(row)
    cal = next(r for r in records if r["load"] == 1 and r["stock"] == "one")
    gain = sum(cal["full"]["edge_clicks"]) / SHOTS / cal["excitation"]
    results = []
    u = guide.visits_before_drain()
    for row in records:
        if (row["load"], row["stock"]) not in ((1, "one"), (2, "top"), (2, "ground")):
            continue
        releases = sum(row["full"]["edge_clicks"])
        visits, ticks = (
            guide.parent.visits(rng, releases)
            if releases
            else (np.zeros_like(u, dtype=int), 0)
        )
        for probe in (1, 2):
            setting = guide.Setting("matched", (6, 5, 5, 0), probe, (2, 2, 2))
            arms = guide.count_panel(rng, visits, setting)
            reward = guide.signed_readout(setting)
            per_release = float(np.sum(u * reward))
            # Conditional on observed clicks, propagate the receiver bound; also
            # cover the rate of those clicks as a source of unconditional uncertainty.
            width = abs(per_release) * SOURCE_WIDTH
            if releases:
                width += releases / SHOTS * guide.bound(reward, u, 6, releases, 0.005)
            shift = (sum(sum(a["right"]) for a in arms) - int(visits.sum()) / 2) / SHOTS
            results.append(
                dict(
                    source_load=row["load"],
                    stock=row["stock"],
                    probe=probe,
                    releases=releases,
                    last_drain=ticks,
                    arms=arms,
                    per_start=shift,
                    predicted=row["excitation"] * per_release / 4,
                    width=width,
                )
            )
    return dict(
        seed=261126,
        starts=SHOTS,
        source_probability_width=SOURCE_WIDTH,
        calibrated_clicks_per_excitation=gain,
        spectra=spectra,
        sources=records,
        downstream=results,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        if any(
            a["full"] != b["full"]
            or a["energy_ports"] != b["energy_ports"]
            or a.get("missing_first_edge") != b.get("missing_first_edge")
            for a, b in zip(output["sources"], saved["sources"], strict=True)
        ):
            raise ValueError("source banks differ")
        if [s["even"] for s in output["spectra"]] != [
            s["even"] for s in saved["spectra"]
        ] or [s["arms"] for s in output["downstream"]] != [
            s["arms"] for s in saved["downstream"]
        ]:
            raise ValueError("beat or downstream banks differ")
