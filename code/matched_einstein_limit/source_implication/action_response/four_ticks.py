"""Compare actual native seam energy and directional image response at the same exposure."""

import argparse
import json
from decimal import Decimal

import numpy as np
import sympy as sp
from mpmath import iv

from ..beam_response import measure as beam
from . import measure


def work_effects(board, route, launches, charge, coeff):
    """Exact rational approximants of the compressed zero-energy effect at t=4."""
    n = len(launches)
    nodes = (*launches, *route[2:])
    assert len(set(nodes)) == len(nodes)
    p = measure.front.PRIOR.sym_matrix(measure.receiving_effect(charge))
    arms = [[] for _ in launches]
    for odd in (False, True):
        h = measure.native.generator(board, route, charge, odd)
        raw = measure.native.kernel_panels(h, beam.indices(nodes), coeff)[:, 7, n:, :n]
        for i, row in enumerate(arms):
            wave = (
                sp.Matrix(raw[0, :, i].reshape(8, 4))
                + sp.I * sp.Matrix(raw[1, :, i].reshape(8, 4))
            ) / measure.native.prior.FRAC
            effect = wave.H * p * wave
            row.append(effect.applyfunc(sp.expand))
    return arms


def parse_interval(s):
    return iv.mpf([str(v) for v in json.loads(s, parse_float=Decimal)])


def raw_effect(e, density, trace_error, kernel_error):
    full = measure.front.full
    center = sum(full.inner(m, x) for m, x in zip(e, density))
    norm = max(measure.front.PRIOR.one_norm(m) for m in e)
    return full.ri_band(center, norm * trace_error + kernel_error)


def volume_ratio(source):
    # Both fixed-volume directions use the same unitary rotation. Its volume
    # derivative drops out of the ratio of the two adjustment angles.
    mid, targets, *_ = source
    nums, scale, error = measure.native.exact_effects()
    full = measure.front.full
    effects = full.rational_effects(nums, scale)
    bounds = [
        raw_effect(e, mid, sp.Rational("6e-7"), sp.Rational(error)) for e in effects
    ]
    _, grad = full.enc.geometry_jet(bounds)
    d = {}
    for name, radius in (("ordinary", "12e-7"), ("unpair", "18e-7")):
        slopes = []
        for e in effects:
            # Traceless source tangents remove the identity part exactly.
            center = sum(x.trace() for x in e) / 8
            centered = [x - center * sp.eye(4) for x in e]
            slopes.append(
                raw_effect(
                    centered, targets[name], sp.Rational(radius), 2 * sp.Rational(error)
                )
            )
        d[name] = sum(a * b for a, b in zip(grad.g, slopes))
    return d["ordinary"] / d["unpair"]


def observed_energy_slopes(counts, shots, log_level, angles):
    probs = {k: np.mean([row[k] / shots for row in counts]) for k in counts[0]}
    errors = {
        k: np.mean(
            [
                np.sqrt(2 * row[k] / shots * log_level / shots) + 2 * log_level / shots
                for row in counts
            ]
        )
        for k in counts[0]
    }
    measured = {}
    for k, angle in angles.items():
        point = -8 * (
            probs[k]
            - probs["equilibrium"]
            + angle * (probs["plus"] - probs["minus"]) / 2
        )
        width = 8 * (
            errors[k]
            + errors["equilibrium"]
            + abs(angle) * (errors["plus"] + errors["minus"]) / 2
        )
        measured[k] = {"energy_slope": point, "sampling_halfwidth": width}
    return measured


def optical_count_ratio(pixel, measured):
    spike = parse_interval(pixel["pairing"]) / parse_interval(pixel["ordinary"])

    def count_interval(name):
        item = measured[name]
        return iv.mpf(item["energy_slope"]) + iv.mpf([-1, 1]) * iv.mpf(
            item["sampling_halfwidth"]
        )

    return spike * count_interval("ordinary") / count_interval("pairing")


def run():
    iv.dps = 30
    source = measure.front.inputs()
    rho = volume_ratio(source)
    states = measure.preparations()
    arms = json.loads((measure.HERE.parent / "signed_transfers.json").read_text())[
        "tests"
    ]
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in arms
    }
    board = measure.native.iso.pairing.formation.federation(1)
    route = measure.native.iso.pairing.formation.region(board)
    launches = beam.screens(board, route)[0]
    assert launches == [0, 5, 9, 8, 7]
    coeff = measure.native.prior.weights()
    optical = json.loads((beam.HERE / "certificate.json").read_text())
    pixel_counts = json.loads((beam.HERE / "counts.json").read_text())["species"]
    panels = {}
    rng = np.random.default_rng(402627)
    shots = 2**58
    # Both signs, five launches and five preparations.
    log_level = np.log(2 * 2 * len(launches) * len(states) / 0.01)
    for q in (-1, 1):
        entries = work_effects(board, route, launches, q, coeff)
        read = {}
        mean_p = []
        for launch, blocks in zip(launches, entries):
            floats = np.array([np.array(b, dtype=complex) for b in blocks])
            # Independent unmonitored floating propagation evaluates the detector at the endpoints.
            for odd in (False, True):
                h = measure.native.generator(board, route, q, odd)
                wave = measure.wave_at(h, launch, 4)[list(route[2:])].reshape(8, 4)
                alternative = wave.conj().T @ measure.receiving_effect(q) @ wave
                np.testing.assert_allclose(
                    alternative, floats[int(odd)], rtol=2e-10, atol=1e-14
                )
            probabilities = {
                k: float(np.einsum("bij,bji->", floats, x).real)
                for k, x in states.items()
            }
            counts = {k: int(rng.binomial(shots, p)) for k, p in probabilities.items()}
            read[str(launch)] = {
                "zero_energy_counts": counts,
                "zero_energy_probabilities": probabilities,
            }
            mean_p.append(counts)
        group = {}
        for key, ids in (
            ("port_0", [0]),
            ("equal_beam_launches", list(range(len(launches)))),
        ):
            effect = [
                sum((entries[i][b] for i in ids), sp.zeros(4)) / len(ids)
                for b in range(2)
            ]
            jets = measure.source_jets(effect, source)
            # Tail for exact propagation and dyadic quantization. A four-column
            # submatrix of U has norm error below 2e-19; ||F - F_tilde||<1e-18.
            # Both matched source tangents have trace norm at most 2+997.
            allowance = iv.mpf(["-1e-15", "1e-15"])
            energy = {k: -8 * (j + allowance) for k, j in jets.items()}
            measured = observed_energy_slopes(
                [mean_p[i] for i in ids], shots, log_level, angles
            )
            record = dict(
                energy={k: str(v) for k, v in energy.items()},
                synthetic_energy_slope=measured,
            )
            if key == "equal_beam_launches":
                sky = optical[str(q)]
                ao, ap = map(parse_interval, (sky["ordinary"], sky["unpairing"]))
                assert (
                    0 not in energy["ordinary"]
                    and 0 not in energy["pairing"]
                    and 0 not in ao
                )
                # E_o/E_p = rho+(e_o-rho*e_p)/E_p, where e's are raw
                # energy tangents and rho=V_o/V_p. The rotation cancels exactly.
                raw_o = -8 * raw_effect(
                    effect,
                    source[1]["ordinary"],
                    sp.Rational("12e-7"),
                    sp.Rational("2e-18"),
                )
                raw_p = -8 * raw_effect(
                    effect,
                    source[1]["unpair"],
                    sp.Rational("18e-7"),
                    sp.Rational("2e-18"),
                )
                direct = energy["ordinary"] / energy["pairing"]
                joint = rho + (raw_o - rho * raw_p) / energy["pairing"]
                shared = iv.mpf([max(direct.a, joint.a), min(direct.b, joint.b)])
                ratio = (ap / ao) * shared
                record["energy_ordinary_over_pairing"] = str(shared)
                record["image_observed_over_energy_prediction"] = str(ratio)
                assert 1 not in ratio
                pixel = pixel_counts[str(q)]["signed_slopes_99pct"]
                count_ratio = optical_count_ratio(pixel, measured)
                record["nominal_count_ratio_at_least_98pct"] = str(count_ratio)
            group[key] = record
        panels[str(q)] = {"settings": read, "surveys": group}
    return {
        "time": 4,
        "ratio_of_volume_adjustment_angles": str(rho),
        "launches": launches,
        "energy_seam": list(route[2:]),
        "seed": 402627,
        "shots_per_setting": shots,
        "charges": panels,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (measure.HERE / "four_ticks.json").write_text(out)
    else:
        print(out, end="")
