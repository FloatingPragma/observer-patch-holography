"""Prepare a contact reference using ordinary mismatch and odd-parity pilots."""

import argparse
from fractions import Fraction as F
import json
from math import log

import numpy as np

from .. import measure as flights
from ...energy_source import instrument
from .measure import HERE, check
from .modular import EQUAL, EVEN, ODD, H, calibrate
from .reservoir import counts


def pilots():
    # Extract parity from the ordinary matter receiver's independent interferometer.
    read = instrument.instruments(1)[1][0][:2, :2]
    np.testing.assert_allclose(read, ODD[1:3, 1:3], atol=3e-14)
    spectral = []
    for sign in (1, -1):
        psi = flights.receiver.receiver.standing_wave(1, 8, sign)
        state = np.outer(psi, psi.conj())
        odd = float(np.trace(read @ state).real)
        mismatch = float(state.trace().real)
        np.testing.assert_allclose([odd, mismatch], [round(odd), 1], atol=3e-14)
        spectral.append((round(odd), 1))
    assert sorted(spectral) == [(0, 1), (1, 1)]
    return [(0, 0), *sorted(spectral), (0, 0)]  # odd result, mismatch result


def run():
    enabled = pilots()
    rng = np.random.default_rng(27111942)
    proposals = 2**21
    starts = rng.multinomial(proposals, [0.25] * 4)
    survive = F(7, 8)
    source = json.loads(HERE.parent.joinpath("formed_heads.json").read_text())["source"]
    output = []
    for name, ncontrast, nodd in (
        ("ordinary energy", 8, 1),
        ("one mismatch", 1, 1),
        ("parity disabled", 8, 0),
    ):
        survivors = []
        stopped = []
        likelihoods = []
        for count, (odd, mismatch) in zip(starts, enabled, strict=True):
            remaining = int(count)
            rejected = []
            switches = [mismatch] * ncontrast + [odd] * nodd
            for active in switches:
                lost = int(rng.binomial(remaining, 1 / 8)) if active else 0
                rejected.append(lost)
                remaining -= lost
            likelihoods.append(survive ** sum(switches))
            stopped.append(rejected)
            survivors.append(remaining)
            assert sum(rejected) + remaining == count
        Z = sum(likelihoods)
        weights = [v / Z for v in likelihoods]
        ref = (
            2 * float(weights[0]) * EQUAL
            + float(weights[1]) * EVEN
            + float(weights[2]) * ODD
        )
        beta = -log(float(survive))
        vals, vec = np.linalg.eigh(ref)
        modular = -(vec * np.log(vals)) @ vec.T
        err = float(
            np.linalg.norm(modular - log(float(Z)) * np.eye(4) - beta * H, ord=2)
        )
        if name == "ordinary energy":
            assert err < 2e-14
        read = calibrate(ref, EVEN - (8 * ODD + EQUAL) / 9)
        held = []
        for row in source:
            if row["orientation"] == "matching":
                mismatch_change = row["independent_unlike"] - row["formed_unlike"]
                signal = calibrate(ref, mismatch_change * ((EVEN + ODD) / 2 - EQUAL))
                held.append(
                    dict(
                        source_sign=row["source_sign"],
                        modular=signal["held_modular"],
                        ordinary_measured_energy=signal["held_energy"],
                        modular_over_ordinary_pilot=signal["held_modular"]
                        / read["beta_equal_to_even"],
                    )
                )
        output.append(
            dict(
                design=name,
                ordinary_inputs=list(map(int, starts)),
                stops_by_gate=stopped,
                survivors=survivors,
                weights=[str(t) for t in weights],
                probability_of_acceptance=float(Z / 4),
                full_generator_residual=err,
                matched_banks_modular=read["held_modular"],
                source=held,
            )
        )
    # Conditioning independent accepted contacts on both total observables
    # gives the same weight to each allowed eigenword. Proposals for local
    # contacts can be refilled before preparing a batch for the total reading.
    Z = 2 + survive**8 + survive**9
    batches = [
        dict(
            contacts=m,
            probability_totals_given_accepted_contacts=float(
                counts(m, 2 * m, m) * survive ** (2 * m) / Z**m
            ),
        )
        for m in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
    ]
    return dict(
        seed=27111942,
        proposals=proposals,
        panels=output,
        microcanonical_selection=batches,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "herald.json"
    result = run()
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
