"""Replay charged spin transport on two periodic interacting material rings."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from ....replay import compare
from ..charged_endpoints.measure import SOURCE_REVISION
from . import certificate, contact_control, ring, spectroscopy


HERE = Path(__file__).resolve().parent


def clock_plan(length):
    # Map the 32-unit spectral circle to a uniform Fourier clock, with its
    # seam in an empty portion of the spectrum. Use the same native H for
    # the pass test and its coherent uncomputation.
    margin = F(1, 200) if length == 4 else F(1, 2500)
    bins = 2**17 if length == 4 else 2**21
    single_error = 1 / (2 * (bins * margin / 32 - 1))
    repeated_error = 2**19 * single_error**10
    assert repeated_error < F(1, 10**10)
    return {
        "energy_circle": 32,
        "raw_energy_at_phase_zero": -8,
        "readout": "return binary majority and require all Fourier registers to erase to blank",
        "bins_per_coherent_clock": bins,
        "repetitions": 19,
        "spectral_margin_for_each_binary_gate": str(margin),
        "single_clock_tail_upper": float(single_error),
        "accepted_Kraus_error_upper": float(repeated_error),
        "three_gate_joint_probability_absolute_error_upper": float(6 * repeated_error),
    }


def cross_check(record, exact):
    for bank, bounds in zip(record["sign_banks"], exact["banks"], strict=True):
        assert bank["sign"] == bounds["sign"]
        window = bank["two_energy_window"]
        for gate in ("parked", "half_beat"):
            for site, bound in enumerate(bounds["local_read_intervals"][gate]):
                lo, hi = map(F, bound["ground_return_per_start"])
                assert (
                    float(lo)
                    <= window["three_waits"][gate]["joint_spin_return"][site]
                    <= float(hi)
                )
                lo, hi = map(F, bound["spin_retention_given_single"])
                assert (
                    float(lo)
                    <= window["three_waits"][gate]["local_spin_retention_given_single"][
                        site
                    ]
                    <= float(hi)
                )
        assert (
            float(
                F(
                    bounds["local_read_intervals"]["half_beat"][2][
                        "local_Bell_absolute_CHSH"
                    ][0]
                )
            )
            > 2.64
        )


def run():
    ring.source()
    exact = certificate.run()
    rings = [spectroscopy.experiment(length) for length in (4, 8)]
    cross_check(rings[0], exact)
    for record in rings:
        length = record["length"]
        plan = clock_plan(length)
        allowance = float(F(plan["spectral_margin_for_each_binary_gate"]))
        assert record["half_filled_ground_gap"] / 2 > allowance
        record["clock_plan"] = plan
        if length == 4:
            finite = []
            for bounds in exact["banks"]:
                gate = bounds["local_read_intervals"]["half_beat"][2]
                # For each setting bound both the raw signed correlation and
                # its admission probability before dividing by that probability.
                phase_error = F(bounds["half_beat_clock"]["relative_phase_error_upper"])
                eps = F(6, 10**10) + 2 * phase_error
                admission = F(gate["branch_and_single_per_original"][0])
                score = F(gate["local_Bell_absolute_CHSH"][0]) - 8 * eps / (
                    admission - eps
                )
                assert score > F(2642, 1000)
                finite.append(
                    {
                        "sign": bounds["sign"],
                        "raw_event_or_signed_expectation_error_upper": float(eps),
                        "rounded_clock_local_CHSH_certified_lower": 2.642,
                    }
                )
            plan["rounded_half_beat_at_four_sites"] = finite
        for bank in record["sign_banks"]:
            window = bank["two_energy_window"]
            assert window["nearest_excluded_above_band"] / 2 > allowance
            parked = window["three_waits"]["parked"]
            moved = window["three_waits"]["half_beat"]
            for key in parked:
                np.testing.assert_allclose(
                    np.roll(parked[key], length // 2), moved[key], atol=3e-10
                )
            # Recoil on both sides of the +pi/2 minimum, at the same
            # |change of covariant wave number| with no mass or speed fit.
            near = length // 4
            lines = bank["lines"]
            p = [
                lines[(near + d) % length]["reservoir_corrected_transfer"]
                for d in (-1, 0, 1)
            ]
            bank["held_equal_offset_recoils"] = [p[0] - p[1], p[2] - p[1]]
    return {
        "schema": "oph.interacting_contact_spin_branch.v1",
        "source_revision": SOURCE_REVISION,
        "archived_source_sha256": ring.ANCESTRY,
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "exact_four_site_enclosures": exact,
        "four_site_fixed_instrument_contact_control": contact_control.run(rings[0]),
        "ring_measurements": rings,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    result = run()
    if options.output:
        with options.output.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare(
            result,
            json.loads((HERE / "results.json").read_bytes()),
            "contact spin branch",
        )
        print(
            "Replayed both charged spin branches, exact four-site translations and local clock/receiver counts."
        )


if __name__ == "__main__":
    main()
