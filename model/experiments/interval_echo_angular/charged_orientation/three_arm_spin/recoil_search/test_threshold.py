"""Source-visible limits, opposite current orders and a receiver-local witness."""

import json
import numpy as np

from .threshold import HERE, SHOTS, WITNESS_SHOTS, run


def test_lowest_mobile_band_and_bound_alternative():
    saved = json.loads((HERE / "threshold.json").read_text())
    fresh = run()
    for data in (saved, fresh):
        assert data["changed_kappa4_homogeneous_bottom"] == 8
        energies = []
        for row in data["lengths"]:
            below = row["below"]
            assert [w["classified_as_visible"] for w in below] == [
                False,
                True,
                True,
                False,
            ]
            assert abs(below[1]["center"] - 8.9381807248) < 1e-9
            assert (
                below[1]["mid99"][1] < 2e-9
                and below[1]["counts_mid_far_other_refusal"][1] == 0
            )
            low, high = row["lowest_and_next_harmonic"]
            energies.append(low["energy"])
            for line in [low, high]:
                assert (
                    line["classified_as_visible"]
                    and line["filtered_mode_weight"] > 0.999999
                )
                assert (
                    line["remote_original99"][0] > 1e-11
                    and line["winding"]["witness99"][1] < 0
                )
                assert all(
                    sum(bank) == WITNESS_SHOTS
                    for bank in line["winding"]["counts_same_different_all_refusals"]
                )
                assert sum(line["counts_mid_far_other_refusal"]) == SHOTS
            assert low["mid99"][0] > 0.55
            a, b, c, d, e = row["free_beat"]
            assert (
                max(a["far_half99"] + e["far_half99"]) < 0.14
                and c["far_half99"][0] > 0.85
            )
            assert b["current99"][0] > 0.01 and d["current99"][1] < -0.01
            for control in row["incoherent_control"]:
                assert 0.47 < control["far_half99"][0] < control["far_half99"][1] < 0.51
            for sample in row["free_beat"]:
                assert (
                    sum(sample["current_counts_forward_backward_else_refused"]) == SHOTS
                )
                assert sum(sample["counts_far_near_refused"]) == SHOTS
        assert 16 < energies[2] < energies[1] < energies[0] < 16.014
    for old, new in zip(saved["lengths"], fresh["lengths"], strict=True):
        np.testing.assert_allclose(
            old["scaled_first_excess"], new["scaled_first_excess"], rtol=5e-10
        )
        for a, b in zip(old["free_beat"], new["free_beat"], strict=True):
            np.testing.assert_allclose(
                [a["current"], a["far_half"]], [b["current"], b["far_half"]], atol=2e-6
            )
