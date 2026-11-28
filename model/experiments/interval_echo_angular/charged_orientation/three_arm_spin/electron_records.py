"""Verify the shared source and independent readouts of the electron studies."""

from fractions import Fraction
from hashlib import sha256
import json
from math import isclose
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def run():
    threshold = read("recoil_search/threshold.json")
    charge = read("charge_probe/threshold_records/results.json")
    recoil = read("recoil_search/electron_identity/records.json")
    magnetic = read("magnetic_response/results.json")
    exchange = read("echo_turn/intact_exchange/results.json")
    authority = sha256((HERE / "recoil_search/threshold.json").read_bytes()).hexdigest()
    assert charge["source_threshold_sha256"] == authority
    assert recoil["threshold_receipt_sha256"] == authority
    assert magnetic["threshold_receipt_sha256"] == authority
    assert charge["receiving_owner"] == magnetic["receiving_owner"] == 13
    assert exchange["intact_turns"]["receiver"] == charge["receiving_tip"] == 160
    checked = []
    for original, c, r, m, x in zip(
        threshold["lengths"],
        charge["lengths"],
        recoil["native"],
        magnetic["cases"],
        exchange["charts"],
        strict=True,
    ):
        length, center = (
            original["length"],
            original["lowest_and_next_harmonic"][0]["energy"],
        )
        assert length == c["length"] == r["length"] == m["arm_links"] == x["length"]
        assert r["source_and_receiving_hubs"] == x["vacancy_hubs"] == [4, 160]
        front = r["readings"][0]
        assert abs(front["prepared_center"] - center) < 1e-11
        lower, upper = m["line_enclosure"]["interval"]
        assert lower <= center <= upper
        source_record = c["middle"][0]
        assert (
            source_record["preparation"] == "lowest line" and source_record["time"] == 0
        )
        accepted = source_record["admission"]
        for check in (
            front["source_admission_before_detector"],
            m["conditional_clock_acceptance"] / 3,
        ):
            assert isclose(accepted, check, rel_tol=1e-9, abs_tol=0)
        assert isclose(
            front["finite_clock_main_line_fraction"],
            m["conditional_filter_weight"],
            rel_tol=0,
            abs_tol=1e-11,
        )
        assert c["cuts"] == original["physical_cuts"]
        for a, b in c["cuts"]:
            assert any(
                [a, b] in (path[j : j + 2] for j in range(length))
                for path in x["arm_paths"]
            )
        negative = c["boundary"]["-1"]
        assert negative["admitted_tip_reading"] == [-1, 1]
        assert negative["boundary_idle_given_signed_tip"] == 1
        assert negative["signed_tip_arrivals"] > 0
        assert source_record["native_multiple_deposits_given_clock"] == 0
        assert r["three_upper_harmonic_secants"]["positive_massive_convexity_rejected"]
        assert all(
            lo <= 0 <= hi for lo, hi in r["held_energy_minus_native_cosine_forecast"]
        )
        assert m["physical_edges"] == 3 * length + 6
        assert (
            max(row["maximum_conditional_score"] for row in m["zero_intrinsic_current"])
            < 1e-10
        )
        assert m["maximum_changed_law_difference"] > 1e-8
        assert len(x["matched_pairs"]) == 6
        assert all(
            pair["head_operator"] == "identity_on_64_words"
            and pair["native_amplitudes"] == [1, 1]
            and pair["antisymmetric_amplitudes"] == [-1, 1]
            for pair in x["matched_pairs"]
        )
        checked.append(
            {
                "physical_arm_links": length,
                "clock_center": center,
                "total_clock_acceptance_per_original_start": accepted,
                "signed_arrivals": negative["signed_tip_arrivals"],
                "observed_boundary_idle_given_arrival": negative[
                    "boundary_idle_given_signed_tip"
                ],
                "same_guide_held_dispersion_tests": len(
                    r["held_energy_minus_native_cosine_forecast"]
                ),
                "stationary_current_contacts": m["physical_edges"],
                "two_vacancy_route_pairs": len(x["matched_pairs"]),
            }
        )
    banks = exchange["banks"]
    assert [Fraction(x) for x in banks["native_phase"][:2]] == [Fraction(1, 81), 0]
    assert [Fraction(x) for x in banks["antisymmetric_changed_rule"][:2]] == [
        0,
        Fraction(1, 81),
    ]
    assert len(exchange["intact_turns"]["forbidden_rotation_indices"]) == 59
    return checked


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
