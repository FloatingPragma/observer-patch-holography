import itertools
import json

import numpy as np

from ..receipt_replay import check_binomial, pinned_receipt
from . import measure, passive_work


def test_native_energy_screen_and_law_control():
    data = measure.run()
    stored = json.loads((measure.HERE / "results.json").read_text())
    for charge in ("-1", "1"):
        record = data["exposures"][charge]
        stored_record = stored["exposures"][charge]
        assert (
            0
            < record["exchange_only_reference_coefficient"]
            < record["active_reference_coefficient"] / 100
        )
        for panel, saved in zip(record["finite"], stored_record["finite"]):
            assert panel["counts"] == saved["counts"]
            for name in ("ordinary", "pairing"):
                observed = panel["count_slopes"][name]
                population = panel["matched_per_time_four"][name]
                assert abs(observed["estimate"] - population) < observed["half_width"]
                assert observed["estimate"] + observed["half_width"] < 0
    assert data["coefficient_slopes"] == stored["coefficient_slopes"]


def test_two_copy_sort_and_collective_limit():
    data = passive_work.run()
    snapshot = json.loads((measure.HERE / "passive_work.json").read_text())
    assert data["two_copy_one_swap_work_lower"] > 0.012
    for n, row in data["finite"].items():
        for key, number in row.items():
            np.testing.assert_allclose(
                number, snapshot["finite"][n][key], rtol=1e-9, atol=1e-12
            )
    assert abs(data["finite"]["1"]["pairing_derivative"]) < 1e-12
    assert 0.98 < data["finite"]["12"]["ordinary_calibrated_pairing_fraction"] < 0.99
    for name in ("ordinary_derivative", "pairing_derivative"):
        np.testing.assert_allclose(
            data["infinite"][name], snapshot["infinite"][name], rtol=1e-10
        )


def test_type_sort_against_explicit_energy_levels():
    state = measure.preparations()["equilibrium"]
    p = np.linalg.eigvalsh(state).ravel()
    perturb = np.stack([np.roll(p, 1) - p, p[::-1] - p]).T
    n = 3
    words = np.array(list(itertools.product(range(8), repeat=n)))
    shell = np.sort((words >= 4).sum(axis=1) / n)

    def direct(populations):
        return shell @ np.sort(np.prod(populations[words], axis=1))[::-1]

    estimate = passive_work.sorted_energy(n, p, perturb)
    assert abs(estimate[0] - direct(p)) < 1e-13
    for offset, direction in enumerate(perturb.T, 1):
        probabilities = np.prod(p[words], axis=1)
        descending = np.argsort(-probabilities)
        derivative = probabilities * (direction[words] / p[words]).sum(axis=1)
        assert abs(shell @ derivative[descending] - estimate[offset]) < 1e-13


def test_same_exposure_ordinary_calibration_on_beam_ensemble():
    from . import four_ticks

    report = four_ticks.run()
    snapshot = pinned_receipt(
        measure.HERE / "four_ticks.json",
        "522a24bea02e07aff83a9b8f8638ae9e83f87ada8b32aa258a3e88133ca20963",
    )
    assert {k: v for k, v in report.items() if k != "charges"} == {
        k: v for k, v in snapshot.items() if k != "charges"
    }
    assert report["time"] == 4 and report["launches"] == [0, 5, 9, 8, 7]
    shots = report["shots_per_setting"]
    log_level = np.log(2 * 2 * 5 * 5 / 0.01)
    both_records_log_level = np.log(2 * 2 * 2 * 5 * 5 / 0.01)
    angles = {
        ("ordinary" if arm["source"] == "ordinary" else "pairing"): arm[
            "angle_per_unit"
        ]
        for arm in json.loads(
            (measure.HERE.parent / "signed_transfers.json").read_text()
        )["tests"]
    }
    pixels = json.loads((four_ticks.beam.HERE / "counts.json").read_text())["species"]
    for charge in ("-1", "1"):
        live = report["charges"][charge]
        frozen = snapshot["charges"][charge]
        for launch in map(str, report["launches"]):
            law = live["settings"][launch]["zero_energy_probabilities"]
            saved_law = frozen["settings"][launch]["zero_energy_probabilities"]
            for label, prob in law.items():
                np.testing.assert_allclose(
                    prob, saved_law[label], atol=1e-14, rtol=2e-12
                )
                for data in (frozen, live):
                    check_binomial(
                        data["settings"][launch]["zero_energy_counts"][label],
                        shots,
                        prob,
                        both_records_log_level,
                    )
        for key, launches in (
            ("port_0", [0]),
            ("equal_beam_launches", report["launches"]),
        ):
            reference = frozen["surveys"][key]
            for data in (frozen, live):
                panel = data["surveys"][key]
                assert panel["energy"] == reference["energy"]
                counts = [
                    data["settings"][str(i)]["zero_energy_counts"] for i in launches
                ]
                slopes = four_ticks.observed_energy_slopes(
                    counts, shots, log_level, angles
                )
                assert slopes == panel["synthetic_energy_slope"]
                if key == "equal_beam_launches":
                    for name in (
                        "energy_ordinary_over_pairing",
                        "image_observed_over_energy_prediction",
                    ):
                        assert panel[name] == reference[name]
                    lower, upper = json.loads(
                        panel["image_observed_over_energy_prediction"]
                    )
                    assert (
                        0.85 < lower < upper < 0.94
                        if charge == "-1"
                        else -0.164 < lower < upper < -0.14
                    )
                    interval = four_ticks.optical_count_ratio(
                        pixels[charge]["signed_slopes_99pct"], slopes
                    )
                    assert str(interval) == panel["nominal_count_ratio_at_least_98pct"]
                    assert interval.b < 1
                    for name, item in slopes.items():
                        lo, hi = json.loads(panel["energy"][name])
                        assert lo < item["energy_slope"] < hi
