"""Independent checks of the phase meter, terminal flights and tidal inference."""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

pytest.importorskip("strict_history", reason="needs the oph-model-sim checkout on PYTHONPATH")

from level202_gravimeter import WideReplay
from matched_einstein_limit import observer_field


HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def run():
    return json.loads((HERE / "level202_gravimeter.json").read_bytes())


def test_home_source_clock_and_strict_releases(run):
    assert (
        run["protocol_sha256"]
        == sha256((HERE / "LEVEL202_GRAVIMETER.md").read_bytes()).hexdigest()
    )
    assert (
        run["clock_archive_sha256"]
        == sha256((HERE / "level202_echo_reach.json").read_bytes()).hexdigest()
    )
    assert abs(run["inside_minus_outside_at_write"]) == 2
    assert 0 < run["archived_write_tick"] < 96
    assert run["strict_one_unit_releases"]["launches"] == {
        "add_one": 25,
        "remove_one": 28,
    }
    ranges = np.array(run["roundtrip_ticks"]) / 2
    for trial in run["strict_one_unit_releases"]["terminal_reads"]:
        assert abs(trial["released_load"] - trial["sham_load"]) == 1
        for sent, received in trial["query_and_return_ticks"]:
            assert (received - sent) / 2 == ranges[0, trial["terminal"]]
    meter = run["receiver_side_energy_meter"]
    for mode, energy in (("even", 32), ("odd", 33)):
        data = meter["readings"][mode]
        expected = np.cos(np.array(meter["holds"]) * energy / 2) ** 2
        np.testing.assert_allclose(
            data["exact_reference_returns"], expected, atol=2e-12
        )
        p = data["short_hold_return_count"] / meter["short_hold_shots"]
        assert np.arccos(2 * p - 1) / meter["holds"][1] == pytest.approx(
            data["estimated_energy"]
        )
        assert abs(data["estimated_energy"] - energy) < 0.05
    for sign in meter["enrollment"]:
        np.testing.assert_allclose(sign["post_read_fidelities"], 1)
        np.testing.assert_allclose(sign["probabilities"], 0.5)
    np.testing.assert_allclose(
        run["measured_exchange_clock"], run["exposure_times"], atol=1e-12
    )


def test_independent_native_flight_on_original_cut(run):
    model = WideReplay(202, 110, run["seed"])
    ports = [12 * int(u) + p for u, p in run["terminal_labels"]]
    index = {p: i for i, p in enumerate(ports)}
    distances = np.array(run["roundtrip_ticks"]) // 2
    source = run["launch_indices"][0]
    starts = run["launch_indices"]
    for bank in run["banks"]:
        if bank["source_mount"] != "settled_exterior":
            continue
        q = bank["tag"]
        for mode, record in bank["source_readings"].items():
            parity = dict(run["exterior_spectator_readouts"])
            parity[source] = mode
            h = np.zeros(distances.shape)
            for i, p in enumerate(ports):
                links = model.links(p)
                h[i, i] = 8 * q * q * (6 - len(links))
                for neighbour, write in links:
                    if not write:
                        h[i, index[neighbour]] = -0.5
                        h[i, i] += 0.5 + 8 * q * q
                    elif neighbour not in index and q == 1:
                        h[i, i] -= 32 + (parity[i] == "odd")
            unit = np.eye(len(ports))[:, starts]
            generator = -1j * sparse.csr_matrix(h - np.eye(len(h)) * h[0, 0])
            mean = []
            for t in (1 - 1 / 16, 1.0, 1 + 1 / 16):
                amplitudes = expm_multiply(t * generator, unit)
                probs = abs(amplitudes) ** 2
                observed = probs[:, 0] @ distances @ probs[:, 1:].mean(axis=1)
                assert observed == pytest.approx(
                    record["mean_ranges"][run["exposure_times"].index(t)], abs=2e-11
                )
                mean.append(observed)
            assert (mean[2] - 2 * mean[1] + mean[0]) * 256 == pytest.approx(
                record["centered_second_differences"][0], abs=2e-10
            )


def test_count_only_acceleration_and_wall_readouts(run):
    for bank in run["banks"]:
        if "measured_acceleration_difference" not in bank:
            assert bank["acceleration_difference"] == 0
            continue
        measurements = []
        for mode in ("even", "odd"):
            repeat = bank["source_readings"][mode]["precision_repeat"]
            means, variances = [], []
            for histogram in repeat["histograms"]:
                counts = np.array(histogram)
                values = np.arange(len(counts))
                n = float(sum(histogram))
                assert n == repeat["shots_per_time"]
                mean = counts @ values / n
                means.append(mean)
                variances.append((counts @ (values - mean) ** 2) / (n * (n - 1)))
            weights = np.array([256, -512, 256])
            measurements.append((weights @ means, (weights**2) @ variances))
        difference = measurements[1][0] - measurements[0][0]
        error = np.sqrt(measurements[1][1] + measurements[0][1])
        assert difference == pytest.approx(bank["measured_acceleration_difference"])
        assert error == pytest.approx(bank["acceleration_difference_se"])
        if bank["tag"] == 1:
            assert difference > 6 * error
        else:
            assert abs(difference) < 4 * error
    for chamber in run["wall_and_volume"]["chambers"]:
        assert chamber["switched_blank_contacts"] > 0
        assert chamber["switch_off_work_each_source"] == 0
        assert chamber["single_loaded_contact_control_work"] == -8.5
    assert len(run["wall_and_volume"]["volume_checks"]) == 10
    for triple in run["wall_and_volume"]["volume_checks"]:
        assert np.linalg.det(triple["twice_gram"]) < 0
    assert run["observer_vacuum_energy_density_J_m3"] is None
    distances = np.array(run["roundtrip_ticks"]) / 2
    for terminals in run["transverse_echoes"]["home_unit_radius_pairwise_diameter_two"]:
        np.testing.assert_array_equal(distances[0, terminals], np.ones(3))
        np.testing.assert_array_equal(
            distances[np.ix_(terminals, terminals)], 2 * (1 - np.eye(3))
        )
    assert run["transverse_echoes"]["home_unit_radius_pairwise_diameter_two"]
    for region in run["transverse_echoes"]["regions"]:
        assert region["flat_min_gram_eigenvalue"] < -1e-5
        assert region["sphere_100_min_gram_eigenvalue"] < -1e-5


def test_vacuum_control_from_held_clock_and_matter_tapes(run):
    for branch in run["controlled_einstein_scalar"]:
        refinements = branch["refinements"]
        errors = [r["max_held_error"] for r in refinements]
        assert errors[0] > errors[1] > errors[2]
        assert errors[2] < 0.0003
        calibration = observer_field.calibrate(branch["fine_calibration_records"])
        for result in refinements[-1]["held"]:
            tape = result["clock_and_scalar_records"]
            dt = tape["tick_unit"]
            field = np.array([r["scalar_meter_ticks"] for r in tape["readings"]]) * dt
            times = np.array([r["proper_clock_ticks"] for r in tape["readings"]]) * dt
            trips = np.array([r["radar_send_return_ticks"] for r in tape["readings"]])
            arms = trips[:, :, 1] - trips[:, :, 0]
            h = (times[2] - times[0]) / 2
            tidal = np.sum((arms[2] + arms[0] - 2 * arms[1]) / arms[1]) / h**2
            v = (field[2] - field[0]) / (2 * h)
            matter = (
                calibration["coupling"]
                * (2 * v * v - calibration["mass_squared"] * field[1] ** 2)
                / 2
            )
            assert tidal + matter == pytest.approx(result["inferred_lambda"], abs=1e-10)
