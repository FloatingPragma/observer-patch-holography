import json
from itertools import product

import numpy as np

from model.experiments.interval_echo_quantum import experiment as ex


def test_interference_archive_reacquired_from_fresh_native_port_evolution():
    archive = json.loads((ex.HERE / 'results.json').read_text())
    measured = ex.interference(ex.federation(0), np.random.default_rng(archive['seed']))
    assert measured['screen_counts_then_absorbed'] == archive['interference']['screen_counts_then_absorbed']
    for counts in measured['screen_counts_then_absorbed'].values():
        assert sum(counts) == measured['shots_per_setting']


def test_native_echo_and_routed_prediction():
    ex.verify_read()
    rho, route = ex.routed_source(ex.federation(0))
    assert route["cell_distance"] == 5
    local = [ex.local_gate(t) for t in (0, np.pi / 2)]
    remote = [ex.local_gate(t) for t in (3 * np.pi / 4, 5 * np.pi / 4)]
    e = [[ex.SCORES @ ex.pointer(rho, a, b) @ ex.SCORES for b in remote] for a in local]
    np.testing.assert_allclose(np.sum(ex.CHSH * e), 1 + np.sqrt(2), atol=1e-12)
    # Any preassigned pair of local results at each independently set device.
    assert {a * (b + d) + c * (b - d) for a, b, c, d in product((-1, 1), repeat=4)} == {
        -2,
        2,
    }


def test_scores_use_both_records_and_all_enrollments():
    record = json.loads((ex.HERE / "results.json").read_text())
    studies = [
        record["free_pair"][k] for k in ("held", "unread_rails", "exchange_only")
    ]
    studies.extend(
        r["measurement"] for r in record["free_pair"]["fixed_axes_nearby_times"]
    )
    studies.extend(record["routed"][k] for k in ("held", "unread_rails"))
    for study in studies:
        merged = np.array(study["counts_xy_positive_negative_quiet"])
        flags = study["capture_then_settings_then_local_records"]
        np.testing.assert_array_equal(merged, sum(np.array(f["counts"]) for f in flags))
        assert merged.sum() == ex.N
        for f in flags:
            c = np.array(f["counts"])
            assert c.sum() == f["enrolled_before_settings"]
            if f["alice_capture"]:
                assert c[:, :, 2, :].sum() == 0
            else:
                assert c[:, :, :2, :].sum() == 0
            if f["bob_capture"]:
                assert c[:, :, :, 2].sum() == 0
            else:
                assert c[:, :, :, :2].sum() == 0
        for name, stats in ex.scores(merged).items():
            assert stats == study["scores"][name]
    free = record["free_pair"]
    score = free["held"]["scores"]["joint_capture"]
    assert score["S"] - score["simultaneous_99_radius"] > 2
    assert ex.all_axes_bound(free["axis_training"])["simultaneous_99_upper_bound"] < 2
    assert abs(free["unread_rails"]["scores"]["joint_capture"]["S"]) < 2


def test_apertures_account_for_absorption():
    record = json.loads((ex.HERE / "results.json").read_text())["interference"]
    screens = record["screen_counts_then_absorbed"]
    for counts in screens.values():
        assert len(counts) == 241 and sum(counts) == record["shots_per_setting"]
    s = record["source"]
    assert screens["phase_0"][s] > 1.8 * (
        screens["first_only"][s] + screens["second_only"][s]
    )
    assert screens["phase_4"][s] < 0.01 * screens["phase_0"][s]


def test_midpoint_includes_all_signed_publications_and_all_capture_outcomes():
    record = json.loads((ex.HERE / "results.json").read_text())
    midpoint = record["free_pair"]["midpoint_publication"]
    assert {tuple(r["gateway_signs"]) for r in midpoint["publications"]} == set(
        product((-1, 0, 1), repeat=2)
    )
    totals = np.zeros((4, 2, 2, 3, 3), dtype=int)
    for publication in midpoint["publications"]:
        flags = np.array(publication["flag_xy_local_counts"])
        assert flags.sum() == publication["launches"]
        totals += flags
    assert totals.sum() == ex.N
    reconstructed = totals.sum(axis=0)
    np.testing.assert_array_equal(reconstructed, midpoint["signed_counts"])
    assert ex.scores(reconstructed) == midpoint["scores"]
    optimized = midpoint["pooled_recalibration"]
    assert ex.capture_axes_bound(optimized["training"]) == optimized["bound"]
    scores = optimized["held"]["scores"]["joint_capture"]
    assert (
        scores["S"] - scores["simultaneous_99_radius"]
        < optimized["bound"]["simultaneous_99_upper"]
    )
