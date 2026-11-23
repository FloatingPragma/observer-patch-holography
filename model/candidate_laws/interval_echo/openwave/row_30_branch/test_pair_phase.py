from copy import deepcopy

import numpy as np
import pytest

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import federation
from tools.strict_json import load
from model.candidate_laws.interval_echo.openwave.row_30_branch import pair_phase as phase
from model.candidate_laws.interval_echo.openwave.row_30_branch import pair_phase_check as check
from model.candidate_laws.interval_echo.openwave.row_30_branch.pair_crossing import crossing


def test_only_pair_coherence_carries_the_sign_and_reader_detects_it():
    sector = PairSector(federation(0))
    inputs, analyzer = phase.apparatus()
    u, v, w, x = phase.rails(0)
    states = np.zeros((sector.P, 2), complex)
    states[[sector.index(u, w), sector.index(v, x)]] = inputs
    # Reduced one-load density (trace 1 after source admission).
    waves = np.zeros((sector.n, sector.n, 2), complex)
    waves[sector.aa, sector.bb] = waves[sector.bb, sector.aa] = states
    first, second = waves.T
    np.testing.assert_allclose(first @ first.conj().T, second @ second.conj().T, atol=1e-13)
    data = phase.observe(sector, 0, states, np.zeros(2), analyzer)
    np.testing.assert_allclose(data, [[0.5, 0, 0, 0, 0, 0.5], [0.5, 0, 0, 0, 0.5, 0]], atol=1e-13)
    src, dest = phase.dephased_controls(data)
    np.testing.assert_allclose(src, dest, atol=1e-13)
    np.testing.assert_allclose(src[:, -2:], 0.25, atol=1e-13)
    assert phase.observe(sector, 1, states, np.zeros(2), analyzer) == [[0.5, 0.5, 0, 0, 0, 0]] * 2


def test_routing_forecast_and_changed_receipts(monkeypatch):
    monkeypatch.setattr(phase, "SOURCES", (2,))
    monkeypatch.setattr(phase, "TIMES", (0, 16))
    actual = phase.forecast()
    archive = load(phase.PREDICTIONS)
    for field in ("receiver", "association"):
        archive[field] = [row for row in archive[field] if row[0] == 2 and row[-3] in (0, 16)]
    archive["all_time_adjacent_given_launch"] = {
        "2": archive["all_time_adjacent_given_launch"]["2"]
    }
    phase.verify_forecast(actual, archive)
    for bad in (True, float("nan"), -1e-12):
        archive["receiver"][0][-1][1] = bad
        with pytest.raises(ValueError, match="bin layout"):
            phase.verify_forecast(actual, archive)
    archive["receiver"][0] = deepcopy(actual["receiver"][0])
    archive["receiver"][0][-1][4:] = archive["receiver"][0][-1][4:][::-1]
    with pytest.raises(AssertionError):
        phase.verify_forecast(actual, archive)
    archive["receiver"] = archive["receiver"][1:]
    with pytest.raises(ValueError):
        phase.verify_forecast(actual, archive)


def test_receiver_refusals_and_binding_crossings():
    assert crossing() == {
        "adjacent_components": {30: 20, 1: 30},
        "ordered_four_move_crossings": 1500,
        "paths_per_crossing": 2,
    }
    n = check.SHOTS
    # Differential admission alone gives the same conditional readout.
    counts = np.array([[0, n - 200, 0, 0, 100, 100], [0, n - 20000, 0, 0, 10000, 10000]])
    assert check.decision(counts, counts / n) == (1.0, False)
    assert check.admission_floor(counts / n) == 0
    counts[1] = [0, n, 0, 0, 0, 0]
    assert check.decision(counts, counts / n) == (1.0, False)
    with pytest.raises(ValueError):
        check.decision(counts[:, 1:], counts / n)


def test_frozen_held_out_counts_and_verdict():
    forecast = {
        (s, d): np.array(rows)
        for s, d, law, t, *rows in load(phase.PREDICTIONS)["receiver"]
        if law == "native" and t == 256 and s != d and s in (2, 3)
    }
    trials = load(check.COUNTS)
    routes = [(r["source"], r["receiver"]) for r in trials]
    assert len(routes) == len(set(routes)) == 6 and set(routes) == forecast.keys()
    for key, row in zip(routes, trials, strict=True):
        p, passed = check.decision(np.array(row["counts"]), forecast[key])
        assert np.isclose(p, row["p_value"], rtol=1e-9, atol=1e-300) and passed == row["passed"]
    exchanged = deepcopy(trials)
    exchanged[0]["counts"][0][-2:] = exchanged[0]["counts"][0][-2:][::-1]
    with pytest.raises(ValueError, match="counts"):
        check.check_trials(trials, exchanged)
