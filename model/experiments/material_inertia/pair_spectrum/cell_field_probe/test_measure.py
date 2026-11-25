from copy import deepcopy
import json
import sys

import numpy as np
import pytest

from . import measure


def test_complete_native_pointer_decoder_and_failed_echo(monkeypatch):
    decoded = measure.pointer_decoder()
    assert decoded["local_words"] == decoded["distinct_complete_records"] == 157
    assert decoded["charge_outcomes"] == [-1, 0, 1]
    law = measure.candidate.native_law
    wrong = law.echo().copy()
    wrong[law.N_SRC, 0] += .25  # corrupt an observed written-memory column
    monkeypatch.setattr(law, "echo", lambda: wrong)
    with pytest.raises(ArithmeticError, match="Native signed pointer"):
        measure.pointer_decoder()


def test_coherent_finite_probe_unitary_gives_the_record_characteristic_function():
    vector = np.array([1, 1j, 2-1j])/np.sqrt(7)
    rho = np.outer(vector, vector.conj())
    reference = np.ones((2, 2))/2
    sx = np.array([[0, 1], [1, 0]])
    sy = np.array([[0, -1j], [1j, 0]])
    for phase in measure.PHASES:
        phase_by_charge = np.exp(1j*phase*np.array([-1, 0, 1]))
        unitary = np.diag(np.ravel([[1, p] for p in phase_by_charge]))
        evolved = unitary @ np.kron(rho, reference) @ unitary.conj().T
        reduced = np.einsum("aiaj->ij", evolved.reshape(3, 2, 3, 2))
        expected = np.diag(rho) @ phase_by_charge
        assert abs(np.trace(reduced@sx)-expected.real) < 2e-15
        assert abs(np.trace(reduced@sy)-expected.imag) < 2e-15


def test_physically_valid_matched_marginals_can_fail_spatial_neutrality():
    n = 24
    x = np.repeat(np.arange(n), n-1)
    offset = np.tile(np.arange(n-1), n)
    y = offset+(offset >= x)
    weights = np.zeros(len(x))
    weights[(x == 0) & (y == 1)] = .5
    weights[(x == 12) & (y == 13)] = .5
    null, report = measure.hard_core_control(weights, x, y, n)
    assert report["same_port_mass"] == 0 and report["maximum_marginal_residual"] < 1e-14
    actual, changed = (measure.cell_records(w, x, y, .25) for w in (weights, null))
    assert actual["different_cell_probability"] == 0
    assert abs(changed["different_cell_probability"]-.5) < 1e-14
    assert actual["summed_cell_charge_second_moment"] == 0
    for a, b in zip(actual["cells"], changed["cells"], strict=True):
        ma = a["minus_zero_plus"][2]-a["minus_zero_plus"][0]
        mb = b["minus_zero_plus"][2]-b["minus_zero_plus"][0]
        assert abs(ma-mb) < 1e-14
    # Same signed mean, different X response: a joint correlation is measured.
    assert actual["cells"][0]["supplied_Ramsey"][0]["X_plus_minus"] == [1., 0.]
    assert abs(changed["cells"][0]["supplied_Ramsey"][0]["X_plus_minus"][0]-.75) < 1e-14


def test_archived_complete_native_and_supplied_probe_partitions():
    receipt = json.loads(measure.ARCHIVE.read_text())
    assert len(receipt["launches"]) == 3
    for launch in receipt["launches"]:
        assert launch["gates"]["native_different_cells_below_001"]
        assert launch["gates"]["exchange_only_different_cells_above_002"]
        assert launch["matched_control"]["maximum_marginal_residual"] < 2e-14
        assert launch["gates"]["matched_model_effect_separation"] > .29
        for panel in launch["panels"].values():
            assert abs(panel["summed_cell_charge_second_moment"]-2*panel["different_cell_probability"]) < 2e-12
            for cell in panel["cells"]:
                partitions = [cell["minus_zero_plus"], cell["original_start_minus_zero_plus_filter_miss"]]
                for reading in cell["supplied_Ramsey"]:
                    partitions += [reading[key] for key in ("X_plus_minus", "Y_plus_minus",
                        "original_start_X_plus_minus_filter_miss", "original_start_Y_plus_minus_filter_miss")]
                for values in partitions:
                    assert min(values) >= 0 and abs(sum(values)-1) < 2e-12


def test_serialization_failure_creates_no_partial_record(tmp_path, monkeypatch):
    destination = tmp_path/"receipt.json"
    monkeypatch.setattr(measure, "run", lambda: {"bad_numpy_scalar": np.bool_(True)})
    monkeypatch.setattr(sys, "argv", ["measure", "--record", str(destination)])
    with pytest.raises(TypeError):
        measure.main()
    assert not destination.exists()


def test_checker_rejects_an_outcome_change_without_modifying_the_receipt():
    original = measure.ARCHIVE.read_bytes()
    actual = json.loads(original)
    changed = deepcopy(actual)
    changed["launches"][0]["panels"]["native_arrival"]["cells"][0]["minus_zero_plus"].reverse()
    with pytest.raises(ArithmeticError):
        measure.compare(actual, changed)
    assert measure.ARCHIVE.read_bytes() == original
