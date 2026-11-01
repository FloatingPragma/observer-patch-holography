"""Complete energy transfer, spectral redistribution and dark-access controls."""
import numpy as np
import pytest
from scipy import sparse

from ..charged_defect import dictionary
from ..common_spin_protection import spectrum, spin
from ..full_fock_field.field import FieldGeometry, FullFock
from . import response as r, algebra


def test_two_dimer_exact_work_spectrum_and_positive_negative_state_controls():
    field = FullFock(8)
    number = sparse.diags([i.bit_count() for i in range(256)], dtype=float)
    k = (field.hamiltonian(FieldGeometry((4, 1, 1)))-number/2).toarray()
    outcomes = dictionary.full_extension(2)
    w = sum((op.conj().T@k@op for _, op in outcomes), start=np.zeros_like(k))-k
    energies, vectors = np.linalg.eigh(w)
    expected = sorted([-np.sqrt(3/512), np.sqrt(3/512)]+[-np.sqrt(1/512)]*3+
                      [np.sqrt(1/512)]*3+[0.]*248)
    np.testing.assert_allclose(energies, expected, atol=2e-15, rtol=0)
    assert abs(np.vdot(w, w).real-3/128) < 1e-15
    for j in (0, -1):
        state = vectors[:, j]
        before = np.vdot(state, k@state).real
        after = sum(np.vdot(op@state, k@(op@state)).real for _, op in outcomes)
        assert abs(after-before-energies[j]) < 2e-15
    result = r.complete_work_identity()
    assert result['transfer_identified_as_heat_or_controller_work'] is False


def test_unital_identity_is_not_claimed_for_an_arbitrary_nonunital_channel():
    # Complete amplitude damping obeys the symmetrized W identity but lacks
    # the unital trace cancellation. This exercises the missing hypothesis.
    k = np.diag([0., 1.]); gamma = .4
    operators = [np.diag([1., np.sqrt(1-gamma)]), np.array([[0., np.sqrt(gamma)], [0., 0.]])]
    w = sum(x.T@k@x for x in operators)-k
    symmetric = sum((x.T@(k@x-x@k)+(x.T@k-k@x.T)@x)/2 for x in operators)
    np.testing.assert_allclose(w, symmetric, atol=2e-16, rtol=0)
    quadratic = sum(np.linalg.norm(k@x-x@k)**2 for x in operators)
    assert abs(np.trace(k@w)+quadratic/2) > .1


def test_complete_active_transition_law_has_resolved_energy_transfer():
    value = r.active_response()
    matrix = np.array(value['complete_energy_transition_matrix'])
    outside = np.array(value['outside_energy_transition_matrix'])
    assert matrix.shape == outside.shape == (90, 90)
    np.testing.assert_allclose(matrix.sum(axis=0), 1, atol=2e-14, rtol=0)
    np.testing.assert_allclose(matrix.sum(axis=1), 1, atol=2e-14, rtol=0)
    assert np.min(value['mean_squared_energy_change_per_input']) > 1e-7
    assert value['lowest_energy_input_mean_transfer'] > .06
    assert np.min(np.diag(outside)) > .1
    # Direct complete Kraus expectation, independent of displayed transition
    # weights, binds the reported energy transfer on the same lowest state.
    columns, _ = r.active_chart(); state = columns[:, 0]
    k = spin.operators().K.matrix(16)
    before = np.vdot(state, k@state).real
    after = sum(np.vdot(op@state, k@(op@state)).real
                for _, op in spin.repair_outcomes('hole', 'E', (0, 1)))
    assert abs(after-before-value['lowest_energy_input_mean_transfer']) < 2e-14


def test_dark_coherent_response_survives_but_simple_ground_probes_are_absent():
    value = r.dark_response()
    assert value['nonzero_beat_frequency'] > .05
    assert value['full_state_reference_error'] < 1e-13
    assert value['native_outcome_probability'] < 1e-24
    records = value['complete_stopping_history_outputs']
    assert len(records) == 41
    assert sum(row['probability'] for row in records) == pytest.approx(1, abs=1e-14)
    assert value['simple_two_annihilator_probes_checked'] == 66
    assert value['largest_ground_probe_dark_overlap'] < 1e-13
    assert value['dark_state_preparation_and_detector_source_executed'] is False
    # A direct supplied eta-raised two-up input is a positive access control;
    # the exclusion concerns simple pair removal from the eta-singlet ground.
    d, _ = algebra.dark_columns(); state = np.asarray(d[:, 0].toarray()).ravel()/2
    assert np.linalg.norm(d.conj().T@state/2) == pytest.approx(1, abs=1e-15)
