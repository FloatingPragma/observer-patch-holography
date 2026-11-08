from fractions import Fraction as F

import numpy as np
import pytest

from model.common_source import slow_geometry as geometry
from model.observer import native
from . import algebra


def _real(value):
    return float(value[0])+np.sqrt(5)*float(value[1])


def _matrices(branch):
    shift = np.zeros((78, 78), dtype=complex)
    for source, target in branch.pairs:
        shift[target, source] = 1
    coordinate = np.array([[_real(x) for x in geometry.load_readout(state)]
                           for state in native.states()])
    return coordinate, shift/np.sqrt(2)


def test_all_actual_columns_and_complete_adaptive_history():
    result = algebra.branch_certificate()
    assert result['source_columns'] == 4680
    assert result['branch_support_size_histogram'] == {'1': 60, '11': 60, '66': 60}
    assert algebra.complete_observer_certificate()['complete_history_columns'] == 312


def test_exact_uniform_mean_is_slow_contraction():
    assert algebra.mean_certificate()['diagonal_operator_entries_checked'] == 234
    # A fixed displacement cannot shrink opposite endpoint preparations.
    states = []
    for port in (0, 3):
        state = [0]*12
        state[port] = 2
        states.append(tuple(state))
    shifts = []
    for state in states:
        before = np.array([_real(x) for x in geometry.load_readout(state)])
        after = sum((np.array([_real(x) for x in geometry.load_readout(
            native.native(state, seam, coin)[0])])
            for seam in geometry.source().seams() for coin in (0, 1)), start=np.zeros(3))/60
        shifts.append(after-before)
    assert np.linalg.norm(shifts[0]) > .01
    np.testing.assert_allclose(shifts[0], -shifts[1], atol=2e-14)


def test_weyl_identity_and_reference_instrument_without_postselection():
    q = np.array([.37, -.41, .29])
    rng = np.random.default_rng(7021)
    psi = rng.normal(size=(78, 3))+1j*rng.normal(size=(78, 3))
    psi /= np.linalg.norm(psi)
    probability = 0
    for branch in algebra.branches(native.F):
        coordinate, kraus = _matrices(branch)
        phases = np.exp(1j*(coordinate@q))
        displacement = np.array([_real(x) for x in branch.displacement])
        np.testing.assert_allclose(phases[:, None]*kraus*phases.conj()[None, :],
                                   np.exp(1j*(q@displacement))*kraus, atol=2e-15)
        output = kraus@psi
        probability += np.linalg.norm(output)**2
        for i in range(3):
            np.testing.assert_allclose((coordinate[:, i]-displacement[i])[:, None]*output,
                                       kraus@(coordinate[:, i, None]*psi), atol=2e-15)
    assert abs(probability-1) < 2e-15


def test_record_acceptance_filter_cannot_be_removed():
    branch = next(r for r in algebra.branches(native.E) if r.coin == 0 and r.current == 1)
    _, kraus = _matrices(branch)
    inside = branch.pairs[0][0]
    outside = next(i for i in range(78) if i not in {a for a, _ in branch.pairs})
    assert abs(np.linalg.norm(kraus[:, inside])**2-.5) < 2e-16
    assert np.linalg.norm(kraus[:, outside]) == 0
    assert np.linalg.matrix_rank(kraus) == 11


def test_branch_quadrature_commutators_and_source_uncertainty():
    branch = next(r for r in algebra.branches(native.E) if r.coin == 0 and r.current == 1)
    coordinate, kraus = _matrices(branch)
    x = np.diag(coordinate[:, 0])
    real = (kraus+kraus.conj().T)/2
    imag = (kraus-kraus.conj().T)/(2j)
    delta = _real(branch.displacement[0])
    np.testing.assert_allclose(x@real-real@x, 1j*delta*imag, atol=1e-15)
    np.testing.assert_allclose(x@imag-imag@x, -1j*delta*real, atol=1e-15)
    source, target = branch.pairs[0]
    psi = np.zeros(78, dtype=complex)
    psi[source], psi[target] = 1/np.sqrt(2), 1j/np.sqrt(2)
    def variance(op):
        mean = np.vdot(psi, op@psi)
        return float(np.vdot(op@psi, op@psi).real-abs(mean)**2)
    rhs = abs(delta*np.vdot(psi, imag@psi))/2
    assert rhs > 0
    assert variance(x)*variance(real) >= rhs**2-1e-14
    # A source position eigenstate makes any universal positive bound false.
    psi[:] = 0
    psi[source] = 1
    assert variance(x) == 0
    assert variance(real) > 0


def test_nonzero_unitary_translation_has_full_operator_error():
    x = np.diag(np.array([-2., -.5, 1., 4.]))
    permutation = np.roll(np.eye(4), 1, axis=0)
    for shift in (F(1, 4), F(-3, 2), F(0)):
        residual = permutation.T@x@permutation-x-float(shift)*np.eye(4)
        assert np.linalg.norm(residual, ord=2) >= float(algebra.finite_translation_residual(shift))
        assert abs(np.trace(residual)+4*float(shift)) < 1e-15
        identity_residual = x-x-float(shift)*np.eye(4)
        assert np.linalg.norm(identity_residual, ord=2) == float(abs(shift))


def test_actual_protected_source_has_no_nonconstant_total_load():
    result = algebra.projected_load_certificate()
    assert result['onsite_occupation_entries_checked'] == 2048
    assert result['original_repair_certificate']['protected_cube_columns'] == 12288
    assert not result['same_total_load_detector_can_locate_protected_spin_excitations']


def test_dressing_quadratic_not_linear_in_charge_leakage():
    theta = .13
    inclusion = np.eye(4)[:, :2]
    charge = np.diag([0., 0., 1., -1.])
    dressed = np.cos(theta)*inclusion+np.sin(theta)*np.eye(4)[:, 2:]
    difference = dressed-inclusion
    np.testing.assert_allclose(dressed.T@charge@dressed, difference.T@charge@difference,
                               atol=2e-17)
    assert np.linalg.norm(dressed.T@charge@dressed, 2) <= np.linalg.norm(difference, 2)**2
    assert algebra.dressed_load_bound(F(3), 10) == F(3, 250000)
    assert algebra.dressed_load_bound(F(3), 20) == algebra.dressed_load_bound(F(3), 10)/64


@pytest.mark.parametrize('edge', [(0, 0), (1, 0), (True, 1), [0, 1], (0, 12)])
def test_invalid_native_seams_rejected(edge):
    with pytest.raises(ValueError):
        algebra.branches(edge)


@pytest.mark.parametrize('weight,edges', [(1, 2), (F(-1), 2), (F(1), 1), (F(1), True)])
def test_invalid_dressed_bound_inputs_rejected(weight, edges):
    with pytest.raises(ValueError):
        algebra.dressed_load_bound(weight, edges)


def test_boolean_shift_is_not_a_coordinate_number():
    with pytest.raises(ValueError):
        algebra.finite_translation_residual(True)
