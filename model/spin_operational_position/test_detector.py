from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as sp

from model.native_spin_observer import circuits as native, execute
from model.protected_particle_readout import execute as local
from . import bounds, circuits, experiment, generic, geometry


def test_exact_source_reflector_certificate_and_failure_control():
    value = geometry.certify()
    assert F(value['mode_isometry_error']) < F(1, 10**10)
    bad = [list(v) for v in geometry.witness()]
    bad[0][0] += 1 << 40
    with pytest.raises(ArithmeticError):
        geometry.certify(tuple(map(tuple, bad)))


def test_complex_elimination_conjugation_is_fixed():
    x = sp.Matrix([sp.I/2, sp.sqrt(3)/2])
    e = sp.Matrix([1, 0])
    v = x-e
    real, imaginary = generic.elimination_denominator(F(), F(1, 2), 1)
    denominator = sp.Rational(real.numerator, real.denominator)+sp.I*sp.Rational(
        imaginary.numerator, imaginary.denominator)
    unitary = sp.eye(2)-v*v.conjugate().T/denominator
    assert sp.simplify(unitary*x-e) == sp.zeros(2, 1)
    assert sp.simplify(unitary.conjugate().T*unitary-sp.eye(2)) == sp.zeros(2)
    wrong = sp.eye(2)-v*v.conjugate().T/sp.conjugate(denominator)
    assert sp.simplify(wrong*x-e) != sp.zeros(2, 1)


@pytest.mark.parametrize('normal', [((1, 2), (0, 0), (-3, 1), (2, -4)),
                                   ((0, -2), (-3, 0), (0, 4), (1, 0))])
def test_refinable_complex_preparation_and_rank_one_source_phase(normal):
    vector = np.array([complex(a, b) for a, b in normal])
    vector /= np.linalg.norm(vector)
    word, record = generic.preparation(normal, 28)
    result = execute.complete(list(native.emit(word)), 2, np.eye(4)[:, :1])
    assert np.linalg.norm(result.success[:, 0]-vector) < float(F(record[
        'preparation_unitary_operator_error_upper']))+1e-12
    word, record = generic.rank_one_phase(normal, F(3, 7), 28)
    result = execute.complete(list(native.emit(word)), 2, np.eye(4))
    target = np.eye(4)+(np.exp(3j/7)-1)*np.outer(vector, vector.conj())
    assert np.linalg.norm(result.success-target, 2) < float(F(record[
        'rank_one_operator_error_upper']))+1e-12
    assert record['pulses'] <= 20*4


def test_exact_rational_phase_diagonal_all_columns():
    phases = tuple(F(j*j-4*j, 9) for j in range(8))
    word = list(native.emit(generic.diagonal(phases)))
    result = execute.complete(word, 3, np.eye(8))
    target = np.diag(np.exp(1j*np.array(phases, float)))
    np.testing.assert_allclose(result.success, target, atol=2e-14)


@pytest.mark.parametrize('bits', [1, 2, 3])
def test_marked_phase_preserves_all_other_columns(bits):
    word = list(native.emit(circuits.diagonal_phase_zero(bits)))
    result = local.compile_block(local.key(word), bits)
    target = np.eye(1 << bits, dtype=complex)
    target[0, 0] = -1
    np.testing.assert_allclose(result.success, target, atol=2e-14)


@pytest.mark.parametrize('pattern', [0, 1, 2, 3])
def test_native_controlled_label_bit_all_columns(pattern):
    word = list(native.emit(circuits.multi_controlled_not((0, 1), 2, pattern)))
    result = execute.complete(word, 3, np.eye(8))
    target = np.eye(8)
    target[:, [pattern, pattern+4]] = target[:, [pattern+4, pattern]]
    np.testing.assert_allclose(result.success, target, atol=2e-14)


def test_complete_coordinate_and_dark_luders_target():
    unitary = geometry.matrix()
    np.testing.assert_allclose(unitary.T@unitary, np.eye(256), atol=2e-15)
    initial, target = experiment.input_and_target()
    np.testing.assert_allclose(np.linalg.norm(initial), 1, atol=2e-15)
    np.testing.assert_allclose(np.linalg.norm(target), 1, atol=3e-15)
    # Dark coherence remains a single vector, not 248 measured fine labels.
    np.testing.assert_allclose(target[:256, 1], initial[:256, 1], atol=2e-15)
    assert np.linalg.norm(target[256:2048]) == 0
    assert all(circuits.label_of(i) == (8+i if i < 8 else 0) for i in range(256))


def test_fixed_source_program_and_arbitrary_reference_error():
    core, physical, report = circuits.plan()
    assert len(core) == 49228 and len(physical) == 49324
    assert F(report['complete_actual_coordinate_halfdiamond_error']) < F(1, 10**6)
    assert report['dark_subspace_luders_backaction_retained']
    assert not report['full_twenty_bit_native_tape_executed']
    # Exposing a caller-owned word must not mutate the cached producer.
    core[0]['angle'] += 1
    assert circuits.plan()[0][0]['angle'] != core[0]['angle']


@pytest.mark.parametrize('spins,outcomes', [(1, 3), (8, 9), (24, 17)])
def test_finite_instrument_budget_covers_all_four_errors(spins, outcomes):
    error = F(1, 10**7)
    record = bounds.generic_instrument_plan(spins, outcomes, error)
    dimension = record['unitary_dimension']
    assert 4*dimension*F(record['normalized_factor_vector_error_budget']) <= error/8
    assert (dimension+1)*F(record['factor_phase_error_budget']) <= error/8
    assert F(record['rank_one_tree_operator_error_upper']) <= error/4
    assert F(record['port_operator_error_upper']) <= error/2
    assert record['native_pulse_upper_bound'] == 20*dimension**2+4*dimension+10*spins


def test_one_coarse_depolarizing_outcome_needs_four_kraus_branches():
    paulis = (np.eye(2), np.array([[0, 1], [1, 0]]),
              np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    kraus = tuple(p/2 for p in paulis)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(2))
    pure = np.diag([1, 0])
    np.testing.assert_allclose(sum(k@pure@k.conj().T for k in kraus), np.eye(2)/2)
    record = bounds.generic_instrument_plan(1, len(kraus), F(1, 100))
    assert record['retained_kraus_branch_count'] == 4
    assert record['coherent_record_bits'] == 2
    assert record['unitary_dimension'] == 8


@pytest.mark.parametrize('normal', [(), ((0, 0), (0, 0)), ((True, 0), (1, 0)),
                                  ((1, 0), (1, 0), (1, 0))])
def test_malformed_complex_normals_rejected(normal):
    with pytest.raises(ValueError):
        generic.preparation(normal, 12)


@pytest.mark.parametrize('args', [(True, 2, F(1, 10)), (1, 0, F(1, 10)),
                                 (1, 2, .1), (1, 2, F(0))])
def test_invalid_generic_budget_rejected(args):
    with pytest.raises(ValueError):
        bounds.generic_instrument_plan(*args)
