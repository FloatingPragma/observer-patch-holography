"""Exact arithmetic and domain controls supporting the dense-effect proof."""
import pytest
import sympy as sp

from . import dense_effects as dense


def test_exact_fixture_and_explicit_source_boundary():
    result = dense.build()
    assert result['regularization']['operator_error_bound'] == '1/51'
    assert result['finite_tomography_counterexample']['additivity_defect'] == '-9/64'
    assert not result['source_dense_effect_menu_realized']
    assert not result['source_equal_effect_record_equivalence_established']
    assert result['canonical_physical_premises_discharged'] == []


@pytest.mark.parametrize('matrix,expected', [
    ([[0, 0], [0, 0]], True),
    ([[0, 0], [0, 1]], True),
    ([[0, 1], [1, 0]], False),
    ([[1, sp.I], [-sp.I, 1]], True),
    ([[1, 2], [2, 1]], False),
    ([[-1, 0], [0, 1]], False),
])
def test_exact_schur_elimination_handles_null_and_negative_directions(matrix, expected):
    assert dense.positive_semidefinite(matrix) is expected


@pytest.mark.parametrize('bad', [0.25, True, sp.sqrt(2), sp.I])
def test_no_float_boolean_or_irrational_radius_is_silently_rationalized(bad):
    with pytest.raises(ValueError):
        dense.rational(bad)


def test_nonhermitian_and_noneffect_inputs_rejected():
    with pytest.raises(ValueError):
        dense.hermitian([[0, 1], [0, 1]])
    with pytest.raises(ValueError):
        dense.effect(sp.eye(2) * 2)
    with pytest.raises(ValueError):
        dense.effect([[sp.Rational(1, 2), 0.0], [0, sp.Rational(1, 2)]])


def test_lipschitz_certificate_keeps_two_noncommuting_order_comparisons():
    first = sp.diag(1, 0)
    second = sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2)], [sp.Rational(1, 2), sp.Rational(1, 2)]])
    assert first * second != second * first
    forward, reverse = dense.lipschitz_order_witness(first, second, sp.Rational(3, 4))
    for left, right in (forward, reverse):
        assert dense.positive_semidefinite(right - left)
    with pytest.raises(ValueError):
        dense.lipschitz_order_witness(first, second, sp.Rational(1, 2))


def test_joint_regularization_preserves_the_sum_at_a_boundary():
    delta = sp.Rational(1, 10)
    first = sp.diag(-delta, 1 + delta)
    second = sp.eye(2) - first
    with pytest.raises(ValueError):
        dense.effect(first)
    left, right = dense.regularize_coexistent(first, second, delta)
    assert dense.effect(left + right) == (1 + 2 * delta) / (1 + 4 * delta) * sp.eye(2)
    with pytest.raises(ValueError):
        dense.regularize_coexistent(sp.eye(2), sp.eye(2), delta)


def test_complex_three_dimensional_reconstruction_has_correct_imaginary_sign():
    factor = sp.Matrix([[1, sp.I, 0], [0, 1, 1], [1, 0, 1]])
    density = factor * factor.H / sp.trace(factor * factor.H)
    result = dense.tomographic_candidate(lambda matrix: sp.simplify(sp.trace(density * matrix)), 3)
    assert result == density
    assert dense.positive_semidefinite(result)
    assert result[0, 1] == sp.I / 6


def test_candidate_tomography_does_not_certify_a_nonlinear_assignment():
    cubic = lambda t: 3 * t - 6 * t**2 + 4 * t**3
    assignment = lambda matrix: cubic(sp.trace(matrix) / 2)
    assert dense.tomographic_candidate(assignment, 2) == sp.eye(2) / 2
    for matrix in (sp.diag(1, 0), sp.diag(0, 1), sp.Matrix([[1, sp.I], [-sp.I, 1]]) / 2):
        assert assignment(matrix) == sp.Rational(1, 2)
        assert assignment(sp.eye(2) - matrix) == 1 - assignment(matrix)
    part = sp.eye(2) / 8
    assert dense.effect(2 * part) == 2 * part
    assert assignment(2 * part) != 2 * assignment(part)
