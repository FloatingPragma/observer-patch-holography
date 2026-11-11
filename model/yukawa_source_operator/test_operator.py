"""Adversarial operator admission and independent normalization checks."""
import numpy as np
import pytest
import sympy as s

from . import operator


def coefficients():
    return {row['origin']: s.Matrix([[1, 2+s.I, 0], [3*s.I, 4, 5], [6, 0, 7-s.I]])
            for row in operator.tensor_lines()}


def test_complex_orientation_and_complete_components():
    original = coefficients()
    derivative = operator.assemble_derivatives(original)
    recovered = operator.extract_coefficients(derivative)
    assert recovered == original
    assert all(value != value.T and value != value.H for value in recovered.values())
    # An entry on a channel is recovered directly, independently of the dual.
    for row in operator.tensor_lines():
        (a, b), weight = next(iter(row['directed'][1].todok().items()))
        assert derivative[1][3*a, 3*b+1]/weight == original[row['origin']][0, 1]


def test_off_channel_kernel_and_inconsistent_Higgs_component_rejected():
    derivative = list(operator.assemble_derivatives(coefficients()))
    bad = s.MutableSparseMatrix(derivative[1])
    bad[0, 0] += 1
    with pytest.raises(ValueError, match='outside'):
        operator.extract_coefficients((derivative[0], bad))
    with pytest.raises(ValueError, match='outside'):
        operator.extract_coefficients((2*derivative[0], derivative[1]))
    with pytest.raises(ValueError):
        operator.extract_coefficients((s.eye(3), s.eye(3)))


def test_canonical_singular_values_are_invariant_under_complex_field_coordinates():
    y = np.array([[1, 2j, 3], [0, 2, 1-1j], [2, 0, 4]], complex)
    left = np.array([[2, 1j, 0], [0, 3, 1], [1j, 0, 2]], complex)
    right = np.array([[1, 0, 1j], [1, 2, 0], [0, 1j, 3]], complex)
    hscale = 7.
    changed = hscale*left.T@y@right
    canonical = operator.canonical_coefficient(changed, left.conj().T@left,
                                               right.conj().T@right, hscale*hscale)
    np.testing.assert_allclose(np.linalg.svd(canonical, compute_uv=False),
                               np.linalg.svd(y, compute_uv=False), rtol=2e-13, atol=2e-13)


def test_kinetic_metric_removes_spurious_family_splitting():
    gram = np.array([[2., 1, 0], [1, 3, 1], [0, 1, 4]])
    result = operator.canonical_coefficient(3*gram, gram, gram, 4)
    np.testing.assert_allclose(result, 1.5*np.eye(3), atol=3e-15)
    with pytest.raises(ValueError):
        operator.canonical_coefficient(gram, np.diag([1, 1, 0]), gram, 1)
    with pytest.raises(ValueError):
        operator.canonical_coefficient(gram, gram, gram, 0)


def test_exact_symbolic_certificate():
    report = operator.certificate()
    assert report['all_27_symbolic_coordinates_recovered_exactly']
    assert sorted(x['directed_rank'] for x in report['tensor_lines']) == [1, 3, 3]


def test_composite_overlap_does_not_supply_a_free_action_Yukawa_vertex():
    result = operator.composite_three_point()
    assert result['free_action_Higgs_derivative_exact'] == '0'
    assert sorted({x['pair_norm_squared_exact'] for x in result['rows']}) == ['2', '6']
    assert all(x['all_435_CAR_pair_contractions_checked'] for x in result['rows'])


def test_existing_trial_action_derivative_is_recovered_with_its_source_ancestry():
    result = operator.declared_action_coefficients()
    assert result['all_three_native_duals_recover_existing_derivative_exactly']
    for row in result['rows']:
        assert len(row['unit_Higgs_coordinate_singular_values_exact']) == 3
        assert float(s.sympify(row['largest_over_smallest_exact'])) < 1.5
    assert result['source_log_odds_exact']['a'] == 'log(25/8)'
    assert result['source_log_odds_exact']['b'] == 'log(8/5)'
