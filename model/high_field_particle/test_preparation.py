from fractions import Fraction as F

import numpy as np
import pytest

from model.native_ordered_ground import source
from model.native_spin_observer import algebra as spin_algebra
from . import algebra, preparation


@pytest.mark.parametrize('sites,edges', [(2, ((0, 1),)),
    (3, ((0, 1), (1, 2))), (3, ((0, 1), (0, 1), (1, 2))),
    (8, source.cube()[0])])
def test_product_residual_from_full_source_columns(sites, edges):
    j, b = F(1), F(400)
    actual, scale = algebra.transformed_source(sites, edges, j, b)
    result = preparation.product_bound(sites, edges, j, b)
    assert F(int(actual[0, 0]), scale) == result['Rayleigh_value']
    residual = sum(F(int(value)**2, scale**2) for value in actual[1:, 0])
    assert residual == result['residual_norm_squared']


def test_cube_short_bound_is_independent_of_eigensolver():
    result = preparation.product_bound(8, source.cube()[0], F(1), F(400))
    assert result['residual_norm_squared'] == 3
    assert result['excited_spectrum_minus_Rayleigh_lower'] == 394
    assert result['product_ground_half_trace_upper'] == F(1, 197)
    assert result['ideal_local_aligned_pointer_success_lower'] == F(155233, 155236)
    assert not result['uses_a_fitted_eigenvector']


def test_nonseparated_product_reference_rejected():
    with pytest.raises(ValueError, match='separation'):
        preparation.product_bound(2, ((0, 1),), F(2), F(1))


def test_actual_complete_native_short_preparation_and_creation():
    value = preparation.execution_certificate()
    assert value['actual_native_total_pulses'] == 167
    assert value['complete_work_rows'] == 131072
    assert value['native_first_failure_slots'] == 608
    assert value['complete_native_isometry_error_bound'] < F(1, 10**13)
    assert value['ideal_exact_ground_success_probability_lower'] == F(155233, 155236)
    assert not value['ideal_exact_ground_success_is_particle_band_weight']
    assert value['full_output_half_trace_error_from_exact_ground_instrument_upper'] < F(1, 100)
    assert value['all_native_first_failures_and_both_pointer_outcomes_retained']
    assert not value['success_renormalized']
    assert not value['literal_pointer_equals_abstract_out_of_code_extension']


@pytest.mark.parametrize('sites', [1, 2, 3])
def test_target_bit_addresses_equal_the_owning_encoding(sites):
    inclusion = spin_algebra.inclusion(sites)
    for mask in range(1 << sites):
        expected = np.flatnonzero(inclusion[:, mask])
        assert expected.tolist() == [preparation._encoded_address(mask, sites)]
        assert inclusion[expected[0], mask] == 1
    assert preparation._encoded_address(0, sites) == sum(2 << (2*x) for x in range(sites))


def test_cube_reversed_encoding_has_the_wrong_exact_phase():
    odd = sum(1 << x for x in source.cube()[1])
    correct, reversed_encoding = {}, {}
    for mask in range(256):
        amplitude = F((-1)**((mask & odd).bit_count()+(mask & 1)), 16)
        correct[preparation._encoded_address(mask, 8)] = amplitude
        reversed_encoding[sum(1 << (2*x+((mask >> x) & 1)) for x in range(8))] = amplitude
    # Fix the source basis convention; never remove the error with a phase fit.
    assert correct != reversed_encoding
    assert all(reversed_encoding[address] == -amplitude for address, amplitude in correct.items())
    assert sum((correct[address]-reversed_encoding[address])**2 for address in correct) == 4


def test_cached_execution_evidence_is_detached_including_fractions():
    value = preparation.execution_certificate()
    rational = value['finite_product_ground_bound']['excited_weight_upper']
    before = (rational.numerator, rational.denominator)
    rational._numerator = 999  # deepcopy(Fraction) alone does not detach its slots.
    value['finite_product_ground_bound']['full_finite_ground_simple'] = False
    value['ideal_exact_ground_success_probability_lower']._denominator = 7
    fresh = preparation.execution_certificate()
    actual = fresh['finite_product_ground_bound']['excited_weight_upper']
    assert (actual.numerator, actual.denominator) == before
    assert fresh['finite_product_ground_bound']['full_finite_ground_simple'] is True
    assert fresh['ideal_exact_ground_success_probability_lower'] == F(155233, 155236)
