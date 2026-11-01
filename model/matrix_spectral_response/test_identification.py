from fractions import Fraction as F

import pytest

from . import identification as i


def test_minimax_bound_is_sharp_for_a_two_endpoint_example():
    # E=(1,3), target masses=(1,6):lambda=4/3 gives both errors1/3.
    bound=i.minimax_relative_error_lower(1,3,6)
    assert bound==F(1,3)
    scale=F(4,3)
    assert abs(scale-1)==abs(scale*3/6-1)==bound
    assert i.minimax_relative_error_lower(1,3,2)==0


def test_invalid_or_floating_dictionary_bounds_are_rejected():
    for values in ((0,1,2),(2,1,2),(1,2,0),(1.,2,3),(True,2,3)):
        with pytest.raises(ValueError):i.minimax_relative_error_lower(*values)


def test_entire_visible_gap_dictionary_fails_even_with_reported_uncertainties():
    result=i.certificate()
    assert F(result['source_spread_upper'])<62
    muon,tau=result['pair_tests']
    assert F(muon['mass_ratio_lower_on_reported_uncertainty_box'])>206
    assert F(tau['mass_ratio_lower_on_reported_uncertainty_box'])>3477
    assert F(muon['common_positive_scale_minimax_relative_error_lower'])>F(54,100)
    assert F(tau['common_positive_scale_minimax_relative_error_lower'])>F(96,100)
    assert result['all_78_down_lines_and_positive_up_line_included'] is True
    assert all(row['exact_assignment_possible'] is False for row in result['pair_tests'])


def test_rejection_does_not_promote_a_mass_claim_or_exclude_other_sectors():
    result=i.certificate()
    assert result['unit_factor_or_source_parameters_fitted'] is False
    assert result['finite_gaps_proved_to_be_relativistic_rest_masses'] is False
    assert result['unobserved_sectors_composites_or_continuum_excluded'] is False
    assert result['whole_OPH_model_falsified'] is False
    assert result['eligible_physical_comparisons_added']==0
