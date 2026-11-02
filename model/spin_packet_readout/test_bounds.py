from fractions import Fraction as F
import pytest
from . import bounds as b


def test_success_is_p_to_four_and_not_p_to_one():
    x=b.packet_success(F(1),F(2));y=b.packet_success(F(1,2),F(2))
    assert y['ideal_success_lower']==x['ideal_success_lower']/16
    assert y['time_LCU_normalization']==29
    assert y['independent_fresh_input_attempts_mean_upper']*y['ideal_success_lower']>=1


def test_zero_success_and_excess_error_cannot_be_conditioned():
    for args in ((0,0,0),(F(1,4),F(1,2),0),(F(1,4),0,F(1,4))):
        with pytest.raises(ValueError):b.conditioned_errors(*args)
    x=b.conditioned_errors(F(1,4),F(1,100),F(1,10000))
    assert 0<x['actual_success_lower']<F(1,4)
    assert x['conditional_trace_error_upper']<F(1,10)


def test_quadrature_has_a_complete_budget():
    x=b.quadrature_plan(1,2,F(1,4),F(1,2),F(1,100))
    assert x['total_time_filter_error_upper']<=F(1,100)
    assert x['coefficient_L1_upper']<5
    assert x['cells']*x['step']==2*x['time_cutoff']


def test_native_plan_counts_control_overhead_and_precision():
    x=b.native_plan(2,1,0,2,5,F(1),F(1),F(1,4),F(1,100),F(1,1000))
    assert x['padded_label_dimension']==8
    assert x['maximum_conjugation_product_error']<=F(1,100)
    assert x['native_complete_isometry_error_upper']<=F(1,1000)
    assert not x['input_ground_preparation_cost_included']
    assert not x['literal_large_tape_executed']


def test_logarithmic_lr_certificate_handles_long_time_without_exponential_allocation():
    x=b.locality_radius(3,1,10**8,F(1,100))
    r=x['radius']
    assert r-3>=x['time_exponent_upper']+x['log2_prefactor_over_error_upper']+2*b.ceil_log2(F(2*r+1))
    assert r<10**11


@pytest.mark.parametrize('bad',[False,1.5,-1])
def test_native_counts_are_exact_nonnegative_integers(bad):
    with pytest.raises(ValueError):b.native_plan(bad,1,0,2,5,1,1,F(1,4),F(1,100),F(1,1000))
