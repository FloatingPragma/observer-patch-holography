from fractions import Fraction as F
import math
from .bounds import exp_upper, infinite_flux_bound


def test_rational_exponential_upper_bounds():
    for x in (F(-25),F(-1,2),F(0),F(1,2),F(3)):
        upper = exp_upper(x)
        assert float(upper) >= math.exp(float(x))*(1-2e-15)
        assert float(upper) < math.exp(float(x))*(1+1e-12)


def test_complete_flux_bound_improves_under_cap_refinement():
    low = infinite_flux_bound(24,3,1,4)
    high = infinite_flux_bound(48,3,1,4)
    key = 'state_norm_error_before_numerical_diagonalization_exact'
    assert F(high[key]) < F(low[key])/1000
    assert F(high[key]) < F(1,10**8)
    assert high['finite_numerical_errors_included'] is False


def test_antiparticle_and_family_domains_retained():
    assert infinite_flux_bound(24,-3,1,4) == infinite_flux_bound(24,3,1,4)
    for m in (1,207,3477):
        assert F(infinite_flux_bound(24,3,m,4)['weighted_root_ratio_upper_exact']) < 1
