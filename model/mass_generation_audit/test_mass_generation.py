from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
import sympy as sy

from . import chain, gap, scaling


def test_original_source_is_bound_without_new_hopping_or_interaction_assignments():
    result = chain.certificate()
    assert result['complete_Fock_witness_dimension'] == 256
    assert result['hopping_after_exact_frame'] == 't=c/2=1/16'
    assert result['contact_changed'] is False
    assert result['native_reference_log_statistics_used_as_new_couplings'] is False
    np.testing.assert_array_equal(chain.frame(1000001), chain.frame(1))


@pytest.mark.parametrize('length', [3, 6, 10, 68, True])
def test_periodic_frame_mismatch_is_not_silently_dropped(length):
    with pytest.raises(ValueError): chain.one_body(length)


@pytest.mark.parametrize('power', [1, 3, 9, 23])
def test_euler_enclosure_against_independent_special_function(power):
    interval = gap.dirichlet_eta(power)
    with mp.workdps(85):
        value = mp.log(2) if power == 1 else (1-mp.mpf(2)**(1-power))*mp.zeta(power)
        lo = mp.mpf(interval.lower.numerator)/interval.lower.denominator
        hi = mp.mpf(interval.upper.numerator)/interval.upper.denominator
        assert lo <= value <= hi


def test_new_strong_gap_enclosure_matches_original_positive_integral():
    interval = gap.strong_gap()
    with mp.workdps(80):
        t, u = mp.mpf(1)/16, mp.mpf(1)
        a = 2*mp.pi*t/u
        integrand = lambda y: 2*mp.sqrt(y*y-1)*mp.exp(-a*y)/(1-mp.exp(-2*a*y))
        value = 16*t*t/u*mp.quad(integrand, [1, 2, 4, 8, 16, 32, mp.inf])
        lo = mp.mpf(interval.lower.numerator)/interval.lower.denominator
        hi = mp.mpf(interval.upper.numerator)/interval.upper.denominator
        assert lo < value < hi
    assert interval.upper-interval.lower < F(1, 10**24)


def test_absolute_tail_bound_tightens_and_rejects_wrong_domain():
    coarse, fine = gap.strong_gap(order=6), gap.strong_gap(order=12)
    assert coarse.lower < fine.lower < fine.upper < coarse.upper
    with pytest.raises(ValueError): gap.strong_gap(F(1), F(1))


def test_original_mu_zero_reference_is_not_mislabeled_a_vacuum():
    result = gap.certificate()
    assert F(result['original_H_removal_energy_interval']['upper_exact']) < 0
    assert F(result['original_H_addition_energy_interval']['lower_exact']) > 0
    assert result['three_plus_one_mass_or_species_identified'] is False


def test_global_bessel_envelope_contains_tight_source_gap():
    coarse, _, _ = scaling.gap_envelope(chain.HOPPING, chain.REPULSION)
    exact = gap.strong_gap()
    assert coarse.lower < exact.lower < exact.upper < coarse.upper
    rows = scaling.certificate()['explicit_fixed_U_enclosures']
    assert F(rows[-1]['charge_gap_envelope']['upper_exact']) < F(1, 10**40)


def test_running_condition_is_an_identity_not_a_mass_fit():
    A, v, M = sy.symbols('A v M', positive=True)
    a = 4*v*sy.exp(-A)/(M*sy.sqrt(2*sy.pi*A))
    t, U = v/(2*a), sy.pi*v/(a*A)
    assert sy.simplify((4/sy.pi)*sy.sqrt(U*t)*sy.exp(-2*sy.pi*t/U)-M) == 0
    assert scaling.certificate()['target_M_determined_by_fixed_source'] is False
