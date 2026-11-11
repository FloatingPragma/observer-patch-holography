from fractions import Fraction as F
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm

from . import spatial as q


def test_complete_spectral_contraction_includes_cross_input_coherence():
    result = q.contractions()
    assert len(result['literal_contractions']) == 9
    assert len(result['actual_source_routes']) == 6
    assert {r['total_effect'] for r in result['literal_contractions']} == {'1/128'}
    assert result['initial_energy_record_changes_the_experiment']
    # Recording the initial line deletes its off-diagonal input Gram terms.
    weights = (F(1,3),F(1,6))
    coherent = sum(weights)**3/16
    measured = sum(w*w for w in weights)*sum(weights)/16
    assert coherent == F(1,128) and measured < coherent


def test_all_literal_roots_shift_by_one_and_no_visible_resonance():
    for family, floor in (('uc',F('0.0749111786')),('dc',F('0.0749111786')),('Q',F('0.0511420752'))):
        data = q.family(family)
        assert F(data['minimum_input_output_detuning_lower']) > floor
        assert F(data['minimum_visible_input_to_one_link_full_H_bin_gap_lower']) > floor-F(1,10**5)
        assert data['disjoint_energy_control_is_operator_identity_not_extra_executed_setting']
        assert data['exact_visible_support_shift_gcd'] == '1'
        assert data['visible_atomic_lines'] == 13
        initial=q.parents()[0]['exact_local_response'][family]['exact_ground_moments']
        shifted=data['one_link_first_hop_exact_moments_in_neutral_root_field']
        assert shifted[0] == ['1/128']
        difference=[64*F(b)-F(a)for a,b in zip(initial[1],shifted[1])]
        assert difference == [F(0)]*(len(difference)-1)+[F(1,2)]
        for row in data['all_atomic_and_full_H_bins']:
            old = tuple(map(F,row['initial_local_energy']))
            shifted = tuple(map(F,row['one_link_atomic_local_energy']))
            # Outward serialization can add one unit in its last decimal.
            assert shifted[0] <= old[0]+1 <= old[1]+1 <= shifted[1]
            w = tuple(map(F,row['one_component_weight']))
            route = tuple(map(F,row['unfiltered_first_hop_weight']))
            assert route[0] <= w[0]/64 <= w[1]/64 <= route[1]
        zeros = [r for r in data['all_atomic_and_full_H_bins'] if r['exact_atomic_probe_zero']]
        assert len(zeros) == (12 if family == 'Q' else 0)
        for row in zeros:
            assert row['unfiltered_first_hop_weight'] == ['0','0']
            assert F(row['original_route_then_energy_probability'][1]) > 0


def test_rational_sine_bounds_and_literal_one_link_two_level_control():
    mp.mp.dps = 100
    low,high = q.sine_interval(F(1,2))
    exact = mp.sin(mp.mpf(1)/2)
    convert = lambda x:mp.mpf(x.numerator)/x.denominator
    assert convert(low) <= exact <= convert(high)
    # One initial line of weight1/2 gives the exact1/128 short-time moment.
    # The finite two-level control checks the time integral and state-norm factor.
    lam = .01; time = 1.; coupling = 1/8
    H0 = np.diag([0.,1.]); T = np.array([[0.,coupling],[coupling,0.]])
    psi = np.array([1/np.sqrt(2),0.])
    actual = expm(-1j*time*(H0+lam*T))@psi
    first = -lam*coupling*(1-np.exp(-1j*time))*psi[0]
    remainder = (lam*coupling*time)**2*np.linalg.norm(psi)/2
    assert abs(actual[1]-first) <= remainder
    kernel = lam**2/16 * .5 * abs(.5*(1-np.exp(-1j*time)))**2
    assert abs(kernel-abs(first)**2)<1e-20
    # Source energy is conserved before the arrival projection. Selecting
    # the second basis vector reads energy1 and requires measurement backaction.
    assert abs(np.vdot(actual,(H0+lam*T)@actual)-np.vdot(psi,(H0+lam*T)@psi))<1e-16
    assert (H0+lam*T)[1,1] == 1


def test_full_H_positive_ordered_records_and_complete_capacity():
    for family in ('uc','dc','Q'):
        data=q.family(family)
        n=data['retained_atomic_lines']
        best=data['all_atomic_and_full_H_bins'][data['best_certified_bin_index']]
        low,high=map(F,best['original_route_then_energy_probability'])
        assert 0 < low < high < F(1,10**21)
        assert F(data['complete_success_nonarrival_lower']) > F('0.4999')
        outside=tuple(map(F,data['original_route_then_outside_all_energy_bins_probability']))
        assert outside[0]==0 < outside[1] < F(data['full_H_original_route_probability_upper'])
        assert data['final_detector_pointer_dimension'] == 1+6*(n+1)
        assert data['coarse_original_attempt_records'] == 4+6*(n+1)
        assert not data['autonomous_energy_gain_or_simultaneous_joint_PVM_claimed']
        assert F(data['full_unitary_Dyson_remainder_per_sqrt_localization_upper']) < F(16,10**17)


def test_domains_and_exact_probability_normalization():
    import sympy as s
    qq,a,p=s.symbols('q a p')
    # Energy bins plus outside sum to each retained route probability.
    ps=s.symbols('p:6')
    assert s.expand(1-qq+qq*(1-a)+qq*a/2+qq*(a/2-sum(ps))+qq*sum(ps)) == 1
    for bad in (-1,0,3):
        with pytest.raises(ValueError):q.family('uc',F(bad))
    with pytest.raises(ValueError):q.sine_interval(F(2))


def test_no_unlisted_source_or_physical_scope_promotion():
    value=q.certificate()
    assert value['scope']['same_H_positive_energy_resolved_ordered_arrival']
    assert not value['scope']['all_motion_proved_virtual_or_globally_confined']
    assert not value['scope']['full_charged_sector_ground_or_string_tension_identified']
    assert value['apparatus']['background_cancelled_explicitly_during_services']
    assert value['apparatus']['additional_recorded_atomic_localization']
    assert not value['apparatus']['exact_interacting_ground_projection_added']
    assert value['apparatus']['additional_elapsed_source_time'].startswith('3 =')
    assert value['exact_idle_or_zero_hopping_arrival_probability'] == '0'
    assert q.PARENT in value['input_sha256']


def test_changed_parent_receipt_is_rejected(monkeypatch):
    original=q.digest
    q.parents.cache_clear()
    monkeypatch.setattr(q,'digest',lambda p:'0'*64 if Path(p)==q.ROOT/q.PARENT else original(p))
    with pytest.raises(ValueError,match='reviewed parent'):q.parents()
    q.parents.cache_clear()
