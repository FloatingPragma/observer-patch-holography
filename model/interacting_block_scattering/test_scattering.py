"""Independent vertex, weighted-bound, interval and complete-record controls."""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import numpy as np
from scipy import sparse

from .model import (inputs,member,menus,exact_two_particle_matrix,
                         independent_bog_vertex,moment_bound,preparation)
from .finite_collision import certificate


def result():
    from .model import calculate
    return calculate()


def test_exact_capacity_vertex_and_free_zero():
    for N in (2,8):
        row=exact_two_particle_matrix(N)
        assert row['exact_capacity_effect_matrix_element']>0
        assert row['free_matrix_element_zero']


def test_independent_quasiparticle_quartic_contraction():
    row=independent_bog_vertex()
    assert row['quartic_Bogoliubov_matrix_element']>0
    assert abs(row['quartic_Bogoliubov_matrix_element']-row['closed_form_matrix_element'])<1e-13


def test_packets_have_six_conserved_elastic_channels():
    for side in (32,64):assert len(menus(side))==6


def test_square_root_remainder_includes_overflow():
    for u in np.concatenate((np.linspace(0,1,201),np.linspace(1,20,201))):
        f=np.sqrt(max(1-u,0))
        assert abs(f-1)<=u+1e-14
        assert abs(f-1+u/2)<=u*u+1e-14


def test_global_weighted_operators_including_overflow():
    # Full six canonical modes up to total occupation six. Inputs stop at
    # four, so all intermediate and output columns are represented exactly.
    basis=[n for n in itertools.product(range(7),repeat=6) if sum(n)<=6]
    ids={n:i for i,n in enumerate(basis)}
    rows=[];columns=[];values=[]
    lowers=[]
    for mode in range(6):
        row=[];col=[];data=[]
        for j,n in enumerate(basis):
            if n[mode]:
                target=list(n);target[mode]-=1
                row.append(ids[tuple(target)]);col.append(j);data.append(np.sqrt(n[mode]))
        lowers.append(sparse.csr_matrix((data,(row,col)),shape=(len(basis),len(basis))))
    total=np.array([sum(n) for n in basis])
    local=[np.array([sum(n[:3]) for n in basis]),np.array([sum(n[3:]) for n in basis])]
    accepted=np.flatnonzero(total<=4)
    t=.125;C=3.
    Hfree=sparse.csr_matrix((len(basis),len(basis)))
    V1=Hfree.copy()
    for color in range(3):
        bx=lowers[color];by=lowers[3+color]
        X=bx+bx.T;Y=by+by.T
        nx=sparse.diags(local[0],dtype=float);ny=sparse.diags(local[1],dtype=float)
        cubicx=bx.T@nx+nx@bx;cubicy=by.T@ny+ny@by
        Hfree-=t*(X@Y)
        V1+=t/2*(cubicx@Y+X@cubicy)
    for N in (1,3):
        HN=sparse.csr_matrix(Hfree.shape)
        for color in range(3):
            ax=sparse.diags(np.sqrt(np.maximum(0,1-local[0]/N)),dtype=float)@lowers[color]
            ay=sparse.diags(np.sqrt(np.maximum(0,1-local[1]/N)),dtype=float)@lowers[3+color]
            HN-=t*((ax+ax.T)@(ay+ay.T))
        difference=(HN-Hfree)[:,accepted].toarray()
        remainder=(HN-Hfree-V1/N)[:,accepted].toarray()
        weight=(total[accepted]+2.)**2
        assert np.linalg.norm(difference/weight,2)<=C/N+1e-13
        assert np.linalg.norm(remainder/(total[accepted]+3.)**3,2)<=C/N**2+1e-13
        weighted=difference*(total[:,None]+4.)**2/(total[accepted]+6.)**4
        assert np.linalg.norm(weighted,2)<=C/N+1e-13


def test_every_error_is_below_collision_amplitude_before_squaring():
    for row in result()['finite_members']:
        b=row['exact_bounds'];signal=F(b['leading_collision_amplitude_lower'])
        assert F(b['second_order_Born_amplitude_error_upper'])<=signal/8
        assert F(b['capacity_projection_and_detector_amplitude_error_upper'])<=signal/8
        p=row['complete_preparation'];N=row['constituent_blocks_per_macrocell_N'];V=row['macro_sites']
        r0=F(p['initial_odds_upper_exact']);c=F(p['odds_contraction_exact']);depth=p['binary_filter_depth']
        assert 2*N*V*r0*c**depth<=(signal/8)**2
        assert F(b['conditional_hit_probability_lower'])==(signal/2)**2>0
        assert F(b['native_prepared_packet_probability_outside_actual_energy_window_upper'])<F(1,1000)
        residual=F(b['exact_admitted_packet_energy_residual_over_Delta_upper'])
        width=F(b['energy_halfwidth_over_Delta'])
        assert F(b['native_prepared_packet_probability_outside_actual_energy_window_upper'])==min(F(1),(residual/width+signal/8)**2)


def test_interval_full_capacity_collision_and_complete_preparation():
    row=certificate()
    assert row['exact_reducing_source_dimension']==1024
    assert F(row['exact_source_amplitude_lower'])>F(26,1000)
    assert F(row['conditional_hit_probability_lower_after_preparation_margin'])>F(17,100000)
    parent,_,_=inputs()
    prep=preparation(parent,1,5,F(row['exact_source_amplitude_lower']))
    assert prep['binary_filter_depth']==0
    assert prep['deterministic_initial_trial']
    assert prep['base_trial_success']=='1'
    assert prep['initial_state']=='unfiltered atomic/cube trial'
    assert prep['filter_depth_counts_all_rounds']
    assert prep['all_preparation_flags_retained']


def test_collision_scope_is_finite_and_nonstationary():
    scope=result()['scope']
    assert not scope['prepared_two_packet_state_is_actual_interacting_ground']
    assert not scope['asymptotic_scattering_matrix_proved']
    assert not scope['nontrivial_interacting_relativistic_UV_limit_proved']
    assert scope['measured_inputs']==[]


def test_complete_time_series_and_native_occupation_intervals():
    rows=result()['interval_certified_five_block_collision']['complete_time_series']
    assert [r['clock_Delta_times_t'] for r in rows]==['0','1/4','1/2','3/4','1']
    for row in rows:
        hit=list(map(F,row['native_preparation_hit_probability_interval_exact']))
        complement=list(map(F,row['native_preparation_complement_probability_interval_exact']))
        assert 0<=hit[0]<=hit[1]<=1
        assert complement==[1-hit[1],1-hit[0]]
        assert row['preparation_success_probability_exact']=='1'
        assert row['preparation_failure_probability_exact']=='0'
        assert row['unconditional_hit_probability_interval_exact']==row['native_preparation_hit_probability_interval_exact']
        assert row['unconditional_complement_probability_interval_exact']==row['native_preparation_complement_probability_interval_exact']
        if row['clock_Delta_times_t']!='0':assert hit[0]>0
        for site in row['occupation_marginals']:
            assert len(site['components'])==3
            for color in site['components']:
                lo,hi=map(F,color['ideal_interval_exact'])
                plo,phi=map(F,color['native_preparation_interval_exact'])
                assert 0<=plo<=lo<=hi<=phi<=1
            twirl=site['reference_traced_source_components']
            assert len(twirl)==3
            mean_lower=sum(F(c['ideal_interval_exact'][0]) for c in site['components'])/3
            mean_upper=sum(F(c['ideal_interval_exact'][1]) for c in site['components'])/3
            for color in twirl:
                assert list(map(F,color['ideal_interval_exact']))==[mean_lower,mean_upper]
                assert list(map(F,color['native_preparation_interval_exact']))==[
                    sum(F(c['native_preparation_interval_exact'][0]) for c in site['components'])/3,
                    sum(F(c['native_preparation_interval_exact'][1]) for c in site['components'])/3]
            assert twirl[0]['native_preparation_interval_exact']==twirl[1]['native_preparation_interval_exact']==twirl[2]['native_preparation_interval_exact']
        assert not row['classical_trajectories_claim']


def test_collision_probability_cannot_be_mutated_in_receipt():
    from copy import deepcopy
    from .certificate import replay_differences
    actual=result()
    broken=deepcopy(actual)
    broken['finite_members'][0]['exact_bounds']['conditional_hit_probability_lower']='0'
    assert replay_differences(broken,actual)==['/finite_members/0/exact_bounds/conditional_hit_probability_lower']


def test_source_custody_rejects_bad_paths_and_hashes():
    import pytest
    from .certificate import validate_pins
    for pins in ({'../elsewhere':'0'},{'/tmp/elsewhere':'0'},
                 {'model/interacting_block_scattering/model.py':'0'}):
        with pytest.raises(ValueError):validate_pins(pins)


def test_first_filter_probability_is_counted_once_when_run():
    import math
    parent,_,_=inputs()
    prep=preparation(parent,1,1,F(1,10000))
    assert prep['binary_filter_depth']==1
    r0=F(prep['initial_odds_upper_exact'])
    p0=1/(1+r0)
    assert abs(math.exp(float(prep['log_success_lower_decimal']))-float(p0*p0))<1e-14


def test_matched_free_baseline_is_nonzero_and_separated_at_final_time():
    direct=result()['interval_certified_five_block_collision']
    assert not direct['same_admitted_free_channel_exactly_zero']
    final=direct['complete_time_series'][-1]
    assert F(final['same_admitted_free_hit_probability_interval_exact'][0])>0
    assert F(final['native_minus_matched_free_hit_difference_interval_exact'][0])>F(1,2000)


def test_full_free_kernel_against_independent_fock_matrix_evolution():
    # This is a numerical cross-check, not the interval certificate's
    # definition. The certified kernel has no intervening Fock cutoff.
    from scipy.sparse.linalg import expm_multiply
    V=5;cutoff=12;D=2.;coupling=.125
    basis=[n for n in itertools.product(range(cutoff+1),repeat=V) if sum(n)<=cutoff]
    ids={n:i for i,n in enumerate(basis)}
    rows=[];cols=[];data=[]
    for origin,n in enumerate(basis):
        rows.append(origin);cols.append(origin);data.append(D*sum(n))
        for x in range(V):
            y=(x+1)%V
            for dx,dy in itertools.product((-1,1),repeat=2):
                if (dx<0 and n[x]==0) or (dy<0 and n[y]==0):continue
                target=list(n);target[x]+=dx;target[y]+=dy
                if sum(target)>cutoff:continue
                amplitude=np.sqrt((n[x]+(dx>0))*(n[y]+(dy>0)))
                rows.append(ids[tuple(target)]);cols.append(origin);data.append(-coupling*amplitude)
    H=sparse.csr_matrix((data,(rows,cols)),shape=(len(basis),len(basis)))
    inputs=[(),(0,)];initial=np.zeros((len(basis),len(inputs)),complex)
    def index(sites):
        n=[0]*V
        for i in sites:n[i]+=1
        return ids[tuple(n)]
    for j,sites in enumerate(inputs):initial[index(sites),j]=1
    evolved=expm_multiply(-1j*H,initial)
    controls=certificate()['complete_time_series'][-1]['independent_one_color_kernel_controls']
    for row in controls:
        real=sum(map(F,row['real_interval_exact']))/2
        imag=sum(map(F,row['imaginary_interval_exact']))/2
        exact=complex(float(real),float(imag))
        numeric=evolved[index(row['outgoing_sites']),inputs.index(tuple(row['incoming_sites']))]
        assert abs(numeric-exact)<2e-9
