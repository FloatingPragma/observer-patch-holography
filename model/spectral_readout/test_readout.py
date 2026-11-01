from fractions import Fraction as F

import numpy as np
import pytest

from model.certified_reference import interval as iv
from . import kernels as k,energy,finite


def contains(interval,value,tolerance=0):
    return float(F(interval.lo,iv.SCALE))-tolerance<=value<=float(F(interval.hi,iv.SCALE))+tolerance


def test_complete_positive_four_outcome_measurement_on_exact_energy_grid():
    # Mixture of two distinct grid eigenstates. Cardinal Fourier phases are exact.
    chars=[]
    for n in range(1,4):
        s1,c1=k.phase(4,n);s3,c3=k.phase(4,3*n)
        chars.append((c1/3+2*c3/3,-s1/3-2*s3/3))
    result=k.probabilities(chars)
    for interval,value in zip(result,(0,F(1,3),0,F(2,3)),strict=True):
        assert F(interval.lo,iv.SCALE)<=value<=F(interval.hi,iv.SCALE)


def test_energy_alias_is_real_and_cannot_be_removed_by_more_fourier_bins():
    # E=0 and E=2*pi/tau give exactly the same sampled evolution for every n.
    zero=[(iv.Interval.point(1),iv.Interval.point(0)) for _ in range(15)]
    output=k.probabilities(zero)
    assert contains(output[0],1)
    assert all(contains(row,0) for row in output[1:])
    assert k.real_energy_leakage(64,1,F(1,2),1)['alias']>0
    assert k.real_energy_leakage(64,1,F(1,2),1,support_upper=2)['alias']==0


@pytest.mark.parametrize('angle',[-1000,-17,-F(1,7),0,F(22,7),137])
def test_argument_halving_trigonometric_enclosures(angle):
    sine,cosine=k.trig_interval(angle,angle)
    assert contains(sine,np.sin(float(angle)),1e-14)
    assert contains(cosine,np.cos(float(angle)),1e-14)


def test_uniform_grid_leakage_bound_on_off_grid_energies():
    m=64;tau=F(1);width=F(1,2)
    bound=float(k.circular_leakage(m,tau,width))
    grid=2*np.pi*np.arange(m)/m
    for theta in np.linspace(-np.pi,np.pi,113):
        distances=np.abs((theta-grid+np.pi)%(2*np.pi)-np.pi)
        amplitudes=np.exp(-1j*np.outer(np.arange(m),theta-grid)).mean(axis=0)
        assert abs(np.sum(abs(amplitudes)**2)-1)<1e-13
        assert np.sum(abs(amplitudes[distances>=float(width)])**2)<=bound+1e-13


def test_exact_operator_moments_include_contact_variance():
    # V embeds into the first two coordinates; U=2,r=1/2 makes W=I exact.
    h=np.array([[0,2],[2,0]],dtype=np.int64);identity=np.eye(2,dtype=np.int64)
    L=np.block([[h+identity,identity],[identity,np.zeros((2,2),dtype=np.int64)]])
    assert np.array_equal((L@L)[:2,:2],h@h+2*h+2*identity)
    assert not np.array_equal((L@L)[:2,:2],(h+identity)@(h+identity))
    assert energy.local_moments()=={
        'first_moment':[F(1,8),F(1,2)],
        'second_moment':[F(19,128),F(67,128)],
        'variance':[F(17,128),F(35,128)]}
    assert energy.energy_tail(4,F(67,128))==F(67,2048)
    assert energy.kinetic_controlled_second_moment(1)==F(5,2)


def test_negative_response_frequency_has_positive_removal_energy():
    result=energy.charge_tag_certificate()
    assert result['two_mode_control']['positive_excitation_costs']==['3','2']
    assert result['two_mode_control']['signed_response_frequencies']==['3','-2']
    assert result['negative_ground_excitation_energy_inferred'] is False


def test_unknown_ground_discrepancy_cannot_silently_be_zero():
    inputs=dict(count_=4,tau=F(1,8),selected_count=1,half_side=4096,
                buffer_radius=64,support_radius=0,support_sites=1)
    unknown=energy.comparison_contract(**inputs)
    assert unknown['state_discrepancy'] is None
    assert unknown['total_infinite_ground_bin_error'] is None
    supplied=energy.comparison_contract(**inputs,state_discrepancy=F(1,100))
    assert supplied['total_infinite_ground_bin_error'] is not None
    assert supplied['state_discrepancy_is_supplied_conditional_premise']
    assert supplied['ground_state_convergence_rate_derived'] is False


def test_complete_bin_event_has_zero_error_but_partial_event_does_not():
    assert k.bin_fourier_error(16,16,F(1,100))==0
    assert k.bin_fourier_error(16,1,F(1,100))==F(3,320)
    assert k.bin_fourier_error(16,15,F(1,100))==F(3,320)
    assert energy.sampling_failure(16,100000,F(1,50))<F(1,10**15)
    assert k.clip_event_estimate(F(-1,4))==0
    assert k.clip_event_estimate(F(7,4))==1
    assert k.clip_event_estimate(F(3,4))==F(3,4)


def test_frozen_78_line_distribution_and_independent_complex_amplitudes():
    m=16;report=finite.readout(m);source=finite.load_source()
    energies=np.array([float((lo+hi)/2) for lo,hi,_,_ in source])
    weights=np.array([float((lo+hi)/2) for _,_,lo,hi in source])
    for row in report['rows']:
        amplitudes=np.exp(-1j*np.outer(np.arange(m),energies-2*np.pi*row['index']/m)).mean(axis=0)
        value=float(weights@abs(amplitudes)**2)
        low,high=map(F,row['probability_interval'])
        assert float(low)-1e-14<=value<=float(high)+1e-14
    assert report['all_outcomes_retained']
    assert report['ambiguous_resolution_pairs']==0
    assert report['joint_error_probability_decimal']==['0.055193158231','0.055193158232']
    assert report['native_full_energy_POVM_executed'] is False


def test_native_time_one_anchor_does_not_claim_higher_harmonics():
    value=finite.native_first_harmonic()
    assert value['previous_complete_native_first_harmonic_agrees']
    assert value['all_higher_harmonics_natively_executed'] is False


def test_weak_convergence_does_not_imply_sharp_atom_probability_convergence():
    # delta_(1/n) -> delta_0, but the mass of {0} jumps from zero to one.
    assert all(1/n!=0 for n in range(1,100))
    assert abs(np.exp(1j/100000)-1)<2e-5
    value=energy.certificate()
    assert value['sharp_interval_convergence_requires_zero_endpoint_mass']
    assert value['common_Hilbert_space_strong_resolvent_convergence_of_dilations_claimed'] is False


@pytest.mark.parametrize('action',[
    lambda:k.circular_leakage(1,1,1),lambda:k.circular_leakage(8,1,4),
    lambda:k.real_energy_leakage(8,1,3,1),lambda:k.exact(0.5,'input'),
    lambda:energy.energy_tail(0,1),lambda:energy.kinetic_controlled_second_moment(-1),
    lambda:energy.sampling_failure(8,0,F(1,10)),
    lambda:k.bin_fourier_error(8,9,F(1,10)),
])
def test_invalid_domains_rejected(action):
    with pytest.raises(ValueError):action()
