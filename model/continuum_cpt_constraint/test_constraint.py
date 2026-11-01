from fractions import Fraction as F

import pytest
import sympy as s

from . import algebra as a, bounds as b


def test_exact_packet_bounds_and_broad_limit():
    rows=[b.gaussian(F(1,8),1,width) for width in (1,2,4,8)]
    assert [row['TV_from_reflection_lower'] for row in rows]==[
        F(32,67),F(128,259),F(512,1027),F(2048,4099)]
    assert all(row['signed_mean']==F(1,2) for row in rows)
    assert all(row['W1_from_reflection_lower']==1 for row in rows)
    assert all(row['clipping_radius']==row['second_moment_upper'] for row in rows)
    assert b.gaussian(1,0,1)['clipping_radius'] is None
    assert b.gaussian(1,0,1)['TV_from_reflection_lower']==0


def test_TV_convention_and_nontrivial_equality_case():
    # The laws share a zero atom: TV is 1/4, not the full norm 1/2.
    law=b.measure(((0,F(3,4)),(4,F(1,4))))
    reflected=b.reflection(law)
    _,mean,second=b.moments(law)
    assert (mean,second)==(1,4)
    assert b.total_variation(law,reflected)==mean*mean/second==F(1,4)
    assert b.wasserstein_one(law,reflected)==2*mean==2
    radius=second/(2*mean)
    assert b.clipped_mean(law,radius)==mean/2


def test_asymmetric_measures_satisfy_bounded_witness_without_poles_assumption():
    for weight in (F(1,7),F(1,3),F(3,5)):
        law=b.measure(((-1,(1-weight)/2),(0,(1-weight)/2),(7,weight)))
        _,mean,second=b.moments(law)
        assert mean>0
        radius=second/(2*mean)
        assert b.clipped_mean(law,radius)>=mean/2
        assert b.total_variation(law,b.reflection(law))>=mean*mean/second
        assert b.wasserstein_one(law,b.reflection(law))>=2*mean


def test_clipping_tail_constant_is_sharp():
    radius=F(7,3)
    for energy in (F(-20),-2*radius,-radius,F(0),radius,2*radius,F(20)):
        clipped=max(-radius,min(radius,energy))
        assert abs(energy-clipped)<=energy*energy/(4*radius)
    assert abs(2*radius-radius)==(2*radius)**2/(4*radius)


def test_full_Fock_CPT_requires_generator_invariance():
    free=a.full_fock_control(2,0)
    shifted=a.full_fock_control(2,s.Rational(1,2))
    assert free['CPT_generator_defect']==s.zeros(4)
    assert shifted['CPT_generator_defect']!=s.zeros(4)
    assert list(shifted['positive_GNS_generator'].diagonal())==[s.Rational(3,2),4,0,s.Rational(5,2)]
    assert shifted['addition_cost']==s.Rational(5,2)
    assert shifted['removal_cost']==s.Rational(3,2)
    law=b.shifted_pair(2,F(1,2))
    assert b.moments(law)==(1,F(1,2),F(17,4))
    assert b.reflection(law)!=law


def test_zero_first_moment_is_not_ground_positivity():
    control=a.centering_failure()
    assert control['original_generator_diagonal']==[0,1,3,4]
    assert control['centered_generator_diagonal']==[0,-1,1,0]
    centered=b.measure(((-1,F(1,2)),(1,F(1,2))))
    assert b.moments(centered)[1]==0 and b.reflection(centered)==centered
    assert control['same_reference_is_centered_ground'] is False


def test_massive_Dirac_trace_and_eight_free_chiralities():
    control=a.clifford_control()
    assert control['Dirac_first_moment_trace']==0
    assert control['mass_term_is_scalar_energy_shift'] is False
    rows=control['free_valleys']
    assert len(rows)==8
    assert sum(row['chirality']==1 for row in rows)==4
    assert sum(row['chirality']==-1 for row in rows)==4
    massive=b.measure(((-5,F(1,2)),(5,F(1,2))))
    assert b.moments(massive)==(1,0,25)
    assert b.reflection(massive)==massive


@pytest.mark.parametrize('call',[
    lambda:b.gaussian(0,1,1),lambda:b.gaussian(1,-1,1),
    lambda:b.gaussian(1,1,0),lambda:b.gaussian(0.5,1,1),
    lambda:b.measure(((0,F(1,2)),)),lambda:b.measure(((0,-1),(1,2))),
    lambda:b.clipped_mean({F(1):F(1)},0),lambda:a.full_fock_control(1,2),
])
def test_invalid_domains_rejected(call):
    with pytest.raises(ValueError):call()
