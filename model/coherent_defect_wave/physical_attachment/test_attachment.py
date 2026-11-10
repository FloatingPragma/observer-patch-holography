"""Physical attachment identities and controls against false inferences."""
from functools import lru_cache
import pytest
import sympy as s

from . import spin,chiral,compensator


@lru_cache(maxsize=1)
def spin_receipt():
    return spin.certificate()


@lru_cache(maxsize=1)
def chiral_receipt():
    return chiral.certificate()


def test_actual_routing_spin_currents_and_finite_symmetry_boundary():
    r=spin_receipt()
    assert r['spin_squared']=='3/4 I4'
    assert r['exact_Lorentz_vector_current_identities']==12
    assert r['finite_ordered_route_quarter_turn_residual_HS_squared']=='12'
    for row in r['finite_complete_position_and_flux_records']:
        before,after,flux=(list(map(s.Rational,row[name])) for name in
                           ('position_before','position_after','directed_edge_flux'))
        assert sum(before)==sum(after)==1
        assert all(after[x]-before[x]==flux[(x-1)%4]-flux[x] for x in range(4))


def test_density_cannot_replace_scalar_mass_bilinear():
    alpha,beta,_,_=spin.spin_data()
    boost=s.Rational(5,3)*s.eye(4)+s.Rational(4,3)*alpha[1]
    assert spin.clean(boost.H*beta*boost-beta)==s.zeros(4)
    assert spin.clean(boost.H*boost-s.eye(4))!=s.zeros(4)


def test_rotation_rejects_nonunit_or_inexact_angle():
    for axis,c,d in ((3,1,0),(0,1,1),(0,.6,.8),(0,s.I,0)):
        with pytest.raises(ValueError):
            spin.rotation(axis,c,d)
    with pytest.raises(ValueError):
        spin.one_axis_shift({(0,0):s.Symbol('x')},0,4)


def test_all_chiral_current_anomalies_and_mass_equations():
    r=chiral_receipt()
    assert r['anomalies']['symmetric_cubic_current_trace_checks']==364
    assert r['mass_and_Yukawa']['unbroken_hypercharge_invariant_Weyl_mass_dimension']==0
    assert r['mass_and_Yukawa']['full_current_Yukawa_identities']==72
    assert r['mass_and_Yukawa']['two_gauge_compatible_mass_square_trace_controls']==['14','48']


def test_removing_a_charged_singlet_breaks_anomaly_cancellation():
    _,_,generators,_=chiral.representation()
    y=generators[-1]
    kept=[i for i,q in enumerate(y.diagonal()) if q!=1]
    reduced=y.extract(kept,kept)
    assert s.trace(y)==s.trace(y**3)==0
    assert s.trace(reduced)==s.trace(reduced**3)==-1


def test_a_constant_charged_mass_violates_actual_hypercharge():
    _,_,generators,_=chiral.representation()
    y=generators[-1]
    trial=s.zeros(15);trial[0,7]=trial[7,0]=1
    assert y.T*trial+trial*y!=s.zeros(15)
    assert len({str(x) for x in y.diagonal()})==5


def test_determinant_charge_and_complete_memory_compensation():
    r=chiral_receipt()['mass_and_Yukawa']
    assert r['full_U3_times_U2_phase_characters']==[1]*5
    assert r['determinant_character_component_checks']==30
    for pair,h in (((2,2),1),((2,4),-1),((4,2),-1)):
        record=compensator.component_control(pair,h)
        p=list(map(s.Rational,record['complete_sixteen_detector_probabilities']))
        assert len(p)==16 and sum(p)==1 and p[0]==p[-1]==s.Rational(1,2)
    with pytest.raises(ValueError):
        compensator.component_control((2,2),-1)
