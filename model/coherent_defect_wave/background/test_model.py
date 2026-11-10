"""Exact background, reversible preparation and complete outcome checks."""
from fractions import Fraction as F
import pytest
import sympy as s
from .model import certificate, compiled_paths, compiled_link_event, background
from ..model import JointObserver, SparseState, combine
from ..gauss import gamma_amplitudes


def test_same_source_background_and_complete_preparation():
    result=certificate()
    assert len(result['fixed_triplet_backgrounds'])==3
    assert sum(F(r['conditional_probability']) for r in result['ten_level_fourier_ancilla_outcomes'])==1
    assert result['joint_success_per_fresh_ambient_trial']=='1/10'
    assert not result['scope']['selected_physical_vacuum']


def test_triplet_domain_and_scratch_reversal():
    observer=JointObserver()
    for internal in (-1,3,4,True,1.0):
        with pytest.raises(ValueError):
            background(observer,internal=internal)
    for edge in observer.law.interior:
        for port in range(12):
            output=compiled_link_event(observer,1,port,0,edge)
            assert compiled_link_event(observer,*output,edge)==(1,port,0)


def test_compiler_rejects_edges_outside_its_interior_domain():
    observer=JointObserver()
    boundary=next(e for e in range(len(observer.law.edges)) if e not in observer.law.interior)
    for edge in (boundary,-1,len(observer.law.edges),True,float(observer.law.interior[0])):
        with pytest.raises(ValueError,match='interior edge'):
            compiled_link_event(observer,1,0,0,edge)


def test_blank_scratch_is_necessary_for_unconditional_link_compilation():
    observer=JointObserver()
    edge=observer.law.edges.index((1,5))
    intended=observer.law.transport(1,edge)[0]
    assert intended != 1
    assert compiled_link_event(observer,1,0,0,edge)==(intended,0,0)
    # A dirty scratch remains reversible but does not implement the promised T_e.
    dirty=compiled_link_event(observer,1,0,1,edge)
    assert dirty==(1,0,1)
    assert compiled_link_event(observer,*dirty,edge)==(1,0,1)


def test_intermediate_port_routing_preserves_gauss_only_in_its_declared_domain():
    observer=JointObserver()

    def route_zero_scratch(state):
        return SparseState({(g,5 if p==0 else 0 if p==5 else p,j,r,a):value
                            for (g,p,j,r,a),value in state.amplitudes},
                           state.norm_denominator)

    def gauge(state):
        return SparseState(gamma_amplitudes(state.as_dict()),state.norm_denominator)

    for internal in (0,3):
        state=observer.basis_state(1,0,internal=internal)
        left=gauge(route_zero_scratch(state))
        right=route_zero_scratch(gauge(state))
        if internal==0:
            assert left==right
        else:
            # Moving a doublet between ports0 and5 changes its gauge sign.
            assert combine(left,right,1,-1).norm_squared()==4
            assert combine(left,right,1,1).norm_squared()==0


def test_faithful_stationarity_does_not_supply_a_positive_energy_vacuum():
    H=s.diag(0,1)
    faithful=s.diag(s.Rational(1,3),s.Rational(2,3))
    lower=s.Matrix([[0,1],[0,0]])
    upper=lower.adjoint()
    assert H*faithful==faithful*H
    assert H*lower-lower*H==-lower
    assert H*upper-upper*H==upper
    assert s.trace(faithful*lower.adjoint()*lower)==s.Rational(2,3)
    assert s.trace(faithful*upper.adjoint()*upper)==s.Rational(1,3)
    # In a pure ground reference the negative-frequency operator has zero norm.
    ground=s.diag(1,0)
    assert s.trace(ground*lower.adjoint()*lower)==0
    assert s.trace(ground*upper.adjoint()*upper)==1
