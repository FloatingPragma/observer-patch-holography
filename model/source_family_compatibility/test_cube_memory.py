from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from model.source_family_compatibility import cube_memory as m,composite


def test_exact_inverse_residual_rejects_a_rounded_guess_as_identity():
    a=(((F(1),F(2)),(F(1,3),F())),((F(1,3),F()),(F(-2),F(1))))
    q,bound,res=m.inverse_proposal(a)
    exact=s.Matrix([[s.Rational(x.numerator,x.denominator)+s.I*s.Rational(y.numerator,y.denominator)for x,y in r]for r in a]).inv()
    approx=np.asarray([[complex(float(x),float(y))for x,y in r]for r in q])
    assert res>0 and res<1
    assert np.linalg.norm(approx-np.array(exact,dtype=complex),2)<=float(bound*res/(1-res))
    assert bound>0


def test_complete_four_filter_records_on_literal_full_source():
    d=composite.source((2,1,1),F(3,2),F(4));num,den,_=composite.ground(d,F(-9,2))
    h=s.Matrix(d['H'])/d['scale'];rho=s.Matrix(num)/den;eye=s.eye(16)
    delta=h+s.Rational(9,2)*eye;eta=s.Rational(1,4)
    filt=eta*(eta*eye-s.I*delta).inv();success=filt.H*filt
    a=s.Matrix(d['a'][0]);effects=[];prob=[]
    for op in(a.H,a):
        effects.extend((op.H*success*op,op.H*(eye-success)*op))
        prob.extend((s.simplify(s.trace(rho*effects[-2])),s.simplify(s.trace(rho*effects[-1]))))
    assert sum(effects,s.zeros(16))==eye
    assert sum(prob)==1 and all(x>=0 for x in prob)
    g=s.trace(rho*(a*(s.I*eta*eye-delta).inv()*a.H+a.H*(s.I*eta*eye+delta).inv()*a))
    assert s.simplify(prob[0]+prob[2]+eta*s.im(g))==0
    assert prob[0]==prob[2] and prob[1]==prob[3]


def test_actual_default_cube_memory_is_not_dropped():
    c=m.certificate()
    assert c['full_Fock_dimension']==65536
    assert c['complete_ground_block_dimension']==4900
    assert c['charged_spin_block_dimension']==3920
    assert c['frequencies'][0]['original_minus_memory_free_probability'][0]==F(23,1000)
    assert c['imaginary_memory_change_strict_lower']==F(15,1000)
    assert c['local_mode_charge_filter_records']==4
    assert not c['arbitrary_reference_channel_for_retarded_columns_claimed']
    assert not c['thermodynamic_charge_gap_proved']
    assert c['physical_mass_assigned']is None


def test_actual_charged_columns_equal_full_Fock_actions():
    d=m._source();base=m.reference.certified_ground();p=d['ground']
    b=d['b'];assert b.shape==(3920,16)
    assert d['symmetry']['complete_source_clock_invariance']
    assert (d['clock_numerator']-d['clock_numerator'].T).nnz==0
    assert p['half_trace_distance_upper']<=m.GROUND_HALFTRACE
    assert base['full_gap_premise']['full_ground_multiplicity']==1
    full=np.zeros(65536,dtype=np.int64);full[d['balanced_indices']]=d['trial']
    operators=m.reference._operators()
    for x in range(8):
        assert np.array_equal((operators['a'][2*x].T@full)[d['charged_indices']],b[:,x])
        assert np.array_equal((operators['b'][2*x].T@full)[d['charged_indices']],b[:,x+8])
        assert np.array_equal(b[:,x]**2,b[:,x+8]**2)


def test_frequency_output_is_detached_from_cached_rational_slots():
    a=m.frequency(F(1,4));target=F(a['difference_interval'][0].numerator,a['difference_interval'][0].denominator)
    a['difference_interval'][0]._numerator=0
    b=m.frequency(F(1,4));assert b['difference_interval'][0]==target


def test_noncanonical_clock_and_solver_inputs_rejected():
    for eta in(0,-1,.25,True):
        with pytest.raises(ValueError):m.frequency(eta)


def test_corrupted_inverse_solver_fails_exact_residual(monkeypatch):
    a=(((F(1),F(1)),),)
    monkeypatch.setattr(np.linalg,'inv',lambda _:np.zeros((1,1),dtype=complex))
    _,_,res=m.inverse_proposal(a)
    assert res==1


def test_bad_source_solve_rejected_by_integer_residual(monkeypatch):
    class Wrong:
        def solve(self,b):return np.zeros(b.shape,dtype=complex)
    monkeypatch.setattr(m,'splu',lambda _:Wrong())
    with pytest.raises(ArithmeticError,match='rational residual'):
        m._rounded_source_solve(m._source(),F(1,4))
