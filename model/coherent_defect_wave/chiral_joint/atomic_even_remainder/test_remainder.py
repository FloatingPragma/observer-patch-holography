"""Actual bipartite source, full-resolvent identity and hostile parity controls."""
import numpy as np
import pytest
import scipy.linalg
import sympy as s

from .bounds import quartic_remainder,uniform_quartic_remainder,choose_quartic_coupling,coefficient_certificate
from .symmetry import certificate,cross_sublattice_terms,validate_matrix_symmetry
from ..periodic_family.model import Configuration
from ..atomic_reference_control.local_bound import positive_pivots


def toy():
    h0=s.diag(0,0,2,3,4);p=s.diag(1,1,0,0,0);u=s.diag(1,1,-1,1,-1)
    t=s.zeros(5)
    for i,j,value in ((0,2,1),(0,4,s.Rational(1,2)),(1,2,s.Rational(-1,3)),
                      (1,4,s.Rational(2,3)),(2,3,s.Rational(1,4)),(3,4,s.Rational(-1,5))):
        t[i,j]=t[j,i]=value
    assert sum(abs(value) for value in t)<=6
    return h0,t,p,u


def test_actual_parity_and_larger_source_counterexample():
    result=certificate()
    assert result['actual_oriented_hopping_terms']==24
    assert result['reference_band_parity_phase']==1
    assert result['reference_band_dimension']==256
    assert result['actual_side_four_counterexample']['same_sublattice']
    assert not result['larger_spectral_lattice_theorem_transferred']
    with pytest.raises(ValueError,match='same-sublattice'):
        cross_sublattice_terms(Configuration(side=4))


def test_exact_symmetric_full_resolvent_identity_and_bound():
    h0,t,p,u=toy();assert validate_matrix_symmetry(h0,t,p,u)==1
    d0=h0[2:,2:];delta=s.Rational(1,100)*t[2:,2:]
    r0=d0.inv();plus=(d0+delta).inv();minus=(d0-delta).inv()
    averaged=(plus+minus)/2-r0
    # An exact rational equivalent of the positive-square-root identity.
    expected=r0*delta*r0*delta*(s.eye(3)-r0*delta*r0*delta).inv()*r0
    assert averaged==expected
    g=s.Integer(2);v=s.Rational(6,100)
    bound=v*v/(g*(g*g-v*v))
    positive_pivots(bound*s.eye(3)-averaged)
    positive_pivots(bound*s.eye(3)+averaged)
    b=s.Rational(1,100)*t[2:,:2]
    assert b.T*plus*b==b.T*minus*b
    assert b.T*r0*delta*r0*b==s.zeros(2)


def test_canonical_graph_evenness_and_quartic_diagnostic():
    h0,t,p,u=toy();h0=np.array(h0,dtype=float);t=np.array(t,dtype=float)
    coupling=s.Rational(1,100);value=float(coupling)
    def canonical(lam):
        h=h0+lam*t;values,vectors=scipy.linalg.eigh(h)
        low=vectors[:,:2];x=low[2:,:]@np.linalg.inv(low[:2,:])
        metric=np.eye(2)+x.T@x;vals,vecs=scipy.linalg.eigh(metric)
        j=np.vstack((np.eye(2),x))@(vecs@np.diag(vals**-.5)@vecs.T)
        assert np.linalg.norm(j.T@j-np.eye(2),2)<2e-12
        assert np.linalg.norm(h@j-j@(j.T@h@j),2)<2e-12
        return j.T@h@j
    plus,minus=canonical(value),canonical(-value)
    assert np.linalg.norm(plus-minus,2)<2e-12
    k2=-t[:2,2:]@np.linalg.inv(h0[2:,2:])@t[2:,:2]
    bound=float(uniform_quartic_remainder(2,6,coupling))
    assert np.linalg.norm(plus-value*value*k2,2)<=bound+2e-12


def test_broken_source_parity_and_nonscalar_reference_rejected():
    h0,t,p,u=toy();wrong=s.MutableDenseMatrix(t);wrong[2,2]=1
    with pytest.raises(ValueError,match='full-source parity'):
        validate_matrix_symmetry(h0,wrong,p,u)
    h0=s.diag(0,0,2,3);p=s.diag(1,1,0,0);u=s.diag(1,-1,-1,1);t=s.zeros(4)
    for i,j in ((0,2),(2,3),(3,1)):t[i,j]=t[j,i]=1
    with pytest.raises(ValueError,match='scalar on the reference'):
        validate_matrix_symmetry(h0,t,p,u)


def test_exact_uniform_coefficient_and_rational_coupling_choice():
    assert coefficient_certificate()['endpoint_coefficient']=='1175/48'
    gap=s.Rational(1,25);tau=360;delta=s.Rational(1,10000)
    coupling=choose_quartic_coupling(gap,tau,delta)
    result=quartic_remainder(gap,tau,coupling)
    assert result['sharp']<=result['uniform']<=coupling**2*delta/4
    assert coupling.p==1 and int(coupling.q)&(int(coupling.q)-1)==0
    with pytest.raises(ValueError,match='graph interval'):
        uniform_quartic_remainder(gap,tau,1)
    with pytest.raises(ValueError,match='exact rational'):
        uniform_quartic_remainder(0.04,tau,coupling)
    with pytest.raises(ValueError,match='second-order gap'):
        choose_quartic_coupling(gap,tau,0)
    assert uniform_quartic_remainder(gap,tau,-coupling)==result['uniform']
    assert uniform_quartic_remainder(gap,tau,0)==0


def test_missing_declared_and_transitive_inputs_fail_closed(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action,'SOURCE_FILES',(*custody.action.SOURCE_FILES,'missing-even-remainder.json'))
        with pytest.raises(FileNotFoundError,match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody,'import_closure',lambda seeds:set(seeds)|{'missing-even-remainder-import.py'})
    with pytest.raises(FileNotFoundError,match='missing even remainder import'):
        custody.input_hashes()
