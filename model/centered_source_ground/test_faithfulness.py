from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from . import faithfulness as f


def test_all_signed_vertex_flips_preserve_actual_cube():
    value=f.signed_vertex_flips()
    assert value['sublattice_sizes']==[4,4]
    assert value['hopping_connected'] and value['hopping_bipartite']
    assert len(value['signed_symmetries'])==8
    # An unsigned vertex flip is not silently substituted for its gauge correction.
    from model.interacting_vacuum import ground
    A=ground.band_data()[0];permutation=np.arange(8)^1
    assert not np.array_equal(A[np.ix_(permutation,permutation)],A)


def test_all_zero_doublon_configurations_leave_the_subspace_under_hopping():
    value=f.no_doublon_kinetic_block()
    assert value['singly_occupied_configurations']==256
    assert value['nonzero_directed_spin_hops_checked']==3072
    assert value['P0_T_P0_zero']


def test_faithfulness_floor_is_strict_and_does_not_require_uniqueness():
    value=f.certificate()
    assert value['double_occupation_strict_lower']=='1/15000'
    assert value['full_four_state_local_CAR_faithful']
    assert value['ground_uniqueness_needed_for_this_proof'] is False
    assert value['infinite_ground_local_faithfulness_inferred'] is False
    assert s.sqrt(17)/2-2>s.Rational(3,50)
    assert s.sqrt(3)/16<s.Rational(1,8)


def test_symbolic_full_site_modular_formula_and_entropy_monotonicity():
    d,a=s.symbols('d a',positive=True)
    N=s.diag(0,1,1,2)
    K=-s.log(a)*s.eye(4)+s.log(a/d)*(N-s.eye(4))**2
    exact=s.diag(-s.log(d),-s.log(a),-s.log(a),-s.log(d))
    assert (K-exact).applyfunc(s.expand_log)==s.zeros(4)
    S=-2*d*s.log(d)-2*(s.Rational(1,2)-d)*s.log(s.Rational(1,2)-d)
    assert s.simplify(s.diff(S,d)-2*(s.log(s.Rational(1,2)-d)-s.log(d)))==0
    assert s.simplify(S.subs(d,s.Rational(1,8))-(3*s.log(2)-s.Rational(3,4)*s.log(3)))==0


def test_full_site_relative_entropy_bound_for_noncommuting_perturbation():
    d=0.03;rho=np.diag([d,0.5-d,0.5-d,d]);epsilon=F(1,7)
    v=np.array([1,2j,3,4j]);v=v/np.linalg.norm(v);sigma=np.outer(v,v.conj())
    delta=sigma-rho;changed=rho+float(epsilon)*delta
    assert np.linalg.norm(rho@sigma-sigma@rho)>0
    eig=np.linalg.eigvalsh(changed)
    relative=float(np.sum(eig*np.log(eig))-np.trace(changed@np.diag(np.log(np.diag(rho)))).real)
    upper=float(f.relative_entropy_upper(epsilon,F(2)))
    assert 0<relative<upper


def test_invalid_entropy_paths_rejected():
    for epsilon in (-1,1,True,0.1):
        with pytest.raises(ValueError):f.relative_entropy_upper(epsilon,F(1))
    with pytest.raises(ValueError):f.relative_entropy_upper(F(1,2),-1)
