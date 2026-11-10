"""Exact controller identities, finite filter law and invalid-premise controls."""
import pytest
import sympy as s
from . import model


def test_actual_ground_root_and_complete_source_join():
    result=model.local_recipe()
    assert result['equitable_dimension']==403 and result['full_H_cyclic_degree']==49
    assert result['full_local_ground_and_403_class_vector_joined_by_same_simple_root']
    assert not result['native_efficient_gate_synthesis_executed']


def test_householder_on_a_nontrivial_retained_spin_reference():
    # Two spin halves give a concrete d=4 instance of the dimension-independent identity.
    j=s.kronecker_product(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]))
    u=s.Matrix([s.Rational(3,5),0,0,s.Rational(4,5)])
    psi=sum((s.kronecker_product(s.eye(4)[:,k],j*s.eye(4)[:,k]) for k in range(4)),s.zeros(16,1))/2
    target=s.kronecker_product(u,j*u)
    assert (psi.T*psi)[0]==(target.T*target)[0]==1 and (psi.T*target)[0]==s.Rational(1,2)
    w=psi-target
    reflection=s.eye(16)-2*w*w.T
    assert reflection.T*reflection==s.eye(16) and reflection*psi==target


def test_cayley_filter_exact_two_level_probability_with_retained_complement():
    h=s.diag(0,1);identity=s.eye(2)
    c=s.simplify((identity-s.I*h)*(identity+s.I*h).inv())
    assert c.H*c==identity
    w=s.kronecker_product(c,c.H);m=8
    filter_matrix=sum((w**k for k in range(m)),s.zeros(4))/m
    p=s.Rational(3,4)
    trial=s.Matrix([s.sqrt(p),s.sqrt(1-p)])
    image=filter_matrix*s.kronecker_product(trial,trial)
    success=(image.H*image)[0]
    assert success==p**2+(1-p)**2 and s.conjugate(image[0])*image[0]/success==s.Rational(9,10)
    assert identity-c.H*c==s.zeros(2)
    control=model.filter_bounds(p,s.Integer(1),s.Integer(1),s.Rational(1,4))
    assert control['length']==8 and control['success_lower']<=success
    assert control['ground_pair_fidelity_lower']<=s.Rational(9,10)
    assert filter_matrix.H*filter_matrix+(s.eye(4)-filter_matrix.H*filter_matrix)==s.eye(4)


def test_full_preparation_has_finite_M_and_positive_readout_bounds():
    result=model.certificate()
    control=result['finite_energy_agreement_filter']
    assert int(control['explicit_integer_M'])>10**20
    assert s.Rational(control['success_probability_lower'])>s.Rational(9999999,10000000)
    assert s.Rational(control['conditional_ground_pair_fidelity_lower'])>1-s.Rational(1,10**12)
    assert s.Rational(control['conditional_one_source_trace_distance_upper'])<s.Rational(1,10**6)
    for low,high in result['same_source_readout_after_approximate_preparation']['conditional_probability_intervals'].values():
        assert 0<s.Rational(low)<=s.Rational(high)<1
    assert not result['same_source_readout_after_approximate_preparation']['prepared_state_asserted_exactly_stationary']


def test_invalid_probability_gap_and_leakage_are_rejected():
    for args in ((s.Rational(1,2),1,1,s.Rational(1,4)),
                 (1,0,1,s.Rational(1,4)), (1,1,0,s.Rational(1,4)),
                 (1,1,3,s.Rational(1,4)), (1,1,1,0), (1,1,1,1), (0.75,1,1,s.Rational(1,4))):
        with pytest.raises(ValueError,match='rational p>1/2'):
            model.filter_bounds(*args)
