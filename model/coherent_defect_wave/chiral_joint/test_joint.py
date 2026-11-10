"""Exact covariance, same-H quantum encounter and complete reference bounds."""
import sympy as s
from .model import (finite_program,path_ward_checks,joint_encounter,joint_column,
                    apply_control,norm2,norm_bound,physical_reference_certificate)


def test_actual_quantum_path_ward():
    result=path_ward_checks()
    assert result['raising_tensor_terms']=={'1':59,'2':279}
    assert result['all_endpoint_and_internal_vertex_Ward_checks']==60
    program=finite_program(4)
    assert program['gauge_edges']==144
    assert max(len(t['path']) for t in program['kinetic_oriented_terms'])==3
    assert not program['quantum_transporter_set_to_identity']


def test_literal_same_H_creation_and_transport():
    result=joint_encounter()
    assert result['first_H_vacuum_norm_squared']=='4'
    assert result['second_H_vacuum_arrival_norm_squared']=='3'
    assert s.sympify(result['complete_arrival_yes_probability_strict_lower'])>0
    assert result['full_active_control_dimension']==122880
    # Remove all transport destinations: the original pair/scalar creation
    # remains but no left-cell arrival can result from two local-only steps.
    first=apply_control({(0,0,0):1})
    local_second={}
    for key,coefficient in first.items():
        for target,value in joint_column(key).items():
            if target[1]==key[1]:local_second[target]=local_second.get(target,0)+coefficient*value
    assert all(not (key[0]&63) for key,value in local_second.items() if value)


def test_positive_complete_physical_reference_filter():
    result=physical_reference_certificate()
    assert result['unit_parameter_full_H_norm_upper_bound']=='1884'
    assert result['prepared_success_state_exactly_rho_beta']
    assert not result['rho_beta_faithful_on_full_unconstrained_carrier']
    assert not result['global_eigensystem_numerically_computed']
    # Exact spectral-function control including a nonphysical sector.
    h=s.diag(-2,1,3);p=s.diag(1,1,0);b=s.Integer(3);beta=s.log(2)
    l=p*s.diag(*[s.exp(-beta*(e+b)/2) for e in h.diagonal()])
    effect=l.H*l
    failure=s.diag(*[s.sqrt(1-x) for x in effect.diagonal()])
    assert s.simplify(effect+failure.H*failure)==s.eye(3)
    dilation=l.row_join(-failure).col_join(failure.row_join(l))
    assert s.simplify(dilation.H*dilation)==s.eye(6)
    post=l*l.H/s.trace(l*l.H)
    expected=p*s.diag(*[s.exp(-beta*e) for e in h.diagonal()])
    expected=expected/s.trace(expected)
    assert s.simplify(post-expected)==s.zeros(3)
    assert s.trace(effect)/s.trace(p)>=s.exp(-2*beta*b)
    assert s.trace(effect)/3>=s.trace(p)*s.exp(-2*beta*b)/3
    assert result['input']=='P/r on the physical Gauss sector, r=rank(P)'
    assert result['physical_input_preserves_total_fermion_parity']
    assert not result['physical_input_uses_unphysical_Gauss_states']
    # A physical sector containing both source parities requires the
    # matched graded ancilla. The same purification is total even.
    parity=s.diag(1,-1)
    purification=s.Matrix([1,0,0,1])/s.sqrt(2)
    total_parity=s.kronecker_product(parity,parity)
    assert total_parity*purification==purification
    reduced=s.zeros(2)
    for a in range(2):
        for b in range(2):
            reduced[a,b]=sum(purification[2*a+j]*s.conjugate(purification[2*b+j]) for j in range(2))
    assert reduced==s.eye(2)/2
