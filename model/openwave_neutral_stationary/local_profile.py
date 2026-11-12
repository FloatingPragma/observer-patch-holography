"""Explicit algebraic projector and gauge-invariant native occupation records."""
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import target
from model.coherent_defect_wave.chiral_joint.local_reference import physical
from model.coherent_defect_wave.chiral_joint.local_reference.algebra import LEPTONS
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate,divide,outward
from .isolation import data


@lru_cache(maxsize=1)
def certificate():
    block=target(2,0);h=block['H'];gram=block['gram'];basis=block['basis'];keys=block['keys']
    assert h.rows==16 and gram*h==h.H*gram
    x=s.Symbol('x');e=s.Symbol('e');z=s.Symbol('z')
    q=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    factor=s.Poly(data()['neutral_factor'].as_expr().subs(x,e),e)
    assert s.Poly(q.as_expr().subs(x,e),e).rem(factor).is_zero
    coefficients=[s.Poly(sum(q.nth(j)*e**(j-i-1) for j in range(i+1,q.degree()+1)),e).rem(factor)
                  for i in range(q.degree())]
    denominator=s.Poly(q.diff().as_expr().subs(x,e),e).rem(factor)
    assert s.gcd(denominator,factor).degree()==0
    numerator=sum(c.as_expr()*z**i for i,c in enumerate(coefficients))
    qp=s.Poly(q.as_expr().subs(x,z),z)
    for expression in (numerator**2-denominator.as_expr()*numerator,(z-e)*numerator):
        remainder=s.Poly(expression,z).rem(qp)
        assert all(s.Poly(c,e).rem(factor).is_zero for c in remainder.all_coeffs())
    powers=[s.eye(h.rows)]
    for _ in range(q.degree()-1):powers.append(h*powers[-1])
    rank=s.Poly(sum(c.as_expr()*s.trace(p) for c,p in zip(coefficients,powers)),e).rem(factor)
    assert rank==denominator
    spin_one=s.Poly(target(2,2)['H'].charpoly(e).as_expr(),e)
    assert s.gcd(factor,spin_one).degree()==0
    metric=physical.metric_on(keys);inverse=gram.inv()
    ec_bits=[i for i,n in enumerate(LEPTONS) if n%15==9]
    rows=[];sum_numerator=s.Poly(0,e);total_effect=s.zeros(h.rows)
    for occupation in range(3):
        for scalar_excited in (False,True):
            indicator=s.diag(*(int(2-sum(mask>>i&1 for i in ec_bits)==occupation and (scalar!=0)==scalar_excited)
                               for _,mask,scalar in keys))
            effect=inverse*basis.H*metric*indicator*basis
            assert gram*effect==effect.H*gram
            total_effect+=effect
            weight=s.Poly(sum(c.as_expr()*s.trace(p*effect) for c,p in zip(coefficients,powers)),e).rem(factor)
            sum_numerator+=weight
            interval=(s.Rational(0),s.Rational(0)) if weight.is_zero else divide(
                evaluate(weight,data()['en']),evaluate(denominator,data()['en']))
            assert 0<=interval[0]<=interval[1]<=1
            rows.append({'native_N_ec':occupation,'native_N_L':occupation,
                'scalar_excited':scalar_excited,'probability_numerator':str(weight.as_expr()),
                'probability_denominator':str(denominator.as_expr()),
                'probability_interval':outward(interval,160),
                'probability_decimal_interval':[float(v) for v in interval]})
    assert total_effect==s.eye(h.rows) and (sum_numerator-denominator).rem(factor).is_zero
    assert all(sum(mask>>i&1 for i in range(6))==2 for _,mask,_ in keys)
    return {'native_local_physical_block_dimension':16,'complete_source_weight_coordinate_count':len(keys),
        'full_local_CAR_scalar_dimension':5*2**30,
        'exact_local_H_matrix':[[str(v) for v in row] for row in h.tolist()],
        'exact_local_Gram_matrix':[[str(v) for v in row] for row in gram.tolist()],
        'complete_block_annihilator':str(q.as_expr()),
        'selected_root_factor':str(factor.as_expr()),'selected_root':'smallest real root',
        'projector':'sum_i c_i(e) H^i / q_prime(e)',
        'projector_numerator_coefficients':[str(c.as_expr()) for c in coefficients],
        'projector_denominator':str(denominator.as_expr()),
        'projector_rank':1,'exact_idempotence_and_eigenvalue_identities':True,
        'physical_embedding_uses_actual_native_CAR_scalar_metric':True,
        'all_twelve_local_Gauss_currents':'0 on this physical neutral eigenstate',
        'raw_L_B3':[0,0],'local_fermion_parity':'even','local_spin':0,
        'neutral_state_is_bare_empty_Fock_vacuum':False,
        'complete_gauge_invariant_native_occupation_records':rows,
        'six_record_probability_sum':'1 exactly in the selected algebraic field',
        'other_local_Fock_scalar_states_retained_as_full_instrument_complement':True}
