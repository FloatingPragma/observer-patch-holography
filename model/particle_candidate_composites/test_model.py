"""Exact source current, complete-instrument and candidate exclusion controls."""
from copy import deepcopy
from fractions import Fraction as F

import pytest
import sympy as s

from . import model as m


def test_all_six_literal_source_polynomials_have_the_claimed_currents():
    rows=m.calculated_rows()
    assert len(rows)==6
    assert [r['physical_gauge_kernel_dimension'] for r in rows]==[6,6,6,6,8,8]
    assert [r['S_increment'] for r in rows]==[1,-1,3,-3,2,-2]
    assert all(r['native_Q6']==0 and r['gauge_generators_annihilating_every_column']==12 for r in rows)
    for row in rows:
        delta=row['native_increment'];local=m.spin.local_multiplet(delta)
        for column in range(local['kernel_basis'].cols):
            for ops,coefficient in m.spin.composite_polynomial(delta,column):
                assert coefficient and m.polynomial_delta(ops)==delta


def test_conditional_baryon_lepton_bookkeeping_and_parity_are_separate():
    rows=m.calculated_rows()
    assert [(r['conditional_B'],r['conditional_L']) for r in rows]==[
        (0,1),(0,-1),(-1,0),(1,0),(-1,1),(1,-1)]
    assert [r['CAR_parity'] for r in rows]==[1,1,1,1,0,0]
    assert all(r['conditional_BL_requires_constituent_field_dictionary'] for r in rows)
    assert not rows[-1]['direct_single_fermion_parity_compatible']


def test_actual_weak_and_color_casimirs_forbid_nonzero_direct_singlet_maps():
    currents=m.spin.sectors.native.native_currents()
    # Weak and color generators act faithfully on their native nontrivial rows.
    for q in (-3,-4):
        indices=next(indices for charge,indices in m.spin.sectors.native_blocks() if charge==q)
        blocks=[A.extract(indices,indices) for A in currents]
        stack=s.Matrix.vstack(*blocks)
        assert stack.nullspace()==[]
        assert sum((A.H*A for A in blocks),s.zeros(len(indices))).det()>0
    assert all(r['direct_unbroken_weak_doublet_current_intertwiner_dimension']==0
               for r in m.calculated_rows())


def test_cube_channels_distinguish_scalar_and_vector_apparatus():
    assert m.tensor_channels('G1','G1')=={'A1':1,'T1':1}
    assert m.tensor_channels('G1','H')=={'E':1,'T1':1,'T2':1}
    assert m.tensor_channels('H','H')=={'A1':1,'A2':1,'E':1,'T1':2,'T2':2}
    for i,row in enumerate(m.calculated_rows()):
        suffix='g' if i<4 and i%2==0 else 'u'
        assert all(label.endswith(suffix) for label in row['uniform_full_cube_irreps'])


def _car(mode,total=2):
    A=s.zeros(2**total)
    for mask in range(2**total):
        if not mask&(1<<mode):
            A[mask|(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return A


def test_ward_rules_on_a_faithful_stationary_polynomial_CAR_control():
    a,b=_car(0),_car(1);N=a*a.H;M=b*b.H;H=N+2*M
    Q=s.eye(4)+(3*s.eye(4)-H)/4;rho=Q*Q/s.trace(Q*Q)
    assert rho.det()>0 and rho*H==H*rho
    assert s.trace(rho*a.H*b)==0 and s.trace(rho*a*b)==0
    assert s.trace(rho*a*a.H)>0
    da=(1,0,0,0,0);db=(0,1,0,0,0)
    assert not m.normal_charge_allowed(da,db)
    assert not m.anomalous_charge_allowed(da,db)
    assert m.anomalous_charge_allowed(da,tuple(-x for x in da))
    assert m.normal_charge_allowed(da,da)
    for left in m.calculated_rows():
        for right in m.calculated_rows():
            assert m.normal_charge_allowed(left['native_increment'],right['native_increment'])==(left==right)
    with pytest.raises(ValueError):m.normal_charge_allowed((1,2),da)


def test_scalar_ward_zero_does_not_erase_reference_operator_or_compensator():
    raising=s.Matrix([[0,0],[1,0]])
    psi=s.Matrix([1,0,0,1])/s.sqrt(2);joint=psi*psi.H
    X=s.kronecker_product(raising,s.eye(2))*joint
    ref=s.Matrix(2,2,lambda i,j:sum(X[2*k+i,2*k+j] for k in range(2)))
    assert s.trace(X)==0 and ref!=s.zeros(2)
    assert s.trace(s.kronecker_product(raising,raising)*joint)==s.Rational(1,2)


def test_complete_meter_cancels_every_initial_charge_with_reference_retained():
    result=m.meter_control();rows=result['rows'];d=7
    assert len(rows)==49
    assert sum((row['K'].H*row['K'] for row in rows),s.zeros(d))==s.eye(d)
    shifts=(-3,-2,-1,1,2,3,0)
    psi=s.zeros(d*d,1)
    for j in range(d):psi[j*d+j]=s.Rational(j+1)/s.sqrt(140)
    assert (psi.H*psi)[0]==1
    mass=0
    for row in rows:
        if row['meter']!=shifts[row['probe']]%7:assert row['K']==s.zeros(d)
        vector=s.kronecker_product(row['K'],s.eye(d))*psi
        mass+=(vector.H*vector)[0]
    assert s.simplify(mass)==1
    assert result['coupling'].H*result['coupling']==s.eye(49)


def test_actual_two_trion_alias_is_nonzero_and_not_a_single_family_decoder():
    row=m.weak_pair_alias()
    assert row['two_G1_column_product_anchor_norm_squared']=='1'
    assert row['native_increment']==(0,4,0,0,2)
    assert row['S_increment']==2 and row['aliases_single_even_meter_record']
    assert row['full_native_increment_differs']
    assert not row['pair_creation_term_preserves_native_U1_five']
    assert not row['physical_Majorana_mass_identified']
    recipe=m.meter_recipe()
    assert recipe['cyclic_radix']==7 and recipe['original_complete_labels']==321
    assert recipe['formal_meter_refined_labels']==2247
    assert recipe['source_controlled_additions_both_couplings']==448
    assert not recipe['multiple_excitation_universal_decoder']


def test_full_CAR_norm_bound_and_complete_original_weights():
    for row in m.calculated_rows():
        for j,t in enumerate(row['normalized_column_monomial_counts']):
            terms=m.spin.composite_polynomial(row['native_increment'],j)
            assert len(terms)==t and sum(s.conjugate(c)*c for _,c in terms)==1
            assert s.simplify(t-sum(abs(c) for _,c in terms)**2)>=0
    member=m.contract.member(0);bound=m.bind_member(member);m.verify_bound(member,bound)
    value=bound['result'];assert value['physical_density_floor_factored']['negative_exponent']==582407133
    for row in value['candidate_rows']:
        for t,p in zip(row['normalized_column_monomial_counts'],row['each_original_row_probability_lower_factored']):
            assert F(p['multiplier_of_density_floor'])==F(1,320*t)
    assert value['mass_predictions_GeV']==[] and value['physical_particle_map']==[]
    wrong=deepcopy(bound);wrong['result']['complete_detector']['cyclic_radix']=6
    with pytest.raises(ValueError):m.verify_bound(member,wrong)


def test_receipt_schema_and_same_source_ids_are_serializable():
    import json
    bound=m.bind_member(m.contract.member(0))
    assert json.loads(json.dumps(bound))==bound
    assert bound['result_sha256']==m.contract.digest(bound['result'])
    for key,value in m.contract.member(0)['ids'].items():assert bound['result'][key]==value
    assert not bound['result']['reference_operator_valued_ward_zero_asserted']
