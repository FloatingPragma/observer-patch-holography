"""Actual gauge kernels, exterior rotations, native charges and complete records."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import model as m


def signed_sectors():
    return [tuple(sign*n for n in d) for _,d,_ in m.sectors.SECTORS for sign in (1,-1)]


def test_complete_native_kernel_dimensions_and_Casimir_are_computed():
    for index,delta in enumerate(signed_sectors()):
        local=m.local_multiplet(delta)
        assert len(local['masks'])==[12,12,90,90,40,40][index]
        assert local['kernel_basis'].cols==[6,6,6,6,8,8][index]
        expected={s.Rational(3,4):2,s.Rational(15,4):4} if index<4 else {s.Integer(2):3,s.Integer(6):5}
        assert local['casimir_multiplicities']==expected
        assert sum(local['casimir_projectors'].values(),s.zeros(local['kernel_basis'].cols))==s.eye(local['kernel_basis'].cols)
        for P in local['casimir_projectors'].values():
            assert m.clean(P*P)==P and P.H==P
    with pytest.raises(ValueError):m.occupation_basis((0,0,0,0,1))


def test_zero_central_charge_is_insufficient_for_a_physical_color_state():
    local=m.local_multiplet((1,0,0,2,0));basis=s.eye(90)[:,0]
    assert any(A*basis!=s.zeros(90,1) for A in local['gauge_matrices'])
    assert all(A*local['kernel_basis']==s.zeros(90,6) for A in local['gauge_matrices'])


def test_actual_exterior_cube_decomposition_and_inversion_labels():
    for i,delta in enumerate(signed_sectors()):
        row=m.cube_restriction(delta)
        expected={'G1':1,'H':1} if i<4 else {'E':1,'T1':1,'T2':1}
        assert row['proper_binary_cube_irreps']==expected
        suffix='g' if i<4 and i%2==0 else 'u'
        assert row['uniform_orbit_full_cube_irreps']=={key+suffix:1 for key in expected}
        local=m.local_multiplet(delta);d=local['kernel_basis'].cols
        for q,D in zip(m.group.elements(),row['matrices']):
            if tuple(map(m.channels._sym,q))==(-1,0,0,0):assert D==(-1)**local['CAR_parity']*s.eye(d)


def test_literal_group_multiplication_on_actual_CAR_kernel():
    delta=(0,2,0,0,1);elements=m.group.elements()
    generator=next(q for q in elements if tuple(map(m.channels._sym,q))==(s.sqrt(2)/2,0,0,s.sqrt(2)/2))
    A=m.exterior_rotation(delta,generator)
    for q in elements:
        assert m.clean(A*m.exterior_rotation(delta,q))==m.exterior_rotation(delta,m.group.multiply(generator,q))


def test_literal_charge_conjugation_matches_every_cube_rotation_and_opposite_sector():
    for _,delta,_ in m.sectors.SECTORS:
        other=tuple(-n for n in delta);C=m.conjugate_multiplet(delta)['kernel_map']
        reverse=m.conjugate_multiplet(other)['kernel_map']
        assert m.clean(reverse*C)==s.eye(C.cols)
        for q in m.group.elements():
            assert m.clean(C*m.exterior_rotation(delta,q)-m.exterior_rotation(other,q)*C)==s.zeros(C.rows,C.cols)
        for mask in m.local_multiplet(delta)['masks']:
            target,a=m.charge_mask(mask);back,b=m.charge_mask(target)
            assert back==mask and a*b==1


def test_composite_operators_reproduce_normalized_kernel_columns():
    for delta in signed_sectors():
        local=m.local_multiplet(delta)
        for col in range(local['kernel_basis'].cols):
            out={}
            for ops,value in m.composite_polynomial(delta,col):
                target,phase=m.apply_operators(m.OMEGA,ops);out[target]=out.get(target,0)+phase*value
            assert s.simplify(sum(s.conjugate(x)*x for x in out.values()))==1
            assert all(s.simplify(out.get(mask,0)-local['kernel_basis'][row,col])==0
                       for row,mask in enumerate(local['masks']))
        control=m.operator_controls(delta)
        assert control['native_current_commutators_exact'] and control['gauge_invariant_composite_operator']
        assert not control['local_canonical_free_CAR_asserted']


def test_distinct_site_exchange_with_adjoint_and_unequal_composite_families():
    first=m.composite_polynomial((0,2,0,0,1),0)[0][0]
    second=m.composite_polynomial((1,0,0,2,0),0)[0][0]
    shifted=tuple((i+60,c) for i,c in second)
    adjoint=tuple((i,not c) for i,c in reversed(shifted))
    anchor=m.OMEGA|(m.OMEGA<<60)
    occupied,_=m.apply_operators(anchor,shifted)
    a,p=m.apply_operators(occupied,first+adjoint)
    b,q=m.apply_operators(occupied,adjoint+first)
    assert a==b and p==-q and p!=0
    # One trilinear composite is not a canonical free creation operator.
    ops=((0,True),(1,True),(2,True));adj=tuple((i,not c) for i,c in reversed(ops))
    def anticommutator(mask):
        out=0
        for sequence in (ops+adj,adj+ops):
            target,phase=m.apply_operators(mask,sequence)
            if target==mask:out+=phase
        return out
    assert anticommutator(0)==1 and anticommutator(1)==0 and anticommutator(7)==1


def test_uniform_cube_orbit_has_no_spectator_exchange_sign():
    cells=list(product(range(2),repeat=3))
    for q in m.group.elements():
        targets=[m.channels.rotate_cell(x,2,m.group.rotation(q)) for x in cells]
        assert sorted(targets)==cells
    assert m.OMEGA.bit_count()==30
    for delta in signed_sectors():
        local=m.local_multiplet(delta)
        assert all((mask.bit_count()-30)%2==local['CAR_parity'] for mask in local['masks'])


def test_bare_mean_and_positive_complete_cube_energy_records():
    receipt=m.bind_member(m.contract.member(0))
    for family in receipt['result']['multiplet_families']:
        for row in family['signed_multiplets']:
            d=sum(abs(int(n)) for n in row['multiplet']['delta'])
            assert F(row['bare_source_energy_mean'])==-720+3*d
            low,high=map(F,row['each_listed_cube_sector_bottom_excitation_interval'])
            assert 0<low<high
            record=row['prepared_joint_cube_energy_readout']
            assert F(record['reported_window_original_weight_lower'])>=F(255,512)
            assert F(record['reported_excitation_window'][0])>0
            assert record['all_cube_irreps_and_energy_records_retained']
            assert record['source_marginal_changes']
            assert record['original_state_and_reference_correlations_retained_in_memory']
            assert record['compensating_opposite_occupation_sector_retained_in_second_memory']
    assert receipt['result']['physical_mass_table_GeV']==[]
    assert not receipt['result']['named_standard_model_species_identified']


def test_binding_rejects_wrong_clock_and_species_relabelling():
    member=m.contract.member(0);value=m.bind_member(member);m.verify_bound(member,value)
    changed=deepcopy(value);changed['result']['named_standard_model_species_identified']=True
    with pytest.raises(ValueError):m.verify_bound(member,changed)
    bad=deepcopy(member);bad['source']['spacing']='2'
    with pytest.raises(ValueError):m.bind_member(bad)
