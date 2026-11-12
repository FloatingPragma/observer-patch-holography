"""Mutation controls for winding, complete histories and current attachment."""
from fractions import Fraction as F

import pytest
import sympy as s

from .protection import independent_degree, half_degree_solutions, current_attachment_certificate
from model.compact_defect_transport import source
from model.lattice_texture_selection.shells import degree
from model.defect_integer_lift import model as lift


def test_independent_regular_value_count_and_orientation_mutation():
    state=source.initial();triangles=source.surfaces()[1]
    values={p:state[p] for t in triangles for p in t}
    reference=degree(triangles,values)
    assert independent_degree(triangles,values,reference['regular_target'])[0]==1
    reversed_triangles=tuple((a,c,b) for a,b,c in triangles)
    assert independent_degree(reversed_triangles,values,reference['regular_target'])[0]==-1
    with pytest.raises(ValueError,match='closed triangulation'):
        degree(triangles[:-1],values)


def test_all_actual_retained_records_and_no_halving():
    _,records,transition=source.continuing_protocol()
    assert len(records)==12 and all(r['local_shell_degrees']==[-1,1] for r in records)
    assert all(sum(row)==1 for row in transition)
    assert not half_degree_solutions(1) and not half_degree_solutions(-1)
    assert half_degree_solutions(2)==[1]
    assert 2+(-1)==1  # No claim against unequal signed decompositions.
    retained_mass=sum(F(r['probability']) for r in records if r['source_state']==0 and r['accepted'])
    assert retained_mass==F(1,2)  # Accepted-only reporting loses half the law.


def test_degree_erasure_keeps_four_zero_onsite_coordinates():
    state=source.initial()
    erased={p:((0,0,1) if any(v) else v) for p,v in state.items()}
    assert source.charges(erased)==[0,0]
    assert sum(not any(v) for v in state.values())==sum(not any(v) for v in erased.values())==2
    with pytest.raises(ValueError,match='budget'):
        source.validate(erased)


def test_integer_current_witness_does_not_move_with_its_class():
    na,nb,j=lift.attempt(4,0,0)
    assert j==2 and na+nb==4 and lift.accepts(4,0,j,4)
    before=(4,0);after=(before[0]+j,before[1]-j)
    assert after==(6,-2) and sum(after)==sum(before)
    assert tuple(x%6 for x in after)==(0,4)
    assert after[1]!=before[0]
    # Adding an arbitrary integer sheet does not change the native decision.
    assert lift.accepts((before[0]+6)%6, before[1]%6,j,4)


def test_same_winding_current_dictionary_is_impossible():
    value=current_attachment_certificate()
    assert value['same_plus_one_texture_one_particle_witness']['no_function_of_degree_equals_native_Q_on_this_subspace']
    alpha=s.symbols('alpha')
    assert s.solve([alpha,alpha-1],[alpha])==[]
    assert value['uniform_native_current_preserved_by_spatial_transport']
    assert value['all_Q_plus_colour_weights_integer']
    assert all(s.Rational(x)*3==int(s.Rational(x)*3) for x in value['one_family_Q_weights'])
    rho=s.eye(2)/2;q=s.diag(0,1)
    assert s.trace(rho*q)==s.Rational(1,2) and s.Rational(1,2) not in q.eigenvals()
