"""Independent finite-incidence, boundary, Gauss and complete-record checks."""
from fractions import Fraction as F
import itertools
import math
import sympy as s
import importlib.util
from pathlib import Path

from . import model as source


def test_cone_matches_independently_rebuilt_simplicial_boundaries():
    import verify_cone_whitney_bridge as independent
    vertices,edges,faces=independent.source_mesh()
    actual=source.geometry();D,C=actual[3:5]
    cone_edges=[(0,i+1) for i in range(12)]+[(a+1,b+1) for a,b in edges]
    cone_faces=[tuple(i+1 for i in f) for f in faces]+[(0,a+1,b+1) for a,b in edges]
    d=s.Matrix(independent.coboundary([(i,) for i in range(13)],cone_edges))
    c=s.Matrix(independent.coboundary(cone_edges,cone_faces))
    assert d==D and c==C
    changed=c.copy();changed[0,12]*=-1
    assert changed*d!=s.zeros(50,13)


def test_static_green_costs_are_exact_and_reference_ground_is_distinguished():
    rows,V=source.charge_geometry()
    assert rows[1]['static_electric_energy']==0
    assert rows[0]['static_electric_energy']==s.Rational(26,217)
    assert rows[2]['static_electric_energy']==s.Rational(33,217)
    assert rows[6]['static_electric_energy']==s.Rational(40,217)
    assert rows[12]['static_electric_energy']==s.Rational(6,31)
    for row in rows:
        assert row['total_energy_above_actual_neutral_ground']==2+row['static_electric_energy']
        assert row['source_energy_excluding_reference']+1==row['total_energy_above_actual_neutral_ground']
    assert V[source.REFERENCE_VERTEX,:]==s.zeros(1,13)


def test_literal_two_mode_total_cutoff_completion_and_heisenberg_equations():
    # Nondegenerate frequencies make the weighted boundary term necessary.
    cutoff=3;omega=(s.Integer(1),s.Integer(2))
    basis=[n for n in itertools.product(range(4),repeat=2) if sum(n)<=cutoff]
    ids={n:i for i,n in enumerate(basis)};dimension=len(basis)
    a=[]
    for j in range(2):
        b=s.zeros(dimension)
        for col,n in enumerate(basis):
            if n[j]:
                after=list(n);after[j]-=1;b[ids[tuple(after)],col]=s.sqrt(n[j])
        a.append(b)
    bare=s.zeros(dimension);pairs=[]
    for j,b in enumerate(a):
        q=(b+b.H)/s.sqrt(2*omega[j]);p=s.I*s.sqrt(omega[j]/2)*(b.H-b)
        bare+=(p*p+omega[j]**2*q*q)/2;pairs.append((q,p))
    exact=s.diag(*[sum(omega[j]*(n[j]+s.Rational(1,2)) for j in range(2)) for n in basis])
    boundary=s.diag(*[(sum(omega[j]*(n[j]+1) for j in range(2))/2 if sum(n)==cutoff else 0) for n in basis])
    assert s.simplify(bare+boundary-exact)==s.zeros(dimension)
    assert bare!=exact
    for j,(q,p) in enumerate(pairs):
        assert s.simplify(exact*q-q*exact+s.I*p)==s.zeros(dimension)
        assert s.simplify(exact*p-p*exact-s.I*omega[j]**2*q)==s.zeros(dimension)
        assert s.trace(q*p-p*q)==0


def test_all_native_size_cutoff_columns_and_bare_ground_obstruction():
    result=source.total_cutoff_control()
    assert result['photon_dimension']==5456
    assert result['mode_boundary_coefficient_checks']==30*5456
    assert result['bare_cutoff_top_state_strictly_below_bare_vacuum']
    assert result['bare_top_state_energy']<result['actual_completed_vacuum_energy']
    assert source.cutoff_control(2)['bare_quadrature_diagonal_at_frequency_one']==[s.Rational(1,2)]*2


def test_full_small_kinematic_gauss_extension_including_unphysical_inputs():
    # Three matter modes, four reference counts and two longitudinal bits.
    # Check every full-carrier nonzero column, not only matched physical inputs.
    count=0;bad_uncompensated=0
    for mask,R,ell in itertools.product(range(8),range(4),range(4)):
        n=[mask>>j&1 for j in range(3)]
        load=[-(ell&1)-((ell>>1)&1),ell&1,(ell>>1)&1]
        gauss=[load[j]-n[j]+R*int(j==0) for j in range(3)]
        for y in range(3):
            if n[y] or R==3 or (y and ell>>(y-1)&1):continue
            target_n=n.copy();target_n[y]=1
            target_ell=ell if y==0 else ell|(1<<(y-1))
            target_load=[-(target_ell&1)-((target_ell>>1)&1),target_ell&1,(target_ell>>1)&1]
            target_gauss=[target_load[j]-target_n[j]+(R+1)*int(j==0) for j in range(3)]
            assert target_gauss==gauss
            bad=[target_load[j]-target_n[j]+R*int(j==0) for j in range(3)]
            assert bad!=gauss
            count+=1;bad_uncompensated+=1
    assert count==bad_uncompensated>0
    native=source.charge_car_check()
    assert native['Gauss_preserving_creation_columns']==13*2**12


def test_complete_field_and_joint_charge_records_keep_all_failures():
    rows=source.field_records()
    for row in rows:
        assert len(row['energy_lines'])==8
        assert row['pump_success_probability']>0 and row['pump_failure_probability']>0
        assert row['probabilities_sum_exactly_one']
        residue=sum(x['unconditional_pump_and_energy_line_probability'] for x in row['energy_lines'])
        assert source.simple(residue+row['pump_failure_probability'])==1
    combined=source.joint_charge_field_records()
    assert len(combined)==26
    for row in combined:
        success=row['elastic_unconditional_probability']+sum(x['unconditional_probability'] for x in row['inelastic_lines'])
        assert source.simple(success+row['charge_pump_failure_probability']+row['charge_success_field_pump_failure_probability'])==1
        assert success>0
        assert row['charge_success_field_pump_failure_probability']>0
        assert row['retained_pointer_dimension']==12


def test_spatial_residue_sum_rules_without_diagonalizing_field_H():
    edges,faces,vertices,D,C,G,K,PT,powers,polynomials=source.geometry()
    for row in source.field_records():
        total=sum(line['spatial_projector_weight'] for line in row['energy_lines'])
        if row['kind']=='transverse_electric_edge':expected=PT[row['index'],row['index']]
        else:expected=(C[row['index'],:]*C[row['index'],:].T)[0]
        assert source.simple(total-expected)==0
    # Boundary edge12 has strictly positive weight in every transverse band.
    row=next(x for x in source.field_records() if x['kind']=='transverse_electric_edge' and x['index']==12)
    assert all(line['vacuum_residue']>0 for line in row['energy_lines'])



def test_complete_longitudinal_ground_gap_and_physical_gap():
    row=source.exact_ground_gaps()
    assert row['complete_longitudinal_patterns_checked']==4096
    assert row['actual_full_kinematic_ground_gap']==s.Rational(26,217)
    assert row['actual_physical_ground_gap']==s.sqrt(4-s.sqrt(5))


def test_eighth_moment_does_not_transfer_past_declared_cutoff_domain():
    row=source.exact_ground_gaps()
    assert row['frequency_one_single_mode_vacuum_eighth_moment']==s.Rational(81,16)
    assert row['unrestricted_Gaussian_eighth_moment']==s.Rational(105,16)
    assert row['higher_Weyl_or_all_moment_transfer_rejected']



def test_fresh_scientific_calculation_matches_preserved_independent_receipt():
    from .certificate import compare_research,scientific_fields
    value=source.calculate()
    assert compare_research(value)['all_scientific_fields_equal']
    import copy
    import pytest
    corrupted=copy.deepcopy(value)
    key='architecture' if 'architecture' in value else 'source'
    corrupted[key]['clock_mutation']='a different source clock'
    with pytest.raises(ValueError):compare_research(corrupted)
