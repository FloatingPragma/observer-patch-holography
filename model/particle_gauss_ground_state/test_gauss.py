"""Exact local constraints, physical cooling and retained charged outcomes."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from . import gauss as g


def test_actual_native_source_has_zero_trace_and_a_physical_blank():
    seed=g.native_spatial_seed()
    assert seed['CAR_modes_per_cell']==60
    assert seed['native_generator_count']==12
    assert seed['normal_ordering_charge_offset']==0
    assert all(s.trace(t)==0 for t in g.native.native_currents())
    assert seed['empty_Fock_constant_link_seed_physical']
    assert seed['filled_60_mode_constant_link_seed_physical']


def test_source_character_loop_survives_Gauss_and_has_nonzero_same_H_coupling():
    row=g.physical_loop_seed()
    assert row['native_irreducible_character_count']==5
    assert row['blank_loop_overlap']==0
    assert row['real_trace_loop_norm_squared']==F(5,2)
    assert row['normalized_blank_loop_H_matrix_element_squared']==F(1,4410)
    assert g.physical_loop_seed(magnetic='0')['normalized_blank_loop_H_matrix_element_squared']==0
    assert not row['protected_topological_defect_or_particle_claimed']


def test_actual_nonabelian_endpoint_projectors_commute_but_components_do_not():
    c=g.endpoint_control();px=c['left_projector'];py=c['right_projector']
    assert c['local_kernel_ranks']==[4,4]
    assert s.trace(px*py)==1
    assert px*py==py*px==c['singlet']*c['singlet'].H
    assert all(x*y==y*x for x in c['left_generators'] for y in c['right_generators'])
    assert any(x*y!=y*x for x in c['left_generators'] for y in c['left_generators'])
    # Removing an endpoint constraint admits states with a nonzero local charge.
    assert any(t*px!=s.zeros(16) for t in c['right_generators'])
    assert not c['Hamiltonian_reducing_sector_claimed']


def test_exact_kernel_method_keeps_full_invariant_multiplicity():
    t=s.diag(s.I,-s.I,0,0)
    r=g.exact_kernel_projector((t,))
    assert r['rank']==2 and r['projector']==s.diag(0,0,1,1)
    assert g.exact_kernel_projector((s.I*s.eye(2),))['rank']==0
    with pytest.raises(ValueError):g.exact_kernel_projector((s.eye(2),))
    with pytest.raises(ValueError):g.exact_kernel_projector(())
    with pytest.raises(ValueError):g.exact_kernel_projector((s.diag(0.0,s.I),))
    p=s.diag(1,0);q=s.ones(2)/2
    with pytest.raises(ValueError):g.joint_projector((p,q))


def test_source_bounds_feed_restricted_cooling_without_replacing_H():
    plan=g.spatial_plan(side=3,link_degree=2,electric='105',target_excess='1/3')
    assert plan['source']['CAR_modes']==60*27
    assert plan['cooling']['mean_ground_energy_excess_upper']<=F(1,3)
    assert plan['source']['dimension_bits_upper']==60*27+9*2*3*27
    assert plan['reference_separately_Gauss_invariant']
    assert not plan['H_or_clock_changed_by_preparation']
    assert not plan['postselection']
    with pytest.raises(ValueError):g.spatial_plan(regulator='smooth')


def test_native_gauss_support_is78_and_monitor_is_complete():
    p=g.native_projector();outside=g.native_gauss_monitor()[1]
    h=g.thermal.exact_control()['H'];gauss=g.thermal.exact_control()['Gauss']
    assert s.trace(p)==78 and s.trace(outside)==130
    assert p+outside==s.eye(208) and p*outside==s.zeros(208)
    assert p*h==h*p and gauss*p==s.zeros(208)
    assert p[6,6]==0  # Empty Fock is charged in the explicitly neutralized control.
    assert p[3*13+6,3*13+6]==1
    for flux_index in range(13):
        assert sum(p[mask*13+flux_index,mask*13+flux_index] for mask in range(16))==6


def test_native_physical_polynomial_state_support_and_stationarity():
    prepared=g.prepared_control();p=g.native_projector();h=g.thermal.real_computation_chart()[1]
    support=[];z=0
    for indices,q in prepared['amplitude_blocks']:
        hb=h.extract(indices,indices);pb=p.extract(indices,indices)
        assert pb*q==q and q*pb==q and hb*q==q*hb
        z+=s.trace(q*q)
        if q!=s.zeros(len(indices)):support.extend(indices)
    assert len(support)==78 and z==prepared['norm_squared']
    assert len(prepared['amplitude_blocks'])==82
    assert not prepared['postselected_full_carrier_run']
    assert prepared['reference_separately_Gauss_invariant']


def test_physical_energy_against_independent78_matrix_Gibbs():
    prepared=g.prepared_control();p=g.native_projector()
    indices=[i for i in range(208) if p[i,i]]
    matrix=np.array(g.thermal.exact_control()['H'].extract(indices,indices),dtype=complex)
    eigenvalues=la.eigvalsh(matrix)
    weights=np.exp(-4*(eigenvalues-eigenvalues[0]));weights/=weights.sum()
    mean=float(prepared['mean_energy'])
    assert abs(mean-weights@eigenvalues)<=26*float(prepared['preparation_trace_error_upper'])
    summary=g.native_summary()
    assert summary['physical_ground_interval'][0]<=eigenvalues[0]<=summary['physical_ground_interval'][1]
    assert 0<summary['mean_physical_ground_excess_upper']<F(11,100)


def test_charged_probe_needs_boundary_reference_and_preserves_original_sum():
    cert=g.compensated_probe_certificate();p=g.native_projector()
    assert len(cert['rows'])==8
    assert cert['compensated_total_charge_zero']
    assert not cert['bare_charged_operator_preserves_source_physical_space']
    plus=s.zeros(208);minus=s.zeros(208)
    for row in g.records.instrument():
        k=row['K']*p
        if row['sign']==1:plus+=k.H*k
        else:minus+=k.H*k
    assert plus==minus==p/2


def test_restricted_ground_inequality_cannot_be_applied_to_bare_charged_operator():
    # Independent counterexample to dropping the physical-domain premise.
    h=s.diag(0,-1);p=s.diag(1,0);a=s.Matrix([[0,0],[1,0]])
    assert s.trace(p*h)==0  # Exact restricted ground; epsilon_physical=0.
    assert s.trace(p*a.H*(h*a-a*h))==-1
    assert a*p!=p*a
    # The enlarged accessible ground is-1, restoring the correct excess1.
    assert s.trace(p*(h+s.eye(2)))==1


def test_complete_clock_records_use_physical_preparation_without_conditioning():
    law=g.records.complete_filter_records(g.prepared_control())
    assert len(law['records'])==16
    assert sum(row['original_probability'] for row in law['records'])==1
    for sign in (-1,1):
        assert sum(row['original_probability'] for row in law['records'] if row['sign']==sign)==F(1,2)
    assert law['all_quantum_references_retained']
    assert not law['negative_bin_is_negative_energy_measurement']


def test_certificate_is_json_safe_and_preserves_physical_boundary():
    c=g.certificate(include_control=False)
    json.dumps(c,allow_nan=False)
    assert c['fixed_cutoff_physical_ground_subsequence']
    assert not c['charged_local_CAR_is_gauge_invariant']
    assert not c['physical_spacing_refinement_or_scattering_proved']
