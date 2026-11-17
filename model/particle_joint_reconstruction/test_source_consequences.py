"""Reject incomplete native evidence and transfer of diagnostic physics."""
from copy import deepcopy
import json
from pathlib import Path
import pytest

from tools.build_ccf_lt2_claims import validate_receipt

ROOT=Path(__file__).resolve().parents[2]


@pytest.fixture(scope='module')
def receipt():
    return json.loads((ROOT/'model/ccf_lt2/RESULTS.json').read_text())


def test_current_complete_native_receipt(receipt):
    validate_receipt(receipt)


@pytest.mark.parametrize('mutation', ['parent','component','polynomial','gauss','state','physical','neutral','column'])
def test_native_receipt_rejects_omissions_and_promotions(receipt,mutation):
    r=deepcopy(receipt)
    if mutation=='parent': r['source_edge']['parent_export_digest']='0'*64
    if mutation=='component': del r['source_triangle']['components']['QQ_nonleptonic']
    if mutation=='polynomial': r['source_edge']['combined_source']['matter_digest']='0'*64
    if mutation=='gauss': r['source_edge']['joint']['all_local_gauss_identities']=False
    if mutation=='state': r['source_triangle']['combined_source']['full_source_state_vector_executed']=True
    if mutation=='physical': r['source_edge']['physical_rows_passed_by_this_check']=[26]
    if mutation=='neutral':
        cols=r['actual_local_capacity_join']['columns']
        cols[1]['diagonal_exact']=cols[0]['diagonal_exact']
    if mutation=='column': r['actual_local_capacity_join']['columns'][0]['no_native_occupation_projection']=False
    with pytest.raises(ValueError): validate_receipt(r,check_bindings=False)


def test_component_results_do_not_discharge_master_physics():
    result=json.loads((Path(__file__).parent/'results.json').read_text())
    assert result['primary_program']=='CCF-LT-2'
    for key in ['full_native_state_propagation_claimed','full_particle_spectrum_derived',
                'interacting_relativistic_primary_limit_proved','proper_geometry_or_Einstein_limit_proved',
                'pure_Yang_Mills_gap_proved']:
        assert result[key] is False
    assert result['source_derived_particle_mass_table_GeV']==[]
    assert result['finite_primary_acceptance_count']==3


def test_receipt_rejects_omitted_scientific_binding(receipt):
    damaged=deepcopy(receipt)
    damaged['scientific_inputs_sha256'].pop(next(iter(damaged['scientific_inputs_sha256'])))
    with pytest.raises(ValueError,match='inventory omitted'):
        validate_receipt(damaged)


@pytest.mark.parametrize('mutation',['omit','sectors','reference','coframe','ward','physical'])
def test_bridge_evidence_cannot_omit_required_domains(mutation):
    from tools.build_ccf_lt2_claims import validate_bridges
    r=json.loads((ROOT/'model/ccf_lt2/bridges/RESULTS.json').read_text())
    validate_bridges(r)
    if mutation=='omit':del r['matter']
    if mutation=='sectors':r['matter']['edge']['conserved'].remove('L')
    if mutation=='reference':r['reference']['actual_source']['members'][1]['geometry_dimension']=48
    if mutation=='coframe':r['gravity']['coframe']['identity_member_exact']=False
    if mutation=='ward':r['gravity']['energy_ward']['omitted_current_residual']=0
    if mutation=='physical':r['kinematics']['physical_rows_claimed_solved']=[26]
    with pytest.raises(ValueError):validate_bridges(r)
