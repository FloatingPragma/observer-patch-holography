from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import pytest
from tools.joint_observer_results import NAMES,project,render,decimal_lower

@pytest.fixture(scope='module')
def data():
    root=Path(__file__).resolve().parents[1]
    return {n:json.loads((root/f'model/{n}/results.json').read_text())for n in NAMES}

def test_joint_state_identity_and_distinct_spectrum_and_refinement_domains(data):
    p=project(data)
    assert len(p['source_packages'])==9
    assert p['proved_finite_interfaces']['tracked_clock_has_complete_same_joint_state_records']
    assert p['physical_premises_closed']==p['physical_identification_ids_closed']==[]
    assert p['new_independent_GeV_rest_masses']==[]
    assert p['official_status']['validated']==0 and p['official_status']['untested']==31
    text='\n'.join(render(p))
    assert '7,320' in text and '4,096' in text and '0.088951675369' in text
    assert 'not flux-dressed particle rest masses' in text
    # Compact copies retain the authenticated table digest and complete scope.
    walk=p['source_packages']['native_gauss_walk']
    assert walk['scope']==data['native_gauss_walk']['scope']
    assert '_zlib_base64' not in json.dumps(p)
    assert data['native_gauss_walk'].get('source_sha256')

@pytest.mark.parametrize('owner,path,value',[
 ('native_gauss_clock',('scope','local_autonomous_clock_proved'),True),
 ('native_gauss_clock',('scope','physical_rest_mass_identified'),True),
 ('native_gauss_clock',('scope','same_charged_source_field_and_clock'),False),
 ('native_gauss_clock',('common_program_sha256',),'wrong'),
 ('native_gauss_walk',('integer_extension','different'),True),
 ('native_gauss_walk',('preparation','link_limit'),5),
 ('native_gauss_walk',('scope','reset_in_inverse_test'),True),
 ('native_gauss_walk',('scope','undressed_source_spectrum_transferred'),True),
 ('native_transport_spectrum',('unit_contract','physical_rest_mass_assignment'),True),
 ('native_transport_spectrum',('unit_contract','source_clock_and_walk_have_proved_common_stationary_generator'),True),
 ('native_transport_spectrum',('scope','spectrum_of_flux_dressed_joint_electric_program'),True),
 ('native_spatial_refinement',('scope','whole_interacting_quantum_field_theory_constructed'),True),
 ('native_spatial_refinement',('scope','four_dimensional_Yang_Mills_mass_gap'),True),
 ('native_spatial_refinement',('scope','Einstein_gravity_derived'),True),
 ('native_criterion_packet',('criterion1','physical_acceptance_complete'),True),
 ('native_criterion_packet',('official_authority','counts','validated'),1),
 ('native_criterion_packet',('official_authority','all31_criteria',0,'upstream_status','M10'),'✅'),
 ('native_gauss_walk',('observed_electric_phase_response','two_step_source_pair_total_variation'),'1/100'),
 ('native_gauss_walk',('observed_electric_phase_response','three_step_response_lower_with_native_event_word_bound'),'0'),
])
def test_false_interface_or_physics_promotion_rejected(data,owner,path,value):
    d=deepcopy(data);row=d[owner]
    for key in path[:-1]:row=row[key]
    row[path[-1]]=value
    with pytest.raises(ValueError):project(d)

def test_retained_matrix_units_cannot_be_replaced_by_matching_labels(data):
    d=deepcopy(data)
    d['native_gauss_clock']['exact_transported_matrix_units'][1]['complete_matrix_units']=12
    with pytest.raises(ValueError,match='matrix-unit'):project(d)

def test_positive_bound_prints_outwards(data):
    q=data['native_gauss_walk']['observed_electric_phase_response']['three_step_response_lower_with_native_event_word_bound']
    assert 0<F(decimal_lower(q))<=F(q)
