"""Scope, complete coverage and reference/rendering regressions for the spread."""
from copy import deepcopy
import hashlib
import json
import xml.etree.ElementTree as ET
import pytest
from . import build, drawing


@pytest.fixture(scope='module')
def data():return build.load()


def test_every_family_and_each_owed_bridge_survives(data):
    build.validate_content(data)
    assert len(data['families']['lanes'])==9
    assert len(data['families']['calibration_paths'])==10
    assert len(data['theorems']['roles'])==5
    assert data['v1']['physical_particle_map']==[]
    assert len(data['candidate_correspondence']['target_rows'])==18
    assert data['candidate_correspondence']['physical_particle_map']==[]
    assert len(data['candidate_identity_tests']['target_rows'])==18
    assert all(row['physical_identity_status']=='unresolved'
               for row in data['candidate_identity_tests']['target_rows'])


@pytest.mark.parametrize('kind',['lost_target','duplicated_target','promoted_identity','invented_mass','changed_source'])
def test_candidate_identity_coverage_and_source_dictionary_are_guarded(data,kind):
    bad=deepcopy(data);value=bad['candidate_identity_tests'];rows=value['target_rows']
    if kind=='lost_target':rows.pop()
    elif kind=='duplicated_target':rows[-1]['target']=rows[0]['target']
    elif kind=='promoted_identity':rows[0]['physical_identity_status']='confirmed'
    elif kind=='invented_mass':rows[0]['source_derived_mass_GeV']='0.000511'
    else:rows[0]['source_component_ids']=['unproved_component']
    with pytest.raises(ValueError):build.validate_content(bad)


@pytest.mark.parametrize('kind',['missing_target','duplicate_target','promoted_identity','promoted_edge','lost_strategy'])
def test_visual_shortcuts_cannot_silently_promote_the_science(data,kind):
    bad=deepcopy(data)
    if kind=='missing_target':bad['families']['lanes'][0]['target_ids'].pop()
    elif kind=='duplicate_target':bad['families']['lanes'][0]['target_ids'][0]='muon'
    elif kind=='promoted_identity':bad['families']['lanes'][0]['targets'][0]['source_physical_identification']=True
    elif kind=='promoted_edge':bad['story']['edges'][2]['status']='exact'
    else:bad['families']['calibration_paths'][0]['removal_strategy']=''
    with pytest.raises(ValueError):build.validate_content(bad)


def test_render_has_real_measurement_links_and_all_derivation_paths(data):
    svg=drawing.render(data)
    root=ET.fromstring(svg)
    assert float(root.get('width'))>float(root.get('height'))
    assert (float(root.get('width')),float(root.get('height'))) == (48000,47865)
    hrefs=[n.get('href') for n in root.iter('{http://www.w3.org/2000/svg}a')]
    assert any('OperationalBornBridge.lean#L214' in h for h in hrefs)
    assert any('PrimitivePortFrameQuotient.lean#L631' in h for h in hrefs)
    assert any('ModularClockComparison.lean#L40' in h for h in hrefs)
    assert sum('.pdf#page=' in h for h in hrefs)>20
    text=' '.join(''.join(root.itertext()).split())
    assert 'no single sharp value' in text
    assert 'one history: 1/8' in text
    assert 'opposite winding does not identify an antiparticle' in text
    assert text.count('Pass:')==10
    for lane in data['families']['lanes']:
        assert lane['title'] in text


def test_render_uses_certified_complete_probabilities(data):
    value=deepcopy(data)
    value['measurement']['position_probabilities']=['1/3','1/3','1/3','0']
    svg=drawing.render(value)
    # A display-only mutation must reach actual bars; it is never a source certificate.
    assert '>1/3</text>' in svg


def test_reader_embeds_the_vector_and_preserves_every_reference_anchor(data):
    text=build.reader(drawing.render(data),data)
    assert '<svg ' in text and 'fetch(' not in text
    assert 'calibrations:[32000,6150]' in text
    assert 'pointermove' in text and 'ArrowRight' in text
    assert 'ModularClockComparison.lean#L40' in text
    assert 'FiniteBuschGleason.lean#L1133' in text
    assert '>Ten dependencies</button>' in text
    assert 'max-height:40vh;overflow:auto' in text
    assert '>Best-guess map</button>' in text
    assert '>Masses & collisions</button>' in text
    assert 'bestguess:[0,30400]' in text and 'comparisons:[38000,30400]' in text
    assert 'model/particle_primary_locality/PROOF.md' in text
    assert 'model/particle_joint_reconstruction/DEPENDENCIES.md' in text
    assert 'model/particle_joint_reconstruction/CANDIDATE.md' in text
    for owner in ('particle_candidate_charged','particle_candidate_composites','particle_candidate_bosons'):
        assert 'model/'+owner+'/PROOF.md' in text
    assert 'model/particle_candidate_confirmation/RESULTS.md' in text
    for owner in ('particle_primary_survival','particle_primary_mass_readout',
                  'particle_primary_spatial_algebra','particle_primary_capacity',
                  'particle_compact_mass_probe', 'particle_primary_momentum',
                  'particle_primary_channel_labels', 'particle_primary_spectral_scale',
                  'particle_primary_endpoint_bounds', 'particle_primary_variational_readout',
                  'particle_primary_joint_spectroscopy', 'particle_primary_occupation_bound',
                  'particle_primary_residual_refinement', 'particle_primary_response_moments',
                  'particle_primary_branch_resolution', 'particle_primary_energy_tightness',
                  'particle_primary_origin_momentum', 'particle_primary_open_coframe',
                  'particle_primary_charged_string', 'particle_primary_odd_residual',
                  'particle_primary_momentum_transfer', 'particle_primary_local_fluctuation',
                  'particle_primary_spectral_resolution', 'particle_primary_charged_readout',
                  'particle_primary_staggered_response', 'particle_primary_current_tail',
                  'particle_primary_sector_threshold', 'particle_primary_populated_momentum',
                  'particle_primary_weighted_survival', 'particle_primary_charge_conjugation',
                  'particle_primary_composite_spin', 'particle_primary_composite_transport',
                  'particle_primary_selected_odd_response', 'particle_primary_covariant_stress',
                  'particle_primary_selected_spin_response', 'particle_primary_odd_resolvent',
                  'particle_primary_odd_low_energy', 'particle_primary_geometric_cut',
                  'particle_primary_coherent_odd_response',
                  'particle_primary_original_shell_moments', 'particle_primary_original_branch', 'particle_primary_critical_refinement', 'particle_primary_boundary_information',
                  'particle_primary_odd_sum_rules', 'particle_primary_branch_certificate', 'particle_primary_work_spectrum', 'particle_primary_boundary_correlation',
                  'particle_primary_reduced_sector', 'particle_primary_regional_fermion', 'particle_primary_weighted_spectral_tail', 'particle_primary_energy_response',
                  'particle_primary_loop_fusion', 'particle_primary_equal_work_control'):
        assert 'model/'+owner+'/PROOF.md' in text


def test_named_candidate_wings_preserve_exact_weights_and_separate_collision_domains(data):
    canvas=drawing.Canvas()
    build.source_connections.render_candidate_wings(canvas,data)
    root=ET.fromstring('<svg xmlns="http://www.w3.org/2000/svg">'+''.join(canvas.items)+'</svg>')
    text=' '.join(' '.join(root.itertext()).split())
    assert 'BEST-GUESS PARTICLE MAP' in text
    assert 'charged regional G1 factor' in text and 'the pair is not automatically a photon' in text
    assert 'Physical spin and named identities are unproved' in text
    assert 'S=N₋₄+N₂+N₆' in text and 'multiple excitations can alias' in text
    assert 'Higgs-conjugate scalar composite gives a dressed charged singlet' in text
    assert 'Actual polar and axial loop channels have nonzero source transitions' in text
    window=data['candidate_identity_tests']['calculated_source_results']['charged']['original_compensated_dressed_energy_window']
    assert 'original probability >'+window['original_dressed_recorded_simple_strict_lower'] in text
    assert '['+','.join(window['reported_transfer_interval'])+'] source-energy units on the same clock' in text
    assert f'All {window["complete_success_outcomes"]} success records and failure are retained' in text
    assert 'Antiparticle and flavor readouts need the same physical field dictionary' in text
    assert 'Y=1/6' in text and 'Y=-2/3' in text and 'Y=-1/2' in text
    assert 'Every family ratio <21.16; odd-family ratios <15.31' in text
    assert 'Measured μ/e =206.7682827(46)' in text
    assert 'excludes only that member-zero assignment' in text
    assert 'Each complete pair row retains original weight >c_g/40' in text
    assert '22,531 native pulses' in text and 'fine-grid incoming-packet program is not executed' in text
    assert '32-row five-qubit controls' in text
    assert 'No fitted GeV mass is added here' in text
    assert max(b for _,b,_,_ in canvas.text_boxes)<47865
    hrefs=[n.get('href') for n in root.iter('{http://www.w3.org/2000/svg}a')]
    for path in ('model/standard_model_source_map/candidates/CANDIDATES.md',
                 'model/standard_model_source_map/candidates/FINITE_SECTOR_RATIOS.md',
                 'model/particle_candidate_confirmation/RESULTS.md',
                 'model/particle_candidate_composites/PROOF.md',
                 'model/particle_candidate_charged/PROOF.md',
                 'model/contact_scattering/README.md','model/scattering_packet_readout/PROOF.md'):
        assert '../../'+path in hrefs


def test_primary_response_audit_authenticates_receipt_and_source_bytes(tmp_path,monkeypatch):
    source=tmp_path/'source.py';source.write_text('exact declared source\n')
    result={'member_id':'control-member','physical_mass_table_GeV':[]}
    digest=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    receipt={'schema':'test.source-bound.v1','source_sha256':{'source.py':build.sha(source)},
             'selected_member_bound':{'result':result,'result_sha256':digest}}
    path=tmp_path/'receipt.json';path.write_text(json.dumps(receipt))
    monkeypatch.setattr(build,'ROOT',tmp_path)
    monkeypatch.setattr(build,'PRIMARY_RESPONSE_RECEIPTS',('receipt.json',))
    rows,paths=build.scientific_receipt_records()
    assert paths=={'source.py','receipt.json'}
    assert rows[0]['member_result_bindings']==[{'member_id':'control-member','result_sha256':digest}]
    assert rows[0]['scientific_experiments_replayed_by_atlas'] is False
    receipt['embedded_duplicate']=deepcopy(receipt['selected_member_bound'])
    path.write_text(json.dumps(receipt))
    duplicate_rows,_=build.scientific_receipt_records()
    assert duplicate_rows[0]['member_result_bindings']==rows[0]['member_result_bindings']
    receipt['embedded_duplicate']['result']['physical_mass_table_GeV']=[1]
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='binding differs'):build.scientific_receipt_records()
    receipt.pop('embedded_duplicate')
    receipt['selected_member_bound']['result']['physical_mass_table_GeV']=[1]
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='binding differs'):build.scientific_receipt_records()
    receipt['selected_member_bound']['result']['physical_mass_table_GeV']=[]
    path.write_text(json.dumps(receipt));source.write_text('different source\n')
    with pytest.raises(ValueError,match='stale primary response input'):build.scientific_receipt_records()


def test_generated_outputs_match_their_exact_source_inventory():
    svg,html,audit=build.build()
    assert build.SVG.read_text()==svg
    assert build.HTML.read_text()==html
    assert build.AUDIT.read_text()==audit
    assert json.loads(audit)['source_physical_identifications']==[]
    identity=json.loads(audit)['candidate_identity_tests']
    assert identity['target_count']==18 and len({row['target'] for row in identity['targets']})==18
    assert all(row['physical_identity_status']=='unresolved' for row in identity['targets'])
    assert identity['physical_particle_map']==[] and identity['mass_predictions_GeV']==[]
    assert not identity['source_scientific_experiments_replayed_by_atlas']
    assert not identity['original_dressed_pair_window']['energy_window_is_a_single_endpoint_or_rest_mass']
