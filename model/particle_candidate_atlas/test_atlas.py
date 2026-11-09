"""Data-integrity, uncertainty and physical-scope controls for the atlas."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pytest

from . import builder as b
from .graphic import render


def test_exact_uncertainty_comparison_retains_numerical_width():
    r=b.compare(['1001/1000','1003/1000'],'1','1/100')
    assert r['numerical_width']=='1/500'
    assert r['residual_over_quoted_uncertainty']==['1/10','3/10']
    assert r['inside_quoted_window'] and not r['joint_statistical_significance']


def test_whole_interval_must_lie_inside_reference_window():
    assert not b.compare(['99/100','102/100'],'1','1/100')['inside_quoted_window']


@pytest.mark.parametrize('interval,error', [(['2','1'],'1'),(['1','2'],'0'),(['1','2'],'-1')])
def test_invalid_numerical_bounds_rejected(interval,error):
    with pytest.raises(ValueError): b.compare(interval,'1',error)


def test_authentication_rejects_changed_and_escaped_inputs(tmp_path):
    p=tmp_path/'input';p.write_text('exact')
    pins={'input':b.digest(p)};b.authenticate(tmp_path,pins)
    p.write_text('changed')
    with pytest.raises(ValueError,match='stale'): b.authenticate(tmp_path,pins)
    with pytest.raises(ValueError): b.authenticate(tmp_path,{'../outside':'0'*64})


def test_reference_conventions_are_not_flattened_to_common_poles():
    refs=b.load_json(b.ROOT/b.REFS);rows={r['target']:r for r in refs['masses']}
    assert len(rows)==18
    for t in ('up','down','strange','charm','bottom'):
        assert rows[t]['confidence_level']=='90%'
        assert 'MSbar' in rows[t]['definition']
    assert 'Breit-Wigner' in rows['Z']['definition']
    assert rows['photon']['upper_bound'] and rows['gluon']['value']=='0'
    for t in ('electron_neutrino','muon_neutrino','tau_neutrino'):
        assert rows[t]['value'] is None
    assert not any(r['established_source_physical_mass_prediction'] for r in rows.values())
    assert all(r['unit'] in ('GeV',None) for r in rows.values())
    assert F(rows['electron']['value'])==F('0.51099895069')/1000
    assert F(rows['electron']['uncertainty'])==F('0.00000000016')/1000
    assert F(rows['tau']['value'])==F('1776.93')/1000
    assert F(rows['photon']['upper_bound'])==F(1,10**27)


@pytest.fixture(scope='module')
def source():
    # Authenticated receipt read, without repeating any scientific producer.
    return b.load_inputs()


@pytest.fixture(scope='module')
def joins(source):
    return b.science_projection(source[0],source[2],tuple(sorted(source[-1].items())))


def test_actual_data_preserve_calibration_and_distinct_actions(source):
    data=b.project(*source)
    assert len(data['lepton_branches'])==5
    assert data['lepton_branches'][0]['positive_root_Koide_ratio']==['2/3','2/3']
    assert all(b.box(r['positive_root_Koide_ratio'])[1]<F(2,3) for r in data['lepton_branches'][1:])
    assert data['independent_fourier_spectrum_diagnostic']['branches_checked']==5
    assert data['physical_particle_map']==[]
    assert not data['canonical_child_export_consumed']  # pure projection claims no authentication


@pytest.mark.parametrize('mutation', ['physical','defect','target','candidate_export','weight','branch','phase_fit'])
def test_source_scope_forgeries_are_rejected(source,mutation):
    child,candidate,parent,refs,pins=copy.deepcopy(source)
    if mutation=='physical':candidate['rows'][0]['physical_mass']='0.511'
    elif mutation=='defect':candidate['separate_defect_inventory']['species_links']=['electron']
    elif mutation=='target':candidate['rows'][-1]['target']='graviton'
    elif mutation=='candidate_export':candidate['rows'][0]['source_object']='unbound'
    elif mutation=='weight':parent['source_weights'][0]['Q']='0'
    elif mutation=='branch':child['certificates']['conditional_P_dependence']['branches'].pop()
    elif mutation=='phase_fit':child['certificates']['calibrated_sanity_test']['phase']['tau_used_to_choose_phase']=True
    # Matched copies must still fail physical/coverage admission, not only equality.
    if mutation in ('physical','defect','target'):
        child['certificates']['candidate_particle_map']=copy.deepcopy(candidate)
    with pytest.raises(ValueError): b.project(child,candidate,parent,refs,pins)


def test_svg_links_and_readable_geometry(source,joins):
    data=b.project(*source)
    data['finite_joins']=joins
    svg=ET.fromstring(render(data));ns={'s':'http://www.w3.org/2000/svg'}
    assert svg.attrib['width']=='4000'
    assert 5500<int(svg.attrib['height'])<12000
    assert all(float(t.attrib['font-size'])>=24 for t in svg.findall('.//s:text',ns))
    metadata=json.loads(svg.find('s:metadata',ns).text)
    assert metadata['input_sha256']==source[-1]
    for a in svg.findall('.//s:a',ns):
        href=a.attrib['href']
        if not href.startswith('https://'):
            assert (b.OUTPUT.parent/href).resolve().is_file(),href
    text=' '.join(t.text or '' for t in svg.findall('.//s:text',ns))
    assert '1 : 1 : (34 + 24√2)' in text
    assert 'all eight valleys' in text
    assert 'not statistical significance' in text
    assert 'Complete conditional Standard Model benchmark' in text
    assert 'Pixel self-consistency selects χ' in text
    assert '1/40000' in text and '8,8,4,4' in text
    assert 'MeV' not in text


def test_original_probabilities_and_no_physical_promotion_in_bounded_joins(joins):
    data=joins
    assert not data['complete_native_experiments_replayed_here']
    assert not data['large_operator_certificates_replayed_here']
    assert data['physical_particle_map']==data['physical_mass_predictions']==[]
    native=data['native_unread_response']
    weights=list(map(F,native['all_four_coarse_record_weights']))
    assert len(weights)==4 and min(weights)>0 and sum(weights)==1
    assert F(native['logical_purity_loss_lower'])>0
    assert F(native['full_CAR_vacuum_purity_loss_lower'])>0
    rows=data['target_admissions']
    assert [r['target'] for r in rows]==list(b.TARGETS)
    assert all(r['physical_mass'] is None and r['defect_link'] is None for r in rows)
    assert len(data['capacity_nonselection']['rows'])==2
    sm=data['conditional_standard_model']
    assert [r['target'] for r in sm['target_rows']]==list(b.TARGETS)
    assert sm['physical_OPH_assignments']==0
    assert all(r['numerical_mass'] is None for r in sm['neutrino_mass_fields'])
    assert sm['anomalies_and_Yukawa_charges']['U1_cubed']=='0'
    assert not sm['all_reference_masses_inverted_as_tree_parameters']
    fit=data['all_three_input_controller']
    assert fit['all_three_masses_used_as_inputs'] and fit['held_out_mass_predictions']==0
    assert fit['complete_retained_flags']==16
    assert not fit['physical_mass_matching_proved'] and not fit['source_selected_controller']
    assert all(F(e)<F(1,10**36) for e in fit['absolute_mass_error_bounds_GeV'])
    for row in data['capacity_nonselection']['rows']:
        assert row['perturbed_fixed_point_existence_certified']
        assert not row['source_program_is_identical_to_original']
        assert not row['explicit_shift_root_uniqueness_certified']
