"""Exact target coverage, source ancestry and rejection of false assignments."""
from copy import deepcopy
import json
import pytest

from tools import strict_json
from . import candidate_map as c
from . import p_dependence


@pytest.fixture(scope='module')
def inputs():
    return (strict_json.load(c.ROOT/c.PARENT),c.calibrated_summary(),
            p_dependence.certificate(),strict_json.load(c.ROOT/c.DEFECTS))


@pytest.fixture(scope='module')
def result(inputs):
    return c.summarize(*inputs)


def test_all_eighteen_targets_and_fifteen_internal_weights(result):
    assert result['target_count']==18 and result['source_internal_weight_count']==15
    assert tuple(row['target'] for row in result['rows'])==c.TARGETS
    assert len({r['target'] for r in result['rows']})==18
    assert [r['target_kind'] for r in result['rows']].count('neutrino_flavor_channel')==3


def test_charged_family_names_are_only_sorted_mode_hypotheses(result):
    rows=result['rows'][:3]
    assert [r['candidate_sorted_charged_mode_index'] for r in rows]==[0,1,2]
    assert all(r['candidate_name_is_a_hypothesis'] for r in rows)
    assert all(r['charged_source_branch_kinds']==['calibrated_coherent','P_MCPR'] for r in rows)
    assert all(r['physical_mass'] is None for r in rows)


def test_quark_neutrino_boson_action_masses_are_not_invented(result):
    assert all(row['source_mass_coefficient'] is None for row in result['rows'][3:])
    assert all(row['candidate_sorted_charged_mode_index'] is None for row in result['rows'][3:])
    assert next(r for r in result['rows'] if r['target']=='gluon')['color_components']==8


def test_every_physical_and_defect_assignment_is_empty(result):
    for row in result['rows']:
        for key in ('physical_mass','physical_spin','physical_lifetime','physical_width',
                    'topological_invariant','identified_topological_defect','defect_creation_map',
                    'physical_particle_or_resonance_map'):
            assert row[key] is None
        assert row['missing_defect_bridge'] and row['missing_physical_links']
        assert row['physical_particle_identified'] is False
    assert result['physical_particle_map']==[] and result['physical_mass_predictions']==[]


def test_calibrated_and_capacity_sources_are_not_collapsed(result):
    coherent,mcpr=result['charged_source_branches']
    assert coherent['kind']=='calibrated_coherent' and not coherent['source_phase_derived']
    assert coherent['exact_positive_root_Koide_ratio']=='2/3'
    assert mcpr['kind']=='P_MCPR' and not mcpr['same_as_calibrated_coherent_action']
    assert tuple(r['branch'] for r in mcpr['variants'])==c.P_BRANCHES
    assert all(r['complete_inverse_memory_branch_count']==2 for r in mcpr['variants'])
    assert not mcpr['individual_projector_is_a_mass_eigenmode']
    assert not mcpr['shape_square_is_absolute_mass_coefficient']


def test_separate_defects_are_not_matched_by_charge_or_count(result):
    inv=result['separate_defect_inventory']
    assert inv['supplied'] and len(inv['entries'])==11
    assert inv['species_links']==[]
    assert not inv['exhaustive_over_all_possible_sources_or_defects']


def test_no_defect_inventory_is_also_an_explicit_valid_scope(inputs):
    value=c.summarize(*inputs[:3])
    assert value['separate_defect_inventory']=={'supplied':False,'entries':[],'species_links':[]}


@pytest.mark.parametrize('kind',['missing_target','duplicate_target','wrong_component','wrong_kind',
                                'physical_parent','gauge_boson_promotion','changed_weight'])
def test_parent_scope_forgeries_rejected(inputs,kind):
    args=deepcopy(inputs);parent=args[0]
    if kind=='missing_target': parent['rows'].pop()
    elif kind=='duplicate_target': parent['rows'][-1]=deepcopy(parent['rows'][0])
    elif kind=='wrong_component': parent['rows'][0]['internal_components']={}
    elif kind=='wrong_kind': parent['rows'][9]['target_kind']='fermion_species'
    elif kind=='physical_parent': parent['rows'][0]['physical_particle_identified']=True
    elif kind=='gauge_boson_promotion': parent['gauge_and_Higgs_readout']['adjoint_directions_identified_with_propagating_bosons']=True
    else: parent['source_weights'][0]['Q']='0'
    with pytest.raises(ValueError): c.summarize(*args)


@pytest.mark.parametrize('kind',['tau_fitted','phase_derived','loop_proved','unsorted','zero_mass'])
def test_calibration_or_physical_promotion_rejected(inputs,kind):
    args=deepcopy(inputs);v=args[1]
    if kind=='tau_fitted': v['phase']['tau_used_to_choose_phase']=True
    elif kind=='phase_derived': v['phase_calibration_is_source_derivation']=True
    elif kind=='loop_proved': v['comparison']['physical_loop_correction_bounds_supplied']=True
    elif kind=='unsorted': v['comparison']['source_action_coefficients'].reverse()
    else: v['comparison']['source_action_coefficients'][0]=['0','0']
    with pytest.raises(ValueError): c.summarize(*args)


@pytest.mark.parametrize('kind',['missing_branch','duplicate_branch','selection','balanced','one_inverse'])
def test_capacity_source_and_balance_forgeries_rejected(inputs,kind):
    args=deepcopy(inputs);v=args[2]
    if kind=='missing_branch': v['branches'].pop()
    elif kind=='duplicate_branch': v['branches'][-1]=deepcopy(v['branches'][0])
    elif kind=='selection': v['primitive_finite_assumptions_force_this_capacity_law']=True
    elif kind=='balanced': v['branches'][0]['source']['positive_root_Koide_ratio']=['2/3','2/3']
    else: v['branches'][0]['admitted_complete_memory_inverse_branches']=1
    with pytest.raises(ValueError): c.summarize(*args)


def test_defect_species_assignment_is_rejected(inputs):
    args=deepcopy(inputs);args[3]['entries'][0]['identified_Standard_Model_species']=['electron']
    with pytest.raises(ValueError): c.summarize(*args)


def test_public_projection_is_json_safe_and_detached(inputs):
    before=json.dumps(inputs,sort_keys=True)
    one=c.summarize(*inputs)
    one['rows'][0]['internal_components']['charged_lepton_left']['source_masks'][0]=-1
    one['charged_source_branches'][1]['variants'][0]['ascending_action_ratios'][0][0]='0'
    fresh=c.summarize(*inputs)
    assert fresh['rows'][0]['internal_components']['charged_lepton_left']['source_masks'][0]>=0
    assert fresh['charged_source_branches'][1]['variants'][0]['ascending_action_ratios'][0]==['1','1']
    assert json.dumps(inputs,sort_keys=True)==before
    json.dumps(fresh,allow_nan=False)
