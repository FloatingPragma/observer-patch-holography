"""Candidate IDs, exact internal weights and guarded physical interpretation."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from . import model as m


@pytest.fixture(scope='module')
def value():return m.build()


def test_exact_fifteen_weights_and_named_family_reuse(value):
    rows={r['id']:r for r in value['source_excitation_rows']}
    internal=[r for r in rows.values() if r['kind']=='native_internal_component']
    masks=[x for r in internal for x in r['exact']['source_masks']]
    assert len(masks)==len(set(masks))==15
    for row in internal:
        e=row['exact'];assert F(e['comparison_Q'])==F(e['T3'])+F(e['Q6'])/6
    targets={r['target']:r for r in value['target_rows']}
    for group in (('electron','muon','tau'),('up','charm','top'),('down','strange','bottom'),
                  ('electron_neutrino','muon_neutrino','tau_neutrino')):
        assert len({tuple(targets[t]['source_component_ids']) for t in group})==1
        assert all(targets[t]['source_flavor_index'] is None for t in group)
    assert rows['charged_lepton_right_conjugate']['exact']['comparison_Q']=='1'
    assert rows['up_right_conjugate']['exact']['SU3']=='3bar'
    weights=value['internal_gauge_multiplets']
    assert [(r['SU3'],r['SU2_dimension'],r['Y']) for r in weights]==[
        ('3',2,'1/6'),('3bar',1,'-2/3'),('3bar',1,'1/3'),('1',2,'-1/2'),('1',1,'1')]
    for r in weights:
        assert all(F(q)==F(t)+F(r['Y']) for q,t in zip(r['component_Q'],r['component_T3'],strict=True))


def test_regional_charged_factor_is_not_the_global_pair_or_neutral_trion(value):
    rows={r['id']:r for r in value['source_excitation_rows']}
    e=rows[m.REGIONAL_ID]['exact'];assert (e['Q6'],e['source_cube_spin'])==(6,'G1')
    assert (e['regional_CAR_parity'],e['complete_pair_CAR_parity'])==(1,0)
    assert rows[m.PAIR_ID]['exact']['internal_probe_irreps']=={'A1':1,'T1':3}
    trions=[r for r in rows.values() if r['kind']=='actual_local_gauge_singlet_composite']
    assert len(trions)==6
    assert all(r['exact']['Q6']==0 and r['exact']['SU2_dimension']==1 for r in trions)
    assert sum(r['exact']['CAR_parity']==1 for r in trions)==4
    assert any('antibaryon' in r['best_guess'] for r in trions)
    assert all(not r['possible_named_targets'] for r in trions)
    line=rows[m.LOOP_ID]
    lo,hi=map(F,line['exact']['transition_interval_source_units'])
    assert lo<63<hi and not line['possible_named_targets']


def test_protected_catalogue_domains_remain_separate(value):
    original=json.loads((m.ROOT/m.INPUTS[0]).read_text())
    assert [r['id'] for r in value['protected_catalogue_rows']]==[r['id'] for r in original['entries']]
    assert all(r['protected_defect_id'] is None for r in value['source_excitation_rows'])
    exclusions={r['target']:r for r in value['rejection_controls'] if r['status']=='excluded'}
    assert 'direct' in exclusions['active_neutrino']['domain']
    assert 'different' not in exclusions['quark']['domain']
    assert value['physical_particle_map']==value['mass_predictions_GeV']==[]


@pytest.mark.parametrize('mutation',['mass','flavor','Lorentz_spin','topological_ID','missing_target','collision','ratio'])
def test_unproved_promotions_are_rejected(value,mutation):
    bad=deepcopy(value)
    if mutation=='mass':bad['target_rows'][0]['source_derived_mass_GeV']=.000511
    elif mutation=='flavor':bad['target_rows'][0]['source_flavor_index']=0
    elif mutation=='Lorentz_spin':bad['target_rows'][0]['physical_Lorentz_spin_proved']=True
    elif mutation=='topological_ID':bad['source_excitation_rows'][0]['protected_defect_id']='ordered_exchange_charge'
    elif mutation=='missing_target':bad['target_rows'][0]['target']='invented_particle'
    elif mutation=='collision':bad['collision_status']['moving_two_packet_collision_executed']=True
    else:bad['finite_sector_ratio_control']['result']['physical_mass_table_GeV']=[1]
    with pytest.raises(ValueError):m.validate(bad)


def test_full_menu_fusion_and_exact_ratio_domains(value):
    c=value['collision_status'];assert c['adjacent_fused_cluster_rank']==44
    assert c['each_complete_pair_row_original_weight_lower']['rational_multiplier']=='1/40'
    assert c['each_complete_pair_row_original_weight_lower']['negative_exponent']==32406376
    assert F(c['adjacent_fused_energy_interval'][1])<F(c['opposite_disjoint_energy_window'][0])
    separate=c['separate_contact_extension']
    assert separate['finite_interaction_echo_executed'] and separate['native_pulses']==22531
    echo=separate['executed_interacting_echo'];free=separate['executed_free_echo']
    assert echo['source_resources']['native_pulses']==22531 and echo['all_outputs_retained']
    assert echo['pointer_probability_interval']==['469/10000','47/1000']
    assert F(free['pointer_probability_upper'])==F(1,10**24)
    assert F(separate['exact_ideal_word_signal_certificate']['ideal_word_interaction_signal_lower'])==F(579197,16000000)
    assert separate['small_32_row_native_control_executed']
    assert all(r['joint_quantum_rows']==32 for r in separate['small_native_controls'])
    assert not separate['fine_spacing_packet_experiment_executed']
    r=value['finite_sector_ratio_control']['result']
    assert F(r['all_three_family_ratio_upper'])<F('21.16')<F(r['observational_comparison']['central_value'])
    assert r['direct_common_scale_identification_with_charged_lepton_hierarchy_excluded_at_member_zero']
    assert not r['exclusion_applies_to_other_capacities_or_other_source_excitations']


def test_nested_candidate_package_does_not_change_frozen_particle_map_closure(value):
    old=json.loads((m.ROOT/m.INPUTS[1]).read_text())
    assert not any('/candidates/' in p for p in old['source_bindings'])
    for path,expected in old['source_bindings'].items():
        assert m.sha(m.ROOT/path)==expected
    for path,expected in value['source_sha256'].items():
        assert m.sha(m.ROOT/path)==expected
    assert 'model/standard_model_source_map/candidates/finite_sector_ratios.py' in value['source_sha256']
