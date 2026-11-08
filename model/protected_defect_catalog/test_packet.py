"""Admission boundaries and live source classification checks."""
import copy
import subprocess
import sys

import pytest

from . import packet, pauli


@pytest.fixture(scope='module')
def produced():
    return packet.build()


def test_retained_catalogue_replays(produced):
    retained = packet.strict_json.load(packet.HERE/'catalog.json')
    assert packet.canonical(retained) == packet.canonical(produced)
    assert (packet.HERE/'CATALOG.md').read_text() == packet.render(produced)
    for law in ('C', 'E'):
        assert produced['bare_complete_branch_classification'][law]['full78']['commutant_complex_dimension'] == 1
        assert produced['bare_complete_branch_classification'][law]['binary66']['commutant_complex_dimension'] == 1
        assert produced['ordered_cover_complete_branch_classification'][law]['commutant_complex_dimension'] == 2


def test_source_matrix_mutation_rejected(monkeypatch):
    original = pauli.source.matrix
    def changed():
        h, scale = original()
        h = h.copy(); h[0, 0] += 1
        return h, scale
    monkeypatch.setattr(pauli.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='actual source Hamiltonian'):
        pauli.source_polynomial()


@pytest.mark.parametrize('forgery', ['physical_mass','global_completeness','species','Pauli_count'])
def test_claim_promotion_rejected(produced, monkeypatch, forgery):
    # Independent immutable reference avoids another expensive full source
    # replay for each pure admission mutation.
    expected = copy.deepcopy(produced)
    monkeypatch.setattr(packet, 'build', lambda: expected)
    changed = copy.deepcopy(produced)
    if forgery == 'physical_mass':
        changed['entries'][2]['physical_mass'] = 1
    elif forgery == 'global_completeness':
        changed['catalogue_is_exhaustive_over_all_OPH_sources_or_defects'] = True
    elif forgery == 'species':
        changed['entries'][2]['identified_Standard_Model_species'] = ['electron']
    else:
        changed['selected_spin_complete_Pauli_classification']['selected_field_charge_words_modulo_phase'].append([0, 1])
    with pytest.raises(ValueError, match='catalogue classification'):
        packet.validate(changed)


def test_no_double_counting_or_physical_assignments(produced):
    entries = produced['entries']
    assert len({row['id'] for row in entries}) == len(entries)
    assert all(not row['identified_Standard_Model_species'] and row['physical_mass'] is None for row in entries)
    assert not produced['physical_particle_map']
    assert len(next(row for row in entries if row['id']=='ordered_exchange_charge')['same_object_realizations']) == 2
    assert len(entries)==10


def test_hole_sector_keeps_position_spin_scope_and_complete_archive(produced):
    entry=next(row for row in produced['entries']if row['id']=='one_hole_spin_readout')
    attachment=produced['one_hole_source_attachment']
    assert attachment['catalogue_entry']==entry['id']
    assert attachment['owning_hole_scientific_producer_replayed']
    assert not attachment['full_owning_program_countercontrol_replayed']
    assert attachment['complete_corrected_position_instrument']=='sum_x P_x rho P_x'
    assert attachment['arbitrary_reference_preserved_in_each_internal_spin_block']
    assert not attachment['unconditional_hole_coherence_preserved_by_complete_scan']
    assert [r['source_one_hole_dimension']for r in attachment['exact_source_cases']]==[4,12,32]
    assert attachment['counted_cube_trace_bound']=='min(2,4/n); side n^2+1 and m<=n^2'
    assert attachment['complete_source_archive_controls'][-1]['complete_source_archive_trace_bound']=='1/16'
    assert attachment['actual_Hubbard_propagation_scope']['general_one_hole_leakage']['T_over_c_half_matrix_element']!=0
    assert not attachment['topological_invariant_assigned']
    assert not attachment['all_source_states_classified']
    assert not attachment['physical_species'] and not attachment['physical_masses']


@pytest.mark.parametrize('forgery',['reference_lost','ground','counting_premises','clock','erased_record','executed_correction','aligned_ground','receipt'])
def test_hole_owner_premise_or_scope_mutation_rejected(monkeypatch,forgery):
    original=packet.strict_json.load(packet.ROOT/'model/source_family_compatibility/results.json')
    changed=copy.deepcopy(original)
    value=changed['actual_complete_one_hole_repair']
    if forgery=='reference_lost':value['arbitrary_reference_preserved_in_each_internal_spin_block']=False
    elif forgery=='ground':value['one_hole_ground_or_relativistic_particle_assigned']=True
    elif forgery=='counting_premises':value['complete_word_simplex_and_counting_reference_premises_discharged']=True
    elif forgery=='clock':value['global_fair_attempt_complete_record_workaround'][0]['physical_clock_or_capacity_trajectory_selected']=True
    elif forgery=='erased_record':value['record_erasure_restores_hole_coherence']=True
    elif forgery=='executed_correction':value['correction_native_pulse_tape_executed']=True
    elif forgery=='aligned_ground':value['actual_full_Hubbard_propagation_scope']['aligned_sector_identified_as_interacting_ground']=True
    else:value['cases'][0]['source_one_hole_dimension']+=1
    monkeypatch.setattr(packet.hole_records,'certificate',lambda:value)
    monkeypatch.setattr(packet.strict_json,'load',lambda path:original if forgery=='receipt' else changed)
    with pytest.raises(ValueError,match='one-hole'):
        packet.hole_source_attachment()


def test_joint_family_has_literal_observables_and_distinct_state_clock(produced):
    entry=next(row for row in produced['entries'] if row['id']=='joint_positive_field_canonical_sector')
    attached=produced['joint_positive_field_source_attachment']
    assert attached['catalogue_entry']==entry['id']
    assert attached['canonical_pairs']==[['Q3','P2'],['Q2','P3']]
    assert attached['source_observables']['P2']=='-p sum G(px) S_x^2'
    assert attached['source_observables']['P3']=='+p sum G(px) S_x^3'
    assert not attached['ordered_zero_field_history_identity_proved']
    assert not attached['repair_identity_generates_the_H_B_clock']
    assert attached['owning_scientific_producer_replayed']
    assert attached['complete_finite_history_example']['all_native_failure_and_pointer_records_retained']
    assert not attached['physical_species'] and not attached['physical_masses']
    uniform=attached['uniform_mode_response']
    assert uniform['weighted_centre_and_variance_tend_to_zero']
    assert not uniform['uniform_original_continuum_residue_lower_proved']
    assert not uniform['nonzero_momentum_shell_proved']


def test_centered_charge_response_has_original_weight_and_local_record_scope(produced):
    attached=produced['centered_full_CAR_source_attachment']
    assert attached['catalogue_entry']=='centered_full_CAR_charge_response'
    assert attached['owning_local_and_charged_producers_replayed']
    assert attached['original_local_addition_and_removal_weights']==['1/2','1/2']
    assert attached['readable_source_charged_band']['each_original_band_probability_lower']=='25/100352'
    assert not attached['global_raw_archive_equality']
    assert not attached['global_sharp_energy_PVM_implemented_by_local_bound']
    assert not attached['odd_transfer_native_service_derived']
    assert not attached['ordered_spin_state_or_aligned_hole_reference_identity']
    assert not attached['physical_species'] and not attached['physical_masses']
    aligned=produced['aligned_hole_source_attachment']
    assert aligned['catalogue_entry']=='one_hole_spin_readout'
    assert aligned['all_hole_numbers_have_exact_exterior_power_dynamics']
    assert not aligned['interacting_scattering_or_relativistic_mass_proved']


@pytest.mark.parametrize('forgery',['global_archive','ground_preparation','physical_mass','original_weight','original_weyl','bare_frame','receipt'])
def test_centered_owner_scope_mutations_rejected(monkeypatch,forgery):
    original=packet.strict_json.load(packet.ROOT/'model/source_family_compatibility/results.json')
    changed=copy.deepcopy(original)
    local=changed['local_interacting_matter_history'];bands=changed['centered_original_charged_bands']
    if forgery=='global_archive':local['global_raw_archive_equality']=True
    elif forgery=='ground_preparation':local['effective_ground_subsequence_preparation_modulus']=True
    elif forgery=='physical_mass':bands['charged_pole_mass_or_Lorentz_covariance_derived']=True
    elif forgery=='original_weight':bands['cubic_readable_parameter_band']['each_original_preparation_weight']='1'
    elif forgery=='original_weyl':local['original_Weyl_bare_repair_local_history_instantiated']=False
    elif forgery=='bare_frame':local['bare_repairs_identified_under_site_dependent_frame']=True
    else:bands['cubic_readable_parameter_band']['source_U']='2'
    monkeypatch.setattr(packet.local_matter,'certificate',lambda:local)
    monkeypatch.setattr(packet.charge_bands,'certificate',lambda:bands)
    monkeypatch.setattr(packet.strict_json,'load',lambda path:original if forgery=='receipt'else changed)
    with pytest.raises(ValueError,match='centered-matter'):
        packet.centered_matter_attachment()


@pytest.mark.parametrize('forgery',['full_ground','selected_clock','particle'])
def test_aligned_reference_cannot_be_promoted_to_full_ground(monkeypatch,forgery):
    parent=copy.deepcopy(packet.strict_json.load(packet.ROOT/'model/source_family_compatibility/results.json'))
    value=parent['aligned_hole_position_dispersion_and_continuum']
    if forgery=='full_ground':value['full_interacting_ground_reference_selected']=True
    elif forgery=='selected_clock':value['source_capacity_or_time_trajectory_selected']=True
    else:value['physical_particle_labels']=['electron']
    monkeypatch.setattr(packet.hole_dispersion,'certificate',lambda:value)
    monkeypatch.setattr(packet.strict_json,'load',lambda path:parent)
    with pytest.raises(ValueError,match='aligned-hole'):
        packet.aligned_hole_attachment()


@pytest.mark.parametrize('forgery',['missing_second_pair','old_limit_identity','physical_species','unread_selection'])
def test_joint_owner_premise_or_scope_mutation_rejected(monkeypatch,forgery):
    original=packet.strict_json.load(packet.joint_field.OUTPUT)
    changed=copy.deepcopy(original)
    key={'missing_second_pair':'second_canonical_pair_on_same_joint_history_subsequence',
         'old_limit_identity':'equality_to_ordered_zero_field_state_clock_limit_proved',
         'physical_species':'physical_particle_assignment',
         'unread_selection':'sequence_forced_by_unread_repair'}[forgery]
    changed['scope'][key]=False if forgery=='missing_second_pair' else ('electron' if forgery=='physical_species' else True)
    # The owner is a separate scientific admission; neither changed output nor
    # matching retained bytes can transfer its unproved physical premises.
    monkeypatch.setattr(packet.joint_field,'build',lambda:changed)
    monkeypatch.setattr(packet.strict_json,'load',lambda path:changed)
    with pytest.raises(ValueError,match='joint-field'):
        packet.joint_field_attachment()


def test_fresh_full_producer_binds_runtime_imports():
    script = '''
from pathlib import Path
import sys
from model.protected_defect_catalog import packet
result = packet.build()
missing = []
for module in tuple(sys.modules.values()):
    filename = getattr(module, '__file__', None)
    if filename:
        path = Path(filename).resolve()
        if path.is_relative_to(packet.ROOT) and path.suffix == '.py':
            name = path.relative_to(packet.ROOT).as_posix()
            if name not in result['source_bindings']:
                missing.append(name)
assert not missing, sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', script], cwd=packet.ROOT,
                   capture_output=True, text=True, check=True)
