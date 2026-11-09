"""Consumer scope guards reject matching forged scientific owner payloads.

The full owner replay is tested separately. These tests isolate the consumer
boundary with genuine retained payloads and explicitly matched producer data,
so a failed test cannot be hidden by a receipt mismatch or an expensive solve.
"""
import copy

import pytest

from model.protected_defect_catalog import packet


def _set(payload, path, value):
    for key in path[:-1]:
        payload = payload[key]
    payload[path[-1]] = value


@pytest.fixture(scope='module')
def kinetic_payload():
    return packet.serialize(packet.kinetic_nodes.certificate())


def test_kinetic_consumer_accepts_matching_scoped_source(kinetic_payload, monkeypatch):
    value = copy.deepcopy(kinetic_payload)
    monkeypatch.setattr(packet.kinetic_nodes, 'certificate', lambda: value)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: {'kinetic_node_to_free_field_map': value})
    assert packet.kinetic_field_attachment() == value


@pytest.mark.parametrize('path,value', [
    (('scope', 'default_interacting_Green_function_winding_derived'), True),
    (('scope', 'physical_Standard_Model_assignment'), True),
    (('scope', 'all_native_repair_protected_charge'), True),
    (('scope', 'real_space_defect_classification'), True),
    (('free_field_scope', 'default_U_over_c_eight_in_this_class'), True),
    (('nodes', 0, 'physical_particle_assignment'), 'electron'),
])
def test_kinetic_scope_forgery_rejected_despite_matching_owner(kinetic_payload, monkeypatch, path, value):
    changed = copy.deepcopy(kinetic_payload)
    _set(changed, path, value)
    monkeypatch.setattr(packet.kinetic_nodes, 'certificate', lambda: changed)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: {'kinetic_node_to_free_field_map': changed})
    with pytest.raises(ValueError, match='kinetic-field map exceeds its physical scope'):
        packet.kinetic_field_attachment()


def test_kinetic_replay_mismatch_is_rejected_separately(kinetic_payload, monkeypatch):
    changed = copy.deepcopy(kinetic_payload)
    changed['nodes'][0]['kinetic_winding'] *= -1
    monkeypatch.setattr(packet.kinetic_nodes, 'certificate', lambda: changed)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: {'kinetic_node_to_free_field_map': kinetic_payload})
    with pytest.raises(ValueError, match='kinetic-field owning source receipt differs'):
        packet.kinetic_field_attachment()


@pytest.fixture(scope='module')
def joint_payload():
    parent = packet.strict_json.load(packet.joint_field.OUTPUT)
    if 'bounded_joint_charged_energy_momentum_bands' not in parent:
        parent['bounded_joint_charged_energy_momentum_bands'] = packet.serialize(
            packet.joint_field.arctan_bands.certificate())
    return parent


def test_joint_consumer_accepts_complete_scoped_bands(joint_payload, monkeypatch):
    value = copy.deepcopy(joint_payload)
    monkeypatch.setattr(packet.joint_field, 'build', lambda: value)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: value)
    result = packet.joint_field_attachment()
    assert result['bounded_charged_joint_bands']['actual_finite_source']['complete_record_count'] == 9
    assert result['physical_species'] == [] and result['physical_masses'] == []


@pytest.mark.parametrize('path,value', [
    (('source_quantum_correction_to_mean_spin_force','liminf_kappa_minus_m_squared_lower'), '1'),
    (('source_quantum_correction_to_mean_spin_force','corrected_kinetic_coefficients_excluded'), True),
    (('source_quantum_correction_to_mean_spin_force','infrared_sharp_component_excluded'), True),
    (('source_quantum_correction_to_mean_spin_force','source_to_history_unbounded_force_norm_preservation_assumed'), True),
    (('source_quantum_correction_to_mean_spin_force','physical_mass_or_species_assigned'), True),
    (('source_quantum_correction_to_mean_spin_force','physical_mass_or_species_assigned'), 0),
    (('source_quantum_correction_to_mean_spin_force','additional_native_detector_execution'), True),
    (('scope','strict_source_bond_minus_squared_order_margin_proved'), False),
    (('scope','unfiltered_microscopic_mean_spin_force_factorization_excluded'), False),
    (('scope','quantum_force_margin_excludes_corrected_infrared_shell'), True),
])
def test_quantum_force_forgery_rejected_despite_matching_owner(joint_payload, monkeypatch, path, value):
    changed = copy.deepcopy(joint_payload)
    _set(changed, path, value)
    monkeypatch.setattr(packet.joint_field, 'build', lambda: changed)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: changed)
    with pytest.raises(ValueError, match='joint-field quantum force correction'):
        packet.joint_field_attachment()


@pytest.mark.parametrize('path,value', [
    (('continuum_band', 'each_original_charge_probability_lower'), '1'),
    (('continuum_band', 'joint_positive_field_clock_used'), False),
    (('continuum_band', 'sharp_shell_or_physical_mass_inferred'), True),
    (('continuum_band', 'old_zero_field_clock_substituted'), True),
    (('continuum_band', 'original_odd_CAR_identification'), True),
    (('full_history_linear_force_closure_proved',), True),
    (('continuum_limit_error_rate_supplied',), True),
    (('physical_species',), ['electron']),
    (('physical_mass_prediction',), '1'),
])
def test_joint_band_forgery_rejected_despite_matching_owner(joint_payload, monkeypatch, path, value):
    changed = copy.deepcopy(joint_payload)
    _set(changed['bounded_joint_charged_energy_momentum_bands'], path, value)
    monkeypatch.setattr(packet.joint_field, 'build', lambda: changed)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: changed)
    with pytest.raises(ValueError, match='joint-field charged response exceeds its owning scope'):
        packet.joint_field_attachment()


_CENTERED_PRODUCERS = {
    'local_matter': 'local_interacting_matter_history',
    'charge_bands': 'centered_original_charged_bands',
    'composite': 'source_generated_composite_CAR_response',
    'cube_response': 'default_cube_composite_response',
    'charge_scaling': 'volume_uniform_charged_scaling',
    'boost_current': 'full_source_energy_dipole',
    'cube_spectrum': 'complete_cube_wavepacket_spectrum',
    'momentum_response': 'actual_cubic_momentum_response',
    'charge_escape': 'default_original_charge_escape',
    'uv_tail': 'uniform_original_energy_tails',
    'local_spectral': 'complete_local_source_energy_readout',
    'nodal_memory': 'actual_nodal_response_memory',
    'ground_marginals': 'ground_marginal_preparation_and_native_records',
    'cube_memory': 'full_cube_frequency_dependent_memory',
    'nodal_channels': 'actual_nodal_scalar_channels',
    'spectral_residue': 'actual_spectral_residue_and_subcontact_band',
    'spectral_window_access': 'complete_Lorentzian_window_access',
    'nodal_locality': 'actual_nodal_locality_and_window_survival',
    'nodal_moment_bounds': 'uniform_full_CAR_star_force_bound',
    'mass_scale_selection': 'fixed_source_scaling_and_mass_shell_scope',
    'source_action_identification': 'complete_source_action_identification',
    'nodal_next_moments': 'actual_local_nodal_higher_moments',
    'native_lorentzian': 'executed_native_dimer_Lorentzian_readout',
    'bulk_action_identification': 'volume_uniform_cubic_action_identification',
}


@pytest.fixture(scope='module')
def centered_payload():
    parent = packet.strict_json.load(packet.ROOT/'model/source_family_compatibility/results.json')
    # The official owner write precedes catalogue admission. Its complete
    # scientific replay is not duplicated by these scope-only mutations.
    assert all(key in parent for key in _CENTERED_PRODUCERS.values())
    return parent


def _stub_centered_producers(monkeypatch, parent):
    for name, key in _CENTERED_PRODUCERS.items():
        monkeypatch.setattr(getattr(packet, name), 'certificate', lambda key=key: copy.deepcopy(parent[key]))
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: copy.deepcopy(parent))


def test_centered_consumer_accepts_original_records_without_thermodynamic_promotion(centered_payload, monkeypatch):
    parent = copy.deepcopy(centered_payload)
    _stub_centered_producers(monkeypatch, parent)
    result = packet.centered_matter_attachment()
    assert result['physical_species'] == [] and result['physical_masses'] == []
    assert result['ground_marginal_preparation_and_native_records']['cut_ground_substituted_for_infinite_ground'] is False


@pytest.mark.parametrize('key,path,value,error', [
    ('volume_uniform_cubic_action_identification', ('known_local_calibration_preparation_consumed',), False, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('interacting_ground_preparation_inferred',), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('complete_history_or_reference_error_from_marginal_inferred',), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('capacity_selects_action_or_clock',), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('large_native_experiment_executed',), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('physical_particles',), ['electron'], 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('physical_masses',), ['1'], 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('finite_time_ratio_witnesses',0,'is_probability_of_the_larger_collar_plan'), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('separate_arbitrary_exterior_collar_plan','pure_two_CAR_collar_input_required'), False, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('separate_arbitrary_exterior_collar_plan','reset_pulse_tape_executed'), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('complete_collision_detector','full_fine_instrument_covariance_asserted'), True, 'local calibration scope'),
    ('volume_uniform_cubic_action_identification', ('uniform_global_attempt_example','other_seam_records_discarded_or_postselected'), True, 'local calibration scope'),
    ('actual_local_nodal_higher_moments', ('source_selected_ground_moments_evaluated',), True, 'selected-ground scope'),
    ('actual_local_nodal_higher_moments', ('actual_zero_atom_decided',), True, 'selected-ground scope'),
    ('actual_local_nodal_higher_moments', ('physical_species_identified',), ['electron'], 'selected-ground scope'),
    ('actual_local_nodal_higher_moments', ('physical_masses',), ['1'], 'selected-ground scope'),
    ('actual_local_nodal_higher_moments', ('actual_source_memory_bound','zero_frequency_atom_inferred'), True, 'selected-ground scope'),
    ('actual_local_nodal_higher_moments', ('actual_local_nodal_moment_readouts',0,'selected_ground_moment'), '328', 'selected-ground scope'),
    ('executed_native_dimer_Lorentzian_readout', ('large_infinite_ground_collar_executed',), True, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('physical_species_or_mass_assigned',), True, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('hardware_shots_sampled',), True, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('postselection_or_output_renormalization_used',), True, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('floating_contraction_diagnostics_used_as_interval_premise',), True, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('dummy_clock_environment_retained',), False, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('raw_pointer_terminal_record_count',), 8194, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('physical_first_failure_records',), 0, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('zero_charge_prefixes_follow_same_declared_controller',), False, 'complete finite-source scope'),
    ('executed_native_dimer_Lorentzian_readout', ('selector_blanks_are_supplied',), False, 'complete finite-source scope'),
    ('complete_source_action_identification', ('capacity_selects_action_or_physical_clock',), True, 'finite observation scope'),
    ('complete_source_action_identification', ('physical_particle_or_mass_identified',), True, 'finite observation scope'),
    ('complete_source_action_identification', ('physical_strong_coupling_or_QCD_alpha_selected',), True, 'finite observation scope'),
    ('complete_source_action_identification', ('thermodynamic_ground_selected_by_this_readout',), True, 'finite observation scope'),
    ('complete_source_action_identification', ('coherent_clock_pulse_tape_executed_here',), True, 'finite observation scope'),
    ('complete_source_action_identification', ('source_clock_and_known_initial_preparation_required',), False, 'finite observation scope'),
    ('complete_source_action_identification', ('source', 'postselection_used'), True, 'finite observation scope'),
    ('fixed_source_scaling_and_mass_shell_scope', ('arbitrary_sparse_subsequence_excluded_by_theorem',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('geometric_limit_has_continuous_dilation_covariance',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('every_unstable_resonance_excluded',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('continuous_timelike_response_excluded',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('actual_full_continuum_limit_proved',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('actual_nodal_atom_positivity_proved',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('running_reference_or_observable_dictionary_excluded',), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('fixed_reference_and_dimensionless_couplings_required',), False, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('full_joint_vague_limit_required',), False, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('physical_species_assigned',), 1, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('physical_mass_assigned',), 'electron', 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('actual_default_source_homogeneity', 'native_large_detector_tape_executed'), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('actual_default_source_homogeneity', 'finite_dimer_is_continuum_momentum_mode'), True, 'dilation exceeds'),
    ('fixed_source_scaling_and_mass_shell_scope', ('actual_default_source_homogeneity', 'conditional_probability_renormalization_applied'), True, 'dilation exceeds'),
    ('actual_nodal_response_memory', ('physical_species_identified',), ['electron'], 'energy readout promotes'),
    ('actual_nodal_response_memory', ('default_source_result', 'zero_frequency_residue_inferred'), True, 'energy readout promotes'),
    ('ground_marginal_preparation_and_native_records', ('thermodynamic_local_marginal_convergence_modulus_proved',), True, 'finite preparation promotes'),
    ('ground_marginal_preparation_and_native_records', ('cut_ground_substituted_for_infinite_ground',), True, 'finite preparation promotes'),
    ('ground_marginal_preparation_and_native_records', ('physical_state_selection_proved',), True, 'finite preparation promotes'),
    ('full_cube_frequency_dependent_memory', ('physical_mass_assigned',), 'electron', 'cube memory exceeds'),
    ('full_cube_frequency_dependent_memory', ('thermodynamic_charge_gap_proved',), True, 'cube memory exceeds'),
    ('full_cube_frequency_dependent_memory', ('native_large_filter_tape_executed',), True, 'cube memory exceeds'),
    ('actual_nodal_scalar_channels', ('uniform_growing_packet_scaling_supplied',), True, 'nodal channels exceed'),
    ('actual_nodal_scalar_channels', ('physical_spin_or_species_identification',), True, 'nodal channels exceed'),
    ('actual_nodal_scalar_channels', ('physical_mass_prediction',), True, 'nodal channels exceed'),
    ('actual_spectral_residue_and_subcontact_band', ('physical_mass_prediction',), '1', 'spectral residue exceeds'),
    ('actual_spectral_residue_and_subcontact_band', ('physical_species',), ['electron'], 'spectral residue exceeds'),
    ('actual_spectral_residue_and_subcontact_band', ('default_shrinking_window_weight_proved',), True, 'spectral residue exceeds'),
    ('actual_spectral_residue_and_subcontact_band', ('uniform_many_body_charge_gap_proved',), True, 'spectral residue exceeds'),
    ('actual_nodal_locality_and_window_survival', ('actual_nodal_atom_positivity_proved',), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('physical_mass_assigned',), 'electron', 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('sharp_dispersion_proved',), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('growing_native_instrument_executed',), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'actual_nodal_atom_lower_proved'), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'native_growing_readout_error_met'), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'positive_energy_annulus_proved'), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'dispersion_or_particle_identification_proved'), True, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'same_normalized_reference_at_every_spacing_required'), False, 'locality promotes'),
    ('actual_nodal_locality_and_window_survival', ('conditional_example', 'zero_energy_included'), False, 'locality promotes'),
    ('uniform_full_CAR_star_force_bound', ('physical_species_identified',), ['electron'], 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('physical_mass_prediction',), '1', 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('spin_only_ground_or_antiferromagnetic_order_assumed',), True, 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('universal_local_bound_attained_in_the_actual_infinite_ground',), True, 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('native_large_tape_executed',), True, 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('source_moment_budget', 'physical_species_identified'), ['electron'], 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('source_moment_budget', 'physical_mass_prediction'), '1', 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('source_moment_budget', 'actual_zero_atom_positive_proved'), True, 'star bounds promote'),
    ('uniform_full_CAR_star_force_bound', ('source_moment_budget', 'source_eighth_moment_computed'), True, 'star bounds promote'),
    ('complete_Lorentzian_window_access', ('physical_species_assigned',), True, 'spectral access exceeds'),
    ('complete_Lorentzian_window_access', ('physical_mass_assigned',), '1', 'spectral access exceeds'),
    ('complete_Lorentzian_window_access', ('uniform_shrinking_window_weight_proved',), True, 'spectral access exceeds'),
    ('complete_Lorentzian_window_access', ('effective_growing_ground_marginal_modulus_proved',), True, 'spectral access exceeds'),
])
def test_centered_memory_and_state_scope_forgery_rejected(centered_payload, monkeypatch, key, path, value, error):
    changed = copy.deepcopy(centered_payload)
    _set(changed[key], path, value)
    _stub_centered_producers(monkeypatch, changed)
    with pytest.raises(ValueError, match=error):
        packet.centered_matter_attachment()


@pytest.mark.parametrize('path,value', [
    (('directed_joint_source_band', 'each_charge_and_direction_original_probability_lower'), '1'),
    (('directed_joint_source_band', 'source_even_collective_internal_charge'), False),
    (('directed_joint_source_band', 'original_odd_CAR_identification'), True),
    (('directed_joint_source_band', 'energy_shell_or_plane_wave_atom_inferred'), True),
    (('directed_joint_source_band', 'effective_finite_source_to_continuum_error_claimed'), True),
    (('finite_periodic_source_control', 'complete_native_pulse_tape_executed'), True),
    (('finite_periodic_source_control', 'finite_control_proves_continuum_floor'), True),
    (('finite_periodic_source_control', 'finite_translation_labels_identified_with_continuum_momentum'), True),
    (('finite_periodic_source_control', 'physical_particle_or_mass_assignment'), True),
])
def test_directed_scope_forgery_rejected_despite_matching_owner(joint_payload, monkeypatch, path, value):
    changed = copy.deepcopy(joint_payload)
    _set(changed['directed_joint_charged_bands'], path, value)
    monkeypatch.setattr(packet.joint_field, 'build', lambda: changed)
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: changed)
    with pytest.raises(ValueError, match='directed bands exceed'):
        packet.joint_field_attachment()
