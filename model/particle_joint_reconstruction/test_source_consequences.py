"""Keep finite source consequences distinct from physical particle claims."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import pytest

from model.particle_common_contract import contract
from tools import strict_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GROUPS = (
    ('primary_extreme_gauge_isolated_transition', 'particle_primary_reduced_sector'),
    ('primary_regional_charged_spin_and_pair_energy_records',
     'particle_primary_regional_fermion'),
    ('primary_selected_weighted_spectral_tail', 'particle_primary_weighted_spectral_tail'),
    ('primary_original_energy_current_response', 'particle_primary_energy_response'),
)


@pytest.fixture(scope='module')
def aggregate():
    # These checks consume the completed replay; they never run a producer.
    return strict_json.load(HERE / 'results.json')


def members(certificate):
    if 'selected_member_bound' in certificate:
        return [certificate['selected_member_bound']]
    return certificate['members']


@pytest.mark.parametrize('group,package', GROUPS)
def test_aggregate_keeps_exact_receipt_and_same_member_zero_contract(aggregate, group, package):
    certificate = aggregate[group]
    assert certificate == strict_json.load(ROOT / 'model' / package / 'results.json')
    first = members(certificate)[0]['result']
    for key in ('member_id', 'source_id', 'state_id', 'reference_id', 'clock_id'):
        assert first[key] == aggregate[key]
    for bound in members(certificate):
        assert bound['result_sha256'] == contract.digest(bound['result'])
        assert bound['result']['physical_mass_table_GeV'] == []
    own_prefix = f'model/{package}/'
    own_pins = {p: h for p, h in certificate['source_sha256'].items()
                if p.startswith(own_prefix)}
    assert own_prefix + 'PROOF.md' in own_pins
    assert own_prefix + 'model.py' in own_pins
    for path, expected in own_pins.items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == expected
        assert aggregate['input_sha256'][path] == expected


def test_scoped_results_do_not_promote_the_physical_vacuum_branch(aggregate):
    for flag in ('primary_extreme_gauge_isolated_transition_proved',
                 'primary_regional_charged_spin_records_proved',
                 'primary_selected_weighted_spectral_tail_proved',
                 'primary_original_energy_current_response_proved',
                 'primary_static_loop_fusion_energy_records_proved',
                 'primary_same_source_equal_work_control_proved'):
        assert aggregate[flag] is True
    assert 'physical-vacuum particle branch' in aggregate['primary_isolated_excitation_branch_domain']
    for flag in ('primary_isolated_excitation_branch_proved',
                 'primary_uniform_fixed_band_positive_energy_weight_proved',
                 'primary_conserved_particle_momentum_proved',
                 'primary_proper_geometry_or_Einstein_limit_proved',
                 'full_particle_spectrum_derived',
                 'interacting_relativistic_primary_limit_proved',
                 'pure_Yang_Mills_gap_proved'):
        assert aggregate[flag] is False
    assert aggregate['finite_primary_acceptance_count'] == 3
    assert aggregate['source_derived_particle_mass_table_GeV'] == []


def test_isolated_gauge_atom_uses_the_rare_reducing_sector_and_full_complement(aggregate):
    result = aggregate[GROUPS[0][0]]['selected_member_bound']['result']
    spectrum = result['interacting_source_spectrum']
    assert result['source_input_is_original_selected_density']
    assert result['same_source_state_reference_clock']
    assert 'occupations zero' in result['reducing_sector']
    assert spectrum['same_source_interaction_coefficient'] == 1
    assert spectrum['isolated_gauge_vacuum_cluster_rank'] == 1
    assert spectrum['isolated_first_electric_cluster_rank'] == 12
    assert spectrum['isolated_A1u_cluster_rank_within_symmetry_sector'] == 1
    ground = list(map(F, spectrum['interacting_gauge_ground_energy_interval']))
    excited = list(map(F, spectrum['interacting_A1u_energy_interval']))
    gap = list(map(F, spectrum['gauge_ground_to_A1u_energy_gap_interval']))
    assert ground == [-F(1, 45675), F(0)]
    assert excited == [F(63) - F(59, 1351350), F(63)]
    assert gap == [excited[0] - ground[1], excited[1] - ground[0]]
    assert F(spectrum['next_A1u_energy_lower']) - excited[1] == F(429, 14)
    assert F(spectrum['actual_ground_to_odd_matrix_element_squared_strict_lower']) == F(4, 49)
    assert not spectrum['A1u_energy_is_globally_nondegenerate_across_other_cube_irreps']
    assert not spectrum['raw_loop_probe_has_no_other_transition_frequencies_claimed']
    assert result['finite_symmetry_resolved_gauge_transition_atom_proved']
    assert result['extreme_sector_supported_contraction_weight_vanishes_on_selected_refinement']
    for key in ('CAR_odd_trion_branch_isolated', 'ambient_primary_vacuum_gap_proved',
                'refinement_uniform_observable_residue_proved',
                'four_dimensional_relativistic_or_Yang_Mills_continuum_proved'):
        assert result[key] is False


def test_gauge_record_clock_loss_multiplies_the_original_population(aggregate):
    result = aggregate[GROUPS[0][0]]['selected_member_bound']['result']
    floor = result['original_energy_local_population_floor']
    event = result['complete_clocked_original_event']
    loss = F(event['one_clock_loss_upper'])
    original = event['original_weight_lower_factored']
    simple = event['simple_original_weight_strict_lower_factored']
    assert original['negative_exponent'] == floor['negative_exponent'] == 32406376
    assert original['base'] == floor['base'] == 2
    assert F(original['rational_multiplier']) == F(4, 49) * (1-loss)**2
    assert F(original['rational_multiplier']) > F(simple['rational_multiplier']) == F(1, 13)
    assert event['all_source_native_sectors_symmetry_and_failure_outcomes_retained']
    assert event['idle_reference_outputs_retained']
    assert not event['empty_sector_success_divided_out']
    assert not event['joint_reference_Loewner_bound_asserted']
    assert not event['all_other_initial_state_clock_leakage_excluded']


def test_regional_spin_menu_preserves_compensator_and_complete_weights(aggregate):
    result = aggregate[GROUPS[1][0]]['members'][0]['result']
    assert result['regional_endpoint_irrep'] == 'G1'
    assert result['regional_endpoint_dimension'] == 2
    assert result['global_probe_CAR_parity'] == 0
    assert result['inside_and_outside_matter_factor_CAR_parity'] == 1
    assert result['endpoint_charge_increment'] + result['compensating_charge_increment'] == 0
    assert result['complete_probe_outcomes'] == result['complete_success_outcomes'] + 1 == 97
    assert len(result['complete_menu']) == 96
    for row in result['complete_menu']:
        assert sorted(value for value in row['charge_increment'] if value) == [-6, 6]
    assert F(result['complete_menu_each_coefficient_recorded_weight_lower']) * 4 == F(
        result['per_chosen_coefficient_recorded_weight_lower'])
    assert F(result['complete_menu_scalar_recorded_simple_lower']) == F(5, 64)
    assert F(result['complete_menu_vector_recorded_simple_lower']) == F(15, 64)
    raw = result['raw_spin_apparatus']
    assert F(raw['each_of_four_joint_spin_rows_recorded_weight_lower']) * 4 == F(
        raw['per_fixed_spin_entry_recorded_weight_lower'])
    assert F(raw['each_of_four_joint_spin_rows_recorded_weight_lower']) > F(1, 16)
    assert raw['raw_doubled_spin_increments'] == [[1, -1], [1, 1], [-1, -1], [-1, 1]]
    assert raw['charge_and_spin_counters_preserve_original_complete_reference_history']
    assert not raw['initial_absolute_spin_or_charge_measured']
    assert not raw['raw_spin_and_Pauli_records_are_same_joint_apparatus']
    assert not raw['spin_increment_is_final_spin_after_energy_clock']
    for key in ('source_carrier_state_reference_clock_changed',
                'spinor_label_is_physical_Lorentz_spin',
                'pair_energy_divided_by_two_is_justified',
                'pair_energy_is_single_endpoint_rest_energy',
                'isolated_charged_particle_or_physical_electron_identified'):
        assert result[key] is False


@pytest.mark.parametrize('index', range(4))
def test_selected_tail_is_a_capacity_dependent_complete_record_bound(aggregate, index):
    result = aggregate[GROUPS[2][0]]['members'][index]['result']
    a = F(result['local_current_caps']['spacing'])
    cap = F(result['selected_second_moment_upper'])
    assert F(result['family_second_moment_upper']) == 4000000/a
    assert cap <= min(F(result['inherited_operator_second_moment_upper']),
                      F(result['selected_weighted_quadratic_second_moment_upper']),
                      4000000/a)
    ideal = result['ideal_original_positive_window']
    lower, upper = map(F, ideal['original_transfer_interval'])
    weight = F(ideal['original_weight_lower'])
    assert lower > 0 and lower * upper == cap
    assert weight == 4 * lower**2 / cap
    recorded = result['finite_reported_positive_window']
    loss = F(recorded['two_complete_clock_probability_loss_upper'])
    assert F(recorded['original_probability_lower']) == weight - loss > weight/2
    eta = F(recorded['energy_resolution'])
    assert list(map(F, recorded['reported_transfer_interval'])) == [lower-2*eta, upper+2*eta]
    assert lower-2*eta > 0
    assert recorded['all_fine_records_retained'] > 0
    assert not recorded['conditioned_or_postselected_probability']
    assert result['source_state_reference_clock_and_action_unchanged']
    assert result['finite_window_floor_order'] == 'a'
    assert result['window_upper_endpoint_order'] == '1/a'
    for key in ('fixed_positive_energy_probability_limit_proved',
                'uniform_selected_second_moment_proved',
                'uniform_selected_energy_weighted_tail_proved'):
        assert result[key] is False


@pytest.mark.parametrize('index,level', enumerate((0, 1, 2, 8)))
def test_original_current_record_retains_exact_factored_pulse_and_work(aggregate, index, level):
    result = aggregate[GROUPS[3][0]]['members'][index]['result']
    a = F(1, 2**level)
    assert F(result['projected_current_norm_squared_lower']) == F(27, 2)/a**4
    witness = result['physical_projected_current_witness']
    assert witness['actual_modes_in_two_cell_chart'] == 120
    assert not witness['witness_is_the_original_preparation']
    assert not witness['full_commutator_replaced_by_projected_commutator']
    q = result['original_double_commutator_lower']
    t = result['pulse_duration_source_clock']
    d = result['complete_binary_total_variation_lower']
    work = result['cyclic_source_work_lower']
    exponent = result['source_density_floor']['negative_exponent']
    assert q['negative_exponent'] == t['negative_exponent'] == exponent
    assert d['negative_exponent'] == 2*exponent
    assert work['negative_exponent'] == 4*exponent
    C = F(result['current_norm_upper'])
    G = F(result['controlled_H_plus_A_norm_upper'])
    assert F(t['rational_multiplier']) == F(q['rational_multiplier'])/(16*G**2*C)
    assert F(d['rational_multiplier']) == F(q['rational_multiplier'])*F(t['rational_multiplier'])/(4*C)
    beta = F(result['selected_modular_slope']['beta'])
    assert F(work['rational_multiplier']) == 2*F(d['rational_multiplier'])**2/beta > 0
    assert list(map(F, result['numeric_record_labels'])) == [-C, C]
    assert F(result['no_pulse_original_plus_probability']) == F(1, 2)
    assert result['same_source_state_reference_clock']
    assert result['applied_control_work_term_retained']
    assert result['finite_Ward_response_commutes_with_spatial_coarsening']
    assert result['complete_lapse_control_and_current_records_constructed']
    for key in ('large_source_histogram_evaluated',
                'pulse_control_is_a_proved_proper_metric',
                'uniform_physical_refinement_response_proved',
                'Lorentz_stress_or_Einstein_dynamics_proved'):
        assert result[key] is False


def test_static_loop_fusion_keeps_full_complement_and_original_menu_weights(aggregate):
    certificate = aggregate['primary_static_loop_fusion_energy_records']
    assert certificate == strict_json.load(ROOT / 'model/particle_primary_loop_fusion/results.json')
    bound = certificate['selected_member_bound']
    result = bound['result']
    assert bound['result_sha256'] == contract.digest(result)
    for key in ('member_id', 'source_id', 'state_id', 'reference_id', 'clock_id'):
        assert result[key] == aggregate[key]
    assert result['same_source_state_reference_clock'] and result['original_selected_density_retained']
    geometry = result['source_loop_geometry']
    assert geometry['combined_electric_cluster_dimension'] == 44
    assert geometry['q6_hexagon_electric_eigenspace_dimension'] == 32
    assert geometry['weak_square_electric_eigenspace_dimension'] == 12
    assert len(geometry['all_ordered_distinct_face_pairs']) == 30
    assert geometry['adjacent_ordered_pairs'] == 24 and geometry['opposite_ordered_pairs'] == 6
    spectrum = result['full_complement_fusion_bounds']
    assert list(map(F, spectrum['interacting_fusion_cluster_interval'])) == [F(1311, 14), F(1335, 14)]
    assert F(spectrum['fused_anchor_weight_outside_cluster_upper']) == F(64, 19321)
    assert list(map(F, spectrum['deficit_relative_to_twice_isolated_A1u_transition'])) == [
        F(429, 14)-F(118, 1351350), F(453, 14)+F(1, 45675)]
    assert not spectrum['small_subspace_assumed_invariant']
    assert not spectrum['moving_two_packet_collision_executed']
    assert not spectrum['bound_state_below_a_two_particle_scattering_threshold_proved']
    delta = F(result['one_clock_loss_upper'])
    for category, probability in zip(result['complete_original_categories'], (F(99, 100), F(4, 5))):
        weight = category['original_weight_lower_factored']
        assert weight['negative_exponent'] == 32406376
        assert F(weight['rational_multiplier']) == probability*(1-delta)**2
        assert F(category['complete_thirty_pair_menu_each_row_multiplier']) == probability*(1-delta)**2/30 > F(1, 40)
    assert result['complete_pair_and_probe_outcomes'] == 90
    assert result['reference_outputs_and_all_failures_retained']
    assert result['extreme_sector_supported_weight_vanishes_on_selected_refinement']
    assert not result['initial_extreme_sector_weight_conditioned_out']
    assert not result['identified_physical_particle_binding_or_mass_proved']
    assert not result['relativistic_dispersion_or_conserved_spatial_momentum_proved']
    assert result['physical_mass_table_GeV'] == []


def test_equal_work_control_is_a_distinct_apparatus_on_the_same_selected_source(aggregate):
    certificate = aggregate['primary_same_source_equal_work_elastic_control']
    assert certificate == strict_json.load(ROOT / 'model/particle_primary_equal_work_control/results.json')
    original = aggregate['primary_critical_profile_work_refinement']['members']
    assert len(certificate['members']) == len(original) == 4
    for level, bound, critical in zip((0, 1, 2, 8), certificate['members'], original):
        result = bound['result']
        source = critical['result']
        assert bound['result_sha256'] == contract.digest(result)
        assert result['source_critical_bound_sha256'] == critical['result_sha256']
        for key in ('member_id', 'source_id', 'state_id', 'reference_id', 'clock_id'):
            assert result[key] == source[key]
        a = F(1, 2**level)
        lo, hi = map(F, result['exact_coin_outward_interval'])
        assert 0 < lo <= hi <= (F(1) if level == 0 else F(45, 32)*a)
        assert F(result['original_nonelastic_probability_upper']) == hi
        assert result['ideal_transfer_limit'] == 'delta_0'
        assert result['source_state_reference_H_and_clock_unchanged']
        assert not result['apparatus_identical_to_critical_phase_menu']
        assert result['original_work_interval'] == [
            source['profiles'][0]['original_signed_work_lower'],
            source['profiles'][0]['original_signed_work_upper']]
        lower_work = F(result['original_work_interval'][0])
        lower_moment, upper_moment = map(F, result['original_second_moment_interval'])
        assert lower_moment == lower_work**2/hi
        assert upper_moment <= 182250/a
        assert result['same_critical_complete_clock_sha256'] == contract.digest(source['complete_two_energy_clock'])
        assert F(result['complete_reported_work_lower']) == lower_work-F(result['complete_clock_mean_error_upper']) > 0
        assert not result['uniform_energy_weighted_integrability']
        assert not result['critical_phase_collapse_proved']
        assert not result['original_critical_current_covariance_evaluated']
        assert result['physical_mass_table_GeV'] == []
