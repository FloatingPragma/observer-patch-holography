"""Inherited-action checks and source/extension boundary controls."""
import copy

import pytest
import sympy as s

from model.common_source import calibrated_matter
from model.particle_joint_reconstruction import candidate


@pytest.fixture(scope='module')
def result():
    return candidate.build_candidate()


def test_actual_mass_kernel_and_calibrated_ancestry(result):
    catalogue = result['inherited_running_coefficient_arm']['catalogue']
    difference = candidate._matrix(catalogue['kernel_checks']['K45']) - calibrated_matter.mass_kernel()
    assert all(s.simplify(value) == 0 for value in difference)
    assert catalogue['calibration'] == calibrated_matter.parameters()
    assert catalogue['kernel_checks']['K45_rank_exact'] == 42
    assert catalogue['kernel_checks']['K45_nullity_exact'] == 3
    assert len(catalogue['fermions']) == 12
    assert len(catalogue['bosons']['rows']) == 6
    assert result['scope']['independent_predicted_reference_mass_count'] == 0


def test_neutrino_extension_keeps_minimal_arm_and_input_ancestry(result):
    minimal = result['inherited_running_coefficient_arm']['catalogue']['fermions'][-3:]
    assert [row['tree_running_mass_exact'] for row in minimal] == ['0']*3
    extended = result['conditional_physical_completion_arm']['neutrino_producer']
    assert extended['calibration'] == calibrated_matter.neutrino_parameters()
    masses = list(map(s.sympify, extended['neutrino_mass_parameters_GeV_exact']))
    packet = extended['calibration']['oscillation_inputs']
    assert s.simplify((masses[1]**2-masses[0]**2)*10**18-s.Rational(packet['delta_m21_squared_eV2'])) == 0
    assert s.simplify((masses[2]**2-masses[0]**2)*10**18-s.Rational(packet['delta_m31_squared_eV2'])) == 0
    assert extended['same_source_current_LH_invariant_identities_exact'] == 12
    assert not extended['low_energy_to_MZ_EFT_running_or_matching_error_certified']
    assert not result['conditional_physical_completion_arm']['all_reference_targets_numerically_matched']


def test_independent_boson_action_checks_and_soft_photon(result):
    control = result['quadratic_cross_checks']
    assert control['independent_inherited_gauge_Hessians_agree']
    assert control['boson_on_shell_degrees'] == 28
    energies = [s.Rational(row['quadratic_energy_GeV']) for row in control['soft_photon_sequence_control']]
    assert all(left > right > 0 for left, right in zip(energies, energies[1:]))
    assert not control['positive_global_gap_of_photon_containing_continuum_scenario']
    damaged = copy.deepcopy(result['inherited_running_coefficient_arm']['catalogue'])
    next(row for row in damaged['bosons']['rows'] if row['name'] == 'Z')['tree_mass_squared_exact'] = '1'
    with pytest.raises(ArithmeticError, match='boson action'):
        candidate.quadratic_consistency(damaged, result['inherited_symbolic_chiral_action'])


def test_primary_carrier_state_clock_are_not_transferred(result):
    boundary = result['extension_boundary']
    assert boundary['primary_CAR_modes_per_cell'] == 60
    assert boundary['extension_static_spin_internal_CAR_modes'] == 90
    assert not boundary['same_carrier_or_basis_change']
    assert not boundary['state_reference_clock_transfer_supplied']
    assert boundary['primary_GeV_conversion'] is None
    assert not result['primary_member']['state']['exact_Gibbs_state']
    assert all(not row['proved_on_primary'] for row in result['assumptions'])
    assert all(not row['witness_supplied_by_candidate'] for row in result['assumptions'])


def test_actual_intertwiner_leakage_and_pure_gap_scope(result):
    control = result['full_rectangular_matching_control']
    assert control['exact']['generator_intertwining_exact']
    assert control['leakage_control_compression_equals_field_generator']
    assert s.Rational(control['leakage_control']['generator_residual_frobenius_squared']) > 0
    assert not control['leakage_control']['generator_intertwining_exact']
    branch = result['pure_Yang_Mills_branch']
    assert branch['finite_transfer_controls']['finite_coupled_source']['full_excitation_span']
    assert not branch['same_as_matter_containing_product_group_primary']
    assert not branch['strong_sector_gap_transferred_to_full_theory']
    assert not branch['Clay_problem_solved']


def test_receipt_replay_rejects_promoted_physical_claim(result):
    damaged = copy.deepcopy(result)
    damaged['scope']['interacting_four_dimensional_QFT_constructed'] = True
    damaged['scope']['independent_predicted_reference_mass_count'] = 18
    damaged['pure_Yang_Mills_branch']['full_photon_containing_theory_has_positive_global_gap'] = True
    with pytest.raises(ValueError, match='fresh inherited calculations'):
        candidate.validate_candidate(damaged)
