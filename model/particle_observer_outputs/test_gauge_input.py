"""Finite operational couplings and exact physical-parameter dependency checks."""
from fractions import Fraction as F
from pathlib import Path

import sympy as s

from model.particle_observer_outputs import gauge_input as g


def test_fixed_link_source_coefficients_are_observer_readable():
    result = g.fixed_link_readout()
    assert result['source_default_is_unambiguous']
    assert result['dimensionless_epsilon_over_abs_J'] == '1'
    assert not result['identified_as_MSbar_alpha_s']
    assert [r['block_dimension'] for r in result['complete_rational_readout_witnesses']] == [1, 2, 3, 6]
    for row in result['complete_rational_readout_witnesses']:
        for key in ('hopping_squared', 'electric_squared', 'electric_over_hopping_squared'):
            lo, hi = map(F, row['intervals'][key])
            assert lo <= 1 <= hi


def test_scalar_program_has_a_real_noncommuting_interaction():
    result = g.fixed_scalar_interface()
    assert result['default_scalar_charge_coupling'] == '1/4'
    assert result['finite_register_dimension_in_control'] == 16
    assert result['nonzero_complete_interaction_hopping_commutator']
    assert result['source_matrix_float_diagnostic_max_difference'] < 1e-14
    assert not result['equal_to_electroweak_g_or_Higgs_lambda_claimed']


def test_complete_Higgs_transition_reads_its_own_vertex_strength():
    result = g.quantum_Higgs_response()
    assert result['inferred_pair_strength_squared'] == '56'
    assert result['inferred_number_cost_squared'] == '1'
    assert result['count_inverse_gap'] == '15'
    for law in result['complete_fixture_joint_probabilities']:
        assert len(law) == 4
        assert sum(map(s.sympify, law.values())) == 1
    assert s.sympify(result['complete_empty_input_probability']) == s.Rational(56, 225)
    coherent = s.sympify(result['coherent_source_pair_strength_squared'])
    assert s.simplify(coherent-(425+100*s.sqrt(2))/1944) == 0
    assert s.simplify(s.sympify(result['coherent_source_gap_squared'])-1-4*coherent) == 0
    assert not result['population_records_fix_cost_sign']
    assert result['coherent_controller_phase_angle'] == '0'
    assert not result['equal_to_Higgs_background_v_or_quartic_lambda_claimed']


def test_same_symmetry_premises_allow_different_dimensionless_observables():
    result = g.action_dependency()
    a, b = result['positive_parameter_counterchoices']
    assert (a['Higgs_squared_over_W_squared'], b['Higgs_squared_over_W_squared']) == ('8', '16')
    assert (a['Z_squared_over_W_squared'], b['Z_squared_over_W_squared']) == ('2', '5')
    assert a['alpha_s_in_declared_curvature_convention'] != b['alpha_s_in_declared_curvature_convention']
    assert a['neutral_photon_mass_squared'] == b['neutral_photon_mass_squared'] == '0'
    assert result['dimensionless_ratios_invariant_under_common_unit_change']
    assert not result['architecture_uniqueness_required']


def test_conditional_gauge_solution_is_reused_with_all_matching_premises():
    result = g.conditional_gauge_ray()
    assert result['producer_and_all_inherited_sources_authenticated']
    assert not result['conditional_solution_discarded']
    assert not result['value_identified_as_alpha_s_at_MZ']
    assert len(result['premises']) == 4
    assert len(result['physical_matching_obligations']) == 5
    assert result['empirical_inputs_in_owner'] == []
    values = result['conditional_values']
    assert 0 < values['alphaY'] < values['alpha2'] < values['alpha3'] < 1
    assert not result['owning_entropy_tail_proof_reexecuted_here']


def test_repair_and_spin_results_keep_their_positive_conclusions_and_domains():
    result = g.finite_response_reuse()
    assert result['scientific_source_bindings_authenticated']
    assert result['static_family_rank'] == 1
    assert result['uniform_repair_family_ratios'] == [1, 1, 1]
    assert result['canonical_normalization_preserves_two_zero_family_coefficients']
    assert result['fixed_graph_example']['spin_exchange'] == '1'
    assert result['finite_spin_limit_is_interacting']
    assert not result['spin_exchange_is_QCD_alpha_s']
    assert not result['relativistic_chiral_continuum_inferred']


def test_no_empirical_mass_inputs_and_all_sources_resolve():
    result = g.certificate()
    assert result['measured_mass_inputs_added'] == []
    assert result['architecture_inputs_allowed']
    assert result['one_common_energy_unit_can_be_declared']
    assert not result['common_unit_determines_dimensionless_mass_ratios']
    assert len(result['eliminated_unknowns']) == 4
    root = Path(__file__).resolve().parents[2]
    for path in result['source_inputs']:
        assert (root/path).is_file(), path
