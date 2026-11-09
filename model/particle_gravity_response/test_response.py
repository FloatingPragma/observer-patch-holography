"""Complete-carrier source energy accounting and action-response controls."""
import pytest
import sympy as s

from model.particle_gravity_response import response as r


def test_actual_quantum_flux_source_energy_partition_and_current():
    source, terms, result, currents = r.rotor_data()
    assert source['H'].rows == 208
    assert result['initial_component_energies'] == [9, 2, 0]
    assert result['initial_component_energy_rates'] == [18, 0, -18]
    assert result['nonzero_unordered_current_pairs'] == [[0, 2], [1, 2]]
    assert result['hopping_coefficient_response_initial_time_slope'] == 6
    assert result['weighted_action_work_identity']
    assert r.clean(sum(terms, s.zeros(208))-source['H']) == s.zeros(208)
    assert r.clean(currents[0, 2]+currents[2, 0]) == s.zeros(208)


def test_partition_rejects_changed_action_nonhermitian_terms_or_lost_input():
    source, terms, _, _ = r.rotor_data()
    with pytest.raises(ValueError, match='same source Hamiltonian'):
        r.partition_certificate(source['H']+s.eye(208), terms, source['state'])
    with pytest.raises(ValueError, match='Hermitian'):
        r.partition_certificate(source['H'], (terms[0]*s.I,*terms[1:]), source['state'])
    with pytest.raises(ValueError, match='normalized'):
        r.partition_certificate(source['H'], terms, source['state']/2)


def test_actual_full_scalar_source_has_complete_energy_current_detector():
    source, state, result, projectors = r.scalar_data()
    assert source['H'].rows == 307
    assert result['current_detector_ranks'] == [1, 1, 305]
    assert result['initial_complete_current_probabilities'] == [1, 0, 0]
    assert result['initial_component_energy_rates'] == [-2*s.sqrt(14),2*s.sqrt(14)]
    assert result['prepared_state_is_gauge_singlet']
    assert not result['full307_carrier_is_all_singlet']
    dark = r.scalar.gauge_invariant_dark_probe(source)['vector']
    assert [r.expect(dark, p) for p in projectors] == [0, 0, 1]
    assert r.expect(state, source['H']) == s.Rational(1,2)


@pytest.mark.parametrize('time',(0,s.pi/30,s.pi/15))
def test_complete_probabilities_and_integrated_flux_use_actual_H(time):
    source, _, _, current_pvm = r.scalar_data()
    state = r.scalar_wavefunction(time)
    law = r.scalar_history(time)
    occupation = [source['vacuum_projector'],source['bright_projector'],source['dark_projector']]
    assert r.expect(state,s.eye(307)) == 1
    for p, expected in zip(occupation,law['occupation_probabilities'],strict=True):
        assert s.simplify(r.expect(state,p)-expected) == 0
    for p, expected in zip(current_pvm,law['current_probabilities'],strict=True):
        assert s.simplify(r.expect(state,p)-expected) == 0
    assert sum(law['occupation_probabilities']) == 1
    assert sum(law['current_probabilities']) == 1
    assert s.simplify(law['number_energy']+law['vertex_energy']-s.Rational(1,2)) == 0
    assert s.simplify(law['number_energy']-s.Rational(1,2)+law['integrated_number_to_vertex_flux']) == 0
    assert r.expect(state,source['H']) == s.Rational(1,2)


def test_pair_cost_response_holds_preparation_and_clock_fixed():
    result = r.action_response_certificate()
    assert result['fixture_occupied_probability_derivative'] == 2*s.sqrt(14)/3375
    assert result['all_occupation_outputs_retained']
    assert result['same_preparation_and_common_time_required']
    assert not result['geometric_metric_identification_of_pair_cost']
    original, state, _, _ = r.scalar_data()
    changed = r.scalar.sector_data(s.diag(1,2,3), s.Rational(7,5))
    assert r.clean(changed['H']-original['H']-s.Rational(2,5)*original['occupation']) == s.zeros(307)
    assert changed['bright_vector'] == original['bright_vector']
    assert r.expect(state,changed['H']) == s.Rational(7,10)


def test_one_program_product_detector_retains_every_flux_and_scalar_output():
    result = r.paired_initial_records()
    assert result['declared_carrier_dimension'] == 63856
    assert len(result['all_outcomes']) == 26
    nonzero = [row for row in result['all_outcomes'] if row['probability']]
    assert len(nonzero) == 4
    assert all(row['probability'] == s.Rational(1,4) for row in nonzero)
    assert {row['electric_flux'] for row in nonzero} == {0,6}
    assert not result['rotor_scalar_interaction_declared']
    assert not result['postselection_or_failure_conditioning']


def test_energy_data_do_not_supply_unmeasured_einstein_tensor_components():
    result = r.gravity_interface()
    assert result['energy_only_tomography_null_direction']['rest_energy_charge'] == '0'
    assert len(result['Einstein_requirements_not_instantiated']) == 6
    assert not result['finite_Einstein_tensor_bound_invoked']
    assert not result['repair_stationary_entropy_identified_with_geometric_entropy']
    assert not result['Einstein_dynamics_or_physical_Newton_constant_derived']
    assert not result['architecture_uniqueness_required']
