from fractions import Fraction as F

import pytest
import sympy as s

from . import preparation as p


def test_actual_native_metric_and_electric_casimirs():
    rows = {row['q6Y']: row for row in p.native_casimirs()}
    assert {q: row['trace_metric_C2'] for q, row in rows.items()} == {
        1: s.Rational(21, 20), -4: s.Rational(4, 5),
        2: s.Rational(7, 10), -3: s.Rational(9, 20), 6: s.Rational(3, 10)}
    assert sum(row['dimension'] for row in rows.values()) == 15


def test_complete_character_unitary_prepares_source_state():
    data = p.one_cell_seed()
    u = data['unitary']
    assert p.clean(u.H*u) == s.eye(11)
    assert p.clean(u*s.eye(11)[:, 0]) == data['state']
    assert p.clean(u*u) == s.eye(11)
    assert data['one_link_energy'] == s.Rational(33, 9100)


def test_all_original_seed_records_and_complements_are_retained():
    data = p.seed_records()
    assert data['complete_record_count'] == 1728
    assert data['positive_record_count'] == 1331
    assert data['zero_record_count'] == 397
    assert sum(row['probability'] for row in data['original_records']) == 1
    assert data['mean_electric_energy'] == s.Rational(99, 9100)
    assert not data['postselection_used']
    assert not data['small_gap_weight_certified_by_L1_seed']


def test_exact_full_native_gap_closure_uses_spatial_beta():
    result = p.native_gap_closure()
    assert result['H_W_squared_spectral_ranks'] == {0: 16, 1: 28, 3: 16}
    assert result['symbol_auxiliary_trace_probe_probabilities'] == {
        0: s.Rational(4, 15), 1: s.Rational(7, 15), 3: s.Rational(4, 15)}
    assert not result['actual_quantum_point_coordinate_preparation_claimed']


def test_finite_degree_planner_and_factored_resources():
    result = p.choose_preparation(1, 20)
    assert result['root_degree'] == 8192
    assert result['per_link_word_degree'] == 67108864
    assert result['total_word_degree'] == 201326592
    assert result['uncapped_binary_exponent'] == -16231
    assert result['code_dimension_upper'] == {'base': 451, 'exponent': 201326592}
    assert not result['large_preparation_coefficients_materialized']


def test_volume_cost_is_not_silently_removed():
    small = p.polynomial_gap_plan(1, 8192)
    larger = p.polynomial_gap_plan(2, 8192)
    assert larger['links'] == 8*small['links']
    assert larger['electric_energy_upper'] == 8*small['electric_energy_upper']
    assert larger['uncapped_binary_exponent'] > small['uncapped_binary_exponent']
    assert p.choose_preparation(2, 20)['root_degree'] > small['root_degree']


def test_same_hamiltonian_clock_controls_gap_and_record_error():
    plan = p.choose_preparation(1, 20)
    result = p.same_clock_interval(plan)
    assert result['common_duration'] > 0
    assert result['common_duration'] < F(1, 10**22)
    assert result['small_gap_probability_upper_through_interval'] == F(1, 2**18)
    assert F(0) < result['complete_final_record_TV_poly_vs_exact_sign_upper'] < F(1, 10**22)
    assert result['same_electric_plaquette_state_and_clock']
    assert not result['intermediate_measurements_covered']


def test_seed_is_not_accepted_as_a_good_gap_preparation():
    with pytest.raises(ValueError):
        p.same_clock_interval(p.polynomial_gap_plan(1, 1))


def test_prepared_plaquette_moment_meets_its_energy_barrier():
    result = p.plaquette_moment_plan(1, 20)
    target = result['requested_all_time_bad_weight']
    assert result['positive_plaquette_sum_mean_upper'] <= target/(2*48000)
    assert result['large_plaquette_uncapped_binary_exponent'] < 0
    assert result['HJL_bad_gap_plaquette_sum_lower'] == F(1, 48000)


def test_energy_confinement_uses_one_distinct_fixed_hamiltonian():
    result = p.energy_confinement_plan()
    assert result['all_time_small_gap_probability_upper'] <= F(1, 2**40)
    assert result['declared_plaquette_coefficient'] > 1
    assert result['kappa_1_100_energy_bound_on_this_preparation'] == 1
    assert result['small_kappa_failure_is_inconclusive_not_a_dynamical_no_go']
    assert result['no_repeated_gap_filtering']
    assert result['same_bound_for_compressed_bad_effect_on_every_containing_word_code']
    assert result['compressed_bad_effect_is_not_truncated_Wilson_gap']


def test_all_time_gap_does_not_make_detector_error_time_independent():
    plan = p.energy_confinement_plan()
    first = p.all_time_final_record_bound(plan, 1)
    second = p.all_time_final_record_bound(plan, 2)
    assert first['complete_final_record_TV_poly_vs_exact_sign_upper'] < F(1, 1000)
    assert second['complete_final_record_TV_poly_vs_exact_sign_upper'] == 2*first[
        'complete_final_record_TV_poly_vs_exact_sign_upper']
    assert not first['spatial_locality_or_cutoff_error_included']
    assert not first['intermediate_measurement_backaction_covered']


@pytest.mark.parametrize('call', [
    lambda: p.polynomial_gap_plan(0, 1),
    lambda: p.polynomial_gap_plan(1, True),
    lambda: p.choose_preparation(1, 0),
    lambda: p.same_clock_interval(p.choose_preparation(1, 20), spacing='0'),
    lambda: p.same_clock_interval(p.choose_preparation(1, 20), occupied_modes=61),
    lambda: p.same_clock_interval(p.choose_preparation(1, 20), plaquette=0.1),
])
def test_theorem_domain_checks(call):
    with pytest.raises(ValueError):
        call()


def test_certificate_keeps_cutoff_and_physical_domains_explicit():
    result = p.certificate()
    assert result['physical_particle_assignments'] == []
    assert result['periodic_boundary_required_for_free_gap_bound']
    assert result['finite_unitary_preparation_not_retry_postselection']
    assert not result['higher_degree_polynomial_coefficients_or_tape_executed']
    cutoff = result['same_source_total_degree_cutoff']
    assert cutoff['initial_degree'] == 201326592
    assert cutoff['first_possible_omitted_Dyson_order'] == 5
