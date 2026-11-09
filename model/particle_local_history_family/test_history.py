from fractions import Fraction as F
import json
from pathlib import Path

import pytest

from model.particle_observer_outputs import field_route
from model.particle_field_causality import bounds
from model.particle_local_gauge_state import state
from . import history as h


def test_exact_binary_ceiling_on_both_sides_of_powers():
    for n in range(-20, 21):
        power = F(2**n) if n >= 0 else F(1, 2**(-n))
        assert h.ceil_log2(power) == n
        assert h.ceil_log2(power*F(999, 1000)) == n
        assert h.ceil_log2(power*F(1001, 1000)) == n+1
    with pytest.raises(ValueError):
        h.ceil_log2(0)


def test_retained_patch_keeps_shared_vertices_and_outgoing_boundary_links():
    patch = h.local_marginal_patch([(0, 0, 0), (1, 0, 0)])
    assert len(patch['retained_outgoing_links']) == 6
    assert len(patch['independent_Haar_vertices']) == 7
    assert ((0, 0, 0), (1, 0, 0)) in patch['retained_outgoing_links']
    assert (2, 0, 0) in patch['independent_Haar_vertices']
    # Translating the patch translates one shared incidence structure,
    # rather than introducing independent gauges for adjacent link ends.
    shifted = h.local_marginal_patch([(3, -1, 7), (4, -1, 7)])
    assert {tuple(v-d for v, d in zip(x, (3, -1, 7)))
            for x in shifted['independent_Haar_vertices']} == set(patch['independent_Haar_vertices'])


def test_fixed_time_volume_plan_has_a_compact_integer_certificate():
    source = state.density_plan(degree=0)
    plan = h.volume_plan(source, '1', 24)
    assert plan['duration'] == 1
    assert plan['local_record_error_binary_exponent_upper'] <= -24
    assert plan['complete_record_TV_upper'] == F(1, 2**24)
    assert plan['box_radius'] > 10**40
    assert plan['locality_constants']['plaquette_norm'] == source['declared_plaquette_coefficient']
    assert not plan['duration_shrunk_with_box_or_volume']
    assert not plan['large_exponential_or_2_to_radius_materialized']


def test_logarithmic_bound_encloses_owning_locality_on_small_control():
    source = {'Legendre_degree': 0, 'spacing': F(1), 'declared_plaquette_coefficient': F(1, 100)}
    plan = h.volume_plan(source, '1/1000', 12, intermediate_instruments=2)
    original = bounds.volume_bound(plan['locality_constants'], F(1, 1000),
                                   plan['box_radius'], interventions=2)
    assert original['complete_local_history_half_trace_upper'] <= plan['complete_record_TV_upper']
    assert plan['periodic_side_sufficient'] >= 2*(plan['box_radius']+5)+1


def test_cutoff_uses_total_history_time_and_support_once():
    source = state.density_plan(degree=0)
    single = h.volume_plan(source, '1', 16)
    plan = h.cutoff_plan(source, single, 16)
    M, x = plan['first_omitted_Dyson_order'], plan['integrated_bounded_action_norm_upper']
    assert M >= 4*x and M >= 3*(16+2)
    assert plan['total_word_cutoff'] == plan['initial_total_word_degree_upper']+plan['word_step']*(M-1)
    assert plan['whole_history_error_binary_exponent_upper'] <= -16
    assert plan['degree_preserving_pointers_are_interleaved_before_Dyson_count']
    # Same box and duration with extra interventions has the same global
    # Dyson budget; the initial degree is not restarted between intervals.
    several = {**single, 'intermediate_instruments': 5}
    other = h.cutoff_plan(source, several, 16)
    assert other['first_omitted_Dyson_order'] == M
    assert other['total_word_cutoff'] == plan['total_word_cutoff']
    assert other['evolution_intervals'] == 6


def test_factorial_tail_bound_against_existing_exact_tail_routine():
    for x in (F(0), F(1, 8), F(2), F(17, 3)):
        M = max(1, h.ceil_fraction(4*x), 18)
        tail = 2*field_route.exponential_tail_upper(x, M, extra=8)
        bound = F(1, 2**(M//3-2))
        assert tail <= 4*F(3, 4)**M <= bound


def test_native_full_Fock_current_is_nonstatic_with_complete_initial_number_law():
    control = h.source_current_control()
    assert control['single_bond_current_variance'] == F(240, 289)
    assert control['origin_number_derivative_mean_square'] == F(1440, 289)
    assert control['two_time_squared_number_increment_second_derivative_at_zero'] == F(2880, 289)
    probabilities = control['initial_complete_number_probabilities']
    assert len(probabilities) == 61 and sum(probabilities) == 1
    assert sum(n*p for n, p in enumerate(probabilities)) == 30
    assert sum((n-30)**2*p for n, p in enumerate(probabilities)) == 15
    assert control['two_time_number_menu_labels'] == 3721


def test_domain_mutations_rejected():
    source = state.density_plan(degree=0)
    for cells in ([], [(0, 1)], [(0, 1, 2.0)]):
        with pytest.raises(ValueError):
            h.local_marginal_patch(cells)
    with pytest.raises(ValueError):
        h.volume_plan(source, duration=0)
    with pytest.raises(ValueError):
        h.volume_plan(source, duration=0.1)
    with pytest.raises(ValueError):
        h.volume_plan(source, observed_cells=2)
    with pytest.raises(ValueError):
        h.volume_plan(source, intermediate_instruments=-1)


def test_certificate_keeps_limits_references_and_resources_distinct():
    cert = h.certificate()
    assert cert['schema'] == 'oph.local_history_family.v1'
    assert cert['total_record_error_upper'] == '1/8388608'
    assert cert['local_polynomial_record_limit_constructed_at_fixed_spacing']
    assert not cert['full_canonical_purifiers_identified_across_volumes']
    assert not cert['full_dynamical_exterior_output_compared']
    assert not cert['roughness_bound_reused_after_local_interventions']
    assert not cert['exact_sign_or_relativistic_limit_claimed']
    assert not cert['large_finite_simulation_executed']
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in cert['source_inputs'])
    assert len(json.dumps(cert)) < 50000
