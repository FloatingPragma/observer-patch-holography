from fractions import Fraction as F
import json
from pathlib import Path

import pytest
import sympy as s

from model.particle_field_causality import bounds
from model.particle_observer_outputs import field_route
from . import experiment as e


def test_actual_native_six_bond_escape_coefficient():
    result = e.native_escape_coefficient()
    assert result['squared_escape_coefficient'] == F(2880, 289)
    assert len(result['bonds']) == 6
    assert all(row['spin_rank'] == 2 and row['internal_dimension'] == 15
               and row['squared_escape_coefficient'] == '480/289'
               for row in result['bonds'])
    assert len({(row['direction'], row['sign']) for row in result['bonds']}) == 6


def test_whole_marginal_degree_and_boundary_reference():
    support = e.marginal_support(5, 7)
    assert support['retained_outgoing_links'] == 3*11**3
    assert support['initial_total_word_degree_upper'] == 3*11**3*7
    assert support['full_Fock_dimension'] == {'base': 2, 'exponent': 60*11**3}
    assert not support['box_only_singlet_projection_performed']
    assert not support['bra_degree_added_to_initial_ket_cutoff']
    # An entangled degree-one ket has a mixed degree-one local marginal.
    # Its nonzero charge is compensated by the conjugate reference.
    psi = s.Matrix([s.sqrt(3)/3, 0, 0, s.sqrt(6)/3])
    rho = psi*psi.H
    marginal = s.Matrix(2, 2, lambda i, j: sum(rho[2*i+k, 2*j+k] for k in range(2)))
    assert marginal == s.diag(s.Rational(1, 3), s.Rational(2, 3))
    assert s.diag(0, -1, 1, 0)*psi == s.zeros(4, 1)
    assert s.trace(marginal*s.diag(0, 1)) == s.Rational(2, 3)


def test_rational_locality_envelope_bounds_owning_theorem():
    constants = bounds.interaction_constants(0, '1', '1/100')
    t = F(1, 100000)
    compact = e.small_clock_volume_bound(constants, t, 5, 21)
    inherited = bounds.volume_bound(constants, t, 5)
    assert inherited['complete_local_history_half_trace_upper'] <= compact['complete_local_record_TV_upper'] < 1
    assert compact['boundary_separation_lower'] == 1
    assert compact['same_plaquette_coefficient'] == F(1, 100)
    assert not compact['dynamical_exterior_state_comparison_claimed']
    assert not compact['open_Wilson_gap_assumed']


def test_boundary_and_clock_domain_mutations_are_rejected():
    constants = bounds.interaction_constants(0, '1', '1/100')
    with pytest.raises(ValueError, match='wraparound'):
        e.small_clock_volume_bound(constants, F(1, 100000), 5, 19)
    with pytest.raises(ValueError, match='wraparound'):
        e.small_clock_volume_bound(constants, F(1, 100000), 5, 22)
    with pytest.raises(ValueError, match='rational exponential'):
        e.small_clock_volume_bound(constants, F(1), 5, 21)
    with pytest.raises(ValueError, match='exact rational'):
        e.small_clock_volume_bound(constants, 0.001, 5, 21)


def test_one_action_and_clock_in_all_composed_steps():
    plan = e.experiment_plan()
    source = plan['periodic_preparation_and_action']
    kappa = source['declared_plaquette_coefficient']
    assert kappa != F(1, 100)
    assert source['all_time_small_gap_probability_upper'] <= F(1, 2**20)
    assert kappa == plan['locality_constants']['plaquette_norm']
    assert kappa == plan['finite_cutoff']['same_plaquette_coefficient']
    assert kappa == plan['periodic_to_open_local_record_comparison']['same_plaquette_coefficient']
    t = plan['clock']['macro_tick']
    assert t == plan['finite_cutoff']['same_duration']
    assert t == plan['periodic_sign_comparison']['duration']
    assert t == plan['periodic_to_open_local_record_comparison']['same_duration']
    assert t == plan['finite_complete_occupation_readout']['duration']


def test_initial_degree_and_first_missing_dyson_order_are_compatible():
    plan = e.experiment_plan()
    cut, support = plan['finite_cutoff'], plan['open_marginal_and_reference']
    assert cut['initial_degree'] == support['initial_total_word_degree_upper']
    assert cut['cutoff']-cut['initial_degree'] == 12
    assert cut['coordinate_step'] == 4
    assert cut['first_possible_omitted_Dyson_order'] == 4
    assert cut['actual_integrated_interaction_norm_upper'] <= F(cut['integrated_interaction_norm_upper'])
    expected = field_route.cutoff_bound(cut['cutoff'], cut['initial_degree'], 4, F(1, 128))
    refined = 2*field_route.exponential_tail_upper(
        cut['actual_integrated_interaction_norm_upper'], 4, extra=1)
    assert F(cut['complete_record_TV_error_upper']) == refined
    assert refined <= F(expected['complete_record_TV_error_upper'])
    assert not cut['total_degree_compression_claimed_exactly_local']


def test_complete_number_menu_and_nonzero_finite_event():
    plan = e.experiment_plan()
    result = plan['finite_complete_occupation_readout']
    rows = result['complete_outcomes']
    assert [row['origin_occupation'] for row in rows] == list(range(61))
    assert sum(row['local_Fock_projector_rank'] for row in rows) == 2**60
    assert all(0 <= row['probability_interval'][0] <= row['probability_interval'][1] <= 1 for row in rows)
    assert sum(row['probability_interval'][0] for row in rows) <= 1 <= sum(row['probability_interval'][1] for row in rows)
    assert rows[59]['probability_interval'][0] > 0
    assert rows[60]['probability_interval'][0] > F(999, 1000)
    assert result['joint_probability_origin_occupation_0_through_58_upper'] == result['state_vector_Taylor_remainder_upper']**2
    assert result['unused_pointer_probabilities'] == [0, 0, 0]
    assert result['sum_of_all61_probabilities'] == 1
    assert not result['conditioning_or_histogram_renormalization']


def test_taylor_enclosure_rejects_unresolved_tick():
    with pytest.raises(ValueError, match='first hopping'):
        e.occupation_enclosure(F(1, 100), F(100))
    with pytest.raises(ValueError, match='wrong sign'):
        e.occupation_enclosure(0, F(100))
    with pytest.raises(ValueError, match='exact rational'):
        e.occupation_enclosure(0.01, F(100))


def test_complete_parent_error_and_rare_event_resolution_are_distinct():
    plan = e.experiment_plan()
    budget = plan['comparison_budget']
    total = sum(budget[key] for key in (
        'periodic_exact_sign_to_periodic_polynomial',
        'periodic_polynomial_to_open_polynomial',
        'open_polynomial_to_finite_word_code'))
    assert total == budget['total_complete_local_record_TV_upper'] < F(1, 10**6)
    assert plan['periodic_exact_sign_p59_interval'][0] == 0
    assert not plan['periodic_exact_sign_rare_event_positivity_resolved']
    assert not plan['interacting_relativistic_continuum_claimed']
    assert not plan['physical_particle_mass_or_species_claimed']


def test_compiler_and_preparation_errors_add_as_complete_channels():
    plan = e.experiment_plan()
    base = plan['comparison_budget']['total_complete_local_record_TV_upper']
    assert e.implementation_budget(plan) == base
    assert e.implementation_budget(plan, preparation_error=F(1, 1000),
                                   step_error=F(2, 1000), detector_error=F(3, 1000),
                                   initial_state_error=F(4, 1000)) == base+F(1, 100)
    assert e.implementation_budget(plan, step_error=F(1)) == 1
    with pytest.raises(ValueError, match='exact complete-channel'):
        e.implementation_budget(plan, step_error=0.001)


def test_certificate_keeps_factored_resources_and_source_bindings():
    cert = e.certificate()
    assert cert['schema'] == 'oph.prepared_field_readout.v1'
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in cert['source_inputs'])
    assert not cert['numerically_materialized_source_or_propagator']
    assert not cert['independent_measured_particle_input_used']
    assert not cert['narrower_native_repair_tape_implementation_claimed']
    assert cert['experiment']['finite_cutoff']['coordinate_dimension_upper']['exponent'] > 10**20
    assert len(json.dumps(cert)) < 100000


def test_mutated_preparation_domain_is_rejected():
    for radius, degree in ((0, 1), (5, 0), (5.0, 2), (5, True)):
        with pytest.raises(ValueError):
            e.marginal_support(radius, degree)
    with pytest.raises(ValueError):
        e.experiment_plan(radius=4)
    with pytest.raises(ValueError):
        e.experiment_plan(target_bad_bits=0)
