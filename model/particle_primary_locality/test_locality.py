from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json

import pytest

from model.particle_common_contract import contract
from . import locality


def test_selected_source_has_actual_scaled_action_and_cutoff_domain():
    for level in (0, 1, 2):
        member = contract.member(level)
        source = locality.selected_source_bound(member)
        a = F(member['source']['spacing'])
        assert source['family']['weighted_rate_times_spacing_upper'] == F(5232, 7)
        assert source['actual_member']['magnetic_coefficient'] == F(1, 7)/a
        assert source['actual_member']['electric_coefficient'] == F(105, 2)/a
        assert source['actual_member']['weighted_interaction_rate_upper'] == F(5232, 7)/a
        assert source['reference_idle_H'] == 0
        assert source['product_link_cutoff_preserves_term_support']
        assert not source['global_Gauss_projection_used_as_local_Hamiltonian']
        assert not source['per_link_word_cutoff_enters_propagation_rate']


def test_native_full_fock_bond_constant_is_reused_without_one_particle_reduction():
    rows = locality.cone.native_wilson_bonds()
    assert len(rows) == 3
    assert all(row['two_site_one_body_dimension'] == 120 for row in rows)
    assert all(row['full_Fock_bond_norm_times_spacing'] == 30 for row in rows)


def test_actual_open_box_incidence_satisfies_bound_at_nonunit_spacing():
    for side in (2, 3, 4):
        output = locality.incidence_control(side, spacing=F(2, 5), magnetic=F(1, 7))
        assert output['bond_count'] == 3*side**2*(side-1)
        assert output['plaquette_count'] == 3*side*(side-1)**2
        assert output['actual_support_weighted_norm_bound_max'] <= F(5232, 7)*F(5, 2)
        assert output['all_open_box_interactions_enumerated']


def test_plaquette_geometry_uses_link_owners_and_real_box_boundaries():
    result = locality.plaquette(2, (1, 0, 0), (1, 2))
    assert set(result['owner_support']) == {(1, 0, 0), (1, 1, 0), (1, 0, 1)}
    assert len(result['oriented_links']) == len(result['vertices']) == 4
    assert sum(row[2] for row in result['oriented_links']) == 0
    with pytest.raises(ValueError):
        locality.plaquette(2, (1, 1, 0), (1, 2))
    with pytest.raises(ValueError):
        locality.plaquette(2, (1, 0, 0), (0, 1))
    with pytest.raises(ValueError):
        locality.plaquette(2, (1, 0, 0), (2, 1))


def test_exact_short_time_exponential_encloses_independent_high_precision_value():
    with localcontext() as ctx:
        ctx.prec = 80
        for x in (F(0), F(1, 7), F(654, 4375), F(2)):
            result = locality.exponential_increment_upper(x)
            upper = result['exp_minus_one_upper']
            exact = (Decimal(x.numerator)/Decimal(x.denominator)).exp()-1
            assert Decimal(upper.numerator)/Decimal(upper.denominator) >= exact
            lower = result['positive_sum']
            assert Decimal(lower.numerator)/Decimal(lower.denominator) <= exact
    with pytest.raises(ValueError):
        locality.exponential_increment_upper('3')


def test_actual_complete_boundary_outputs_have_nontrivial_finite_bounds():
    for level, expected_bits in enumerate((2, 7, 30, 124)):
        member = contract.member(level)
        experiment = locality.complete_boundary_experiment(member)
        assert experiment['complete_output_half_trace_upper']['negative_exponent'] == expected_bits
        assert experiment['graph_support_distance'] == member['source']['side']-1
        assert experiment['observer_support_cell_count'] == 1
        assert experiment['all_observer_records_and_quantum_reference_retained']
        assert experiment['remote_field_and_private_record_traced']
        assert not experiment['number_of_record_labels_multiplies_error_bound']
        assert not experiment['outcome_postselection']
        assert not experiment['full_spatial_experiment_executed']
    first = locality.complete_boundary_experiment(contract.member(0))
    assert first['observer_instrument']['outcome_count'] == 1331


def test_fixed_physical_separation_and_time_use_actual_member_geometry():
    member = contract.member(6)
    experiment = locality.complete_boundary_experiment(member, remote_anchor=(64, 0, 0))
    assert experiment['assigned_physical_support_separation'] == 1
    assert experiment['duration'] == F(1, 10000)
    assert experiment['rational_exponential_margin'] == F(3067, 8750)
    assert experiment['complete_output_half_trace_upper']['negative_exponent'] == 21
    assert experiment['source_bound']['actual_member']['spacing'] == F(1, 64)


def test_exact_equal_time_and_inside_cone_domains_are_not_confused():
    member = contract.member(0)
    zero = locality.complete_boundary_experiment(member, duration='0')
    assert zero['complete_output_half_trace_upper']['exact_zero']
    interior = locality.complete_boundary_experiment(member, duration='1')
    assert interior['rational_exponential_margin'] < 0
    assert interior['complete_output_half_trace_upper']['negative_exponent'] == 0
    assert not interior['complete_output_half_trace_upper']['exact_zero']
    with pytest.raises(ValueError):
        locality.complete_boundary_experiment(member, duration='-1')
    with pytest.raises(ArithmeticError):
        locality.complete_boundary_experiment(member, remote_anchor=(0, 0, 0))


def test_typed_result_rejects_changed_clock_cutoff_and_forged_claim():
    member = contract.member(0)
    bound = locality.bind_member(member)
    locality.verify_bound(member, bound)
    for section, field, value in (
        ('source', 'word_projection', 'global_total_degree'),
        ('source', 'electric_coefficient', '1'),
        ('reference', 'idle_H', 1),
        ('clock', 'Schrodinger', 'exp(-2 i t H)'),
    ):
        changed = deepcopy(member)
        changed[section][field] = value
        with pytest.raises(ValueError):
            locality.common_contract_payload(changed)
    changed = deepcopy(bound)
    changed['result']['joint_state_observable_refinement_proved'] = True
    changed['result_sha256'] = contract.digest(changed['result'])
    with pytest.raises(ValueError):
        locality.verify_bound(member, changed)


def test_complete_certificate_is_json_safe_and_keeps_physical_scope():
    report = locality.certificate()
    encoded = json.dumps(report, allow_nan=False)
    assert len(encoded) < 60000
    assert report['same_source_locality_premise_proved']
    assert not report['interacting_relativistic_continuum_proved']
    assert not report['full_spatial_field_run_executed']
    payload = report['selected_member']['result']
    assert not payload['Lorentz_covariance_or_measured_light_speed_proved']
    assert not payload['particle_mass_dispersion_or_scattering_identified']
