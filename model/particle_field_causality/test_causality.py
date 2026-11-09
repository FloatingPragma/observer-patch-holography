from fractions import Fraction as F
import itertools

import pytest
import sympy as s

from . import bounds
from model.particle_observer_outputs import field_route
from model.particle_candidate_atlas.interacting_chiral import creator


def test_ball_counts_against_literal_cubic_geometry():
    for radius in range(6):
        actual = sum(sum(abs(x) for x in point) <= radius
                     for point in itertools.product(range(-radius, radius+1), repeat=3))
        assert actual == bounds.ball_size(radius)


def test_actual_native_full_fock_bond_norm():
    result = bounds.native_bond_check()
    assert len(result) == 3
    assert all(row['full_Dirac_M15_Fock_bond_norm'] == '120/17' for row in result)
    assert bounds.interaction_constants(0, '1', '1/100')[
        'weighted_local_interaction_upper'] == F(72816, 425)


def test_generic_bound_and_sharp_source_bound_are_distinct():
    generic = bounds.interaction_constants(0, sharpen_zero=False)
    sharp = bounds.interaction_constants(0)
    assert generic['weighted_local_interaction_upper'] > sharp['weighted_local_interaction_upper']
    assert bounds.interaction_constants(2)['coordinate_word_degree'] == 5
    assert not generic['outer_volume_in_bound']


def test_polynomial_support_uses_actual_source_gamma():
    gamma, _ = field_route.source_data.load()['gamma']
    beta = gamma[3]
    gamma = s.Matrix(gamma[0]).applyfunc(lambda x: s.nsimplify(x, rational=True))
    beta = s.Matrix(beta).applyfunc(lambda x: s.nsimplify(x, rational=True))
    sites = 7
    kernel = s.zeros(4*sites)
    for x in range(sites):
        kernel[4*x:4*x+4, 4*x:4*x+4] = 2*beta
    hop = -beta*(s.eye(4)-gamma)/2
    for x in range(sites-1):
        kernel[4*x:4*x+4, 4*(x+1):4*(x+2)] = hop
        kernel[4*(x+1):4*(x+2), 4*x:4*x+4] = hop.H
    polynomial = field_route.sign_polynomial(kernel, 1)
    for x in range(sites):
        for y in range(sites):
            if abs(x-y) > 3:
                assert polynomial[4*x:4*x+4, 4*y:4*y+4].is_zero_matrix
    assert not polynomial[0:4, 8:12].is_zero_matrix


def test_odd_car_fields_require_anticommutator():
    a, b = creator(2, 0).H, creator(2, 1).H
    assert a*b+b*a == s.zeros(4)
    assert a*b-b*a != s.zeros(4)
    n = a.H*a
    assert n*b-b*n == s.zeros(4)


def test_remote_complete_detector_bound_is_nontrivial():
    constants = bounds.interaction_constants(0, '1', '1/100')
    result = bounds.commutator_bound(constants, '1/64', 32)
    assert F(0) < result['complete_local_record_TV_upper'] < F(1, 10**7)
    assert bounds.commutator_bound(constants, 0, 32)['complete_local_record_TV_upper'] == 0
    assert bounds.commutator_bound(constants, '1/64', 33)[
        'uncapped_normalized_bound'] == result['uncapped_normalized_bound']/2


def test_clock_and_common_energy_scale_cancel_in_bound():
    first = bounds.interaction_constants(0, '1', '1/100')
    second = bounds.interaction_constants(0, '1/3', '3/100')
    left = bounds.commutator_bound(first, '1/64', 32)
    right = bounds.commutator_bound(second, '1/192', 32)
    assert left['uncapped_normalized_bound'] == right['uncapped_normalized_bound']


def test_complete_history_volume_search_is_uniform_in_exterior():
    constants = bounds.interaction_constants(0, '1', '1/100')
    result = bounds.choose_box(constants, '1/64', '1/1000000', interventions=2)
    assert result['complete_local_history_half_trace_upper'] <= F(1, 10**6)
    assert result['source_radius'] >= 5
    assert not result['unbounded_exterior_volume_in_bound']
    assert 'entire radius-R box' in result['initial_state_requirement']
    assert not result['large_Fock_allocation_used']
    no_interventions = bounds.volume_bound(constants, '1/64', result['source_radius'])
    assert result['complete_local_history_half_trace_upper'] == 5*no_interventions[
        'complete_local_history_half_trace_upper']


def test_total_degree_projection_is_not_local():
    result = bounds.global_cutoff_counterexample()
    assert result['parent_disjoint_shift_commutator'] == '0'
    assert result['compressed_commutator_norm'] == '1'
    assert result['witness_input'] == [0, 1]


def test_representation_error_is_added_twice_for_two_experiments():
    assert bounds.cutoff_record_bound('1/10', '1/20', '1/30') == F(11, 60)
    assert bounds.cutoff_record_bound('9/10', '1/10', '1/10') == 1


@pytest.mark.parametrize('call', [
    lambda: bounds.ball_size(-1),
    lambda: bounds.interaction_constants(True),
    lambda: bounds.interaction_constants(0, 1.0),
    lambda: bounds.interaction_constants(0, '0'),
    lambda: bounds.commutator_bound(bounds.interaction_constants(0), 1, 0),
    lambda: bounds.choose_box(bounds.interaction_constants(0), 1, 1),
    lambda: bounds.volume_bound(bounds.interaction_constants(0), 1, 3),
    lambda: bounds.cutoff_record_bound('2', '0', '0'),
])
def test_invalid_theorem_domain_is_rejected(call):
    with pytest.raises(ValueError):
        call()


def test_certificate_preserves_physical_boundaries():
    result = bounds.certificate()
    assert result['physical_particle_assignments'] == []
    assert not result['uniform_refinement_or_Lorentz_boost_claimed']
    assert not result['strict_Hamiltonian_light_cone_claimed']
    assert not result['norm_continuous_dynamics_on_all_link_BH_claimed']
