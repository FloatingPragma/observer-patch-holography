"""Exact source matrices, cube orbits, CAR mutations and retained band bounds."""
from copy import deepcopy
from fractions import Fraction as F

import pytest
import sympy as s

from model.particle_common_contract import contract
from model.cubic_particle_channels import group
from model.particle_primary_survival import survival
from . import channels as c


@pytest.fixture(scope='module')
def member():
    return contract.member(1)


def test_actual_native_blocks_and_both_internal_charts():
    value = c.native_number_control()
    assert [(x['Q6'], x['CAR_modes_per_cell']) for x in value['blocks']] == [(-4, 12), (-3, 8), (1, 24), (2, 12), (6, 4)]
    assert value['native_generators_checked'] == 12
    assert value['both_mask_and_packet_charts_checked']


def test_actual_clifford_covariance_for_every_binary_cube_element():
    value = c.rotation_control()
    assert value['binary_rotations'] == 48
    assert value['proper_spatial_rotations'] == 24
    assert c.spin_lift(group.negative(group.ONE)) == -s.eye(4)


def test_spin_lift_product_on_independent_cube_generators():
    generators = [q for q in group.elements() if q[0] == group.Q(0, F(1, 2)) and
                  sum(bool(x) for x in q[1:]) == 1][:3]
    assert len(generators) == 3
    for q in group.elements():
        for r in generators:
            assert c.apparatus.clean(c.spin_lift(q)*c.spin_lift(r)) == c.spin_lift(group.multiply(q, r))


@pytest.mark.parametrize('side', [4, 8, 32])
def test_selected_plaquette_has_regular_24_element_orbit(side):
    value = c.orbit_control(side)
    assert value['orbit_size'] == 24
    assert value['coefficient_representation_irreps'] == {'A1': 1, 'A2': 1, 'E': 2, 'T1': 3, 'T2': 3}
    assert value['A1_average_coefficient'] == '1/24'
    assert not value['original_single_plaquette_is_cubic_scalar']
    # Independent Haar-edge condition for distinct loop orthogonality.
    def edges(vertices):
        return {tuple(sorted((x, y))) for x in vertices for y in vertices
                if sum(abs(a-b) for a, b in zip(x, y)) == 1}
    sets = [edges(p) for p in value['orbit']]
    assert all(len(x) == 4 for x in sets)
    assert all(left-right and right-left for i, left in enumerate(sets) for right in sets[i+1:])


def test_ghost_boundary_reduction_is_needed_and_interior_graph_rotates():
    value = c.boundary_graph_control(3)
    assert value['owned_links'] == 81
    assert value['one_link_ghost_vertices'] == 27
    assert value['interior_links'] == 54
    assert not value['ambient_owned_link_allocation_claimed_cube_covariant']


def test_exact_CAR_control_detects_pairing_and_interblock_mutations():
    value = c.car_control()
    assert value['separate_block_commutators_zero'] == [True, True]
    assert value['cross_block_mutation_preserves_total']
    assert value['cross_block_mutation_breaks_individual']
    assert value['pair_mutation_breaks_total']


def test_scalar_blank_norm_and_band_coefficients_are_source_derived(member):
    value = c.selected_channel(member)
    old = survival.selected_witness(member)
    a = F(member['source']['spacing']); width = F(old['source_spectral_width_upper'])
    assert F(value['cubic_scalar_blank_commutator_norm_squared']) == F(old['blank_commutator_norm_squared'])/24
    k = F(4949, 480); r = F(23, 24)**4
    coefficient = F(value['finite_apparatus_original_band_weight_lower']['coefficient'])
    assert coefficient == (F(37, 96)*r*k-F(1, 2))/(a*a*width*width)
    assert coefficient > 0
    assert value['original_positive_band']['lower']['coefficient'] == str(1/a)
    assert value['delta_N_q'] == {'-4': 0, '-3': 0, '1': 0, '2': 0, '6': 0}
    assert value['rest_mass_or_GeV'] is None
    assert not value['massless_spectrum_excluded']


def test_same_source_binding_rejects_changed_observable_and_clock(member):
    bound = c.bind_member(member)
    changed = deepcopy(bound)
    changed['result']['cubic_scalar_norm_upper'] = '2'
    with pytest.raises(ValueError):
        c.verify_bound(member, changed)
    changed_member = deepcopy(member)
    changed_member['reference']['idle_H'] = 1
    with pytest.raises(ValueError):
        c.selected_channel(changed_member)


def test_cubic_scalar_does_not_determine_a_unique_continuum_spin():
    assert group.decompose(lambda q: group.spin_character(q, 0)) == {'A1': 1}
    spin4 = group.decompose(lambda q: group.spin_character(q, 8))
    assert spin4['A1'] == 1
    assert sum(group.DIMENSIONS[group.IRREPS.index(k)]*v for k, v in spin4.items()) == 9


def test_inactive_primary_member_does_not_gain_a_scalar_response():
    with pytest.raises(ValueError, match='inactive'):
        c.selected_channel(contract.member(0))
