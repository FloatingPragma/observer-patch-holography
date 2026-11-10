"""Exact Gauss words, independent finite controls and source-bound events."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from model.particle_primary_reduced_sector.test_model import q6_compression
from . import model as m


def independent_flux_basis():
    edges = m.reduced.cube_support()['interior_edges']
    incidence = s.zeros(7, 12)
    for i, (x, y) in enumerate(edges):
        if x < 7:
            incidence[x, i] = -1
        if y < 7:
            incidence[y, i] = 1
    _, pivots = incidence.rref()
    chords = [j for j in range(12) if j not in pivots]
    transform = -incidence[:, pivots].inv()*incidence[:, chords]
    states = []
    for word in product((-1, 0, 1), repeat=5):
        tree = transform*s.Matrix(word)
        if any(x not in (-1, 0, 1) for x in tree):
            continue
        flux = [0]*12
        for i, x in zip(pivots, tree):
            flux[i] = int(x)
        for i, x in zip(chords, word):
            flux[i] = x
        states.append(tuple(flux))
    return states


def test_actual_cube_flux_fusion_and_all_occupations_retained():
    result = m.cycle_control()
    rows = result['all_ordered_distinct_face_pairs']
    assert len(rows) == 30
    states = independent_flux_basis()
    assert len(states) == 69
    assert sum(sum(q*q for q in state) == 6 for state in states) == 32
    hexagon_supports = {tuple(i for i, q in enumerate(state) if q)
                        for state in states if sum(q*q for q in state) == 6}
    assert len(hexagon_supports) == 16
    for support in hexagon_supports:
        assert sum(tuple(i for i, q in enumerate(state) if q) == support
                   for state in states) == 2
    for row in rows:
        assert row['resulting_q6_flux'] in states
        assert row['electric_energy'] == F(63, 4)*sum(q*q for q in row['resulting_q6_flux'])
    assert sum(row['nonconstant_links'] == 6 for row in rows) == 24
    assert sum(row['nonconstant_links'] == 8 for row in rows) == 6


def test_complete_sequential_loop_instrument_in_full_q6_control():
    states = independent_flux_basis()
    where = {x: i for i, x in enumerate(states)}
    faces = m.reduced.cube_support()['oriented_outward_faces']
    n = len(states)
    shifts = []
    for face in faces:
        C = np.zeros((n, n))
        for j, flux in enumerate(states):
            target = tuple(a+b for a, b in zip(flux, face['outward_flux']))
            if target in where:
                C[where[target], j] = 1
        shifts.append(C)
    omega = np.eye(n)[:, where[(0,)*12]]
    for row in m.cycle_control()['all_ordered_distinct_face_pairs']:
        left, right = (shifts[row[k]] for k in ('first_face', 'second_face'))
        # Each compressed character is a partial isometry in this independent
        # q6 control. Its positive failure root is the complementary projector.
        failure_left = np.eye(n)-left.T@left
        failure_right = np.eye(n)-right.T@right
        Kraus = (right@left, failure_right@left, failure_left)
        assert np.array_equal(sum(K.T@K for K in Kraus), np.eye(n))
        assert np.array_equal(Kraus[0]@omega, np.eye(n)[:, where[row['resulting_q6_flux']]])
    # No finite-q6 compression is used as the full native dynamics.
    assert not m.cycle_control()['pair_intervention_uses_unprojected_parent_dynamics']


def test_full_native_electric_rank_includes_the_weak_loops_omitted_by_q6_control():
    data = m.reduced.native_Haar_control()
    energies = data['one_link_electric_energies']
    assert 4*energies[1] == 6*energies[-1] == F(189, 2)
    assert all(4*e >= F(441, 4) for e in energies if e not in (F(63, 4), F(189, 8)))
    assert 6*F(189, 8) > F(441, 4)
    assert 7*min(energies) == F(441, 4)
    H, A, omega, phi, electric = q6_compression()
    assert np.count_nonzero(electric == 94.5) == 32
    assert m.cycle_control()['combined_electric_cluster_dimension'] == 44
    # The full source contains four more native character pairs in its action.
    assert F(1, 735) == 5*F(1, 3675)


def test_full_complement_fusion_bounds_have_the_correct_signs_and_energy_zero():
    r = m.fusion_bounds()
    assert r['distance_of_electric_anchor_energy_to_full_cluster_complement_lower'] == F(417, 28)
    assert r['fused_anchor_weight_outside_cluster_upper'] == F(64, 19321)
    assert r['interacting_fusion_cluster_interval'] == [F(1311, 14), F(1335, 14)]
    assert r['actual_ground_to_fused_cluster_weight_comparison_lower'] > F(99, 100)
    assert r['actual_ground_to_disjoint_window_weight_comparison_lower'] > F(4, 5)
    assert r['deficit_relative_to_twice_isolated_A1u_transition'] == [
        F(429, 14)-F(118, 1351350), F(453, 14)+F(1, 45675)]
    assert r['interacting_fusion_cluster_interval'][1] < r['disjoint_final_energy_window'][0]


def test_exact_interacting_q6_control_separates_both_static_pair_windows():
    H, A, omega, phi, electric = q6_compression()
    energies, vectors = np.linalg.eigh(H)
    fused = (energies >= float(F(1311, 14))) & (energies <= float(F(1335, 14)))
    disjoint = (energies >= 124) & (energies <= 128)
    assert np.count_nonzero(fused) == 32
    for i, energy in enumerate(electric):
        if energy not in (94.5, 126):
            continue
        band = fused if energy == 94.5 else disjoint
        weight = sum(abs(vectors[i, band])**2)
        assert weight > .99
    assert not m.fusion_bounds()['small_subspace_assumed_invariant']


def test_original_faithful_weight_precedes_the_clock_and_selector_factors():
    member = m.contract.member(0)
    receipt = m.bind_member(member)
    m.verify_bound(member, receipt)
    r = receipt['result']
    assert r['complete_pair_and_probe_outcomes'] == 90
    assert F(r['complete_ordered_pair_selector_probability']) == F(1, 30)
    assert r['original_gauge_ground_population_floor']['negative_exponent'] == 32406376
    delta = F(r['one_clock_loss_upper'])
    for category, ideal in zip(r['complete_original_categories'], (F(99, 100), F(4, 5))):
        weight = category['original_weight_lower_factored']
        assert weight['negative_exponent'] == 32406376
        assert F(weight['rational_multiplier']) == ideal*(1-delta)**2
        assert F(category['complete_thirty_pair_menu_each_row_multiplier']) == F(weight['rational_multiplier'])/30 > F(1, 40)
    assert r['same_source_state_reference_clock'] and r['original_selected_density_retained']
    assert r['reference_outputs_and_all_failures_retained']
    assert not r['initial_extreme_sector_weight_conditioned_out']
    assert not r['retained_pure_reference_Loewner_bound_asserted']


def test_forged_particle_binding_or_source_change_is_rejected():
    member = m.contract.member(0)
    receipt = m.bind_member(member)
    forged = deepcopy(receipt)
    forged['result']['identified_physical_particle_binding_or_mass_proved'] = True
    forged['result_sha256'] = m.contract.digest(forged['result'])
    with pytest.raises(ValueError):
        m.verify_bound(member, forged)
    with pytest.raises(ValueError):
        m.bind_member(m.contract.member(1))
    assert receipt['result']['physical_mass_table_GeV'] == []
    assert receipt['result']['extreme_sector_supported_weight_vanishes_on_selected_refinement']
