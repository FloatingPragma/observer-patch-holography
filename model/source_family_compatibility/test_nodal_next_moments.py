"""Independent occupation controls for actual local nodal moment readouts."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product

import numpy as np
import pytest

from . import composite, momentum_response as momentum, nodal_next_moments as n
from .local_matter import annihilator


@lru_cache(maxsize=None)
def gammas(sites):
    matrices = []
    for mode in range(2*sites):
        a = np.array(annihilator(sites, mode).tolist(), dtype=np.complex128)
        matrices.extend((a+a.conj().T, -1j*(a-a.conj().T)))
    return tuple(matrices)


def mask_matrix(mask, sites, hermitian=False):
    out = np.eye(1 << (2*sites), dtype=np.complex128)
    basis = gammas(sites)
    for j in range(4*sites):
        if mask >> j & 1:
            out = out@basis[j]
    return out*((1j)**(mask.bit_count()*(mask.bit_count()-1)//2) if hermitian else 1)


def polynomial_matrix(operator, sites):
    out = np.zeros((1 << (2*sites),)*2, dtype=np.complex128)
    for mask, (re, im) in operator.items():
        out += complex(re, im)*mask_matrix(mask, sites)
    return out


def word_matrix(word, sites):
    index = {point: j for j, point in enumerate(sorted(sites))}
    mask = sum(1 << (4*index[row[:3]]+row[3]) for row in word)
    return mask_matrix(mask, len(sites), hermitian=True)


def literal_source_matrix(sites):
    """Independent creation/annihilation formula for the declared source."""
    sites = tuple(sorted(sites))
    count = len(sites)
    a = [np.array(annihilator(count, j).tolist(), dtype=np.complex128)
         for j in range(2*count)]
    identity = np.eye(1 << (2*count))
    h = np.zeros_like(identity, dtype=np.complex128)
    for j, x in enumerate(sites):
        nu, nd = a[2*j].conj().T@a[2*j], a[2*j+1].conj().T@a[2*j+1]
        h += 8*(nu-identity/2)@(nd-identity/2)
        for k, y in enumerate(sites):
            if k <= j:
                continue
            difference = tuple(y[l]-x[l] for l in range(3))
            if sum(map(abs, difference)) != 1:
                continue
            axis = next(l for l in range(3) if difference[l])
            start = x if difference[axis] > 0 else y
            coefficient = momentum.signed_link(start, axis)/2
            for spin in (0, 1):
                term = a[2*j+spin].conj().T@a[2*k+spin]
                h += coefficient*(term+term.conj().T)
    return h


def test_sparse_source_and_four_derivatives_match_full_occupation_matrices():
    sites = ((0, 0, 0), (1, 0, 0), (2, 0, 0))
    _, _, adjacency, terms = n._source(sites)
    source = composite.source((3, 1, 1), F(1), F(8))
    h4 = polynomial_matrix(terms, 3)
    assert np.array_equal(h4, 4*source['H']/source['scale']+24*np.eye(64))
    assert np.array_equal(h4, 4*literal_source_matrix(sites))
    operator = {13: (0, 1), 14: (-1, 0)}
    matrix = 2*source['b'][0]
    assert np.array_equal(matrix, polynomial_matrix(operator, 3))
    for _ in range(4):
        matrix = matrix@h4-h4@matrix
        operator = n._derivative(operator, adjacency, terms)
        assert np.array_equal(matrix, polynomial_matrix(operator, 3))
    contracted = n._b_contraction(operator, sites)
    got = sum((c*mask_matrix(mask, 3, hermitian=True) for mask, c in contracted.items()),
              np.zeros((64, 64), dtype=np.complex128))
    b = source['b'][0]-source['b'][4]
    want = matrix@b.T+b.T@matrix
    assert np.array_equal(got, (want+want.conj().T)/2)


def occupation_transport(sites, perm, signs, shift):
    moved = tuple(tuple(momentum.rotate(x, perm, signs)[j]+shift[j] for j in range(3))
                  for x in sites)
    target = tuple(sorted(moved))
    index = {x: j for j, x in enumerate(target)}
    mode_map, mode_sign = [], []
    for x, y in zip(sites, moved, strict=True):
        rotation = momentum.rotation_phase(x, perm)
        z = momentum.rotate(x, perm, signs)
        translation = (-1 if (shift[1]*z[0]+shift[2]*(z[0]+z[1])) % 2 else 1)
        for spin in (0, 1):
            mode_map.append(2*index[y]+spin)
            mode_sign.append(rotation*translation)
    destinations, phases = [], []
    for state in range(1 << (2*len(sites))):
        occupied = [j for j in range(2*len(sites)) if state >> j & 1]
        permutation = sum(mode_map[i] > mode_map[j] for i in occupied for j in occupied if i < j)
        phase = -1 if permutation % 2 else 1
        for j in occupied:
            phase *= mode_sign[j]
        destinations.append(sum(1 << mode_map[j] for j in occupied))
        phases.append(phase)
    return target, np.array(destinations), np.array(phases)


def transported(matrix, destinations, phases):
    out = np.empty_like(matrix)
    out[destinations[:, None], destinations[None, :]] = phases[:, None]*matrix*phases[None, :]
    return out


def test_all_spatial_lifts_and_odd_unit_shifts_preserve_actual_source_and_words():
    sites = ((0, 0, 0), (1, 0, 0), (1, 1, 0))
    word = ((0, 0, 0, 0), (1, 0, 0, 2), (1, 1, 0, 0), (1, 1, 0, 3))
    h = literal_source_matrix(sites)
    observable = word_matrix(word, sites)
    vector = np.zeros(64)
    vector[[1, 9, 16]] = (1, 2, 3)
    density_numerator = np.outer(vector, vector)
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            shift = (1, -1, 2)
            target, destinations, phases = occupation_transport(sites, perm, signs, shift)
            transformed_word, sign = n.transform_word(word, perm, signs, shift)
            want = sign*word_matrix(transformed_word, target)
            assert np.array_equal(transported(h, destinations, phases), literal_source_matrix(target))
            assert np.array_equal(transported(observable, destinations, phases), want)
            transformed_density = transported(density_numerator, destinations, phases)
            assert np.trace(density_numerator@observable) == np.trace(transformed_density@want)
            inverse = tuple(perm.index(j) for j in range(3))
            inverse_signs = tuple(signs[inverse[j]] for j in range(3))
            inverse_shift = momentum.rotate(tuple(-j for j in shift), inverse, inverse_signs)
            restored, inverse_sign = n.transform_word(transformed_word, inverse, inverse_signs, inverse_shift)
            assert restored == word and sign*inverse_sign == 1


def test_naked_coordinate_permutation_loses_the_actual_gauge_sign():
    word = ((0, 0, 0, 0), (1, 1, 0, 0))
    moved, sign = n.transform_word(word, (1, 0, 2))
    assert moved == word and sign == -1
    assert not np.array_equal(word_matrix(word, ((0, 0, 0), (1, 1, 0))),
                              sign*word_matrix(moved, ((0, 0, 0), (1, 1, 0))))


@pytest.fixture(scope='module')
def receipt():
    return n.certificate()


def test_actual_sparse_source_readouts_and_ground_scope(receipt):
    rows = receipt['actual_local_nodal_moment_readouts']
    assert [row['invariant_Hermitian_words'] for row in rows] == [19, 493, 27965]
    assert [row['maximum_occupied_sites_per_word'] for row in rows] == [3, 4, 5]
    assert all(row['selected_ground_moment'] is None for row in rows)
    assert all(row['quotient_is_invariant_state_expectation_equality'] for row in rows)
    assert all(not row['quotient_is_operator_equality'] for row in rows)
    assert not receipt['source_selected_ground_moments_evaluated']
    assert not receipt['actual_zero_atom_decided']
    assert receipt['physical_species_identified'] == receipt['physical_masses'] == []


@pytest.mark.parametrize('order', (4, 6, 8))
def test_complete_weighted_score_retains_every_setting_and_both_outcomes(receipt, order):
    plan = n.score_plan(order)
    assert sum(row['probability'] for row in plan['settings']) == 1
    assert all(row['probability'] > 0 for row in plan['settings'])
    assert all(abs(row['plus_score']) == plan['score_absolute_bound'] for row in plan['settings'])
    assert plan['both_pointer_outcomes_per_setting_retained']
    assert plan['identity_word_negative_outcome_is_zero']
    assert not plan['large_native_readout_executed']
    assert not plan['state_or_service_error_implied_zero']
    # The tracial functional annihilates every nonempty distinct Majorana word.
    mean = sum(row['probability']*row['plus_score'] for row in plan['settings'] if not row['word'])
    assert mean == n.nodal_moment(order)['tracial_control_constant']
    for row in plan['settings']:
        assert len(row['word']) % 2 == 0 and tuple(sorted(set(row['word']))) == row['word']


def test_onsite_force_identity_does_not_supply_a_ground_norm_identity(receipt):
    assert receipt['onsite_force_identity']['nonzero_anticommutator_norm_difference_words'] == 300
    control = receipt['independent_full_source_control']
    assert control['whole_Fock_dimension'] == 64 and control['full_ground_rank'] == 2
    assert control['onsite_norm_difference_from_alpha_squared'] == F(-216, 805)
    assert control['onsite_center_force_ground_norm_squared']-(F(3, 14)**2)*control['center_force_ground_norm_squared'] == F(-216, 805)
    assert not control['finite_graph_zero_hopping_mode_is_infinite_lattice_momentum']
    assert not control['finite_ground_replaces_default_infinite_ground']


@pytest.mark.parametrize('s,t,u', [(1, 1, 1), (1, 1, 4), (6, F(1203, 98), F(492987, 19208))])
def test_christoffel_bound_satisfies_independent_normal_equations(s, t, u):
    row = n.christoffel(s, t, u)
    _, m2, m4, m6, m8 = row['original_moments_0_2_4_6_8']
    one, b, d = row['polynomial_coefficients_1_E2_E4']
    assert one == 1
    assert m2+b*m4+d*m6 == 0
    assert m4+b*m6+d*m8 == 0
    norm = 1+2*b*m2+2*d*m4+b*b*m4+2*b*d*m6+d*d*m8
    assert norm == row['original_zero_atom_upper']
    assert 0 <= norm < row['degree_two_upper']
    assert not row['selected_source_values_certified_by_helper']


@pytest.mark.parametrize('s', (0, F(3, 2)))
def test_singular_memory_kernel_case_has_no_zero_over_zero(s):
    row = n.christoffel(s, 0, 0)
    assert row['original_zero_atom_upper'] == F(s)/(16+F(s))
    assert row['case'] == 'memory_in_kernel'


@pytest.mark.parametrize('values', [(1, 0, 1), (0, 0, 1), (0, 1, 1), (1, 2, 1),
                                   (-1, 1, 1), (True, 1, 1), (1.0, 1, 1)])
def test_impossible_or_inexact_source_moments_are_rejected(values):
    with pytest.raises(ValueError):
        n.christoffel(*values)


def test_actual_memory_floor_improves_bound_without_a_positive_zero_limit(receipt):
    row = receipt['actual_source_memory_bound']
    assert row['source_t_upper'] == 17082
    assert row['source_u_upper'] == F(33378201, 4)
    eta, lower = row['eta'], row['source_s_strict_lower']
    old = eta*eta*lower*lower/(eta*eta*lower+5175**2)
    assert row['M_eta_lower'] > 1500*old
    assert n.memory_budget(eta/100)['M_eta_lower'] < row['M_eta_lower']/9000
    assert row['lower_bound_is_order_eta_squared'] and not row['zero_frequency_atom_inferred']


def test_public_rationals_do_not_alias_cached_source_results(receipt):
    expected = repr(n.nodal_moment(6))
    changed = n.nodal_moment(6)
    changed['source_observable_norm_upper']._numerator += 1
    assert repr(n.nodal_moment(6)) == expected
    control = n.finite_control()
    expected = repr(n.finite_control())
    control['memory_s_t_u'][0]._numerator += 100
    assert repr(n.finite_control()) == expected


@pytest.mark.parametrize('order', (True, 2, 5, 10, 6.0))
def test_source_readout_order_is_explicit(order):
    with pytest.raises(ValueError):
        n.nodal_moment(order)


def test_symmetry_transform_rejects_odd_repeated_or_inexact_words():
    for word in (((0, 0, 0, 0),), ((0, 0, 0, 0), (0, 0, 0, 0))):
        with pytest.raises(ValueError):
            n.transform_word(word)
    with pytest.raises(ValueError):
        n.transform_word(((0, 0, 0, 0), (1, 0, 0, 0)), (0, True, 2))
