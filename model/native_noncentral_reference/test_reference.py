"""Independent controls for exact fair preparation, path ranks and modular gluing."""
from fractions import Fraction as F
from itertools import product, combinations
from math import prod

import numpy as np
import pytest

from . import reference as r, packet


@pytest.mark.parametrize('total', range(1, 12))
def test_all_columns_faithful_noncentral_and_complete_acceptance(total):
    matrix = r.p.block_numerators(total)
    weights = r.occupation_weights(total)
    budget = r.filter_budget(total)
    target = r.filtered_reference(total)
    assert sum(target) == 1
    assert min(target) >= F(1, 8192)
    assert budget['acceptance'] >= F(1, 12)
    assert budget['acceptance'] + budget['rejection'] == 1
    assert int(matrix.min()) >= r.DENOMINATOR // 8192
    assert np.all(matrix.max(axis=0) > matrix.min(axis=0))
    for column in range(matrix.shape[1]):
        accepted = tuple(F(budget['path_multiplicity'] * w, r.DENOMINATOR) for w in weights)
        assert all(int(matrix[i, column]) >= budget['path_multiplicity'] * w for i, w in enumerate(weights))
        assert sum(accepted) == budget['acceptance']
        assert tuple(x / budget['acceptance'] for x in accepted) == target
        assert sum(F(int(matrix[i, column]), r.DENOMINATOR) - x for i, x in enumerate(accepted)) == budget['rejection']


def test_tree_weight_protocol_is_literal_and_fixed():
    assert r.tree_weights() == (1, 2, 3, 4, 3, 2, 3, 2, 3, 3, 2, 2)
    assert r.filter_budget(2)['partition_weight'] == 409
    assert r.filter_budget(2)['acceptance'] == F(7771, 65536)
    assert r.filter_budget(2)['rejection'] == F(57765, 65536)
    assert r.filter_budget(2)['source_eigenvalue_floor'] == F(2, 409)


@pytest.mark.parametrize('initial', [1, 3, 63, 2047])
def test_native_forward_coin_count_matches_reference_and_backward_suffix(initial):
    population = {initial: 1}
    for edge in r.MIXING:
        following = {}
        for state, count in population.items():
            for coin in (0, 1):
                output = r.p.readout.native.image(state, edge, coin, r.p.PORTS)[0]
                following[output] = following.get(output, 0) + count
        population = following
    states = r.p.masks(initial.bit_count())
    expected = tuple(F(population.get(state, 0), r.DENOMINATOR) for state in states)
    assert expected == r.direct_reference(initial)
    final = states[len(states) // 2]
    counts = r.suffix_counts(initial.bit_count(), final)
    assert counts[0][states.index(initial)] == population[final]
    assert all(0 <= value <= 1 << (60-t) for t, row in enumerate(counts) for value in row)


@pytest.mark.parametrize('initial,final', [(1, 2), (3, 288), (63, 126), (2047, 4094)])
def test_rank_unrank_endpoint_and_boundaries(initial, final):
    states = r.p.masks(initial.bit_count())
    count = r.suffix_counts(initial.bit_count(), final)[0][states.index(initial)]
    for rank in (0, 1, count // 2, count-1):
        coins = r.unrank(initial, final, rank)
        assert r.path_rank(initial, coins) == (final, rank, count)
    with pytest.raises(ValueError):
        r.unrank(initial, final, count)


def test_lexicographic_path_rank_by_independent_small_enumeration():
    # Vary the final eight coins after a common52-coin prefix. Each endpoint
    # list is a consecutive final block in its lexicographic endpoint order.
    initial = 3
    grouped = {}
    for suffix in product((0, 1), repeat=8):
        coins = (0,)*52 + suffix
        final, rank, count = r.path_rank(initial, coins)
        grouped.setdefault(final, []).append(rank)
    for ranks in grouped.values():
        assert ranks == list(range(min(ranks), max(ranks)+1))


def test_actual_tree_branch_and_reference_factorization_with_external_coherence():
    trace = r.p.local_trace(3, 0, (0,)*11)
    rho = {(3, 0, 3, 0): F(1, 2), (3, 1, 3, 1): F(1, 2),
           (3, 0, 3, 1): F(1, 2), (3, 1, 3, 0): F(1, 2),
           (3, 0, 5, 1): F(1, 7)}
    branch = r.p.readout.apply_branch(r.p.FAMILY, trace['slots'], rho)
    endpoint = trace['terminal_mask']
    assert branch == {(endpoint, a, endpoint, b): F(1, 4096) for a in (0, 1) for b in (0, 1)}
    actual = r.p.unread(branch, r.MIXING)
    sigma = r.direct_reference(endpoint)
    expected = {(mask, a, mask, b): value / 4096
                for mask, value in zip(r.p.masks(2), sigma, strict=True) for a in (0, 1) for b in (0, 1)}
    assert actual == expected


def test_complete_native_acceptance_and_rejection_boundary_words():
    values = packet.example_words()
    assert [value['flag']['status'] for value in values] == ['ACCEPT', 'ACCEPT', 'REJECT']
    for value in values:
        actual = value['actual_native_word']
        assert actual['attempts'] == 71
        assert actual['exported_records'] == 70
        assert actual['retained_records'] == 1
        assert actual['raw_weight'] == str(F(1, 2**71))


def test_constant_code_and_zero_words_are_not_accepted():
    trace = r.p.local_trace(0, 1, (0,)*71)
    result = r.classify(trace['slots'][:11], trace['slots'][11:])
    assert result['status'] == 'CONSTANT_CODE'
    assert not result['local_total_measured']
    assert not result['accepted']
    full = r.p.local_trace(r.p.FULL, 1, (0,)*71)
    assert full['slots'] == trace['slots']
    assert r.classify((0,)*11, trace['slots'][11:])['status'] == 'ZERO_RAW_WORD'
    assert r.classify(trace['slots'][:11], (0,)*60)['status'] == 'ZERO_RAW_WORD'
    # A syntactically live but impossible mixing contrast remains zero.
    wrong = list(trace['slots'][11:])
    wrong[0] = 14
    assert r.classify(trace['slots'][:11], tuple(wrong))['status'] == 'ZERO_RAW_WORD'


@pytest.mark.parametrize('total', [1, 2, 6, 11])
def test_regional_marginals_compose_without_independent_state_replacement(total):
    sigma = r.filtered_reference(total)
    large = (0, 1, 2, 5, 7)
    small = (0, 2, 7)
    first = r.marginal(total, sigma, large)
    second = [F()] * 8
    positions = tuple(large.index(i) for i in small)
    for mask, value in enumerate(first):
        second[sum(((mask >> bit) & 1) << i for i, bit in enumerate(positions))] += value
    assert tuple(second) == r.marginal(total, sigma, small)
    assert sum(first) == 1


def test_modular_ratios_are_context_independent_and_nontrivial():
    # Port0/1 local swap has a factor1/2 in every q2 exterior context.
    for outside in range(2, 12):
        assert r.modular_ratio(2, 1, 2, 1 << outside) == F(1, 2)
    sigma = r.filtered_reference(2)
    marginal = r.marginal(2, sigma, (0, 1))
    assert marginal[1] / marginal[2] == F(1, 2)
    bigger = r.marginal(2, sigma, (0, 1, 5))
    assert bigger[1] / bigger[2] == bigger[5] / bigger[6] == F(1, 2)


def test_disjoint_number_preserving_regions_commute_on_same_sector():
    states = r.p.masks(2)
    index = {state: i for i, state in enumerate(states)}
    def exchange(a, b):
        matrix = np.zeros((66, 66), dtype=np.int64)
        for state in states:
            if ((state >> a) & 1) != ((state >> b) & 1):
                matrix[index[state ^ ((1 << a) | (1 << b))], index[state]] = 1
        return matrix
    left, right = exchange(0, 1), exchange(5, 7)
    assert np.array_equal(left @ right, right @ left)
    assert np.any(left @ right)
    weights = np.array(r.occupation_weights(2), dtype=np.int64)
    assert np.any(weights[:, None]*left != left*weights[None, :])


def test_raw_program_budget_and_mutation_isolation():
    program = r.p.program(blocks=1)
    assert len(program['commands']) == 71
    assert program['declared_resources']['attempts_upper'] == 71
    first = r.filter_budget(2)
    first['acceptance'] = 0
    assert r.filter_budget(2)['acceptance'] == F(7771, 65536)


@pytest.mark.parametrize('bad', [True, 0, 12, -1, 2.0, '2'])
def test_nonconstant_sector_domain_fails_closed(bad):
    with pytest.raises(ValueError):
        r.filter_budget(bad)


def test_rank_validation_rejects_boolean_coins_and_out_of_sector_endpoint():
    with pytest.raises(ValueError):
        r.path_rank(3, (True,)+(0,)*59)
    with pytest.raises(ValueError):
        r.suffix_counts(2, 1)
    with pytest.raises(ValueError):
        r.unrank(3, 1, 0)


def test_exact_receipt_and_geometric_scope():
    packet.check()
    value = packet.produce()
    assert value['scope']['exact_faithful_noncentral_reference_per_nonconstant_tree_branch']
    assert not value['scope']['observer_tree_weight_protocol_is_unique_vacuum_selection']
    assert not value['scope']['full_archive_source_factorization']
    assert not value['scope']['small_ball_geometric_stress_or_Einstein_premises_derived']
