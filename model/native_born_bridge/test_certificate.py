from dataclasses import replace
from fractions import Fraction as Q
from itertools import product
from unittest.mock import patch

import pytest

from . import certificate as c


def uniform():
    return tuple(Q(1, 78) for _ in c.BASIS)


def test_complete_native_tree_matches_original_adaptive_execution():
    tree = c.source_tree()
    actual = {(p.initial, p.final, tuple((r[2], r[3]) for r in p.record), p.conditional_weight)
              for p in tree.paths}
    expected = set()
    for i, state in enumerate(c.BASIS):
        for first_weight, first in c.native.tick(c.native.State(state)):
            for last_weight, last in c.native.tick(first):
                expected.add((i, c.INDEX[last.loads], last.prefix, first_weight*last_weight))
    assert actual == expected
    assert len(tree.paths) == 312
    assert len({(p.record, p.final) for p in tree.paths}) == 312


def test_zero_mass_outcomes_and_joint_terminal_states_are_retained():
    tree = c.source_tree()
    weights = tuple(Q(i == 0) for i in range(78))
    distribution = c.record_law(tree, weights)
    joint = c.execute(tree, weights)
    assert set(distribution) == set(tree.records)
    assert any(w == 0 for w in distribution.values())
    assert len(joint) == len(tree.records)*78
    assert sum(joint.values()) == sum(distribution.values()) == 1


@pytest.mark.parametrize('damage', ['omit', 'duplicate', 'current', 'policy'])
def test_complete_tree_custody_mutations_fail(damage):
    tree = c.source_tree()
    if damage == 'omit':
        bad = replace(tree, paths=tree.paths[1:])
    elif damage == 'duplicate':
        bad = replace(tree, paths=tree.paths+(tree.paths[0],))
    elif damage == 'current':
        p = tree.paths[0]
        a, b, coin, current = p.record[0]
        bad = replace(tree, paths=(replace(p, record=((a, b, coin, current+1),)+p.record[1:]),)+tree.paths[1:])
    else:
        rules = dict(tree.rules)
        rules[()] = (0, 2)
        bad = replace(tree, rules=tuple(sorted(rules.items())))
    with pytest.raises(ValueError):
        c.validate_tree(bad)


def test_source_primitive_mutation_is_caught_by_independent_endpoint_formula():
    with patch.object(c.native, 'native', return_value=(c.BASIS[0], (0, 99))):
        with pytest.raises(ValueError, match='independent balanced repair'):
            c.endpoint(c.BASIS[0], c.native.E, 0)


def test_coarse_graining_pullback_and_coexistent_effect_additivity():
    tree = c.source_tree()
    law = tuple(Q(i+1, 3081) for i in range(78))
    dist = c.record_law(tree, law)
    parts = [frozenset(r for i, r in enumerate(tree.records) if i % 3 == k) for k in range(3)]
    effects = [c.event_effect(tree, p) for p in parts]
    assert all(sum(e[i] for e in effects) == 1 for i in range(78))
    for selected, effect in zip(parts, effects):
        assert sum(dist[r] for r in selected) == sum(law[i]*effect[i] for i in range(78))
    first = {r: Q((i % 3), 7) for i, r in enumerate(tree.records)}
    second = {r: Q((i % 4), 9) for i, r in enumerate(tree.records)}
    summed = {r: first[r]+second[r] for r in tree.records}
    assert c.valuation(dist, summed) == c.valuation(dist, first)+c.valuation(dist, second)
    assert c.valuation(dist, {r: 1 for r in tree.records}) == 1
    assert c.valuation(dist, {r: 0 for r in tree.records}) == 0
    # Point indicators recover every diagonal entry uniquely.
    assert all(c.valuation(dist, {s: int(s == r) for s in tree.records}) == dist[r]
               for r in tree.records)


def test_initial_mixing_and_correlated_classical_reference_preserve_weights():
    tree = c.source_tree()
    first = tuple(Q(i == 0) for i in range(78))
    second = tuple(Q(i == 17) for i in range(78))
    a, b = c.record_law(tree, first), c.record_law(tree, second)
    mixed = tuple(Q(2, 7)*x+Q(5, 7)*y for x, y in zip(first, second))
    actual = c.record_law(tree, mixed)
    reference_joint = {(r, 0): Q(2, 7)*a[r] for r in tree.records}
    reference_joint.update({(r, 1): Q(5, 7)*b[r] for r in tree.records})
    assert all(actual[r] == reference_joint[r, 0]+reference_joint[r, 1] for r in tree.records)
    assert sum(w for (_, ref), w in reference_joint.items() if ref == 0) == Q(2, 7)
    assert sum(w for (_, ref), w in reference_joint.items() if ref == 1) == Q(5, 7)


def test_fair_native_suffix_implements_dyadic_acceptance_without_postselection():
    tree = c.source_tree()
    numerators = {r: i % 5 for i, r in enumerate(tree.records)}
    result = c.dyadic_acceptance(tree, uniform(), numerators, 2)
    dist = c.record_law(tree, uniform())
    for r in tree.records:
        assert result['record_masses'][r, True] == dist[r]*Q(numerators[r], 4)
        assert result['record_masses'][r, False] == dist[r]*Q(4-numerators[r], 4)
    assert result['executed_suffixes'] == 1248
    assert result['accept']+result['reject'] == 1


def test_dyadic_sandwich_bounds_every_exact_record_effect():
    dist = c.record_law(c.source_tree(), uniform())
    values = {r: Q(i % 13, 13) for i, r in enumerate(dist)}
    previous = None
    for bits in (0, 1, 3, 8, 32):
        result = c.dyadic_enclosure(dist, values, bits)
        assert result['lower'] <= result['value'] <= result['upper']
        assert result['upper']-result['lower'] <= Q(1, 2**bits)
        if previous:
            assert previous['lower'] <= result['lower'] <= result['upper'] <= previous['upper']
        previous = result


def test_two_coexistent_effects_share_one_disjoint_native_refinement():
    tree = c.source_tree()
    dist = c.record_law(tree, uniform())
    counts = {r: (i % 3, 1, 3-i % 3) for i, r in enumerate(tree.records)}
    result = c.dyadic_partition(tree, uniform(), counts, 2)
    A, B, reject = result['outcome_masses']
    assert A+B+reject == 1
    union = c.dyadic_acceptance(tree, uniform(), {r: row[0]+row[1] for r, row in counts.items()}, 2)
    assert union['accept'] == A+B
    assert set(label for _, label in result['record_masses']) == {0, 1, 2}
    for (original, suffix, label), weight in result['raw_record_masses'].items():
        assert len(suffix) == 2 and original in dist and label in (0, 1, 2) and weight >= 0
    assert sum(result['raw_record_masses'].values()) == 1
    with pytest.raises(ValueError, match='partition the whole coin domain'):
        c.dyadic_partition(tree, uniform(), {r: (3, 2, -1) for r in tree.records}, 2)


def test_prefix_pushforward_and_explicit_halts():
    short = c.compile_tree(1, {(): c.native.E})
    long = c.source_tree()
    assert c.refinement_certificate(short, long)['all_initial_laws_and_classical_references']
    select = lambda p: None if p and p[0][2] == 0 else c.native.E
    halt = c.compile_tree(2, c.policy(2, select))
    longer = c.compile_tree(3, c.policy(3, select))
    assert len(halt.paths) == 234
    assert {len(p.record) for p in halt.paths} == {1, 2}
    assert c.refinement_certificate(halt, longer)['input_columns_checked'] == 78
    bad = c.compile_tree(2, c.policy(2, lambda p: c.native.F))
    with pytest.raises(ValueError, match='changes an earlier policy'):
        c.refinement_certificate(short, bad)


@pytest.mark.parametrize('bad', [True, .5, Q(-1), 2])
def test_initial_law_invalid_entries_and_mass_rejected(bad):
    with pytest.raises(ValueError):
        c.initial_law((bad,)+tuple(Q() for _ in range(77)))


def test_missing_policy_or_record_effect_cannot_silently_drop_outcomes():
    with pytest.raises(ValueError, match='omits'):
        c.compile_tree(1, {})
    dist = c.record_law(c.source_tree(), uniform())
    with pytest.raises(ValueError, match='every public outcome'):
        c.valuation(dist, {})
    with pytest.raises(ValueError, match='unit interval'):
        c.valuation(dist, {r: 2 for r in dist})
    assert c.record_law(c.compile_tree(0, {}), uniform()) == {(): 1}
    nonedge = next((a, b) for a in range(12) for b in range(a+1, 12)
                   if (a, b) not in c.carrier_seams())
    with pytest.raises(ValueError, match='not a committed carrier seam'):
        c.compile_tree(1, {(): nonedge})


def test_diagonal_weights_do_not_select_off_diagonal_entries():
    # Two rank-one density matrices: all principal minors nonnegative, trace one.
    states = [((Q(1, 2), off), (off, Q(1, 2))) for off in (Q(1, 2), Q(-1, 2))]
    assert states[0] != states[1]
    for rho in states:
        assert rho[0][0]*rho[1][1]-rho[0][1]*rho[1][0] == 0
        assert rho[0][0]+rho[1][1] == 1
    for a, b in product((Q(), Q(1, 3), Q(1)), repeat=2):
        assert states[0][0][0]*a+states[0][1][1]*b == states[1][0][0]*a+states[1][1][1]*b


def test_quantum_diagonal_restriction_is_checked_after_classical_execution():
    result = c.quantum_diagonal_comparison()
    assert result['endpoint_coin_branches_checked'] == 4680
    assert result['exact_diagonal_blocks_equal']
    assert not result['matrix_comparison_used_to_define_classical_weights']


def test_conditional_complex_tomography_is_retained_without_probability_promotion():
    result = c.tomography_attachment()
    assert result['fresh_exact_design']['real_hermitian_rank'] == 4356
    assert result['coherent_code_services_and_quantum_semantics_supplied']
    assert not result['laboratory_data_or_universal_operational_additivity']
