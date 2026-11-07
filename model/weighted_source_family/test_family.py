"""Independent native attachment, weighted restriction and complete-record checks."""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
from math import prod

import numpy as np
import pytest
import sympy as s

from . import family as f
from ..native_capacity_restriction import restriction as fair


def path(weights):
    return f.Family(weights, tuple((i, i+1) for i in range(len(weights)-1)),
                    tuple(i+1 for i in range(len(weights)-1)))


def masks_in(mask, indices):
    return sum(((mask >> i) & 1) << j for j, i in enumerate(indices))


@pytest.mark.parametrize('weights', [(1, 1), (2, 3, 5), (25, 8, 25, 8), (1, 4, 2, 3, 2)])
def test_all_number_sectors_native_detailed_balance_and_unique_population(weights):
    family = path(weights)
    for m in range(len(weights)+1):
        pi = f.canonical(weights, m)
        columns = {x: f.population_column(family, x) for x in pi}
        assert all(sum(row.values()) == 1 for row in columns.values())
        for x, y in product(pi, repeat=2):
            assert pi[x]*columns[x].get(y, 0) == pi[y]*columns[y].get(x, 0)
        visited, pending = {next(iter(pi))}, [next(iter(pi))]
        while pending:
            for y in columns[pending.pop()]:
                if y not in visited:
                    visited.add(y); pending.append(y)
        assert visited == set(pi)
        density = {(x, 0, x, 0): p for x, p in pi.items()}
        total = defaultdict(F)
        for _, output in f.attempt(family, density, c_weight=F(2, 5), repair=F(3, 4)):
            for key, value in output.items():
                total[key] += value
        assert dict(total) == density


@pytest.mark.parametrize('law', ['C', 'E', 'mixture'])
def test_literal_H66_attachment_including_all_zero_raw_slots(law):
    family = f.local_family()
    states = f.local.basis()
    masks = [sum(value << i for i, value in enumerate(x)) for x in states]
    marker = {(i, i, i, i): F(1, 66) for i in range(66)}
    initial = {(masks[i], i, masks[i], i): value for (i, _, _, _), value in marker.items()}
    command = {'op':'repair', 'selector':{'kind':'uniform'}, 'law':law, 'repeat':1}
    native = f.interpreter.program(commands=[command])
    expected = f.interpreter._outcomes(native, command, marker)
    actual = f.attempt(family, initial, law=law)
    assert len(expected) == len(actual) == 721
    for outcome, (record, density) in zip(expected, actual, strict=True):
        translated = {(masks[i], a, masks[j], b): value
                      for (i, a, j, b), value in outcome.density.items()}
        assert record == outcome.record
        assert density == translated
    # Each positive branch is a scalar partial injection. The marker comparison
    # fixes every domain/image and therefore every complex matrix-unit action.
    assert all(len({i for i, _, _, _ in d}) == len(d) for _, d in actual)


@pytest.mark.parametrize('weights', [(1, 1, 1, 1), (25, 8, 25, 8), (2, 3, 4, 2, 5)])
def test_every_window_full_mixture_and_nested_actual_partial_trace(weights):
    n = len(weights)
    for m in range(n+1):
        full = f.canonical(weights, m)
        for bits in range(1 << n):
            window = tuple(i for i in range(n) if bits >> i & 1)
            traced = defaultdict(F)
            for mask, value in full.items():
                traced[masks_in(mask, window)] += value
            exact = f.restricted_state(weights, m, window)
            assert exact == dict(traced)
            assert sum(f.window_law(weights, m, window)) == 1
            inner = tuple(i for i in range(len(window)) if i % 2 == 0)
            nested = defaultdict(F)
            for mask, value in exact.items():
                nested[masks_in(mask, inner)] += value
            assert dict(nested) == f.restricted_state(weights, m, tuple(window[i] for i in inner))


@pytest.mark.parametrize('n,m,v', [(4, 2, 2), (8, 3, 4), (12, 2, 12), (12, 0, 4), (12, 12, 3)])
def test_equal_weight_reduction_reuses_fair_owner(n, m, v):
    assert f.window_law((7,)*n, m, tuple(range(v))) == fair.number_law(n, m, v)


@pytest.mark.parametrize('weights', [(2, 3, 5), (25, 25)+(8,)*4])
@pytest.mark.parametrize('z', [F(1, 8), F(2)])
def test_fugacity_mixture_is_product_and_exactly_restriction_compatible(weights, z):
    law = f.grand_canonical(weights, z)
    assert sum(law) == 1
    mixture = {mask: p*law[m] for m in range(len(weights)+1)
               for mask, p in f.canonical(weights, m).items()}
    direct = {mask: prod((z*w if mask >> i & 1 else 1)/(1+z*w)
                        for i, w in enumerate(weights)) for mask in range(1 << len(weights))}
    assert mixture == direct
    chosen = (0, len(weights)-1)
    traced = defaultdict(F)
    for mask, p in direct.items():
        traced[masks_in(mask, chosen)] += p
    selected = tuple(weights[i] for i in chosen)
    assert dict(traced) == {mask: prod((z*w if mask >> i & 1 else 1)/(1+z*w)
                                     for i, w in enumerate(selected)) for mask in range(4)}


@pytest.mark.parametrize('first,second', [(1, 1), (25, 8), (3, 7)])
def test_local_generator_is_exact_actual_symmetrization_and_positive_projector(first, second):
    family = path((first, second))
    states = (0, 1, 2, 3)
    column = s.Matrix([[f.population_column(family, x).get(y, 0) for x in states] for y in states])
    omega = s.diag(1, s.sqrt(first), s.sqrt(second), s.sqrt(first*second))
    H = f.local_generator(first, second)
    assert s.simplify(H-(s.eye(4)-omega.inv()*column*omega)) == s.zeros(4)
    assert s.simplify(H*H-H) == s.zeros(4)
    assert s.simplify(H*omega*s.ones(4, 1)) == s.zeros(4, 1)
    assert H == H.T


@pytest.mark.parametrize('weights', [(25, 8), (25, 25, 8, 8), (2, 3, 5)])
def test_native_modular_gap_cancels_capacity_occupation_and_fugacity(weights):
    family = path(weights)
    z = F(2, 7)
    product_law = {mask: prod((z*w if mask >> i & 1 else 1)/(1+z*w)
                             for i, w in enumerate(weights)) for mask in range(1 << len(weights))}
    for m in range(len(weights)+1):
        pi = f.canonical(weights, m)
        for mask, edge, coin in product(pi, family.seams, (0, 1)):
            target, current, _ = f.image(mask, edge, coin, len(weights))
            assert pi[target]/pi[mask] == f.modular_ratio(family, edge, current)
            assert product_law[target]/product_law[mask] == pi[target]/pi[mask]


@pytest.mark.parametrize('z', [F(1, 8), F(3, 5)])
def test_product_modular_covariance_on_all_partial_trace_matrix_units(z):
    weights = (25, 25, 8, 8)
    global_pi = {mask: prod((z*w if mask >> i & 1 else 1)/(1+z*w)
                           for i, w in enumerate(weights)) for mask in range(16)}
    for subset in range(16):
        inside = tuple(i for i in range(4) if subset >> i & 1)
        outside = tuple(i for i in range(4) if not subset >> i & 1)
        local = {mask: prod((z*weights[i] if mask >> j & 1 else 1)/(1+z*weights[i])
                           for j, i in enumerate(inside)) for mask in range(1 << len(inside))}
        for left, right in product(range(16), repeat=2):
            if masks_in(left, outside) == masks_in(right, outside):
                assert global_pi[left]/global_pi[right] == local[masks_in(left, inside)]/local[masks_in(right, inside)]
            # Other units have zero partial trace on both sides because the
            # diagonal modular action preserves their distinct exterior labels.


@pytest.mark.parametrize('level', [0, 1, 2])
def test_existing_source_graph_and_conditional_H66_reference(level):
    family = f.federation_family(level)
    graph = f.geometry.federation(level)
    assert len(family.weights) == 240*4**level
    assert len(family.seams) == 630*4**level
    assert Counter(i//12 for edge in graph.inter_seams for i in edge) == {
        cell:3 for cell in range(graph.cells)}
    cell = tuple(range(12))
    assert len(f.boundary(family, cell)) == 3
    assert len([e for e in family.seams if set(e) <= set(cell)]) == 30
    law = f.window_law(family.weights, 2, cell)
    assert 0 < law[2] < 1
    conditional = f.canonical(family.weights[:12], 2)
    target = {sum(value << i for i, value in enumerate(x)): p
              for x, p in zip(f.local.basis(), f.local.reference(f.local.Parameters(
                  port_weights=(25, 25)+(8,)*10)), strict=True)}
    assert conditional == target
    bounds = f.cut_bounds(family, cell, attempts=4, internal_repairs=3)
    assert bounds['boundary_probability'] == F(3, len(family.seams))
    assert bounds['no_boundary_through_kth_internal'] == F(10, 11)**3
    assert graph.inter_seams and all(family.coin(e, 0)+family.coin(e, 1) == 1 for e in graph.inter_seams)


def test_declared_weighted_extension_is_rejected_by_unchanged_interpreter():
    native = f.interpreter.program(federation={'level':0, 'total':2})
    native['parameters']['port_weights'] = [25, 25]+[8]*10
    native['declared_resources'] = f.interpreter._native_cost(native)
    with pytest.raises(ValueError, match='not an invented bias lift'):
        f.interpreter.admit(native)


@pytest.mark.parametrize('law', ['C', 'E', 'mixture'])
def test_complete_boundary_control_bound_with_entangled_reference(law):
    family = path((25, 8, 3))
    # All eight configurations are correlated with an untouched reference.
    density = {(i, i, j, j): F(1, 8) for i, j in product(range(8), repeat=2)}
    actual = f.attempt(family, density, law=law)
    control = f.attempt(family, density, law=law, cut=(0, 1))
    assert tuple(r for r, _ in actual) == tuple(r for r, _ in control)
    distance = 0.
    reference = [defaultdict(F), defaultdict(F)]
    for (record, a), (_, b) in zip(actual, control, strict=True):
        difference = np.zeros((64, 64))
        for sign, d in ((1, a), (-1, b)):
            for (i, r, j, q), value in d.items():
                difference[8*i+r, 8*j+q] += sign*float(value)
        distance += np.abs(np.linalg.eigvalsh(difference)).sum()/2
        for k, d in enumerate((a, b)):
            for (i, r, j, q), value in d.items():
                if i == j:
                    reference[k][r, q] += value
        if record.kind == 'repair' and record.seam != (1, 2):
            assert a == b
    assert reference[0] == reference[1] == {(i, i): F(1, 8) for i in range(8)}
    assert 0 < distance <= float(f.cut_bounds(family, (0, 1), attempts=1,
                                             internal_repairs=1)['global_attempt_error_upper'])+1e-14


def test_independent_clock_event_count_and_wait_summation():
    family = path((2, 3, 4))
    bound = f.cut_bounds(family, (0, 1), attempts=4, internal_repairs=3)
    # There is no exterior edge in this independent enumerator.
    counts = defaultdict(F)
    for choices in product(range(2), repeat=4):
        p = prod(family.selection(family.seams[i]) for i in choices)
        counts[any(i == 1 for i in choices)] += p
    assert counts[True] == bound['global_attempt_error_upper']
    assert bound['no_boundary_through_kth_internal'] == F(1, 3)**3
    # Sum geometric exterior waits exactly: relevant-event probability a/(a+b).
    a, b, exterior = F(30, 630), F(3, 630), F(597, 630)
    assert a/(1-exterior) == a/(a+b) == F(10, 11)


@pytest.mark.parametrize('mutation', ['zero_weight', 'disconnected', 'unordered', 'list_edge'])
def test_invalid_families_rejected(mutation):
    args = {'weights':(2, 3, 4), 'seams':((0, 1), (1, 2)), 'seam_weights':(1, 1)}
    if mutation == 'zero_weight': args['weights'] = (2, 0, 4)
    if mutation == 'disconnected': args.update(seams=((0, 1),), seam_weights=(1,))
    if mutation == 'unordered': args['seams'] = ((1, 2), (0, 1))
    if mutation == 'list_edge': args['seams'] = ([0, 1], [1, 2])
    with pytest.raises(ValueError): f.Family(**args)


@pytest.mark.parametrize('value', [False, 0, (), {}])
def test_explicit_invalid_local_parameters_do_not_select_default(value):
    with pytest.raises(ValueError, match='Parameters required'):
        f.local_family(value)


@pytest.mark.parametrize('edge,coin', [((0, 1), 2), ((0, 1), True), ((1, 0), 0), ((0, 2), 0)])
def test_native_image_domain_rejects_undeclared_addresses_or_coins(edge, coin):
    with pytest.raises(ValueError):
        f.image(1, edge, coin, 2)
