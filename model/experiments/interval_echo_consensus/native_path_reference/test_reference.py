"""Regression tests for complete native records, the quotient and exact sampling bounds."""

from collections import Counter
from fractions import Fraction
import json
from math import log

import numpy as np
import pytest
import sympy as sp

from . import measure


@pytest.fixture(scope='module')
def receipt():
    return measure.run()


def test_complete_replay_and_retained_marker_tapes(receipt):
    measure.compare(receipt, json.loads((measure.HERE / 'results.json').read_text()))
    for panel in receipt['levels']:
        depth = panel['depth']
        raw = panel['source_records']
        assert len({tuple(row['tape']) for row in raw}) == len(raw) == 4**depth
        counts = Counter(tuple(row['output']) for row in raw)
        for row in panel['reference']:
            assert Fraction(row['probability']) == Fraction(counts[tuple(row['record'])], len(raw))
        for row in raw:
            for stage in row['records']:
                assert sum(stage['complete_flags']) == 1
                assert stage['complete_flags'][stage['position']]
        for row in panel['biased_native_source_records']:
            ones = sum(row['tape'])
            assert Fraction(row['tape_probability']) == Fraction(3**(2 * depth - ones), 4**(2 * depth))
            assert row['gates'].get('source_resets', 0) == ones


def test_refinement_and_cover_dimensions(receipt):
    prior = None
    for panel, path_rank, single_rank in zip(receipt['levels'], (3, 8, 21), (3, 5, 7), strict=True):
        assert panel['prefix_cover']['rank'] == path_rank
        assert panel['individual_time_cover']['rank'] == single_rank
        assert panel['path_only_cover_on_full_source']['invisible_classical_dimension'] == panel['complete_tapes'] - path_rank
        assert panel['weighted_reference_cost_nats'] == pytest.approx(log(2))
        optimum = {tuple(row['record']): Fraction(row['probability']) for row in panel['conditional_optimizer']}
        assert all(h[0] != 0 for h in optimum)
        if prior is not None:
            assert measure.push(optimum, lambda h: h[:-1]) == prior
        prior = optimum


def test_finite_exact_sampling_obstruction(receipt):
    result = receipt['exact_mean_projection']
    a, b, c = map(sp.sympify, result['optimizer'])
    assert a == c and sp.simplify(b + 2 * c) == 1
    assert a.is_rational is False
    x = sp.Symbol('x')
    assert sp.minimal_polynomial(a, x) == 7 * x**2 - 8 * x + 2
    # All m/k in [0,1] have nonzero integer polynomial numerator.
    for k in range(1, 129):
        assert all(abs(7 * m * m - 8 * m * k + 2 * k * k) >= 1 for m in range(k + 1))
    assert all(p['full_source_biased_pulse_excess_nats'] > .13 for p in receipt['levels'])


def test_dyadic_stopping_partition_and_bound(receipt):
    sampler = receipt['adaptive_sampler']
    cap = sampler['cutoff_reads']
    leaves = sampler['stopping_leaves']
    intervals = sorted(
        (int(''.join(map(str, r['tape'])), 2) << (cap - len(r['tape'])), 2**(cap - len(r['tape'])))
        for r in leaves
    )
    cursor = 0
    for lo, width in intervals:
        assert lo == cursor
        cursor += width
    assert cursor == 2**cap
    assert sum(Fraction(r['probability']) for r in leaves if r['cutoff']) == Fraction(2, 2**cap)
    assert Fraction(sampler['finite_expected_reads']) == 3 - Fraction(4, 2**cap)
    a = sp.sympify(receipt['exact_mean_projection']['optimizer'][0])
    target = sp.floor(a * 2**cap) / 2**cap
    assert sp.Rational(sampler['finite_output_law'][0]) == target
    assert 0 < 2 * (a - target) < sp.Rational(sampler['finite_total_variation_upper_bound'])
    # The same entry point also runs without a cap on a normally seeded native register.
    registers = measure.prepare(measure.federation(0))
    for _ in range(16):
        record = measure.sample_mean(registers)
        assert not record['cutoff'] and record['complete_flags'][record['terminal']]


def test_native_coherent_pair_and_changed_law(receipt):
    quantum = receipt['quantum_boundary']
    assert np.allclose(quantum['native_sharp_diagonal_probabilities'], .5)
    assert np.allclose(quantum['actual_native_half_contact_read'], [[0, 1], [1, 0]])
    assert np.allclose(quantum['contrast_zero_same_half_contact_read'], quantum['actual_native_half_contact_read'])
