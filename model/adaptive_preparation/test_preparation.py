from fractions import Fraction as F

import pytest

from . import certificate as c
from model.completion import adaptive_reference
from tools import strict_json


@pytest.fixture(scope='module')
def kernel():
    return c.native_counts()


def test_receipt_replays_and_bound_covers_transient_inputs():
    result = c.report()
    assert result == strict_json.load(c.HERE / 'results.json')
    assert result['kernel']['states'] == 132
    assert result['kernel']['stationary_support'] == 129
    assert result['complete_record_basis_inputs_checked'] == 264


def test_independent_event_executor_reproduces_complete_population(kernel):
    states, counts = kernel
    # Execute all basis inputs together under distinct rational weights,
    # then compare their accumulated exact image with the native event path.
    initial = {state: F(i + 1, sum(range(1, len(states) + 1))) for i, state in enumerate(states)}
    actual = adaptive_reference.advance(initial)
    expected = {state: sum(initial[s] * F(counts[i][j], 120) for i, s in enumerate(states))
                for j, state in enumerate(states)}
    assert actual == {s: p for s, p in expected.items() if p}


def test_periodic_and_disconnected_controls_do_not_certify_mixing():
    for matrix in ([[0, 1], [1, 0]], [[1, 0], [0, 1]]):
        with pytest.raises(ValueError, match='minorization'):
            c.block_minorization(matrix, 1, 64, F(1, 10))


def test_false_minorization_rejected(kernel):
    _, counts = kernel
    with pytest.raises(ValueError, match='minorization'):
        c.block_minorization(counts, 120, 1, F(1, 1000000))


def test_budget_uses_exact_one_sided_error_and_extra_full_record_attempt():
    for target in (F(1, 1000), F(1, 10**6), F(1, 10**9)):
        row = c.budget(target)
        bound = F(row['total_variation_bound'])
        assert bound <= target < bound / F(3, 5)
        assert row['full_last_record_attempts'] == row['source_bit_attempts'] + 1
    assert c.budget(F(1, 10), steps=2, minorization=F(1))['source_bit_attempts'] == 2


@pytest.mark.parametrize('bad', (F(0), F(1), 0.01, True))
def test_invalid_accuracy_rejected(bad):
    with pytest.raises(ValueError):
        c.budget(bad)


def test_negative_or_unnormalized_counts_rejected():
    for matrix in ([[2, -1], [0, 1]], [[1, 1], [0, 1]], [[True, 0], [0, 1]]):
        with pytest.raises(ValueError, match='stochastic counts'):
            c.block_minorization(matrix, 1, 2, F(1, 2))


def test_exact_two_state_chain_has_known_common_measure():
    power, minima, mass = c.block_minorization([[3, 1], [1, 3]], 4, 2, F(3, 4))
    assert power.tolist() == [[10, 6], [6, 10]]
    assert minima == [6, 6]
    assert mass == F(3, 4)


def test_mutated_stationary_reference_fails_native_residual(kernel, monkeypatch):
    states, counts = kernel
    data = strict_json.load(c.ROOT / 'model/completion/stationary_reference_results.json')
    rows = data['stationary_joint']
    pairs = [(i, j) for i in range(len(rows)) for j in range(i + 1, len(rows))
             if rows[i]['probability'] != rows[j]['probability']]
    i, j = pairs[0]
    rows[i]['probability'], rows[j]['probability'] = rows[j]['probability'], rows[i]['probability']
    monkeypatch.setattr(c.strict_json, 'load', lambda path: data)
    with pytest.raises(ValueError, match='native kernel residual'):
        c.reference(states, counts)
