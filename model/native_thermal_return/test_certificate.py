from fractions import Fraction as Q
from itertools import repeat
import random
from unittest.mock import patch

import pytest

from . import certificate as c


def streams(seed):
    coin_rng, seam_rng = random.Random(seed), random.Random(seed+10000)
    return (coin_rng.randrange(2) for _ in repeat(None)), (seam_rng.randrange(30) for _ in repeat(None))


def test_two_census_source_tables_and_attempt_budgets():
    first, last = c.success_tables()
    assert len(first) == len(last) == 66
    assert 1-min(first) <= Q(35, 2**31)
    assert 1-min(last) <= Q(33, 2**31)
    assert c.bounds(warmup=512, absorb=True)['native_attempts_max'] == 5774
    assert c.bounds(clock='attempt', warmup=512, absorb=True)['native_attempts_max'] == 5712


def test_complete_joint_endpoint_body_replayed_independently():
    assert c.thermal.loop_law(cap=2) == c.thermal.native_oracle(cap=2)


def test_all_census_and_waiting_failures_are_present_before_conditioning():
    law = c.record_law(seed_scans=1, retries=2)
    assert sum(law.values()) == 1
    for flag in ('initial_census_failure', 'terminal_census_failure', 'waiting_timeout', 'return', 'different'):
        assert sum(p for (kind, _), p in law.items() if kind == flag) > 0
    assert sum(c.conditioned(law).values()) == 1


def test_uniform_failure_bound_covers_nonuniform_input_and_zero_support():
    first, last = c.success_tables(1, 3)
    limit = 2-min(first)-min(last)
    for p in (tuple(Q(i == 0) for i in range(66)), tuple(Q(i+1, 2211) for i in range(66))):
        actual = c.record_law(p, seed_scans=1, retries=3)
        ideal = c.ideal_census_law(p)
        assert c.tv(actual, ideal) <= limit


def test_finite_census_success_is_not_an_unbiased_free_measurement():
    actual = c.record_law(seed_scans=1, retries=2)
    ideal = c.ideal_census_law()
    assert c.tv(c.conditioned(actual), c.conditioned(ideal)) > 0


def test_mass_loss_in_a_native_body_column_is_rejected():
    body = c.thermal.loop_law(cap=1)
    matrices = [[list(row) for row in matrix] for matrix in body.counts]
    matrices[0][0][0] += 1
    changed = body._replace(counts=tuple(tuple(map(tuple, matrix)) for matrix in matrices))
    with pytest.raises(ValueError, match='lost an input column'):
        c._compose(changed, c.population(), *c.success_tables())


@pytest.mark.parametrize('law', ['C', 'E'])
def test_executed_native_censuses_restore_then_retrodict_without_initial_disclosure(law):
    for i in (0, 12, 37, 65):
        coins, seams = streams(400+i)
        result = c.execute(c.BASIS[i], coins, seams, cap=32, law=law)
        assert result.initial_census.outcome.initial == c.BASIS[i]
        assert result.initial_census.output == c.BASIS[i]
        assert result.initial_census.outcome.success
        assert result.flag in ('return', 'different')
        body_endpoint = list(c.BASIS[i])
        fieldword = []
        for edge, _, current, _ in result.body_records:
            body_endpoint[edge[0]] -= current
            body_endpoint[edge[1]] += current
            if set(edge).intersection(c.core.source().region):
                fieldword.append(body_endpoint[(0, 1)[len(fieldword)]])
        assert result.word == tuple(fieldword)
        assert result.terminal_census.outcome.initial == tuple(body_endpoint)
        assert result.flag == ('return' if tuple(body_endpoint) == c.BASIS[i] else 'different')
        assert result.attempts <= c.bounds()['native_attempts_max']


def test_wait_timeout_is_not_a_return_and_does_not_execute_terminal_census():
    coins, _ = streams(123)
    remote = next(i for i, edge in enumerate(c.census.quantum.seams())
                  if not set(edge).intersection(c.core.source().region))
    result = c.execute(c.BASIS[0], coins, repeat(remote), cap=1)
    assert result.flag == 'waiting_timeout'
    assert result.word == () and result.terminal_census is None
    assert len(result.body_records) == 1


def test_full_startup_uses_a_retained_covering_sweep_before_restoration():
    initial = (2,)+(0,)*11
    coins, seams = streams(87)
    result = c.execute(initial, coins, seams, absorb=True, warmup=3, clock='attempt')
    assert len(result.preparation_records) == 33
    assert result.initial_census.outcome.success
    assert max(result.initial_census.outcome.initial) == 1
    with pytest.raises(ValueError, match='configuration'):
        c.execute(initial, *streams(87))


def test_attempt_clock_full_words_and_trace_probability():
    actual, ideal = c.attempt_record_law(length=1), c.attempt_record_law(length=1, ideal=True)
    assert len(actual) == 3*66+1 and sum(actual.values()) == 1
    expected = sum(dict(row).get(i, 0) for i, row in enumerate(c.attempt_rows()))
    assert sum(p for (flag, _), p in ideal.items() if flag == 'return') == Q(expected, 66*60)
    assert not any(flag == 'waiting_timeout' for flag, _ in actual)
    result = c.execute(c.BASIS[11], *streams(99), clock='attempt', length=2)
    assert result.flag in ('return', 'different') and len(result.word) == 2
    state = list(result.initial_census.outcome.initial)
    word = []
    for edge, _, current, _ in result.body_records:
        state[edge[0]] -= current; state[edge[1]] += current
        word.append(c.BASIS.index(tuple(state)))
    assert result.word == tuple(word)


def test_conditional_bounds_charge_rare_returns_and_degenerate_denominators():
    loose = c.bounds(warmup=0)
    assert loose['conditional_return_error_upper'] == 1
    tight = c.bounds(clock='attempt', warmup=512, absorb=True)
    assert tight['waiting_error_upper'] == 0
    assert 0 < tight['conditional_return_error_upper'] < Q(1, 1000)
    assert c.bounds(warmup=None)['exact_uniform_preparation_assumed']
    assert not tight['exact_uniform_preparation_assumed']


def test_thermal_energy_is_reconstructed_before_the_loop_length():
    spectral = c.core.spectral_certificate()
    assert spectral['effective_transfer_lower'] == Q(3, 38)
    assert spectral['effective_transfer_upper'] == 1
    assert c.thermal.kms_defect_lower(Q(1, 2)) == Q(25, 1188)
    assert c.thermal.local_closure_control()


@pytest.mark.parametrize('bad', [True, 2.0, -1])
def test_exact_public_protocol_parameters_reject_coercions(bad):
    with pytest.raises(ValueError):
        c.bounds(length=bad)
    with pytest.raises(ValueError):
        c.attempt_record_law(length=bad)


def test_frozen_source_pin_change_is_not_silently_resealed():
    old = c.strict_json.load
    def corrupted(path):
        result = old(path)
        if str(path).endswith('configuration_results.json'):
            result['source_sha256']['model/native_clock/configuration.py'] = '0'*64
        return result
    with patch.object(c.strict_json, 'load', side_effect=corrupted):
        with pytest.raises(ValueError, match='inherited source binding changed'):
            c.inherited_pins()
