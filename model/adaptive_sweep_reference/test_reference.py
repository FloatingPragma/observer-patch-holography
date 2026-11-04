"""Independent endpoint propagation and complete quantum-output controls."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from model.native_equilibration import sweeps
from model.observer import native
from model.observer.scalar_response.experiment import compare_receipt
from . import reference as ref


def _independent_advance(blocks, only_c=False):
    """Propagate exact populations event by event, without channel powers.

    Local outputs use endpoint floor/ceiling averages instead of native().
    Arrays have an arbitrary last column index, including an original-state
    index for the complete transition matrix.
    """
    index = {state: i for i, state in enumerate(ref.STATES)}
    result = np.zeros_like(blocks)
    for c in (0, 1):
        if only_c and c:
            continue
        word = ref.SEAMS if c == 0 else ref.SEAMS[::-1]
        routed = np.zeros_like(blocks)
        for d in (0, 1):
            for i, state in enumerate(ref.STATES):
                a, b = word[0]
                total = state[a]+state[b]
                for coin in (0, 1):
                    out = list(state)
                    out[a], out[b] = ((total//2, (total+1)//2) if coin == 0
                                      else ((total+1)//2, total//2))
                    current = state[a]-out[a]
                    raw = (0 if current == 0 else 1 if current == 1 else 2)
                    raw += 4*(total % 2) if c else 0
                    after = 2*(d ^ coin ^ (raw & 1))+c
                    routed[after, index[tuple(out)]] += blocks[2*c+d, i]
        for a, b in word[1:]:
            updated = np.zeros_like(routed)
            for i, state in enumerate(ref.STATES):
                total = state[a]+state[b]
                for coin in (0, 1):
                    out = list(state)
                    out[a], out[b] = ((total//2, (total+1)//2) if coin == 0
                                      else ((total+1)//2, total//2))
                    updated[:, index[tuple(out)]] += routed[:, i]
            routed = updated
        result += routed
    return result


def test_local_attachment_and_rank_one_E_history_certificate():
    result = ref.structural_certificate()
    assert result['checked_raw_maps'] == 2040
    assert result['E_distinct_pair_survivors_after_each_attempt'][0] == 66*65
    assert result['E_distinct_pair_survivors_after_each_attempt'][-1] == 0
    assert result['full_carrier_support_after_C_and_E_sweeps'] == [66, 66]


def test_observer_controller_continues_and_replays_every_raw_record():
    checkpoint = ref.Checkpoint((0, 1))
    initial_memory = checkpoint.memory
    state = ref.STATES[17]
    for event in range(150):
        before = checkpoint
        law, seam = checkpoint.selected_event()
        selector = 1 if law == 'E' else 0
        coin = (event//7+event) % 2
        out, (_, current) = native.native(state, seam, coin)
        raw = {0: 0, 1: 1, -1: 2}[current]
        raw += 4*((state[seam[0]]+state[seam[1]]) % 2) if selector else 0
        record = 'native', coin, raw
        checkpoint = checkpoint.commit(record)
        state = out
        assert checkpoint.history[:-1] == before.history
        assert checkpoint.history[-1] == record
        assert ref.Checkpoint.restore(checkpoint.export()) == checkpoint
        if checkpoint.phase:
            assert checkpoint.selected_event() == (law, ref.order(selector)[checkpoint.phase])
    assert checkpoint.phase == 0 and len(checkpoint.history) == 150
    assert checkpoint.recover_initial_control() == (initial_memory, 0)


def test_four_sweep_minorization_against_independent_elementary_propagation():
    value = np.eye(264, dtype=object).reshape(4, 66, 264)
    only_c = value.copy()
    minima = []
    for _ in range(4):
        value = _independent_advance(value)
        only_c = _independent_advance(only_c, only_c=True)
        minima.append(int((value-only_c).min()))
    measured = (value-only_c).reshape(264, 264)
    certificate = ref.minorization()
    assert minima == certificate['minimum_measured_path_counts_by_sweep']
    assert minima == [0, 0, 0, 1586432651050269168489450503667712]
    serialized = json.dumps([[int(x) for x in row] for row in measured], separators=(',', ':'))
    assert hashlib.sha256(serialized.encode()).hexdigest() == certificate['rank_one_path_count_sha256']
    assert F(3, 10) < F(certificate['depolarizing_weight']) < F(1, 3)


def test_boundary_population_uses_actual_sweep_order_and_current_feedback():
    actual = ref.population_channel()
    independent = _independent_advance(np.eye(264, dtype=object).reshape(4, 66, 264))
    np.testing.assert_array_equal(actual, independent.reshape(264, 264))
    for memory in range(4):
        selector = memory//2
        indices = tuple(range(30)) if selector == 0 else tuple(reversed(range(30)))
        expected = sweeps.classical_power(1, order=indices)
        marginal = actual[:, 66*memory:66*(memory+1)].reshape(4, 66, 66).sum(axis=0)
        np.testing.assert_array_equal(marginal, np.array(expected.counts, dtype=object))
    # Memory is initially fixed; source configurations affect its next value.
    conditional = actual[:, :66].reshape(4, 66, 66).sum(axis=1)
    assert len({tuple(conditional[:, i]) for i in range(66)}) > 1


def _entangled_input(full=False):
    d = 78 if full else 66
    basis = native.states(binary=not full)
    vector = np.zeros(2*d, dtype=complex)
    left = next(i for i, state in enumerate(basis) if max(state) == (2 if full else 1))
    right = next(i for i, state in enumerate(basis) if max(state) == 1 and i != left)
    vector[2*left], vector[2*right+1] = .6, .8j
    blocks = np.zeros((4, 2*d, 2*d), dtype=complex)
    blocks[0] = np.outer(vector, vector.conj())
    return blocks


def _reference_marginal(blocks, dimension):
    return sum(np.trace(b.reshape(dimension, 2, dimension, 2), axis1=0, axis2=2)
               for b in blocks)


def test_complete_source_execution_preserves_stationary_joint_state():
    blocks = np.array([np.eye(66)/264]*4, dtype=complex)
    actual = ref.archive_traced_tick(blocks)
    np.testing.assert_allclose(actual, blocks, rtol=0, atol=2e-17)


@pytest.mark.parametrize('full', [False, True])
def test_source_execution_includes_untouched_reference_and_absorption(full):
    blocks = _entangled_input(full)
    actual = ref.archive_traced_tick(blocks, reference_dimension=2, full=full)
    d = 78 if full else 66
    np.testing.assert_allclose(_reference_marginal(actual, d),
                               _reference_marginal(blocks, d), rtol=0, atol=2e-15)
    assert all(np.linalg.eigvalsh(b).min() > -2e-15 for b in actual)
    if full:
        bad = [2*i+j for i, s in enumerate(native.states()) if max(s) == 2 for j in (0, 1)]
        np.testing.assert_array_equal(actual[:, bad, :], 0)
        np.testing.assert_array_equal(actual[:, :, bad], 0)


def test_execution_matches_classical_integer_channel_and_retains_subnormalization():
    rng = np.random.default_rng(1138)
    population = rng.random(264)
    population /= population.sum()
    blocks = np.array([np.diag(row) for row in population.reshape(4, 66)], dtype=complex)
    actual = ref.archive_traced_tick(blocks)
    expected = np.asarray(ref.population_channel(), dtype=float)@population/ref.DENOMINATOR
    np.testing.assert_allclose(np.array([np.diag(b).real for b in actual]).ravel(), expected,
                               rtol=0, atol=5e-18)
    small = ref.archive_traced_tick(blocks/8)
    np.testing.assert_allclose(small, actual/8, rtol=0, atol=1e-18)
    np.testing.assert_array_equal(ref.archive_traced_tick(blocks*0), 0)


def test_four_sweep_quantum_contraction_on_entangled_input():
    blocks = _entangled_input()
    marginal = _reference_marginal(blocks, 66)
    stationary = np.array([np.kron(np.eye(66)/264, marginal)]*4)
    initial = sum(np.linalg.norm(b-s, 'nuc') for b, s in zip(blocks, stationary))/2
    for _ in range(4):
        blocks = ref.archive_traced_tick(blocks, reference_dimension=2)
    actual = sum(np.linalg.norm(b-s, 'nuc') for b, s in zip(blocks, stationary))/2
    assert actual <= float(ref.preparation_bound(4))*initial+2e-14
    assert actual < initial


def test_exact_budget_covers_full_input_and_has_no_uniform_phase_claim():
    binary = ref.preparation_budget(F(1, 10**6))
    full = ref.preparation_budget(F(1, 10**6), full=True)
    assert full['sweeps'] == binary['sweeps']+1
    assert F(full['half_trace_bound']) <= F(1, 10**6)
    assert ref.preparation_bound(full['sweeps']-1, full=True) > F(1, 10**6)
    assert ref.minorization()['arbitrary_phase_inputs_mix_to_uniform_phase'] is False
    assert ref.minorization()['stationary_reference_modular_automorphism'] == 'identity'


@pytest.mark.parametrize('value', [True, -1, 1.0, F(1), np.int64(1), 100001])
def test_invalid_sweep_counts_fail_closed(value):
    with pytest.raises(ValueError):
        ref.preparation_bound(value)


@pytest.mark.parametrize('value', [True, 0, -1, .1, float('nan')])
def test_invalid_tolerance_fails_closed(value):
    with pytest.raises(ValueError):
        ref.preparation_budget(value)


def test_invalid_quantum_inputs_fail_closed():
    for bad in (np.full((4, 66, 66), np.nan), np.zeros((4, 66, 66), dtype=int),
                np.array([-np.eye(66)/264]*4), np.array([np.eye(66)]*4)):
        with pytest.raises(ValueError):
            ref.archive_traced_tick(bad)


def test_retained_certificate_replays_and_detects_forged_bound():
    saved = json.loads(Path(ref.OUTPUT).read_text())
    actual = ref.run()
    compare_receipt(saved, actual)
    actual['full_carrier_mixing'] = dict(actual['full_carrier_mixing'], depolarizing_weight='1')
    with pytest.raises((AssertionError, ValueError)):
        compare_receipt(saved, actual)
