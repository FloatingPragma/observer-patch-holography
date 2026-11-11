from fractions import Fraction as F

import pytest

from .record_observer import Observer, encode, memory_bound, mixing_horizon, report, tables


def test_every_native_branch_is_reconstructed_and_exact_relaxation_is_short():
    result = report()
    assert result['exhaustive_native_branch_checks'] == 3960
    assert result['dobrushin_one_attempt_contraction'] == '14/15'
    horizon = result['one_percent_relaxation']
    assert horizon['attempts'] == 39
    assert F(horizon['worst_TV_exact']) <= F(1, 100)
    assert F(horizon['previous_worst_TV_exact']) > F(1, 100)
    assert result['memory_at_one_percent_horizon']['retained_register_bits_upper'] == 298
    assert result['memory_at_one_percent_horizon']['total_including_incoming_packet_bits_upper'] == 311


def test_observers_agree_after_complete_native_histories_without_hidden_load_reads():
    from model.observer import native
    states, seams, distances = tables()
    for initial in range(66):
        source = states[initial]
        alice = Observer.prepare(initial, 39)
        bob = Observer.prepare(initial, 39)
        for sequence in range(39):
            address, coin = (7 * sequence + initial) % 30, sequence % 2
            source, (_, current) = native.native(source, seams[address], coin)
            packet = encode(address, coin, current)
            alice = alice.consume(sequence, packet)
            bob = bob.consume(sequence, packet)
            assert alice == bob
            assert states[alice.state_index] == source
            assert alice.distance() == distances[states.index(source)]
        assert alice.status == 'full' and len(alice.tape) == 39
        assert alice.consume(39, 0) == alice


def test_missing_or_duplicate_record_halts_inference():
    observer = Observer.prepare(0, 10)
    skipped = observer.consume(1, 0)
    assert skipped.status == 'sequence_error' and skipped.distance() is None
    assert skipped.consume(0, 0) == skipped


def test_impossible_native_record_halts_without_inventing_distance():
    states, seams, _ = tables()
    a, b = seams[0]
    initial = next(i for i, state in enumerate(states) if state[a] == state[b] == 0)
    observer = Observer.prepare(initial, 5).consume(0, encode(0, 0, 1))
    assert observer.status == 'invalid_record'
    assert observer.distance() is None
    assert observer.tape == ()


def test_no_operation_is_retained_and_counts_against_capacity():
    states, seams, _ = tables()
    a, b = seams[0]
    initial = next(i for i, state in enumerate(states) if state[a] == state[b])
    observer = Observer.prepare(initial, 1).consume(0, encode(0, 0, 0))
    assert observer.state_index == initial
    assert observer.tape == (encode(0, 0, 0),) and observer.status == 'full'


@pytest.mark.parametrize('letter', [-1, 120, 128, 1.0, True])
def test_invalid_record_alphabet_is_not_silently_truncated(letter):
    assert Observer.prepare(0, 1).consume(0, letter).status == 'invalid_record'


def test_wrong_initial_preparation_is_not_discovered_by_a_noop_record():
    states, seams, distances = tables()
    from model.observer import native
    rows = {}
    for i, state in enumerate(states):
        _, record = native.native(state, seams[0], 0)
        rows.setdefault(record, []).append(i)
    pair = next((a, b, record) for record, ids in rows.items() for a in ids for b in ids
                if distances[a] != distances[b])
    a, b, (coin, current) = pair
    first = Observer.prepare(a, 1).consume(0, encode(0, coin, current))
    second = Observer.prepare(b, 1).consume(0, encode(0, coin, current))
    assert first.status == second.status == 'full'
    assert first.distance() != second.distance()


def test_zero_capacity_and_uncertified_relaxation_are_explicit():
    assert Observer.prepare(0, 0).status == 'full'
    assert memory_bound(0)['retained_native_tape_bits'] == 0
    with pytest.raises(ValueError, match='not certified'):
        mixing_horizon([[1, 0], [0, 1]], [1, 1], F(1, 100), 10)
