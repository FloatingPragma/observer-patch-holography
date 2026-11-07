"""Exact operator, complete-record and counterexample checks."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F

import pytest

from . import causality as c, packet
from ..observer import native


@pytest.fixture(scope='module')
def receipt():
    return packet.produce()


def test_exhaustive_full_carrier_operator_identity(receipt):
    rows = receipt['exact_operator_census']['full_carriers']
    assert [row['dimension'] for row in rows] == [66, 78]
    assert [row['raw_operator_pairs'] for row in rows] == [196875, 196875]
    assert [row['zero_composed_operators_retained'] for row in rows] == [158760, 139860]


@pytest.mark.parametrize('shell', ('H66', 'H78'))
def test_raw_operators_independently_match_native_endpoint_equation(shell):
    states = native.states(shell == 'H66')
    index = {value: i for i, value in enumerate(states)}
    for edge in c.source.local.SEAMS:
        slots = c.raw_slots(shell, edge)
        assert slots[0].pairs == tuple((i, i) for i in range(len(states)))
        for slot in slots[1:]:
            law, coin, raw = slot.label
            label = c.source.circuit.record_values(law, coin).get(raw)
            expected = []
            for i, loads in enumerate(states):
                after, (_, current) = native.native(loads, edge, coin)
                measured = current if law == 'C' else loads[edge[0]]-loads[edge[1]]
                if measured == label:
                    expected.append((i, index[after]))
            assert slot.pairs == tuple(expected)


def test_entangled_reference_and_exported_records(receipt):
    rows = receipt['actual_interpreter_comparisons']
    assert len(rows) == 6
    for row in rows:
        assert row['complete_blocks'] == 625
        assert row['zero_blocks_retained'] == 576
        assert row['same_after_reversible_record_permutation']
        assert not row['raw_chronological_outputs_equal']
        assert row['archive_slots_per_output'] == max(0, 2-row['FIFO_capacity'])
    assert all(row['zero_blocks_retained'] == 624
               for row in receipt['complex_linearity_controls'])


@pytest.mark.parametrize('law', ('C', 'E'))
def test_fixed_record_law_keeps_zero_weight_other_law(law):
    row = c.complete_pair('H78', 1, law=law)
    assert row['complete_blocks'] == 625
    assert row['zero_blocks_retained'] > 576


def test_complete_tape_permutation_is_invertible_across_fifo_boundary():
    tape = tuple(range(7));order = (3, 1, 4, 0, 6, 2, 5)
    inverse = tuple(order.index(i) for i in range(len(order)))
    for capacity in (1, 2, 7, 8):
        archive, memory = c.permute_tape(tape, capacity, order)
        assert len(memory) == min(capacity, len(tape))
        assert c.permute_tape(archive+memory, capacity, inverse) == c.split_tape(tape, capacity)
    for bad in ((0,), (0,)*7, tuple(range(6)), (True, 1, 2, 3, 4, 5, 6)):
        with pytest.raises(ValueError):
            c.permute_tape(tape, 1, bad)


@pytest.fixture(scope='module')
def history():
    value, _ = c.paired_programs()
    x = c._index('H66', (0, 1))
    return c.source.execute(value, {(x, 0, x, 0): F(1)})


def test_comparison_rejects_stops_missing_zero_slots_and_duplicate_records(history):
    zero = next(i for i, row in enumerate(history.branches) if not row.density)
    missing = history.branches[:zero]+history.branches[zero+1:]
    duplicate = history.branches[:-1]+(history.branches[0],)
    bads = (replace(history, status='RESOURCE_STOP'),
            replace(history, attempted_steps=1), replace(history, branches=missing),
            replace(history, branches=duplicate))
    for bad in bads:
        with pytest.raises((ValueError, ArithmeticError)):
            c.comparison_blocks(bad, 1)
    first = history.branches[0]
    with pytest.raises((ValueError, ArithmeticError)):
        c.comparison_blocks(replace(history, branches=(replace(first, archive=()),)+history.branches[1:]), 1)


def test_literal_record_comparison_without_permutation_fails(history):
    _, value = c.paired_programs()
    x = c._index('H66', (0, 1))
    reverse = c.source.execute(value, {(x, 0, x, 0): F(1)})
    assert c.comparison_blocks(history, 1) != c.comparison_blocks(reverse, 1)
    assert c.comparison_blocks(history, 1) == c.comparison_blocks(reverse, 1, swap=True)
    assert history.program_digest != reverse.program_digest


def test_feedback_timeout_and_overlap_boundaries(receipt):
    value = receipt['boundaries']
    assert value['adaptive_controller']['receiver_probabilities'] == ['1/4', '0']
    assert value['intrinsic_FIFO_timeout']['receiver_probabilities'] == ['8/33', '0']
    assert value['overlapping_fixed_seams']['receiver_probabilities'] == ['8/33', '4/33']


@pytest.mark.parametrize('kwargs', ({'shell':'H6'}, {'capacity':True},
    {'edges':((0, 5), (0, 1))}, {'edges':((5, 0), (1, 7))},
    {'law':'adaptive'}, {'repair_probability':0.5}))
def test_rejects_outside_fixed_exact_domain(kwargs):
    with pytest.raises(ValueError):
        c.paired_programs(**kwargs)


def test_receipt_binds_own_proof_and_transitive_actual_source(receipt):
    bindings = receipt['source_bindings']
    for name in (*packet.SOURCE_FILES, 'model/native_program/interpreter.py',
                 'model/observer/native.py', 'model/local_repair_family/local_bias.py'):
        assert name in bindings
    assert all(not name.startswith('docs/') for name in bindings)
    assert bindings == packet.source_bindings()


def test_receipt_verifier_rejects_changed_census_scope_and_source(receipt, tmp_path, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: receipt)
    path = tmp_path/'receipt.json'
    path.write_text(packet.canonical(receipt));assert packet.verify(path)
    changes = [deepcopy(receipt) for _ in range(3)]
    changes[0]['exact_operator_census']['full_carriers'][0]['raw_operator_pairs'] -= 1
    changes[1]['scope']['three_spatial_dimensions_or_Lorentz_limit_derived'] = True
    name = next(iter(changes[2]['source_bindings']))
    changes[2]['source_bindings'][name] = '0'*64
    for change in changes:
        path.write_text(packet.canonical(change))
        with pytest.raises(ValueError):
            packet.verify(path)
    path.write_text('{"schema": 1, "schema": 2}')
    with pytest.raises(ValueError):
        packet.verify(path)
