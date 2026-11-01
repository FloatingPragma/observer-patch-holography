from fractions import Fraction as F
from itertools import product
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy import sparse

from ..charged_defect import dictionary
from . import controller as c, history, instrument as ins, packet


def test_recordwise_memory_inverse_and_incoming_selector():
    for memory, record in product(product(range(2), repeat=2), c.LABELS):
        assert c.retreat(c.advance(memory, record), record) == memory
        assert c.advance(c.retreat(memory, record), record) == memory
        assert c.retreat(memory, record)[0] == memory[1]
    assert c.permutation_certificate()['all_68_memory_record_pairs_checked']


def test_native_balanced_seed_and_unchanged_spectator():
    seed = c.seed_certificate()
    assert sum(row['probability'] for row in seed['uniform_memory']) == 1
    assert all(row['probability'] == F(1, 4) for row in seed['uniform_memory'])
    assert len({tuple(row['memory']) for row in seed['uniform_memory']}) == 4
    assert seed['full_two_event_raw_archive_slots'] == 64


@pytest.mark.parametrize('memory', [(), (0,), (0, 1, 0), [0, 1], (True, 0), (0, 2)])
def test_invalid_memory(memory):
    with pytest.raises(ValueError): c.advance(memory, ('outside',))


@pytest.mark.parametrize('record', [('native', True, 0), ('native', 1, 8), ('native', 2, 0),
                                     ('outside', 0), ['outside'], ('unknown',)])
def test_invalid_record(record):
    with pytest.raises(ValueError): c.advance((0, 0), record)


def test_checkpoint_roundtrip_and_explicit_feedback():
    before = c.Checkpoint((1, 0))
    a = before.commit(('native', 0, 4))
    b = before.commit(('native', 0, 5))
    assert a.memory[0] != b.memory[0]
    assert c.selector(a.memory) != c.selector(b.memory)
    for result in (a, b):
        assert c.Checkpoint.restore(json.loads(json.dumps(result.export()))) == result
        assert result.recover_initial_memory() == before.memory


def test_literal_signed_source_attachment():
    cert = ins.attachment_certificate()
    assert cert['all_E_corner_records_rank_at_most_one']
    assert [r['input_and_output_column_counts'] for r in cert['selected_instruments']] == [2, 2]


def test_complete_native_gram_on_all_corner_inputs():
    for bit in (0, 1):
        ops = ins.corner(bit)
        assert np.array_equal(sum(a.T@a for a in ops), 2*np.eye(3))
        assert np.array_equal(sum(a@a.T for a in ops), 2*np.eye(3))
    assert ins.balance_certificate()['complete_and_unital_exact_integer_check']


def test_full_Fock_controller_marginal_without_corner_restriction():
    # Independent dictionary implementation includes the nonzero full-Fock
    # outside outcome, unlike the invariant three-state corner.
    rho = sparse.diags([3**j.bit_count()/4**12 for j in range(4096)], format='csr')
    outputs = [sparse.csr_matrix(rho.shape, dtype=complex) for _ in range(4)]
    source = {bit: dictionary.full_extension(3, law=('C', 'E')[bit],
                                             seam=((0, 1), (0, 2))[bit]) for bit in (0, 1)}
    for first, second in product(range(2), repeat=2):
        for record, op in source[first]:
            after = c.advance((first, second), record)
            outputs[2*after[0]+after[1]] += op@rho@op.conj().T/4
    assert max(sparse.linalg.norm(out-rho/4) for out in outputs) < 1e-15


def test_all_second_order_zero_and_nonzero_records_are_addressable():
    for memory in product(range(2), repeat=2):
        gram = np.zeros((3, 3), dtype=object)
        for records in product(c.LABELS, repeat=2):
            checkpoint, matrix, denom = history.branch(memory, records)
            assert denom == 4
            assert checkpoint.recover_initial_memory() == memory
            gram += matrix.T@matrix
        assert np.array_equal(gram, 4*np.eye(3, dtype=int))
    _, outside, _ = history.branch((0, 0), [('outside',)])
    assert not np.any(outside)


def test_independent_path_products_include_field_coherences():
    expected, _ = ins.exact_channels()
    for memory in product(range(2), repeat=2):
        paths = [(c.Checkpoint(memory), np.eye(3, dtype=np.int64))]
        for _ in range(8):
            following = []
            for checkpoint, before in paths:
                for record, op in zip(c.LABELS, ins.corner(checkpoint.memory[0]), strict=True):
                    value = op@before
                    if np.any(value): following.append((checkpoint.commit(record), value))
            paths = following
        result = np.zeros((36, 9), dtype=np.int64)
        rank_one_counts = np.zeros((12, 3), dtype=np.int64)
        for checkpoint, op in paths:
            after = 2*checkpoint.memory[0]+checkpoint.memory[1]
            result[9*after:9*after+9] += np.kron(op, op)
            initial = memory
            saw_e = False
            for record in checkpoint.history:
                saw_e |= initial[0] == 1
                initial = c.advance(initial, record)
            if saw_e:
                assert np.count_nonzero(op) == 1
                rank_one_counts[3*after:3*after+3] += op*op
        before = 2*memory[0]+memory[1]
        assert np.array_equal(result, np.linalg.matrix_power(expected, 8)[:, 9*before:9*before+9])
        assert np.min(rank_one_counts) >= 1


def test_minorization_is_CP_not_just_population_mixing():
    channel, _ = ins.exact_channels()
    power = np.linalg.matrix_power(channel, 8)/256
    for before, after in product(range(4), repeat=2):
        # Choi matrix is an independent matrix-unit reshuffle, including all
        # off-diagonal input columns. Subtract the claimed depolarizing map.
        choi = np.zeros((9, 9))
        for i, j in product(range(3), repeat=2):
            out = power[9*after:9*after+9, 9*before+3*i+j].reshape(3, 3)
            choi[3*i:3*i+3, 3*j:3*j+3] = out
        choi -= np.eye(9)/256
        assert np.linalg.eigvalsh(choi).min() >= -1e-14
    assert ins.mixing_certificate()['depolarizing_weight'] == F(3, 64)


def test_arbitrary_reference_contraction_on_CQ_domain():
    rng = np.random.default_rng(1138)
    rho = []
    for _ in range(4):
        a = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
        rho.append(a@a.conj().T)
    norm = sum(np.trace(x).real for x in rho)
    rho = [x/norm for x in rho]
    reference = sum(x.reshape(3, 2, 3, 2).trace(axis1=0, axis2=2) for x in rho)
    stationary = np.kron(np.eye(3)/12, reference)
    distance = lambda blocks: sum(np.abs(np.linalg.eigvalsh(x-stationary)).sum() for x in blocks)
    initial = distance(rho)
    for _ in range(8):
        output = [np.zeros((6, 6), complex) for _ in range(4)]
        for first, second in product(range(2), repeat=2):
            for record, op in zip(c.LABELS, ins.corner(first), strict=True):
                after = c.advance((first, second), record)
                full = np.kron(op, np.eye(2))
                output[2*after[0]+after[1]] += full@rho[2*first+second]@full.T/2
        rho = output
    assert distance(rho) <= F(61, 64)*initial+1e-13


def test_archive_retains_initial_information_and_constant_marginal():
    cert = history.enumerate_positive(8)
    assert cert['positive_basis_histories'] == 3072
    assert cert['exact_output_population'] == [F(1, 12)]*12
    assert cert['initial_memory_and_field_basis_recovered_from_final_state_and_history']
    assert not cert['full_history_channel_forgets_initial_information']


def test_one_bit_reversibility_is_not_a_balance_proof():
    control = ins.negative_controls()
    assert control['uniform_input_trace_distance_after_one_tick'] == F(1, 6)
    assert control['one_bit_map_reversible_for_each_fixed_record']


def test_nonuniform_working_memory_is_not_stationary():
    channel, _ = ins.exact_channels()
    initial = np.zeros(36); initial[:9] = np.eye(3).ravel()/3
    assert np.linalg.norm(channel@initial/2-initial) > .1


def test_volume_uniform_local_error_and_growing_archive_cost():
    one = history.resource_certificate(4, 0, 4)
    many = history.resource_certificate(10000, 1000000, 4)
    assert one['local_reference_log_error'] == many['local_reference_log_error']
    assert many['record_storage_bits_upper'] == 5000006
    assert one['native_seed_or_memory_rounding_error'] == 0
    assert many['maximum_tick_record_alphabet_with_abort_bypass'] == 18


@pytest.mark.parametrize('ticks', [-1, True, 1.5, 1000001])
def test_invalid_resources(ticks):
    with pytest.raises(ValueError): history.resource_certificate(ticks=ticks)


def test_scopes_and_no_success_renormalization():
    row = history.continuing_theorem()
    assert row['no_success_normalization']
    assert not row['full_growing_history_is_stationary']
    assert not row['arbitrary_feedback_policies_are_covered']
    assert not ins.mixing_certificate()['full_Fock_or_interleaved_H_mixing_claimed']


def test_retained_receipt_and_every_source_hash():
    old = json.loads(packet.OUTPUT.read_text())
    assert old == packet.run()
    for path, value in old['source_sha256'].items():
        assert hashlib.sha256((packet.ROOT/path).read_bytes()).hexdigest() == value


def test_cli_rejects_modified_receipt(tmp_path):
    old = json.loads(packet.OUTPUT.read_text())
    old['finite_quantum_mixing']['depolarizing_weight'] = '1/2'
    corrupt = tmp_path/'corrupt.json'; corrupt.write_text(json.dumps(old))
    result = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.adaptive_reference',
                             '--check', str(corrupt)], cwd=packet.ROOT, capture_output=True, text=True)
    assert result.returncode != 0
