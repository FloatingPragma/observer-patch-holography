"""Complete recovery, actual noisy controls and hostile certificate rejection."""

import copy
import fnmatch
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest

from . import algebra, algebra_check, build, executor, instrument_check, recovery, recovery_check, verify


@pytest.fixture(scope='module')
def evidence():
    return build.candidate()


@pytest.fixture(scope='module')
def responses():
    tape, flags, count = recovery.program()
    columns, spans = recovery.responses(tape, count)
    return tape, flags, columns, spans


def test_final_receipt_and_all_parent_evidence():
    verify.verify()


def test_rebuilt_candidate_has_same_meaning(evidence):
    stored = verify.parent.parent.load(verify.HERE/'receipt.json')['evidence']
    verify.verify_evidence(evidence)
    verify.verify_evidence(stored)
    for field in ('decoder', 'recovery', 'scaling'):
        assert evidence[field] == stored[field]


def test_all_noisy_recovery_locations_checked_both_directions(evidence):
    result = recovery_check.replay(recovery.CONFIG)
    assert result == evidence['recovery']
    assert result['locations'] == 19664
    assert result['one_location_faults'] == 65328
    assert sum(result['residual_histogram'].values()) == 65328
    assert not result['failures']


@pytest.mark.parametrize('config,count', [
    (dict(recovery.CONFIG, verification=[]), 192),
    (dict(recovery.CONFIG, rounds=1), 1565),
])
def test_bare_or_unrepeated_recovery_really_fails(config, count):
    produced, independent = recovery.audit(config), recovery_check.replay(config)
    assert produced == independent
    assert len(independent['failures']) == count
    for index, labels, _ in independent['failures'][:4]:
        actual = executor.execute((index, labels), config=config)['residual']
        assert actual not in recovery.correctable()


def test_actual_adaptive_execution_matches_symbolic_census(responses):
    tape, flags, columns, spans = responses
    examples = list(recovery.single_faults(columns, spans))
    chosen = examples[::1999]
    # Explicitly exercise rejection of the first cat, a data idle, reads,
    # and final conditional corrections rather than only default success.
    first_reject = next(row for row in examples if any(row[2] >> b & 1 for fs in flags for b in fs[0]))
    chosen.append(first_reject)
    for op in ('i', 'mz', 'r', 'cx', 'cz', 'correct_x', 'correct_z'):
        chosen.append(next(row for row in reversed(examples) if tape[row[0]][0] == op))
    for index, labels, signature in chosen:
        actual = executor.execute((index, labels))
        assert actual['residual'] == recovery.interpret(signature, flags, 4)
        assert actual['locations'] == len(tape)
    branch = executor.execute((first_reject[0], first_reject[1]))
    assert any(first and not second for first, second in branch['candidate_rejections'])


def test_exhausted_local_supply_keeps_output_and_all_later_slots(responses):
    tape, flags, _, _ = responses
    first, second = flags[0]
    # One read fault rejects each cat; exhaustion needs both faults.
    indices = [next(i for i, (op, _, b) in enumerate(tape) if op == 'mz' and b == flag)
               for flag in (first[0], second[0])]
    result = executor.execute(faults=[(i, [1]) for i in indices])
    assert result['candidate_rejections'][0] == [True, True]
    assert result['local_failure'] and result['residual'] is None
    assert result['continued_frame'] == (0, 0)
    assert result['locations'] == len(tape)
    assert len(result['syndrome_history']) == 4


@pytest.mark.parametrize('x,z', [(0, 0)]+[(x << q, z << q) for q in range(7) for x, z in ((1, 0), (0, 1), (1, 1))])
def test_actual_recovery_corrects_incoming_error(x, z):
    result = executor.execute(incoming=(x, z))
    assert result['residual'] in recovery.stabilizers()
    assert result['syndrome_history'] == [recovery.syndrome(x, z)]*4


def test_initial_syndrome_translation_preserves_selection(responses):
    _, flags, columns, spans = responses
    for _, _, signature in list(recovery.single_faults(columns, spans))[::31]:
        if any(signature >> b & 1 for fs in flags for b in fs[0]):
            continue  # actual replacement is clean; incoming syndrome is constant
        history = [(signature >> (14+6*r)) & 63 for r in range(4)]
        selected = next(i for i in range(1, 4) if history[i] == history[i-1])
        for initial in range(64):
            changed = [x ^ initial for x in history]
            assert next(i for i in range(1, 4) if changed[i] == changed[i-1]) == selected


def test_ideal_recovery_on_every_syndrome_subspace():
    for s in range(64):
        result = executor.execute(incoming=recovery.correction(s))
        assert result['residual'] == (0, 0)
        assert result['syndrome_history'] == [s]*4


def test_full_quantum_syndrome_instrument(responses):
    instrument_check.check_program(responses[0])


@pytest.mark.parametrize('change', ['cat_h', 'coupling', 'basis', 'cat_read'])
def test_quantum_instrument_rejects_wrong_ideal_circuit(responses, change):
    tape = copy.deepcopy(responses[0])
    if change == 'cat_h':
        index = next(i for i, g in enumerate(tape) if g[:2] == ['h', [7]])
        tape.pop(index)
    elif change in ('coupling', 'basis'):
        index = next(i for i, (op, qs, _) in enumerate(tape)
                     if op == 'cx' and qs[0] == 7 and qs[1] < 7)
        if change == 'coupling':
            tape.pop(index)
        else:
            tape[index][0] = 'cz'
    else:
        index = next(i for i, g in enumerate(tape) if g[:2] == ['mz', [7]])
        tape.pop(index)
    with pytest.raises(ValueError):
        instrument_check.check_program(tape)


@pytest.mark.parametrize('kwargs', [
    {'fault': (-1, [1])}, {'fault': (999999, [1])}, {'fault': (True, [1])},
    {'fault': (0., [1])}, {'fault': (0, [])}, {'fault': (0, [1, 2])},
    {'fault': (0, [True])}, {'fault': (0, [-1])}, {'fault': (0, [4])},
    {'fault': (0, [0])}, {'fault': (0, [1.])}, {'fault': [0]},
    {'faults': [(0, [1]), (0, [2])]}, {'fault': (0, [1]), 'faults': [(0, [2])]},
    {'faults': {0: [1]}}, {'incoming': (128, 0)}, {'incoming': (0, -1)},
    {'incoming': (True, 0)}, {'incoming': (0,)}, {'config': dict(recovery.CONFIG, rounds=True)},
])
def test_adaptive_replay_rejects_invalid_fault_requests(kwargs):
    with pytest.raises(ValueError):
        executor.execute(**kwargs)


def test_decoder_corrects_nonpauli_channel_with_spectator():
    code, maps = algebra_check.check_decoder(algebra.decoder())
    gamma, bit = .23, 1 << 4
    identity = np.eye(128)
    z, x, xz = [algebra_check.pauli(a, b) for a, b in ((0, bit), (bit, 0), (bit, bit))]
    errors = [((1+np.sqrt(1-gamma))*identity+(1-np.sqrt(1-gamma))*z)/2,
              np.sqrt(gamma)*(x-xz)/2]
    assert np.allclose(sum(e.T@e for e in errors), identity)
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    result = np.zeros((4, 4))
    for e in errors:
        for k in maps:
            v = np.kron(k@e@code, np.eye(2))@bell
            result += np.outer(v, v)
    assert np.allclose(result, np.outer(bell, bell))


def test_full_encoded_basis_transversal_clifford_identities():
    code, _ = algebra_check.check_decoder(algebra.decoder())
    h = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
    tensor_h = np.array([[1.]])
    for _ in range(7):
        tensor_h = np.kron(tensor_h, h)
    assert np.allclose(tensor_h@code, code@h)
    physical_s_dagger = np.array([(-1j)**i.bit_count() for i in range(128)])
    assert np.allclose(physical_s_dagger[:, None]*code, code@np.diag([1., 1j]))
    pair_code = np.kron(code, code)
    # First block controls the second at every one of the seven positions.
    indices = np.arange(128*128)
    physical_cx = ((indices >> 7) << 7) | ((indices & 127) ^ (indices >> 7))
    logical_cx = [0, 1, 3, 2]
    assert np.allclose(pair_code[physical_cx], pair_code[:, logical_cx])


def test_nonclifford_reference_keeps_both_outcomes_and_controlled_phase(evidence):
    raw = np.array(evidence['interfaces']['t'])
    t = raw[..., 0]+1j*raw[..., 1]
    x, z, s = np.array([[0., 1.], [1., 0.]]), np.diag([1., -1.]), np.diag([1., 1j])
    q = t@x@t.conj().T
    assert np.allclose(q, q.conj().T) and np.allclose(q@q, np.eye(2))
    assert np.allclose(q, np.exp(-1j*np.pi/4)*s@x)
    assert not np.allclose(s@x, (s@x).conj().T)
    plus = t@np.ones(2)/np.sqrt(2)
    minus = z@plus
    assert np.allclose(q@plus, plus) and np.allclose(q@minus, -minus)
    assert np.allclose(z@minus, plus)
    cx = np.eye(4)[[0, 1, 3, 2]]
    output = cx@np.kron(np.eye(2), plus[:, None])
    first, second = output[[0, 2]], output[[1, 3]]
    branches = [first, s@second]
    assert np.allclose(sum(k.conj().T@k for k in branches), np.eye(2))
    for k in branches:
        relative = t.conj().T@k
        assert np.allclose(relative, np.eye(2)*np.trace(relative)/2)
    # Omitting feed-forward leaves the same read probabilities but corrupts
    # the unconditional quantum channel.
    wrong_fidelity = sum(abs(np.trace(t.conj().T@k)/2)**2 for k in (first, second))
    assert np.isclose(wrong_fidelity, .75)


def test_same_syndrome_probabilities_do_not_certify_logical_output(evidence):
    bad = copy.deepcopy(evidence['decoder'])
    # Logical Z after every branch leaves all branch probabilities and
    # completeness unchanged, but damages arbitrary quantum inputs.
    for row in bad:
        for entry in row['kraus_rows'][1]:
            entry[1] *= -1
    with pytest.raises(ValueError, match='decoder branch'):
        algebra_check.check_decoder(bad)


def test_native_gate_and_memory_boundaries(evidence):
    assert evidence['memory']['logical_error_by_weight'] == [0, 0, 21, 7, 28, 0, 7, 1]
    assert evidence['accounting']['raw_input_entangled_fidelity'] == .875
    assert evidence['accounting']['after_encoding_corrected_fidelity'] == 1.
    assert evidence['scaling']['levels'][2]['accounting_q_degree'] == 0
    assert evidence['scaling']['levels'][3]['accounting_q_degree'] < 0


def test_full_native_carrier_extension_charges_complement_without_postselection(evidence):
    raw = np.array(evidence['accounting']['logical_effect'])
    effect = raw[..., 0]+1j*raw[..., 1]
    # Three output dimensions: logical qubit plus an orthogonal failure.
    keep = np.zeros((3, 6))
    keep[:2, :2] = np.eye(2)
    branches = [keep]
    for q in range(2, 6):
        failure = np.zeros((3, 6))
        failure[2, q] = 1
        branches.append(failure)
    assert np.allclose(sum(k.T@k for k in branches), np.eye(6))
    output_effect = np.eye(3, dtype=complex)
    output_effect[:2, :2] = effect
    physical = sum(k.T@output_effect@k for k in branches)
    expected = np.eye(6, dtype=complex)
    expected[:2, :2] = effect
    assert np.allclose(physical, expected)
    state = np.array([1., 0., 1., 0., 0., 0.])/np.sqrt(2)
    actual = np.vdot(state, physical@state).real
    assert np.isclose(actual, .5*effect[0, 0].real+.5)
    assert not np.isclose(actual, effect[0, 0].real)


MUTATIONS = [
    (('decoder',), []),
    (('decoder', 3, 'syndrome'), 0),
    (('decoder', 0, 'kraus_rows'), []),
    (('decoder', 0, 'kraus_rows', 0, 0, 1), 0),
    (('decoder', 0, 'denominator_squared'), 16),
    (('interfaces', 't', 1, 1), [1., 0.]),
    (('interfaces', 'cx', 1, 3), [0., 0.]),
    (('recovery', 'config', 'rounds'), 1),
    (('recovery', 'config', 'candidates'), True),
    (('recovery', 'config', 'verification'), []),
    (('recovery', 'one_location_faults'), 1),
    (('recovery', 'idle_locations'), 0),
    (('recovery', 'residual_histogram'), {}),
    (('recovery', 'failures'), [[1, [1], 'abort']]),
    (('recovery', 'response_sha256'), '0'*64),
    (('memory', 'logical_error_by_weight'), [0]*8),
    (('memory', 'raw_entangled_fidelity_by_weight'), [1]*8),
    (('memory', 'single_phase_branches'), []),
    (('continuous_noise',), []),
    (('continuous_noise', 0, 'rate'), 0.),
    (('continuous_noise', 0, 'zero_jump_probability'), 1.),
    (('continuous_noise', 0, 'channel_choi', 0, 0, 0), .9),
    (('continuous_noise', 0, 'fault_choi', 0, 0, 0), -1.),
    (('accounting', 'cases'), []),
    (('accounting', 'cases', 1, 'decoded_expectation'), 1.),
    (('accounting', 'cases', 1, 'raw_code_penalty'), 0.),
    (('accounting', 'cases', 0, 'raw_code_penalty'), False),
    (('accounting', 'abort_value'), 0.),
    (('accounting', 'raw_input_entangled_fidelity'), 1.),
    (('scaling', 'selected_level'), 3),
    (('scaling', 'levels', 2, 'accounting_q_degree'), -1),
    (('scaling', 'local_fault_q_degree'), -8),
    (('scaling', 'full_library_constant'), 19664),
    (('scaling', 'auxiliary_records'), 'discard failures'),
    (('scaling', 'local_failure'), 'abort experiment on every lower-level flag'),
]


@pytest.mark.parametrize('path,value', MUTATIONS, ids=[str(p) for p, _ in MUTATIONS])
def test_hostile_physics_rejected_with_hashes_bypassed(evidence, path, value):
    bad = copy.deepcopy(evidence)
    cursor = bad
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError):
        verify.verify_evidence(bad)


@pytest.mark.parametrize('config', [dict(recovery.CONFIG, rounds=0), dict(recovery.CONFIG, rounds=True),
    dict(recovery.CONFIG, candidates=1), dict(recovery.CONFIG, verification=[[0, 0]]),
    dict(recovery.CONFIG, verification=[[0, 4]]), dict(recovery.CONFIG, verification=[[0, False]])])
def test_invalid_recovery_configuration(config):
    with pytest.raises(ValueError):
        recovery.program(config)


@pytest.mark.parametrize('text', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}'])
def test_strict_json(tmp_path, text):
    path = tmp_path/'invalid.json'
    path.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError):
        verify.parent.parent.load(path)


def test_optimized_checker_rejects_forged_physics_and_custody(tmp_path):
    packet = verify.parent.parent.load(verify.HERE/'receipt.json')
    for part in ('sources', 'evidence'):
        bad = copy.deepcopy(packet)
        if part == 'sources':
            bad[part][next(iter(bad[part]))] = '0'*64
        else:
            bad[part]['decoder'] = []
        path = tmp_path/(part+'.json')
        path.write_text(json.dumps(bad), encoding='utf-8')
        result = subprocess.run([sys.executable, '-O', '-m', 'm1_fixed_noise.verify', str(path)],
                                env=dict(os.environ, PYTHONPATH=str(verify.ROOT/'code')),
                                cwd=verify.ROOT, capture_output=True, text=True)
        assert result.returncode != 0 and 'ValueError' in result.stderr


def test_verifier_with_producers_disabled():
    script = 'import sys\n'
    for name in ('algebra', 'build', 'executor', 'noise', 'recovery', 'scaling'):
        script += f"sys.modules['m1_fixed_noise.{name}'] = None\n"
    script += 'from m1_fixed_noise import verify\nverify.verify()\n'
    result = subprocess.run([sys.executable, '-O', '-c', script],
                            env=dict(os.environ, PYTHONPATH=str(verify.ROOT/'code')),
                            cwd=verify.ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


