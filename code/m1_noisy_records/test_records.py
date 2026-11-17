"""Channel-level, fault-level and hostile verifier controls."""

import copy
import fnmatch
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import archive, archive_check, circuits, independent as check, verify
from . import computation_code, computation_code_check
from . import instruments, instrument_check


@pytest.fixture(scope='module')
def evidence():
    return circuits.candidate()


def test_receipt_and_all_parent_custody():
    verify.verify()


def test_rebuild_has_same_discrete_evidence(evidence):
    assert evidence == verify.load(verify.HERE/'receipt.json')['evidence']
    check.verify_evidence(evidence)


def test_explicit_computation_code_meets_published_spread_hypothesis(evidence):
    computation_code_check.verify(evidence['computation_code'])
    assert evidence['computation_code']['corrects'] == 31
    assert evidence['computation_code']['correction_then_gate_spread_bound'] == 16
    assert computation_code.rows(5, 2) == computation_code_check.monomials(5, 2)
    assert computation_code_check.small_complete_distances() == [(3, 1, 16, 4), (4, 1, 32, 8), (5, 2, 65536, 8)]


@pytest.mark.parametrize('key', ['corrects', 'quantum_distance', 'quantum_dimension',
                                'procedure_spread_bound', 'correction_then_gate_spread_bound',
                                'tolerated_faults_per_rectangle', 'generator_sha256'])
def test_code_margin_cannot_be_asserted_by_receipt(key, evidence):
    row = dict(evidence['computation_code'])
    row[key] = '0'*64 if key == 'generator_sha256' else row[key]+1
    with pytest.raises(ValueError):
        computation_code_check.verify(row)


@pytest.mark.parametrize('angles', circuits.CASES)
def test_two_independent_full_instrument_computations(angles):
    tape = circuits.tape(*angles)
    v = circuits.execute(tape)
    np.testing.assert_allclose(v, check.isometry(tape), atol=1e-14)
    check.close(check.choi_from_isometry(v), check.reference(angles), 'all Choi entries')
    # A four-dimensional reference witnesses all input coherences, not just
    # computational input distributions. The Choi state has trace four.
    choi = check.choi_from_isometry(v)
    assert abs(np.trace(choi)-4) < 1e-12
    assert np.linalg.eigvalsh(choi).min() > -1e-12


@pytest.mark.parametrize('change', ['environment', 'feedback', 'abort_phase', 'accept_rotation'])
def test_physically_wrong_tapes_fail_without_hash_checks(change):
    angles = circuits.CASES[0]
    row = dict(angles=list(angles), tape=circuits.tape(*angles))
    ops = row['tape']
    if change == 'environment':
        ops.remove(['cx', 2, 5])
    elif change == 'feedback':
        ops.remove(['cz', 3, 0])
    elif change == 'abort_phase':
        ops.remove(['tdg', 4])
    else:
        ops[-1][2] *= -1
    with pytest.raises(ValueError, match='complete history instrument'):
        check.verify_case(row, angles)


def test_toffoli_phase_error_hidden_from_truth_tables_is_rejected():
    tape = circuits.toffoli(0, 1, 2)+[['t', 0]]
    actual = circuits.execute(tape, qubits=3, inputs=3)
    ideal = circuits.execute(circuits.toffoli(0, 1, 2), qubits=3, inputs=3)
    np.testing.assert_allclose(abs(actual)**2, abs(ideal)**2, atol=1e-14)
    with pytest.raises(ValueError, match='coherent phases'):
        check.verify_toffoli(tape)


def test_grouped_native_complement_retains_exact_full_channel(evidence):
    instrument_check.verify(evidence['native_instruments'])


@pytest.mark.parametrize('d', [2, 4])
@pytest.mark.parametrize('mutation', ['omit_failure', 'coherent_failures',
                                     'public_coherence', 'wrong_reset'])
def test_wrong_native_complements_fail_without_custody(d, mutation):
    row = instruments.candidate()[0 if d == 2 else 1]
    v = instrument_check.decode_matrix(row).reshape(7-d, 2, d+1, 6)
    if mutation == 'omit_failure':
        v[1:] = 0
    elif mutation == 'coherent_failures':
        v[1] += v[2]
        v[2] = 0
    elif mutation == 'public_coherence':
        # Still an isometry, but the public success/failure register is no
        # longer classical on inputs coherent across code and complement.
        v[0] += v[1]
        v[1] = 0
    else:
        # A TP instrument with the correct failure probability but wrong
        # conditional output must also be rejected.
        v[1, 1, 0, d] = v[1, 1, d, d]
        v[1, 1, d, d] = 0
    flat = v.reshape(-1, 6)
    if mutation in ('public_coherence', 'wrong_reset'):
        np.testing.assert_allclose(flat.conj().T @ flat, np.eye(6), atol=1e-14)
    with pytest.raises(ValueError, match='native'):
        instrument_check.check_dilation(flat, d)


def test_grouped_reset_and_public_feedback_on_all_input_coherences():
    # One reset outcome contains two private Kraus branches. The other outcome
    # leaves the input unchanged. A later X is controlled only by public y.
    reset = [np.array([[1, 0], [0, 0]]), np.array([[0, 1], [0, 0]])]
    groups = [[np.sqrt(.3)*k for k in reset], [np.sqrt(.7)*np.eye(2)]]
    v = instruments.dilate(groups).reshape(3, 4, 2)
    feedback = np.block([[np.array([[0, 1], [1, 0]]), np.zeros((2, 2))],
                         [np.zeros((2, 2)), np.eye(2)]])
    for i in range(2):
        for j in range(2):
            actual = sum(np.outer(k[:, i], k[:, j].conj()) for k in v)
            actual = feedback @ actual @ feedback.conj().T
            expected = np.zeros((4, 4))
            if i == j:
                expected[1, 1] = .3
            expected[2+i, 2+j] = .7
            np.testing.assert_allclose(actual, expected, atol=1e-14)
    # Reset branches may be rotated within their private environment without
    # changing the observable instrument. The verifier checks the channel,
    # not an arbitrary choice of Kraus representation.
    native = instruments.candidate()[0]
    w = instrument_check.decode_matrix(native).reshape(5, 6, 6)
    a, b = w[1].copy(), w[2].copy()
    w[1], w[2] = (a+b)/np.sqrt(2), (a-b)/np.sqrt(2)
    instrument_check.check_dilation(w.reshape(-1, 6), 2)


@pytest.mark.parametrize('groups', [None, [], [[]], ([],), [[np.eye(2)*.9]],
                                    [[np.eye(2), np.eye(3)]], [[np.ones(2)]],
                                    [[np.eye(2, dtype=bool)]], [[np.array([[np.nan]])]],
                                    [[np.array([[np.inf]])]], [[np.array([['1']])]],
                                    [[np.eye(17)]], [[np.eye(1)/np.sqrt(33)]*33]])
def test_incomplete_or_malformed_instruments_rejected(groups):
    with pytest.raises(ValueError):
        instruments.dilate(groups)


@pytest.mark.parametrize('field,value', [('code_dimension', True), ('code_dimension', 3),
                                       ('environment', 4), ('real', []), ('extra', 1)])
def test_native_instrument_schema_is_strict(field, value):
    rows = instruments.candidate()
    rows[0][field] = value
    with pytest.raises(ValueError):
        instrument_check.verify(rows)


@pytest.mark.parametrize('value', [True, '0', None, float('nan'), float('inf')])
def test_native_instrument_matrix_rejects_non_numeric_entries(value):
    rows = instruments.candidate()
    rows[0]['real'][0][0] = value
    with pytest.raises(ValueError):
        instrument_check.verify(rows)


@pytest.mark.parametrize('indices', [[], [0], [1, 0], [0, 0], [0, 1, 1]])
def test_native_instrument_catalog_is_complete(indices):
    rows = instruments.candidate()
    with pytest.raises(ValueError):
        instrument_check.verify([rows[i] for i in indices])


@pytest.mark.parametrize('op', [[], ['unknown', 0], ['cx', 1, 1], ['cx', -1, 2],
                               ['h', True], ['h', 0, 1], ['ry', 0, float('nan')],
                               ['ry', 0, float('inf')], ['ry', 0, True]])
def test_malformed_quantum_tapes(op):
    with pytest.raises(ValueError):
        check.isometry([op])
    with pytest.raises(ValueError):
        circuits.execute([op])


def test_archive_independent_full_fault_census(evidence):
    archive_check.verify(evidence['archive'])
    assert evidence['archive']['census'] == dict(
        locations=515, single_fault_cases=8100,
        output_error_histogram=[6365, 1735, 0, 0, 0, 0],
        live_input_error_histogram=[1130, 6090, 880, 0, 0, 0])


@pytest.mark.parametrize('n', [5, 9, 13])
def test_general_refresh_with_multiple_old_errors_and_multiple_faults(n):
    p = archive.schedule(n)
    ops = archive_check.validate(p)
    r = (n-1)//4
    rng = np.random.default_rng(812+n)
    for _ in range(24):
        bit = int(rng.integers(2))
        inputs = [bit]*n
        for i in rng.choice(n, r, replace=False):
            inputs[int(i)] ^= 1
        faults = [(int(i), int(rng.integers(1, 1 << (len(ops[int(i)])-1))))
                  for i in rng.choice(len(ops), r, replace=False)]
        result, live_wrong = archive_check.run(p, inputs, faults, include_live=True)
        assert sum(x != bit for x in result) <= r
        assert live_wrong <= 2*r


def test_two_fault_witness_is_not_vacuously_accepted():
    p = archive.schedule(5)
    ops = archive.locations(p)
    # The initial input already has one error. Two different initial idles
    # corrupt two more source bits before any voter copies them.
    faults = [(ops.index(['idle', bit]), 1) for bit in (1, 2)]
    result, live_wrong = archive_check.run(p, [1, 0, 0, 0, 0], faults, include_live=True)
    assert result == [1]*5
    assert live_wrong == 3


@pytest.mark.parametrize('option', [0, 1, None, 'yes'])
def test_live_trace_option_is_not_silently_coerced(option):
    p = archive.schedule(5)
    with pytest.raises(ValueError):
        archive.run(p, 0, include_live=option)
    with pytest.raises(ValueError):
        archive_check.run(p, [0]*5, include_live=option)


def test_false_live_census_fails_independent_replay(evidence):
    row = copy.deepcopy(evidence['archive'])
    row['census']['live_input_error_histogram'][0] += 1
    row['census']['live_input_error_histogram'][1] -= 1
    with pytest.raises(ValueError, match='exact fault census'):
        archive_check.verify(row)


@pytest.mark.parametrize('mutation', ['idle', 'shared', 'missing_voter', 'duplicate_wire', 'output', 'sort'])
def test_archive_structural_corruptions_fail_before_census(mutation):
    p = archive.schedule(5)
    if mutation == 'idle':
        p['layers'][0].remove(['idle', 0])
    elif mutation == 'shared':
        p['layers'][1][0][2] = 1
    elif mutation == 'missing_voter':
        p['layers'][-1][0] = ['idle', p['outputs'][0]]
    elif mutation == 'duplicate_wire':
        p['layers'][1].append(['idle', 0])
    elif mutation == 'output':
        p['outputs'][0] = 0
    else:
        op = next(op for layer in p['layers'] for op in layer if op[0] == 'sort')
        op[0] = 'copy'
    with pytest.raises(ValueError):
        archive_check.validate(p)


@pytest.mark.parametrize('faults', [[(-1, 1)], [(10000, 1)], [(0, 0)], [(0, 2)],
                                   [(0, True)], [(0, 1), (0, 1)], [('0', 1)], [(0,)]])
def test_no_silent_discard_of_trash_faults(faults):
    with pytest.raises(ValueError):
        archive_check.run(archive.schedule(5), [0]*5, faults)
    with pytest.raises(ValueError):
        archive.run(archive.schedule(5), 0, faults[0] if len(faults) == 1 else faults)


@pytest.mark.parametrize('word', [-1, 32, True, 1.0])
def test_producer_rejects_invalid_initial_words(word):
    with pytest.raises(ValueError):
        archive.run(archive.schedule(5), word)


def test_producer_rejects_numerically_equal_malformed_schedule():
    p = archive.schedule(5)
    p['width'] = float(p['width'])
    with pytest.raises(ValueError):
        archive.run(p, 0)
    p = archive.schedule(5)
    p['layers'][1][0][1] = False
    with pytest.raises(ValueError):
        archive.run(p, 0)


def test_all_resource_entries_are_load_bearing(evidence):
    for field in evidence['budget']:
        row = dict(evidence['budget'])
        row[field] += 1
        with pytest.raises(ValueError):
            check.verify_budget(row)
        row[field] = True
        with pytest.raises(ValueError):
            check.verify_budget(row)


def test_superpolynomial_archive_bound_has_proved_asymptotics():
    # Constants are illustrative, not the unknown physical library constants.
    q = 2**200
    value = archive_check.log_failure_majorant(q, rate_constant=17.)
    assert value < -20*math.log(q)
    q2 = 2**300
    assert archive_check.log_failure_majorant(q2, 17.)/math.log(q2) < value/math.log(q)


def test_terminal_noise_law_and_last_bit_contraction(evidence):
    check.verify_terminal(evidence['terminal'])
    for rate, horizon in ((.3, .5), (2., .17), (1.1, 2.)):
        generator = rate*np.array([[-1, 1], [1, -1]])
        actual = expm(horizon*generator)
        p = (1-math.exp(-2*rate*horizon))/2
        np.testing.assert_allclose(actual, [[1-p, p], [p, 1-p]], atol=1e-14)
        assert abs(np.linalg.norm(actual[:, 0]-actual[:, 1], 1)/2-(1-2*p)) < 1e-14
        failure = sum(math.comb(5, j)*p**j*(1-p)**(5-j) for j in range(3, 6))
        assert abs(failure-(10*p**3-15*p**4+6*p**5)) < 1e-14
        assert 0 < (1-failure)**100 < (1-failure)**10 < 1


def test_export_preserves_all_classical_quantum_correlations(evidence):
    check.verify_export(evidence['export'])
    bad = copy.deepcopy(evidence['export'])
    bad['tape'].append(['cx', 2, 0])  # leaking the public copy back into the application
    with pytest.raises(ValueError, match='public export'):
        check.verify_export(bad)
    bad = copy.deepcopy(evidence['export'])
    bad['tape'].pop()
    with pytest.raises(ValueError, match='public export'):
        check.verify_export(bad)


def test_passive_diagnostic_noise_is_harmless_only_without_feedback():
    bell = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    rho = np.outer(bell, bell.conj())
    flip_aux = np.kron(np.array([[0, 1], [1, 0]]), np.eye(2))
    noisy = .7*rho+.3*flip_aux@rho@flip_aux.conj().T
    partial = lambda state: np.einsum('aiaj->ij', state.reshape(2, 2, 2, 2))
    np.testing.assert_allclose(partial(noisy), partial(rho), atol=1e-14)
    # A passive flag becoming a control would invalidate that exclusion.
    clean = np.diag([1., 0., 0., 0.])
    corrupted = .7*clean+.3*flip_aux@clean@flip_aux.conj().T
    feedback = check.elementary(['cx', 1, 0], 2)
    assert np.linalg.norm(partial(feedback@corrupted@feedback.conj().T)-partial(clean)) > .1


def test_archive_majority_is_not_a_perfect_bare_gate():
    # A correct public codeword still gives a corrupt naked result if the
    # reader's last bit flips. This fault is not distributed over code lanes.
    word = [1]*5
    perfect_majority = int(sum(word) > len(word)//2)
    noisy_last_bit = perfect_majority ^ 1
    assert noisy_last_bit != 1


def test_classical_comparator_has_a_finite_native_boolean_realization():
    tape = circuits.toffoli(0, 1, 2)+[['cx', 0, 3], ['cx', 1, 3], ['cx', 2, 3]]
    actual = circuits.execute(tape, qubits=4, inputs=2)
    for j in range(4):
        a, b = j & 1, j >> 1
        output = j+4*(a & b)+8*(a | b)
        expected = np.zeros(16)
        expected[output] = 1
        np.testing.assert_allclose(actual[:, j], expected, atol=1e-14)
    # Keeping the old a,b as local garbage and exposing c,d is the entire
    # irreversible classical comparator. No state-dependent postselection.


def test_strict_json_rejects_duplicate_and_nonfinite(tmp_path):
    for raw in ('{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'):
        path = tmp_path/'bad.json'
        path.write_text(raw, encoding='utf-8')
        with pytest.raises(ValueError):
            verify.load(path)


def test_optimized_replay_with_producers_disabled(tmp_path):
    script = '''
import sys
class Block:
    def find_spec(self, fullname, path=None, target=None):
        if fullname in ('m1_noisy_records.circuits', 'm1_noisy_records.archive', 'm1_noisy_records.build', 'm1_noisy_records.computation_code', 'm1_noisy_records.instruments'):
            raise RuntimeError('producer called')
sys.meta_path.insert(0, Block())
from m1_noisy_records import verify
packet = verify.load(verify.HERE/'receipt.json')
verify.verify_evidence(packet['evidence'])
packet['evidence']['budget']['fault_power'] = 1
try:
    from m1_noisy_records.independent import verify_budget
    verify_budget(packet['evidence']['budget'])
except ValueError:
    pass
else:
    raise RuntimeError('optimized verifier accepted trash')
'''
    env = dict(os.environ, PYTHONPATH=str(verify.ROOT/'code'))
    result = subprocess.run([sys.executable, '-O', '-c', script], env=env,
                            cwd=verify.ROOT, text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout+result.stderr


def test_source_and_claim_tampering_is_rejected(tmp_path):
    original = verify.load(verify.HERE/'receipt.json')
    for section in ('sources', 'claims'):
        row = copy.deepcopy(original)
        row[section][next(iter(row[section]))] = '0'*64
        path = tmp_path/(section+'.json')
        path.write_text(json.dumps(row), encoding='utf-8')
        with pytest.raises(ValueError, match='custody'):
            verify.verify(path)
