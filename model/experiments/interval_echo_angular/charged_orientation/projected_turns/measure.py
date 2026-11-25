"""A complete native half-turn test with a saturating postselection code."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.symmetry import rotations
from ..measure import compile_turn, rotated

HERE = Path(__file__).resolve().parent
SHOTS, SEED = 2**18, 707153
I2 = np.eye(2)
X = np.array([[0., 1], [1, 0]])
Y = np.array([[0., -1j], [1j, 0]])
Z = np.diag([1., -1])
INPUTS = [(label, sign, np.linalg.eigh(axis)[1][:, (sign + 1) // 2])
          for label, axis in (('X', X), ('Y', Y), ('Z', Z)) for sign in (-1, 1)]


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def close(a, b, message):
    check(np.max(abs(np.asarray(a) - np.asarray(b))) < 8e-13, message)


def compare(a, b, path='projected_turns'):
    if isinstance(b, dict):
        check(isinstance(a, dict) and a.keys() == b.keys(), path)
        for key in b:
            compare(a[key], b[key], path + '.' + key)
    elif isinstance(b, list):
        check(isinstance(a, list) and len(a) == len(b), path)
        for i, (left, right) in enumerate(zip(a, b, strict=True)):
            compare(left, right, f'{path}[{i}]')
    elif isinstance(b, float):
        check(np.isfinite(a) and np.isclose(a, b, rtol=1e-10, atol=2e-12), path)
    else:
        check(type(a) is type(b) and a == b, path)


def setup():
    identity = tuple(range(12))
    halves = [tuple(r) for r in sorted(rotations())
              if tuple(r) != identity and all(r[r[i]] == i for i in identity)]
    r, s = next((a, b) for a, b in combinations(halves, 2)
                if all(a[b[i]] == b[a[i]] for i in identity))
    R, S = np.eye(12)[:, r], np.eye(12)[:, s]
    close(R @ S, S @ R, 'full physical half-turns do not commute')
    root = 0
    orbit = [root, r[root], s[root], r[s[root]]]
    check(len(set(orbit)) == 4, 'source orbit is not regular Klein four')
    signs = list(product((1, -1), repeat=2))
    chars = np.zeros((12, 4))
    rows = []
    for i, (a, b) in enumerate(signs):
        chars[orbit, i] = np.array([1, a, b, a * b]) / 2
        angle = np.arctan2(a, b) % (2 * np.pi)
        rows.append(np.array([np.cos(angle / 2), np.sin(angle / 2)]) / np.sqrt(2))
    code = chars @ np.asarray(rows)
    close(code.conj().T @ code, I2, 'source isometry fails')
    A, B = code.T @ R @ code, code.T @ S @ code
    close(A, X / np.sqrt(2), 'first compressed turn differs')
    close(B, Z / np.sqrt(2), 'second compressed turn differs')
    P = code @ code.T
    close(P @ P, P, 'code filter is not projective')
    pulse_words = [compile_turn(perm) for perm in (r, s)]
    for matrix, word in zip((R, S), pulse_words, strict=True):
        close(rotated(np.eye(12), word), matrix, 'native pulse word does not implement spatial turn')
    # Native pi contacts are positive load exchanges, including source/host words.
    swap = np.eye(9)[[law.SRC.index((b, a)) for a, b in law.SRC]]
    close(expm(-1j * np.pi * law.seam_term(8)), swap, 'native contact has another phase')
    return R, S, code, P, orbit, pulse_words


def branch(operations, P, history, state, fair_reference=False):
    v = state.copy()
    for operation, hit in zip(operations, history, strict=True):
        v = operation @ v
        v = v / np.sqrt(2) if fair_reference else (P if hit else np.eye(len(P)) - P) @ v
    return v


def complete(left, right, P, state, filtered, rng, fair_reference=False):
    histories = list(product((True, False), repeat=len(left))) if filtered else [()]
    probabilities, labels = [], []
    for history in histories:
        if filtered:
            a = branch(left, P, history, state)
            b = branch(right, P, history, state, fair_reference)
        else:
            a, b = state.copy(), state.copy()
            for operation in left:
                a = operation @ a
            for operation in right:
                b = operation @ b
        probabilities += [float(np.linalg.norm((a + b) / 2)**2),
                          float(np.linalg.norm((a - b) / 2)**2)]
        labels += [dict(code_acceptance=list(history), interference=name) for name in ('plus', 'minus')]
    check(min(probabilities) >= 0 and abs(sum(probabilities) - 1) < 2e-13,
          'complete instrument lost a rejection history')
    counts, complement = sample_acquisition(probabilities, rng)
    return dict(labels=labels, probabilities=probabilities, counts=counts.tolist(),
                sampling_l1_complement_bound=complement)


def sample_acquisition(probabilities, rng):
    raw = np.asarray(probabilities)
    order = np.argsort(raw, kind='stable')
    sampling = raw[order].copy()
    sampling[-1] = 1 - sum(sampling[:-1])
    check(sampling[-1] >= 0, 'negative explicit complement')
    counts = np.empty(len(raw), np.int64)
    counts[order] = rng.multinomial(SHOTS, sampling)
    return counts, float(abs(sampling[-1] - raw[order[-1]]))


def calibrated_turns(strength, exchange_rate, words, R, S):
    kinetic = exchange_rate * np.array([[0.5, -0.5], [-0.5, 0.5]])
    pulse = expm(-1j * np.pi * (kinetic + strength * np.eye(2)))
    implemented = []
    for word in words:
        op = np.eye(12, dtype=complex)
        for a, b in word:
            step = np.eye(12, dtype=complex)
            step[np.ix_([a, b], [a, b])] = pulse
            op = step @ op
        implemented.append(op)
    error = [float(np.linalg.norm(op - target, 2))
             for op, target in zip(implemented, (R, S), strict=True)]
    return implemented, dict(contrast=strength, exchange_rate=exchange_rate,
                             spatial_pi_calibration_errors=error,
                             full_turn_calibration_errors=[float(np.linalg.norm(op @ op - np.eye(12), 2))
                                                           for op in implemented],
                             raw_order_commutator_norm=float(np.linalg.norm(implemented[0] @ implemented[1]
                                                                           - implemented[1] @ implemented[0], 2)))


def run():
    R, S, V, P, orbit, words = setup()
    A, B = V.T @ R @ V, V.T @ S @ V
    leaks = [float(np.linalg.norm((np.eye(12) - P) @ op @ V, 2)) for op in (R, S)]
    commutator = float(np.linalg.norm(A @ B - B @ A, 2))
    close(commutator, 2 * leaks[0] * leaks[1], 'leakage bound is not saturated')
    close(A @ A, I2 / 2, 'compressed literal full turn has another sign')
    rng = np.random.default_rng(SEED)
    records = []
    for charge in (-1, 1):
        for axis, orientation, logical in INPUTS:
            state = V @ logical
            panels = dict(
                unfiltered_order=complete([R, S], [S, R], P, state, False, rng),
                filtered_order=complete([R, S], [S, R], P, state, True, rng),
                unfiltered_full_turn=complete([R, R], [np.eye(12)] * 2, P, state, False, rng),
                filtered_full_turn=complete([R, R], [np.eye(12)] * 2, P, state, True, rng, True),
                unfiltered_double_full_turn=complete([R] * 4, [np.eye(12)] * 4, P, state, False, rng),
                filtered_double_full_turn=complete([R] * 4, [np.eye(12)] * 4, P, state, True, rng, True))
            close(panels['filtered_order']['probabilities'][:2], [0, 0.25],
                  'selected opposite-order phase is not negative')
            close(panels['filtered_full_turn']['probabilities'][:2], [0.25, 0],
                  'selected literal full-turn phase is not positive')
            close(panels['filtered_double_full_turn']['probabilities'][:2], [1 / 16, 0],
                  'selected double-full-turn phase is not positive')
            records.append(dict(charge=charge, axis=axis, orientation=orientation, panels=panels))
    spinor = dict(half_turn_x_real=(-1j * X).real.tolist(),
                  half_turn_x_imag=(-1j * X).imag.tolist())
    sx, sz = -1j * X, -1j * Z
    close(sx @ sz, -sz @ sx, 'supplied spinor order control fails')
    close(sx @ sx, -I2, 'supplied spinor full-turn control fails')
    spinor['interference'] = dict(
        order=complete([sx, sz], [sz, sx], I2, INPUTS[0][2], False, rng),
        full_turn=complete([sx, sx], [I2, I2], I2, INPUTS[0][2], False, rng),
        double_full_turn=complete([sx] * 4, [I2] * 4, I2, INPUTS[0][2], False, rng))
    controls = []
    for strength, rate in ((4, 1), (8, 0.5)):
        (r, s), calibration = calibrated_turns(strength, rate, words, R, S)
        if rate == 1:
            check(max(calibration['spatial_pi_calibration_errors']) < 1e-12,
                  'even-contrast control does not preserve exact spatial turns')
        else:
            check(max(calibration['spatial_pi_calibration_errors']) > 0.1,
                  'half-exchange calibration control failed to detect a changed pulse')
        bank = []
        for axis, orientation, logical in INPUTS:
            state = V @ logical
            bank.append(dict(axis=axis, orientation=orientation,
                             unfiltered_order=complete([r, s], [s, r], P, state, False, rng),
                             filtered_order=complete([r, s], [s, r], P, state, True, rng),
                             unfiltered_full_turn=complete([r, r], [np.eye(12)] * 2, P, state, False, rng),
                             filtered_full_turn=complete([r, r], [np.eye(12)] * 2, P, state, True, rng, True)))
        controls.append(dict(calibration=calibration, records=bank))
    return dict(shots=SHOTS, seed=SEED, source_orbit=orbit, native_pi_words=[[list(edge) for edge in word] for word in words],
                code_real=V.real.tolist(), code_imag=V.imag.tolist(),
                ideal_source_acceptance=1, ideal_source_refusal=0,
                per_turn_code_success=[0.5, 0.5], leakage_norms=leaks,
                compressed_commutator_norm=commutator,
                exact_success_sum_upper_bound_for_anticommuting_scaled_unitaries=1,
                records=records, changed_law_controls=controls, supplied_spinor_reference=spinor,
                scope='Guarded spatial half-turns with supplied coherent code/filter/reference; no two-core exchange',
                official_status='OW-07 untested at 0df73dc6; full local criterion open')


def acquisition_banks(data):
    for source in data['records']:
        for key in ('unfiltered_order', 'filtered_order', 'unfiltered_full_turn',
                    'filtered_full_turn', 'unfiltered_double_full_turn', 'filtered_double_full_turn'):
            yield source['panels'][key]
    for key in ('order', 'full_turn', 'double_full_turn'):
        yield data['supplied_spinor_reference']['interference'][key]
    for control in data['changed_law_controls']:
        for source in control['records']:
            for key in ('unfiltered_order', 'filtered_order', 'unfiltered_full_turn', 'filtered_full_turn'):
                yield source[key]


def verify_acquisition(actual, archived):
    # The probability-sorted sampler can interchange exactly degenerate bins
    # when a matrix library changes the last bit. Validate the independently
    # recomputed probabilities first, then replay each original random draw.
    rng = np.random.default_rng(SEED)
    for fresh, saved in zip(acquisition_banks(actual), acquisition_banks(archived), strict=True):
        compare(fresh['probabilities'], saved['probabilities'], 'detector probabilities')
        counts, complement = sample_acquisition(saved['probabilities'], rng)
        fresh['counts'] = counts.tolist()
        fresh['sampling_l1_complement_bound'] = complement
    compare(actual, archived)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--receipt', type=Path, default=HERE / 'results.json')
    args = parser.parse_args()
    result = run()
    if args.output:
        text = json.dumps(result, indent=2, allow_nan=False) + '\n'
        with args.output.open('x') as stream:
            stream.write(text)
    else:
        verify_acquisition(result, json.loads(args.receipt.read_text()))
        print('Replayed complete projected native half-turn histories and literal full-turn controls.')
