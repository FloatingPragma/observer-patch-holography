"""Attach source-compiled phase tomography to actual native C/E records.

The phase controller remains an explicit coherent service. Both detector
resets and the terminal pointer are executed by the literal native carrier.
All source export complements, reset timeouts and fine readouts are retained.
"""
from collections import defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[3]
WORKSPACE = BASE.parent
RER = WORKSPACE/'reverse-engineering-reality'
sys.path.insert(0, str(RER))
from model.observer import native
from model.observer.source_scalar import compiler
from model.observer.source_discrimination import circuit
from model.native_interface import pointer
from model.native_clock import configuration
from model.reference_recovery import recovery


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gauge():
    path = HERE.parent/'gauge_loop/gauge_loop.py'
    spec = importlib.util.spec_from_file_location('actual_row1_gauge_loop', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def complete_word(word, bits, initial):
    """Full input-column extension of the existing complete local executor.

    The source matrices, contraction and failure extraction are the existing
    ones. All input columns are contracted together, preserving cross-column
    coherences and arbitrary external references.
    """
    tape = compiler._validate(word, bits)
    state = initial.copy()
    gram = np.zeros((initial.shape[1],)*2, complex)
    first_failures = []
    for event, step in enumerate(tape):
        full = compiler.native_pulse(step['kind'], step['angle'])
        targets = step['targets']
        failures = recovery._failure_vectors(state, full, targets, bits)
        for complement, branch in enumerate(failures):
            effect = branch.conj().T@branch
            gram += effect
            first_failures.append({'event': event, 'complement': complement,
                                   'all_column_effect_trace': float(np.trace(effect).real)})
        dimension = 1 << len(targets)
        state = compiler._apply(state, full[:dimension, :dimension], targets, bits)
    residual = float(np.max(abs(state.conj().T@state+gram-initial.conj().T@initial)))
    assert residual < 2e-10
    return state, gram, first_failures, residual


def actual_pointer():
    """Fine literal readout after the actual six-pulse coupling."""
    initial = np.eye(8, dtype=complex)[:, [pointer.READY_CODE, pointer.READY_CODE+4]]
    executed, failure, slots, residual = complete_word(pointer._local_word(), 3, initial)
    result = {}
    for law in ('C', 'E'):
        grouped = defaultdict(lambda: defaultdict(lambda: np.zeros((2, 2), complex)))
        metadata = {}
        for z in range(4):
            source = circuit.decode(z)
            coupling = executed[[z, z+4]]
            for coin in (0, 1):
                output, (_, current) = native.native(source, pointer.SEAM, coin)
                label = current if law == 'C' else source[0]-source[1]
                key = coin, label
                bit = int(2*current+2*coin-1 < 0)
                grouped[key][output] += coupling/math.sqrt(2)
                metadata[key] = {'law': law, 'coin': coin, 'label': label,
                                 'current': current, 'decoded_pointer_bit': bit}
        result[law] = [(metadata[key], list(environment.values()))
                       for key, environment in sorted(grouped.items())]
        effects = sum(k.conj().T@k for _, ks in result[law] for k in ks)
        assert np.max(abs(effects+failure-np.eye(2))) < 2e-10
    # Complete reduced superoperator, including arbitrary references.
    worst = 0.
    for law, outcomes in result.items():
        for bit in (0, 1):
            actual = sum(np.kron(k.conj(), k) for meta, ks in outcomes
                         if meta['decoded_pointer_bit'] == bit for k in ks)
            projector = np.diag([int(bit == 0), int(bit == 1)])
            worst = max(worst, float(np.max(abs(actual-np.kron(projector, projector)))))
    return result, {'source_pulses': len(pointer._local_word()),
                    'full_matrix_unit_reduced_instrument_error': worst,
                    'export_failure_effect_trace': float(np.trace(failure).real),
                    'export_first_failure_slots': slots, 'gram_error': residual,
                    'exact_reference_certificate': pointer.exact_instrument_audit()}


def reset_contract():
    """Actual finite reset of two initially nonblank native carriers."""
    initial = native.states()[0]
    entries = []
    success = Fraction(1)
    for code in (0, pointer.READY_CODE):
        target = circuit.decode(code)
        values = configuration.success_probabilities(4, 32, full=True, mode='reset', target=target)
        cert = configuration.certificate(4, 32, full=True, mode='reset', target=target)
        p = values[0]
        entries.append({'target_code': code, 'initial_native_loads': list(initial),
                        'target_native_loads': list(target), 'success_probability': str(p),
                        'certificate': cert,
                        'all_fine_histories': {'policy': 'model.native_clock.configuration.policy',
                            'executor': 'model.native_clock.configuration.execute',
                            'amplitude': '2^(-attempts/2) on each column of word_support',
                            'failure_retains': 'timeout reason, every seam/coin/current, final native carrier'},
                        'complete_input_success_probabilities': list(map(str, values))})
        success *= p
    return entries, success


def run():
    target, target_path = load_gauge()
    source_before = sha(target_path)
    outcomes, pointer_receipt = actual_pointer()
    resets, p_reset = reset_contract()
    angles = np.linspace(0., 2*math.pi, 33)
    source_charge = np.array([mask.bit_count()-3 for mask in range(64)])
    initial = np.zeros((128, 64), complex)
    initial[:64] = np.eye(64)
    contrasts = {law: {q: np.zeros(33, complex) for q in range(-3, 4)} for law in ('C', 'E')}
    fine_rows = []
    export_rows = []
    max_error = max_gram = max_failure = 0.
    for angle_index, theta in enumerate(angles):
        for axis in ('X', 'Y'):
            prepare = [compiler.pulse('H', math.pi/2, [6])]
            rotate = [] if axis == 'X' else compiler.remap(compiler.compile_diagonal([0., -math.pi/2]), [6])
            word = prepare+target.phase_word(float(theta))+rotate+compiler.inverse(prepare)
            matrix, fail, slots, gram_error = complete_word(word, 7, initial)
            max_gram = max(max_gram, gram_error)
            max_failure = max(max_failure, float(np.linalg.norm(fail, 2)))
            export_rows.append({'angle_index': angle_index, 'axis': axis,
                                'pulses': len(word), 'first_failure_slots': slots,
                                'complete_failure_effect_trace': float(np.trace(fail).real)})
            # Trace the final logical pointer and auxiliary environment only.
            # Source outputs remain complete subnormalized 64-column operators.
            for law, native_records in outcomes.items():
                probabilities = {bit: np.zeros(64) for bit in (0, 1)}
                for meta, ks in native_records:
                    rates = np.zeros(64)
                    for k in ks:
                        continuing = np.kron(k, np.eye(64))@matrix
                        rates += np.sum(abs(continuing)**2, axis=0)
                    probabilities[meta['decoded_pointer_bit']] += rates
                    fine_rows.append({'angle_index': angle_index, 'axis': axis, **meta,
                                      'probability_by_occupation_mask_after_successful_reset': rates.tolist(),
                                      'unconditional_weight_multiplier': str(p_reset/2),
                                      'poststate_recipe': 'sum_env (k_env tensor I64) A_theta |input><input| A_theta* (k_env* tensor I64); encoded record fibre is unchanged'})
                expected = np.cos(source_charge*theta) if axis == 'X' else np.sin(source_charge*theta)
                contrast = probabilities[0]-probabilities[1]
                max_error = max(max_error, float(np.max(abs(contrast-expected))))
                for q in range(-3, 4):
                    contrasts[law][q][angle_index] += np.mean(contrast[source_charge == q])*(1 if axis == 'X' else 1j)
    degrees = {law: {str(q): target.degree(curve) for q, curve in families.items()}
               for law, families in contrasts.items()}
    assert all(degrees[law] == {str(q): q for q in range(-3, 4)} for law in degrees)
    assert max_error < 2e-10 and source_before == sha(target_path)
    # Exact one-to-one C10 encoding, with two unused fibre labels retained.
    with np.load(BASE/'next/C10/joint_source_states.npz', allow_pickle=False) as archived:
        encoded = [int(mask)+64*int(f[0]) for mask, f in zip(archived['occupation_masks'], archived['half_records'])]
        assert len(set(encoded)) == 384 and all(0 <= value < 384 for value in encoded)
    pins = [Path(__file__), target_path, RER/'model/native_interface/pointer.py',
            RER/'model/native_clock/configuration.py', RER/'model/observer/native.py',
            RER/'model/observer/source_scalar/compiler.py', RER/'model/reference_recovery/recovery.py',
            BASE/'next/C10/joint_source_states.npz']
    return {'schema': 'oph.row1.native_pointer_bridge.v1',
            'checks': {'actual_C_and_E_readout_closes_phase_loop': True,
                       'all_source_occupation_columns_and_reference_coherences_retained': True,
                       'export_failure_effects_retained': True,
                       'both_nonblank_auxiliary_native_resets_have_exact_complete_laws': True,
                       'C10_mask_fibre_encoding_is_injective': True,
                       'target_unchanged': True},
            'character_degrees_from_actual_native_records': degrees,
            'max_phase_contrast_error': max_error, 'max_complete_gram_error': max_gram,
            'max_meter_export_failure_operator_norm': max_failure,
            'pointer': pointer_receipt, 'reset': resets,
            'joint_reset_success_probability': str(p_reset),
            'joint_reset_timeout_probability': str(1-p_reset),
            'native_record_rows': fine_rows, 'source_export_rows': export_rows,
            'encoding': {'isometry': '|mask,F> -> |mask> tensor |F0>; F1..F5 are recovered by the actual Gauss equations',
                         'source_dimension': 384, 'padded_memory_dimension': 512,
                         'invalid_fibre_labels': [6, 7], 'invalid_complement': 'retained as rejected source admission, not charge states'},
            'closed_interfaces': ['logical pointer readout replaced by literal C/E coin/current instrument',
                                  'pointer and ready auxiliary initialized by finite native reset, with all timeouts retained',
                                  'all full source pulse export effects included; no success normalization'],
            'remaining_premises': ['coherent reference-preserving proper-code transfer and memory routing',
                                   'signed analog source-response pulses with supplied angle settings',
                                   'source-operation schedule pauses C10 between native evolution intervals; no unproved physical-time equivalence',
                                   'C10 CAR/gauge source is not derived here from committed repair',
                                   'gauge character winding has no demonstrated spatial-defect or physical electric identity'],
            'source_sha256': {str(p.relative_to(WORKSPACE)): sha(p) for p in pins},
            'full_criterion_met': False, 'official_rows_changed': []}


if __name__ == '__main__':
    result = run()
    (HERE/'NATIVE_POINTER_BRIDGE.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['checks', 'character_degrees_from_actual_native_records',
                                                'max_phase_contrast_error', 'max_complete_gram_error',
                                                'max_meter_export_failure_operator_norm']}))
