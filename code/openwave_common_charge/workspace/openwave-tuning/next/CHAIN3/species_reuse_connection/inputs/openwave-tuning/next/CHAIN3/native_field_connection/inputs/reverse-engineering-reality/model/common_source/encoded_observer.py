"""Active native quantum memory in the canonical continued observer.

The native mode is part of quantum_native.JointHistory: the work program,
record menu, output APIs and complete-map auditor are shared. This module
supplies the encoded protocol's resource ledger and reproducible controls.
Memory is assessed before the optional destructive terminal logical read.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import continuation_bounds as bounds, encoded_memory as memory
from . import native_routing as routing, protocol, quantum_native as q


def passive_obstruction(dimension):
    """Exact uniform lower bounds proved in ENCODED_MEMORY.md, not fitted.

    Q_2 has eigenvalues +/-sqrt(1/2). Q_4 has eigenvalues
    +/-sqrt((3+/-sqrt(6))/2). Their distance from odd pi/2 gives the
    stated strictly positive minimum of (I+sin Q_d)/2. The passive
    record-selected 01/23 experiment therefore cannot preserve the qubit.
    """
    if type(dimension) is not int or dimension not in (2, 4):
        raise ValueError('the proved finite pointer cutoffs are 2 and 4')
    lower = F(1, 8) if dimension == 2 else F(1, 1536)
    return {'uniform_second_smooth_plus_lower_exact': str(lower),
            'uniform_passive_worst_case_decoder_error_lower_exact': str(lower/2),
            'scope': 'ideal finite pointer; arbitrary incoming work state; before terminal read',
            'infinite_occupation_uniform_lower_bound_claimed': False}


def _encoded(history):
    q._validate_menu(history)
    if history.native_mode not in ('protected_qubit', 'passive_qubit'):
        raise ValueError('the declared encoded observer mode is required')


def phase_distance(history):
    """Complete CQ output distance for logical +/- with declared work input.

    Each native difference is A X A*, where A has two logical columns.
    A thin QR reduces its nonzero spectrum to a 2x2 Hermitian matrix;
    there is no near-equal probability subtraction or code projection.
    Every raw branch, including aborts and numerical padding, contributes.
    This floating evaluation is a diagnostic, not an outward certificate.
    """
    _encoded(history)
    flip = np.array([[0., 1.], [1., 0.]])
    terms, plus, minus = [], [], []
    for branch in history.branches:
        a, v = branch.native_map, branch.work_state
        weight = float(np.vdot(v, v).real)
        _, triangular = np.linalg.qr(a, mode='reduced')
        difference = triangular@flip@triangular.conj().T
        terms.append(weight*float(np.abs(np.linalg.eigvalsh(difference)).sum())/2)
        plus.append(weight*float(np.vdot(a[:, 0]+a[:, 1], a[:, 0]+a[:, 1]).real)/2)
        minus.append(weight*float(np.vdot(a[:, 0]-a[:, 1], a[:, 0]-a[:, 1]).real)/2)
    result = math.fsum(terms)
    if (not math.isfinite(result) or not 0 <= result <= 1+q.TOLERANCE
            or any(not math.isfinite(total) or abs(total-1) > q.TOLERANCE
                   for total in (math.fsum(plus), math.fsum(minus)))):
        raise ArithmeticError('phase comparison requires complete finite normalized outputs')
    return result


def logical_entanglement_fidelity(history):
    """Bell-input fidelity with the retained encoded qubit, before readout."""
    _encoded(history)
    bell = np.eye(2, dtype=complex)/math.sqrt(2)
    expected = q._logical_code()@bell
    terms, masses = [], []
    for a, v in q.output(history, bell).values():
        weight = float(np.vdot(v, v).real)
        terms.append(weight*abs(np.vdot(expected, a))**2)
        masses.append(weight*float(np.vdot(a, a).real))
    result, mass = math.fsum(terms), math.fsum(masses)
    if (not math.isfinite(result) or not 0 <= result <= 1+q.TOLERANCE
            or not math.isfinite(mass) or abs(mass-1) > q.TOLERANCE):
        raise ArithmeticError('logical/reference comparison discarded an output')
    return result


def memory_witness(active, passive):
    """Compare complete corrected and uncorrected observer outputs."""
    _encoded(active)
    _encoded(passive)
    if (active.native_mode != 'protected_qubit' or passive.native_mode != 'passive_qubit'
            or (active.law, active.dimension, active.steps, active.field_backend)
            != (passive.law, passive.dimension, passive.steps, passive.field_backend)
            or active.field_backend != 'source'):
        raise ValueError('matched source-executed active and passive observers required')
    targets = q._programs(active.dimension, active.steps, 'product')
    q_bar = math.fsum(float(weight)*math.fsum(
        branch.probability for branch in targets[q._record_current(active.law, coin, raw)].branches
        if branch.record[2] == 0)
        for (coin, raw), weight in q._code_weights(active.law, (0, 1)).items())
    distances = {'active': phase_distance(active), 'passive': phase_distance(passive)}
    fidelity = logical_entanglement_fidelity(active)
    if (abs(distances['active']-1) > q.TOLERANCE
            or abs(distances['passive']-(1-q_bar)) > q.TOLERANCE
            or abs(fidelity-1) > q.TOLERANCE):
        raise ArithmeticError('source memory differs from the independent logical channel')
    terminal = {setting: q.readout(active, q.prepared_native(setting, native_mode='protected_qubit'))
                for setting in q.circuit.SETTINGS[:3]}
    for setting, expected in zip(q.circuit.SETTINGS[:3], (1., 0., .5), strict=True):
        if abs(terminal[setting]['plus_code_zero_probability']-expected) > q.TOLERANCE:
            raise ArithmeticError('decoded logical read does not preserve its source phase')
    return {'prepared_second_smooth_plus_probability_diagnostic': q_bar,
            'complete_pre_read_logical_phase_distances_diagnostic': distances,
            'active_encoded_bell_fidelity_diagnostic': fidelity,
            'passive_distance_formula': '1 - q_bar',
            'passive_any_decoder_worst_case_error_lower_formula': 'q_bar / 2',
            'evaluated_passive_decoder_error_lower_diagnostic': q_bar/2,
            'exact_uniform_passive_obstruction': passive_obstruction(active.dimension),
            'source_decoded_complete_terminal_reads': terminal,
            'terminal_read_is_destructive': True,
            'logical_memory_is_assessed_before_terminal_read': True}


def resources(law, dimension, steps, *, native_mode='protected_qubit'):
    """Enumerate compatible realized-path schedules, including early halts.

    Each full-register source permutation pays for every macro invocation;
    cached host contractions do not reduce the physical source-word cost.
    These are scheduled upper bounds including raw zero-probability paths.
    """
    q._law(law)
    if q._mode(native_mode) == 'full_carrier':
        raise ValueError('an encoded native mode required')
    histories = q._programs(dimension, steps, 'source')
    protected = native_mode == 'protected_qubit'
    encoder, decoder = memory.encoding_resources(), memory.decoding_resources()
    first = routing.routed_repair_resources(law, (0, 1))
    last = {bit: routing.routed_repair_resources(law, q.final_seam(native_mode, bit))
            for bit in (0, 1)}
    prefix = 3+encoder['native_pulses']+first['native_pulses']
    candidates = [(prefix+dict(histories[0].preparation_resources)['native_pulses'],
                   {'halt': 'first_raw_abort'})]
    record_count = 1 << q.circuit.RECORD_BITS[law]
    # Compile each declared recovery schedule once on the host. Each realized
    # invocation still pays the complete source-word cost in its path below.
    recovery_cost = {(seam, coin, raw): memory.recovery_resources(
        law, seam, coin, raw)['native_pulses']
        for seam in ((0, 1), (2, 3)) for coin in (0, 1)
        for raw in range(record_count)
        if protected and q._record_current(law, coin, raw) is not None}
    recovery_max = [0, 0]
    for coin1 in (0, 1):
        for raw1 in range(record_count):
            current = q._record_current(law, coin1, raw1)
            if current is None:
                continue
            first_recovery = recovery_cost.get(((0, 1), coin1, raw1), 0)
            recovery_max[0] = max(recovery_max[0], first_recovery)
            for work in histories[current].branches:
                seam = q.final_seam(native_mode, work.record[2])
                cost = prefix+first_recovery+dict(work.resources)['native_pulses']+last[work.record[2]]['native_pulses']
                for coin2 in (0, 1):
                    for raw2 in range(record_count):
                        aborted = q._record_current(law, coin2, raw2) is None
                        recover = recovery_cost.get((seam, coin2, raw2), 0)
                        recovery_max[1] = max(recovery_max[1], recover)
                        total = cost+recover+(0 if aborted else decoder['native_pulses']+1)
                        candidates.append((total, {'first_coin': coin1, 'first_raw': raw1,
                            'work_record': list(work.record), 'final_seam': list(seam),
                            'last_coin': coin2, 'last_raw': raw2, 'last_abort': aborted}))
    maximum, witness = max(candidates, key=lambda pair: pair[0])
    return {'maximum_native_pulses_per_attempt': maximum,
            'maximum_proper_code_load_events_per_attempt': maximum,
            'maximum_proper_code_export_events_per_attempt': maximum,
            'maximum_proper_code_service_events_per_attempt': 2*maximum,
            'maximizing_scheduled_record': witness,
            'logical_input_preparation_pulses': 3,
            'encoder': encoder, 'terminal_decoder': decoder,
            'maximum_first_recovery_native_pulses': recovery_max[0],
            'maximum_last_recovery_native_pulses': recovery_max[1],
            'native_input_qubits': 9, 'supplied_logical_qubits': 1,
            'initial_blank_native_code_qubits': 8,
            'blank_quantum_qubits_added_during_recovery': 0,
            'archive_record_rewrites_during_recovery': 0,
            'native_coin_initializations': 2, 'native_record_initializations': 2,
            'native_classical_readout_registers': 5,
            'field_preparations': 1, 'field_pointer_initializations_and_reads': 4,
            'peak_code_qubits_excluding_six_level_processor': 9+2+2*(dimension.bit_length()-1)+q.circuit.RECORD_BITS[law],
            'fixed_width_terminal_record_bits_including_two_status_bits': 2*(1+q.circuit.RECORD_BITS[law])+4+9+2,
            'source_word_counts_include_all_repeated_invocations': True,
            'physical_services_or_host_tree_storage_counted': False}


def required_service_keys(program):
    return (*bounds.required_service_keys(program.operations), 'logical_input_preparation',
            'native_encoder', 'native_logical_storage_through_work',
            'first_native_coin_and_recorded_repair', 'first_native_recovery',
            'last_native_coin_and_recorded_routed_repair', 'last_native_recovery',
            'native_terminal_decoder_rotation_and_readout')


def run(steps=q.STEPS):
    rows = []
    for dimension in (2, 4):
        program = protocol.continuation_program(dimension, steps)
        budget = bounds.program_bound(program)
        laws = {}
        for law in ('C', 'E'):
            histories = {mode: q.execute(law, dimension, steps, native_mode=mode)
                         for mode in ('protected_qubit', 'passive_qubit')}
            laws[law] = {'complete_channel_audits': {mode: q.audit(history) for mode, history in histories.items()},
                         'memory_witness': memory_witness(histories['protected_qubit'], histories['passive_qubit']),
                         'resources': {mode: resources(law, dimension, steps, native_mode=mode) for mode in histories}}
        rows.append({'dimension': dimension, 'laws': laws, 'ideal_complete_work_output_budget': budget,
                     'physical_service_coverage': bounds.service_error_budget(
                         F(budget['complete_ideal_product_to_infinite_instrument_error_upper']),
                         required_stage_keys=required_service_keys(program))})
    return {'schema': 'oph.common_source.encoded_observer.v1', 'steps': steps, 'cutoffs': rows,
            'scope': {'canonical_runtime': 'quantum_native.JointHistory with explicit native_mode',
                      'logical_input_dimension': 2, 'native_output_rows_retained': 512,
                      'encoded_interventions_added_to_native_law': True,
                      'arbitrary_erased_full_carrier_coherence_recovered': False,
                      'logical_memory_requires_scalar_continuum_limit': False,
                      'unknown_quantum_memory_reconstructed_from_classical_tape': False,
                      'record_fault_tolerance_or_physical_service_accuracy_established': False,
                      'physical_clock_calibration_established': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(q.canonical(result))
    if args.check:
        q.compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'source_corrected_native_memory': True, 'cutoffs': [2, 4]}))


if __name__ == '__main__':
    main()
