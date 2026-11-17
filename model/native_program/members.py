"""Literal admission of retained native observers and their counted services.

The long record-clock receipt is prior execution evidence, checked against
its source hashes. Its native schedule and generated coherent words are
admitted afresh. No readout, clock success, physical time or target field is
substituted for an operation in the source grammar.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

from tools import strict_json
from . import interpreter as p
from ..local_repair_family import acquisition, inference, learned_flow, record_clock, modular_clock

ROOT = Path(__file__).resolve().parents[2]
RECORD_RECEIPT = 'model/local_repair_family/record_clock_results.json'


def _fixed(edge, law='mixture', repeat=1):
    return {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(edge)},
            'law': law, 'repeat': repeat}


def checked_record_receipt():
    receipt = strict_json.load(ROOT/RECORD_RECEIPT)
    record_clock.contract()
    if receipt['schema'] != 'oph.local_repair_family.record_clock.v1':
        raise ValueError('the retained record-clock schema changed')
    for name, expected in receipt['source_sha256'].items():
        source = (ROOT/name).resolve()
        if not source.is_relative_to(ROOT) or not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise ValueError('retained record-clock source binding changed: '+name)
    return receipt


def record_clock_member():
    """Fresh grammar admission and resource recomputation for one named member.

    This constructor has no caller-selected function or state. Inference
    consumes only complete retained seam/coin counts; the simulator weights
    enter only the declared native apparatus. The existing ``learned_flow``
    encoding/proper-code service is explicit, not a generic tensor map.
    """
    retained = checked_record_receipt()
    counts = {tuple(row['seam']): (row['trials'], row['coin1_count']) for row in retained['acquisition']['blocks']}
    inferred = inference.infer_counts(counts)
    if inferred['status'] != 'accept':
        return {'status': 'STOP', 'reason': inferred['stop_reason'], 'native_source_unchanged': True}
    training = p.program(commands=[_fixed(edge, 'alternating', acquisition.SAMPLES_PER_SEAM) for edge in p.local.SEAMS])
    preparation = p.program(commands=[{'op': 'repair', 'selector': {'kind': 'uniform'},
                                      'law': 'mixture', 'repeat': record_clock.PREPARATION_ATTEMPTS}])
    continuations = {law: p.program(commands=[_fixed(edge, law) for edge in
                     (*modular_clock.REQUESTED_SEAMS, modular_clock.PROBE)]) for law in ('C', 'E')}
    source = learned_flow.source_program(inferred['relative_weights'])
    words = {name: p.source_word(word, 9) for name, word in source.items()}
    source_pulses = len(source['probe'])+max(2*len(source['half_clock']), len(source['clock']))+len(source['analyzer'])
    native_cost = acquisition.native_word_costs()
    training_pulses = sum(sum(row.values())*acquisition.SAMPLES_PER_SEAM//2 for row in native_cost.values())
    preparation_pulses = record_clock.PREPARATION_ATTEMPTS*max(v for row in native_cost.values() for v in row.values())
    future_cost = acquisition.native_word_costs((*modular_clock.REQUESTED_SEAMS, modular_clock.PROBE))
    future_pulses = max(sum(row[law] for row in future_cost.values()) for law in ('C', 'E'))
    total_attempts = (training['declared_resources']['attempts_upper']+
                      preparation['declared_resources']['attempts_upper']+3)
    total_pulses = training_pulses+preparation_pulses+future_pulses+source_pulses
    if (total_attempts != retained['resources']['total_native_attempts']
            or total_pulses != retained['resources']['complete_supplied_realization_pulses_upper']
            or source_pulses != retained['resources']['coherent_clock_pulses']
            or total_pulses > 2**48 or total_attempts > 2**32):
        raise ArithmeticError('literal schedule/service counts disagree with retained execution')
    # Bind both complete readout laws and the actual native-to-code injection.
    states, codes, pairs = learned_flow._source_chart()
    if len(states) != 66 or len(set(codes)) != 66 or len(pairs) != 10:
        raise ArithmeticError('the named member lost its explicit proper-code embedding')
    terminal = retained['acquisition']['private_replay']['terminal_state_index']
    terminal_state = retained['acquisition']['private_replay']['terminal_state']
    if list(states[terminal]) != terminal_state:
        raise ArithmeticError('retained terminal source disagrees with the shared source chart')
    if (retained['acquisition']['resources']['source_resets_during_training'] != 0
            or retained['acquisition']['resources']['FIFO_capacity'] != 1
            or retained['resources']['training_count_storage_bits'] != 60*acquisition.SAMPLES_PER_SEAM.bit_length()
            or retained['acquisition']['resources']['exported_records'] != 30*acquisition.SAMPLES_PER_SEAM-1):
        raise ArithmeticError('continuing source/FIFO/count-register custody changed')
    return {'status': 'ADMITTED', 'prior_execution_receipt': RECORD_RECEIPT,
            'prior_execution_sha256': hashlib.sha256((ROOT/RECORD_RECEIPT).read_bytes()).hexdigest(),
            'training_program': training, 'preparation_program': preparation,
            'native_readout_programs': continuations,
            'source_words': {name: {'digest': p.digest(word), 'resources': word['declared_resources']}
                             for name, word in words.items()},
            'source_word_literal_payloads_admitted': True,
            'source_encoding': {'source_dimension': 66, 'work_dimension': 512,
                                'basis_addresses': list(codes), 'fixed_probe_pairs': len(pairs),
                                'owner': 'model/local_repair_family/learned_flow.py'},
            'composition_order': ['native acquisition', 'record-only inference or retained stop',
                                  '1600 native preparation attempts on the terminal source', 'probe',
                                  'native repair', 'half_clock', 'native repair', 'half_clock',
                                  'analyzer', 'complete native detector'],
            'composition_custody': {
                'executor': 'model/local_repair_family/record_clock.py: report and continue_from_records',
                'preparation': 'continued_source_density takes the actual acquisition terminal index; no source reset',
                'native_terminal_index': terminal,
                'native_terminal_state': terminal_state,
                'clock_input': 'the complete prepared density is passed to learned_flow.execute on the same H66 ordering',
                'literal_word_owner': 'learned_flow.source_program receives record-inferred relative weights only',
                'source_chart_injection': 'circuit.encode on each of the 66 states; full512 unitary word then complete source-address export',
                'training_archive': 'all emitted records; prior receipt stream hash is evidence, not a replacement for storage',
                'training_FIFO_slot': retained['acquisition']['retained_FIFO_slot'],
                'training_exports': retained['acquisition']['resources']['exported_records'],
                'count_register_inputs': 'complete trials/coin-one counts, independently checked against full raw histograms',
                'preparation_archive_comparison': 'exports and final FIFO traced only for the declared product-reference comparison',
                'success_conditioning': False,
                'long_composition_execution_evidence': 'the source-bound prior receipt; this constructor freshly validates membership and accounting'},
            'resources': {'total_native_attempts': total_attempts,
                          'complete_supplied_realization_pulses_upper': total_pulses,
                          'complete_supplied_realization_load_export_calls_upper': 2*total_pulses,
                          'clock_pulses_per_realization_upper': source_pulses,
                          'record_storage_bits_upper': (total_attempts+1)*40,
                          'FIFO_slots': 1, 'FIFO_bits_upper': 40,
                          'export_archive_bits_upper': (total_attempts-1)*40,
                          'terminal_status_bits_upper': 40,
                          'literal_program_bytes_stored': sum(v['declared_resources']['program_description_bytes']
                              for v in (training, preparation, *continuations.values(), *words.values())),
                          'native_controller_bits': 1,
                          'global_event_counter_bits': total_attempts.bit_length(),
                          'source_status_bits': 2,
                          'source_position_counter_bits': source_pulses.bit_length(),
                          'processor_complement_index_bits': 3,
                          'address_register_qubits': 9,
                          'initial_address_register_preparations': 1,
                          'initial_six_level_processor_preparations_upper': total_pulses,
                          'count_register_bits': 60*acquisition.SAMPLES_PER_SEAM.bit_length(),
                          'estimated_rational_weights_storage_bits': sum(
                              w.numerator.bit_length()+w.denominator.bit_length() for w in inferred['relative_weights'])},
            'retained_native_training_events_executed': retained['acquisition']['protocol']['total_events'],
            'coherent_realization_pulses_reexecuted_here': False,
            'complete_output_semantics': {
                'native_raw_zero_slots': 'complete typed native instrument',
                'address_export_aborts': 'retained by learned_flow.execute',
                'six_level_complement_first_aborts': 'complete export channel; retained proof bound in learned_flow.source_budget',
                'all_trillion_pulse_abort_matrices_materialized': False,
                'untouched_reference': 'full operator identity/certified complete-channel bound',
                'physical_apparatus_error_supplied': True},
            'clock_parameter': 'dimensionless modular parameter s',
            'source_selected_physical_clock_or_field_claimed': False}
