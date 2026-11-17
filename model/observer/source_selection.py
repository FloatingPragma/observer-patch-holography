"""Source-executed record-phase ambiguity with a fixed coherent read interface.

The fixed-seam C/E dilation and its record-phase variant induce the same
complete committed instrument. An X-basis read before commitment distinguishes
them; committing first removes that distinction. The new read codes are not
native currents. Every coin, read code and native/reference amplitude remains.
This checks a source-selection obligation, without selecting either rule.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from .source_discrimination import circuit
from .source_scalar import compiler
from .scalar_response.experiment import canonical, compare_receipt


ROOT = Path(__file__).resolve().parents[2]
SIZE = 512
TOLERANCE = 2e-10
SOURCE_PATHS = (
    'model/observer/source_selection.py', 'model/observer/native.py',
    'model/observer/source_discrimination/circuit.py',
    'model/observer/source_scalar/compiler.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json',
)


def _law(value):
    if type(value) is not str or value not in ('C', 'E'):
        raise ValueError('native instrument C or E required')
    return value


def _boolean(value):
    if type(value) is not bool:
        raise ValueError('an explicit Boolean experiment switch required')
    return value


def _original(value):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('amplitudes cannot be Boolean, masked or missing')
    if isinstance(value, (tuple, list)):
        for item in value:
            _original(item)
    elif isinstance(value, np.ndarray) and value.dtype.kind not in 'iufc':
        raise ValueError('original numeric amplitudes required')


def _state(value):
    _original(value)
    raw = np.asarray(value)
    if (raw.ndim not in (1, 2) or raw.shape[0] != SIZE
            or (raw.ndim == 2 and raw.shape[1] < 1)
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite 512-component amplitudes or reference columns required')
    result = np.array(raw, dtype=complex, copy=True).reshape(SIZE, -1)
    if not np.isfinite(result).all():
        raise ValueError('amplitudes must remain finite during conversion')
    return result


def _immutable(value):
    raw = np.ascontiguousarray(value, dtype=complex)
    if not np.isfinite(raw).all():
        raise ArithmeticError('nonfinite source output')
    return np.frombuffer(raw.tobytes(), dtype=complex).reshape(raw.shape)


def preparation_word():
    """Prepare (|1>+|2>)/sqrt(2) from blank native code, with its exact phase."""
    return ([compiler.pulse('Y', -math.pi/4, [0])]+compiler.cnot(0, 1)
            +[compiler.pulse('X', math.pi/2, [0]),
              compiler.pulse('I', -math.pi/2, [0])])


def phase_word(twisted):
    """Phase (-1)**low_record_bit; zero-angle controls retain both services."""
    angle = math.pi/2 if _boolean(twisted) else 0.
    return [compiler.pulse('Z', angle, [3]), compiler.pulse('I', -angle, [3])]


def read_word():
    return [compiler.pulse('Y', math.pi/4, [3])]


def prepared_state():
    blank = np.zeros(SIZE)
    blank[0] = 1
    return _immutable(compiler.execute(preparation_word(), 9, blank))


@lru_cache(maxsize=1)
def _coin():
    return _immutable(compiler.execute(
        [compiler.pulse('Y', -math.pi/4, [0])], 1, [1., 0.]))


@lru_cache(maxsize=8)
def _repair(law, coin, twisted):
    return _immutable(compiler.execute(
        circuit.repair_word(law, coin)+phase_word(twisted), 3+circuit.RECORD_BITS[law]))


@lru_cache(maxsize=2)
def _reader(law):
    return _immutable(compiler.execute(read_word(), 3+circuit.RECORD_BITS[law]))


def _local_inputs(value, records):
    # Literal spectator tensor factors, including every native padding row.
    local = value.reshape(64, 8, -1).transpose(1, 0, 2).reshape(8, -1)
    blank = np.zeros((8*records, local.shape[1]), complex)
    blank[:8] = local
    return blank


def _native_output(local, columns):
    return _immutable(local.reshape(8, 64, columns).transpose(1, 0, 2).reshape(SIZE, columns))


def committed_branches(law, state, *, twisted=False):
    """Baseline computational-record instrument, including reserved raw codes.

    Full native/reference arrays may be unnormalized. Scalar record phases
    cancel in each CP block; no assertion identifies their coherent lifts.
    """
    law, twisted, value = _law(law), _boolean(twisted), _state(state)
    records = 1 << circuit.RECORD_BITS[law]
    blank = _local_inputs(value, records)
    result = []
    for coin in (0, 1):
        blocks = (_repair(law, coin, twisted)@blank).reshape(records, 8, -1)*_coin()[coin]
        for raw, block in enumerate(blocks):
            result.append({'coin': coin, 'raw_record': raw,
                           'abort': raw not in circuit.record_values(law, coin),
                           'amplitude': _native_output(block, value.shape[1])})
    return tuple(result)


def read_branches(law, state, *, twisted=False, commit_before_read=False):
    """Complete output after X reading the low record bit and Z reading others.

    The low read bit zero denotes the plus outcome. These read codes have no
    native-current/abort interpretation. In the dephased control only, the
    earlier computational record is separately retained; every raw value,
    including reserved values, receives the specified later read operation.
    This comparator is not the native observer's halt-on-abort program.
    """
    law, twisted = _law(law), _boolean(twisted)
    committed, value = _boolean(commit_before_read), _state(state)
    records = 1 << circuit.RECORD_BITS[law]
    blank = _local_inputs(value, records)
    result = []
    for coin in (0, 1):
        coherent = _repair(law, coin, twisted)@blank*_coin()[coin]
        for old in (range(records) if committed else (None,)):
            incoming = coherent
            if old is not None:
                incoming = np.zeros_like(coherent)
                incoming[8*old:8*(old+1)] = coherent[8*old:8*(old+1)]
            blocks = (_reader(law)@incoming).reshape(records, 8, -1)
            for code, block in enumerate(blocks):
                result.append({'coin': coin, 'committed_record': old, 'read_code': code,
                               'amplitude': _native_output(block, value.shape[1])})
    return tuple(result)


def resources(law):
    law = _law(law)
    variants = {}
    counts = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')
    for coin in (0, 1):
        for twisted in (False, True):
            word = ([compiler.pulse('Y', -math.pi/4, [0])]+preparation_word()
                    +circuit.repair_word(law, coin)+phase_word(twisted)+read_word())
            variants[f'coin_{coin}_twisted_{str(twisted).lower()}'] = compiler.resources(word)
    first = next(iter(variants.values()))
    if any(any(row[key] != first[key] for key in counts) for row in variants.values()):
        raise ArithmeticError('declared equal-pulse source schedules differ')
    return {**{key: first[key] for key in counts}, 'source_word_variants': variants,
            'native_data_qubits': 9, 'quantum_record_qubits': circuit.RECORD_BITS[law],
            'native_blank_qubits_in_prepared_witness': 9,
            'native_register_preparations': 1, 'record_register_preparations': 1,
            'fresh_fair_coin_qubits': 1, 'new_read_pointer_qubits': 0,
            'coin_initializations_and_reads': 1, 'final_record_register_reads': 1,
            'six_level_processor_additional': True,
            'coherent_final_classical_bits': 1+circuit.RECORD_BITS[law],
            'dephased_control_final_classical_bits': 1+2*circuit.RECORD_BITS[law],
            'dephased_control_extra_record_measurements': 1,
            'physical_initialization_storage_readout_costs_included': False}


def _committed_audit(law):
    original = committed_branches(law, np.eye(SIZE))
    twisted = committed_branches(law, np.eye(SIZE), twisted=True)
    errors, gram = [], np.zeros((SIZE, SIZE), complex)
    for left, right in zip(original, twisted, strict=True):
        if {k: v for k, v in left.items() if k != 'amplitude'} != {
                k: v for k, v in right.items() if k != 'amplitude'}:
            raise ArithmeticError('phase variant changed a committed record')
        sign = (-1)**(left['raw_record'] & 1)
        errors.append(float(np.linalg.norm(right['amplitude']-sign*left['amplitude'], 'fro')))
        gram += left['amplitude'].conj().T@left['amplitude']
    completeness = float(np.linalg.norm(gram-np.eye(SIZE), 'fro'))
    if max(*errors, completeness) > TOLERANCE:
        raise ArithmeticError('complete source committed-map phase identity failed')
    return {'all_input_columns': SIZE, 'raw_outcomes_including_reserved': len(original),
            'maximum_scalar_phase_operator_error_diagnostic': max(errors),
            'completeness_frobenius_error_diagnostic': completeness,
            'arbitrary_reference_identity': 'K_twisted[c,r]=(-1)**(r mod 2) K_original[c,r]',
            'diagnostics_are_not_outward_roundoff_certificates': True}


def run():
    initial = prepared_state()
    laws = {}
    for law in ('C', 'E'):
        runs = []
        for committed in (False, True):
            for twisted in (False, True):
                branches = read_branches(law, initial, twisted=twisted, commit_before_read=committed)
                rows = [{**{key: value for key, value in branch.items() if key != 'amplitude'},
                         'probability': float(np.vdot(branch['amplitude'], branch['amplitude']).real)}
                        for branch in branches]
                mass = math.fsum(row['probability'] for row in rows)
                plus = math.fsum(row['probability'] for row in rows if row['read_code'] % 2 == 0)
                expected = .5 if committed else .25 if twisted else .75
                if (not math.isfinite(mass) or not math.isfinite(plus)
                        or abs(mass-1) > TOLERANCE or abs(plus-expected) > TOLERANCE):
                    raise ArithmeticError('source phase/interface contrast failed')
                runs.append({'commit_before_read': committed, 'twisted': twisted,
                             'probability_mass_diagnostic': mass, 'plus_probability_diagnostic': plus,
                             'ideal_plus_probability_exact': '1/2' if committed else '1/4' if twisted else '3/4',
                             'complete_read_outcomes': rows})
        laws[law] = {'complete_committed_map_audit': _committed_audit(law),
                     'runs': runs, 'resources_per_attempt': resources(law)}
    return {'schema': 'oph.observer.source_selection.v1',
            'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'laws': laws,
            'contract': {'seam': [0, 1], 'preparation': '(|native code 1>+|native code 2>)/sqrt(2)',
                         'phase': 'source Z(pi/2) and I(-pi/2) on the low unread record bit; both angles zero in control',
                         'read': 'source Y(pi/4) then complete record-register computational read',
                         'read_dictionary': 'low bit 0 is X-plus; higher bits are unchanged Z reads, not native currents',
                         'outputs': 'every branch API returns all 512 unnormalized native rows and reference columns'},
            'scope': {'same_complete_committed_CE_instrument': True,
                      'same_coherent_record_instrument': False,
                      'original_current_retained_in_coherent_read': False,
                      'source_primitives_and_code_services_unchanged': True,
                      'phase_selected_by_this_experiment': False,
                      'coherent_access_required_for_distinction': True,
                      'physical_calibration_or_empirical_comparison': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'source_record_phase_check': True, 'selected_repair_rule': False}))


if __name__ == '__main__':
    main()
