"""A source-coupled work-bit instrument read through native C/E currents.

The ready auxiliary is the native configuration of source code 2.  Its
preparation is a separate, complete native configuration-reset instrument;
this module never supplies that pure state by an unreported input reset.
On a conclusive reset its purity factors it from every work/reference state.

Two existing source CNOTs copy a chosen work bit into auxiliary orientation,
then the committed native coin/current pointer supplies the outcome.  Only
three source bits are active.  All four auxiliary codes reached even by the
raw numerical word, 0..3 with spectator label zero, are real full-carrier
configurations.  No auxiliary membership measurement or padding extension
is needed.  Arbitrary work padding and untouched references are retained.

The ideal reduced instrument is exactly the two Lueders maps P_b rho P_b.
The raw source word is certified separately against its full 8D permutation
target, without projection or normalization.  Existing proper-code
load/export, coherent two-bit routing, and control precision remain supplied;
the logical pointer blank and its projective read are the interfaces closed
by native reset and this construction, not by assuming a new ideal detector.

For initially correlated auxiliary/work/reference, a successful reset can
steer the work marginal.  The exact read acts on that POST-RESET marginal.
After discarding private reset records, success and failure marginals sum
to the original work state.  A uniform reset timeout t therefore contributes
at most t to the complete original-work comparison, with timeout retained.
The raw-word error adds on the successful subnormalized state, giving
t + delta_word; no conditional-success independence is claimed.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from model.compiled_native_interface import construction as compiled
from model.observer import native
from model.observer.scalar_response.experiment import canonical, compare_receipt as _portable_compare


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('pointer_results.json')
WORD_ERROR_GATE = F(1, 10**11)
DENSITY_ROUNDOFF_GATE = F(1, 10**13)
SOURCE_PATHS = tuple(sorted(set((*compiled.SOURCE_PATHS,
    'model/native_interface/pointer.py', 'model/observer/scalar_response/experiment.py',
    'tools/strict_json.py'))))
SEAM = (0, 1)
READY_CODE = 2


def ready_state():
    """Target for the separate native reset; no state is prepared by this call."""
    return compiled.circuit.decode(READY_CODE)


def _law(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('the committed native C or E instrument is required')
    return law


def _address(bit, work_qubits):
    if (type(work_qubits) is not int or not 1 <= work_qubits <= 10
            or type(bit) is not int or not 0 <= bit < work_qubits):
        raise ValueError('one to ten work bits and an in-range bit address required')
    return bit, work_qubits


def _local_word():
    # Local bits 0,1 are the auxiliary seam code; bit2 is the work control.
    # Starting at code2, a control one gives 2 -> 3 -> 1.
    return compiled.compiler.cnot(2, 0)+compiled.compiler.cnot(0, 1)


def source_word(bit=7, work_qubits=9):
    """Six source pulses; auxiliary source bits follow the work register.

    Seven other auxiliary encoding bits stay zero.  This is a factorized
    execution of a fixed wire-address lift, not a dense two-carrier matrix.
    """
    bit, work_qubits = _address(bit, work_qubits)
    return compiled.compiler.remap(_local_word(), [work_qubits, work_qubits+1, bit])


@lru_cache(maxsize=1)
def _raw_payload():
    matrix = compiled.compiler.execute(_local_word(), 3)
    return matrix.dtype.str, matrix.tobytes()


def local_matrix():
    dtype, payload = _raw_payload()
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(8, 8)


def target():
    """Independent complete integer permutation target for the emitted word."""
    result = np.zeros((8, 8), dtype=int)
    for source in range(8):
        a, b, control = source & 1, (source >> 1) & 1, source >> 2
        a ^= control
        b ^= a
        result[a+2*b+4*control, source] = 1
    return result


@lru_cache(maxsize=2)
def _effect_counts(law):
    """Exact native coarse-pointer effects on the four actually used codes."""
    states = tuple(compiled.circuit.decode(code) for code in range(4))
    groups, outcomes = defaultdict(list), {}
    for code, state in enumerate(states):
        for coin in (0, 1):
            output, record = native.native(state, SEAM, coin)
            if record != (coin, state[0]-output[0]):
                raise ArithmeticError('native current/coin grammar changed')
            current = record[1]
            label = current if law == 'C' else state[0]-state[1]
            key = coin, label
            bit = int(2*current+2*coin-1 < 0)
            if key in outcomes and outcomes[key] != bit:
                raise ArithmeticError('the native record does not determine its pointer bit')
            outcomes[key] = bit
            groups[key].append((code, output))
    effects = [[[0]*4 for _ in range(4)] for _ in (0, 1)]
    for key, mapping in groups.items():
        if len({output for _, output in mapping}) != len(mapping):
            raise ArithmeticError('a native recorded branch is not a partial injection')
        for left, out_left in mapping:
            for right, out_right in mapping:
                effects[outcomes[key]][left][right] += int(out_left == out_right)
    for i in range(4):
        for j in range(4):
            if effects[0][i][j]+effects[1][i][j] != 2*int(i == j):
                raise ArithmeticError('native pointer effects are not complete')
            if i != j and (effects[0][i][j] or effects[1][i][j]):
                raise ArithmeticError('native pointer effects are not diagonal')
    return tuple(tuple(effects[bit][code][code] for code in range(4)) for bit in (0, 1))


def native_effects(law):
    """Integer numerators over two, derived from committed source records."""
    return _effect_counts(_law(law))


def exact_instrument_audit():
    """Check every logical matrix unit, not only the two outcome probabilities."""
    ideal = target()
    if not np.array_equal(ideal.T@ideal, np.eye(8, dtype=int)):
        raise ArithmeticError('the complete word target is not a permutation')
    checked = 0
    for law in ('C', 'E'):
        effects = native_effects(law)
        for bit in (0, 1):
            for a in (0, 1):
                for b in (0, 1):
                    for c in (0, 1):
                        for d in (0, 1):
                            numerator = sum(effects[bit][z]*int(ideal[z+4*a, READY_CODE+4*c])
                                            *int(ideal[z+4*b, READY_CODE+4*d]) for z in range(4))
                            wanted = 2*int(a == b == c == d == bit)
                            if numerator != wanted:
                                raise ArithmeticError('native pointer is not the complete Lueders instrument')
                            checked += 1
    if tuple(ready_state()[:2]) != (1, 0):
        raise ArithmeticError('the native reset target is not the ready orientation')
    return {'ready_native_configuration': list(ready_state()), 'ready_source_code': READY_CODE,
            'auxiliary_source_codes_checked': list(range(4)),
            'logical_superoperator_entries_checked': checked,
            'native_effect_numerators_over_two': {law: [list(row) for row in native_effects(law)]
                                                for law in ('C', 'E')},
            'ideal_reduced_instrument': 'rho -> (P0 rho P0, P1 rho P1)',
            'arbitrary_work_padding_and_untouched_references': True,
            'native_record_readout_not_a_supplied_projective_work_read': True}


def certificate():
    """Fresh original-binary source comparison with portable fixed envelopes."""
    audit = exact_instrument_audit()
    actual = compiled.exact._array(local_matrix(), (8, 8))
    ideal = target()
    squared = F(0)
    for row in range(8):
        for column in range(8):
            real, imag = compiled.exact._component(actual[row, column])
            squared += (real-int(ideal[row, column]))**2+imag**2
    if squared > WORD_ERROR_GATE**2:
        raise ArithmeticError('raw native-pointer source word exceeds its operator envelope')
    delta = WORD_ERROR_GATE*(2+WORD_ERROR_GATE)/2
    return {'operator_error_upper': WORD_ERROR_GATE, 'half_diamond_upper': delta,
            'channel_norm_upper': (1+WORD_ERROR_GATE)**2,
            'source_pulses': len(_local_word()), 'proper_code_service_events': 2*len(_local_word()),
            'native_pointer_attempts': 1, 'exact_instrument': audit,
            'raw_word_is_not_normalized_or_projected': True,
            'logical_auxiliary_reset_is_separate_and_must_retain_timeout': True,
            'proper_code_load_export_and_coherent_routing_remain_supplied': True,
            'physical_service_error_authenticated': False}


def comparison_bound(auxiliary_timeout):
    """Complete original-work bound, including the retained reset timeout."""
    if type(auxiliary_timeout) not in (int, F) or not 0 <= auxiliary_timeout <= 1:
        raise ValueError('an exact uniform auxiliary reset timeout in [0,1] required')
    return F(auxiliary_timeout)+certificate()['half_diamond_upper']


def _superoperators(law):
    effects = native_effects(law)
    actual = local_matrix()
    result = np.zeros((2, 2, 2, 2, 2), complex)
    for bit in (0, 1):
        for z in range(4):
            block = actual[np.ix_([z, z+4], [READY_CODE, READY_CODE+4])]
            result[bit] += (effects[bit][z]/2)*np.einsum('ac,bd->abcd', block, block.conj())
    return result


def read_bit(rho, law, bit=7, work_qubits=9, *, reference_dimension=1):
    """Raw source/native reduced output blocks, on the post-reset work state.

    Returns (outcome_zero, outcome_one); both are subnormalized and retain
    every work/reference coordinate.  The separately retained auxiliary
    reset-timeout branch is NOT consumed or synthesized here.  This API is
    invoked only after the canonical reset has conclusively produced the
    pure ready configuration.  Its private classical reset tape is forgotten.

    Numerical admission screens the aggregate negative spectral mass,
    Hermitian defect and trace excess; it is not an outward PSD certificate.
    The channel theorem assumes a genuinely positive subnormalized input.
    Entries within the small roundoff gate are used unchanged, so numerical
    output diagnostics can retain roundoff rather than exact positivity.
    """
    law, (bit, work_qubits) = _law(law), _address(bit, work_qubits)
    if type(reference_dimension) is not int or not 1 <= reference_dimension <= 3:
        raise ValueError('one to three numerical reference levels required')
    dimension = (1 << work_qubits)*reference_dimension
    if (type(rho) is not np.ndarray or rho.dtype.kind not in 'fc'
            or rho.shape != (dimension, dimension) or not np.isfinite(rho).all()):
        raise ValueError('a plain finite work/reference density on the declared register required')
    with np.errstate(over='ignore', invalid='ignore'):
        mass = np.trace(rho)
        hermitian_error = np.sqrt(dimension)*np.linalg.norm(rho-rho.conj().T, 'fro')/2
    if (not np.isfinite(mass) or not np.isfinite(hermitian_error)
            or hermitian_error > float(DENSITY_ROUNDOFF_GATE)):
        raise ValueError('a Hermitian subnormalized work/reference density required')
    # Symmetrize only the validation matrix, never the propagated input.
    eigenvalues = np.linalg.eigvalsh(rho/2+rho.conj().T/2)
    with np.errstate(over='ignore', invalid='ignore'):
        defect = (hermitian_error+max(0., -mass.real, mass.real-1)
                  -np.minimum(eigenvalues, 0).sum())
    if not np.isfinite(eigenvalues).all() or not np.isfinite(defect) or defect > float(DENSITY_ROUNDOFF_GATE):
        raise ValueError('a positive subnormalized density within aggregate roundoff required')
    # No spectral projection or density renormalization enters the execution.
    certificate()
    matrices = _superoperators(law)
    bases = np.array([i for i in range(1 << work_qubits) if not (i >> bit) & 1])
    indices = tuple(((bases+(value << bit))[:, None]*reference_dimension
                     +np.arange(reference_dimension)).reshape(-1) for value in (0, 1))
    outputs = []
    for coefficients in matrices:
        output = np.zeros_like(rho, dtype=complex)
        for a in (0, 1):
            for b in (0, 1):
                for c in (0, 1):
                    for d in (0, 1):
                        output[np.ix_(indices[a], indices[b])] += coefficients[a, b, c, d]*rho[np.ix_(indices[c], indices[d])]
        if not np.isfinite(output).all():
            raise ArithmeticError('raw native-pointer output is nonfinite')
        outputs.append(output)
    return tuple(outputs)


def _json_exact(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    return value


def run():
    report = certificate()
    vector = np.zeros(1024, complex)
    for index, amplitude in ((2*16, .5), (2*144+1, .5j), (2*511, -.5), (2*383+1, -.5j)):
        vector[index] = amplitude
    rho = np.outer(vector, vector.conj())
    diagonal = np.repeat([(i >> 7) & 1 for i in range(512)], 2)
    witnesses = []
    for law in ('C', 'E'):
        outputs = read_bit(rho, law, reference_dimension=2)
        wanted = tuple(rho*(diagonal[:, None] == bit)*(diagonal[None, :] == bit) for bit in (0, 1))
        error = max(float(np.max(abs(a-b))) for a, b in zip(outputs, wanted, strict=True))
        mass_error = abs(sum(np.trace(value) for value in outputs)-1)
        if error > 2e-12 or mass_error > 2e-12:
            raise ArithmeticError('complete native-pointer reference witness failed')
        witnesses.append({'law': law, 'work_dimension': 512, 'reference_dimension': 2,
            'native_and_padding_input_addresses': [16, 144, 511, 383],
            'complete_output_max_entry_error_diagnostic': error,
            'mass_error_diagnostic': float(mass_error),
            'probabilities_diagnostic': [float(np.trace(value).real) for value in outputs]})
    return {'schema': 'oph.native_current_pointer.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'certificate': _json_exact(report), 'source_word_for_work_bit7': source_word(),
        'full_register_reference_executions': witnesses,
        'scope': {'auxiliary_ready_state_must_come_from_complete_native_reset': True,
            'conditional_herald_or_reset_can_steer_initial_work': True,
            'complete_original_work_bound': 'reset_timeout + half_diamond_word_error',
            'work_projective_readout_and_logical_pointer_blank_not_separately_supplied': True,
            'proper_code_service_blanks_not_derived': True,
            'numerical_density_admission_gate': str(DENSITY_ROUNDOFF_GATE),
            'numerical_density_screen_is_outward_psd_certificate': False,
            'new_physical_readout_or_clock_accuracy_claimed': False}}


def compare_receipt(retained, actual):
    """Portable replay never relaxes a diagnostic's domain or execution gate."""
    _portable_compare(retained, actual)
    for report in (retained, actual):
        for row in report['full_register_reference_executions']:
            if any(not 0 <= row[key] <= 2e-12 for key in (
                    'complete_output_max_entry_error_diagnostic', 'mass_error_diagnostic')):
                raise ValueError('native pointer diagnostic violates its execution domain')
            if (any(value < 0 for value in row['probabilities_diagnostic'])
                    or abs(sum(row['probabilities_diagnostic'])-1) > 2e-12):
                raise ValueError('native pointer probabilities violate their execution domain')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', nargs='?', const=str(OUTPUT), type=Path)
    mode.add_argument('--write', nargs='?', const=str(OUTPUT), type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native current pointer receipt verified.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
