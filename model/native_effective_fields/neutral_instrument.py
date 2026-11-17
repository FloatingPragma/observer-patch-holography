"""Source-coded preparation and complete one-bit-pointer neutral readout.

Two actual fixed-seam currents herald a pure orientation at common port 2.
The existing mixer creates its phase; a bounded first-global-movement run
retains timeout. Nine source CNOT/copy/read/reset operations analyze every
encoded address. The resulting POVM agrees with the spatial phase analyzer;
its complement instrument additionally destroys address coherence.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import numpy as np

from . import neutral_transport as neutral
from ..compiled_native_interface import construction as compiled
from ..common_source import exact_output as exact
from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..native_relational_readout import construction as stopped
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..primitive_observer import runtime
from ..record_laws import burn_in
from ..predictions import repair_response
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PAIR, COMMON, FIRST_SEAM = (0, 1), 2, (2, 3)
SOURCE_PATHS = tuple(sorted(set((*neutral.SOURCE_PATHS, *compiled.SOURCE_PATHS, *repair_response.SOURCE_PATHS,
    'model/record_laws/burn_in.py', 'model/predictions/repair_response.py',
    'model/native_effective_fields/neutral_instrument.py'))))


def herald(law):
    """Coarse-grain TWO actual repair records; all other words are failure.

    First current -1 on (2,3) places a token at 2. Then either nonzero
    current on (0,1) fixes the other token's orientation. Total charge two
    makes the successful output one-dimensional, even for coherent inputs.
    """
    groups = [[], [], []]
    for first in interface.fixed_maps(law, FIRST_SEAM):
        a = dict(first.mapping)
        for second in interface.fixed_maps(law, PAIR):
            b = dict(second.mapping)
            mapping = tuple((i, b[j]) for i, j in a.items() if j in b)
            flag = (1 if second.current == -1 else 2) if first.current == -1 and second.current else 0
            k = np.zeros((66, 66), complex)
            for i, j in mapping:
                k[j, i] = .5
            groups[flag].append((mapping, k))
    return tuple(tuple(group) for group in groups)


def herald_audit(law):
    groups = herald(law)
    counts, successful = [F(0)]*3, {}
    for flag, group in enumerate(groups):
        outputs = set()
        for mapping, _ in group:
            outputs.update(j for _, j in mapping)
            counts[flag] += F(len(mapping), 264)
        if flag:
            wanted = neutral.addresses(PAIR, COMMON)[flag-1]
            if outputs != {wanted} or counts[flag] != F(1, 264):
                raise ArithmeticError('two native records do not herald the stated pure orientation')
            successful[flag] = wanted
    if sum(counts) != 1:
        raise ArithmeticError('failed native preparation records were removed')
    return tuple(counts), successful


def copy_word(bit):
    if type(bit) is not int or not 0 <= bit < 9:
        raise ValueError('one of the nine encoded work bits required')
    return compiler.cnot(bit, 9)


@lru_cache(maxsize=1)
def copy_matrix():
    """Execute original source arrays on work bit 0 and pointer bit 1."""
    a = compiler.execute(compiler.cnot(0, 1), 2)
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(4, 4)


def copy_certificate(matrix=None):
    """Exact rational enclosure of the RAW binary CNOT, all four columns."""
    matrix = copy_matrix() if matrix is None else matrix
    if type(matrix) is not np.ndarray:
        raise ValueError('an original base ndarray required')
    a = exact._array(matrix, (4, 4))
    error2 = F(0)
    for i in range(4):
        for j in range(4):
            real, imag = exact._component(a[i, j])
            real -= int(i == (j ^ (2 if j & 1 else 0)))
            error2 += real*real+imag*imag
    eps = exact._sqrt_interval(error2)[1]
    return {'operator_error_upper': eps, 'half_diamond_upper': eps*(2+eps)/2,
            'channel_norm_upper': (1+eps)**2}


def pointer_kraus(matrix=None):
    """Blank pointer, source copy, binary pointer read, then pointer reset.

    Reset supplies zero before each copy: one initial blank and eight
    subsequent resets. The final pointer need not be reset. These operations
    retain all nine classical readout bits. No work projection
    or normalization is performed on the raw copy arrays.
    """
    a = copy_matrix() if matrix is None else exact._array(matrix, (4, 4))
    return tuple(a[np.ix_([2*bit, 2*bit+1], [0, 1])] for bit in (0, 1))


def coarse_address(address):
    if type(address) is not int or not 0 <= address < 512:
        raise ValueError('one complete nine-bit readout address required')
    for n, common in enumerate(neutral.common_ports(PAIR)):
        for sign, index in enumerate(neutral.addresses(PAIR, common)):
            if address == circuit.encode(quantum.STATES[index]):
                return 2*n+sign
    return 20


def ideal_readout():
    """Twenty singleton effects and a 46-map dephasing complement.

    This is the native-code restriction of the complete 512-address
    source pointer instrument. It has the old POVM, not its Lueders update.
    """
    b = compiled.ideal_native_matrix()
    groups = [[] for _ in range(21)]
    for i, state in enumerate(quantum.STATES):
        k = np.zeros((66, 66), complex); k[i, :] = b[i, :]
        groups[coarse_address(circuit.encode(state))].append(k)
    return runtime.Instrument(kraus=tuple(tuple(group) for group in groups), visible=True,
        resources=(('source_native_pulses', 61), ('proper_code_service_events', 122),
                   ('pointer_blank_resets', 9), ('pointer_bit_reads', 9)),
        binding={'raw_record_alphabet': 512, 'coarse_outcomes': 21,
                 'raw_bits_retained_by_pointer_execution': True,
                 'ideal_native_surrogate_not_raw_array_projection': True,
                 'complement_instrument_is_address_dephasing': True})


def raw_readout_probabilities(native_density):
    """Execute full source mixer and nine actual copy/read instruments.

    Tensor contraction keeps all original binary entries and all 512 raw
    records, including non-binary carrier and padding addresses. This is
    floating execution; its arithmetic residual is a diagnostic only.
    """
    rho = exact._array(native_density, (66, 66))
    embedding = compiled.embedding()
    encoded = embedding@rho@embedding.T
    b = compiled.local_matrix()
    raw = np.einsum('ab,ibjd,cd->iajc', b, encoded.reshape(64, 8, 64, 8), b.conj(), optimize=True)
    tensor = raw.reshape((2,)*18)
    labels = list(range(18))
    kraus = pointer_kraus()
    effects = np.array([k.conj().T@k for k in kraus])
    for axis in range(9):
        out = [x for x in labels if x not in (axis, axis+9)]+[18+axis]
        tensor = np.einsum(tensor, labels, effects, [18+axis, axis+9, axis], out, optimize=True)
        labels = out
    probabilities = tensor.transpose([labels.index(18+i) for i in range(9)]).reshape(512)
    if np.max(np.abs(probabilities.imag)) > 1e-11 or min(probabilities.real) < -1e-11:
        raise ArithmeticError('raw source pointer produced invalid diagnostic masses')
    # No clipping, normalization or removal of a raw outcome.
    return probabilities.real


def analyzer_audit():
    read, old = ideal_readout(), neutral.readout(PAIR)
    effect_error = 0.
    for group, old_group in zip(read.kraus, old.kraus, strict=True):
        effect = sum((k.conj().T@k for k in group), np.zeros((66, 66), complex))
        target = sum((k.conj().T@k for k in old_group), np.zeros((66, 66), complex))
        effect_error = max(effect_error, float(np.linalg.norm(effect-target, 'fro')))
    if effect_error > 1e-12:
        raise ArithmeticError('pointer copies fail to implement the spatial phase POVM')
    active = {i for common in neutral.common_ports(PAIR) for i in neutral.addresses(PAIR, common)}
    i, j = [i for i in range(66) if i not in active][:2]
    v = np.zeros(66); v[i] = v[j] = 1/np.sqrt(2)
    rho = np.outer(v, v)
    actual = sum((k@rho@k.conj().T for k in read.kraus[20]), np.zeros_like(rho, dtype=complex))
    difference = float(np.linalg.norm(actual-rho, 'nuc')/2)
    if abs(difference-.5) > 1e-12:
        raise ArithmeticError('complement disturbance countercontrol changed')
    raw = raw_readout_probabilities(stopped.dense(neutral.preparation(PAIR, COMMON)))
    predicted = np.zeros(512); predicted[circuit.encode(quantum.STATES[neutral.addresses(PAIR, COMMON)[0]])] = 1
    return {'same_complete_21_outcome_POVM_residual_diagnostic': effect_error,
            'complement_update_trace_distance_countercontrol': difference,
            'raw_pointer_readout_count': len(raw), 'raw_mass_sum_diagnostic': float(sum(raw)),
            'raw_pointer_prepared_plus_L1_error_diagnostic': float(sum(abs(raw-predicted))),
            'raw_padding_mass_diagnostic': float(sum(p for code, p in enumerate(raw) if code not in set(circuit.SOURCE_CODES))),
            'full_Lueders_instrument_claimed': False}


def source_budget(warmup=369):
    mix, copy = compiled.certify(), copy_certificate()
    # Chronological preparation mixer, analyzer mixer, then nine copies.
    error, gain = F(0), F(1)
    for step in [mix, mix]+[copy]*9:
        error = step['channel_norm_upper']*error+step['half_diamond_upper']
        gain *= step['channel_norm_upper']
    delta = burn_in.binary_bound(warmup)
    return {'warmup_attempts': warmup, 'warmup_half_trace_error_upper': delta,
            'raw_source_half_diamond_error_upper': error,
            'raw_source_channel_norm_upper': gain,
            'complete_output_error_upper': error+gain*delta,
            'two_mixers_source_pulses': 68, 'nine_pointer_copies_source_pulses': 27,
            'maximum_source_pulses_per_trial': 95, 'maximum_proper_code_service_events': 190,
            'pointer_blank_resets_per_success': 9, 'pointer_bit_reads_per_success': 9,
            'native_herald_attempts': 2, 'movement_attempt_cap': 2,
            'maximum_total_native_attempts': warmup+4,
            'raw_source_comparison_padding_extension': 'dephase native-code/complement before native repair; act natively on code and retain complement in an absorbing failure record',
            'physical_control_reset_read_error_included': False}


def herald_instrument(law):
    groups = herald(law)
    return runtime.Instrument(kraus=tuple(tuple(k for _, k in group) for group in groups),
        visible=True, resources=(('native_attempts', 2), ('native_record_readouts', 2)),
        binding={'two_actual_seams': [list(FIRST_SEAM), list(PAIR)],
                 'outcomes': ['failed', 'plus_orientation', 'minus_orientation']})


def readout_factors():
    """Numerically stable factorization of the IDEAL source readout target.

    P_y then B then address dephasing equals sum_i |i><i| B rho B* |i><i|
    for addresses with label y. The factors are a simulation identity, not
    three separately assumed physical operations. Exact half-entry P_y
    avoids subtractive cancellation at mathematically impossible outcomes.
    """
    old = neutral.readout(PAIR)
    phase = runtime.Instrument(kraus=old.kraus, visible=True,
        resources=ideal_readout().resources,
        binding={'factor_of_ideal_source_pointer_target_only': True})
    b = runtime.Instrument(kraus=((compiled.ideal_native_matrix(),),))
    basis = tuple(np.diag(np.eye(66)[i]).astype(complex) for i in range(66))
    diagonal = runtime.Instrument(kraus=(basis,))
    return phase, b, diagonal


def exact_law(law, cap=2):
    cap = stopped.horizon(cap, maximum=2)
    masses, _ = herald_audit(law)
    source = stopped.source(law)
    rows = {(0, 0, 20): masses[0]}
    for flag, sign in ((1, 1), (2, -1)):
        initial = neutral.preparation(PAIR, COMMON, sign)
        for first in range(cap+1):
            state = (neutral.exact_first_current(source, initial, first=first) if first
                     else neutral.exact_first_current(source, initial, tail=cap))
            for outcome, p in enumerate(neutral.exact_probabilities(state, PAIR)):
                rows[flag, first, outcome] = masses[flag]*p
    if sum(rows.values()) != 1 or min(rows.values()) < 0:
        raise ArithmeticError('complete exact protocol lost failure or timeout')
    return rows


def execution(law):
    """Execute finite continuations with native herald masses retained.

    Each sign continuation has a certain native orientation input proved
    by the actual two-current herald. Conditional output is multiplied by
    its separately executed unnormalized herald mass in the full law.
    """
    prep, marked = herald_instrument(law), stopped.source(law)
    herald_program = runtime.Program('native_two_current_herald_'+law, 66, np.eye(66)/66,
        (runtime.Event('two_current_records', instruments=(prep,)),))
    prepared = runtime.execute(herald_program)
    expected, rows = exact_law(law), {(0, 0, 20): prepared.branches[0].probability}
    residual, checkpoint_error, raw_error = 0., 0., 0.
    read, analysis_basis, address_dephasing = readout_factors()
    for flag, sign in ((1, 1), (2, -1)):
        branch = next(b for b in prepared.branches if b.record == (flag,))
        i = neutral.addresses(PAIR, COMMON)[flag-1]
        pure = np.zeros((66, 66), complex); pure[i, i] = 1
        residual = max(residual, float(np.linalg.norm(branch.density-float(F(1, 264))*pure, 'fro')))
        # The initial basis state is exactly the native successful output.
        hold = runtime.Instrument(kraus=((np.zeros((66, 66)),), (np.eye(66),)), visible=True)
        b = runtime.Instrument(kraus=((compiled.ideal_native_matrix(),),), resources=(('source_native_pulses', 34),))
        events = (runtime.Event('source_phase_preparation', instruments=(b,)),
                  runtime.Event('first_global_current', instruments=(marked.instrument,)),
                  runtime.Event('second_current_or_latch', instruments=(marked.instrument, hold), control=0),
                  runtime.Event('source_spatial_pointer', instruments=(read,)),
                  runtime.Event('ideal_pointer_target_factorization', sequence=(analysis_basis, address_dephasing)))
        p = runtime.Program('source_neutral_'+law+'_'+str(sign), 66, pure, events)
        result = runtime.execute(p)
        prefix = runtime.execute(p, stop_after=2)
        replay = runtime.execute(p, checkpoint=runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix)))
        checkpoint_error = max(checkpoint_error, runtime.complete_distance(result, replay))
        for output in result.branches:
            flags, outcome = output.record[:2], output.record[2]
            possible = flags != (1, 0)
            first = flags.index(1)+1 if 1 in flags else 0
            wanted = expected[flag, first, outcome] if possible else F(0)
            residual = max(residual, abs(branch.probability*output.probability-float(wanted)))
            wanted_density = np.zeros((66, 66), complex)
            if possible:
                state = (neutral.exact_first_current(marked, neutral.preparation(PAIR, COMMON, sign), first=first)
                         if first else neutral.exact_first_current(marked, neutral.preparation(PAIR, COMMON, sign), tail=2))
                if outcome < 20:
                    common = neutral.common_ports(PAIR)[outcome//2]
                    address = neutral.addresses(PAIR, common)[outcome % 2]
                    wanted_density[address, address] = float(wanted)
                else:
                    for address in range(66):
                        if coarse_address(circuit.encode(quantum.STATES[address])) == 20:
                            wanted_density[address, address] = float(F(1, 264)*state.get((address, address), F(0)))
                rows[flag, first, outcome] = branch.probability*output.probability
                expected_native_attempts = first if first else 2
                if dict(output.resources).get('native_attempts') != expected_native_attempts:
                    raise ArithmeticError('a halted branch incurred a native repair')
            residual = max(residual, float(np.linalg.norm(branch.probability*output.density-wanted_density, 'fro')))
        # Raw source address readout, applied to the independently summed
        # completed and unresolved ideal source states; all addresses kept.
        initial = neutral.preparation(PAIR, COMMON, sign)
        for first in (0, 1, 2):
            state = (neutral.exact_first_current(marked, initial, first=first) if first
                     else neutral.exact_first_current(marked, initial, tail=2))
            raw = raw_readout_probabilities(stopped.dense(state))
            coarse = [sum(p for code, p in enumerate(raw) if coarse_address(code) == out) for out in range(21)]
            target = neutral.exact_probabilities(state, PAIR)
            raw_error = max(raw_error, sum(abs(p-float(q)) for p, q in zip(coarse, target, strict=True)))
    if rows.keys() != expected.keys():
        raise ArithmeticError('finite execution lost a declared complete output label')
    if residual > 2e-11 or checkpoint_error > 2e-11 or raw_error > 2e-11:
        raise ArithmeticError('source-coded finite instrument failed its exact-law comparison')
    return {'candidate': law, 'complete_coarse_record_count': len(expected),
            'complete_probability_sum': str(sum(expected.values())),
            'failed_herald_probability': str(expected[0, 0, 20]),
            'herald_probability_per_sign': '1/264',
            'success_and_timeout_probability': str(sum(p for (flag, first, _), p in expected.items() if flag and not first)),
            'exact_complete_law': [[list(key), str(value)] for key, value in sorted(expected.items())],
            'runtime_probability_residual_diagnostic': residual,
            'checkpoint_complete_distance_diagnostic': checkpoint_error,
            'raw_source_analyzer_L1_residual_diagnostic': float(raw_error),
            'surrogate_instrument_uses_ideal_target_not_compressed_raw_arrays': True}


def finite_warmup_audit(warmup=369):
    """Actual finite diagonal preparation; quantum/reference bound is uniform.

    Every success is pure independently of preparation. Warm-up controls
    the success weights and full failed/successful experiment without
    dividing by the small herald probability.
    """
    if tuple(repair_response.source().states) != tuple(quantum.STATES):
        raise ArithmeticError('warm-up and neutral source coordinate ordering differ')
    prior = repair_response.warmup_prior(quantum.STATES[0], warmup)
    mass = []
    for group in herald('C')[1:]:
        mass.append(sum((prior[i]/4 for mapping, _ in group for i, _ in mapping), F(0)))
    prior_error = sum((abs(p-F(1, 66)) for p in prior), F(0))/2
    full_error = (sum((abs(p-F(1, 264)) for p in mass), F(0))
                  + abs(sum(mass)-F(1, 132)))/2
    if not full_error <= prior_error <= burn_in.binary_bound(warmup):
        raise ArithmeticError('actual native warm-up violates complete-law contraction')
    return {'actual_warmup_attempts': warmup,
            'actual_full_flagged_law_error_diagnostic': float(full_error),
            'actual_prior_error_diagnostic': float(prior_error),
            'all_failures_retained_and_no_success_probability_division': True}


def run():
    laws = {law: exact_law(law) for law in ('C', 'E')}
    full = sum((abs(p-laws['E'][key]) for key, p in laws['C'].items()), F(0))/2
    complete = sum((abs(p-laws['E'][key]) for key, p in laws['C'].items() if key[0] and key[1]), F(0))/2
    budget = source_budget()
    if budget['raw_source_half_diamond_error_upper'] > F(1, 10**10):
        raise ArithmeticError('source implementation exceeds its fixed error envelope')
    compact = {k: str(F(-(-(v*10**12).numerator//(v*10**12).denominator), 10**12)) if type(v) is F else v
               for k, v in budget.items()}
    return {'schema': 'oph.native_neutral_instrument.v1', 'all_pass': True,
        'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'analyzer': analyzer_audit(), 'source_budget': compact,
        'finite_native_warmup': finite_warmup_audit(),
        'copy_pointer_source_words': [copy_word(bit) for bit in range(9)],
        'finite_executions': [execution(law) for law in ('C', 'E')],
        'finite_complete_C_E_total_variation': str(full),
        'finite_completed_movement_subdistribution_half_L1': str(complete),
        'ideal_uncapped_pure_preparation_distance_not_finite_trial_score': '4/69',
        'scope': {'phase_preparation_supplied_as_arbitrary_state': False,
            'spatial_analyzer_supplied_as_arbitrary_projectors': False,
            'existing_coherent_proper_code_and_pointer_reset_read_services_required': True,
            'uniform_random_seams_fair_coins_global_current_access_required': True,
            'warmup_tapes_discarded_reference_safe_bound': True,
            'native_or_source_binary_array_errors_are_not_physical_fidelity': True,
            'passive_native_repairs_create_coherence_from_diagonal_inputs': False,
            'particle_pole_or_physical_clock_identification': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path); modes.add_argument('--check', type=Path)
    args = parser.parse_args(); receipt = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), receipt)
        print('Source-coded neutral instrument receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(receipt))
        print(f'Wrote {args.write}')
    else:
        print(canonical(receipt), end='')


if __name__ == '__main__':
    main()
