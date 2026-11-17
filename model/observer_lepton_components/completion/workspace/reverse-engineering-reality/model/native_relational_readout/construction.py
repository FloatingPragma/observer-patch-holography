"""Join the native current counter to a coherent measurement at its first tick.

The first-current channel is an exactly summed source instrument, not a
fractional channel or a supplied-duration evolution. The bounded runtime
retains its complete unresolved no-current tail and deterministic halt flags.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from ..observer import native
from ..native_quantum_history import construction as quantum
from ..native_clock import clock
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/native_relational_readout/__init__.py',
    'model/native_relational_readout/__main__.py',
    'model/native_relational_readout/construction.py',
    *quantum.SOURCE_PATHS,
    *clock.SOURCE_PATHS,
)))
TOLERANCE = 2e-11
HORIZON = 6


def bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= TOLERANCE:
        raise ArithmeticError(f'{name} exceeds the relational readout tolerance')
    return value


def horizon(value, maximum=256):
    if type(value) is not int or not 1 <= value <= maximum:
        raise ValueError(f'an integer attempt bound from one to {maximum} required')
    return value


@dataclass(frozen=True)
class MarkedSource:
    law: str
    no_current_counts: tuple
    moving_maps: tuple
    instrument: runtime.Instrument


def source(law):
    """Enumerate actual committed labels and classify their actual current."""
    law = quantum.candidate(law)
    groups = {}
    bits = {}
    for seam_number, seam in enumerate(quantum.mixing.seams()):
        for coin in (0, 1):
            for i, state in enumerate(quantum.STATES):
                output, record = native.native(state, seam, coin)
                if record != (coin, state[seam[0]]-output[seam[0]]) or output not in quantum.INDEX:
                    raise ArithmeticError('native current grammar or binary support changed')
                label = record[1] if law == 'C' else state[seam[0]]-state[seam[1]]
                key = (seam_number, coin, label)
                bit = abs(record[1])
                if bit not in (0, 1) or key in bits and bits[key] != bit:
                    raise ArithmeticError('candidate record does not determine a current bit')
                bits[key] = bit
                groups.setdefault(key, {})[i] = quantum.INDEX[output]
    no_counts = [[0]*66 for _ in range(66)]
    moving = []
    matrices = [[], []]
    for key in sorted(groups):
        mapping, bit = groups[key], bits[key]
        if len(set(mapping.values())) != len(mapping):
            raise ArithmeticError('a native recorded Kraus map is not a partial injection')
        matrix = np.zeros((66, 66), complex)
        for i, output in mapping.items():
            matrix[output, i] = 1/math.sqrt(60)
        matrices[bit].append(matrix)
        if bit:
            moving.append(tuple(sorted(mapping.items())))
        else:
            if any(output != i for i, output in mapping.items()):
                raise ArithmeticError('a zero-current branch moves a native configuration')
            for i in mapping:
                for j in mapping:
                    no_counts[i][j] += 1
    if any(no_counts[i][j] > 52 or no_counts[i][j] < 0 for i in range(66) for j in range(66)):
        raise ArithmeticError('source zero-current coherence lacks the uniform escape bound')
    if set(no_counts[i][i] for i in range(66)) != {50, 52}:
        raise ArithmeticError('the native minimum escape probability changed')
    instrument = runtime.Instrument(kraus=tuple(tuple(group) for group in matrices), visible=True,
        resources=(('native_attempts', 1), ('native_current_readouts', 1)),
        binding={'candidate': law, 'scheduler': 'uniform thirty seams and fair placement coin',
                 'outcome': 'absolute current bit',
                 'forgotten_labels': ['seam', 'coin', 'committed candidate repair label'],
                 'native_channel_semantics': True, 'native_pulse_words_emitted': False})
    return MarkedSource(law, tuple(map(tuple, no_counts)), tuple(moving), instrument)


@lru_cache(maxsize=1)
def phase_pairs():
    return quantum.phase_apparatus()[3]


def exact_preparation(sign=1):
    """Exact rational version of the existing declared phase-pair service."""
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('phase sign must be minus or plus one')
    rho = {}
    for i, j in phase_pairs():
        rho[i, i] = rho[j, j] = F(1, 20)
        rho[i, j] = rho[j, i] = F(sign, 20)
    return rho


def exact_stopped_state(marked, *, sign=1, first_at=None, through=None):
    """Sum moving CP branches after zero-current Schur evolution, exactly.

    Neither bound means the infinite first-current instrument. ``first_at``
    retains exactly that first-event count; ``through`` sums events up to it.
    """
    if first_at is not None and through is not None:
        raise ValueError('choose an exact first count or a cumulative count bound')
    if first_at is not None:
        first_at = horizon(first_at)
    if through is not None:
        through = horizon(through)
    output, initial = {}, exact_preparation(sign)
    for mapping in marked.moving_maps:
        mapping = dict(mapping)
        for (i, j), value in initial.items():
            if i not in mapping or j not in mapping:
                continue
            count = marked.no_current_counts[i][j]
            if first_at is not None:
                weight = F(count, 60)**(first_at-1)/60
            else:
                weight = F(1, 60-count)
                if through is not None:
                    weight *= 1-F(count, 60)**through
            key = (mapping[i], mapping[j])
            output[key] = output.get(key, F(0))+value*weight
    return output


def exact_tail(marked, count, sign=1):
    count = horizon(count)
    return {(i, j): value*F(marked.no_current_counts[i][j], 60)**count
            for (i, j), value in exact_preparation(sign).items()}


def dense(state):
    result = np.zeros((66, 66), complex)
    for (i, j), value in state.items():
        result[i, j] = float(value)
    return result


def trace(state):
    return sum((value for (i, j), value in state.items() if i == j), F(0))


def phase_probabilities(state):
    coherence = sum((state.get((i, j), F(0))+state.get((j, i), F(0)) for i, j in phase_pairs()), F(0))
    return ((trace(state)+coherence)/2, (trace(state)-coherence)/2)


def program(marked, count=HORIZON, sign=1):
    count = horizon(count, 6)
    rho, _, read, _ = quantum.phase_apparatus(sign)
    latch = runtime.Instrument(kraus=((np.zeros((66, 66), complex),), (np.eye(66),)), visible=True,
        resources=(('latched_flag_copies', 1),),
        binding={'deterministic_outcome': 1, 'native_updates': 0,
                 'declared_logical_halt_controller': True})
    events = [runtime.Event('first_native_current', instruments=(marked.instrument,))]
    for n in range(1, count):
        events.append(runtime.Event(f'continue_or_keep_halted_{n+1}',
            instruments=(marked.instrument, latch), control=n-1))
    events.append(runtime.Event('complete_phase_readout', instruments=(read,)))
    return runtime.Program(f'first_current_phase_{marked.law}_{count}_{sign}', 66, rho, tuple(events),
        binding={'candidate': marked.law, 'phase_sign': sign,
                 'rule': 'stop native attempts at first moving current; retain unresolved tail',
                 'attempt_budget': count,
                 'flag_semantics': 'zero until first current; then repeated deterministic one',
                 'physical_duration_supplied': False, 'physical_halt_device_claimed': False})


def runtime_check(marked, count=HORIZON, sign=1):
    p = program(marked, count, sign)
    history = runtime.execute(p)
    read = p.events[-1].instruments[0]
    error, probability_error, impossible_mass = 0., 0., 0.
    rows = []
    finite, tail = exact_stopped_state(marked, sign=sign, through=count), exact_tail(marked, count, sign)
    if trace(finite)+trace(tail) != 1:
        raise ArithmeticError('the exact first-current instrument discarded its unresolved tail')
    moving_probabilities = [0., 0.]
    no_current_probabilities = [0., 0.]
    stopped_states = {n: exact_stopped_state(marked, sign=sign, first_at=n) for n in range(1, count+1)}
    stopped_states[None] = tail
    for branch in history.branches:
        flags, phase = branch.record[:-1], branch.record[-1]
        possible = all(a <= b for a, b in zip(flags, flags[1:]))
        first = flags.index(1)+1 if 1 in flags else None
        state = stopped_states[first] if possible else {}
        expected_probability = phase_probabilities(state)[phase]
        k = read.kraus[phase][0]
        wanted = k@dense(state)@k.conj().T
        error = max(error, float(np.linalg.norm(branch.density-wanted, 'fro')))
        probability_error = max(probability_error, abs(branch.probability-float(expected_probability)))
        if not possible:
            impossible_mass += branch.probability
        else:
            expected_cost = count if first is None else first
            resources = dict(branch.resources)
            if (resources.get('native_attempts', 0) != expected_cost
                    or resources.get('native_current_readouts', 0) != expected_cost
                    or resources.get('latched_flag_copies', 0) != count-expected_cost
                    or resources.get('phase_readouts', 0) != 1):
                raise ArithmeticError('halt feedback changed the retained native resource cost')
            target = no_current_probabilities if first is None else moving_probabilities
            target[phase] += branch.probability
        rows.append({'flags': list(flags), 'phase': phase,
                     'exact_probability': str(expected_probability), 'runtime_probability': branch.probability})
    prefix = runtime.execute(p, stop_after=2)
    restored = runtime.checkpoint_from_json(p, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(p, checkpoint=restored)
    return {'candidate': marked.law, 'phase_sign': sign, 'attempt_budget': count,
            'complete_flag_and_phase_records': rows,
            'full_terminal_density_dimension': 66,
            'complete_terminal_density_residual_diagnostic': bounded(error, 'stopped state'),
            'complete_probability_residual_diagnostic': bounded(probability_error, 'stopped phase law'),
            'impossible_flag_words_mass_diagnostic': bounded(impossible_mass, 'irreversible flag history'),
            'first_current_and_phase_probabilities': list(map(str, phase_probabilities(finite))),
            'no_current_and_phase_probabilities': list(map(str, phase_probabilities(tail))),
            'runtime_first_current_and_phase_probabilities': moving_probabilities,
            'runtime_no_current_and_phase_probabilities': no_current_probabilities,
            'unresolved_tail_mass': str(trace(tail)),
            'uniform_all_input_tail_bound': str(F(13, 15)**count),
            'actual_branch_native_costs_checked': True,
            'all_zero_probability_and_unresolved_tail_records_retained': True,
            'checkpoint_complete_cq_distance_diagnostic': bounded(runtime.complete_distance(continued, history), 'latched checkpoint'),
            'expanded_tree_resources': dict(history.expanded_tree_resources)}


def exact_audit(marked):
    state = exact_stopped_state(marked)
    if trace(state) != 1:
        raise ArithmeticError('the infinite first-current instrument is not normalized')
    degrees = []
    kappa = F(1, 6) if marked.law == 'C' else F(3, 10)
    for i, j in phase_pairs():
        degree = 60*(1-kappa)-marked.no_current_counts[i][j]
        if degree.denominator != 1:
            raise ArithmeticError('the phase-sector native moving degree is not integral')
        degrees.append(int(degree))
    histogram = {str(d): degrees.count(d) for d in sorted(set(degrees))}
    if histogram != {'3': 2, '4': 4, '5': 4}:
        raise ArithmeticError('the native phase-sector degree multiplicities changed')
    predicted_coherence = sum((F(degree, 60*kappa+degree) for degree in degrees), F(0))/10
    probabilities = phase_probabilities(state)
    if probabilities != ((1+predicted_coherence)/2, (1-predicted_coherence)/2):
        raise ArithmeticError('native first-event sum differs from the sector formula')
    # The infinite map is CP as a convergent sum of T1 T0^n. Its diagonal
    # action is stochastic: sum outgoing moving multiplicities/(60-W_ii)=1.
    for i in range(66):
        moving_count = sum(i in dict(mapping) for mapping in marked.moving_maps)
        if F(moving_count, 60-marked.no_current_counts[i][i]) != 1:
            raise ArithmeticError('the eventual-current instrument loses a basis-state trace')
    # Global inserted idling cancels in every resolvent multiplier. This is
    # an exact all-operator-input identity, not a trajectory resampling claim.
    q = F(1, 2)
    for row in marked.no_current_counts:
        for count in row:
            original = F(1, 60-count)
            lazy = q/F(60)/(1-((1-q)+q*F(count, 60)))
            if lazy != original:
                raise ArithmeticError('global idle rescaling changed the stopped instrument')
    return {'candidate': marked.law, 'complete_phase_probabilities': list(map(str, probabilities)),
            'remaining_phase_coherence': str(predicted_coherence),
            'phase_sector_moving_degree_histogram': histogram,
            'zero_current_Kraus_maps_are_diagonal': True,
            'zero_current_Schur_denominator': 60,
            'zero_current_Schur_counts_sha256': hashlib.sha256(canonical(marked.no_current_counts).encode()).hexdigest(),
            'resolvent': 'R_xy=60/(60-Wcount_xy); stopped map=T1(R entrywise times rho)',
            'all_66_basis_input_trace_identities_checked': True,
            'uniform_tail_bound_every_density_input': '(13/15)^N',
            'positive_global_idle_rescaling_preserves_complete_stopped_quantum_instrument': True,
            'global_idle_control_probability': '1/2',
            'global_idle_invariance_preserves_attempt_cost_or_zero_records': False}


def clock_attachment(sources):
    """The coherent preparation realizes the earlier clock's Palm class law."""
    kernel = clock.transfer.native_counts()
    if kernel.states != quantum.STATES:
        raise ArithmeticError('the current clock and phase observer use different native bases')
    reduction = clock.lumping(kernel)
    prepared = exact_preparation()
    classes = tuple(sum((value for (i, j), value in prepared.items()
                        if i == j and reduction['classes'][i] == c), F(0)) for c in range(3))
    if classes != reduction['palm']:
        raise ArithmeticError('the phase preparation does not have the derived clock class law')
    waits, tail = clock.waiting_law(reduction, 16, start='palm')
    rows = []
    for law, marked in sources.items():
        joint = []
        kappa = F(1, 6) if law == 'C' else F(3, 10)
        for n in range(1, 17):
            state = exact_stopped_state(marked, first_at=n)
            probabilities = phase_probabilities(state)
            coherence = sum((F(multiplicity*degree, 600)*(1-kappa-F(degree, 60))**(n-1)
                             for degree, multiplicity in ((3, 2), (4, 4), (5, 4))), F(0))
            expected = ((waits[n-1]+coherence)/2, (waits[n-1]-coherence)/2)
            if probabilities != expected:
                raise ArithmeticError('the joint current-clock/phase law differs from its source formula')
            joint.append({'first_current_attempt': n, 'joint_phase_probabilities': list(map(str, probabilities))})
        if trace(exact_tail(marked, 16)) != tail:
            raise ArithmeticError('the coherent experiment lost the native clock no-event tail')
        rows.append({'candidate': law, 'joint_first_current_phase_law': joint,
                     'no_current_after_sixteen_attempts': str(tail)})
    return {'prepared_class_law': list(map(str, classes)),
            'equals_existing_native_clock_Palm_class_law': True,
            'same_first_current_waiting_law_for_C_and_E': True,
            'mean_first_current_attempts': '33/5',
            'joint_law': rows,
            'fresh_preparation_was_not_conditioned_on_an_unexecuted_event': True}


def run():
    sources = {law: source(law) for law in ('C', 'E')}
    if sources['C'].moving_maps != sources['E'].moving_maps:
        raise ArithmeticError('the native moving-current CP maps differ between candidates')
    exact = {law: exact_audit(marked) for law, marked in sources.items()}
    delta = F(exact['C']['complete_phase_probabilities'][0])-F(exact['E']['complete_phase_probabilities'][0])
    if delta != F(3644, 69069):
        raise ArithmeticError('the first-native-current phase discriminator changed')
    return {'schema': 'oph.native_relational_readout.v1', 'all_pass': True,
            'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'exact_first_current_measurement': exact,
            'same_moving_current_CP_map_for_both_candidates': True,
            'C_minus_E_plus_probability': str(delta),
            'native_clock_attachment': clock_attachment(sources),
            'finite_executions': [runtime_check(marked, sign=sign) for marked in sources.values() for sign in (1, -1)],
            'scope': {'native_current_and_coherent_readout_share_one_executed_history': True,
                      'phase_preparation_and_effects_supplied': True,
                      'logical_halt_and_record_latch_supplied': True,
                      'physical_halt_device_or_clock_calibration_derived': False,
                      'uniform_scheduler_and_native_candidate_supplied': True,
                      'independent_clock_or_Hamiltonian_supplied': False,
                      'infinite_first_current_law_is_exactly_summed_not_finitely_executed': True,
                      'finite_horizon_failures_or_no_current_tail_discarded': False,
                      'native_pulse_words_emitted': False,
                      'laboratory_measurement_identified_or_empirical_agreement_claimed': False,
                      'new_Lean_proof_or_outward_roundoff_bound': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if result.get('all_pass') is not True:
        raise ArithmeticError('native relational readout failed')
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native relational readout receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')
