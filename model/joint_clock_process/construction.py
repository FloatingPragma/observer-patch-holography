"""Full quantum record laws, exact memory, and two source-local event clocks.

All algebraic identities use rational matrices. The actual primitive runtime
uses binary64 Kraus execution; its comparison residual is only a diagnostic.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp

from ..federated_quantum_sector import construction as q
from ..observer.scalar_response.experiment import canonical
from ..primitive_observer import runtime
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*q.SOURCE_PATHS,
    'model/federated_quantum_sector/construction.py',
    'model/federated_observer/MEMORY.md',
    'model/relational_event_clock/kernel.py',
    'model/primitive_observer/runtime.py',
    'model/common_source/open_protocol.py',
    'model/joint_clock_process/__init__.py',
    'model/joint_clock_process/__main__.py',
    'model/joint_clock_process/construction.py',
    'model/joint_clock_process/clock_change.py',
    'Lean/Time/TimeOrderLedger.lean',
    'Lean/Time/ClockReadout.lean',
    'Lean/Screen/RegionalContinuity.lean'))))
MARKS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))


@lru_cache(maxsize=1)
def geometry():
    graph, vertices, tape, _, _ = q.source.tape_data()
    edges = tuple(sorted(set(tape)))
    if vertices != (0, 1, 72, 82) or edges != ((0, 1), (1, 82), (72, 82)):
        raise ArithmeticError('source crossing-path addresses changed')
    if any(edge not in graph.seams for edge in edges):
        raise ArithmeticError('active edge is not a source federation seam')
    return vertices, edges


def _law(law):
    return q.candidate(law)


def _cap(value):
    if type(value) is not int or not 1 <= value <= 8:
        raise ValueError('an exact attempt cap from1 through8 required')
    return value


def _matrix(matrix, size):
    if not isinstance(matrix, sp.MatrixBase) or matrix.shape != (size, size):
        raise ValueError('an exact matrix of the declared size required')
    if any(value.has(sp.Float) or value.is_finite is not True for value in matrix):
        raise ValueError('finite exact quantum matrix entries required')
    return sp.ImmutableMatrix(matrix)


def vector(matrix):
    return sp.Matrix(list(_matrix(matrix, 4)))


def unvector(column):
    if not isinstance(column, sp.MatrixBase) or column.shape != (16, 1):
        raise ValueError('a complete sixteen-entry operator column required')
    return sp.Matrix(4, 4, list(column))


def superoperator(group):
    result = sp.zeros(16)
    if not group:
        raise ValueError('a nonempty exact Kraus group required')
    for weight, matrix in group:
        if type(weight) is not F or weight < 0:
            raise ValueError('nonnegative rational Kraus weight required')
        matrix = sp.Matrix(matrix)
        _matrix(matrix, 4)
        result += sp.Rational(weight.numerator, weight.denominator)*sp.kronecker_product(matrix, matrix.conjugate())
    return sp.ImmutableMatrix(result)


@lru_cache(maxsize=2)
def source_groups(law):
    """Current-derived joint marks; hidden seams and coins are CP alternatives."""
    _law(law)
    vertices, edges = geometry()
    clocks = (vertices[0], vertices[2])
    groups = [[] for _ in MARKS]
    for edge in edges:
        for current, group in zip(q.CURRENT, q.native_groups(law, vertices, edge), strict=True):
            mark = tuple(current*(int(edge[1] == p)-int(edge[0] == p)) for p in clocks)
            for weight, matrix in group:
                groups[MARKS.index(mark)].append((weight/F(len(edges)), matrix))
    return tuple(tuple(group) for group in groups)


@lru_cache(maxsize=2)
def source_maps(law):
    return tuple(superoperator(group) for group in source_groups(_law(law)))


def source_instrument(law):
    return q._runtime_instrument(source_groups(_law(law)), visible=True,
        resources=(('native_attempts', 1), ('joint_clock_reads', 1)),
        binding={'candidate': law, 'source_edges': geometry()[1],
                 'uniform_controlled_three_edge_schedule': True,
                 'clock_ports': [0, 72], 'marks': MARKS,
                 'coin_and_seam_are_hidden_CP_alternatives': True})


@lru_cache(maxsize=1)
def field_maps():
    """R traces the within-carrier address; J adjoins its maximally mixed state."""
    reduce = sp.zeros(4, 16)
    lift = sp.zeros(16, 4)
    for a, b, h in product(range(2), repeat=3):
        reduce[2*a+b, 4*(2*a+h)+2*b+h] = 1
        lift[4*(2*a+h)+2*b+h, 2*a+b] = sp.Rational(1, 2)
    if reduce*lift != sp.eye(4):
        raise ArithmeticError('field restriction and lift do not split')
    return sp.ImmutableMatrix(reduce), sp.ImmutableMatrix(lift)


def reduced(matrix):
    reduce, _ = field_maps()
    return sp.Matrix(2, 2, list(reduce*vector(matrix)))


def memory_forecast(initial, maps):
    """Exact per-record Schur recurrence, without querying future fine states.

    Maps can be outcome maps, interventions or readouts. Their complete menus,
    not this linear routine, must supply CP and trace-preservation hypotheses.
    The hidden forcing is retained; individual memory terms are not CP maps.
    """
    initial = _matrix(initial, 4)
    if type(maps) is not tuple or len(maps) > 8:
        raise ValueError('at most eight exact record-conditioned maps required')
    for matrix in maps:
        _matrix(matrix, 16)
    reduce, lift = field_maps()
    hidden = sp.eye(16)-lift*reduce
    x = [reduce*vector(initial)]
    y0 = hidden*vector(initial)
    blocks = [memory_blocks(sp.ImmutableMatrix(m)) for m in maps]
    for n, (a, b, _, _) in enumerate(blocks):
        forcing = y0
        for previous in range(n):
            forcing = blocks[previous][3]*forcing
        current = a*x[n]+b*forcing
        for j in range(n):
            memory = blocks[j][2]*x[j]
            for middle in range(j+1, n):
                memory = blocks[middle][3]*memory
            current += b*memory
        x.append(current)
    return tuple(sp.ImmutableMatrix(2, 2, list(v)) for v in x)


@lru_cache(maxsize=64)
def memory_blocks(matrix):
    matrix = _matrix(matrix, 16)
    reduce, lift = field_maps()
    hidden = sp.eye(16)-lift*reduce
    return tuple(sp.ImmutableMatrix(m) for m in
                 (reduce*matrix*lift, reduce*matrix*hidden,
                  hidden*matrix*lift, hidden*matrix*hidden))


def field_unitary(mode):
    if type(mode) is not str or mode not in ('identity', 'hadamard', 'phase'):
        raise ValueError('a declared field intervention required')
    logical = (sp.eye(2) if mode == 'identity' else
               sp.Matrix([[1, 1], [1, -1]]) if mode == 'hadamard' else sp.diag(1, -1))
    weight = F(1, 2) if mode == 'hadamard' else F(1)
    return ((weight, sp.kronecker_product(logical, sp.eye(2)).tolist()),)


def field_readout():
    result = []
    for sign in (1, -1):
        projector = sp.kronecker_product(sp.Matrix([[1, sign], [sign, 1]])/2, sp.eye(2))
        result.append(((F(1), projector.tolist()),))
    return tuple(result)


def intervention_groups(feedback):
    if type(feedback) is not bool:
        raise ValueError('Boolean record feedback selection required')
    # The declared observer acts after every recorded exit from clock port0.
    return tuple(field_unitary('hadamard' if feedback and mark == (-1, 0) else 'identity')
                 for mark in MARKS)


@lru_cache(maxsize=2)
def event_maps(law):
    maps = source_maps(_law(law))
    q0 = maps[0]
    resolvent = (sp.eye(16)-q0).inv()
    return tuple(sp.ImmutableMatrix(mark*resolvent) for mark in maps[1:])


def event_certificate(law):
    maps, stopped = source_maps(_law(law)), event_maps(law)
    trace = sp.Matrix([[int(i//4 == i % 4) for i in range(16)]])
    if trace*maps[0] != sp.Rational(5, 6)*trace:
        raise ArithmeticError('joint no-event survival is not5/6 on every input')
    if trace*sum(stopped, sp.zeros(16)) != trace:
        raise ArithmeticError('completed joint event instrument lost mass')
    destinations = (1, 0, 3, 2)
    effects = []
    for event, destination in zip(stopped, destinations, strict=True):
        effect = trace*event
        expected = sp.zeros(16)
        expected[5*destination, :] = effect
        if event != expected or any(effect[i] for i in range(16) if i//4 != i % 4):
            raise ArithmeticError('source completed mark lacks its exact measure-prepare form')
        effects.append(tuple(effect[5*i] for i in range(4)))
    transition = sp.Matrix(4, 4, lambda i, j: effects[i][destinations[j]])
    if any(sum(transition[:, j]) != 1 for j in range(4)):
        raise ArithmeticError('joint event chain incomplete')
    stationary = sp.ones(4, 1)/4
    if transition*stationary != stationary:
        raise ArithmeticError('source event reference is not stationary')
    switch = sum(transition[i, j]/4 for i, j in product(range(4), repeat=2)
                 if (i < 2) != (j < 2))
    local_survival = []
    for clock in (0, 1):
        mask = sp.diag(*(int((i < 2) != (clock == 0)) for i in range(4)))
        survival = sp.ones(1, 4)*(mask*transition)**2
        local_survival.append(max(survival))
    activity = sp.Rational(1, 7)
    held = (1-activity)*sp.eye(16)+activity*maps[0]
    if any(activity*mark*(sp.eye(16)-held).inv() != event
           for mark, event in zip(maps[1:], stopped, strict=True)):
        raise ArithmeticError('irrelevant identity holds changed the stopped instrument')
    return {'candidate': law, 'one_attempt_no_event_effect': '5/6 I',
        'cap_no_event_probability': '(5/6)^cap', 'mean_attempts_per_union_event': '6',
        'mark_order': [list(v) for v in MARKS[1:]],
        'post_event_native_addresses': [geometry()[0][i] for i in destinations],
        'event_effects_on_initial_basis_states': [[str(v) for v in row] for row in effects],
        'transition_columns_given_previous_mark': [[str(v) for v in row] for row in transition.tolist()],
        'stationary_mark_probabilities': ['1/4']*4,
        'stationary_probability_consecutive_union_ticks_use_different_clocks': str(switch),
        'largest_no_local_tick_probability_over_two_following_union_events': list(map(str, local_survival)),
        'independent_identity_hold_activity_one_seventh_exactly_cancelled': True,
        'every_event_map_is_rank_one_output_measure_prepare': True,
        'arbitrary_reference_outputs_require_retained_conditional_reference_operators': True,
        'completed_event_instrument_preserves_initial_quantum_coherence': False,
        'complete_map_sha256': _digest([list(m) for m in stopped])}


def _digest(value):
    def encode(item):
        if isinstance(item, (list, tuple)):
            return [encode(v) for v in item]
        return str(item)
    return hashlib.sha256(canonical(encode(value)).encode()).hexdigest()


def _dense(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def _apply_exact(matrix, state, reference=1):
    if type(reference) is not int or not 1 <= reference <= 4:
        raise ValueError('a bounded exact reference dimension required')
    _matrix(matrix, 16); _matrix(state, 4*reference)
    result = sp.zeros(4*reference)
    for r, s in product(range(reference), repeat=2):
        block = sp.Matrix(4, 4, lambda i, j: state[reference*i+r, reference*j+s])
        output = unvector(matrix*vector(block))
        for i, j in product(range(4), repeat=2):
            result[reference*i+r, reference*j+s] = output[i, j]
    return result


def _reduce_reference(state, reference):
    _matrix(state, 4*reference)
    return sp.Matrix(2*reference, 2*reference, lambda i, j:
        sum(state[(2*(i//reference)+h)*reference+i % reference,
                  (2*(j//reference)+h)*reference+j % reference] for h in range(2)))


def initial_state(reference=False, phase=1):
    if type(reference) is not bool or type(phase) is not int or phase not in (-1, 1):
        raise ValueError('Boolean reference and exact preparation phase required')
    rd = 2 if reference else 1
    state = sp.zeros(4*rd)
    first, second = 0, 2*rd+int(reference)
    state[first, first] = state[second, second] = sp.Rational(1, 2)
    state[first, second] = state[second, first] = sp.Rational(phase, 2)
    return state


def experiment(law, *, feedback=False, reference=False, cap=3, phase=1):
    """Complete fixed-cap execution; unresolved event waits remain in its law.

    Native attempts are a resource cap. The two clock counters are derived
    only from retained signed repair records. An intervention is applied on
    each clock-A exit; it does not reset or replace a branch's density.
    """
    _law(law); _cap(cap)
    if cap > 3:
        raise ValueError('actual complete runtime allocation allows cap at most3')
    initial = initial_state(reference, phase)
    rd = 2 if reference else 1
    interventions = intervention_groups(feedback)
    maps = source_maps(law)
    combined = tuple(superoperator(group)*matrix for group, matrix in zip(interventions, maps, strict=True))
    reads = tuple(superoperator(group) for group in field_readout())
    events = []
    for attempt in range(cap):
        events.append(runtime.Event(f'native_joint_clock_{attempt}', instruments=(source_instrument(law),)))
        events.append(runtime.Event(f'local_exit_feedback_{attempt}', control=attempt,
            instruments=tuple(q._runtime_instrument((group,),
                resources=(('declared_field_interventions', int(feedback and index == 1)),),
                binding={'mode': 'field Hadamard' if feedback and index == 1 else 'identity',
                         'trigger': list(MARKS[index])}) for index, group in enumerate(interventions))))
    events.append(runtime.Event('terminal_field_X', instruments=(q._runtime_instrument(field_readout(),
        visible=True, resources=(('declared_field_readouts', 1),),
        binding={'algebra': 'carrier qubit tensored with within-carrier identity',
                 'outcome_signs': [1, -1]}),)))
    default = sp.Matrix(4, 4, lambda i, j: sum(initial[rd*i+r, rd*j+r] for r in range(rd)))
    program = runtime.Program(f'joint_clock_{law}_{feedback}_{reference}_{cap}_{phase}', 4, _dense(default),
        tuple(events), preparation_resources=(('declared_coherent_preparations', 1),),
        binding={'candidate': law, 'active_one_particle_ports': geometry()[0],
                 'all_other236_source_ports_vacuum': True, 'native_attempt_cap': cap,
                 'timeout_retained_when_fewer_than_requested_ticks': True,
                 'native_clock_events_exclude_declared_external_apparatus': True,
                 'field_and_hidden_factorization': 'ordered active addresses=(carrier bit, internal bit)'})
    exact = {(): initial}
    for _ in range(cap):
        exact = {word+(mark,): _apply_exact(matrix, state, rd)
                 for word, state in exact.items() for mark, matrix in enumerate(combined)}
    exact = {word+(outcome,): _apply_exact(matrix, state, rd)
             for word, state in exact.items() for outcome, matrix in enumerate(reads)}
    return program, initial, exact, combined, reads


def execute_experiment(law, *, feedback=False, reference=False, cap=3, phase=1):
    program, initial, exact, maps, reads = experiment(law, feedback=feedback,
        reference=reference, cap=cap, phase=phase)
    supplied = _dense(initial) if reference else None
    actual = runtime.execute(program, initial=supplied)
    prefix = runtime.execute(program, initial=supplied, stop_after=2)
    checkpoint = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program, checkpoint=checkpoint)
    residual = runtime.complete_distance(actual, continued)
    rd, rows, matrices = actual.reference_dimension, [], []
    plus, first_timeout, second_timeout, total = (sp.Rational(0) for _ in range(4))
    for branch in actual.branches:
        expected = exact[branch.record]
        residual = max(residual, float(np.linalg.norm(branch.density-_dense(expected), 'fro')))
        probability = sp.trace(expected)
        total += probability
        plus += probability*int(branch.record[-1] == 0)
        marks = branch.record[:-1]
        events = sum(mark != 0 for mark in marks)
        first_timeout += probability*int(events < 1)
        second_timeout += probability*int(events < 2)
        counts = tuple(sum(MARKS[mark][clock] != 0 for mark in marks) for clock in (0, 1))
        chosen = tuple(maps[mark] for mark in marks)+(reads[branch.record[-1]],)
        predicted = sp.zeros(2*rd)
        for r, s in product(range(rd), repeat=2):
            block = sp.Matrix(4, 4, lambda i, j: initial[rd*i+r, rd*j+s])
            reduced_block = memory_forecast(block, chosen)[-1]
            for i, j in product(range(2), repeat=2):
                predicted[rd*i+r, rd*j+s] = reduced_block[i, j]
        if predicted != _reduce_reference(expected, rd):
            raise ArithmeticError('word-dependent field memory lost a quantum output')
        expected_resources = {'native_attempts': cap, 'joint_clock_reads': cap,
            'declared_field_readouts': 1,
            'declared_field_interventions': int(feedback)*marks.count(1)}
        if not reference:
            expected_resources['declared_coherent_preparations'] = 1
        if dict(branch.resources) != expected_resources:
            raise ArithmeticError('complete branch resource accounting changed')
        rows.append({'attempt_marks_and_readout': list(branch.record),
                     'clock_tick_counts': list(counts), 'probability': str(probability),
                     'first_event_timeout': events < 1, 'second_event_timeout': events < 2})
        matrices.append([list(branch.record), list(expected)])
    if total != 1 or residual > 2e-11 or len(rows) != 2*5**cap:
        raise ArithmeticError('complete runtime history failed its exact oracle')
    return {'candidate': law, 'record_controlled_field_Hadamard': feedback,
        'initial_coherent_phase': phase, 'reference_dimension': rd,
        'attempt_cap': cap, 'complete_record_count': len(rows),
        'complete_field_plus_probability': str(plus),
        'first_event_timeout_probability': str(first_timeout),
        'second_event_timeout_probability': str(second_timeout),
        'complete_probability': '1', 'complete_record_law': rows,
        'exact_all_terminal_work_reference_matrices_sha256': _digest(matrices),
        'exact_memory_matches_every_reduced_terminal_matrix': True,
        'checkpoint_after_first_attempt_and_feedback_replayed': True,
        'full_state_and_checkpoint_residual_diagnostic': residual,
        'preparation_resource_scope': ('external joint reference preparation unpriced' if reference
            else 'one declared coherent preparation'),
        'compiled_pulse_or_hardware_service_error_certified': False}


def memory_witnesses():
    maps = source_maps('C')
    reduce, lift = field_maps()
    first, second = sp.zeros(4), sp.zeros(4)
    first[0, 0], second[1, 1] = 1, 1
    if reduced(first) != reduced(second):
        raise ArithmeticError('hidden-preparation witness no longer has equal fields')
    event_laws = [[sp.trace(unvector(event*vector(initial))) for event in event_maps('C')]
                  for initial in (first, second)]
    initial = unvector(lift*sp.Matrix([1, 0, 0, 0]))
    chosen = (maps[0], maps[1])
    exact = reduced(unvector(chosen[1]*chosen[0]*vector(initial)))
    memory = memory_forecast(initial, chosen)[-1]
    memoryless = sp.Matrix(2, 2, list((reduce*chosen[1]*lift)*(reduce*chosen[0]*lift)*(reduce*vector(initial))))
    if memory != exact or exact == memoryless:
        raise ArithmeticError('the nonzero field memory witness failed')
    return {'same_initial_reduced_field': [[str(v) for v in row] for row in reduced(first).tolist()],
        'different_native_preparations': [0, 1],
        'completed_first_mark_laws': [[str(v) for v in row] for row in event_laws],
        'first_event_law_total_variation': str(sum(abs(a-b) for a, b in zip(*event_laws, strict=True))/2),
        'zero_initial_hidden_term_control_word': ['no event', 'clock_A_exit'],
        'exact_word_probability': str(sp.trace(exact)),
        'memoryless_word_probability': str(sp.trace(memoryless)),
        'exact_memory_restores_full_reduced_matrix': True,
        'memory_is_not_an_assumed_finite_window_or_CP_kernel': True}


def run():
    from .clock_change import run as clock_change_run
    executed = [execute_experiment(law, feedback=feedback)
                for law, feedback in product(('C', 'E'), (False, True))]
    executed.extend(execute_experiment(law, reference=True, cap=2) for law in ('C', 'E'))
    return {'schema': 'oph.joint_clock_process.v1',
        'event_closure': [event_certificate(law) for law in ('C', 'E')],
        'field_memory_witnesses': memory_witnesses(),
        'observer_clock_conversion': clock_change_run(),
        'complete_runtime_experiments': executed,
        'scope': {'finite_observer_chosen_field_and_two_local_record_clocks': True,
                  'time_parameter': 'native attempt resource cap and local signed-event counts',
                  'common_physical_clock_or_modular_identification': False,
                  'new_relativistic_spacetime_or_empirical_confirmation': False,
                  'stationary_state_conserved_sector_schedule_and_apparatus_selected_from_OPH': False},
        'source_bindings': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS})}
