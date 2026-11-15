"""Exact finite incident-clock control, with a supplied global-attempt cutoff.

The private scheduler chooses uniformly among three declared path seams. The
observer sees only signed currents at port zero and complete/timeout status.
Remote native maps keep acting; an incident zero current is still a tick.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import numpy as np

from . import construction as local

q, runtime = local.quantum, local.runtime
EDGES = ((0, 1), (1, 2), (2, 3))
SOURCE_PATHS = tuple(sorted(set((*local.SOURCE_PATHS,
    'model/native_action/stopped_control.py', 'model/primitive_observer/runtime.py'))))


def _cap(cap):
    local._integer(cap, 1, 4)
    return cap


def _schedule(schedule):
    if (type(schedule) is not tuple or not 1 <= len(schedule) <= 6
            or any(type(edge) is not int or edge not in (0, 1, 2) for edge in schedule)
            or schedule.count(0) > 2
            or (schedule.count(0) == 2 and schedule[-1] != 0)):
        raise ValueError('a bounded path tape stopped no later than its second incident attempt required')
    return schedule


def schedules(cap=3):
    """Private terminal scheduler tapes; their exact weights are 3**(-length)."""
    _cap(cap)
    waiting, finished = [()], []
    while waiting:
        word = waiting.pop()
        if len(word) == cap or word.count(0) == 2:
            finished.append(word)
        else:
            waiting.extend((*word, edge) for edge in (2, 1, 0))
    return tuple(finished)


def outcomes():
    """Complete observer alphabet, retaining impossible zero-output branches."""
    return tuple((status, word) for status, lengths in (('complete', (2,)), ('timeout', (0, 1)))
                 for length in lengths for word in product((-1, 0, 1), repeat=length))


def schedule_column(law, left, right, schedule):
    """Exact full-algebra instrument column for one private tape, including Z0 feedback."""
    q.candidate(law); local._integer(left, 0, 15); local._integer(right, 0, 15)
    _schedule(schedule)
    states = {((), left, right): F(1)}
    for edge in schedule:
        following = defaultdict(F)
        for (word, a, b), weight in states.items():
            for current in ((-1, 0, 1) if edge == 0 else (None,)):
                for (i, j), value in q.native_unit(law, a, b, EDGES[edge], current).items():
                    newword = (*word, current) if edge == 0 else word
                    sign = (-1)**((i & 1)+(j & 1)) if edge == 0 and not word and current == 0 else 1
                    following[newword, i, j] += weight*value*sign
        states = {key: value for key, value in following.items() if value}
    result = {word: {} for word in product((-1, 0, 1), repeat=schedule.count(0))}
    for (word, a, b), value in states.items():
        result[word][a, b] = value
    return result


def exact_column(law, left, right, cap=3):
    """All output matrix units; equality on these columns extends to every reference."""
    q.candidate(law); local._integer(left, 0, 15); local._integer(right, 0, 15)
    result = {record: defaultdict(F) for record in outcomes()}
    for tape in schedules(cap):
        status = 'complete' if tape.count(0) == 2 else 'timeout'
        for word, matrix in schedule_column(law, left, right, tape).items():
            for address, value in matrix.items():
                result[status, word][address] += value/F(3**len(tape))
    return {record: {key: value for key, value in matrix.items() if value}
            for record, matrix in result.items()}


def exact_probe(law, cap=3):
    """Root port and reference output for the exact existing Bell probe (four 1/2 entries)."""
    result = {record: defaultdict(F) for record in outcomes()}
    for (a, r), (b, s) in product(((8, 0), (9, 1)), repeat=2):
        for record, matrix in exact_column(law, a, b, cap).items():
            for (i, j), value in matrix.items():
                if i//2 == j//2:
                    result[record][2*(i & 1)+r, 2*(j & 1)+s] += value/2
    return {record: {key: value for key, value in matrix.items() if value}
            for record, matrix in result.items()}


def program(law, schedule):
    """Compile actual native instruments, then local feedback from the first current only."""
    q.candidate(law); _schedule(schedule)
    identity = runtime.Instrument(kraus=((np.eye(16),),))
    phase = runtime.Instrument(kraus=((np.diag([(-1.)**(mask & 1) for mask in range(16)]),),))
    events, hits = [], 0
    for index, edge in enumerate(schedule):
        instrument = local._instrument(local.full_groups(law, EDGES[edge]), visible=edge == 0)
        events.append(runtime.Event(f'native_{index}', instruments=(instrument,)))
        if edge == 0:
            hits += 1
            if hits == 1:
                events.append(runtime.Event('first_current_feedback', control=0,
                                            instruments=(identity, phase, identity)))
    initial = local._initial()
    reduced = np.trace(initial.reshape(16, 2, 16, 2), axis1=1, axis2=3)
    return runtime.Program(f'stopped_{law}_{schedule}', 16, reduced, tuple(events),
        binding={'private_path_tape': list(schedule), 'uniform_seam_count': 3,
                 'stop_after_incident_attempts': 2, 'root_port': 0,
                 'clock_counts_zero_currents': True,
                 'preparation_and_local_phase_service_supplied': True}), initial


def execute(law, cap=3):
    """Numerical complete local outputs, summed over the private uniform scheduler."""
    q.candidate(law); _cap(cap)
    result = {record: np.zeros((4, 4), complex) for record in outcomes()}
    for tape in schedules(cap):
        declared, initial = program(law, tape)
        history = runtime.execute(declared, initial=initial)
        status = 'complete' if tape.count(0) == 2 else 'timeout'
        for branch in history.branches:
            record = status, tuple(value-1 for value in branch.record)
            result[record] += local.partial_trace(branch.density, (0,), 4, 2)/3**len(tape)
    return result


def continuation(law):
    """Serialize an actual timeout, then extend its retained state and current record."""
    prefix, initial = program(law, (2, 1, 0))
    history = runtime.execute(prefix, initial=initial)
    restored = runtime.checkpoint_from_json(prefix, runtime.checkpoint_to_json(history))
    full, _ = program(law, (2, 1, 0, 2, 0))
    extended = runtime.branch_experiment(restored, full.events[len(prefix.events):],
        name='incident_timeout_continued', binding={'additional_attempts': 2})
    direct = runtime.execute(full, initial=initial)
    if direct.records != extended.records:
        raise ArithmeticError('continuation changed the complete record alphabet')
    # Sum entrywise absolute real/imaginary differences BEFORE rounding a
    # subtraction. Each matrix unit has nuclear norm one, so half this sum
    # is a rigorous complete-output trace-distance upper bound for the two
    # represented histories, not an eigensolver diagnostic.
    return sum((abs(F(float(x.real))-F(float(y.real)))
                +abs(F(float(x.imag))-F(float(y.imag)))
                for a, b in zip(direct.branches, extended.branches, strict=True)
                for x, y in zip(a.density.flat, b.density.flat, strict=True)), F(0))/2


def run():
    """Reproducible exact controls; numerical residual bounds apply to this fixed probe."""
    rows = []
    for law in ('C', 'E'):
        for left, right in product(range(16), repeat=2):
            column = exact_column(law, left, right)
            if sum(value for matrix in column.values() for (a, b), value in matrix.items()
                   if a == b) != int(left == right):
                raise ArithmeticError('complete stopped instrument lost trace')
        exact, actual = exact_probe(law), execute(law)
        error = sum(abs(F(float(actual[key][a, b].real))-exact[key].get((a, b), F(0)))
                    +abs(F(float(actual[key][a, b].imag)))
                    for key in outcomes() for a, b in product(range(4), repeat=2))/2
        timeout = sum(value for (status, _), matrix in exact.items() if status == 'timeout'
                      for (a, b), value in matrix.items() if a == b)
        if timeout != F(20, 27) or error > F(1, 10**12):
            raise ArithmeticError('complete probe output differs from its exact stopped law')
        resumed = continuation(law)
        if type(resumed) is not F or not 0 <= resumed <= F(1, 10**12):
            raise ArithmeticError('serialized timeout continuation changed the retained output')
        rows.append({'candidate': law, 'cap_global_attempts': 3, 'incident_ticks_required': 2,
            'carrier': 'declared four-port path; uniform scheduler on its three seams only',
            'observer_outcomes': len(outcomes()), 'private_terminal_tapes': len(schedules()),
            'all_input_matrix_units': 256, 'arbitrary_reference_extension': 'linearity of exact CP maps',
            'exact_timeout_mass': str(timeout), 'fixed_probe_half_nuclear_error_upper': '1/1000000000000',
            'serialized_timeout_continuation_error_upper': '1/1000000000000',
            'observer_sees_private_schedule': False, 'remote_native_updates_retained': True,
            'external_attempt_cutoff_supplied': True, 'local_phase_service_supplied': True,
            'phase_feedback_changes_local_root_output': law == 'C',
            'maximum_native_attempts': 3, 'maximum_feedback_interventions': 1,
            'maximum_retained_signed_currents': 2, 'untouched_reference_dimension': 2,
            'complete_local_output': [{'status': status, 'currents': list(word),
                'matrix': [[str(exact[status, word].get((a, b), F(0))) for b in range(4)]
                           for a in range(4)]} for status, word in outcomes()]})
    return rows
