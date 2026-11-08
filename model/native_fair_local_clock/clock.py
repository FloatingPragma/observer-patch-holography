"""Admitted fixed interventions and uniform epochs with exact attempt budgets."""
from dataclasses import fields
from fractions import Fraction as F
from math import comb, factorial
import hashlib
import json

from model.native_program import interpreter as native
from model.native_local_history_limit import limit
from model.native_exchange_clock import exchange
from model.native_rooted_refinement import routing, attachment


def typed(value):
    """Type-sensitive serialization, including exact sparse query operators."""
    if type(value) is F:
        return ['fraction', str(value)]
    if type(value) is tuple:
        return ['tuple', [typed(x) for x in value]]
    if type(value) is list:
        return ['list', [typed(x) for x in value]]
    if type(value) is dict:
        rows = [(typed(k), typed(v)) for k, v in value.items()]
        return ['dict', sorted(rows, key=lambda x: json.dumps(x[0], sort_keys=True))]
    if type(value) is native.Record:
        return ['Record', [[f.name, typed(getattr(value, f.name))] for f in fields(value)]]
    if value is None or type(value) in (str, int, bool):
        return [type(value).__name__, value]
    raise ValueError('an exact typed source value required')


def digest(value):
    return hashlib.sha256(json.dumps(typed(value), sort_keys=True).encode()).hexdigest()


def _word(word):
    if type(word) is not tuple:
        raise ValueError('a literal tuple of fixed native repairs required')
    for item in word:
        if type(item) is not tuple or len(item) != 2:
            raise ValueError('each fixed operation is its ordered seam and law')
        edge, law = item
        if (type(edge) is not tuple or len(edge) != 2 or
                any(type(v) is not int or v < 0 for v in edge) or edge[0] >= edge[1]
                or type(law) is not str or law not in ('C', 'E', 'mixture')):
            raise ValueError('literal ordered endpoints and native record partition required')
    return word


def compile_plan(level, total, partition, words, laws):
    """Each bin starts with a fixed word; all remaining ticks are uniform.

    Every fixed operation consumes one native attempt. No source operation
    is compressed into a zero-duration native service.
    """
    native._integer(level); native._integer(total)
    if type(words) is not tuple or type(laws) is not tuple:
        raise ValueError('one literal fixed word and uniform law per bin required')
    seams = 630 * 4**level
    ticks = limit.aligned_clock(seams, partition)
    if len(words) != len(ticks)-1 or len(laws) != len(words):
        raise ValueError('clock bins, interventions and laws must have equal lengths')
    commands, intervals, epochs = [], [], []
    for index, (left, right, word, law) in enumerate(zip(ticks[:-1], ticks[1:], words, laws, strict=True)):
        _word(word)
        if type(law) is not str or law not in ('C', 'E', 'mixture'):
            raise ValueError('a fixed native uniform-epoch record partition required')
        if len(word) > right-left:
            raise ValueError('fixed interventions exceed their counted clock bin')
        position = left
        for ordinal, (edge, fixed_law) in enumerate(word):
            commands.append({'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(edge)},
                             'law': fixed_law, 'repeat': 1})
            intervals.append((position, position+1, index, 'fixed', ordinal))
            position += 1
        background = right-position
        if background:
            commands.append({'op': 'repair', 'selector': {'kind': 'uniform'},
                             'law': law, 'repeat': background})
            intervals.append((position, right, index, 'uniform', None))
        epochs.append({'bin': index, 'total_attempts': right-left,
                       'fixed_attempts': len(word), 'uniform_attempts': background,
                       'fixed_duration': F(len(word), seams),
                       'uniform_duration': F(background, seams),
                       'background_exchange_time': exchange.native_time(seams, background)})
    program = native.program(federation={'level': level, 'total': total}, biased=False,
                             commands=commands, capacity=1, overflow='export_oldest')
    if program['parameters']['repair_probability'] != '1':
        raise ArithmeticError('the declared clock requires a repair at every attempt')
    return {'schema': 'oph.native_fair_local_clock.plan.v1', 'level': level, 'total': total,
            'seams': seams, 'partition': partition, 'words': words, 'laws': laws,
            'program': program, 'program_digest': native.digest(program),
            'command_intervals': tuple(intervals), 'epochs': tuple(epochs),
            'total_attempts': ticks[-1], 'fixed_attempts': sum(map(len, words)),
            'uniform_attempts': sum(e['uniform_attempts'] for e in epochs)}


def authenticate_plan(value):
    if type(value) is not dict:
        raise ValueError('an authenticated native clock plan required')
    try:
        expected = compile_plan(*(value[k] for k in ('level', 'total', 'partition', 'words', 'laws')))
    except (KeyError, TypeError) as error:
        raise ValueError('the complete native plan fields required') from error
    if typed(value) != typed(expected):
        raise ValueError('clock plan differs from its literal source compilation')
    return True


def attempt_label(value, attempt):
    """Counted tick, bin and operation role; the public observation drops tick."""
    authenticate_plan(value); native._integer(attempt, 1, value['total_attempts'])
    for start, stop, index, role, ordinal in value['command_intervals']:
        if start < attempt <= stop:
            return {'attempt': attempt, 'time': F(attempt, value['seams']),
                    'bin': index, 'role': role, 'fixed_ordinal': ordinal}
    raise ArithmeticError('an admitted attempt lost its source operation')


def budget(value, window_seams, root_size, degree, radius, initial_error=F(0)):
    """Full selected-local CQ comparison, conditional on stated support data.

    Every fixed operation and every observed endpoint must lie in the root
    set. The window contains its graph radius-radius ball. Input error is a
    joint window/reference half-trace bound. No archive factorization follows.
    """
    authenticate_plan(value)
    native._integer(window_seams, 0, value['seams']); native._integer(root_size, 1)
    native._integer(degree, 1); native._integer(radius)
    limit.rational(initial_error, upper=F(1))
    M, U, L = value['seams'], value['uniform_attempts'], value['fixed_attempts']
    length = radius+1
    cone = min(F(1), root_size*comb(U, length)*F(degree, M)**length) if U >= length else F(0)
    categorical = min(F(1), F(U*window_seams**2, M**2))
    reservation = min(F(1), F(window_seams*L, M))
    infinite_cone = min(F(1), root_size*(degree*value['partition'][-1])**length/factorial(length))
    return {'initial_joint_half_trace_upper': initial_error,
            'uniform_attempt_cone_upper': cone,
            'complete_window_categorical_Poisson_upper': categorical,
            'fixed_attempt_clock_reservation_upper': reservation,
            'native_to_common_finite_window_upper': min(F(1), initial_error+cone+categorical+reservation),
            'additional_infinite_process_cone_if_supplied': infinite_cone,
            'native_total_attempts': value['total_attempts'],
            'native_uniform_attempts': U, 'native_fixed_attempts': L,
            'per_seam_time': value['partition'][-1],
            'uniform_exchange_times': tuple(e['background_exchange_time'] for e in value['epochs'])}


def window_plan(window, total, partition, words, laws):
    """Compile canonical window addresses to literal source seam addresses."""
    routing.authenticate_window(window)
    converted = []
    for word in words:
        _word(word); row = []
        for edge, law in word:
            if edge not in window['internal_seams']:
                raise ValueError('fixed word contains an edge outside the authenticated window')
            row.append((tuple(window['source_ports'][v] for v in edge), law))
        converted.append(tuple(row))
    return compile_plan(window['level'], total, partition, tuple(converted), laws)


def attach_window(value, window, root_cell_radius, observed_seams, final_ports,
                  *, initial_error=F(0), initial_reference=None):
    """Authenticate finite geometry and the exact observation support premise.

    Root cells and dynamical window are concentric addressed corner balls.
    Each port path crosses at most one cell per step, so their radius
    difference supplies the graph-ball premise without a physical metric.
    """
    authenticate_plan(value); routing.authenticate_window(window)
    native._integer(root_cell_radius, 0, window['cell_radius'])
    if value['level'] != window['level']:
        raise ValueError('window and programme name different source levels')
    if initial_reference is not None and (type(initial_reference) is not str or not initial_reference):
        raise ValueError('an explicit common joint-reference identifier required')
    inner = routing.routed_window(window['level'], root_cell_radius)
    names = {v: i for i, v in enumerate(window['source_ports'])}
    root = tuple(sorted(names[v] for v in inner['source_ports']))
    if (type(observed_seams) is not tuple or observed_seams != tuple(sorted(set(observed_seams)))
            or any(type(e) is not tuple or len(e) != 2 or any(type(v) is not int for v in e)
                   or e not in window['internal_seams'] or not set(e) <= set(root) for e in observed_seams)):
        raise ValueError('observed seams must have both endpoints in the declared root')
    if (type(final_ports) is not tuple or final_ports != tuple(sorted(set(final_ports)))
            or any(type(v) is not int or v not in root for v in final_ports)):
        raise ValueError('final ports must lie in the declared root')
    source_root = set(inner['source_ports'])
    if any(not set(edge) <= source_root for word in value['words'] for edge, _ in word):
        raise ValueError('every intervention endpoint must lie in the declared root')
    radius = window['cell_radius']-root_cell_radius
    return {'plan_digest': digest(value), 'window_digest': attachment.window_digest(window),
            'root_cell_radius': root_cell_radius, 'root_ports': root,
            'declared_initial_reference': initial_reference,
            'observed_uniform_seams': observed_seams, 'final_ports': final_ports,
            'graph_radius_margin': radius,
            'budget': budget(value, len(window['internal_seams']), len(root), 6, radius, initial_error)}


def canonical_error(value, window, density):
    authenticate_plan(value); routing.authenticate_window(window)
    if value['level'] != window['level']:
        raise ValueError('one source level required')
    return limit.canonical_state_bound(240*4**value['level'], value['total'],
                                       len(window['canonical_ports']), density)


def compare_attachments(left_plan, left_window, left_attachment,
                        right_plan, right_window, right_attachment):
    """A shared finite-window target requires the same transported instruments."""
    for plan, window, declared in ((left_plan, left_window, left_attachment),
                                   (right_plan, right_window, right_attachment)):
        expected = attach_window(plan, window, declared['root_cell_radius'],
            declared['observed_uniform_seams'], declared['final_ports'],
            initial_error=declared['budget']['initial_joint_half_trace_upper'],
            initial_reference=declared['declared_initial_reference'])
        if typed(expected) != typed(declared):
            raise ValueError('the observation attachment differs from exact replay')
    if not routing.same_window(left_window, right_window):
        raise ValueError('a common authenticated routed window required')
    if (left_attachment['declared_initial_reference'] is None or
            left_attachment['declared_initial_reference'] != right_attachment['declared_initial_reference']):
        raise ValueError('both initial errors must refer to the same explicitly declared joint reference')
    for key in ('partition', 'laws'):
        if typed(left_plan[key]) != typed(right_plan[key]):
            raise ValueError('the common clock and uniform record laws must coincide')
    def canonical_words(plan, window):
        inverse = {v: i for i, v in enumerate(window['source_ports'])}
        return tuple(tuple((tuple(inverse[v] for v in e), law) for e, law in word)
                     for word in plan['words'])
    if canonical_words(left_plan, left_window) != canonical_words(right_plan, right_window):
        raise ValueError('the fixed intervention instruments must transport identically')
    for key in ('root_cell_radius', 'root_ports', 'observed_uniform_seams', 'final_ports'):
        if typed(left_attachment[key]) != typed(right_attachment[key]):
            raise ValueError('the two observation maps must agree')
    return min(F(1), left_attachment['budget']['native_to_common_finite_window_upper']+
                    right_attachment['budget']['native_to_common_finite_window_upper'])


def prepared_readout_plan(level, total, blocks=32, preparation_time=F(4), readout_time=F(8)):
    """Actual fixed preparation followed by a full uniform E readout epoch.

    The preparation bin also includes its remaining uniform background
    attempts. No readout accuracy or local number distribution is inferred
    from declaring this common programme.
    """
    from model.native_fair_local_preparation import preparation
    limit.rational(preparation_time); limit.rational(readout_time)
    if not preparation_time or not readout_time:
        raise ValueError('positive preparation and readout bins required')
    word = tuple((edge, 'E') for edge in preparation.edges(blocks))
    plan = compile_plan(level, total, (F(0), preparation_time, preparation_time+readout_time),
                        (word, ()), ('mixture', 'E'))
    if typed(plan['program']['commands'][:len(word)]) != typed(preparation.commands(blocks)):
        raise ArithmeticError('the clock prefix differs from the owning preparation word')
    return plan


def readout_coverage_budget(value, epoch, tree_edges, erasure):
    """Coupon bound for a declared E-tree decoder on this uniform epoch.

    Its complete instrument attachment and initial constant-corner weight
    belong to the readout theorem. This routine does not infer those data
    from a source number or preparation error.
    """
    authenticate_plan(value); native._integer(epoch, 0, len(value['epochs'])-1)
    native._integer(tree_edges, 0, value['seams']); limit.rational(erasure, upper=F(1))
    if value['laws'][epoch] != 'E':
        raise ValueError('the contrast coverage decoder requires a uniform E epoch')
    U = value['epochs'][epoch]['uniform_attempts']
    whole_seam_times = U//value['seams']
    if whole_seam_times > 4096:
        raise ValueError('displayed exact coverage bound exceeds its finite arithmetic budget')
    uncovered = min(F(1), F(tree_edges, 2**whole_seam_times))
    return {'uniform_attempts': U, 'per_seam_duration': F(U, value['seams']),
            'whole_seam_times': whole_seam_times,
            'coverage_failure_upper': uncovered, 'initial_constant_corner_weight': erasure,
            'coverage_or_erasure_upper': uncovered+(1-uncovered)*erasure}
