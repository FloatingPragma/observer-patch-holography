"""All-record native memory and an executed nonlinear observation process.

Fine labels are grouped only when their new memory and field operation agree.
Their exact conditional label probabilities are retained, so classical splitting
recovers the full record law without changing its trace distance.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s
from scipy import linalg as la

from ..common_source import native_family_record_memory as memory
from ..common_source import native_repair_mass as repair
from ..compact_interacting_response.finite import numerical_block
from . import bounds


P, Q = F(1, 32), F(1, 16)
DURATION = F(1, 32)


@lru_cache(maxsize=1)
def law():
    """Literal 128 source records times four independent declared binary filters."""
    records = memory._validated_records(repair.native_records())
    _, (gc, ge) = memory.geometry()
    binary = ((1-12*P)*(1-6*Q), (1-12*P)*6*Q, 12*P*(1-6*Q), 72*P*Q)
    transition = s.zeros(12)
    enabled = s.zeros(12)
    classes = []
    for old in range(12):
        grouped = {}
        for index, (cc, cr, ec, er, _, raw) in enumerate(records):
            for left, right in product((False, True), repeat=2):
                active = cc == 0 and cr == 1 and left
                event_e = ec == 0 and er in (4, 5) and right
                new = memory._destination(old, active, event_e, 2, gc, ge)
                weight = s.Rational(raw)*s.Rational(binary[2*int(left)+int(right)])
                label = 4*index+2*int(left)+int(right)
                grouped.setdefault((new, active), []).append((label, weight))
                transition[new, old] += weight
                if active:
                    enabled[new, old] += weight
        rows = []
        for (new, active), labels in sorted(grouped.items()):
            weight = sum((weight for _, weight in labels), s.Integer(0))
            rows.append({'new': new, 'enabled': active, 'weight': weight, 'labels': tuple(labels)})
        if sorted(label for row in rows for label, _ in row['labels']) != list(range(512)):
            raise ArithmeticError('the complete 512-label alphabet was not retained')
        if sum(row['weight'] for row in rows) != 1:
            raise ArithmeticError('native record weights do not sum to one')
        classes.append(tuple(rows))
    declared = memory.architecture_from_factors(s.Rational(P), s.Rational(Q))
    if transition != declared['T']:
        raise ArithmeticError('literal native records differ from the declared memory recursion')
    stationary = declared['stationary']
    reset = s.zeros(12)
    reset[2, :] = s.ones(1, 12)
    remainder = (transition-(1-s.Rational(P))*reset)/s.Rational(P)
    if (any(x < 0 for x in remainder) or any(sum(remainder[:, j]) != 1 for j in range(12))
            or transition*stationary != stationary):
        raise ArithmeticError('the actual reset decomposition or invariant law failed')
    return {'classes': tuple(classes), 'transition': s.ImmutableMatrix(transition),
            'enabled': s.ImmutableMatrix(enabled), 'stationary': stationary,
            'projectors': declared['projectors']}


def law_certificate():
    data = law()
    pi = data['stationary']
    variation = sum(abs(x) for x in data['transition'][:, 2]-pi)/2
    if variation == 0:
        raise ArithmeticError('resetting a continuing memory lost its countercontrol')
    return {'old_memory_inputs': 12, 'fine_labels_per_input': 512,
            'native_raw_labels': 128, 'independent_binary_filter_labels': 4,
            'event_factors': [str(P), str(Q)],
            'event_factor_origin': 'unchanged native records plus declared independent binary filters',
            'binary_filter_probabilities': ['3/8', '3/8'],
            'transition': [[str(x) for x in data['transition'][i, :]] for i in range(12)],
            'enabled_transition': [[str(x) for x in data['enabled'][i, :]] for i in range(12)],
            'stationary': [str(x) for x in pi],
            'record_class_counts': [len(rows) for rows in data['classes']],
            'zero_weight_fine_labels_per_input':
                [sum(weight == 0 for row in rows for _, weight in row['labels']) for rows in data['classes']],
            'conditional_reset_port_next_step_vs_stationary_TV': str(variation),
            'fine_labels_reconstructible_without_conditioning': True}


def probe():
    """A complete four-mode nonlinear probe, not a selected compact-core slice."""
    fixed, interaction, annihilators = numerical_block(1/4, 1/2, 1/64, 3/128)
    pairing, _, _ = numerical_block(1., 0., 0., 0.)
    data = law()
    coefficients = [float(projector[0, 0]) for projector in data['projectors']]
    gain = .5
    pulses = [fixed+gain*x*pairing for x in coefficients]
    weights = data['enabled']*data['stationary']
    mean = fixed+gain*sum(float(x)*y for x, y in zip(weights, coefficients, strict=True))*pairing
    n0 = annihilators[0].T@annihilators[0]
    pair01 = annihilators[0].T@annihilators[1].T+annihilators[1]@annihilators[0]
    rotate = la.expm(-1j*np.pi*pair01/4)
    feedback = rotate
    final = rotate@n0@rotate.conj().T
    if (la.norm(fixed@interaction-interaction@fixed) == 0
            or la.norm(n0@mean-mean@n0) == 0 or la.norm(n0@final-final@n0) < .1):
        raise ArithmeticError('a required interaction or noncommuting detector was lost')
    # Exact triangle estimate: pairing norm <=2, interaction norm <=4,
    # kinetics norm <=1/64+3/128, gain<=1/2, projector coefficient<=1.
    norm = F(1, 2)+2+F(5, 128)+1
    if norm >= 4 or max(la.norm(h, 2) for h in [fixed, mean, *pulses]) > float(norm)+1e-12:
        raise ArithmeticError('the independently bounded probe norm failed')
    return {'fixed': fixed, 'pulses': pulses, 'mean': mean,
            'middle': (np.eye(16)-n0, n0), 'final': (np.eye(16)-final, final),
            'feedback': feedback, 'norm': F(4)}


def initial_vectors():
    """Classically correlated memory plus a retained entangled two-row reference."""
    a = np.zeros((2, 32), complex)
    a[0, 0], a[0, 2*5+1] = .5, .5
    a[1, 2*3], a[1, 2*12+1] = .5, .5j
    return a, np.array([2, 10], dtype=np.int64)


def apply(matrix, vectors):
    return np.einsum('ij,bjr->bir', matrix, vectors.reshape(-1, 16, 2)).reshape(-1, 32)


def advance(actual, ideal, labels, active, outcomes, step, operators):
    """One complete memory instrument; never reset memory at a detector boundary."""
    off = la.expm(-1j*step*operators['fixed'])
    on = [la.expm(-1j*step*h) for h in operators['pulses']]
    target = la.expm(-1j*step*operators['mean'])
    arrays = [[], [], [], [], []]
    for old, classes in enumerate(law()['classes']):
        keep = labels == old
        if not np.any(keep):
            continue
        for row in classes:
            # Zero-probability classes are defined by the exact label grammar.
            # Their zero output does not require allocation of an amplitude row.
            if row['weight'] == 0:
                continue
            factor = np.sqrt(float(row['weight']))
            arrays[0].append(factor*apply(on[row['new']] if row['enabled'] else off, actual[keep]))
            arrays[1].append(factor*apply(target, ideal[keep]))
            arrays[2].append(np.full(np.sum(keep), row['new'], dtype=np.int64))
            arrays[3].append(np.full(np.sum(keep), row['enabled'], dtype=bool))
            arrays[4].append(outcomes[keep])
    return tuple(np.concatenate(array) for array in arrays)


def detector(actual, ideal, labels, active, outcomes, operators, *, final=False):
    arrays = [[], [], [], [], []]
    for bit, projector in enumerate(operators['final'] if final else operators['middle']):
        left, right = apply(projector, actual), apply(projector, ideal)
        if not final:
            # Feedforward uses a retained detector bit and the most recent
            # actual native enabling record. The same branch policy acts on
            # the ideal comparison, including branches of zero probability.
            selected = active ^ bool(bit)
            left[selected] = apply(operators['feedback'], left[selected])
            right[selected] = apply(operators['feedback'], right[selected])
        arrays[0].append(left)
        arrays[1].append(right)
        arrays[2].append(labels.copy())
        arrays[3].append(active.copy())
        arrays[4].append(2*outcomes+bit)
    return tuple(np.concatenate(array) for array in arrays)


def execute(steps=2):
    if type(steps) is not int or not 1 <= steps <= 4:
        raise ValueError('the explicit complete-record evaluator admits one through four steps per interval')
    operators = probe()
    actual, labels = initial_vectors()
    ideal = actual.copy()
    active = np.zeros(len(labels), bool)
    outcomes = np.zeros(len(labels), np.int64)
    for interval in range(2):
        for _ in range(steps):
            actual, ideal, labels, active, outcomes = advance(
                actual, ideal, labels, active, outcomes, float(DURATION/steps), operators)
        actual, ideal, labels, active, outcomes = detector(
            actual, ideal, labels, active, outcomes, operators, final=interval == 1)
    norms_a = np.sum(abs(actual)**2, axis=1)
    norms_b = np.sum(abs(ideal)**2, axis=1)
    overlaps = np.sum(actual.conj()*ideal, axis=1)
    discriminants = (norms_a+norms_b)**2-4*abs(overlaps)**2
    if np.min(discriminants) < -1e-12:
        raise ArithmeticError('the unnormalized pure-output trace distance is inconsistent')
    error = float(np.sum(np.sqrt(np.maximum(discriminants, 0))))
    if max(abs(norms_a.sum()-1), abs(norms_b.sum()-1)) > 1e-11:
        raise ArithmeticError('complete intervention history lost probability')
    one = bounds.segment(operators['norm'], DURATION, P, steps)
    bound_squared = min(F(4), 4*one['complete_record_trace_norm_upper_squared'])
    if error > np.sqrt(float(bound_squared))+1e-10:
        raise ArithmeticError('complete adaptive observation exceeded the analytic bound')
    probabilities = [float(np.sum(norms_a[outcomes == label])) for label in range(4)]
    if min(probabilities) <= 0:
        raise ArithmeticError('an intervention outcome was lost')
    return {'steps_per_interval': steps, 'duration_per_interval': str(DURATION),
            'retained_class_history_rows': len(labels),
            'complete_record_trace_distance': error,
            'analytic_trace_bound_squared': str(bound_squared),
            'detector_probabilities': probabilities,
            'normalization_error': float(max(abs(norms_a.sum()-1), abs(norms_b.sum()-1)))}


def certificate():
    rows = [execute(steps) for steps in (1, 2, 4)]
    errors = [row['complete_record_trace_distance'] for row in rows]
    if not errors[0] > errors[1] > errors[2] > 0:
        raise ArithmeticError('the complete-history refinement control failed')
    # Store outward decimal bins, not platform-dependent raw eigensolver bits.
    # Each bin is deliberately wider than binary64 arithmetic discrepancies.
    summaries = []
    for row in rows:
        error = row['complete_record_trace_distance']
        lower = int(np.floor(error*10000))
        summaries.append({'steps_per_interval': row['steps_per_interval'],
                          'duration_per_interval': row['duration_per_interval'],
                          'retained_class_history_rows': row['retained_class_history_rows'],
                          'complete_record_trace_distance_observed_interval':
                              [str(F(lower, 10000)), str(F(lower+1, 10000))],
                          'normalization_error_below': '1e-11',
                          'all_four_detector_outcomes_positive': True,
                          'analytic_trace_bound_squared': row['analytic_trace_bound_squared']})
    return {'law': law_certificate(), 'executions': summaries,
            'field_rows': 16, 'untouched_reference_rows': 2,
            'initial_memory_classically_correlated': True,
            'native_memory_reset_between_intervals': False,
            'middle_detector_and_final_detector_commute': False,
            'feedback_depends_on_native_and_detector_record': True,
            'unread_bound_can_replace_complete_history_bound': False,
            'probe_is_a_reducing_compact_core_sector': False,
            'source_native_detector_compilation_claimed': False,
            'fine_record_factorization_is_lossless': True,
            'numerical_outputs_are_interval_proofs': False}


def archive_countercontrol():
    """The unread bound and complete-archive bound differ even for two exact branches."""
    target = s.diag(1, 0)
    vectors = [s.Matrix([s.Rational(4, 5), sign*s.Rational(3, 5)]) for sign in (-1, 1)]
    densities = [vector*vector.H for vector in vectors]
    unread = sum(densities, s.zeros(2))/2-target
    trace_norm = lambda matrix: sum(abs(value)*multiplicity
                                    for value, multiplicity in matrix.eigenvals().items())
    delta = trace_norm(unread)
    full = sum(trace_norm(rho-target) for rho in densities)/2
    if delta != s.Rational(18, 25) or full != s.Rational(6, 5) or full*full != 2*delta:
        raise ArithmeticError('the complete-record square-root control failed')
    return {'unread_trace_distance': str(delta), 'complete_record_trace_distance': str(full),
            'square_root_bound_saturated': True, 'unread_rate_is_complete_rate': False}
