"""Bind a finite adaptive instrument tree to the actual source operations.

All returned maps act on the full field input. Abort maps retain their own
output spaces and are never coherently combined or success-normalized.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import json
import math
import numpy as np
from scipy.linalg import expm

from model.full_fock_field.field import FieldGeometry, FullFock, source_word
from model.full_fock_field.local_simulation import split_bound
from model.field_thermal.thermal import execute_with_exports
from model.observer.source_scalar import compiler
from model.reference_observation.observation import frozen


@dataclass(frozen=True)
class Read:
    duration: F
    mode: int
    steps: int = 4

    def __post_init__(self):
        if (type(self.duration) not in (int, F) or type(self.mode) is not int
                or self.mode < 0 or type(self.steps) is not int or not 1 <= self.steps <= 16
                or abs(self.duration) > 4):
            raise ValueError('bounded rational duration, mode and split count required')


def _geometry(geometry):
    if not isinstance(geometry, FieldGeometry) or geometry.volume not in (1, 2):
        raise ValueError('one or two complete field sites required')


def instrument(geometry, read):
    """Execute on all input matrix columns, thereby covering any reference."""
    _geometry(geometry)
    if not isinstance(read, Read) or read.mode >= geometry.modes:
        raise ValueError('read mode outside the complete field')
    d = 1 << geometry.modes
    word = list(source_word(geometry, float(read.duration), read.steps))
    word += compiler.cnot(read.mode, geometry.modes)
    blank_pointer = np.zeros((2*d, d), complex)
    blank_pointer[:d] = np.eye(d)
    output, failures = execute_with_exports(word, geometry.modes+1, blank_pointer)
    space = FullFock(geometry.modes)
    unitary = expm(-1j*float(read.duration)*space.hamiltonian(geometry).toarray())
    number = space.number(read.mode).toarray()
    actual = (frozen(output[:d]), frozen(output[d:]))
    target = (frozen((np.eye(d)-number)@unitary), frozen(number@unitary))
    # Each event and complement basis value has a separate Kraus map.
    aborts = tuple((event, tuple(branches)) for event, branches in enumerate(failures))
    return actual, target, aborts, len(word)


def adaptive(geometry=FieldGeometry((2, 1, 1)), first=Read(F(1,4), 0),
             second=(Read(F(1,4), 1), Read(F(-1,4), 1)), *, diagnostic_only=False):
    """Read once, select the next signed evolution using that actual bit, read again.

    Two durable record bits suffice. A fresh blank pointer is admitted for
    each read; the first is retained as a classical record, never erased.
    Every first-stage abort terminates; second-stage aborts retain the first
    record. Execution cost is one branch, not the sum of counterfactual trees.
    """
    _geometry(geometry)
    if type(diagnostic_only) is not bool:
        raise ValueError('explicit Boolean diagnostic mode required')
    if type(second) is not tuple or len(second) != 2 or any(not isinstance(r, Read) for r in second):
        raise ValueError('exactly one declared continuation per first outcome required')
    a, b, failures, initial_cost = instrument(geometry, first)
    aborts = {((), event): branches for event, branches in failures}
    actual, target, costs = {}, {}, []
    for bit, read in enumerate(second):
        c, d, stopped, cost = instrument(geometry, read)
        costs.append(initial_cost+cost)
        for outcome in range(2):
            actual[bit, outcome] = frozen(c[outcome]@a[bit])
            target[bit, outcome] = frozen(d[outcome]@b[bit])
        for event, branches in stopped:
            aborts[((bit,), event)] = tuple(frozen(k@a[bit]) for k in branches)
    dim = 1 << geometry.modes
    failure_gram = sum((k.conj().T@k for branches in aborts.values() for k in branches),
                       np.zeros((dim, dim), complex))
    success_gram = sum(k.conj().T@k for k in actual.values())
    difference_gram = failure_gram+sum((actual[key]-target[key]).conj().T@
                                      (actual[key]-target[key]) for key in actual)
    tp_residual = float(np.linalg.norm(success_gram+failure_gram-np.eye(dim), 2))
    if not math.isfinite(tp_residual) or tp_residual > 2e-9:
        raise ArithmeticError('complete source instrument lost input/reference probability')
    error = math.sqrt(float(np.linalg.norm(difference_gram, 2)))
    # At each adaptive layer, direct-sum record control costs the maximum
    # branch error, rather than adding mutually exclusive continuations.
    def bound(read):
        return split_bound(geometry.sides, read.duration, read.steps,
                           F(geometry.speed), F(geometry.interaction), refined=True)
    dilation_bound = bound(first)+max(map(bound, second))
    analytic = min(F(1), dilation_bound)
    admitted = math.isfinite(error) and error <= float(dilation_bound)+2e-9
    if not math.isfinite(error) or (not admitted and not diagnostic_only):
        raise ArithmeticError('complete adaptive source exceeds its ideal-word error budget')
    return {'success': actual, 'ideal': target, 'aborts': aborts,
            'failure_gram': frozen(failure_gram), 'trace_preservation_residual': tp_residual,
            'numerical_dilation_difference_norm': error,
            'ideal_word_half_diamond_bound': analytic,
            'ideal_source_error_budget_satisfied': admitted,
            'diagnostic_only': diagnostic_only,
            'maximum_native_pulses_per_history': max(costs),
            'maximum_load_export_events_per_history': 2*max(costs),
            'field_modes': geometry.modes, 'record_bits': 2,
            'initialized_pointer_bits_per_history': 2}


def report():
    tree = adaptive()
    d = 1 << tree['field_modes']
    # Maximally entangled field/reference input; total norms include every abort.
    probabilities = {str(key): float(np.linalg.norm(k, 'fro')**2/d)
                     for key, k in tree['success'].items()}
    return {key: value for key, value in tree.items() if key not in
            {'success', 'ideal', 'aborts', 'failure_gram', 'ideal_word_half_diamond_bound'}} | {
        'schema': 'oph.complete_adaptive_field.v1',
        'ideal_word_half_diamond_bound': str(tree['ideal_word_half_diamond_bound']),
        'entangled_input_joint_probabilities': probabilities,
        'entangled_input_abort_probability': float(np.trace(tree['failure_gram']).real/d),
        'complete_abort_slots': len(tree['aborts']),
        'outcomes': 4,
        'domain': 'complete two-site CAR algebra, arbitrary input and untouched reference',
        'clock': 'Hamiltonian model coordinate; first +1/4, then +1/4 or -1/4 from the first record',
        'physical_3_plus_1_field_identified': False,
        'numerical_error_is_outward_certificate': False,
    }


if __name__ == '__main__':
    print(json.dumps(report(), indent=2))
