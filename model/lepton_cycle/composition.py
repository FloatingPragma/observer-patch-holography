"""Complete repeated cycles, finite reference errors and declared service counts."""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from .instrument import LABELS, P, TRACE_EFFECT


def _rounds(n):
    if type(n) is not int or n < 1:
        raise ValueError('positive integer number of rounds required')
    return n


def history_multipliers(grams, n):
    """All formal histories, including zero branches, for a bounded replay."""
    _rounds(n)
    if n > 6 or set(grams) != set(LABELS):
        raise ValueError('complete four-outcome instrument and at most six explicit rounds required')
    return {word: np.prod(np.stack([grams[label] for label in word]), axis=0)
            for word in product(LABELS, repeat=n)}


def complete_multiplier(grams, n):
    _rounds(n)
    if set(grams) != set(LABELS):
        raise ValueError('complete four-outcome instrument required')
    return sum(grams.values()) ** n


def error_budget(n, t=F(1), *, cycle='weak', per_cycle_channel_error=F(0), preparation_error=F(0)):
    """Upper bounds in half diamond norm, hence trace distance with any reference.

Implementation inputs are authenticated complete-channel errors, not inferred hardware accuracies.
Preparation error includes independence from the continuing register/reference, not only a source marginal.
"""
    _rounds(n)
    if any(not isinstance(v, F) or v < 0 for v in (t, per_cycle_channel_error, preparation_error)):
        raise ValueError('nonnegative exact time and error budgets required')
    if cycle == 'weak':
        ideal = P * (1 - P) * t * t / n
    elif cycle == 'mass':
        ideal = (1 + TRACE_EFFECT ** 2) * t * t / (4 * n)
    else:
        raise ValueError('weak or mass cycle required')
    implemented = min(F(1), ideal + n * (per_cycle_channel_error + preparation_error))
    return {'cycle': cycle, 'rounds': n, 'declared_time': str(t),
            'ideal_half_diamond_bound': str(min(F(1), ideal)),
            'implemented_half_diamond_bound': str(implemented),
            'per_cycle_channel_error': str(per_cycle_channel_error),
            'per_fresh_pair_preparation_error': str(preparation_error)}


def required_rounds(accuracy, t=F(1), *, cycle='weak'):
    if not isinstance(accuracy, F) or not 0 < accuracy < 1:
        raise ValueError('exact accuracy between zero and one required')
    coefficient = F(error_budget(1, t, cycle=cycle)['ideal_half_diamond_bound'])
    # Use the untruncated coefficient for long declared times.
    coefficient = (P * (1 - P) if cycle == 'weak' else (1 + TRACE_EFFECT**2) / 4) * t*t
    return max(1, math.ceil(coefficient / accuracy))


def resources(n):
    """Literal finite service counts; no compiled elementary-gate or energy claim."""
    _rounds(n)
    return {'cycles': n, 'fresh_prepared_pairs_including_terminal_active_pair': n + 1,
            'blank_double_loads_prepared': 2 * (n + 1),
            'declared_uniform_address_preparation_gates': 2 * (n + 1),
            'forward_native_cell_ticks': 2 * n, 'inverse_native_cell_ticks': 2 * n,
            'controlled_analyzer_or_exclusive_kicks': n, 'binary_return_tests': 2 * n,
            'cell_register_swaps': 2 * n, 'classical_record_bits': 2 * n,
            'exported_address_pairs': n, 'exported_address_hilbert_dimension_per_pair': 3600,
            'binary_embedding_qubits_per_exported_pair': 12,
            'formal_complete_history_count': f'4^{n}',
            'compiled_pulse_count_claimed': False, 'physical_duration_or_energy_cost_claimed': False}


def fixed_entry_mask(multiplier, tolerance=1e-12):
    """Numerical diagnostic of the exact entrywise fixed-algebra condition."""
    return np.abs(np.asarray(multiplier) - 1) <= tolerance


def implementation_design(accuracy, *, cycle='weak'):
    """Split a target error equally between ideal finite rounds and service errors."""
    n = required_rounds(accuracy / 2, cycle=cycle)
    allowance = accuracy / (2 * n)
    row = error_budget(n, cycle=cycle, per_cycle_channel_error=allowance)
    return {'target_trace_distance': str(accuracy), 'rounds': n,
            'maximum_combined_channel_and_fresh_pair_error_per_cycle': str(allowance),
            'certified_total_error': row['implemented_half_diamond_bound']}
