"""Finite transfer reconstruction of the actual native stationary history.

The positive configuration kernel determines H=-log(K) in Euclidean native
attempt units. Marked transition insertions reproduce the observer's two
current bits and are checked against complete CQ runtime output. The analytic
continuation exp(-i*t*H) is not identified with native updates or physical time.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..observer import native
from ..common_source import slow_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json
from . import runtime


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    'model/primitive_observer/transfer.py', 'model/primitive_observer/runtime.py',
    'model/common_source/open_protocol.py', 'model/common_source/open_environment.py',
    'model/observer/native.py', 'model/common_source/slow_geometry.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/requirements.txt',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/SourceLaw/ReflectionPositivity.lean',
)
TOLERANCE = 2e-11
GAPS = (0, 1, 4)


@dataclass(frozen=True)
class Kernel:
    states: tuple
    seams: tuple
    counts: tuple
    marked_counts: tuple


def native_counts():
    """Count all native outcomes; no transfer matrix or field action is input."""
    states = native.states(binary=True)
    seams = tuple(slow_geometry.source().seams())
    if len(states) != 66 or len(seams) != 30:
        raise ArithmeticError('the declared binary carrier or scheduler changed')
    index = {state: i for i, state in enumerate(states)}
    counts = [[0]*66 for _ in states]
    marked = [[[0]*66 for _ in states] for _ in (0, 1)]
    for i, state in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                output, record = native.native(state, seam, coin)
                if (record[0] != coin or type(record[1]) is not int or record[1] not in (-1, 0, 1)
                        or output not in index
                        or output != tuple(state[k]-(int(k == seam[0])-int(k == seam[1]))*record[1]
                                           for k in range(12))):
                    raise ArithmeticError('native current, conservation or binary support changed')
                j = index[output]
                bit = int(abs(record[1]) == 1)
                counts[i][j] += 1
                marked[bit][i][j] += 1
    return Kernel(states, seams, tuple(map(tuple, counts)),
                  tuple(tuple(map(tuple, matrix)) for matrix in marked))


def exact_audit(kernel):
    """Check the finite hypotheses of positivity and the Markov reconstruction."""
    n, counts, marked = len(kernel.states), kernel.counts, kernel.marked_counts
    if (n != 66 or len(kernel.seams) != 30 or len(counts) != n or len(marked) != 2
            or any(len(matrix) != n or any(len(row) != n for row in matrix) for matrix in (counts, *marked))
            or any(type(value) is not int or value < 0 for matrix in (counts, *marked)
                   for row in matrix for value in row)):
        raise ValueError('complete nonnegative integer native transition counts required')
    if (any(sum(row) != 60 for row in counts)
            or any(counts[i][j] != counts[j][i] for i in range(n) for j in range(n))
            or any(marked[0][i][j]+marked[1][i][j] != counts[i][j]
                   for i in range(n) for j in range(n))):
        raise ArithmeticError('the complete native kernel must be symmetric and stochastic')
    for i in range(n):
        for j in range(n):
            if (marked[0][i][j] != (counts[i][i] if i == j else 0)
                    or marked[1][i][j] != (0 if i == j else counts[i][j])
                    or (i != j and counts[i][j] not in (0, 1))):
                raise ArithmeticError('native current bits do not label wait versus token exchange')
    degrees = tuple(60-counts[i][i] for i in range(n))
    if Counter(degrees) != {8: 30, 10: 36}:
        raise ArithmeticError('the actual two-token exchange degrees changed')
    reached, pending, tree = {0}, [0], []
    while pending:
        i = pending.pop(0)
        for j in range(n):
            if i != j and counts[i][j] and j not in reached:
                reached.add(j)
                pending.append(j)
                tree.append((i, j))
    if len(reached) != n:
        raise ArithmeticError('native transition support is disconnected')
    lower = F(min(counts[i][i]-sum(counts[i][j] for j in range(n) if j != i)
                  for i in range(n)), 60)
    if lower != F(2, 3):
        raise ArithmeticError('the exact positive transfer lower bound changed')
    return {'state_dimension': n, 'primitive_branches_checked': n*60,
            'transition_denominator': 60,
            'integer_counts': [list(row) for row in counts],
            'symmetric': True, 'row_and_column_sum': '1',
            'stationary_probability_per_configuration': '1/66',
            'connected_spanning_tree': [list(edge) for edge in tree],
            'exchange_degree_histogram': {'8': 30, '10': 36},
            'minimum_wait_probability': '5/6',
            'transfer_identity': 'K=I-L_token/60',
            'positive_operator_lower_bound': 'K >= (2/3)*I',
            'upper_operator_bound': 'K <= I',
            'positivity_proof': 'symmetric diagonal-dominance bound, exact integer counts',
            'constant_ground_direction_unique': True,
            'hidden_native_choices_are_probability_summed': True}


def _row_times(row, counts):
    result = [F(0)]*len(row)
    for i, weight in enumerate(row):
        if weight:
            for j, count in enumerate(counts[i]):
                if count:
                    result[j] += weight*F(count, 60)
    return tuple(result)


def _times_column(counts, vector):
    return tuple(sum((F(value, 60)*vector[j] for j, value in enumerate(row) if value), F(0))
                 for row in counts)


def _gap(value):
    if type(value) is not int or not 0 <= value <= 16:
        raise ValueError('an integer gap between zero and sixteen native attempts is required')
    return value


def exact_record_law(kernel, gap):
    """Full marked two-event law and each unnormalized terminal population."""
    gap = _gap(gap)
    prior = (F(1, 66),)*66
    result = {}
    for a in (0, 1):
        intermediate = _row_times(prior, kernel.marked_counts[a])
        for _ in range(gap):
            intermediate = _row_times(intermediate, kernel.counts)
        for b in (0, 1):
            terminal = _row_times(intermediate, kernel.marked_counts[b])
            result[a, b] = (sum(terminal), terminal)
    if sum(probability for probability, _ in result.values()) != 1:
        raise ArithmeticError('the complete current-bit law is not normalized')
    return result


def _bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or value < 0 or value > TOLERANCE:
        raise ArithmeticError(f'{name} exceeds the finite transfer replay tolerance')
    return value


def spectral_reconstruction(kernel):
    """Numerical evaluation of the unique self-adjoint log of positive K."""
    exact_audit(kernel)
    matrix = np.asarray(kernel.counts, dtype=float)/60
    values, vectors = np.linalg.eigh(matrix)
    if (not np.isfinite(values).all() or values[0] < 2/3-TOLERANCE
            or values[-1] > 1+TOLERANCE or values[-2] >= 1-TOLERANCE):
        raise ArithmeticError('the exact positive connected kernel has an inconsistent numerical spectrum')
    energies = -np.log(values)
    hamiltonian = (vectors*energies)@vectors.T
    vacuum = np.ones(66)/math.sqrt(66)
    diagnostics = {
        'minimum_transfer_eigenvalue': float(values[0]),
        'maximum_transfer_eigenvalue': float(values[-1]),
        'gap_per_Euclidean_native_attempt': float(-math.log(values[-2])),
        'largest_energy_per_Euclidean_native_attempt': float(energies[0]),
        'exp_minus_H_vs_native_K_residual': _bounded(np.linalg.norm(expm(-hamiltonian)-matrix, 'fro'),
                                                    'exponential reconstruction'),
        'H_vacuum_residual': _bounded(np.linalg.norm(hamiltonian@vacuum), 'ground direction'),
        'orthonormal_eigenbasis_residual': _bounded(np.linalg.norm(vectors.T@vectors-np.eye(66), 'fro'),
                                                   'orthonormal spectral basis'),
    }
    return hamiltonian, vacuum, diagnostics


def instruments(kernel):
    """Compile actual current-record C Kraus branches, retaining |j| only."""
    index = {state: i for i, state in enumerate(kernel.states)}
    groups = {}
    for seam_number, seam in enumerate(kernel.seams):
        for coin in (0, 1):
            for i, state in enumerate(kernel.states):
                output, record = native.native(state, seam, coin)
                key = (seam_number, *record)
                groups.setdefault(key, np.zeros((66, 66), complex))[index[output], i] = 1/math.sqrt(60)
    outcomes = tuple(tuple(groups[key] for key in sorted(groups) if abs(key[2]) == bit)
                     for bit in (0, 1))
    if tuple(map(len, outcomes)) != (60, 60):
        raise ArithmeticError('the complete native current Kraus dictionary changed')
    binding = {'native_rule': 'C', 'carrier': 'binary total-two, dimension 66',
               'scheduler': 'uniform thirty seams and fair placement coin',
               'record': 'bit = abs(native committed current)',
               'discarded_labels': ['seam', 'coin', 'current direction'],
               'native_channel_semantics': True, 'native_pulse_words_emitted': False}
    visible = runtime.Instrument(kraus=outcomes, visible=True,
        resources=(('native_attempts', 1), ('retained_current_bits', 1)), binding=binding)
    hidden = runtime.Instrument(kraus=(outcomes[0]+outcomes[1],),
        resources=(('native_attempts', 1),), binding={**binding, 'record': 'all labels forgotten'})
    return visible, hidden


def program(kernel, gap, compiled=None):
    gap = _gap(gap)
    visible, hidden = instruments(kernel) if compiled is None else compiled
    events = [runtime.Event('first_native_current_bit', instruments=(visible,))]
    if gap:
        events.append(runtime.Event('unread_native_attempts', sequence=(hidden,), repeat=gap))
    events.append(runtime.Event('second_native_current_bit', instruments=(visible,)))
    return runtime.Program(f'native_current_transfer_gap_{gap}', 66, np.eye(66)/66, tuple(events),
        binding={'preparation': 'stationary classical mixture I66/66',
                 'native_attempt_gap_between_recorded_events': gap,
                 'reconstruction_vacuum_is_not_this_prepared_density': True,
                 'complete_native_current_bits_retained': True})


def _reflection_check(kernel):
    """Exact cylinder factorization fixture for the general proof in the note."""
    # F reads two future current bits; its reflected past reads their reverse.
    values = {(0, 0): F(1), (0, 1): F(2), (1, 0): F(-1), (1, 1): F(3)}
    future = [F(0)]*66
    ones = (F(1),)*66
    for (a, b), value in values.items():
        conditional = _times_column(kernel.marked_counts[a],
                                    _times_column(kernel.marked_counts[b], ones))
        future = [x+value*y for x, y in zip(future, conditional, strict=True)]
    propagated = _times_column(kernel.counts, future)
    gram = sum((x*y/F(66) for x, y in zip(future, propagated, strict=True)), F(0))
    lower = F(2, 3)*sum((x*x/F(66) for x in future), F(0))
    direct = F(0)
    for (a, b), left_value in values.items():
        row = _row_times(_row_times((F(1, 66),)*66, kernel.marked_counts[b]), kernel.marked_counts[a])
        row = _row_times(row, kernel.counts)
        for (c, d), right_value in values.items():
            final = _row_times(_row_times(row, kernel.marked_counts[c]), kernel.marked_counts[d])
            direct += left_value*right_value*sum(final)
    if direct != gram or gram < lower or lower <= 0:
        raise ArithmeticError('the native marked-cylinder reflection factorization failed')
    return {'fixture': 'two current bits on each side of one unmarked reflection bond',
            'future_function_values_00_01_10_11': ['1', '2', '-1', '3'],
            'direct_native_marked_path_form': str(direct),
            'conditional_future_kernel_form': str(gram),
            'exact_positive_lower_bound_for_fixture': str(lower),
            'all_cylinder_argument': 'Markov conditional independence and reversibility reduce the form to f* K f/66',
            'non_Gaussian_stationary_history': True,
            'full_cylinder_reflection_proof_is_analytic': True,
            'new_Lean_proof_claimed': False}


def run():
    kernel = native_counts()
    exact = exact_audit(kernel)
    hamiltonian, vacuum, spectrum = spectral_reconstruction(kernel)
    compiled = instruments(kernel)
    rows = []
    for gap in GAPS:
        source_program = program(kernel, gap, compiled)
        history = runtime.execute(source_program)
        expected = exact_record_law(kernel, gap)
        if history.records != ((0, 0), (0, 1), (1, 0), (1, 1)):
            raise ArithmeticError('the full two-current record dictionary changed')
        transfer = expm(-gap*hamiltonian)
        marked = tuple(np.asarray(matrix, dtype=float)/60 for matrix in kernel.marked_counts)
        probability_error, output_error, spectral_error = 0., 0., 0.
        for branch in history.branches:
            probability, terminal = expected[branch.record]
            probability_error = max(probability_error, abs(branch.probability-float(probability)))
            output_error += float(np.linalg.norm(branch.density-np.diag(list(map(float, terminal))), 'fro'))
            a, b = branch.record
            predicted = float(vacuum@marked[a]@transfer@marked[b]@vacuum)
            spectral_error = max(spectral_error, abs(predicted-float(probability)))
        if any(dict(branch.resources) != {'native_attempts': gap+2, 'retained_current_bits': 2}
               for branch in history.branches):
            raise ArithmeticError('native event or complete record accounting changed')
        both = expected[1, 1][0]
        first = sum(expected[1, b][0] for b in (0, 1))
        second = sum(expected[a, 1][0] for a in (0, 1))
        rows.append({'unread_native_attempts': gap, 'total_native_attempts': gap+2,
                     'complete_records': [[0, 0], [0, 1], [1, 0], [1, 1]],
                     'exact_probabilities': [str(expected[a, b][0]) for a in (0, 1) for b in (0, 1)],
                     'runtime_probabilities': [branch.probability for branch in history.branches],
                     'single_event_moving_probability': str(first),
                     'moving_bit_covariance': str(both-first*second),
                     'quantum_output_dimension_per_record': 66,
                     'expanded_tree_resources': dict(history.expanded_tree_resources),
                     'runtime_probability_residual_diagnostic': _bounded(probability_error, 'current law'),
                     'complete_terminal_density_residual_diagnostic': _bounded(output_error, 'terminal CQ output'),
                     'spectral_transfer_probability_residual_diagnostic': _bounded(spectral_error, 'spectral current law')})
    checkpoint_program = program(kernel, 1, compiled)
    prefix = runtime.execute(checkpoint_program, stop_after=1)
    restored = runtime.checkpoint_from_json(checkpoint_program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(checkpoint_program, checkpoint=restored)
    uninterrupted = runtime.execute(checkpoint_program)
    checkpoint_error = _bounded(runtime.complete_distance(continued, uninterrupted), 'checkpoint continuation')
    return {'schema': 'oph.primitive_observer.native_transfer.v1', 'all_pass': True,
            'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'exact_kernel_attachment': exact,
            'spectral_reconstruction': {'Hamiltonian': 'H=-log(K)',
                'vacuum': 'constant normalized vector on native configurations',
                'vacuum_selected_as': 'unique ground direction representing the stationary Markov history',
                'self_adjoint_log_unique': True, 'parameter': 'Euclidean native attempt spacing set to one',
                'finite_dimension': 66, 'diagnostics': spectrum},
            'observer_current_records': rows,
            'reflection_positivity': _reflection_check(kernel),
            'checkpoint': {'restored_after_recorded_events': 1,
                           'complete_quantum_state_retained': True,
                           'serialized_state_independently_replayed': True,
                           'complete_cq_distance_diagnostic': checkpoint_error},
            'scope': {'native_architecture_and_uniform_scheduler_supplied': True,
                      'all_hidden_seam_coin_and_direction_branches_summed': True,
                      'retained_bit_comes_from_existing_native_current': True,
                      'retained_bit_is_coarse_absolute_current': True,
                      'complete_raw_seam_coin_signed_current_dictionary_retained': False,
                      'configuration_access_needed_for_reported_bit_law': False,
                      'initial_stationary_mixture_supplied_for_runtime_experiment': True,
                      'transfer_ground_vector_is_actual_native_preparation': False,
                      'marked_T0_plus_T1': 'K, not the identity',
                      'marked_insertions_are_claimed_POVM_effects': False,
                      'new_oscillator_action_or_coefficients_received': False,
                      'native_pulse_words_emitted': False,
                      'real_time_unitary_identified_with_native_or_physical_dynamics': False,
                      'classical_transfer_selects_complete_quantum_instrument': False,
                      'physical_clock_local_field_limit_or_empirical_prediction': False,
                      'roundoff_or_physical_service_certificate': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if result.get('all_pass') is not True:
        raise ArithmeticError('native transfer construction failed')
    if args.check:
        expected = strict_json.loads(args.check.read_text(encoding='utf-8'))
        compare_receipt(expected, result)
        print('Native transfer receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
