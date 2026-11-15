"""Exact single-carrier passage laws, with physical acceptance left explicit.

The positive result uses a disclosed total-one calibration preparation of the
unchanged native primitive. It is not the 66-state total-two experiment. A
second causal witness executes directly in that existing total-two sector.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as sp

from ..observer import native
from ..observer_reconstruction import native_geometry as ng
from ..relative_observer_geometry import intrinsic_geometry as intrinsic
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
PORTS = 12
SOURCE_PATHS = tuple(sorted(set((*intrinsic.SOURCE_PATHS,
    'model/physical_space_time/construction.py'))))


def _port(value):
    if type(value) is not int or value not in range(PORTS):
        raise ValueError('a native integer port from zero through eleven required')
    return value


def _matrix_rows(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


@lru_cache(maxsize=1)
def recovered_carrier():
    """Use the established current-record calibration, without coordinates."""
    records = ng.execute_calibration()
    result = intrinsic.recover_structure(records)
    if not result['rank_three_structural_gram_certified']:
        raise ArithmeticError('the record-derived carrier certificate failed')
    seams = tuple(tuple(edge) for edge in result['seams'])
    laplacian = sp.ImmutableMatrix([[sp.Rational(x) for x in row]
                                   for row in result['laplacian_exact']])
    return seams, laplacian


def kernel(receiver, *, holding=F(1), seam_weights=None):
    """Enumerate actual native branches; split remote versus incident attempts.

    Matrices are column-stochastic. Holding slots do not count as an incident
    repair attempt. Positive seam weights and fair independent coins are
    declared source services, never inferred from a seeded trajectory.
    """
    receiver = _port(receiver)
    if type(holding) is not F or not 0 < holding <= 1:
        raise ValueError('an exact positive holding acceptance at most one required')
    seams, laplacian = recovered_carrier()
    if seam_weights is None:
        weights = {edge: F(1, len(seams)) for edge in seams}
    elif (type(seam_weights) is dict and set(seam_weights) == set(seams)
          and all(type(w) is F and w > 0 for w in seam_weights.values())
          and sum(seam_weights.values(), F(0)) == 1):
        weights = seam_weights
    else:
        raise ValueError('positive exact normalized weights on every native seam required')
    remote, incident = (1-sp.Rational(holding))*sp.eye(PORTS), sp.zeros(PORTS)
    checks = 0
    for sender in range(PORTS):
        before = tuple(int(i == sender) for i in range(PORTS))
        for seam in seams:
            for coin in (0, 1):
                after, (_, current) = native.native(before, seam, coin)
                if sum(after) != 1 or any(value not in (0, 1) for value in after):
                    raise ArithmeticError('the unchanged primitive lost its total-one invariant')
                # The decoder sees only the disclosed initial receiver bit and
                # its incident signed current; the hidden output is an audit.
                decoded = before[receiver]+current*(int(receiver == seam[1])-int(receiver == seam[0]))
                if decoded != after[receiver]:
                    raise ArithmeticError('receiver occupation is not recovered by local currents')
                target = after.index(1)
                matrix = incident if receiver in seam else remote
                matrix[target, sender] += sp.Rational(holding*weights[seam]/2)
                checks += 1
    attempt = remote+incident
    p = sum((sp.Rational(holding*weights[e]) for e in seams if receiver in e), sp.Rational(0))
    ones = sp.ones(1, PORTS)
    if ones*attempt != ones or ones*incident != p*ones or ones*remote != (1-p)*ones:
        raise ArithmeticError('state-independent incident exposure failed')
    completed = incident*(sp.eye(PORTS)-remote).inv()
    if ones*completed != ones or any(value < 0 for value in completed):
        raise ArithmeticError('the completed incident kernel lost probability mass')
    if seam_weights is None and attempt != sp.eye(PORTS)-sp.Rational(holding)*laplacian/60:
        raise ArithmeticError('actual passage dynamics differ from current-reconstructed response')
    return {'attempt': sp.ImmutableMatrix(attempt), 'remote': sp.ImmutableMatrix(remote),
            'incident': sp.ImmutableMatrix(incident), 'completed': sp.ImmutableMatrix(completed),
            'tick_probability': p, 'native_branch_checks': checks}


def mean_arrival(matrix, receiver):
    """Exact Dirichlet problem, h(receiver)=0 and h=1+P^T h elsewhere."""
    receiver = _port(receiver)
    if not isinstance(matrix, sp.MatrixBase) or matrix.shape != (PORTS, PORTS):
        raise ValueError('a twelve-port transition matrix required')
    if (any(value.is_Rational is not True for value in matrix)
            or any(value < 0 for value in matrix)
            or sp.ones(1, PORTS)*matrix != sp.ones(1, PORTS)):
        raise ValueError('a nonnegative column-stochastic transition required')
    keep = [i for i in range(PORTS) if i != receiver]
    operator = sp.eye(PORTS-1)-matrix.extract(keep, keep).T
    if not operator.det():
        raise ValueError('the target is not reached with finite mean from every port')
    reduced = operator.inv()*sp.ones(PORTS-1, 1)
    result = sp.zeros(PORTS, 1)
    for i, value in zip(keep, reduced, strict=True):
        result[i] = value
    if any(result[i] <= 0 for i in keep):
        raise ArithmeticError('a first-arrival mean must be positive off its target')
    if any(result[i] != 1+(matrix.T*result)[i] for i in keep):
        raise ArithmeticError('the exact first-step equations failed')
    return sp.ImmutableMatrix(result)


def arrival_law(matrix, sender, receiver, cap):
    """All arrival bins through cap plus the remaining, unconditioned timeout."""
    sender, receiver = _port(sender), _port(receiver)
    if sender == receiver or type(cap) is not int or not 1 <= cap <= 64:
        raise ValueError('distinct ports and an integer cap from one through 64 required')
    mean_arrival(matrix, receiver)  # Authenticate the admitted finite stochastic domain.
    vector = sp.zeros(PORTS, 1)
    vector[sender] = 1
    probabilities = []
    for _ in range(cap):
        vector = sp.Matrix(matrix*vector)
        probabilities.append(vector[receiver])
        vector[receiver] = 0
    timeout = sum(vector)
    if sum(probabilities)+timeout != 1:
        raise ArithmeticError('arrival records and timeout must exhaust probability one')
    return {'first_arrival_probabilities': [str(x) for x in probabilities],
            'timeout_probability': str(timeout), 'cap': cap,
            'postselection_or_timeout_renormalization': False}


def graph_distances():
    seams, _ = recovered_carrier()
    neighbors = {i: set() for i in range(PORTS)}
    for a, b in seams:
        neighbors[a].add(b)
        neighbors[b].add(a)
    distances = []
    for start in range(PORTS):
        row = [None]*PORTS
        row[start] = 0
        pending = deque([start])
        while pending:
            point = pending.popleft()
            for neighbor in sorted(neighbors[point]):
                if row[neighbor] is None:
                    row[neighbor] = row[point]+1
                    pending.append(neighbor)
        distances.append(row)
    if any(value is None for row in distances for value in row):
        raise ArithmeticError('the recovered carrier must be connected')
    return distances


def native_witness(initial, receiver, schedule, coins):
    """Execute a fixed word and recover receiver state from its local tape."""
    receiver = _port(receiver)
    if (type(initial) is not tuple or len(initial) != PORTS
            or any(type(v) is not int or v not in (0, 1) for v in initial)
            or sum(initial) not in (1, 2)):
        raise ValueError('a disclosed total-one or total-two binary preparation required')
    if len(schedule) != len(coins):
        raise ValueError('one coin for each declared native seam required')
    seams, _ = recovered_carrier()
    state, field, records, trajectory = initial, initial[receiver], [], [initial[receiver]]
    for seam, coin in zip(schedule, coins, strict=True):
        if seam not in seams or type(coin) is not int or coin not in (0, 1):
            raise ValueError('a native oriented seam and binary coin required')
        state, (_, current) = native.native(state, seam, coin)
        if receiver in seam:
            field += current*(int(receiver == seam[1])-int(receiver == seam[0]))
            records.append({'seam': list(seam), 'current': current, 'decoded_field': field})
        if field != state[receiver]:
            raise ArithmeticError('the incident-current decoder failed')
        trajectory.append(field)
    return {'initial_occupied_ports': [i for i, value in enumerate(initial) if value],
            'schedule': [list(edge) for edge in schedule], 'coins': list(coins),
            'receiver': receiver, 'receiver_field_by_native_attempt': trajectory,
            'accessible_incident_records': records, 'incident_tick_count': len(records),
            'native_attempt_count': len(schedule)}


def run():
    seams, laplacian = recovered_carrier()
    distance = graph_distances()
    attempt_means, incident_means = sp.zeros(PORTS), sp.zeros(PORTS)
    kernels = []
    for receiver in range(PORTS):
        data = kernel(receiver)
        h = mean_arrival(data['attempt'], receiver)
        g = mean_arrival(data['completed'], receiver)
        if data['tick_probability'] != sp.Rational(1, 6) or g != h/6:
            raise ArithmeticError('the stopped incident clock and attempt means disagree')
        arrival_by_shell = {0: sp.Rational(1, 2), 1: sp.Rational(19, 222),
                            2: sp.Rational(1, 74), 3: sp.Rational(1, 222)}
        if any(data['completed'][receiver, sender] != arrival_by_shell[distance[receiver][sender]]
               for sender in range(PORTS)):
            raise ArithmeticError('first-tick arrival does not respect the declared shells')
        attempt_means[:, receiver], incident_means[:, receiver] = h, g
        kernels.append(data)
    if incident_means != incident_means.T:
        raise ArithmeticError('the uniform carrier passage metric must be symmetric')
    if any(incident_means[i, j] > incident_means[i, k]+incident_means[k, j]
           for i in range(PORTS) for j in range(PORTS) for k in range(PORTS)):
        raise ArithmeticError('the passage metric failed a triangle inequality')
    projector = sp.eye(PORTS)-sp.ones(PORTS)/PORTS
    centered = -projector*incident_means*projector/2
    inverse_response = (laplacian+sp.ones(PORTS)/PORTS).inv()-sp.ones(PORTS)/PORTS
    if centered != 60*inverse_response or centered.rank() != 11:
        raise ArithmeticError('passage geometry is not the exact all-mode response metric')
    receiver, farthest = 0, next(i for i in range(PORTS) if distance[0][i] == 3)
    clock = kernels[receiver]
    power = sp.eye(PORTS)
    for attempts in range(4):
        if any((power[i, j] > 0) != (distance[i][j] <= attempts)
               for i in range(PORTS) for j in range(PORTS)):
            raise ArithmeticError('native attempt response lost exact shortest-path support')
        power = clock['attempt']*power
    minimum_arrival = min(clock['completed'][0, j] for j in range(1, PORTS))
    if minimum_arrival != sp.Rational(1, 222):
        raise ArithmeticError('the uniform first-arrival tail bound changed')
    laws = []
    for shell in (1, 2, 3):
        sender = next(i for i in range(PORTS) if distance[receiver][i] == shell)
        first = clock['completed'][receiver, sender]
        if first <= 0:
            raise ArithmeticError('remote waiting must permit a first-tick arrival in each shell')
        laws.append({'sender': sender, 'receiver': receiver, 'graph_distance': shell,
            'mean_native_attempts': str(attempt_means[sender, receiver]),
            'mean_receiver_incident_ticks': str(incident_means[sender, receiver]),
            'arrival_law_in_receiver_incident_ticks': arrival_law(clock['completed'], sender, receiver, 8),
            'arrival_law_in_native_attempts': arrival_law(clock['attempt'], sender, receiver, 8)})
    held = kernel(0, holding=F(2, 7))
    if held['completed'] != clock['completed'] or mean_arrival(held['attempt'], 0) != sp.Rational(7, 2)*attempt_means[:, 0]:
        raise ArithmeticError('independent identity holding changed the incident experiment')
    # The same causal obstruction also occurs in the established 66-state
    # sector. The initial difference is supported at ports two and three.
    words = tuple(native_witness(tuple(int(i in occupied) for i in range(PORTS)), 0,
                                ((2, 10), (0, 10)), (0, 1))
                  for occupied in ((1, 2), (1, 3)))
    if words[0]['receiver_field_by_native_attempt'] != [0, 0, 1] or words[1]['receiver_field_by_native_attempt'] != [0, 0, 0]:
        raise ArithmeticError('the two-particle local-clock causal witness changed')
    if min(distance[0][i] for i in (2, 3)) != 2 or any(w['incident_tick_count'] != 1 for w in words):
        raise ArithmeticError('the clock counterexample lost its support or tick count')
    # A differing operational metric is not a contradiction: its readout
    # dictionary differs. It prevents identifying physical dimension by rank alone.
    shells = Counter((distance[i][j], str(attempt_means[i, j]), str(incident_means[i, j]))
                     for i in range(PORTS) for j in range(i+1, PORTS))
    return {
        'schema': 'oph.physical_space_time.single_carrier.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'contract': {'carriers': 1, 'ports': PORTS, 'seams': len(seams),
            'passage_preparation': 'one known native load at sender, every other port empty',
            'calibration_sector_total_load': 1,
            'same_preparation_as_retained_binary_total_two_experiment': False,
            'scheduler': 'independent uniform native seams and fair placement coins',
            'receiver_clock': 'every seam attempt incident to receiver, including zero current',
            'arrival_record': 'first incident current that changes disclosed receiver load from zero to one',
            'initialization_reset_and_faithful_local_record_services_are_supplied': True,
            'remote_attempt_counter_used_by_receiver': False,
            'source_coordinates_used_by_estimator': False,
            'federation_required': False},
        'native_branch_readback_checks': sum(data['native_branch_checks'] for data in kernels),
        'attempt_kernel': _matrix_rows(clock['attempt']),
        'receiver_zero_completed_incident_kernel': _matrix_rows(clock['completed']),
        'mean_arrival_native_attempts': _matrix_rows(attempt_means),
        'mean_arrival_receiver_incident_ticks': _matrix_rows(incident_means),
        'passage_shells': [{'graph_distance': d, 'mean_native_attempts': h,
                           'mean_receiver_incident_ticks': g, 'unordered_pairs': count}
                          for (d, h, g), count in sorted(shells.items())],
        'completed_laws_and_retained_timeouts': laws,
        'finite_completion_bound': {
            'minimum_conditional_arrival_probability_next_incident_tick': str(minimum_arrival),
            'tail_after_n_receiver_incident_ticks_upper': '(221/222)^n',
            'native_attempt_tail_after_3n_attempts_upper': '(215999/216000)^n',
            'timeout_records_included_in_each_capped_law': True},
        'derived_equalities': {
            'native_attempt_kernel': 'P=I-L/60',
            'completed_receiver_kernel': 'K_v=B_v(I-A_v)^(-1)',
            'all_132_distinct_port_means_satisfy_first_step_equations': True,
            'mean_receiver_ticks_equal_mean_native_attempts_over_six': True,
            'incident_passage_metric_symmetric_and_all_triangles_checked': True,
            'squared_euclidean_passage_distance_matrix': 'H_incident=60 R_effective',
            'centered_passage_gram': '60 L^+',
            'centered_passage_gram_rank': int(centered.rank()),
            'existing_slow_band_gram_rank': 3,
            'rank_three_selected_as_unique_physical_dimension': False},
        'clock_controls': {
            'independent_identity_holding_acceptance': '2/7',
            'completed_incident_kernel_unchanged_exact': True,
            'native_attempt_means_multiplied_by': '7/2',
            'movement_only_clock_first_arrival_tick_for_every_distinct_pair': 1,
            'movement_only_clock_can_distinguish_passage_distances': False},
        'causal_boundary': {
            'first_tick_arrival_at_graph_distance_three_probability': str(clock['completed'][0, farthest]),
            'native_attempt_propagation_respects_graph_path_support': True,
            'one_graph_edge_per_receiver_incident_tick_bound': False,
            'total_two_witnesses': list(words),
            'witness_initial_difference_support': [2, 3],
            'witness_support_graph_distance_to_receiver': 2,
            'disproves_every_possible_physical_clock_or_spacetime_construction': False},
        'acceptance': {
            'finite_operational_signal_and_local_clock_laws_established': True,
            'finite_operational_passage_metric_established': True,
            'physical_three_plus_one_dimensional_spacetime_established': False,
            'lorentzian_signal_law_or_spatial_continuum_established': False,
            'lab_clock_identification_or_empirical_comparison_established': False,
            'specific_federation_architecture_needed_for_these_results': False,
            'compatible_extended_observable_domain_needed_for_continuum_claim': True},
        'all_pass': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Single-carrier passage geometry and local-clock receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
