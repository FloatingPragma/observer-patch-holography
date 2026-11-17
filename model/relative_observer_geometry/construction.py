"""Reconstruct local fields from complete incident native-current records.

Known preparation, incident-record access and common event identifiers are
declared. Event identifiers give source order, never a physical duration.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from pathlib import Path

import numpy as np

from ..observer import native
from ..observer_reconstruction import native_geometry
from ..federated_observer import refinement as observer_refinement
from ..native_action import screen_compatibility
from ..primitive_observer import algebra, runtime, transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
LETTERS = tuple(product((0, 1), (-1, 0, 1)))
TAPE = ((0, 1), (0, 1), (0, 5), (0, 5), (0, 1))
COVER = (tuple(range(7)), tuple(range(5, 12)))
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/relative_observer_geometry/__init__.py',
    'model/relative_observer_geometry/__main__.py',
    'model/relative_observer_geometry/construction.py',
    'model/observer_reconstruction/native_geometry.py',
    *observer_refinement.SOURCE_PATHS,
    *screen_compatibility.SOURCE_PATHS,
    'model/primitive_observer/algebra.py', 'model/common_source/operational.py',
    *transfer.SOURCE_PATHS,
    'Lean/Time/ClockComparison.lean', 'Lean/QFT/FiniteCausalObserverNet.lean',
    'Lean/QFT/ObserverNetDescent.lean', 'Lean/QFT/SourceHistoryGNSDynamics.lean',
    'Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean',
)))


@lru_cache(maxsize=1)
def source():
    kernel = transfer.native_counts()
    transfer.exact_audit(kernel)
    return kernel


def _ports(ports):
    if (type(ports) is not tuple or not ports
            or any(type(i) is not int or i not in range(12) for i in ports)
            or len(set(ports)) != len(ports)):
        raise ValueError('distinct nonempty native port indices required')
    return tuple(sorted(ports))


def _record(record):
    if type(record) is not tuple or len(record) != 5:
        raise ValueError('record must be (event id, first port, second port, coin, current)')
    event, a, b, coin, current = record
    if (any(type(i) is not int for i in record) or event < 1
            or (a, b) not in source().seams
            or (coin, current) not in ((0, 0), (0, 1), (1, -1), (1, 0))):
        raise ValueError('an ordered native seam record with positive event id required')
    return record


@dataclass(frozen=True)
class LocalCounter:
    """Only disclosed local preparation and incident public records are stored.

    Completeness/authenticity of the incident stream is an access premise;
    a local observer cannot detect an omitted external boundary current.
    """
    ports: tuple[int, ...]
    initial: tuple[int, ...]
    records: tuple[tuple[int, ...], ...] = ()

    def __post_init__(self):
        if _ports(self.ports) != self.ports:
            raise ValueError('local ports must be in canonical increasing order')
        if (type(self.initial) is not tuple or len(self.initial) != len(self.ports)
                or any(type(v) is not int or v not in (0, 1) for v in self.initial)
                or sum(self.initial) > 2 or type(self.records) is not tuple):
            raise ValueError('a disclosed binary local preparation is required')
        previous = 0
        values = list(self.initial)
        for record in self.records:
            event, a, b, coin, current = _record(record)
            if event <= previous or not (a in self.ports or b in self.ports):
                raise ValueError('distinct increasing incident event records required')
            previous = event
            if a in self.ports and b in self.ports:
                difference = values[self.ports.index(a)]-values[self.ports.index(b)]
                expected = (difference-(2*coin-1)*abs(difference))//2
                if current != expected:
                    raise ValueError('committed current contradicts known endpoint loads and coin')
            for j, port in enumerate(self.ports):
                values[j] += current*(int(port == b)-int(port == a))
            if any(v not in (0, 1) for v in values) or sum(values) > 2:
                raise ValueError('current record contradicts the disclosed binary preparation')

    @property
    def values(self):
        values = list(self.initial)
        for _, a, b, _, current in self.records:
            for j, port in enumerate(self.ports):
                values[j] += current*(int(port == b)-int(port == a))
        return tuple(values)

    @property
    def ticks(self):
        return sum(current != 0 for _, _, _, _, current in self.records)

    def continue_with(self, records):
        if type(records) is not tuple:
            raise ValueError('immutable incident record continuation required')
        return LocalCounter(self.ports, self.initial, self.records+records)

    def checkpoint(self):
        return {'schema': 'oph.local_current_counter.v1', 'ports': list(self.ports),
                'initial': list(self.initial), 'records': [list(r) for r in self.records]}

    @staticmethod
    def restore(value):
        if (type(value) is not dict or set(value) != {'schema', 'ports', 'initial', 'records'}
                or value['schema'] != 'oph.local_current_counter.v1'
                or any(type(value[k]) is not list for k in ('ports', 'initial', 'records'))
                or any(type(r) is not list for r in value['records'])):
            raise ValueError('complete local counter checkpoint required')
        return LocalCounter(tuple(value['ports']), tuple(value['initial']),
                            tuple(tuple(r) for r in value['records']))


def local(initial, records, ports):
    """Apparatus hands each observer only its declared local inputs."""
    ports = _ports(ports)
    if (type(initial) is not tuple or len(initial) != 12
            or any(type(v) is not int or v not in (0, 1) for v in initial)
            or sum(initial) != 2 or type(records) is not tuple):
        raise ValueError('binary total-two preparation and immutable record history required')
    for record in records:
        _record(record)
    if any(a[0] >= b[0] for a, b in zip(records, records[1:])):
        raise ValueError('public event ids must be strictly increasing')
    return LocalCounter(ports, tuple(initial[i] for i in ports),
                        tuple(r for r in records if r[1] in ports or r[2] in ports))


def glue(views, ports=tuple(range(12))):
    """Glue values at one declared common history cut, not marginal laws.

    Simultaneity here means a shared retained prefix; cut agreement and
    source authenticity must accompany the supplied views.
    """
    ports = _ports(ports)
    if type(views) is not tuple or not views or any(type(v) is not LocalCounter for v in views):
        raise ValueError('nonempty tuple of local counters required')
    joined = {}
    for view in views:
        for port, value in zip(view.ports, view.values, strict=True):
            if port in joined and joined[port] != value:
                raise ValueError('local sections disagree on their overlap')
            joined[port] = value
    if not set(ports) <= set(joined):
        raise ValueError('declared regions do not cover the requested ports')
    return tuple(joined[i] for i in ports)


def execute(initial, tape, coins):
    """Exact diagonal native execution; hidden states are audit output only."""
    if (type(initial) is not tuple or len(initial) != 12
            or any(type(v) is not int for v in initial) or initial not in source().states
            or type(tape) is not tuple or type(coins) is not tuple):
        raise ValueError('source preparation, immutable seam tape and coins required')
    if len(tape) != len(coins) or len(tape) > 128:
        raise ValueError('equal bounded tape and coin lengths required')
    state, records = initial, []
    for event, (seam, coin) in enumerate(zip(tape, coins, strict=True), 1):
        if (type(seam) is not tuple or len(seam) != 2
                or any(type(i) is not int for i in seam) or seam not in source().seams
                or type(coin) is not int or coin not in (0, 1)):
            raise ValueError('admitted source seam and native coin required')
        state, (_, current) = native.native(state, seam, coin)
        records.append((event, *seam, coin, current))
    return state, tuple(records)


def affine_comparison(shared):
    """Two shared-event anchors calibrate; all further pairs independently test.

    The common event identity is a premise. Readings are exact rationals,
    strictly increasing on this comparison domain, not physical seconds.
    """
    if type(shared) is not tuple or len(shared) < 2:
        raise ValueError('at least two shared events required')
    if any(type(row) is not tuple or len(row) != 3 for row in shared):
        raise ValueError('shared event and two exact readings required')
    if (any(type(e) is not int or e < 1 or type(a) not in (int, F)
            or type(b) not in (int, F) for e, a, b in shared)
            or any(x[0] >= y[0] or x[1] >= y[1] or x[2] >= y[2]
                   for x, y in zip(shared, shared[1:]))):
        raise ValueError('ordered shared events and strictly increasing exact clocks required')
    _, a, b = shared[0]; _, c, d = shared[1]
    slope = F(d-b)/F(c-a)
    offset = F(b)-slope*a
    residuals = tuple(F(y)-slope*x-offset for _, x, y in shared[2:])
    return slope, offset, residuals


def branch_audit():
    checks = 0
    for initial in source().states:
        for seam in source().seams:
            for coin in (0, 1):
                target, records = execute(initial, (seam,), (coin,))
                views = tuple(local(initial, records, region) for region in COVER)
                if glue(views) != target:
                    raise ArithmeticError('incident current reconstruction failed')
                checks += 1
    return {'native_branches_checked_exactly': checks,
            'identity': 'n_i(k)=n_i(0)+sum_{e<=k} j_e (1[i=b_e]-1[i=a_e])',
            'extension_to_every_finite_prefix': 'sum the exact one-step identities',
            'boundary_currents_retained': True}


def seam_instrument(seam):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(i) is not int for i in seam) or seam not in source().seams):
        raise ValueError('admitted source seam required')
    matrices = [np.zeros((66, 66), complex) for _ in LETTERS]
    index = {state: i for i, state in enumerate(source().states)}
    for i, state in enumerate(source().states):
        for coin in (0, 1):
            target, letter = native.native(state, seam, coin)
            matrices[LETTERS.index(letter)][index[target], i] = 1/np.sqrt(2)
    return runtime.Instrument(kraus=tuple((m,) for m in matrices), visible=True,
        resources=(('native_attempts', 1), ('signed_current_reads', 1)),
        binding={'seam': list(seam), 'outcomes': [list(x) for x in LETTERS],
                 'scheduler': 'declared fixed seam tape', 'coin': 'fair'})


def runtime_check():
    initial = tuple(int(i in (0, 11)) for i in range(12))
    index = {state: i for i, state in enumerate(source().states)}
    density = np.zeros((66, 66)); density[index[initial], index[initial]] = 1
    program = runtime.Program('relative_current_observers', 66, density,
        tuple(runtime.Event(f'read_{i+1}', instruments=(seam_instrument(seam),))
              for i, seam in enumerate(TAPE[:3])),
        preparation_resources=(('disclosed_native_preparations', 1),),
        binding={'preparation': list(initial), 'complete_incident_record_access': True})
    full = runtime.execute(program)
    prefix = runtime.execute(program, stop_after=1)
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program, checkpoint=restored)
    exact = {}
    for coins in product((0, 1), repeat=3):
        target, records = execute(initial, TAPE[:3], coins)
        word = tuple(LETTERS.index((r[3], r[4])) for r in records)
        exact[word] = target, records
    error = 0.
    values = {}
    for branch in full.branches:
        expected = np.zeros((66, 66))
        values[branch.record] = 0
        if branch.record in exact:
            target, records = exact[branch.record]
            expected[index[target], index[target]] = 1/8
            views = tuple(local(initial, records, region) for region in COVER)
            if glue(views) != target:
                raise ArithmeticError('runtime record observer disagrees with terminal state')
            values[branch.record] = views[0].values[0]
        error = max(error, float(np.linalg.norm(branch.density-expected, 'fro')))
    view = algebra.observer_view(full)
    recorded = algebra.public_expectation(view, values).real
    hidden = algebra.expectation(view, np.diag([state[0] for state in source().states])).real
    distance = runtime.complete_distance(full, continued)
    if max(error, abs(recorded-hidden), distance) > 2e-11:
        raise ArithmeticError('complete runtime or public field readback changed')
    return {'retained_record_words': len(full.branches), 'positive_record_words': len(exact),
            'zero_probability_record_words': len(full.branches)-len(exact),
            'complete_density_max_frobenius_diagnostic': error,
            'checkpoint_complete_distance_diagnostic': distance,
            'public_n0_expectation': float(recorded), 'hidden_n0_audit_expectation': float(hidden),
            'native_attempts_per_trial': 3, 'initial_preparation_disclosed': list(initial),
            'roundoff_enclosure': False}


def history_check():
    """Complete finite diagonal ensemble; no path or failed read is discarded."""
    joint = Counter(); left_fibers = Counter(); right_fibers = Counter()
    retained = []
    for initial in source().states:
        for coins in product((0, 1), repeat=len(TAPE)):
            target, records = execute(initial, TAPE, coins)
            views = tuple(local(initial, records, region) for region in COVER)
            if glue(views) != target:
                raise ArithmeticError('complete history gluing failed')
            for view in views:
                split = len(view.records)//2
                prefix = LocalCounter(view.ports, view.initial, view.records[:split])
                continued = LocalCounter.restore(prefix.checkpoint()).continue_with(view.records[split:])
                if continued != view:
                    raise ArithmeticError('observer-local checkpoint failed')
            a = local(initial, records, (0,)); b = local(initial, records, (1,))
            joint[a.ticks, b.ticks, target[0], target[1]] += 1
            left_fibers[views[0].initial, views[0].records] += 1
            right_fibers[views[1].initial, views[1].records] += 1
            retained.append([list(initial), [list(r) for r in records], list(target)])
    denominator = len(retained)
    if denominator != 2112 or sum(joint.values()) != denominator:
        raise ArithmeticError('complete finite history mass changed')
    return {'declared_tape': [list(x) for x in TAPE], 'complete_histories': denominator,
            'history_probability': str(F(1, denominator)),
            'initial_uniform_binary_label_disclosed_before_execution': True,
            'overlapping_regions': [list(x) for x in COVER],
            'all_pathwise_overlaps_and_full_coverage_exact': True,
            'all_local_checkpoint_continuations_exact': True,
            'joint_records': [{'ticks_port0': a, 'ticks_port1': b, 'n0': x, 'n1': y,
                               'probability': str(F(count, denominator))}
                              for (a,b,x,y), count in sorted(joint.items())],
            'complete_history_sha256': hashlib.sha256(canonical(retained).encode()).hexdigest(),
            'local_record_fiber_counts': [len(left_fibers), len(right_fibers)],
            'GNS_pullback_identity': 'sum_h p(h)|f(view(h))|^2=sum_v p(v)|f(v)|^2',
            'pushforward_probability_mass_exact': sum(left_fibers.values()) == denominator
                and sum(right_fibers.values()) == denominator,
            'homomorphic_matrix_algebra_retractions_constructed': False}


def clock_checks():
    initial = tuple(int(i in (0, 11)) for i in range(12))
    _, records = execute(initial, TAPE, (0, 1, 0, 1, 0))
    readings = []
    for k, record in enumerate(records, 1):
        if record[1:3] == (0, 1) and record[4]:
            a = local(initial, records[:k], (0,)).ticks
            b = local(initial, records[:k], (1,)).ticks
            readings.append((record[0], a, b))
    slope, offset, residual = affine_comparison(tuple(readings))
    # Same shared seam counter, with a declared affine display convention.
    success = tuple((event, i, 2*i+5) for i, (event, _, _) in enumerate(readings, 1))
    good = affine_comparison(success)
    if residual != (F(-2),) or good != (F(2), F(5), (F(0),)):
        raise ArithmeticError('shared-event comparison controls changed')
    return {'record_history': [list(r) for r in records], 'history_probability_given_preparation': '1/32',
            'local_incident_tick_pairs': [list(r) for r in readings],
            'first_two_anchors_affine_slope': str(slope), 'offset': str(offset),
            'predesignated_third_event_residual': str(residual[0]),
            'affine_shared_seam_display_pairs': [list(r) for r in success],
            'affine_shared_seam_third_residual': str(good[2][0]),
            'parameter': 'integer local nonzero incident movements or shared-seam movements',
            'common_event_identity_and_record_access_supplied': True,
            'physical_time_or_Lorentz_frame_identification': False}


def controls():
    initial = tuple(int(i in (0, 11)) for i in range(12))
    target, record = execute(initial, ((0, 5),), (0,))
    complete = local(initial, record, (0,))
    omitted = local(initial, (), (0,))
    # Both priors generate the same no-current record on seam(0,1), coin1.
    second = tuple(int(i in (2, 11)) for i in range(12))
    _, first_records = execute(initial, ((0, 1),), (1,))
    _, second_records = execute(second, ((0, 1),), (1,))
    if first_records != second_records or initial[0] == second[0]:
        raise ArithmeticError('hidden preparation ambiguity disappeared')
    if complete.values == omitted.values or complete.values != (target[0],):
        raise ArithmeticError('boundary-current omission control disappeared')
    # Full regional marginal laws do not determine cross-region correlations.
    laws = (((0, 10), (1, 11)), ((0, 11), (1, 10)))
    regions = (tuple(range(6)), tuple(range(6, 12)))
    marginals = [[Counter(tuple(int(i in state) for i in region) for state in law)
                  for region in regions] for law in laws]
    if marginals[0] != marginals[1]:
        raise ArithmeticError('marginal-law descent control is invalid')
    return {'omitted_boundary_current': {'correct_n0': complete.values[0], 'incorrect_n0': omitted.values[0]},
            'same_record_hidden_preparations': {'record': [list(r) for r in first_records],
                'possible_initial_n0': [initial[0], second[0]]},
            'same_disjoint_regional_marginals_different_global_laws': {
                'configurations_law_A': [list(x) for x in laws[0]],
                'configurations_law_B': [list(x) for x in laws[1]],
                'probability_per_configuration': '1/2', 'E_n0_n10_A': '1/2', 'E_n0_n10_B': '0'},
            'local_state_or_record_marginals_imply_unique_joint_law': False}


def spatial_shells(gram_qsqrt5, seams):
    """Check the exact Euclidean shells encoded by the record-recovered Gram.

    Each entry [a,b] denotes a+b*sqrt(5). This is a finite carrier metric,
    with no length scale or continuum interpretation attached.
    """
    if (type(gram_qsqrt5) is not list or len(gram_qsqrt5) != 12
            or any(type(row) is not list or len(row) != 12 for row in gram_qsqrt5)
            or type(seams) not in (tuple, list)):
        raise ValueError('a complete twelve-port exact Gram and seam list are required')
    gram = []
    for row in gram_qsqrt5:
        converted = []
        for entry in row:
            if (type(entry) is not list or len(entry) != 2
                    or any(type(value) is not str for value in entry)):
                raise ValueError('Gram entries must be exact rational coefficient pairs')
            try:
                converted.append((F(entry[0]), F(entry[1])))
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError('Gram coefficients must be exact rationals') from error
        gram.append(converted)
    edge_set = set()
    for edge in seams:
        if (type(edge) not in (tuple, list) or len(edge) != 2
                or any(type(port) is not int or not 0 <= port < 12 for port in edge)
                or edge[0] == edge[1]):
            raise ValueError('distinct native seam endpoints are required')
        edge_set.add(tuple(sorted(edge)))
    if len(edge_set) != 30:
        raise ValueError('the complete thirty-seam carrier is required')
    if any(gram[i][j] != gram[j][i] for i in range(12) for j in range(12)):
        raise ValueError('the recovered Gram must be symmetric')
    if any(gram[i][i] != (F(1), F(0)) for i in range(12)):
        raise ValueError('the recovered Gram must have unit diagonal')
    adjacency = [set() for _ in range(12)]
    for a, b in edge_set:
        adjacency[a].add(b); adjacency[b].add(a)
    if any(len(neighbours) != 5 for neighbours in adjacency):
        raise ValueError('the recovered carrier graph must be five-regular')
    graph_distances = []
    for origin in range(12):
        distances = {origin: 0}
        frontier = [origin]
        while frontier:
            current = frontier.pop(0)
            for neighbour in adjacency[current]:
                if neighbour not in distances:
                    distances[neighbour] = distances[current]+1
                    frontier.append(neighbour)
        if len(distances) != 12 or max(distances.values()) != 3:
            raise ValueError('the source seam graph must have diameter three')
        graph_distances.append(distances)

    # Classify pairs directly from the exact Gram, then independently compare
    # its shortest shell with the executed native seam interface.
    shell_counts = Counter()
    squared_distance_by_shell = {}
    graph_shell = {1: ('edge', (F(0), F(1, 5))),
                   2: ('two_edge', (F(0), F(-1, 5))),
                   3: ('antipode', (F(-1), F(0)))}
    for i in range(12):
        for j in range(i+1, 12):
            dot = gram[i][j]
            squared = (F(2)-2*dot[0], -2*dot[1])
            distance = graph_distances[i][j]
            if distance not in graph_shell or dot != graph_shell[distance][1]:
                raise ValueError('the recovered Gram disagrees with its graph-distance shell')
            shell = graph_shell[distance][0]
            if (shell == 'edge') != ((i, j) in edge_set):
                raise ValueError('native seams must be exactly the shortest Gram shell')
            shell_counts[shell] += 1
            if shell in squared_distance_by_shell and squared_distance_by_shell[shell] != squared:
                raise ArithmeticError('a graph shell has inconsistent exact Euclidean distances')
            squared_distance_by_shell[shell] = squared
    expected_counts = {'edge': 30, 'two_edge': 30, 'antipode': 6}
    expected_distances = {'edge': [F(2), F(-2, 5)],
                          'two_edge': [F(2), F(2, 5)],
                          'antipode': [F(4), F(0)]}
    if dict(shell_counts) != expected_counts or squared_distance_by_shell != {
            key: tuple(value) for key, value in expected_distances.items()}:
        raise ArithmeticError('the exact spatial shell census changed')
    return {'exact_pair_counts': expected_counts,
            'squared_chord_lengths_qsqrt5': {
                'edge': '2-2*sqrt(5)/5', 'two_edge': '2+2*sqrt(5)/5',
                'antipode': '4'},
            'native_seams_are_exactly_shortest_shell': True,
            'metric_is_unitless_and_finite': True,
            'length_calibration_or_continuum_claimed': False}


def geometry_connection():
    recovery = native_geometry.recover(native_geometry.execute_calibration())
    laplacian = recovery['laplacian_exact']
    recovered = tuple((a, b) for a in range(12) for b in range(a+1, 12)
                      if F(laplacian[a][b]) == -1)
    if recovered != source().seams or not recovery['rank_three_unit_gram_recovered']:
        raise ArithmeticError('record-recovered graph does not match counter interface')
    shells = spatial_shells(recovery['gram_qsqrt5'], recovered)
    permutation = (7, 2, 10, 0, 5, 11, 8, 4, 1, 9, 3, 6)
    relabelled_gram = [[None]*12 for _ in range(12)]
    for i in range(12):
        for j in range(12):
            relabelled_gram[permutation[i]][permutation[j]] = recovery['gram_qsqrt5'][i][j]
    relabelled_seams = tuple(tuple(sorted((permutation[a], permutation[b]))) for a, b in recovered)
    if spatial_shells(relabelled_gram, relabelled_seams) != shells:
        raise ArithmeticError('the exact spatial metric shell census is not port-relabeling covariant')
    return {'response_record_count': recovery['record_count'],
            'response_recovered_seams_exactly_match_current_interfaces': True,
            'recovered_seams': [list(x) for x in recovered],
            'rank_of_certified_slow_response_Gram': 3,
            'exact_finite_spatial_metric': shells,
            'metric_shell_census_covariant_under_port_relabeling': True,
            'source_screen_refinement': refinement_metric_bridge(),
            'distinct_total_two_calibration_preparations_at_single_ports': 12,
            'calibration_sector_is_same_as_binary_history_preparation': False,
            'physical_spatial_dimension_inferred': False,
            'new_spatial_coordinates_or_Lorentz_cones_supplied': False}


def _ball_distances(ball):
    """Exact port-hop distances in a complete finite source ball."""
    adjacency = dict(ball.adjacency)
    distance = {port: 0 for port in range(12)}
    queue = list(range(12))
    for port in queue:
        for neighbour in adjacency[port]:
            if neighbour not in distance:
                distance[neighbour] = distance[port]+1
                queue.append(neighbour)
    return distance


def _ball_view(ball, distance, radius):
    """Restrict a complete ball to a rooted induced sub-ball."""
    ports = tuple(sorted(port for port, hops in distance.items() if hops <= radius))
    allowed = set(ports)
    rows = dict(ball.adjacency)
    adjacency = tuple((port, tuple(neighbour for neighbour in rows[port]
                                   if neighbour in allowed)) for port in ports)
    addresses = tuple(row for row in ball.addresses if row[0] in allowed)
    return ports, adjacency, addresses


def _full_screen_hop_profiles():
    """Exact port- and carrier-ball volumes on every retained source level."""
    profiles = []
    by_level = {}
    for level in range(6):
        graph = screen_compatibility.screen(level)
        root = 4**level-1
        distance = {port: 0 for port in range(12*root, 12*root+12)}
        queue = list(distance)
        for port in queue:
            for neighbour in graph.graph.neighbour_lists[port]:
                if neighbour not in distance:
                    distance[neighbour] = distance[port]+1
                    queue.append(neighbour)
        if len(distance) != graph.ports:
            raise ArithmeticError('native seam graph ball failed to cover its source federation')
        diameter = max(distance.values())
        port_spheres = [0]*(diameter+1)
        carrier_first_hop = {}
        for port, hops in distance.items():
            port_spheres[hops] += 1
            cell = port//12
            carrier_first_hop[cell] = min(hops, carrier_first_hop.get(cell, hops))
        carrier_spheres = [0]*(diameter+1)
        for hops in carrier_first_hop.values():
            carrier_spheres[hops] += 1
        port_balls, carrier_balls = [], []
        port_total = carrier_total = 0
        for port_count, carrier_count in zip(port_spheres, carrier_spheres, strict=True):
            port_total += port_count
            carrier_total += carrier_count
            port_balls.append(port_total)
            carrier_balls.append(carrier_total)
        row = {'level': level, 'root_carrier': root,
            'carriers': graph.carriers, 'ports': graph.ports,
            'native_seams': graph.seams,
            'root_port_eccentricity_in_native_seam_hops': diameter,
            'port_ball_counts_by_radius': port_balls,
            'carrier_ball_counts_by_radius': carrier_balls}
        profiles.append(row)
        by_level[level] = row
    doubling = []
    for level in range(3, 6):
        profile = by_level[level]
        for lower in (16, 32, 64, 128, 256):
            upper = 2*lower
            if upper < len(profile['port_ball_counts_by_radius']):
                port_counts = profile['port_ball_counts_by_radius']
                carrier_counts = profile['carrier_ball_counts_by_radius']
                doubling.append({'level': level, 'radius_pair': [lower, upper],
                    'port_ball_growth_ratio': str(F(port_counts[upper], port_counts[lower])),
                    'carrier_ball_growth_ratio': str(F(carrier_counts[upper], carrier_counts[lower]))})
    return {'distance': 'unweighted shortest-path count of native repair seams',
        'profiles': profiles, 'exact_doubling_ratios': doubling,
        'not_a_calibrated_physical_length_or_asymptotic_dimension': True}


def _rotation_between(source_vector, target_vector):
    """Return an orthogonal matrix sending one unit vector to another."""
    source_vector = np.asarray(source_vector, dtype=float)
    target_vector = np.asarray(target_vector, dtype=float)
    cosine = float(np.dot(source_vector, target_vector))
    cross = np.cross(source_vector, target_vector)
    skew = np.array([[0., -cross[2], cross[1]],
                     [cross[2], 0., -cross[0]],
                     [-cross[1], cross[0], 0.]])
    if cosine <= -1.+1e-12:
        raise ArithmeticError('the selected endpoint alignment is antipodal')
    return np.eye(3)+skew+skew@skew/(1.+cosine)


def refinement_metric_bridge():
    """Transport the finite carrier metric through the actual source screens.

    Port-hop counts are the dimensionless graph metric of the native seam
    network.  Only the root carrier has the separately recovered Euclidean
    Gram metric; this routine checks whether the current construction fixes
    its frame relative to neighboring carriers.
    """
    levels = (3, 5)
    balls = tuple(screen_compatibility.root_ball(level, 10) for level in levels)
    distances = tuple(_ball_distances(ball) for ball in balls)
    views9 = tuple(_ball_view(ball, distance, 9) for ball, distance in zip(balls, distances))
    views10 = tuple(_ball_view(ball, distance, 10) for ball, distance in zip(balls, distances))
    if views9[0][:2] != views9[1][:2]:
        raise ArithmeticError('the declared source screens disagree inside the matched radius-nine ball')
    if views10[0][1] == views10[1][1]:
        raise ArithmeticError('the source-screen continuation control no longer distinguishes radius ten')

    recovery = native_geometry.recover(native_geometry.execute_calibration())
    gram = np.asarray(recovery['gram_numeric'], dtype=float)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    coordinates = eigenvectors[:, -3:] @ np.diag(np.sqrt(eigenvalues[-3:]))
    # Select an actual inter-carrier seam at the matched root.  Its endpoint
    # label supplies no cross-carrier Gram entry or local frame transition.
    root_cell = balls[0].root
    graph = screen_compatibility.screen(levels[0]).graph.neighbour_lists
    cross_edge = next((tuple(sorted((port, neighbour)))
                       for port in range(12*root_cell, 12*root_cell+12)
                       for neighbour in graph[port] if neighbour//12 != root_cell))
    root_endpoint = next(port for port in cross_edge if port//12 == root_cell)
    neighbour_endpoint = next(port for port in cross_edge if port//12 != root_cell)
    root_port, neighbour_port = root_endpoint % 12, neighbour_endpoint % 12
    u, v = coordinates[root_port], coordinates[neighbour_port]
    # Two relative orientations both put this cross-seam pair in the already
    # recovered shortest shell.  A free rotation about the aligned endpoint
    # changes the other cross-carrier Gram entries while preserving both
    # carriers' complete internal Gram matrices.
    target_dot = 1./np.sqrt(5.)
    reference_axis = np.eye(3)[int(np.argmin(np.abs(u)))]
    tangent = reference_axis-np.dot(reference_axis, u)*u
    tangent /= np.linalg.norm(tangent)
    target = target_dot*u+np.sqrt(1.-target_dot**2)*tangent
    first_rotation = _rotation_between(v, target)
    angle = np.pi/2
    cross_target = np.array([[0., -target[2], target[1]],
                             [target[2], 0., -target[0]],
                             [-target[1], target[0], 0.]])
    about_target = (np.cos(angle)*np.eye(3)+(1.-np.cos(angle))*np.outer(target, target)
                    +np.sin(angle)*cross_target)
    second_rotation = about_target@first_rotation
    first_frame = coordinates@first_rotation.T
    second_frame = coordinates@second_rotation.T
    local_gram_residual = max(float(np.max(np.abs(coordinates@coordinates.T-gram))),
                              float(np.max(np.abs(first_frame@first_frame.T-gram))),
                              float(np.max(np.abs(second_frame@second_frame.T-gram))))
    first_link_dot = float(np.dot(u, first_frame[neighbour_port]))
    second_link_dot = float(np.dot(u, second_frame[neighbour_port]))
    cross_block_difference = float(np.max(np.abs(coordinates@first_frame.T
                                                 -coordinates@second_frame.T)))
    if (local_gram_residual > 2e-12
            or abs(first_link_dot-target_dot) > 2e-12
            or abs(second_link_dot-target_dot) > 2e-12
            or cross_block_difference < 1e-3):
        raise ArithmeticError('the independent-frame gluing control is not discriminating')

    common_ports, common_edges, common_addresses = views9[0]
    counts = [sum(hops <= radius for hops in distances[0].values()) for radius in range(11)]
    if counts != [sum(hops <= radius for hops in distances[1].values())
                  for radius in range(11)]:
        raise ArithmeticError('finite root graph-hop shells disagree across the two screens')
    first_rows, second_rows = dict(views10[0][1]), dict(views10[1][1])
    witness = next(port for port in sorted(first_rows.keys() | second_rows.keys())
                   if first_rows.get(port) != second_rows.get(port))
    root_maps = []
    for ball in balls:
        address_map = dict(ball.addresses)
        root_maps.append([address_map[port] for port in range(12)])
    address_maps = tuple(dict(view[2]) for view in views9)
    if (set(address_maps[0]) != set(address_maps[1])
            or any(address_maps[0][port] % 12 != port % 12
                   or address_maps[1][port] % 12 != port % 12
                   for port in common_ports)):
        raise ArithmeticError('a matched screen address changed its native carrier port label')
    tower_certificates = []
    for coarse, fine in ((0, 1), (1, 2)):
        graph_map = observer_refinement.graph_parent_certificate(coarse, fine)
        defect = observer_refinement.defect_norm(coarse, fine)
        if graph_map['cell_parent_graph_morphism'] is not True or defect != 1:
            raise ArithmeticError('the retained cell-ancestry obstruction changed')
        tower_certificates.append({**graph_map,
            'uniform_per_edge_laplacian_intertwining_defect_inf': str(defect),
            'coordinate_frame_or_metric_scaling_derived': False})
    source_routing = []
    for level in levels:
        federation = screen_compatibility.screen(level)
        receipt = federation.gluing_receipt
        source_routing.append({'level': level,
            'routing_mode': receipt['routing_mode'],
            'source_derived': receipt['source_derived'],
            'uses_global_cell_centres_and_local_tangent_frame_convention': True,
            'intercarrier_seams': federation.inter_count,
            'antipodal_consistent_pairs': receipt['antipodal_consistent_pairs'],
            'antipodal_consistent_fraction_exact': str(F(
                receipt['antipodal_consistent_pairs'], federation.inter_count))})
    return {
        'source_levels': list(levels),
        'root_cells': [ball.root for ball in balls],
        'root_port_global_addresses': root_maps,
        'recovered_root_carrier_gram_equal_across_levels': True,
        'all_matched_source_address_ports_preserved': True,
        'root_internal_seams_match_the_record_recovered_icosahedral_shell': True,
        'matched_root_graph_ball': {'radius_in_native_seam_hops': 9,
            'ports': len(common_ports), 'induced_seams': sum(map(len, (row for _, row in common_edges)))//2,
            'canonical_address_and_adjacency_equal': True,
            'root_port_metric_is_exactly_the_same_unitless_chord_metric': True},
        'root_graph_hop_ball_port_counts_radii_zero_through_ten': counts,
        'full_federation_native_hop_growth': _full_screen_hop_profiles(),
        'first_induced_graph_difference': {'radius': 10, 'canonical_port': witness,
            'source_level_3_neighbours': list(first_rows.get(witness, ())),
            'source_level_5_neighbours': list(second_rows.get(witness, ()))},
        'existing_federation_refinement_transport': {
            'map': 'cell ancestry with the same port label copied to every child',
            'cell_ancestry_map_record_derived': False,
            'level_certificates': tower_certificates,
            'path_independent_port_frame_or_metric_scaling_recovered': False},
        'existing_intercarrier_routing': {
            'level_receipts': source_routing,
            'routing_is_derived_from_native_observer_records': False,
            'record_derived_path_independent_frame_transition_available': False},
        'intercarrier_frame_nonuniqueness_control': {
            'conditional_equal_shortest_shell_length_used_only_as_control': True,
            'actual_level_3_cross_seam_global_ports_root_then_neighbor': [
                root_endpoint, neighbour_endpoint],
            'local_port_labels': [root_port, neighbour_port],
            'two_candidate_squared_link_lengths': [
                '2-2/sqrt(5)', '2-2/sqrt(5)'],
            'max_local_gram_residual_diagnostic': local_gram_residual,
            'cross_carrier_gram_entries_differ_between_candidates': cross_block_difference > 1e-3,
            'same_local_gram_and_same_link_length_leave_relative_frame_underdetermined': True},
        'operational_readout': 'exact native seam-hop neighborhoods plus an exact finite Euclidean metric within each root carrier',
        'physical_global_space_or_continuum_recovered': False,
        'physical_length_or_intercarrier_frame_calibration_supplied': False,
    }


def run():
    return {'schema': 'oph.relative_observer_geometry.v1', 'all_pass': True,
            'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            'exact_continuity': branch_audit(), 'complete_execution': runtime_check(),
            'complete_local_histories': history_check(), 'shared_event_clocks': clock_checks(),
            'record_geometry_connection': geometry_connection(), 'negative_controls': controls(),
            'scope': {'source_native_rule_carrier_and_record_services_supplied': True,
                'absolute_fields_observer_readable_given_known_initial_local_preparation': True,
                'hidden_initial_fields_determined_by_current_records_alone': False,
                'shared_prefix_is_a_record_cut_not_spacetime_simultaneity': True,
                'fundamental_global_clock_or_field_assumed': False,
                'observer_modular_flow_attached_to_these_local_ticks': False,
            'finite_carrier_spatial_metric_recovered': True,
            'matched_screen_root_metric_and_native_hop_neighborhoods_certified': True,
            'intercarrier_frame_gluing_or_physical_continuum_recovered': False,
            'continuum_spacetime_or_physical_length_calibration': False,
                'new_physical_prediction_or_empirical_comparison': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Relative observer geometry receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')
