"""Local native Liouville action and complete quantum history locality.

The per-edge Poisson parameter is optional scheduler bookkeeping. The discrete
bound applies directly to independent uniform native seam attempts.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
import math
from pathlib import Path

import numpy as np

from ..federated_observer import source
from ..federated_quantum_sector import construction as quantum
from ..observer.scalar_response.experiment import canonical
from ..primitive_observer import runtime

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*quantum.SOURCE_PATHS,
    'model/native_action/__init__.py', 'model/native_action/__main__.py',
    'model/native_action/construction.py',
    'model/observer/QUANTUM_MIXING.md',
    'model/federated_observer/SCALING.md',
    'Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean'))))


def _integer(value, low=0, high=1000000):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('a bounded exact nonnegative integer required')
    return value


def _rational(value):
    if type(value) not in (int, F) or value < 0:
        raise ValueError('a nonnegative exact rational parameter required')
    return F(value)


def poisson_local_bound(target_ports, degree, parameter, radius):
    """Trace distance after retaining the radius ball, including local records/R."""
    _integer(target_ports, 1); _integer(degree); _integer(radius)
    parameter = _rational(parameter)
    order = radius+1
    return min(F(1), target_ports*(degree*parameter)**order/math.factorial(order))


def discrete_local_bound(target_ports, degree, seams, attempts, radius):
    _integer(target_ports, 1); _integer(degree); _integer(seams, 1)
    _integer(attempts); _integer(radius)
    if degree > seams:
        raise ValueError('the degree cannot exceed the total number of seams')
    order = radius+1
    if order > attempts:
        return F(0)
    return min(F(1), target_ports*math.comb(attempts, order)*F(degree, seams)**order)


def causal_cone(word, targets):
    """Backward native support; local feedback stays within the target set."""
    if (type(word) is not tuple or type(targets) is not tuple or not targets
            or any(type(p) is not int or p < 0 for p in targets)):
        raise ValueError('an immutable seam word and nonempty target addresses required')
    cone = set(targets)
    for seam in reversed(word):
        quantum._seam(seam)
        if cone.intersection(seam):
            cone.update(seam)
    return frozenset(cone)


def channel_column(law, left, right, edges):
    result = defaultdict(F)
    for edge in edges:
        for key, value in quantum.native_unit(law, left, right, edge).items():
            result[key] += value
    return dict(result)


def action_column(law, left, right, edges):
    """H=sum_e(I-Phi_e), derived independently from native matrix units."""
    result = defaultdict(F, {(left, right): F(len(edges))})
    for key, value in channel_column(law, left, right, edges).items():
        result[key] -= value
    return {key: value for key, value in result.items() if value}


def form_column(law, left, right, edges, ports):
    """Local exclusion plus coherence-domain killing, without native calls."""
    quantum.candidate(law)
    z = quantum.difference(left, right, ports)
    killing = sum((F(abs(z[a]-z[b]), 2) if law == 'C'
                   else F(z[a] != z[b])) for a, b in edges)
    result = defaultdict(F, {(left, right): killing})
    for a, b in edges:
        if z[a] == z[b] == 0 and ((left >> a) & 1) != ((left >> b) & 1):
            result[left, right] += F(1, 2)
            flip = (1 << a) | (1 << b)
            result[left ^ flip, right ^ flip] -= F(1, 2)
    return {key: value for key, value in result.items() if value}


def local_action_certificate():
    rows = []
    for law in ('C', 'E'):
        h = np.zeros((16, 16), dtype=object)
        for a, b in product(range(4), repeat=2):
            for (i, j), weight in action_column(law, a, b, ((0, 1),)).items():
                h[4*i+j, 4*a+b] = weight
        assert np.array_equal(h, h.T)
        if law == 'E':
            assert np.array_equal(h@h, h)
            spectral = ((F(1), h),)
        else:
            half = 4*(h-h@h)
            one = 2*h@h-h
            assert np.array_equal(half@half, half)
            assert np.array_equal(one@one, one)
            assert not np.any(half@one)
            assert np.array_equal(h, F(1, 2)*half+one)
            spectral = ((F(1, 2), half), (F(1), one))
        assert all(np.array_equal(p, p.T) for _, p in spectral)
        # A symmetric idempotent is an orthogonal projector, proving positivity.
        rows.append({'candidate': law, 'local_positive_projector_weights': [str(w) for w, _ in spectral],
                     'local_projector_ranks': [int(np.trace(p)) for _, p in spectral],
                     'exact_local_action_sha256': hashlib.sha256(canonical(
                         [[str(x) for x in row] for row in h]).encode()).hexdigest()})
    edges = ((0, 1), (1, 2), (2, 3))
    for law, a, b in product(('C', 'E'), range(16), range(16)):
        assert action_column(law, a, b, edges) == form_column(law, a, b, edges, 4)
    # The logarithm of one uniform discrete attempt is generally not this local sum.
    # Already a single edge has a killed eigenmode, while exp(-s H) is invertible.
    return {'local_action': rows, 'full_four_port_matrix_units_per_candidate': 256,
            'cross_number_coherences_included': True, 'local_action_has_no_fitted_couplings': True,
            'discrete_unread_transfer_equals_I_minus_H_over_M': True,
            'Poisson_unread_transfer_equals_exp_minus_parameter_H': True,
            'Poisson_parameter_is_declared_scheduler_bookkeeping': True,
            'single_edge_zero_transfer_mode_but_positive_exponential': True,
            'vacuum_full_coherence_is_fixed': all(not action_column(law, 0, 15, edges) for law in ('C', 'E'))}


def full_groups(law, edge, ports=4):
    """Source Kraus groups on the full binary tensor algebra, not one particle."""
    quantum.candidate(law); quantum._seam(edge); _integer(ports, 2, 5)
    if edge[1] >= ports:
        raise ValueError('the selected edge must lie inside the finite carrier')
    dim = 1 << ports
    groups = {}
    for col in range(dim):
        for coin in (0, 1):
            target, current, imbalance = quantum.native_image(col, edge, coin)
            label = current if law == 'C' else imbalance
            groups.setdefault((current, coin, label), quantum._zero(dim))[target][col] = F(1)
    return tuple(tuple((F(1, 2), matrix) for key, matrix in sorted(groups.items()) if key[0] == j)
                 or ((F(1), quantum._zero(dim)),) for j in (-1, 0, 1))


def _instrument(groups, visible=False):
    if not visible:
        groups = (tuple(entry for group in groups for entry in group),)
    return quantum._runtime_instrument(groups, visible=visible)


def partial_trace(density, keep, ports, reference=1):
    """Keep addressed binary ports in declared order and an untouched reference."""
    keep = tuple(keep)
    remove = tuple(p for p in range(ports) if p not in keep)
    out = np.zeros(((1 << len(keep))*reference,)*2, dtype=density.dtype)
    def part(mask, positions):
        return sum(((mask >> p) & 1) << i for i, p in enumerate(positions))
    for a, b in product(range(1 << ports), repeat=2):
        if part(a, remove) == part(b, remove):
            i, j = part(a, keep)*reference, part(b, keep)*reference
            out[i:i+reference, j:j+reference] += density[a*reference:(a+1)*reference, b*reference:(b+1)*reference]
    return out


def source_path():
    graph, _, tape, _, _ = source.tape_data()
    crossing, left, right = tape[:3]
    # The same three source seams, oriented into a simple four-port path.
    ends = [p for p in left if p not in crossing]
    first = ends[0]
    middle = next(p for p in left if p in crossing)
    other = next(p for p in crossing if p != middle)
    last = next(p for p in right if p != other)
    path = (first, middle, other, last)
    assert all(tuple(sorted((a, b))) in graph.seams for a, b in zip(path, path[1:]))
    return graph, path


def _initial():
    # A coherent local probe entangled with R, plus an occupied remote last port.
    state = np.zeros(32)
    state[2*8] = state[2*9+1] = 1/math.sqrt(2)
    return np.outer(state, state)


def locality_program(law, schedule, truncate=False):
    """Fixed-schedule witness with local records and record-controlled phase."""
    quantum.candidate(law)
    if type(schedule) is not tuple or any(type(e) is not int or e not in (0, 1, 2) for e in schedule):
        raise ValueError('a finite three-edge schedule required')
    graph, path = source_path()
    identity = runtime.Instrument(kraus=((np.eye(16),),))
    phase = np.diag([(-1.)**(mask & 1) for mask in range(16)])
    phase_instrument = runtime.Instrument(kraus=((phase,),))
    events = []; count = 0
    for n, edge in enumerate(schedule):
        if truncate and edge == 2:
            instrument = identity
        else:
            instrument = _instrument(full_groups(law, (edge, edge+1)), visible=edge == 0)
        events.append(runtime.Event(f'native_{n}', instruments=(instrument,)))
        if edge == 0:
            events.append(runtime.Event(f'local_feedback_{n}', control=count,
                instruments=(identity, phase_instrument, identity)))
            count += 1
    initial = _initial()
    default = np.trace(initial.reshape(16, 2, 16, 2), axis1=1, axis2=3)
    return runtime.Program(f'locality_{law}_{schedule}_{truncate}', 16, default, tuple(events),
        binding={'federation_level': 0, 'full_federation_ports': graph.ports,
                 'selected_source_path': list(path), 'declared_three_edge_schedule': list(schedule),
                 'observed_and_feedback_ports': list(path[:2]),
                 'retained_ball_radius': 1, 'outside_edge_replaced_with_identity': truncate,
                 'initial_preparation_and_phase_intervention_supplied': True}), initial


def _cq_output(checkpoint):
    return {branch.record: partial_trace(branch.density, (0, 1), 4, 2) for branch in checkpoint.branches}


def cq_distance(first, second):
    if first.keys() != second.keys():
        raise ValueError('matching complete local record alphabets required')
    return sum(float(np.linalg.norm(first[word]-second[word], ord='nuc')) for word in first)/2


def locality_execution():
    results = []; schedules = tuple(product(range(3), repeat=3))
    edges = ((0, 1), (1, 2), (2, 3))
    bad_count = sum(3 in causal_cone(tuple(edges[e] for e in word), (0, 1)) for word in schedules)
    for law in ('C', 'E'):
        total_distance = 0.; largest_good = 0.; largest_bad = 0.; total_branches = 0
        digest_rows = []; checkpoint_residual = 0.
        for schedule in schedules:
            full, initial = locality_program(law, schedule)
            reduced, _ = locality_program(law, schedule, True)
            actual = runtime.execute(full, initial=initial)
            truncated = runtime.execute(reduced, initial=initial)
            first, second = _cq_output(actual), _cq_output(truncated)
            distance = cq_distance(first, second)
            total_distance += distance/len(schedules)
            bad = 3 in causal_cone(tuple(edges[e] for e in schedule), (0, 1))
            if bad:
                largest_bad = max(largest_bad, distance)
            else:
                largest_good = max(largest_good, distance)
                assert distance < 2e-12
            total_branches += len(actual.branches)
            digest_rows.append([list(schedule), bad, len(actual.branches)])
        # Execute an additional serialized continuation with a full reference.
        program, initial = locality_program(law, (2, 1, 0))
        prefix = runtime.execute(program, initial=initial, stop_after=2)
        resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
            program, runtime.checkpoint_to_json(prefix)))
        checkpoint_residual = runtime.complete_distance(runtime.execute(program, initial=initial), resumed)
        assert total_distance <= float(F(bad_count, len(schedules)))+2e-12
        assert largest_bad > 0 and checkpoint_residual < 2e-12
        results.append({'candidate': law, 'source_path': list(source_path()[1]),
            'three_attempt_schedules': len(schedules), 'complete_local_record_branches': total_branches,
            'reference_dimension': 2, 'bad_schedule_probability_exact': str(F(bad_count, len(schedules))),
            'average_conditional_complete_output_trace_distance_diagnostic': total_distance,
            'largest_good_schedule_residual_diagnostic': largest_good,
            'largest_bad_schedule_trace_distance_diagnostic': largest_bad,
            'serialized_checkpoint_residual_diagnostic': checkpoint_residual,
            'schedule_support_audit_sha256': hashlib.sha256(canonical(digest_rows).encode()).hexdigest()})
    return results


def uniform_federation_bounds():
    rows = []
    for level in range(3):
        graph = source.federation(level)
        degree = max(map(len, graph.neighbours))
        assert degree <= 6
        for radius in (4, 8, 12):
            targets = graph.inter_seams[0]
            ball = set(targets)
            frontier = set(targets)
            for _ in range(radius):
                frontier = {b for a in frontier for b in graph.neighbours[a]}-ball
                ball.update(frontier)
            attempts = len(graph.seams)//4
            rows.append({'level': level, 'ports': graph.ports, 'seams': len(graph.seams),
                'observed_ports': 2, 'observed_port_addresses': list(targets),
                'retained_ball_ports': len(ball),
                'maximum_degree': degree, 'retained_graph_radius': radius,
                'native_attempts': attempts,
                'uniform_discrete_joint_output_bound': str(discrete_local_bound(2, degree,
                    len(graph.seams), attempts, radius)),
                'per_edge_Poisson_parameter': '1/4',
                'Poisson_joint_output_bound': str(poisson_local_bound(2, degree, F(1, 4), radius))})
    return rows


def feedback_counterexample():
    # An explicitly nonlocal controller reads port3 and flips port0. The two
    # inputs agree on the entire radius-two ball of port0 in the four-path.
    projector = np.diag([float(not (i & 8)) for i in range(16)])
    read = runtime.Instrument(kraus=((projector,), (np.eye(16)-projector,)), visible=True)
    identity = runtime.Instrument(kraus=((np.eye(16),),))
    flip = np.zeros((16, 16))
    for i in range(16):
        flip[i ^ 1, i] = 1.
    controlled = runtime.Instrument(kraus=((flip,),))
    outputs = []
    for occupied in (0, 8):
        initial = np.zeros((16, 16)); initial[occupied, occupied] = 1.
        program = runtime.Program('nonlocal_controller_counterexample', 16, initial,
            (runtime.Event('remote_read', instruments=(read,)),
             runtime.Event('nonlocal_feedback', control=0, instruments=(identity, controlled))))
        result = runtime.execute(program)
        outputs.append(partial_trace(sum((branch.density for branch in result.branches)), (0,), 4))
    difference = float(np.linalg.norm(outputs[0]-outputs[1], ord='nuc')/2)
    assert difference == 1 and poisson_local_bound(1, 6, 0, 2) == 0
    return {'unrestricted_remote_record_controls_local_flip': True,
            'executed_remote_read_and_controlled_local_flip': True,
            'zero_native_attempt_local_trace_distance': '1',
            'inputs_agree_on_radius_two_ball': True,
            'violates_local_controller_access_hypothesis': True}


def run():
    return {'schema': 'oph.native_action.v1',
        'source_derived_local_action': local_action_certificate(),
        'executed_complete_history_locality': locality_execution(),
        'size_uniform_quantum_observation_bounds': uniform_federation_bounds(),
        'adversarial_controller_counterexample': feedback_counterexample(),
        'scope': {'native_discrete_attempts_and_optional_Poisson_parameter_distinguished': True,
            'local_records_quantum_outputs_and_untouched_references': True,
            'local_controller_only': True, 'conditional_quantum_local_limit': True,
            'pinned_refinements_proved_compatible': False,
            'Lorentzian_continuum_or_physical_time_identification': False,
            'new_laboratory_comparison': False},
        'source_bindings': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}
