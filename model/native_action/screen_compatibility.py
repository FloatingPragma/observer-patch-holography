"""Exact observation neighborhoods of the retained screen federation.

Addresses are transported by a rooted graph isomorphism, preserving every
local port label. The comparison is finite and uses per-seam attempt intensity,
not equal global attempt counts or a physical clock identification.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..federated_observer import source
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..primitive_observer import runtime
from . import construction as local

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*local.SOURCE_PATHS,
    'model/native_action/screen_compatibility.py',
    'Lean/QFT/CPRestrictionNet.lean',
    'Lean/Tower/ConsensusTower.lean'))))


@lru_cache(maxsize=6, typed=True)
def screen(level):
    """Use the pinned builder unchanged, extending the executed level range."""
    if type(level) is not int or level not in range(6):
        raise ValueError('an exact screen level between zero and five required')
    graph = source._source().build_federation(level, 'port_pair')
    assert graph.carriers == 20*4**level and graph.ports == 12*graph.carriers
    assert graph.seams == 630*4**level
    assert graph.graph.degree.max() == 6 and graph.graph.degree.min() == 5
    return graph


@dataclass(frozen=True)
class Ball:
    level: int
    radius: int
    root: int
    # Canonical port -> original global port; local port labels are unchanged.
    addresses: tuple[tuple[int, int], ...]
    adjacency: tuple[tuple[int, tuple[int, ...]], ...]


def root_ball(level, radius):
    if type(radius) is not int or not 0 <= radius <= 24:
        raise ValueError('an exact graph radius between zero and twenty-four required')
    graph = screen(level)
    root = 4**level-1  # central child at each subdivision of source face zero
    cells = {root: 0}
    distances = {12*root+p: 0 for p in range(12)}
    queue = list(distances)
    adjacency = graph.graph.neighbour_lists
    for port in queue:
        if distances[port] >= radius:
            continue
        for other in sorted(adjacency[port], key=lambda w: (w//12 != port//12, w % 12)):
            if other not in distances:
                if other//12 not in cells:
                    cells[other//12] = len(cells)
                distances[other] = distances[port]+1
                queue.append(other)
    encode = lambda p: 12*cells[p//12]+p % 12
    addresses = tuple(sorted((encode(p), int(p)) for p in queue))
    rows = tuple(sorted((encode(p), tuple(sorted(encode(q) for q in adjacency[p]
                                                if q in distances))) for p in queue))
    return Ball(level, radius, root, addresses, rows)


def identify(first, second):
    """Check every induced edge; no graph hash substitutes for this proof."""
    if (type(first) is not Ball or type(second) is not Ball
            or first.radius != second.radius or first.adjacency != second.adjacency):
        raise ValueError('the declared rooted native-port balls do not agree')
    if any(ball != root_ball(ball.level, ball.radius) for ball in (first, second)):
        raise ValueError('a ball certificate must contain its complete source neighborhood')
    left, right = dict(first.addresses), dict(second.addresses)
    if left.keys() != right.keys():
        raise ArithmeticError('canonical address domains differ')
    mapping = {left[p]: right[p] for p in left}
    assert len(set(mapping.values())) == len(mapping)
    for a, b in mapping.items():
        assert a % 12 == b % 12
        original = {mapping[q] for q in screen(first.level).graph.neighbour_lists[a] if q in mapping}
        target = {q for q in screen(second.level).graph.neighbour_lists[b] if q in mapping.values()}
        assert original == target
    return mapping


def comparison_certificate(first, second, radius):
    left, right = root_ball(first, radius), root_ball(second, radius)
    mapping = identify(left, right)
    next_left, next_right = root_ball(first, radius+1), root_ball(second, radius+1)
    assert next_left.adjacency != next_right.adjacency
    a, b = dict(next_left.adjacency), dict(next_right.adjacency)
    witness = next(p for p in sorted(a.keys() | b.keys()) if a.get(p) != b.get(p))
    return {'source_levels': [first, second],
        'source_root_cells': [left.root, right.root],
        'global_ports': [screen(first).ports, screen(second).ports],
        'global_seams': [screen(first).seams, screen(second).seams],
        'matched_graph_radius': radius, 'matched_ports': len(mapping),
        'matched_induced_seams': sum(len(row) for _, row in left.adjacency)//2,
        'local_port_labels_preserved': True,
        'root_is_central_child_of_source_face_zero_at_each_level': True,
        'canonical_port_adjacency': [[p, list(row)] for p, row in left.adjacency],
        'source_address_map': [[a, b] for a, b in sorted(mapping.items())],
        'next_radius_identity_extension_fails': {
            'radius': radius+1, 'canonical_port': witness,
            'first_neighbours': list(a[witness]) if witness in a else None,
            'second_neighbours': list(b[witness]) if witness in b else None,
            'other_port_frame_isomorphisms_not_excluded': True},
        'unconditional_complete_output_bounds': [
            {'per_seam_Poisson_parameter': str(s),
             'trace_distance_allowance': str(min(F(1), 2*local.poisson_local_bound(12, 6, s, radius)))}
            for s in (F(1, 8), F(1, 4))]}


def density_return_moments(level, steps):
    """Exact root block of Q^n, Q=I-L/12; sparse integer arithmetic only.

    Q is a common density uniformization, not a native global-attempt map.
    The Poisson mean is six Q steps per unit per-seam parameter. These steps
    are not deterministic native global attempts.
    """
    graph = screen(level)
    if type(steps) is not int or not 0 <= steps <= 24:
        raise ValueError('a bounded exact return expansion order required')
    root = 12*(4**level-1)
    result = [[[0]*12 for _ in range(12)] for _ in range(steps+1)]
    for start in range(12):
        weights = {root+start: 1}
        for step in range(steps+1):
            for stop in range(12):
                result[step][stop][start] = weights.get(root+stop, 0)
            following = defaultdict(int)
            for p, value in weights.items():
                neighbours = graph.graph.neighbour_lists[p]
                following[p] += (12-len(neighbours))*value
                for q in neighbours:
                    following[q] += value
            weights = following
    return result


def density_certificate():
    # A port-ball of radius nine makes every returning path of length<=18
    # internal; degrees at the radius-nine boundary first matter at order19.
    left = density_return_moments(3, 24)
    right = density_return_moments(5, 24)
    first_difference = next((n for n, (a, b) in enumerate(zip(left, right, strict=True)) if a != b), None)
    assert first_difference is None or first_difference >= 19
    order = 24 if first_difference is None else first_difference-1
    witness = None
    if first_difference is not None:
        i, j = next((i, j) for i, j in product(range(12), repeat=2)
                    if left[first_difference][i][j] != right[first_difference][i][j])
        witness = {'local_readout_port': i, 'local_initial_port': j,
                   'first_numerator': left[first_difference][i][j],
                   'second_numerator': right[first_difference][i][j],
                   'response_coefficient_difference':
                       str(F(left[first_difference][i][j]-right[first_difference][i][j], 12**first_difference))}
    parameter = F(1, 4)
    mean = 6*parameter
    tail = mean**(order+1)/math.factorial(order+1)/(1-mean/(order+2))
    return {'source_levels': [3, 5], 'density_uniformization': 'Q=I-L/12',
        'equivalent_Poisson_density_semigroup': 'exp(-s L/2)=exp[6s(Q-I)]',
        'first_moment_domain': 'every binary occupation sector; not the full quantum process',
        'root_block_dimension': 12,
        'common_exact_return_numerators_through_order': order,
        'common_exact_return_numerators': left[:order+1],
        'denominator_at_order_n': '12**n',
        'first_difference_order_within_executed_domain': first_difference,
        'first_difference_witness': witness,
        'single_initial_root_port_to_single_final_root_port_response_bound': {
            'per_seam_Poisson_parameter': str(parameter),
            'absolute_probability_difference_bound': str(tail),
            'method': 'common return coefficients, then Poisson factorial-tail majorant',
            'root_prepared_particle_and_other_ports_empty': True,
            'not_a_complete_quantum_history_bound': True},
        'common_response_fitted_to_data': False,
        'root_occupation_measurement_or_initial_preparation_service_supplied': True}


def initial_state_certificate():
    """Nontrivial compatible binary states and the distinct dilute-two limit."""
    rho = F(1, 3)
    sigma = {i: rho**i.bit_count()*(1-rho)**(2-i.bit_count()) for i in range(4)}
    for law in ('C', 'E'):
        output = defaultdict(F)
        for i, weight in sigma.items():
            for (a, b), value in local.quantum.native_unit(law, i, i, (0, 1)).items():
                output[a, b] += weight*value
        assert dict(output) == {(i, i): weight for i, weight in sigma.items()}
        # Matrix-unit charge preservation is the branchwise modular covariance
        # identity, valid at every rho in (0,1), not just at the tested value.
        for i, j, current in product(range(4), range(4), (-1, 0, 1)):
            for a, b in local.quantum.native_unit(law, i, j, (0, 1), current):
                assert a.bit_count()-b.bit_count() == i.bit_count()-j.bit_count()
    m = 222
    sizes = [screen(level).ports for level in (3, 5)]
    distributions = [[F(math.comb(m, j)*math.comb(v-m, 2-j), math.comb(v, 2))
                      for j in range(3)] for v in sizes]
    distance = sum(abs(a-b) for a, b in zip(*distributions, strict=True))/2
    return {'compatible_reference': 'tensor product diag(1-rho,rho) on every binary port',
        'executed_reference_density': str(rho),
        'arbitrary_density_domain': '0<rho<1; density is a declared preparation',
        'endpoint_stationarity_verified_for_both_candidates': True,
        'full_tensor_domain_includes_all_global_number_sectors': True,
        'per_carrier_total_two_constraint_imposed': False,
        'local_embeddings': 'A -> A tensor I under the exact port-address map',
        'state_preserving_restriction': 'identity tensor product-reference expectation',
        'branchwise_modular_charge_covariance_checked': True,
        'modular_flow': 'sigma_u(A)=exp(i*u*log(rho/(1-rho))*N) A exp(-i*u*log(rho/(1-rho))*N)',
        'modular_parameter_identified_with_attempt_intensity': False,
        'current_POVM_and_fixed_total_observables_have_trivial_modular_flow': True,
        'uniform_global_total_two_control': {
            'global_ports': sizes, 'matched_ball_ports': m,
            'local_total_probabilities': [[str(p) for p in row] for row in distributions],
            'local_restriction_trace_distance': str(distance),
            'vacuum_trace_distance': [str(1-row[0]) for row in distributions],
            'fixed_m_limit_is_local_vacuum': True,
            'complete_output_bound_at_s_one_quarter_with_this_state_mismatch':
                str(min(F(1), distance+2*local.poisson_local_bound(12, 6, F(1, 4), 9)))}}


def _probe_paths():
    mapping = identify(root_ball(3, 9), root_ball(5, 9))
    graph = screen(3)
    root = 12*(4**3-1)
    a = next(p for p in range(root, root+12)
             if any(q//12 != p//12 for q in graph.graph.neighbour_lists[p]))
    b = next(q for q in graph.graph.neighbour_lists[a] if q//12 != a//12)
    c = min(q for q in graph.graph.neighbour_lists[a] if q//12 == a//12)
    d = min(q for q in graph.graph.neighbour_lists[b] if q//12 == b//12)
    path = (c, a, b, d)
    return path, tuple(mapping[p] for p in path)


def _groups(law, path, edge):
    """Compile actual sorted-global-endpoint calls and transport signed records."""
    local.quantum.candidate(law)
    first, second = path[edge], path[edge+1]
    seam = tuple(sorted((first, second)))
    sign = 1 if first < second else -1
    groups = {}
    for col in range(16):
        occupied = tuple(sorted(path[i] for i in range(4) if (col >> i) & 1))
        for coin in (0, 1):
            target, current = source.binary_step(occupied, seam, coin)
            target_mask = sum((1 << i) for i in range(4) if path[i] in target)
            label = current if law == 'C' else int(seam[0] in occupied)-int(seam[1] in occupied)
            # Reversing the endpoints complements the placement coin and
            # reverses both the signed current and imbalance labels. Coin and
            # imbalance remain hidden, but this bijection preserves Kraus groups.
            mapped_coin = coin if sign == 1 else 1-coin
            groups.setdefault((sign*current, mapped_coin, sign*label), local.quantum._zero(16))[target_mask][col] = F(1)
    return tuple(tuple((F(1, 2), matrix) for key, matrix in sorted(groups.items()) if key[0] == j)
                 or ((F(1), local.quantum._zero(16)),) for j in (-1, 0, 1))


def probe_program(level, law):
    if type(level) is not int or level not in (3, 5):
        raise ValueError('one of the two matched probe levels required')
    path = _probe_paths()[(3, 5).index(level)]
    identity = runtime.Instrument(kraus=((np.eye(16),),))
    phase = runtime.Instrument(kraus=((np.diag([(-1.)**(i & 1) for i in range(16)]),),))
    events = []
    for n, edge in enumerate((1, 0, 2, 1)):
        instrument = local._instrument(_groups(law, path, edge), visible=True)
        events.append(runtime.Event(f'edge_{n}', instruments=(instrument,)))
        if n == 1:
            events.append(runtime.Event('local_feedback', control=1,
                                        instruments=(identity, phase, identity)))
    initial = local._initial()
    reduced = np.trace(initial.reshape(16, 2, 16, 2), axis1=1, axis2=3)
    return runtime.Program(f'screen_{level}_{law}', 16, reduced, tuple(events),
        binding={'source_level': level, 'source_path': list(path),
            'native_seam_tape': [list(sorted((path[e], path[e+1]))) for e in (1, 0, 2, 1)],
            'record_orientation': 'positive current goes from path[e] to path[e+1]',
            'schedule': 'four disclosed seams, all other seams held; not the full uniform scheduler',
            'entangled_preparation_and_local_feedback_service_supplied': True}), initial


def probe_execution():
    paths = _probe_paths()
    rows = []
    for law in ('C', 'E'):
        programs, outputs, residuals = [], [], []
        for level in (3, 5):
            program, initial = probe_program(level, law)
            programs.append(program)
            output = runtime.execute(program, initial=initial)
            prefix = runtime.execute(program, initial=initial, stop_after=2)
            continued = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
                program, runtime.checkpoint_to_json(prefix)))
            residuals.append(runtime.complete_distance(output, continued))
            outputs.append(output)
        first = {b.record: b.density for b in outputs[0].branches}
        second = {b.record: b.density for b in outputs[1].branches}
        distance = local.cq_distance(first, second)
        assert distance < 1e-12 and max(residuals) < 1e-12
        assert len(first) == 81
        rows.append({'candidate': law, 'source_levels': [3, 5],
            'source_paths': [list(path) for path in paths],
            'complete_current_words': len(first), 'untouched_reference_dimension': 2,
            'complete_record_and_terminal_state_distance_diagnostic': distance,
            'checkpoint_continuation_residual_diagnostic': max(residuals),
            'all_unselected_seams_held': True})
    return rows


def run():
    source.verify_source_pins()
    return {'schema': 'oph.native_action.screen_compatibility.v1',
        'finite_rooted_screen_comparisons': [comparison_certificate(a, b, r)
                                            for a, b, r in ((1, 3, 3), (2, 4, 6), (3, 5, 9))],
        'compatible_nontrivial_reference_and_dilution_control': initial_state_certificate(),
        'exact_density_return_prediction': density_certificate(),
        'executed_native_record_transport': probe_execution(),
        'scope': {'actual_pinned_port_pair_screen_source_used': True,
            'new_larger_elementary_carrier_required': False,
            'native_instruments_and_current_signs_transported': True,
            'matching_ball_initial_state_and_reference_restrictions_required': True,
            'fixed_total_two_preparation_is_not_compatible_across_screen_sizes': True,
            'all_level_compatible_exhaustion_or_refinement_limit_proved': False,
            'spherical_cell_parent_coarse_graining_repaired': False,
            'physical_time_dimension_or_relativistic_QFT_recovered': False,
            'independent_laboratory_comparison': False},
        'source_bindings': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                            for path in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native screen observation-neighborhood comparison replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
