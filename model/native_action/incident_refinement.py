"""Source-screen comparison at an observer's incident-event count.

All reported allowances are rational. No fitted time coordinate or sampled
confidence interval enters the certificate. See INCIDENT_REFINEMENT.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import screen_compatibility as screen
from . import stopped_control
from . import clock_histories

ROOT = Path(__file__).resolve().parents[2]


def integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('a bounded exact integer required')
    return value


def word_expectation(ticks, incident_seams, length, early_incident=0):
    """Expected ordered subsequences of ONE word before the kth incident hit.

    early_incident counts incident letters except for the word's final letter.
    The original total seam count cancels from the exact negative-binomial sum.
    """
    integer(ticks, 1, 128); integer(incident_seams, 1, 1000000)
    integer(length, 1, 256); integer(early_incident, 0, length-1)
    top = ticks+length-1-early_incident
    return F(math.comb(top, length) if top >= length else 0, incident_seams**length)


def slot_timeout_bound(slots, ticks, incident_seams, total_seams):
    """Rational upper bound for fewer than ticks successes in slots trials.

    Binomial probability is q^slots times a finite sum. q^slots<=exp(-slots*p),
    and the reciprocal of a positive Taylor partial sum bounds that exponential.
    This avoids enormous exact integers for a half-million comparison slots.
    """
    integer(slots, 1, 1000000); integer(ticks, 1, 128)
    integer(total_seams, 1, 1000000); integer(incident_seams, 1, total_seams)
    if ticks > slots:
        return F(1)
    if incident_seams == total_seams:
        return F(0)
    p = F(incident_seams, total_seams)
    polynomial = sum((math.comb(slots, j)*(p/(1-p))**j for j in range(ticks)), F(0))
    exponential_lower = sum(((slots*p)**j/math.factorial(j) for j in range(65)), F(0))
    return min(F(1), polynomial/exponential_lower)


def geometry():
    """Rebuild and check the actual complete source neighborhoods, not a hash."""
    left, right = (screen.root_ball(level, 9) for level in (3, 5))
    mapping = screen.identify(left, right)
    rows = dict(left.adjacency)
    ports = tuple(sorted(rows))
    root = set(range(12))
    observed = tuple(sorted(root | {q for p in root for q in rows[p]}))
    edges = tuple((p, q) for p in ports for q in rows[p] if p < q)
    incident = tuple(e for e in edges if root.intersection(e))
    exits = []
    for ball in (left, right):
        address = dict(ball.addresses)
        graph = screen.screen(ball.level).graph.neighbour_lists
        exits.append(tuple(len(graph[address[p]])-len(rows[p]) for p in ports))
    if exits[0] != exits[1]:
        raise ValueError('the two source boundary degree rows differ')
    if (len(ports), len(edges), len(observed), len(incident)) != (222, 561, 15, 33):
        raise ValueError('the retained source observation contract changed')
    return ports, edges, incident, observed, exits[0], tuple(sorted(mapping.items()))


def exit_counts(ports, edges, observed, exits, order):
    """Count all walks whose FIRST exit occurs at each length (repeats allowed)."""
    integer(order, 1, 128)
    # Internal API: only canonical immutable source data reaches this routine.
    indices = {p: i for i, p in enumerate(ports)}
    weight = [int(p in observed) for p in ports]
    counts = []
    for _ in range(order):
        counts.append(sum(a*b for a, b in zip(weight, exits, strict=True)))
        following = [0]*len(ports)
        for a, b in edges:
            i, j = indices[a], indices[b]
            following[i] += weight[j]
            following[j] += weight[i]
        weight = following
    return tuple(counts)


def exit_bound(counts, ticks, incident_seams, target_ports=15, degree=6):
    """One full-to-truncated complete quantum-output allowance, including R."""
    integer(ticks, 1, 128); integer(incident_seams, 1, 1000000)
    integer(target_ports, 1, 1000000); integer(degree, 1, 1000000)
    if (type(counts) is not tuple or not 1 <= len(counts) <= 128
            or any(type(c) is not int or c < 0 for c in counts)):
        raise ValueError('a nonempty bounded tuple of nonnegative walk counts required')
    order = len(counts)
    partial = sum((c*word_expectation(ticks, incident_seams, length)
                   for length, c in enumerate(counts, 1)), F(0))
    ratio = F(degree*(ticks+order+1), incident_seams*(order+2))
    if ratio >= 1:
        raise ValueError('the retained order does not certify a geometric tail')
    first = target_ports*degree**(order+1)*word_expectation(ticks, incident_seams, order+1)
    remainder = first/(1-ratio)
    return min(F(1), partial+remainder), remainder


def _density_step(vector, edges, ports, seam_count):
    """Directed floor of sum_e P_e / seam_count on a fixed integer grid.

    For the binary native law, mean density uses (I+swap_e)/2. Positive
    integer coefficients avoid subtraction and yield componentwise lower
    bounds. The caller establishes int64's worst-case coefficient bound.
    """
    index = {p: i for i, p in enumerate(ports)}
    degree = np.zeros(len(ports), dtype=np.int64)
    neighbours = np.zeros(len(ports), dtype=np.int64)
    for a, b in edges:
        i, j = index[a], index[b]
        degree[i] += 1; degree[j] += 1
        neighbours[i] += vector[j]; neighbours[j] += vector[i]
    return ((2*len(edges)-degree)*vector+neighbours)//(2*seam_count)


def density_response(ports, edges, incident, ticks, cap=512, bits=52):
    """Certified column of the stopped mean-density transfer K^ticks.

    This is a first-moment reduction, NOT a replacement for the quantum
    instrument or its complete record distribution.
    """
    integer(ticks, 1, 16); integer(cap, 1, 2048); integer(bits, 8, 52)
    if (type(ports) is not tuple or not ports or len(set(ports)) != len(ports)
            or any(type(p) is not int or p < 0 for p in ports)
            or ports != tuple(sorted(ports)) or ports[0] != 0
            or type(edges) is not tuple or type(incident) is not tuple
            or not edges or not incident):
        raise ValueError('a canonical finite density graph with root zero required')
    for collection in (edges, incident):
        if (any(type(e) is not tuple or len(e) != 2
                or any(type(p) is not int for p in e)
                or e[0] >= e[1] or not set(e) <= set(ports) for e in collection)
                or collection != tuple(sorted(set(collection)))):
            raise ValueError('distinct ordered graph seams required')
    if not set(incident) <= set(edges):
        raise ValueError('incident seams must belong to the source graph')
    grid = 1 << bits
    if 2*len(edges)*grid > np.iinfo(np.int64).max:
        raise ValueError('the exact integer coefficient may overflow int64')
    remote = tuple(e for e in edges if e not in set(incident))
    vector = np.zeros(len(ports), dtype=np.int64)
    vector[0] = grid
    responses = []
    for tick in range(1, ticks+1):
        waiting, completed = vector, np.zeros(len(ports), dtype=np.int64)
        for _ in range(cap):
            completed += _density_step(waiting, incident, ports, len(edges))
            waiting = _density_step(waiting, remote, ports, len(edges))
        vector = completed
        deficit = 1-F(sum(map(int, vector)), grid)
        timeout = 1-(1-F(len(remote), len(edges))**cap)**tick
        if not 0 <= timeout <= deficit <= 1:
            raise ArithmeticError('directed density certificate lost its mass bound')
        low = F(int(vector[0]), grid)
        responses.append({'scheduled_ticks': tick,
            'return_probability_interval': [str(low), str(min(F(1), low+deficit))],
            'unresolved_wait_allowance': str(timeout),
            'rounding_allowance': str(deficit-timeout),
            'total_column_l1_allowance': str(deficit)})
    return responses


def run():
    ports, edges, incident, observed, exits, mapping = geometry()
    counts = exit_counts(ports, edges, observed, exits, 48)
    rows = []
    for tick in (1, 2, 4, 8, 16):
        bound, remainder = exit_bound(counts, tick, len(incident))
        rows.append({'scheduled_ticks': tick, 'one_screen_to_ball_allowance': str(bound),
            'two_screen_complete_output_allowance': str(min(F(1), 2*bound)),
            'walk_series_remainder_allowance': str(remainder)})
    responses = density_response(ports, edges, incident, 4)
    rho = F(1, 3)
    for row in responses:
        tick = row['scheduled_ticks']
        spatial, _ = exit_bound(counts, tick, len(incident))
        lo, hi = map(F, row['return_probability_interval'])
        # On the product input, <n_0(0)n_0(k)> = rho^2+rho(1-rho)K^k[0,0].
        # The complete-history locality allowance is charged AFTER this
        # observable map, not multiplied by a guessed density prefactor.
        covariance = [rho*(1-rho)*lo, rho*(1-rho)*hi]
        row['ball_connected_two_time_field_interval'] = list(map(str, covariance))
        row['either_full_screen_connected_field_interval'] = [
            str(max(-rho*rho, covariance[0]-spatial)),
            str(min(rho*(1-rho), covariance[1]+spatial))]
        row['either_full_screen_joint_11_interval'] = [
            str(max(F(0), rho*rho+covariance[0]-spatial)),
            str(min(rho, rho*rho+covariance[1]+spatial))]
    histories = clock_histories.run()
    for example in (*histories['examples'], histories['forced_timeout_control']):
        example['unbounded_completion_vs_external_cutoff_allowance'] = str(slot_timeout_bound(
            example['external_comparison_slot_cap'], example['requested_scheduled_hits'],
            histories['scheduled_root_incident_seams'], histories['comparison_union_seams']))
    paths = (*screen.SOURCE_PATHS, *stopped_control.SOURCE_PATHS, *clock_histories.SOURCE_PATHS,
             'model/native_action/incident_refinement.py')
    return {'schema': 'oph.native_action.incident_refinement.v1',
        'geometry': {'source_levels': [3, 5], 'root_ports': 12,
            'record_and_field_support_ports': len(observed), 'incident_seams': len(incident),
            'matched_ports': len(ports), 'retained_seams': len(edges),
            'source_map_sha256': hashlib.sha256(canonical(mapping).encode()).hexdigest(),
            'boundary_exits': sum(exits), 'first_exit_walk_counts': list(counts)},
        'complete_stopped_quantum_history_comparison': rows,
        'finite_density_fields': {'preparation': 'independent Bernoulli occupation 1/3 at every port',
            'local_variance_at_every_size': '2/9', 'density_grid_bits': 52,
            'retained_attempt_cap_per_interval': 512,
            'maximum_retained_native_attempts_in_four_successful_intervals': 2048,
            'exact_integer_grid_units': 1 << 52, 'responses': responses,
            'expected_retained_native_attempts_for_four_ticks': str(F(4*len(edges), len(incident))),
            'expected_full_source_attempts_for_four_ticks': [
                str(F(4*screen.screen(level).seams, len(incident))) for level in (3, 5)]},
        'complete_quantum_clock_controls': stopped_control.run(),
        'full_source_clock_histories': histories,
        'scope': {'same_complete_local_initial_restriction_including_reference': True,
            'local_record_feedback_and_final_quantum_output_covered': True,
            'scheduler_attempt_count_not_exposed_to_observer': True,
            'field_covariance_execution_is_diagonal_sector': True,
            'uniform_full_source_scheduler_not_selected_subgraph_scheduler': True,
            'all_level_refinement_or_relativistic_field_limit_proved': False,
            'physical_duration_calibration_or_density_selection_derived': False},
        'source_bindings': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                            for p in sorted(set(paths))}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Incident-clock screen comparison and certified fields replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
