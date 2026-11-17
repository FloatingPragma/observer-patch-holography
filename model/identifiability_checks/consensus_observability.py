"""What completed-history repair can identify in the existing laboratory model.

This attaches the finite-DAG result in common_source/HISTORY_REPAIR.md to the
actual mixed-state instrument runtime. Compatibility is rho_v=T_v(rho_parent)
with the prepared boundary fixed. T_v and the causal tree are supplied by the
declared apparatus, before solving any equations. A local update changes only
rho_v. Induction on depth proves uniqueness and closure after at most depth
covering rounds, from every finite initial assignment. No contraction or
physical-time hypothesis is needed for this finite acyclic statement.

Therefore every terminal measurement, conditional readout and data score is
identical for all completed schedules, with the SAME apparatus. This theorem
does not choose the apparatus maps: two different actions can both close.
Intermediate inconsistent matrices are solver variables, not physical states.
The returned floating diagnostics do not certify physical service fidelity.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from tools import strict_json
from model.anharmonic_transfer import calibration, data, experiment as transfer
from model.anharmonic_observations.replay import _compare
from model.primitive_observer import anharmonic, runtime


ROOT = Path(__file__).resolve().parents[2]
TOLERANCE = 2e-9


@dataclass(frozen=True)
class Node:
    parent: int
    stage: int
    record: tuple[int, ...]
    outcome: int
    event: runtime.Event


@dataclass(frozen=True)
class Graph:
    program: runtime.Program
    nodes: tuple[Node, ...]
    terminals: tuple[int, ...]
    depth: int


def graph(program):
    """Read only declarations; no runtime trajectory or empirical data are used."""
    program = runtime.validate_program(program)
    nodes, frontier = [], [(0, ())]
    for stage, event in enumerate(program.events):
        following = []
        for parent, record in frontier:
            for outcome in range(event.outcomes):
                child_record = (*record, outcome) if event.visible else record
                nodes.append(Node(parent, stage, child_record, outcome, event))
                following.append((len(nodes), child_record))
        frontier = following
    return Graph(program, tuple(nodes), tuple(index for index, _ in frontier),
                 len(program.events))


def _channel(instrument, density, outcome):
    # Independent of runtime._apply: hidden alternatives remain a CP sum.
    if instrument.schur is not None:
        return instrument.schur*density
    result = np.zeros_like(density)
    for operator in instrument.kraus[outcome]:
        result += operator@density@operator.conj().T
    return result


def edge(system, index, density):
    node = system.nodes[index-1]
    event = node.event
    if event.sequence:
        result = density.copy()
        for _ in range(event.repeat):
            for instrument in event.sequence:
                result = _channel(instrument, result, 0)
        return result
    # The parent's record, not the new outcome, selects any feedback map.
    record = () if node.parent == 0 else system.nodes[node.parent-1].record
    choice = 0 if event.control is None else record[event.control]
    return _channel(event.instruments[choice], density, node.outcome)


def schedule(system, kind, *, covers=None, seed=1907):
    if kind not in ('forward', 'reverse', 'shuffled'):
        raise ValueError('declared schedule kind required')
    covers = (1 if kind == 'forward' else system.depth) if covers is None else covers
    if type(covers) is not int or covers < 0:
        raise ValueError('nonnegative integer cover count required')
    rng = np.random.default_rng(seed)
    result = []
    for _ in range(covers):
        order = list(range(1, len(system.nodes)+1))
        if kind == 'reverse':
            order.reverse()
        elif kind == 'shuffled':
            rng.shuffle(order)
        result.extend(order)
    return tuple(result)


def repair(system, updates, initial=None):
    """Apply local equalities; arbitrary finite matrices are admitted guesses.

    Memoization only reuses an edge on byte-identical input. It never accepts
    solved states, empirical responses or a direct-execution oracle.
    """
    expected = set(range(1, len(system.nodes)+1))
    if initial is not None and set(initial) != expected:
        raise ValueError('initial guesses must contain every nonboundary node only')
    values = {0: system.program.initial_density.copy()}
    for index in expected:
        value = (np.zeros_like(values[0]) if initial is None else
                 np.asarray(initial[index], dtype=complex))
        if value.shape != values[0].shape or not np.isfinite(value).all():
            raise ValueError('finite common-size initial matrices required')
        values[index] = value.copy()
    cache = {}

    def target(index):
        parent = values[system.nodes[index-1].parent]
        key = index, parent.tobytes()
        if key not in cache:
            # Exact linearity shortcut; no physical state is renormalized.
            cache[key] = (edge(system, index, parent) if np.any(parent) else
                          np.zeros_like(parent))
        return cache[key]

    count = 0
    for index in updates:
        if type(index) is not int or index not in expected:
            raise ValueError('a local update must name a nonboundary node')
        values[index] = target(index).copy()
        count += 1
    residual = max((float(np.linalg.norm(values[i]-target(i), 'fro'))
                    for i in expected), default=0.)
    if not math.isfinite(residual):
        raise ArithmeticError('nonfinite constraint residual')
    return {'assignment': values, 'updates': count,
            'maximum_constraint_residual': residual,
            'closed_within_tolerance': residual <= TOLERANCE}


def terminal_law(system, result):
    return {(system.nodes[index-1].record if index else ()):
            float(np.trace(result['assignment'][index]).real)
            for index in system.terminals}


def apparatus_inputs(*, dimension=32, subdivisions=128, harmonic=False, case_label='a'):
    """Existing primary preparation and apparatus; no positive-time observations."""
    protocol = strict_json.load(ROOT/'model/anharmonic_transfer/protocol.json')
    settings = transfer._configuration(protocol)
    parameters = calibration.transferred_parameters(protocol)
    case = next(case for case in protocol['cases'] if case['label'] == case_label)
    apparatus = anharmonic.Apparatus(dimension=dimension, subdivisions=subdivisions,
                                     phi=math.pi*case['executor_phi_pi'],
                                     h=protocol['readout']['h'], **settings)
    if harmonic:
        apparatus = replace(apparatus, omega=0., theta=0.)
    return apparatus, case['preparation_sign']*abs(parameters['alpha_real']), parameters['thermal_occupation']


def selected_program(*, dimension=32, subdivisions=128, harmonic=False):
    """First positive endpoint of existing case a; no new fitted settings."""
    return anharmonic.program(*apparatus_inputs(dimension=dimension,
        subdivisions=subdivisions, harmonic=harmonic), 2)


def endpoint_graph(apparatus, alpha, thermal, steps):
    """Share unread prefixes of separate experiments, never measured states.

    Each endpoint has its own four-outcome read on its declared prefix. These
    endpoints are distinct prepared experiments, not sequential observations
    of a copied quantum state. Parent nodes are computed matrices, as in the
    existing snapshots API. Normalization is checked separately per endpoint.
    """
    if (type(steps) is not tuple or not steps or any(type(s) is not int or s < 0 for s in steps)
            or tuple(sorted(set(steps))) != steps):
        raise ValueError('distinct increasing nonnegative endpoint steps required')
    program = anharmonic.program(apparatus, alpha, thermal, steps[-1], read_steps=())
    base = graph(program)
    nodes = list(base.nodes)
    read = runtime.Event('separate_endpoint_read', instruments=(anharmonic.read_instrument(apparatus),))
    groups = {}
    for step in steps:
        parent = 5*step
        group = []
        for label in range(4):
            nodes.append(Node(parent, 5*step, (step, label), label, read))
            group.append(len(nodes))
        groups[step] = tuple(group)
    return Graph(program, tuple(nodes), tuple(i for group in groups.values() for i in group),
                 len(program.events)+1), groups


def _readout(probabilities, h):
    probabilities = np.asarray(probabilities)
    acceptance = float(probabilities[0]+probabilities[1])
    mass = float(probabilities.sum())
    if min(probabilities) < -TOLERANCE or abs(mass-1) > TOLERANCE or acceptance <= 0:
        raise ArithmeticError('complete endpoint law is not a normalized four-outcome law')
    mean = float((probabilities[0]-probabilities[1])/acceptance)
    return {'probabilities': probabilities.tolist(), 'acceptance_probability': acceptance,
            'total_probability_mass': mass, 'predicted_position_estimator': mean/(math.sqrt(2)*h)}


def all_primary_endpoints(progress=None):
    """Locally repair all 130 existing primary apparatus experiments.

    No positive-duration comparison responses or previously solved terminal
    states enter the graph or repair. The prepared boundary retains the prior
    empirical preparation calibration. Direct runtime snapshots validate full
    output matrices only after the repaired endpoint collection is computed.
    """
    rows, case_checks = [], []
    steps = tuple(range(0, 51, 2))
    for case in 'abcde':
        if progress is not None:
            progress(f'Repairing and checking all 26 primary endpoints in case {case}.')
        apparatus, alpha, thermal = apparatus_inputs(case_label=case)
        system, groups = endpoint_graph(apparatus, alpha, thermal, steps)
        result = repair(system, schedule(system, 'forward'))
        if not result['closed_within_tolerance']:
            raise ArithmeticError('full primary endpoint constraints did not close')
        for step, indices in groups.items():
            probabilities = [float(np.trace(result['assignment'][i]).real) for i in indices]
            rows.append({'case': case, 'time_ms': step/5,
                         **_readout(probabilities, apparatus.h)})
        difference = 0.
        for step, direct in anharmonic.snapshots(apparatus, alpha, thermal, steps):
            branches = {branch.record[-1]: branch.density for branch in direct.branches}
            for label, index in enumerate(groups[step]):
                difference = max(difference, float(np.linalg.norm(
                    result['assignment'][index]-branches[label], 'fro')))
        if difference > TOLERANCE:
            raise ArithmeticError('full primary repaired output differs from direct matrices')
        case_checks.append({'case': case, 'nodes_including_boundary': len(system.nodes)+1,
                            'updates': result['updates'], 'endpoints': len(groups),
                            'maximum_constraint_residual': result['maximum_constraint_residual'],
                            'maximum_complete_endpoint_matrix_difference': difference})
    if len(rows) != 130:
        raise ArithmeticError('the complete 130-endpoint collection is required')
    # Seal every prediction before positive-duration comparison data or old
    # endpoint predictions are opened. Prior preparation calibration is supplied.
    # This is reproducible reuse, not a new holdout.
    seal = hashlib.sha256(json.dumps(rows, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    old = strict_json.load(ROOT/'model/anharmonic_transfer/results.json')
    old_primary = next(run for run in old['prediction']['runs'] if run['name'] == 'primary')
    oracle = {(row['case'], row['time_ms']): row
              for row in old_primary['predictions']+old_primary['zero_time_predictions']}
    residual = max(float(np.max(np.abs(np.asarray(row['probabilities'])
        -np.asarray(oracle[row['case'], row['time_ms']]['spin_detector_joint_probabilities']).ravel())))
        for row in rows)
    if residual > TOLERANCE:
        raise ArithmeticError('completed repair changed the existing empirical endpoint law')
    # Positive-duration comparison data enter only after every prediction and
    # complete-state check; the separate preparation calibration remains supplied.
    scores = transfer.score([row for row in rows if row['time_ms'] > 0], data.load_validation())
    return {'endpoint_semantics': 'separate prepared experiments with shared unread computed prefixes',
            'endpoint_count': 130, 'scored_positive_duration_count': 125,
            'primary_occupation_dimension': 32, 'subdivisions': 128,
            'case_checks': case_checks, 'predictions': rows, 'scores': scores,
            'prediction_payload_sha256_before_scoring': seal,
            'maximum_retained_primary_probability_difference': residual,
            'new_independent_data': False, 'positive_duration_parameters_fitted': False,
            'rejected_branch_statistics_empirically_validated': False}


def audit_program(program):
    system = graph(program)
    results = {kind: repair(system, schedule(system, kind))
               for kind in ('forward', 'reverse', 'shuffled')}
    # The validation-only oracle is unavailable to graph extraction and repair.
    direct = runtime.execute(program)
    reference = {branch.record: branch.density for branch in direct.branches}
    rows = []
    for kind, result in results.items():
        difference = max(float(np.linalg.norm(result['assignment'][index]
                            -reference[system.nodes[index-1].record], 'fro'))
                         for index in system.terminals)
        law = terminal_law(system, result)
        if (not result['closed_within_tolerance'] or difference > TOLERANCE
                or abs(sum(law.values())-1) > TOLERANCE):
            raise ArithmeticError('completed repair disagrees with complete source execution')
        rows.append({'schedule': kind, 'updates': result['updates'],
                     'maximum_constraint_residual': result['maximum_constraint_residual'],
                     'maximum_terminal_matrix_difference': difference,
                     'complete_record_law': [{'record': list(key), 'probability': value}
                                              for key, value in sorted(law.items())]})
    unfair = repair(system, (system.terminals[-1],)*len(system.nodes)*system.depth)
    if unfair['closed_within_tolerance']:
        raise ArithmeticError('unfair unfinished assignment earned a closure certificate')
    readout = anharmonic.readout(direct)
    return {'node_count_including_boundary': len(system.nodes)+1,
            'causal_depth': system.depth, 'complete_terminal_records': len(reference),
            'schedules': rows, 'position_estimator': readout['position_estimator'],
            'unfair_control': {key: value for key, value in unfair.items() if key != 'assignment'}}, results


def run(progress=None):
    if progress is not None:
        progress('Checking fair schedules and the changed-action and unfair controls.')
    programs = {name: selected_program(harmonic=name == 'harmonic_control')
                for name in ('primary', 'harmonic_control')}
    audits = {name: audit_program(program)[0] for name, program in programs.items()}
    laws = [np.array([row['probability'] for row in audits[name]['schedules'][0]['complete_record_law']])
            for name in programs]
    tv = float(np.abs(laws[0]-laws[1]).sum()/2)
    if tv <= 1e-4:
        raise ArithmeticError('different supplied actions did not expose the retained contrast')
    full_primary = all_primary_endpoints(progress)
    # Every endpoint is now predicted. Reusing this one observed value adds
    # neither independent data nor a new fit.
    observed = next(row['observed_position_estimator'] for row in data.load_validation()
                    if row['case'] == 'a' and row['time_ms'] == .4)
    sources = sorted({Path(__file__).resolve().relative_to(ROOT).as_posix(),
        'model/primitive_observer/runtime.py', 'model/primitive_observer/anharmonic.py',
        'model/common_source/open_protocol.py', 'tools/strict_json.py',
        'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
        'model/anharmonic_transfer/protocol.json',
        'model/anharmonic_transfer/results.json',
        'model/anharmonic_observations/replay.py',
        'model/anharmonic_transfer/data/calibration_receipt.json.gz',
        'model/anharmonic_transfer/data/asymmetric_wells_expectation_x.csv',
        *transfer.SOURCE_PATHS})
    return {'schema': 'oph.consensus_observability.v1',
        'apparatus': {'case': 'a', 'published_time_ms': .4, 'occupation_dimension': 32,
                      'subdivisions': 128, 'complete_outcome': '2*spin+detector',
                      'actual_existing_empirical_endpoint': True},
        'runs': audits,
        'full_primary_empirical_comparison': full_primary,
        'changed_action_complete_law_total_variation': tv,
        'reused_observation': {'position_estimator': observed,
            'residuals': {name: audit['position_estimator']-observed
                          for name, audit in audits.items()},
            'new_independent_data': False, 'statistical_rejection_claimed': False},
        'theorem': {'nearest_result': 'model/common_source/HISTORY_REPAIR.md',
            'compatibility': 'fixed prepared root and rho_v=T_v(rho_parent)',
            'local_update': 'replace only rho_v by T_v(rho_parent)',
            'domain': 'finite acyclic instrument tree and arbitrary finite initial solver matrices',
            'schedule': 'node fairness; depth covering rounds suffice',
            'parameter': 'auxiliary repair count, not observer or laboratory duration',
            'endpoint_identification': 'all completed schedules give the same complete law and every common data score',
            'action_selection': 'uniqueness conditional on maps does not choose those maps'},
        'scope': {'all_solver_intermediate_matrices_are_physical_states': False,
                  'native_hardware_realizes_the_global_repair_solver': False,
                  'observer_meaning_maps_derived': False, 'physical_action_selected': False,
                  'new_empirical_confirmation': False, 'floating_error_certified': False,
                  'prior_empirical_preparation_calibration_supplied': True,
                  'positive_duration_comparison_outcomes_used_to_build_or_solve_constraints': False},
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    options = parser.parse_args()
    actual = run(progress=lambda message: print(message, file=sys.stderr, flush=True))
    if options.check:
        compare_receipt(strict_json.load(options.check), actual)
        print('Completed apparatus-history repair and empirical endpoint equivalence replayed.')
    elif options.write:
        with options.write.open('x') as handle:
            handle.write(json.dumps(actual, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        print(json.dumps(actual, indent=2, sort_keys=True, allow_nan=False))


def compare_receipt(expected, actual):
    """Validate each pre-score seal, then compare portable numerical outputs."""
    copies = []
    for receipt in (expected, actual):
        receipt = deepcopy(receipt)
        result = receipt['full_primary_empirical_comparison']
        seal = hashlib.sha256(json.dumps(result['predictions'], sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        if result.pop('prediction_payload_sha256_before_scoring') != seal:
            raise ValueError('complete endpoint prediction seal differs from its payload')
        copies.append(receipt)
    _compare(*copies)


if __name__ == '__main__':
    main()
