"""Replay complete repair identification and continuing source histories."""
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from . import bridge, calibration, core, source
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..observer.source_scalar import compiler
from ..primitive_observer import runtime


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS, *calibration.SOURCE_PATHS,
    'model/repair_identification/bridge.py', 'model/repair_identification/experiment.py',
    'model/primitive_observer/runtime.py', 'model/common_source/open_protocol.py',
    'model/common_source/protocol.py',
    'tools/strict_json.py'))))
TOLERANCE = bridge.TOLERANCE


def candidate_tables():
    controls = core.controls()
    return {name: controls['E_n_quarter' if name == 'E_collision_half' else name]
            for name in source.CANDIDATES}


def _mass(amplitude):
    return float(np.vdot(amplitude, amplitude).real)


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _probe(candidate, descriptor, table):
    execution = source.execute_probe(candidate, descriptor.name)
    bridge.validate_probe_execution(execution, candidate, descriptor.name)
    actual, diagnostics = bridge.coarse_probe(execution)
    exact = bridge.exact_probe(table, descriptor.name)
    distance = bridge.distance(actual, exact)+diagnostics['full_output_projection_bound_numeric']
    probabilities = [math.fsum(_mass(branch.amplitude) for branch in execution.branches
        if (None, None) == label and branch.abort or
        not branch.abort and (branch.current, branch.record[-1]) == label) for label in bridge.MENU]
    observed = probabilities[bridge.outcome(descriptor.current, 1)]
    expected = bridge._complex(core.probability(table, descriptor)).real
    error = abs(observed-expected)
    mass = math.fsum(probabilities)
    _require(max(error, distance, abs(mass-1)) <= TOLERANCE,
             'emitted probe does not match the complete quantum output law')
    return {'probe': descriptor.name, 'predicted_joint_probability': str(core.probability(table, descriptor)),
            'executed_joint_probability': observed, 'joint_probability_residual': error,
            'complete_outcome_probabilities': probabilities, **diagnostics,
            'complete_source_to_exact_distance_bound_numeric': distance,
            'maximum_branch_resources': dict(bridge.maximum_cost(branch.resources for branch in execution.branches))}


def _fields(blocks):
    rows = []
    for record, rho in blocks.items():
        probability = float(np.trace(rho).real)
        currents = [bridge.MENU[index][0] or 0 for index in record]
        pointers = [bridge.MENU[index][1] or 0 for index in record]
        rows.append((probability, 1+int(record[0] != 6), sum(currents), pointers[0]*pointers[1]))
    return {'mean_attempt_clock': math.fsum(p*t for p, t, _, _ in rows),
            'mean_integrated_current': math.fsum(p*j for p, _, j, _ in rows),
            'mean_squared_integrated_current': math.fsum(p*j*j for p, _, j, _ in rows),
            'pointer_correlation': math.fsum(p*v for p, _, _, v in rows)}


def _history(candidate, probe, reference, table):
    raw = source.continue_history(candidate, probe, reference=reference)
    _require(bridge._original_equal((raw.candidate, raw.probe, raw.reference_dimension, raw.feedback),
                                   (candidate, probe, 2 if reference else 1, True)),
             'the continuing source history has a different requested protocol')
    actual, raw_diagnostics = bridge.coarse_raw(raw)
    exact = bridge.exact_history(table, probe, reference)
    canonical_history, canonical_diagnostics = bridge.run_history(candidate, probe, reference)
    canonical_output = {branch.record: branch.density for branch in canonical_history.branches}
    projection = raw_diagnostics['full_output_projection_bound_numeric']
    exact_distance = bridge.distance(actual, exact)+projection
    runtime_distance = bridge.distance(actual, canonical_output)+projection
    _require(max(exact_distance, runtime_distance, abs(raw_diagnostics['full_probability_mass']-1)) <= TOLERANCE,
             'a complete continuing source history disagrees with its exact or canonical output')
    # The canonical runtime intentionally omits the cost of a supplied input.
    # Charge that independently executed joint source preparation explicitly.
    supplied = source.preparation_resources(probe, True) if reference else ()
    canonical_maximum = bridge.maximum_cost(runtime.add_resources(branch.resources, supplied)
                                            for branch in canonical_history.branches)
    raw_maximum = bridge.maximum_cost(branch.resources for branch in raw.branches)
    _require(all(count <= dict(canonical_maximum).get(key, 0) for key, count in raw_maximum),
             'the canonical service envelope omitted a raw source resource')
    # Checkpoint byte lengths and floating matrix digests vary across BLAS
    # implementations. Validate locally, retain a portable boundedness claim.
    checkpoint_size = canonical_diagnostics.pop('checkpoint_bytes')
    return {'candidate': candidate, 'probe': probe, 'untouched_reference_dimension': 2 if reference else 1,
            'complete_visible_histories': len(actual), **raw_diagnostics,
            'canonical_runtime': {**canonical_diagnostics,
                'checkpoint_within_runtime_byte_limit': checkpoint_size <= runtime.MAX_CHECKPOINT_BYTES,
                'supplied_preparation_cost': dict(supplied),
                'maximum_branch_resource_envelope': dict(canonical_maximum)},
            'maximum_raw_branch_resources': dict(raw_maximum),
            'complete_source_to_exact_distance_bound_numeric': exact_distance,
            'complete_source_to_runtime_distance_bound_numeric': runtime_distance,
            'complete_visible_probabilities': [float(np.trace(value).real) for value in actual.values()],
            'exact_clock_and_field_moments': _fields(exact),
            'executed_clock_and_field_moments': _fields(actual)}


def run():
    tables = candidate_tables()
    rank = core.rank_certificate()
    candidates = {}
    for candidate, table in tables.items():
        probabilities = core.predictions(table)
        reconstructed = core.reconstruct(probabilities)
        _require(reconstructed == table, 'the complete quadrature inverse changed the table')
        candidates[candidate] = {'table': {key: str(value) for key, value in table.items()},
            'exact_inverse_matches': True,
            'source_probes': [_probe(candidate, probe, table) for probe in core.PROBES]}
    histories = [_history(candidate, probe, reference, tables[candidate])
                 for candidate, probe in (('C', 'g_re'), ('E', 'g_re'), ('C_phase_quarter', 'g_im'),
                                          ('E', 'n_re'), ('E_collision_half', 'n_re'))
                 for reference in (False, True)]
    centers = tuple(value.a for value in core.predictions(tables['E_collision_half']))
    tubes = calibration.analyze_tubes(centers, (F(1, 1000),)*9, tables, 2)
    _require(tubes['compatible_locked_candidates'] == ['E_collision_half'], 'collision control was not separated')
    return {'schema': 'oph.repair_identification.results.v1', 'all_pass': True,
        'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'source_control': compiler.source_summary(),
        'probe_word_sha256': {probe.name: dict(source.word_hashes(probe.name)) for probe in core.PROBES},
        'reference_preparation_word_sha256': {probe.name: dict(source.word_hashes(probe.name, True))['preparation']
                                             for probe in core.PROBES},
        'common_visible_outcomes': [list(item) for item in bridge.MENU],
        'rank': {'real_parameters': list(core.COORDINATES), 'rank': rank['rank'],
                 'jacobian': [[str(entry) for entry in row] for row in rank['jacobian']],
                 'determinant': str(rank['determinant']),
                 'admissible_interior_radius': str(rank['interior_coordinate_radius']),
                 'diagonal_dominance_margin': str(rank['interior_diagonal_dominance_margin'])},
        'candidates': candidates, 'continuing_histories': histories,
        'exact_probability_tube_control': tubes,
        'error_budget': {'probability_delta_to_instrument': 'min(1,7*delta)',
                         'common_tube_radius_to_two_call_history': 'min(1,28*radius)',
                         'two_call_example_radius': '1/1000', 'history_diameter': '7/250',
                         'numerical_guard': TOLERANCE,
                         'machine_roundoff_certified': False, 'laboratory_error_certified': False},
        'scope': {'identified_object': 'the classified coarse definite-current instrument on all 78 total-two states',
                  'executed_scope': '36 source probes and ten complete two-event histories, five with untouched references',
                  'empirical_confirmation': False, 'physical_repair_law_selected': False,
                  'physical_seconds_selected': False,
                  'assumptions': ['the record-law classification premises and same law at each calibrated repair',
                                  'declared source controls, coherent proper-code loading/export, blanks and classical reads'],
                  'not_identified': ['arbitrary finer coin/contrast records or coherent dilations',
                                     'physical apparatus calibration, physical time conversion and continuum field dictionary']}}
