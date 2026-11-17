"""Join finite observer operations and retain the unchanged empirical test."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..anharmonic_transfer import calibration, data
from ..anharmonic_transfer import experiment as prior
from ..common_source.open_observer import SOURCE_PATHS as COMMON_SOURCES
from . import anharmonic, runtime


FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]
PRIMARY_STEPS = tuple(range(0, 51, 2))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def predict_laboratory(*, progress=None):
    """Execute every primary endpoint before loading any later observations."""
    protocol_bytes = (ROOT/'model/anharmonic_transfer/protocol.json').read_bytes()
    protocol = strict_json.loads(protocol_bytes.decode())
    settings = prior._configuration(protocol)
    preparation = calibration.transferred_parameters(protocol)
    rows, zero, costs = [], [], {}
    for case in protocol['cases']:
        label = case['label']
        if progress:
            progress(f'Primitive observer: all 26 independent endpoints, case {label}', flush=True)
        apparatus = anharmonic.Apparatus(
            dimension=32, subdivisions=128, phi=math.pi*case['executor_phi_pi'],
            h=protocol['readout']['h'], **settings)
        alpha = case['preparation_sign']*abs(preparation['alpha_real'])
        for step, history in anharmonic.snapshots(
                apparatus, alpha, preparation['thermal_occupation'], PRIMARY_STEPS):
            reading = anharmonic.readout(history)
            # Full output states survive in history; the compact receipt
            # stores their complete classical law and numerical state checks.
            purities = [float(np.trace(branch.density@branch.density).real)
                        for branch in history.branches]
            row = {
                'case': label, 'time_ms': step/5, 'steps': step,
                'published_coordinate_seconds': step*apparatus.published_step_duration,
                'drive_elapsed_seconds': step*apparatus.drive_step_duration,
                'predicted_position_estimator': reading['position_estimator'],
                'predicted_spin_mean': reading['selected_mean'],
                'spin_detector_joint_probabilities': reading['spin_detector_joint_probabilities'],
                'acceptance_probability': reading['acceptance_probability'],
                'rejection_probability': reading['rejection_probability'],
                'total_probability_mass': history.probability_mass,
                'complete_quantum_output_dimension': history.work_dimension,
                'visible_record_count': len(history.branches),
                'unnormalized_quantum_branch_purities': purities,
                'retained_quantum_state_reset_count': 0,
            }
            (rows if step else zero).append(row)
            if step == PRIMARY_STEPS[-1]:
                costs[label] = {
                    'resources_per_trial_including_preparation_upper_bounds': dict(history.branches[0].resources),
                    'preparation_cost_semantics': 'componentwise upper bound over the hidden capped thermal draw',
                    'final_spin_commitment_extra_native_pulses': 3,
                    'pointer_blank_and_commitment_physical_costs_included': False}
    expected = {(case, step/5) for case in 'abcde' for step in PRIMARY_STEPS[1:]}
    if len(rows) != 125 or {(row['case'], row['time_ms']) for row in rows} != expected:
        raise ArithmeticError('every declared positive-time target must be predicted once')
    return {
        'schema': 'oph.primitive_observer.laboratory_prediction.v1',
        'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
        'dimension': 32, 'subdivisions': 128, 'preparation': preparation,
        'predictions': rows, 'zero_time_predictions': zero, 'costs': costs,
        'positive_time_outcomes_used': False,
        'earlier_endpoint_measurements_applied_to_later_trials': False,
        'host_prefix_reuse': 'separate declared experiments share unread computed prefixes only',
        'initial_state_action_noise_detector_and_clock_calibration_supplied': True,
        'rejected_spin_pointer_statistics_empirically_validated': False,
        'complete_output_states_available_from_runtime': True,
        'compact_receipt_contains_full_density_matrices': False,
    }


def validate_laboratory_prediction(prediction):
    """Check the complete declared endpoint menu before opening observations."""
    if (prediction['positive_time_outcomes_used'] is not False
            or type(prediction['dimension']) is not int or prediction['dimension'] != 32
            or type(prediction['subdivisions']) is not int or prediction['subdivisions'] != 128):
        raise ValueError('the declared unfitted primary prediction is required')
    protocol_bytes = (ROOT/'model/anharmonic_transfer/protocol.json').read_bytes()
    if prediction['protocol_sha256'] != hashlib.sha256(protocol_bytes).hexdigest():
        raise ValueError('prediction is detached from the frozen apparatus protocol')
    protocol = strict_json.loads(protocol_bytes.decode())
    settings = prior._configuration(protocol)
    apparatus_by_case = {case['label']: anharmonic.Apparatus(
        phi=math.pi*case['executor_phi_pi'], h=protocol['readout']['h'], **settings)
        for case in protocol['cases']}
    for key, steps in (('predictions', PRIMARY_STEPS[1:]), ('zero_time_predictions', (0,))):
        rows = prediction[key]
        expected = {(case, step) for case in 'abcde' for step in steps}
        if (len(rows) != len(expected)
                or any(type(row['steps']) is not int for row in rows)
                or {(row['case'], row['steps']) for row in rows} != expected):
            raise ValueError('every declared endpoint must appear exactly once')
        for row in rows:
            apparatus = apparatus_by_case[row['case']]
            law = np.asarray(row['spin_detector_joint_probabilities'], dtype=float)
            if (law.shape != (2, 2) or not np.isfinite(law).all()
                    or np.min(law) < -2e-10 or abs(float(law.sum())-1) > 2e-9):
                raise ValueError('a complete normalized four-outcome law is required')
            accept, reject = map(float, law.sum(axis=1))
            if accept <= 0:
                raise ValueError('the scored conditional detector requires positive acceptance')
            mean = float((law[0, 0]-law[0, 1])/accept)
            fields = {'time_ms': row['steps']/5,
                      'published_coordinate_seconds': row['steps']*apparatus.published_step_duration,
                      'drive_elapsed_seconds': row['steps']*apparatus.drive_step_duration,
                      'acceptance_probability': accept, 'rejection_probability': reject,
                      'total_probability_mass': accept+reject,
                      'predicted_spin_mean': mean,
                      'predicted_position_estimator': mean/(math.sqrt(2)*.4)}
            for name, expected_value in fields.items():
                value = row[name]
                if (isinstance(value, bool) or not isinstance(value, (float, int))
                        or not math.isfinite(value) or abs(value-expected_value) > 2e-10):
                    raise ValueError(f'endpoint field {name} disagrees with its executed law')
            if (type(row['complete_quantum_output_dimension']) is not int
                    or row['complete_quantum_output_dimension'] != 64
                    or type(row['visible_record_count']) is not int or row['visible_record_count'] != 4
                    or type(row['retained_quantum_state_reset_count']) is not int
                    or row['retained_quantum_state_reset_count'] != 0):
                raise ValueError('the complete retained apparatus output is required')
            purities = np.asarray(row['unnormalized_quantum_branch_purities'], dtype=float)
            masses = law.ravel()
            if (purities.shape != (4,) or not np.isfinite(purities).all()
                    or np.any(purities < masses*masses/64-2e-9)
                    or np.any(purities > masses*masses+2e-9)):
                raise ValueError('quantum branch purities contradict their retained masses')


def evaluate_laboratory(prediction):
    """Seal first, score second; old predictions are an implementation oracle."""
    validate_laboratory_prediction(prediction)
    payload = deepcopy(prediction)
    seal = digest(payload)
    observations = data.load_validation()
    scores = prior.score(payload['predictions'], observations)
    # Read the old live receipt only AFTER the new complete prediction is
    # sealed. It neither initializes nor supplies an intermediate trajectory.
    old = strict_json.load(ROOT/'model/anharmonic_transfer/results.json')
    primary = next(run for run in old['prediction']['runs'] if run['name'] == 'primary')
    reference = {(row['case'], row['time_ms']): row
                 for row in primary['predictions']+primary['zero_time_predictions']}
    probability_residuals, estimator_residuals = [], []
    for row in payload['predictions']+payload['zero_time_predictions']:
        expected = reference[row['case'], row['time_ms']]
        probability_residuals.append(float(np.max(np.abs(
            np.asarray(row['spin_detector_joint_probabilities'])
            -np.asarray(expected['spin_detector_joint_probabilities'])))))
        estimator_residuals.append(abs(row['predicted_position_estimator']
                                      -expected['predicted_position_estimator']))
    maximum_probability = max(probability_residuals)
    maximum_estimator = max(estimator_residuals)
    if (not all(math.isfinite(value) for value in (*probability_residuals, *estimator_residuals))
            or maximum_probability > 2e-10 or maximum_estimator > 2e-9):
        raise ArithmeticError('the joined observer changed the frozen apparatus prediction')
    if digest(prediction) != seal or digest(payload) != seal:
        raise ArithmeticError('prediction changed during scoring')
    return {'prediction': payload, 'prediction_sha256_before_scoring': seal,
            'scores': scores, 'provenance': data.verify_provenance(),
            'independent_existing_executor_comparison': {
                'all_130_endpoint_complete_classical_laws_compared': True,
                'maximum_joint_probability_residual': maximum_probability,
                'maximum_position_estimator_residual': maximum_estimator,
                'old_receipt_loaded_only_after_prediction_seal': True,
                'full_quantum_state_compatibility_tested_separately': True},
            'empirical_exact_agreement_established': False,
            'physical_discrimination_from_conventional_quantum_mechanics': False}


def run(*, progress=None):
    from . import action, completion_report, observer_report, selection
    native_action = action.run()
    pointer_selection = selection.run()
    completion = completion_report.run()
    observer = observer_report.run()
    laboratory = evaluate_laboratory(predict_laboratory(progress=progress))
    sources = sorted({*prior.SOURCE_PATHS, *COMMON_SOURCES, *action.SOURCE_PATHS,
                      *(path.relative_to(ROOT).as_posix() for path in FOLDER.glob('*.py')
                        if not path.name.startswith('test_')),
                      'model/common_source/protocol.py', 'model/common_source/open_protocol.py',
                      'model/common_source/open_environment.py', 'model/common_source/quantum_native.py',
                      'model/common_source/operational.py', 'model/observer/native.py',
                      'code/source_selection_model/response.json',
                      'code/m1_source_realization/receipt.json',
                      'Lean/Geometry/SourceActionMeasure.lean',
                      'Lean/SourceLaw/WeightedSeams.lean',
                      'Lean/SourceLaw/VacuumIdentity.lean',
                      'Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean',
                      'model/anharmonic_transfer/protocol.json',
                      'model/anharmonic_transfer/results.json'})
    checks = {
        'native_relative_mean_action_recovered': native_action['all_pass'],
        'projective_pointer_selection_attached': pointer_selection['all_pass'],
        'fixed_seam_native_completion_attached': completion['all_pass'],
        'same_history_observer_interfaces_attached': observer['all_pass'],
        'all_125_later_observations_retained': laboratory['scores']['metrics']['combined']['observations'] == 125,
        'existing_empirical_law_preserved':
            laboratory['independent_existing_executor_comparison']['maximum_joint_probability_residual'] <= 2e-10,
        'no_later_data_refit': not laboratory['prediction']['positive_time_outcomes_used'],
    }
    return {'schema': 'oph.primitive_observer.results.v1',
            'source_bindings': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
            'native_action': native_action, 'pointer_selection': pointer_selection,
            'completion': completion, 'observer': observer, 'laboratory': laboratory,
            'checks': checks, 'all_pass': all(checks.values()),
            'all_pass_meaning': 'declared finite composition and replay contracts, not exact empirical agreement',
            'scope': {
                'one_shared_complete_instrument_runtime': True,
                'common_source_and_laboratory_are_distinct_supplied_apparatus_programs': True,
                'all_finite_OPH_theorem_hypotheses_discharged': False,
                'formal_OperationalObserver_tower_instantiated': False,
                'source_selected_physical_action_or_spacetime_recovered': False,
                'native_relative_mean_action_derived_under_declared_schedule': True,
                'all_universe_observations_reproduced': False}}
