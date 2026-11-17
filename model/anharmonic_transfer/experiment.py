"""Seal all transferred predictions before opening later experimental values."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..anharmonic_observations.experiment import _compact_operations, _check_law
from . import source, certified_comparison
from .calibration import transferred_parameters, zero_time_displacement
from .data import load_metadata, load_zero, load_zero_tokens, load_validation

FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]
RUN_NAMES = ('primary', 'step_refinement', 'occupation_refinement',
             'preparation_adjusted', 'instantaneous_carrier_control',
             'harmonic_control', 'frozen_control')
SUPPORTED_PROTOCOL_DIGEST = 'bc32e46a581c225348f6d48fb0690730c3f7d0d6e8b9c7a7d7a6ae020c30dbc6'
SOURCE_PATHS = tuple('model/anharmonic_transfer/'+name+'.py' for name in
                    ('source', 'environment', 'bounds', 'calibration', 'data', 'experiment',
                     'certified_comparison')) + tuple(certified_comparison.CERTIFICATE_PINS) + (
    'model/anharmonic_observations/source.py',
    'model/anharmonic_observations/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'model/relativistic_scalar/exact.py')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def _configuration(protocol):
    # Any scientific or interpretive change requires an explicitly revised
    # supported protocol, rather than a receipt claiming unexecuted settings.
    if digest(protocol) != SUPPORTED_PROTOCOL_DIGEST:
        raise ValueError('Unsupported protocol revision')
    expected = {
        'apparatus': {'delta_radians_per_second': '2*pi*500',
                      'omega_radians_per_second': 'pi/0.000150',
                      'omega_carrier_radians_per_second': 'pi/0.000035',
                      'pulse_duration_seconds': .00005, 'nominal_step_seconds': .0002,
                      'theta': .8, 'dephasing_per_second': 18.},
        'numerics': {'primary_dimension': 32, 'primary_subdivisions': 128,
                     'refinement_dimension': 64, 'refinement_subdivisions': 256,
                     'quadrature_order': 16},
        'readout': {'h': .4, 'theta': '-pi/2'},
    }
    for section, entries in expected.items():
        for key, value in entries.items():
            actual = protocol[section].get(key)
            if isinstance(actual, bool) or actual != value:
                raise ValueError(f'Unsupported protocol field: {section}.{key}')
    cases = [{'label': label, 'potential_phase_pi': -phase,
              'executor_phi_pi': phase, 'preparation_sign': sign}
             for label, phase, sign in zip('abcde', (0., .05, .05, .05, .1), (-1, -1, 1, 0, -1), strict=True)]
    if protocol['cases'] != cases or protocol['runs'] != list(RUN_NAMES):
        raise ValueError('Unsupported case mapping or run inventory')
    # Preserve the previously reviewed calibration source, not a later fit.
    if (protocol['preparation']['prior_receipt'] != 'model/anharmonic_transfer/data/calibration_receipt.json.gz'
            or protocol['preparation']['prior_receipt_sha256'] !=
            '80504ae1e92e94de0e5f50484da4a5dbbdb0981389c5733066468b36e961f215'):
        raise ValueError('Unsupported preparation calibration source')
    return {'delta': 2*math.pi*500, 'omega': math.pi/.000150,
            'omega_carrier': math.pi/.000035, 'pulse_duration': .00005,
            'theta': .8, 'dephasing_rate': 18., 'quadrature_order': 16}


def _key(row):
    return row['case'], row['time_ms']


def _record(case, record, h):
    read = record['reads'][0]
    if len(record['reads']) != 1 or (read['u'], read['v'], read['theta']) != (0., h, -math.pi/2):
        raise ValueError('Only the declared finite-difference instrument is scored')
    _check_law(read, record['spin_probabilities'])
    return {'case': case, 'time_ms': record['steps']/5, 'steps': record['steps'],
            'published_coordinate_seconds': record['simulated_duration_seconds'],
            'drive_elapsed_seconds': record['drive_elapsed_duration_seconds'],
            'predicted_position_estimator': read['selected_mean']/(math.sqrt(2)*h),
            'predicted_spin_mean': read['selected_mean'],
            'predicted_detector_probabilities': deepcopy(read['selected_detector_probabilities']),
            'spin_detector_joint_probabilities': deepcopy(read['spin_detector_joint_probabilities']),
            'acceptance_probability': record['spin_probabilities'][0],
            'rejection_probability': record['spin_probabilities'][1],
            'exact_position_mean_diagnostic': record['selected_position_mean'],
            'ideal_finite_model_error': deepcopy(record['execution_error_budget'])}


def predict(*, progress=None):
    protocol_bytes = (FOLDER/'protocol.json').read_bytes()
    protocol = json.loads(protocol_bytes)
    apparatus = _configuration(protocol)
    metadata, zero = load_metadata(), load_zero()
    if metadata['labels'] != list('abcde'):
        raise ValueError('Case metadata differs from the protocol')
    steps = [round(time/.2) for time in metadata['times_ms']]
    if steps != list(range(0, 51, 2)):
        raise ValueError('Expected all 26 published coordinates')
    previous = transferred_parameters(protocol)
    occupation = certified_comparison.certified_inputs(protocol, previous, apparatus)
    h = protocol['readout']['h']
    preparations = {case['label']: {
        'transferred_alpha': case['preparation_sign']*abs(previous['alpha_real']),
        'zero_time_adjusted_alpha': zero_time_displacement(zero[case['label']], previous['thermal_occupation'], h),
        'thermal_occupation': previous['thermal_occupation'],
        'reported_zero_time_estimator': zero[case['label']],
    } for case in protocol['cases']}
    configurations = [
        ('primary', 32, 128, {}), ('step_refinement', 32, 256, {}),
        ('occupation_refinement', 64, 256, {}),
        ('preparation_adjusted', 32, 128, {}),
        ('instantaneous_carrier_control', 32, 128, {'finite_spin_duration': False}),
        ('harmonic_control', 32, 128, {'omega': 0., 'theta': 0.}),
    ]
    words, runs = {}, []
    for name, dimension, subdivisions, extra in configurations:
        rows, zero_rows, ledgers = [], [], {}
        for case in protocol['cases']:
            label = case['label']
            if progress:
                progress(f'Executing {name}/{label}: d={dimension}, subdivisions={subdivisions}', flush=True)
            prep = preparations[label]
            alpha = prep['zero_time_adjusted_alpha' if name == 'preparation_adjusted' else 'transferred_alpha']
            execution = source.run(dimension, alpha, prep['thermal_occupation'], steps,
                                   [source.position_setting(h)], subdivisions=subdivisions,
                                   phi=math.pi*case['executor_phi_pi'], **(apparatus | extra))
            for snapshot in execution['records']:
                row = _record(label, snapshot, h)
                (rows if snapshot['steps'] else zero_rows).append(row)
            ledgers[label] = _compact_operations(
                {key: value for key, value in execution.items() if key != 'records'}, words)
        runs.append({'name': name, 'dimension': dimension, 'subdivisions': subdivisions,
                     'predictions': rows, 'zero_time_predictions': zero_rows,
                     'execution': ledgers})
    primary = runs[0]
    frozen = []
    for row in primary['predictions']:
        zero_row = next(item for item in primary['zero_time_predictions'] if item['case'] == row['case'])
        frozen.append({**deepcopy(zero_row), 'time_ms': row['time_ms'],
                       'requested_steps': row['steps'], 'steps': 0,
                       'published_coordinate_seconds': row['published_coordinate_seconds']})
    runs.append({'name': 'frozen_control', 'predictions': frozen,
                 'zero_time_predictions': deepcopy(primary['zero_time_predictions']),
                 'execution': 'primary zero-time source instrument, without evolution'})
    return {'schema': 'oph.anharmonic_transfer.predictions.v1',
            'protocol': protocol, 'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
            'source_bindings': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
            'provenance': metadata['provenance'], 'prior_parameters': previous,
            'preparations': preparations, 'runs': runs, 'contracted_words': words,
            'certified_unrestricted_apparatus': occupation,
            'positive_time_outcomes_used': False}


def score(rows, observations):
    keys, actual_keys = [_key(r) for r in rows], [_key(r) for r in observations]
    if (len(keys) != len(set(keys)) or len(actual_keys) != len(set(actual_keys))
            or set(keys) != set(actual_keys) or not keys):
        raise ValueError('Each later observation must be predicted exactly once')
    lookup = {_key(r): r['observed_position_estimator'] for r in observations}
    scored = []
    for row in rows:
        observed = lookup[_key(row)]
        predicted = row['predicted_position_estimator']
        if row['time_ms'] <= 0 or not math.isfinite(observed) or not math.isfinite(predicted):
            raise ValueError('Finite predictions and positive-time observations required')
        scored.append({'case': row['case'], 'time_ms': row['time_ms'],
                       'observed_position_estimator': observed, 'residual': predicted-observed})
    metrics = {}
    for case in (*'abcde', 'combined'):
        subset = [r['residual'] for r in scored if case == 'combined' or r['case'] == case]
        if not subset:
            raise ValueError('Missing a declared case')
        metrics[case] = {'observations': len(subset),
                         'position_estimator_rmse': float(np.sqrt(np.mean(np.square(subset)))),
                         'maximum_absolute_residual': max(abs(v) for v in subset)}
    return {'metrics': metrics, 'rows': scored, 'statistical_significance_claimed': False}


def evaluate(prediction):
    apparatus = _configuration(prediction['protocol'])
    raw_protocol = (FOLDER/'protocol.json').read_bytes()
    if (prediction['protocol_sha256'] != hashlib.sha256(raw_protocol).hexdigest()
            or prediction['protocol'] != json.loads(raw_protocol)):
        raise ValueError('Prediction protocol and frozen protocol byte binding differ')
    if [r['name'] for r in prediction['runs']] != list(RUN_NAMES):
        raise ValueError('Prediction run inventory must match the complete protocol')
    if prediction['positive_time_outcomes_used']:
        raise ValueError('Later observations must not enter predictions')
    occupation = certified_comparison.certified_inputs(prediction['protocol'], prediction['prior_parameters'], apparatus)
    if digest(prediction.get('certified_unrestricted_apparatus')) != digest(occupation):
        raise ValueError('Sealed occupation certificate differs from its authenticated fixed apparatus')
    expected_targets = {(case, step/5) for case in 'abcde' for step in range(2, 51, 2)}
    expected_axes = {'primary': (32, 128), 'step_refinement': (32, 256),
                     'occupation_refinement': (64, 256), 'preparation_adjusted': (32, 128),
                     'instantaneous_carrier_control': (32, 128), 'harmonic_control': (32, 128)}
    for run in prediction['runs']:
        keys = [_key(row) for row in run['predictions']]
        if len(keys) != 125 or set(keys) != expected_targets:
            raise ValueError('Complete fixed prediction target inventory required before scoring')
        if run['name'] in expected_axes and (run['dimension'], run['subdivisions']) != expected_axes[run['name']]:
            raise ValueError('Refinement and control axes differ from the fixed protocol')
    seal = digest(prediction)
    observations = load_validation()
    certified_residuals = certified_comparison.compare(occupation, observations, load_zero_tokens())
    scores = {run['name']: score(run['predictions'], observations) for run in prediction['runs']}
    lookup = {r['name']: r for r in prediction['runs']}
    refinements = {}
    for name, left, right in (('step', 'primary', 'step_refinement'),
                              ('occupation', 'step_refinement', 'occupation_refinement')):
        a, b = lookup[left]['predictions'], lookup[right]['predictions']
        if [_key(r) for r in a] != [_key(r) for r in b]:
            raise ValueError('Refinement targets changed')
        differences = [abs(x['predicted_spin_mean']-y['predicted_spin_mean']) for x, y in zip(a, b, strict=True)]
        refinements[name] = {'runs': [left, right], 'maximum_binary_probability_change': max(differences)/2,
                             'maximum_position_estimator_change': max(differences)/(math.sqrt(2)*.4),
                             'rigorous_occupation_error_bound': False}
    checks = {'all_125_later_observations_retained': all(s['metrics']['combined']['observations'] == 125 for s in scores.values()),
              'no_positive_time_fit': not prediction['positive_time_outcomes_used']}
    result = {'schema': 'oph.anharmonic_transfer.results.v1', 'prediction': deepcopy(prediction),
              'prediction_payload_sha256_before_scoring': seal, 'scores': scores,
              'numerical_refinements': refinements, 'checks': checks,
              'certified_unrestricted_apparatus_comparison': certified_residuals,
              'all_pass': all(checks.values()),
              'all_pass_meaning': 'Execution and comparison contracts only, not empirical success.'}
    if not result['all_pass']:
        raise ArithmeticError('Transfer comparison contract failed')
    return result
