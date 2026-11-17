"""Freeze source predictions before reading the positive-duration responses."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import source
from .calibration import calibrate
from .data import load_calibration, load_metadata, load_validation


FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]
SOURCE_PATHS = (
    'model/anharmonic_observations/source.py',
    'model/anharmonic_observations/calibration.py',
    'model/anharmonic_observations/data.py',
    'model/anharmonic_observations/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'model/scalar_observations/experimental_readout/executable_calibration.py',
)
RUN_NAMES = frozenset({'primary', 'occupation_refinement', 'dephasing_refinement',
                       'harmonic_control', 'coherent_control', 'frozen_control'})


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def _apparatus(protocol):
    """Resolve the fixed published apparatus without evaluating expressions.

    Changing the experiment requires a newly reviewed protocol; its numeric
    fields cannot silently diverge from defaults in the execution module.
    """
    expected = {'delta_radians_per_second': '2*pi*500',
                'omega_radians_per_second': 'pi/0.000150',
                'pulse_duration_seconds': .00005, 'step_duration_seconds': .0002,
                'spin_rotation_angle': .8, 'asymmetry_phase': 0.,
                'dephasing_per_second': 18., 'dephasing_operator': 'sqrt(18) * N'}
    apparatus = protocol['apparatus']
    for key, value in expected.items():
        if isinstance(apparatus.get(key), bool) or apparatus.get(key) != value:
            raise ValueError(f'Unsupported apparatus protocol value: {key}')
    if apparatus['step_duration_seconds'] != 4*apparatus['pulse_duration_seconds']:
        raise ValueError('Apparatus protocol step duration must equal four pulse durations')
    return {'delta': 2*math.pi*500, 'omega': math.pi/.000150,
            'pulse_duration': apparatus['pulse_duration_seconds'],
            'theta': apparatus['spin_rotation_angle'], 'phi': apparatus['asymmetry_phase'],
            'dephasing_rate': apparatus['dephasing_per_second']}


def _step(time_ms, step_duration_seconds=.0002):
    if (not math.isfinite(float(time_ms)) or not math.isfinite(step_duration_seconds)
            or step_duration_seconds <= 0):
        raise ValueError('finite time and positive step duration required')
    value = float(time_ms)/(1000*step_duration_seconds)
    result = round(value)
    if not math.isfinite(value) or abs(value-result) > 1e-10 or result < 0:
        raise ValueError('time must be a nonnegative multiple of 0.2 ms')
    return result


def _targets(metadata, *, positive_only=True):
    rows = []
    for family, components in (('characteristic', ('real', 'imag')), ('position', ('imag',))):
        item = metadata[family]
        for time in item['times_ms']:
            if positive_only and time <= 0:
                continue
            for beta in item['probe_coordinates']:
                for component in components:
                    rows.append({'family': family, 'time_ms': float(time),
                                 'beta_imag': float(beta), 'component': component})
    return rows


def _key(row):
    return row['family'], row['time_ms'], row['beta_imag'], row['component']


def _compact_operations(value, words):
    """Store identical contracted words once; retain their complete metadata."""
    if isinstance(value, dict):
        if 'word_sha256' in value:
            # Names distinguish descriptive roles; float diagnostics must not
            # determine exact dictionary keys across otherwise valid replays.
            identifier = _digest({key: value.get(key) for key in
                                  ('name', 'word_sha256', 'register_qubits')})
            if identifier in words and words[identifier] != value:
                raise ArithmeticError('Conflicting metadata for the same contracted source word')
            words.setdefault(identifier, deepcopy(value))
            return {'operation_metadata_sha256': identifier}
        return {key: _compact_operations(item, words) for key, item in value.items()}
    if isinstance(value, list):
        return [_compact_operations(item, words) for item in value]
    return value


def _check_law(read, spin_probabilities):
    probabilities = np.asarray(read['selected_detector_probabilities'], dtype=float)
    joint = np.asarray(read['spin_detector_joint_probabilities'], dtype=float)
    spin = np.asarray(spin_probabilities, dtype=float)
    if (probabilities.shape != (2,) or joint.shape != (2, 2) or spin.shape != (2,)
            or not np.isfinite(joint).all() or not np.isfinite(probabilities).all()
            or not np.isfinite(spin).all() or not math.isfinite(read['selected_mean'])
            or spin[0] <= 0 or np.min(spin) < -2e-10 or np.max(spin) > 1+2e-10
            or np.min(joint) < -2e-10 or np.max(joint) > 1+2e-10
            or np.min(probabilities) < -2e-10 or np.max(probabilities) > 1+2e-10
            or abs(sum(probabilities)-1) > 2e-10 or abs(np.sum(joint)-1) > 2e-10
            or abs(read['selected_mean']-(probabilities[0]-probabilities[1])) > 2e-10
            or not np.allclose(np.sum(joint, axis=1), spin, rtol=0, atol=2e-10)
            or not np.allclose(joint[0], probabilities*spin[0], rtol=0, atol=2e-10)):
        raise ArithmeticError('invalid complete selected-detector law')


def _extract(execution, targets, *, frozen=False, step_duration_seconds=.0002):
    snapshots, lookup = {}, {}
    for snapshot in execution['records']:
        step = snapshot['steps']
        if type(step) is not int or step < 0 or step in snapshots:
            raise ValueError('source snapshots require unique nonnegative integer steps')
        snapshots[step] = snapshot
        for read in snapshot['reads']:
            if (read['u'] != 0. or not math.isfinite(read['v'])
                    or read['theta'] not in (0., -math.pi/2)):
                raise ValueError('source readout requires u=0 and declared real/imaginary phases')
            key = step, read['v'], 'real' if read['theta'] == 0. else 'imag'
            if key in lookup:
                raise ValueError('duplicate source detector setting')
            lookup[key] = read
    target_keys = [_key(target) for target in targets]
    if len(target_keys) != len(set(target_keys)):
        raise ValueError('duplicate target statistic')
    rows = []
    for target in targets:
        step = 0 if frozen else _step(target['time_ms'], step_duration_seconds)
        snapshot = snapshots[step]
        read = lookup[(step, target['beta_imag'], target['component'])]
        probabilities = read['selected_detector_probabilities']
        joint = read['spin_detector_joint_probabilities']
        _check_law(read, snapshot['spin_probabilities'])
        rows.append({**target, 'predicted_mean': read['selected_mean'],
                     'predicted_probabilities': deepcopy(probabilities),
                     'spin_detector_joint_probabilities': deepcopy(joint),
                     'acceptance_probability': snapshot['spin_probabilities'][0],
                     'rejection_probability': snapshot['spin_probabilities'][1]})
    return rows


def predict(*, progress=None):
    """No call to load_validation occurs anywhere on this prediction path."""
    protocol_bytes = (FOLDER/'protocol.json').read_bytes()
    protocol = json.loads(protocol_bytes)
    apparatus = _apparatus(protocol)
    step_duration = 4*apparatus['pulse_duration']
    metadata = load_metadata()
    zero = load_calibration()
    scan = zero['characteristic']
    preparation = calibrate(scan['probe_coordinates'], scan['spin_real'][0], scan['spin_imag'][0])
    parameters = preparation['parameters']
    targets = _targets(metadata)
    zero_targets = [row for row in _targets(metadata, positive_only=False)
                    if row['time_ms'] == 0. and row['family'] == 'position']
    settings = sorted({(row['beta_imag'], 0. if row['component'] == 'real' else -math.pi/2)
                       for row in targets+zero_targets})
    settings = [dict(u=0., v=beta, theta=theta) for beta, theta in settings]
    steps = sorted({0, *(_step(row['time_ms'], step_duration) for row in targets)})
    numeric = protocol['numerics']
    dimension = numeric['primary_occupation_dimension']
    subdivision = numeric['dephasing_subdivisions_per_50us']
    configurations = [
        ('primary', dimension, subdivision, {}),
        ('occupation_refinement', numeric['occupation_check_dimension'], subdivision, {}),
        ('dephasing_refinement', dimension, numeric['dephasing_refinement_subdivisions_per_50us'], {}),
        ('harmonic_control', dimension, subdivision, {'omega': 0., 'theta': 0.}),
        ('coherent_control', dimension, subdivision, {'dephasing_rate': 0.}),
    ]
    runs = []
    for name, size, subdivisions, extra in configurations:
        if progress:
            progress(f'Executing {name}: occupation {size}, channel subdivisions {subdivisions}', flush=True)
        execution = source.run(size, parameters['alpha_real'], parameters['thermal_occupation'],
                               steps, settings, subdivisions=subdivisions, **(apparatus | extra))
        words = {}
        ledger = _compact_operations({key: value for key, value in execution.items()
                                      if key != 'records'}, words)
        ledger['contracted_words'] = words
        ledger['snapshots'] = [{key: value for key, value in snapshot.items() if key != 'reads'}
                               for snapshot in execution['records']]
        ledger['executed_detector_laws'] = sum(len(row['reads']) for row in execution['records'])
        rows = _extract(execution, targets, step_duration_seconds=step_duration)
        runs.append({'name': name, 'predictions': rows, 'source_execution': ledger})
        if name == 'primary':
            zero_prediction = _extract(execution, zero_targets, step_duration_seconds=step_duration)
            frozen = _extract(execution, targets, frozen=True, step_duration_seconds=step_duration)
    runs.append({'name': 'frozen_control', 'predictions': frozen,
                 'source_execution': {'origin': 'primary source preparation and t=0 detector laws',
                                      'evolution_executed': False}})
    result = {'schema': 'oph.anharmonic_observations.predictions.v1',
              'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
              'protocol': protocol, 'calibration': preparation,
              'source_binding_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                        for name in SOURCE_PATHS},
              'provenance': metadata['provenance'], 'runs': runs,
              'zero_duration_position_predictions': zero_prediction,
              'positive_duration_values_used': False}
    return result


def score(predictions, observations):
    expected = {}
    for family, components in (('characteristic', ('real', 'imag')), ('position', ('imag',))):
        data = observations[family]
        for name in ('times_ms', 'probe_coordinates'):
            values = np.asarray(data[name], dtype=float)
            if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
                    or len(set(values)) != len(values)):
                raise ValueError(f'Validation {name} must be finite and unique; duplicate coordinates prohibited')
        for component in components:
            if np.asarray(data['spin_'+component]).shape != (len(data['times_ms']), len(data['probe_coordinates'])):
                raise ValueError('Validation spin-mean shape does not match its coordinate grid')
        for index, time in enumerate(data['times_ms']):
            if time <= 0:
                raise ValueError('validation must contain only positive durations')
            for probe, beta in enumerate(data['probe_coordinates']):
                for component in components:
                    key = family, float(time), float(beta), component
                    expected[key] = float(data['spin_'+component][index, probe])
    keys = [_key(row) for row in predictions]
    if len(keys) != len(set(keys)) or set(keys) != set(expected):
        raise ValueError('predictions must cover every positive-duration statistic exactly once')
    rows = []
    for prediction in predictions:
        observed = expected[_key(prediction)]
        predicted = prediction['predicted_mean']
        if not math.isfinite(observed) or abs(observed) > 1 or not math.isfinite(predicted) or abs(predicted) > 1+2e-10:
            raise ValueError('finite binary spin means required')
        _check_law({'selected_mean': predicted,
                    'selected_detector_probabilities': prediction['predicted_probabilities'],
                    'spin_detector_joint_probabilities': prediction['spin_detector_joint_probabilities']},
                   [prediction['acceptance_probability'], prediction['rejection_probability']])
        rows.append({**prediction, 'observed_mean': observed,
                     'observed_probability_zero': (1+observed)/2,
                     'residual': predicted-observed})
    result = {}
    for family in ('characteristic', 'position', 'combined'):
        selected = [row for row in rows if family == 'combined' or row['family'] == family]
        rmse = float(np.sqrt(np.mean([row['residual']**2 for row in selected])))
        result[family] = {'scalar_means': len(selected), 'spin_mean_rmse': rmse,
                          'binary_probability_rmse': rmse/2,
                          'maximum_absolute_residual': max(abs(row['residual']) for row in selected)}
        if family != 'combined':
            result[family]['rows'] = selected
    positions = []
    for time in sorted({row['time_ms'] for row in result['position']['rows']}):
        selected = [row for row in result['position']['rows'] if row['time_ms'] == time]
        beta = np.asarray([row['beta_imag'] for row in selected])
        centered = beta-beta.mean()
        positions.append({'time_ms': time, **{
            label+'_position': float(centered@np.asarray([row[key] for row in selected])/
                                      (centered@centered)/math.sqrt(2))
            for label, key in (('predicted', 'predicted_mean'), ('observed', 'observed_mean'))}})
    result['derived_position'] = positions
    result['statistical_significance_claimed'] = False
    return result


def evaluate(prediction):
    """Unseal raw positive-duration responses only after all predictions exist."""
    names = [run['name'] for run in prediction['runs']]
    if len(names) != len(set(names)) or set(names) != RUN_NAMES:
        raise ValueError('Prediction run inventory must contain each declared run exactly once')
    prediction_hash = _digest(prediction)
    observations = load_validation()
    runs = []
    for run in prediction['runs']:
        comparison = score(run['predictions'], observations)
        runs.append({key: value for key, value in run.items() if key != 'predictions'} |
                    {'comparison': comparison})
    primary = next(run['predictions'] for run in prediction['runs'] if run['name'] == 'primary')
    refinements = {}
    for run in prediction['runs']:
        if run['name'] in ('occupation_refinement', 'dephasing_refinement'):
            if [_key(row) for row in run['predictions']] != [_key(row) for row in primary]:
                raise ValueError('refinement settings changed')
            differences = [abs(a['predicted_mean']-b['predicted_mean'])
                           for a, b in zip(primary, run['predictions'], strict=True)]
            refinements[run['name']] = {
                'maximum_spin_mean_change': max(differences),
                'maximum_binary_probability_change': max(differences)/2,
                'spin_mean_rms_change': float(np.sqrt(np.mean(np.asarray(differences)**2))),
                'rigorous_error_bound': False}
    zeros = load_calibration()['position']
    zero_rows = prediction['zero_duration_position_predictions']
    zero_errors = [row['predicted_mean']-float(zeros['spin_imag'][0, i])
                   for i, row in enumerate(zero_rows)]
    result = {key: value for key, value in prediction.items() if key != 'runs'}
    result.update({'schema': 'oph.anharmonic_observations.results.v1',
                   'prediction_payload_sha256_before_scoring': prediction_hash,
                   'runs': runs, 'numerical_refinements': refinements,
                   'independent_zero_duration_position_check': {
                       'scalar_means': len(zero_rows),
                       'spin_mean_rmse': float(np.sqrt(np.mean(np.asarray(zero_errors)**2))),
                       'used_in_calibration': False},
                   'checks': {'all_707_later_statistics_retained': all(
                       run['comparison']['combined']['scalar_means'] == 707 for run in runs),
                       'no_positive_duration_fit': not prediction['positive_duration_values_used']},
                   'all_pass_meaning': 'Execution and comparison contracts only; not empirical agreement or significance.'})
    result['all_pass'] = all(result['checks'].values())
    if not result['all_pass']:
        raise ArithmeticError('anharmonic comparison contract failed')
    return result
