"""Predict reserved oscillator readouts with native source operations, offline."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from ...observer.scalar_response.experiment import canonical, compare_receipt
from .data import FOLDER, load_observations
from .executable_bounds import cutoff_certificate
from .executable_calibration import calibrate_state
from .executable_source import source_run


def _json_value(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def predict(axis_means, *, extent, dimension):
    """Produce all off-axis predictions without accepting their observations."""
    calibration = calibrate_state(axis_means, extent=extent)
    parameters = calibration['parameters']
    alpha, thermal = parameters['alpha_real'], parameters['thermal_occupation']
    indices = [(i, j) for i in range(1, 11) for j in range(1, 11)]
    settings = [{'u': extent*i/10, 'v': extent*j/10, 'theta': theta}
                for i, j in indices for theta in (0., -math.pi/2)]
    execution = source_run(dimension, alpha, thermal, settings)
    if len(execution['records']) != len(settings):
        raise ArithmeticError('Source execution did not return every requested readout')
    for setting, record in zip(settings, execution['records'], strict=True):
        if any(record.get(key) != value for key, value in setting.items()):
            raise ArithmeticError('Source readout settings or order changed')
        probabilities = record['probabilities']
        if (len(probabilities) != 2 or
                not all(math.isfinite(value) and -2e-10 <= value <= 1+2e-10
                        for value in probabilities) or
                abs(sum(probabilities)-1) > 2e-10 or
                not math.isfinite(record['mean']) or
                abs(record['mean']-(probabilities[0]-probabilities[1])) > 2e-10):
            raise ArithmeticError('Invalid source binary readout')
    predictions = []
    for index, real, imag in zip(indices, execution['records'][::2],
                                 execution['records'][1::2], strict=True):
        u, v = real['u'], real['v']
        certificate = cutoff_certificate(
            dimension, thermal=F(thermal), alpha_real=F(alpha),
            beta_real=F(u), beta_imag=F(v))
        predictions.append({
            'index': list(index), 'u': u, 'v': v,
            'source_real': real['mean'], 'source_imag': imag['mean'],
            'source_real_probabilities': real['probabilities'],
            'source_imag_probabilities': imag['probabilities'],
            'occupation_certificate': _json_value(certificate),
        })
    return {'coordinate_extent': extent, 'calibration': calibration,
            'source_execution': execution, 'predictions': predictions}


def score(prediction, observed):
    """Descriptive comparison; calibrated-parameter uncertainty is unbounded here."""
    expected = {(i, j) for i in range(1, 11) for j in range(1, 11)}
    if set(observed) != expected:
        raise ValueError('Exactly the 100 off-axis observations must be scored')
    indices = [tuple(point['index']) for point in prediction['predictions']]
    if len(indices) != len(expected) or set(indices) != expected:
        raise ValueError('Exactly one source prediction per off-axis point is required')
    rows, squared, zero_squared = [], [], []
    for point in prediction['predictions']:
        key = tuple(point['index'])
        value = complex(observed[key])
        if not all(math.isfinite(x) and abs(x) <= 1 for x in (value.real, value.imag)):
            raise ValueError('Observed spin components must be finite and in [-1,1]')
        if not all(math.isfinite(point[key]) and abs(point[key]) <= 1+2e-10
                   for key in ('source_real', 'source_imag')):
            raise ValueError('Source spin components must be finite and in [-1,1] up to rounding')
        errors = [point['source_real']-value.real, point['source_imag']-value.imag]
        squared.extend(error**2 for error in errors)
        zero_squared.extend([value.real**2, value.imag**2])
        rows.append({'index': point['index'], 'observed_real': value.real,
                     'observed_imag': value.imag,
                     'observed_real_probability_zero': (1+value.real)/2,
                     'observed_imag_probability_zero': (1+value.imag)/2,
                     'real_residual': errors[0], 'imag_residual': errors[1]})
    rmse = math.sqrt(sum(squared)/len(squared))
    return {'withheld_complex_points': len(rows), 'withheld_real_statistics': len(squared),
            'spin_mean_rmse': rmse, 'binary_probability_rmse': rmse/2,
            'zero_response_spin_mean_rmse': math.sqrt(sum(zero_squared)/len(zero_squared)),
            'maximum_component_residual': math.sqrt(max(squared)),
            'parameter_uncertainty_propagated': False,
            'sampling_or_apparatus_error_certified': False,
            'statistical_significance_claimed': False, 'rows': rows}


def _audit(prediction):
    """Independent closed-form reference check, after source predictions exist."""
    parameters = prediction['calibration']['parameters']
    alpha, thermal = parameters['alpha_real'], parameters['thermal_occupation']
    errors, bounds, valid = [], [], []
    for point in prediction['predictions']:
        u, v = point['u'], point['v']
        amplitude = math.exp(-(thermal+.5)*(u*u+v*v))
        reference = complex(amplitude*math.cos(2*alpha*v), amplitude*math.sin(2*alpha*v))
        error = max(abs(point['source_real']-reference.real),
                    abs(point['source_imag']-reference.imag))/2
        bound = F(point['occupation_certificate']['probability_error_upper'])
        errors.append(error)
        bounds.append(bound)
        # This is a numerical consistency check, not the rational proof itself.
        valid.append(error <= float(bound)+2e-10)
    execution = prediction['source_execution']
    operations = [execution['preparation']['displacement'],
                  *execution['shared_detector_operations'].values()]
    for record in execution['records']:
        operations.extend([record['rotation'], record['conditional_detector_phase']])
    return {
        'reference_used_to_generate_predictions': False,
        'maximum_source_vs_unrestricted_probability_difference': max(errors),
        'maximum_certified_ideal_probability_error_exact': str(max(bounds)),
        'maximum_certified_ideal_probability_error': float(max(bounds)),
        'cutoff_consistent_with_independent_numeric_reference': all(valid),
        'maximum_full_phase_operation_residual': max(
            operation['full_phase_target_residual'] for operation in operations),
        'maximum_probability_mass_residual': max(
            abs(record['probability_mass']-1) for record in execution['records']),
        'numeric_reference_check_allowance': 2e-10,
        'numerical_or_physical_errors_in_rational_cutoff_bound': False,
    }


def run(*, progress=None):
    protocol_bytes = (FOLDER/'executable_protocol.json').read_bytes()
    protocol = json.loads(protocol_bytes)
    observations = load_observations()
    runs = []
    for extent in protocol['coordinate_extents']:
        if progress:
            progress(f'Executing dimension {protocol["occupation_dimension"]}, coordinate extent {extent}', flush=True)
        prediction = predict(observations['calibration'], extent=extent,
                             dimension=protocol['occupation_dimension'])
        prediction['comparison'] = score(prediction, observations['validation'])
        prediction['numerical_audit'] = _audit(prediction)
        runs.append(prediction)
    checks = {
        'axis_and_validation_sets_disjoint': set(observations['calibration']).isdisjoint(observations['validation']),
        'both_published_coordinate_extents_retained': [item['coordinate_extent'] for item in runs] == [4, 5],
        'all_off_axis_points_and_both_outcomes_retained': all(
            len(item['predictions']) == 100 and
            all(record['both_detector_outcomes_retained'] for record in item['source_execution']['records'])
            for item in runs),
        'source_probability_mass_preserved': all(
            item['numerical_audit']['maximum_probability_mass_residual'] < 2e-10 for item in runs),
        'source_agrees_with_independent_reference_within_cutoff_and_numeric_allowance': all(
            item['numerical_audit']['cutoff_consistent_with_independent_numeric_reference'] for item in runs),
        'predictions_computed_from_source_branches': all(
            not item['source_execution']['analytic_characteristic_used_for_prediction'] for item in runs),
    }
    bindings = ['executable.py', 'executable_source.py', 'canonical_source.py', 'executable_bounds.py',
                'executable_calibration.py', 'data.py', '../source.py',
                '../../observer/source_scalar/compiler.py']
    receipt = {
        'schema': 'oph.scalar_observations.executable_experimental_readout.v1',
        'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
        'implementation_sha256': {name: hashlib.sha256((FOLDER/name).read_bytes()).hexdigest()
                                  for name in bindings},
        'selection_note': protocol['selection_order'], 'scope': protocol['scope'],
        'provenance': observations['provenance'], 'orientation': observations['orientation'],
        'runs': runs, 'checks': checks, 'all_pass': all(checks.values()),
        'all_pass_meaning': 'Source execution and comparison contracts only; not a scientific significance verdict.',
    }
    if not receipt['all_pass']:
        raise ArithmeticError('Executable comparison contract failed')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run(progress=print)
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': result['all_pass'], 'runs': [
        {'extent': item['coordinate_extent'],
         'parameters': item['calibration']['parameters'],
         'spin_mean_rmse': item['comparison']['spin_mean_rmse'],
         'binary_probability_rmse': item['comparison']['binary_probability_rmse'],
         'occupation_probability_error_upper': item['numerical_audit']['maximum_certified_ideal_probability_error']}
        for item in result['runs']]}, indent=2))


if __name__ == '__main__':
    main()
