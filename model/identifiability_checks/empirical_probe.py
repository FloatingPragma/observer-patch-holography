"""Retrospective empirical comparison; no fitting and no replacement predictions.

Run from the repository root with single-threaded BLAS. This probe reads the
already scored observations. It evaluates a deliberately simpler, supplied
cosine Hamiltonian by an independent unsplit Lindblad integration, retaining
all observations. No result here is a new excluded-data test.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import zlib

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, expm
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

FOLDER = Path(__file__).resolve().parent
MODEL = FOLDER.parent
FROZEN_TRANSFER = FOLDER/'data/transfer_reference.json.gz'
FROZEN_TRANSFER_SHA256 = '71403f474c39705d20a6c30dc2f36b7556943dd4f1734a5a1a395184744d2398'


def frozen_transfer_bytes(path):
    """Read the fixed reference input independently of live replay."""
    compressed = path.read_bytes()
    try:
        raw = gzip.decompress(compressed)
    except (OSError, EOFError, zlib.error) as error:
        raise ValueError('Frozen transfer receipt has invalid gzip data') from error
    if hashlib.sha256(raw).hexdigest() != FROZEN_TRANSFER_SHA256:
        raise ValueError('Frozen transfer receipt differs from the authenticated reference')
    return raw


def target_law(dimension, alpha, thermal, phase, times):
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    number = np.arange(dimension, dtype=float)
    coordinate = a + a.T
    eigenvalues, eigenvectors = eigh(coordinate)
    cosine = (eigenvectors * np.cos(np.pi / 3 * eigenvalues + phase)) @ eigenvectors.T
    hamiltonian = 2 * np.pi * 500 * np.diag(number) + 4000 * cosine
    ratio = thermal / (1 + thermal)
    weights = (1 - ratio) * ratio**number
    weights[-1] = ratio**(dimension - 1)
    displacement = expm(alpha * (a.T - a))
    initial = (displacement * weights) @ displacement.conj().T
    difference_squared = (number[:, None] - number[None, :])**2

    def rhs(_time, vector):
        rho = vector.reshape((dimension, dimension))
        return (-1j * (hamiltonian @ rho - rho @ hamiltonian)
                - 9 * difference_squared * rho).ravel()

    solution = solve_ivp(rhs, (0., float(times[-1])), initial.astype(complex).ravel(),
                         t_eval=times, method='DOP853', rtol=2e-10, atol=2e-12)
    if not solution.success:
        raise ArithmeticError(solution.message)
    detector = expm(.4j * coordinate)
    states = solution.y.T.reshape((-1, dimension, dimension))
    predictions = [float(np.trace(detector @ rho).imag / (np.sqrt(2) * .4))
                   for rho in states]
    return predictions, max(float(abs(np.trace(rho) - 1)) for rho in states)


def run():
    paths = {
        'transfer_receipt': FROZEN_TRANSFER,
        'transfer_data': MODEL / 'anharmonic_transfer/data/asymmetric_wells_expectation_x.csv',
        'audit_probe': Path(__file__).resolve(),
    }
    transfer_bytes = frozen_transfer_bytes(paths['transfer_receipt'])
    receipt = json.loads(transfer_bytes)
    data = np.genfromtxt(paths['transfer_data'], delimiter=',', skip_header=1)
    times, observations = data[:, 0] / 1000, data[:, 1:]
    parameters = receipt['prediction']['prior_parameters']
    cases = receipt['prediction']['protocol']['cases']
    refinements = {}
    for dimension in (32, 48):
        predictions = []
        trace_errors = []
        for case in cases:
            values, trace_error = target_law(
                dimension, case['preparation_sign'] * abs(parameters['alpha_real']),
                parameters['thermal_occupation'], np.pi * case['potential_phase_pi'], times)
            predictions.append(values)
            trace_errors.append(trace_error)
        array = np.asarray(predictions).T
        refinements[str(dimension)] = {
            'dimension': dimension,
            'predictions_in_csv_order': array.tolist(),
            'rmse_positive_by_case': dict(zip('abcde', np.sqrt(np.mean(
                (array[1:] - observations[1:])**2, axis=0)).tolist(), strict=True)),
            'rmse_positive_combined': float(np.sqrt(np.mean((array[1:] - observations[1:])**2))),
            'maximum_trace_residual': max(trace_errors),
        }
    quantized_residual = 63.6 * observations - np.round(63.6 * observations)
    primary = receipt['prediction']['runs'][0]
    zero_errors = {
        row['case']: row['predicted_position_estimator'] - observations[0, index]
        for index, row in enumerate(primary['zero_time_predictions'])
    }
    return {
        'schema': 'oph.empirical_audit.retrospective.v1',
        'status': 'Already scored outcomes inspected; diagnostic comparison only. No fits.',
        'input_sha256': bound_pins({name: hashlib.sha256(
            transfer_bytes if name == 'transfer_receipt' else path.read_bytes()).hexdigest()
                         for name, path in paths.items()}),
        'fixed_target_model': {
            'hamiltonian': 'delta*N+4000*cos((pi/3)*(a+a_dagger)+potential_phase)',
            'delta_per_second': float(2 * np.pi * 500),
            'lindblad': 'sqrt(18)*N',
            'time_parameter': 'published simulated duration',
            'detector': 'Im Tr[rho*D(0.4i)]/(sqrt(2)*0.4)',
            'spin': 'No ancilla; this is the target Hamiltonian, not the physical pulse sequence.',
            'preparation': parameters,
            'integration': 'DOP853 rtol=2e-10 atol=2e-12; no trace/positivity repair',
        },
        'simpler_target_model': refinements,
        'target_cutoff_maximum_estimator_difference': float(np.max(np.abs(
            np.asarray(refinements['32']['predictions_in_csv_order'])
            - np.asarray(refinements['48']['predictions_in_csv_order'])))),
        'zero_prediction_rmse_positive': float(np.sqrt(np.mean(observations[1:]**2))),
        'quantization_rows': int(observations.size),
        'maximum_distance_of_63point6_times_data_from_integer': float(np.max(abs(quantized_residual))),
        'primary_zero_time_estimator_residuals': zero_errors,
        'original_primary_rmse': receipt['scores']['primary']['metrics']['combined']['position_estimator_rmse'],
        'normalization_inferred_or_corrected': False,
    }


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True, allow_nan=False))
