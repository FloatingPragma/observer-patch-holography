"""Attach retained unrestricted-apparatus intervals to unchanged observations.

This uses the independent occupation evaluator, not the native finite-source
prediction as an exact centre. It supplies no experimental error distribution.
"""
from __future__ import annotations

from decimal import Decimal, InvalidOperation
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = 'model/occupation_control/results.json'
CERTIFICATE_PINS = {
    'model/occupation_control/DERIVATION.md': '8e58825368295e9fa815bea4742bf52c4dac33c664961bf9098deed8b59eb5e6',
    'model/occupation_control/__main__.py': 'b442422dbc00e0873a934c8c9348ea26af75e851c36721784744dd9f88887cfc',
    'model/occupation_control/experiment.py': '4073a2d90e53004a0eb0c0a1fc185e4c94768798b07d3a934b161dd896920285',
    'model/occupation_control/intervals.py': '7c11fbde0360e2153fa8d3c154b128a0f238a85b2da1b345ea769b5df4c66bbc',
    CERTIFICATE: 'c07df1b9984e2d24a07916712880b09b3e567e60bc6f1f73d52dd07326557316',
    'model/occupation_control/solver.py': '00df02931eea0513303b070186c90a746931561326839bee8b23eb650f454c53',
}


def decimal_token(value):
    """Preserve a deposited finite decimal, rather than its binary conversion."""
    if type(value) is not str:
        raise ValueError('An original decimal string is required')
    try:
        number = Decimal(value)
    except InvalidOperation as error:
        raise ValueError('A finite decimal string is required') from error
    if not number.is_finite():
        raise ValueError('A finite decimal string is required')
    return F(number)


def unavoidable_allowance(observed, centre, radius):
    """Exact distance from a measured central value to a certified interval.

    If |truth-centre| <= radius and |observed-truth| <= remaining, then
    remaining >= this result. No stochastic interpretation is supplied.
    """
    if any(type(value) is not F for value in (observed, centre, radius)) or radius < 0:
        raise ValueError('Exact rationals and a nonnegative certified radius required')
    return max(F(0), abs(observed-centre)-radius)


def _binary(value, label, *, nonnegative=False):
    if type(value) is not float or not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f'Finite original binary value required for {label}')
    return F(value)


def certified_inputs(protocol, prior, apparatus, *, root=ROOT):
    """Authenticate and seal both retained refinements before later data open."""
    root = Path(root)
    for name, expected in CERTIFICATE_PINS.items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Occupation certificate source or artifact changed: {name}')
    retained = strict_json.load(root/CERTIFICATE)
    expected_parameters = {
        'alpha_real': -abs(prior['alpha_real']), 'thermal': prior['thermal_occupation'],
        'delta': apparatus['delta'], 'omega': apparatus['omega'],
        'pulse_duration': apparatus['pulse_duration'], 'theta': apparatus['theta'],
        'gamma': apparatus['dephasing_rate'], 'omega_carrier': apparatus['omega_carrier'],
        'probe_h': protocol['readout']['h'], 'phi': 0., 'gap_duration': 0.,
    }
    if (retained['schema'] != 'oph.occupation_control.results.v1'
            or [run['dimension'] for run in retained['runs']] != [32, 64]):
        raise ValueError('Both declared occupation refinements are required')
    runs = []
    for run in retained['runs']:
        if (run['case'] != 'a' or run['parameters'].keys() != expected_parameters.keys()
                or any(type(value) is not float or value != expected_parameters[key]
                       for key, value in run['parameters'].items())):
            raise ValueError('Occupation certificate and transferred apparatus parameters differ')
        flags = run['certificate']
        for key, expected in {
            'state_error_includes_occupation': True,
            'state_error_includes_time_discretization': True,
            'state_error_includes_evaluator_roundoff': True,
            'preparation_and_detector_displacements_included': True,
            'native_source_words_executed': False, 'native_compiler_error_included': False,
            'apparatus_or_calibration_uncertainty_included': False,
            'empirical_outcomes_loaded': False,
        }.items():
            if flags.get(key) is not expected:
                raise ValueError('Occupation certificate scope changed')
        if [row['steps'] for row in run['records']] != [0, 2, 10, 25, 50]:
            raise ValueError('Complete retained occupation endpoint inventory required')
        records = []
        for row in run['records']:
            step = row['steps']
            if (type(step) is not int or row['published_coordinate_ms'] != step/5
                    or row['drive_elapsed_seconds'] != step*4*expected_parameters['pulse_duration']):
                raise ValueError('Occupation certificate uses another endpoint or clock')
            centre = _binary(row['selected_estimator'], 'estimator')
            radius = row['estimator_error_upper']
            if (type(row['acceptance_positive_certified']) is not bool
                    or row['acceptance_positive_certified'] != (radius is not None)):
                raise ValueError('Conditional certificate availability changed')
            radius = None if radius is None else _binary(radius, 'radius', nonnegative=True)
            records.append({'case': 'a', 'time_ms': str(F(step, 5)),
                            'centre_exact_binary': str(centre),
                            'radius_exact_binary': None if radius is None else str(radius),
                            'conditional_certificate_available': radius is not None})
        runs.append({'dimension': run['dimension'], 'records': records})
    return {'certificate': CERTIFICATE, 'source_pins': dict(CERTIFICATE_PINS),
            'source_and_artifact_hashes_verified': True, 'certificate_recomputed_here': False,
            'parameters': expected_parameters, 'runs': runs,
            'target': 'fixed unrestricted spin-oscillator apparatus; independent semantic evaluator',
            'native_compiler_or_physical_service_error_included': False,
            'calibration_and_experimental_uncertainty_included': False,
            'new_parameter_fit': False}


def compare(attachment, observations, zero_tokens):
    """Retain every held-out row and explain every unavailable comparison."""
    expected = {(case, F(step, 5)) for case in 'abcde' for step in range(2, 51, 2)}
    observed = {}
    for row in observations:
        if (type(row['case']) is not str or type(row['time_ms']) is not float
                or not math.isfinite(row['time_ms'])
                or type(row['observed_position_estimator']) is not float
                or not math.isfinite(row['observed_position_estimator'])):
            raise ValueError('Original finite observation scalars required')
        key = row['case'], F(str(row['time_ms']))
        if key in observed or key not in expected:
            raise ValueError('Duplicate or foreign held-out observation')
        token = row['observed_position_estimator_decimal']
        value = decimal_token(token)
        if float(value) != row['observed_position_estimator']:
            raise ValueError('Observed decimal token and scored floating value differ')
        observed[key] = token, value
    if set(observed) != expected or set(zero_tokens) != set('abcde'):
        raise ValueError('All 125 held-out observations and five separate zero-time values required')
    zero_values = {case: decimal_token(token) for case, token in zero_tokens.items()}

    def scored(record, token, value):
        centre = F(record['centre_exact_binary'])
        radius = F(record['radius_exact_binary'])
        required = unavoidable_allowance(value, centre, radius)
        return {'observed_decimal': token, 'certified_interval': [str(centre-radius), str(centre+radius)],
                'minimum_uncovered_absolute_allowance': str(required),
                'zero_additional_allowance_compatible': required == 0}

    runs = []
    for run in attachment['runs']:
        by_key = {(r['case'], F(r['time_ms'])): r for r in run['records']}
        if len(by_key) != len(run['records']):
            raise ValueError('Duplicate certified coordinate')
        rows, floors = [], []
        for key in sorted(expected):
            token, value = observed[key]
            row = {'case': key[0], 'time_ms': str(key[1])}
            record = by_key.get(key)
            if record is None:
                row['status'] = 'no_retained_occupation_endpoint'
            elif not record['conditional_certificate_available']:
                row['status'] = 'reference_acceptance_not_certified_positive'
            else:
                row.update(status='certified_interval', **scored(record, token, value))
                floors.append(F(row['minimum_uncovered_absolute_allowance']))
            rows.append(row)
        zero, unmatched = [], []
        for key, record in by_key.items():
            if key[1] == 0:
                zero.append({'case': key[0], 'time_ms': '0',
                             **scored(record, zero_tokens[key[0]], zero_values[key[0]])})
            elif key not in observed:
                unmatched.append({'case': key[0], 'time_ms': str(key[1]),
                                  'reason': 'no observation at this exact published coordinate; no interpolation'})
        runs.append({'dimension': run['dimension'], 'held_out_rows': rows,
                     'held_out_count': len(rows), 'certified_overlap_count': len(floors),
                     'unavailable_count': len(rows)-len(floors),
                     'minimum_common_uncovered_allowance_on_certified_held_out_rows': str(max(floors, default=F(0))),
                     'zero_time_diagnostics_not_held_out': zero,
                     'certificate_endpoints_without_observations': unmatched})
    return {'runs': runs,
            'relation': 'remaining combined experimental/calibration/model discrepancy >= distance(observed, certified interval)',
            'units': 'dimensionless two-point position estimator',
            'experimental_uncertainty_applied': None, 'statistical_rejection_claimed': False,
            'full_uncertainty_agreement_established': False,
            'published_measurements_per_position_point': 200,
            'accepted_counts_or_raw_trials_authenticated': False,
            'source': 'https://arxiv.org/html/2603.04744v1#S3.SS2',
            'measurement_boundary': 'Paper reports 200 measurements per point and h=0.4; the deposit omits accepted counts, raw trials and calibration covariance. The unresolved estimator normalization is not corrected.',
            'overlap_boundary': 'These are existing observations, not additional independent evidence. Zero-time values are separate, and neither missing endpoints nor potential cases are interpolated.'}
