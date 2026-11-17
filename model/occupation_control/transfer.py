"""Certify every deposited changed-potential endpoint without fitting outcomes.

The unchanged independent evaluator encloses the unrestricted apparatus.
Predictions are sealed before measurements open; this ordering is an audit
boundary within a retrospective test, not prospective or blinded evidence.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from model.anharmonic_transfer import certified_comparison as comparison, data
from model.anharmonic_transfer.calibration import transferred_parameters
from model.anharmonic_transfer.experiment import _configuration
from . import intervals as iv
from .__main__ import compare as numerical_compare
from .experiment import PARAMETERS, CASES, run

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
STEPS = tuple(range(0, 51, 2))
CASES_ORDER = tuple('abcde')
SOURCE_PATHS = tuple(sorted({
    *comparison.CERTIFICATE_PINS,
    'model/occupation_control/transfer.py',
    'model/anharmonic_transfer/calibration.py',
    'model/anharmonic_transfer/data.py',
    'model/anharmonic_transfer/certified_comparison.py',
    'model/anharmonic_transfer/experiment.py',
    'model/anharmonic_transfer/protocol.json',
    'model/anharmonic_transfer/data/calibration_receipt.json.gz',
    'tools/strict_json.py',
}))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def settings():
    """Verify identity with the frozen supplied apparatus and calibration."""
    protocol = strict_json.load(ROOT/'model/anharmonic_transfer/protocol.json')
    apparatus = _configuration(protocol)
    prior = transferred_parameters(protocol)
    # Also authenticates the original proof sources and calibration contract.
    comparison.certified_inputs(protocol, prior, apparatus)
    expected = {
        'alpha_real': -abs(prior['alpha_real']), 'thermal': prior['thermal_occupation'],
        'delta': apparatus['delta'], 'omega': apparatus['omega'],
        'pulse_duration': apparatus['pulse_duration'], 'theta': apparatus['theta'],
        'gamma': apparatus['dephasing_rate'], 'omega_carrier': apparatus['omega_carrier'],
        'probe_h': protocol['readout']['h'],
    }
    if PARAMETERS != expected or any(type(v) is not float for v in PARAMETERS.values()):
        raise ValueError('The fixed apparatus or preparation parameters changed')
    for item in protocol['cases']:
        phase, side = CASES[item['label']]
        # Match the original executor's arithmetic, including binary pi/20.
        expected_phase = {'a': 0., 'b': math.pi/20, 'c': math.pi/20,
                          'd': math.pi/20, 'e': math.pi/10}[item['label']]
        if phase != expected_phase or {'left': -1, 'right': 1, 'center': 0}[side] != item['preparation_sign']:
            raise ValueError('The potential/preparation case mapping changed')
    return {'dimension': 64, 'steps': 50, 'pulse_cells': 5, 'order': 12,
            'record_steps': list(STEPS), 'cases': list(CASES_ORDER),
            'parameters': expected, 'new_parameter_fit': False,
            'target': 'unchanged unrestricted supplied spin-oscillator apparatus',
            'selection': 'all five cases and every deposited coordinate; no outcome-dependent selection',
            'classification': 'retrospective supplied-apparatus postdiction; not an OPH-specific prediction',
            'physical_uncertainty_included': False, 'native_compiler_error_included': False}


def conditional_error(acceptance, half_trace_error, trace_defect):
    """Exact trace-aware conditioning bound for a Hermitian approximate path.

    For H=A-pQ, 0<=A<=Q<=I and physical p in [0,1], the spectral interval
    [-p,1-p] has width one. Subtracting its midpoint gives
    |Tr(H Delta)| <= ||Delta||_1/2 + |Tr Delta|/2.
    No positivity of the approximate state or independence is assumed.
    """
    if (any(type(x) is not F for x in (acceptance, half_trace_error, trace_defect))
            or acceptance <= 0 or half_trace_error < 0 or trace_defect < 0):
        raise ValueError('Positive exact acceptance and nonnegative rational errors required')
    effect_error = half_trace_error+trace_defect/2
    return None if acceptance <= effect_error else effect_error/acceptance


def trace_aware_readout(row):
    """Enclose the full-readout trace defect and final scalar conversions."""
    drift = row['mass_drift']
    if type(drift) is not float or not math.isfinite(drift) or drift < 0:
        raise ValueError('Original finite trace-defect scalar required')
    # detector() rounds the exact sum of binary diagonal entries once.
    # Zero is exact: any nonzero binary sum is at least the least subnormal.
    drift_upper = 0. if drift == 0. else math.nextafter(drift, math.inf)
    error = conditional_error(F(row['acceptance_exact_binary_sum']),
                              F(row['half_trace_error_upper']), F(drift_upper))
    probability = mean = estimator = None
    if error is not None:
        lower = iv.mul(iv.square_root(iv.exact(2.)), iv.exact(PARAMETERS['probe_h']))[0]
        # No probability-radius cap: the theorem permits signed approximate
        # states, and the separately rounded readout need not lie in [0,1].
        probability = iv.rational(error+F(row['selected_probability_evaluation_error_upper']))[1]
        mean = iv.rational(2*error+F(row['selected_mean_evaluation_error_upper']))[1]
        estimator = iv.rational(2*error/F(float(lower))+F(row['estimator_evaluation_error_upper']))[1]
    return {'trace_defect_upper': drift_upper,
            'acceptance_positive_certified': error is not None,
            'selected_probability_error_upper': probability,
            'selected_mean_error_upper': mean, 'estimator_error_upper': estimator}


def seal_runs(runs):
    """Bind independently computed paths before opening comparison values."""
    runs = deepcopy(runs)
    for result in runs:
        for row in result['records']:
            row['trace_aware_readout'] = trace_aware_readout(row)
    prediction = {'schema': 'oph.occupation_control.transfer_prediction.v1',
                  'settings': settings(),
                  'source_bindings': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                      for name in SOURCE_PATHS},
                  'runs': runs}
    prediction['seal'] = digest(prediction)
    validate_prediction(prediction)
    return prediction


def predict(progress=None):
    configuration = settings()
    runs = []
    for case in CASES_ORDER:
        if progress:
            progress(case)
        runs.append(run(dimension=configuration['dimension'], case=case,
                        steps=configuration['steps'], pulse_cells=configuration['pulse_cells'],
                        order=configuration['order'], record_steps=configuration['record_steps']))
    return seal_runs(runs)


def validate_prediction(prediction):
    """Check scope, exact endpoint coverage and own-path radius accounting.

    This validates a receipt's structure, not its differential equation.
    Full replay recomputes every path and enclosure through ``predict``.
    """
    numerical_compare(prediction, prediction)  # Reject all nonfinite floats.
    if (set(prediction) != {'schema', 'settings', 'source_bindings', 'runs', 'seal'}
            or prediction.get('schema') != 'oph.occupation_control.transfer_prediction.v1'):
        raise ValueError('Transfer prediction schema changed')
    unsigned = {k: v for k, v in prediction.items() if k != 'seal'}
    if prediction.get('seal') != digest(unsigned):
        raise ValueError('Prediction seal changed before comparison')
    if digest(prediction['settings']) != digest(settings()):
        raise ValueError('Fixed transfer settings changed')
    bindings = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS}
    if prediction['source_bindings'] != bindings:
        raise ValueError('Certificate source bindings changed')
    if [r['case'] for r in prediction['runs']] != list(CASES_ORDER):
        raise ValueError('All five cases in their fixed order are required')
    for result in prediction['runs']:
        case = result['case']
        for key, expected in [('dimension', 64), ('steps', 50), ('pulse_cells', 5), ('order', 12)]:
            if type(result[key]) is not int or result[key] != expected:
                raise ValueError('Fixed numerical menu changed')
        phase, side = CASES[case]
        alpha = {'left': PARAMETERS['alpha_real'], 'right': -PARAMETERS['alpha_real'], 'center': 0.}[side]
        expected = {**PARAMETERS, 'phi': phase, 'alpha_real': alpha,
                    'gap_duration': abs(phase)/PARAMETERS['omega_carrier']}
        if result['parameters'] != expected or any(type(v) is not float for v in result['parameters'].values()):
            raise ValueError('Case parameters differ from the same transferred apparatus')
        if [r['steps'] for r in result['records']] != list(STEPS):
            raise ValueError('Every deposited endpoint is required, including zero time')
        flags = result['certificate']
        for key, value in {
            'state_error_includes_occupation': True, 'state_error_includes_time_discretization': True,
            'state_error_includes_evaluator_roundoff': True, 'preparation_and_detector_displacements_included': True,
            'native_source_words_executed': False, 'native_compiler_error_included': False,
            'apparatus_or_calibration_uncertainty_included': False, 'empirical_outcomes_loaded': False,
        }.items():
            if flags.get(key) is not value:
                raise ValueError('Certificate scope changed')
        for row in result['records']:
            step = row['steps']
            if (type(step) is not int or row['published_coordinate_ms'] != step/5
                    or row['drive_elapsed_seconds'] != step*4*(PARAMETERS['pulse_duration']+expected['gap_duration'])):
                raise ValueError('Published coordinate or drive clock changed')
            for key in ('selected_estimator', 'half_trace_error_upper', 'estimator_evaluation_error_upper'):
                if type(row[key]) is not float:
                    raise ValueError('Original binary readout and error scalars required')
            if row['estimator_error_upper'] is not None and type(row['estimator_error_upper']) is not float:
                raise ValueError('Original binary selected radius required')
            budget = row['budget']
            if any(type(v) is not float or v < 0 for v in budget.values()):
                raise ValueError('Nonnegative original binary error budgets required')
            epsilon = iv.sum_up(list(budget.values()))
            if row['half_trace_error_upper'] != epsilon:
                raise ValueError('Half-trace radius differs from its outward budget')
            q = F(row['acceptance_exact_binary_sum'])
            available = q > 2*F(epsilon)
            if row['acceptance_positive_certified'] is not available:
                raise ValueError('Uncertified conditioning cannot be promoted')
            radius = None
            if available:
                lower = iv.mul(iv.square_root(iv.exact(2.)), iv.exact(PARAMETERS['probe_h']))[0]
                radius = iv.rational(4*F(epsilon)/q/F(float(lower))+
                                     F(row['estimator_evaluation_error_upper']))[1]
            if row['estimator_error_upper'] != radius:
                raise ValueError('Selected radius differs from its own full-path error')
            if digest(row['trace_aware_readout']) != digest(trace_aware_readout(row)):
                raise ValueError('Trace-aware radius differs from its full-path error and trace defect')


def attachment(prediction):
    validate_prediction(prediction)
    records = []
    for result in prediction['runs']:
        for row in result['records']:
            radius = row['trace_aware_readout']['estimator_error_upper']
            records.append({'case': result['case'], 'time_ms': str(F(row['steps'], 5)),
                            'centre_exact_binary': str(F(row['selected_estimator'])),
                            'radius_exact_binary': None if radius is None else str(F(radius)),
                            'conditional_certificate_available': radius is not None})
    return {'runs': [{'dimension': 64, 'records': records}]}


def evaluate(prediction):
    attached = attachment(prediction)  # Before any empirical loader.
    scored = comparison.compare(attached, data.load_validation(), data.load_zero_tokens())
    scored['overlap_boundary'] = ('Every deposited case/time is retained; zero-time diagnostics are separate. '
                                  'These are existing observations, not additional independent evidence.')
    rows = scored['runs'][0]['held_out_rows']
    summaries = []
    for case in CASES_ORDER:
        selected = [row for row in rows if row['case'] == case]
        available = [row for row in selected if row['status'] == 'certified_interval']
        floors = [F(row['minimum_uncovered_absolute_allowance']) for row in available]
        largest = max(available, key=lambda r: F(r['minimum_uncovered_absolute_allowance']), default=None)
        summaries.append({'case': case, 'observations': len(selected), 'certified': len(available),
                          'strictly_positive_discrepancy_floors': sum(value > 0 for value in floors),
                          'minimum_common_remaining_allowance': str(max(floors, default=F(0))),
                          'largest_floor_time_ms': None if largest is None else largest['time_ms']})
    scored['per_case'] = summaries
    return {'schema': 'oph.occupation_control.transfer_results.v1', 'prediction': prediction,
            'data_provenance': data.verify_provenance(), 'comparison': scored}


def verify_receipt(receipt):
    """Recompute exact comparisons from this receipt's own certified path."""
    numerical_compare(receipt, receipt)
    expected = evaluate(receipt['prediction'])
    if digest(receipt) != digest(expected):
        raise ValueError('Receipt comparisons or provenance differ from its own prediction')


def replay(actual, expected):
    # Exact fractions are checked against each receipt's own binary path.
    # Floating replay tolerance never becomes a physical uncertainty band.
    verify_receipt(actual)
    verify_receipt(expected)
    a, b = deepcopy(actual['prediction']), deepcopy(expected['prediction'])
    a.pop('seal'); b.pop('seal')
    numerical_compare(a, b)


def read_receipt(path):
    raw = Path(path).read_bytes()
    if Path(path).suffix == '.gz':
        raw = gzip.decompress(raw)
    receipt = strict_json.loads(raw.decode('utf-8'))
    verify_receipt(receipt)
    return receipt


def write_receipt(path, result):
    verify_receipt(result)
    raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    if Path(path).suffix == '.gz':
        raw = gzip.compress(raw, mtime=0)
    with Path(path).open('xb') as stream:
        stream.write(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = evaluate(predict(progress=lambda case: print(f'Certifying case {case}', flush=True)))
    if args.output:
        write_receipt(args.output, result)
    else:
        replay(result, read_receipt(args.check))
        print('All 125 fixed-apparatus comparisons and unrestricted bounds replayed.')


if __name__ == '__main__':
    main()
