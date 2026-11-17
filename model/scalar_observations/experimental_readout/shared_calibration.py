"""One axis-calibrated oscillator state against both reserved quadratures.

Exact interval cells cover the grid predictions of all displaced-thermal
states compatible with the axis confidence boxes, modulo phase aliases.
A cover intersection is not a feasibility certificate.
A separate, axis-only selected state must satisfy every box simultaneously
to supply a common witness. Neither calculation authenticates experimental
shots, preparation selection, or a physical service-error allowance.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import factorial, pi
from pathlib import Path

from tools import strict_json
from ...calibrated_measurement.bounds import log_interval, root_interval
from ...interacting_band_extension.thermal_spectrum import exponential_interval
from .data import FOLDER, load_observations
from .executable_calibration import AXIS_KEYS, _axes, calibrate_state
from .stats import concentration_certificate
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[3]
OUTPUT = FOLDER/'shared_calibration_results.json'
VALIDATION_KEYS = frozenset((i, j) for i in range(1, 11) for j in range(1, 11))
PI_UPPER = F(22, 7)
ROUNDING_GRID = 10**30
REPORTED_MEAN_ROUNDING = F(1, 10**12)
SOURCE_PATHS = (
    'tools/strict_json.py',
    'model/scalar_observations/experimental_readout/shared_calibration.py',
    'model/scalar_observations/experimental_readout/test_shared_calibration.py',
    'model/scalar_observations/experimental_readout/executable_calibration.py',
    'model/scalar_observations/experimental_readout/data.py',
    'model/scalar_observations/experimental_readout/stats.py',
    'model/scalar_observations/experimental_readout/executable_results.json',
    'model/scalar_observations/experimental_readout/provenance.json',
    'model/scalar_observations/experimental_readout/data/cf_0ms_real.csv',
    'model/scalar_observations/experimental_readout/data/cf_0ms_imag.csv',
    'model/scalar_observations/experimental_readout/data/SOURCE_README.md',
    'model/calibrated_measurement/bounds.py',
    'model/interacting_band_extension/thermal_spectrum.py',
)


def _canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def _hash(value):
    return hashlib.sha256(_canonical(value).encode()).hexdigest()


def _rounded(interval):
    lo, hi = interval
    return F((lo*ROUNDING_GRID).__floor__(), ROUNDING_GRID), F((hi*ROUNDING_GRID).__ceil__(), ROUNDING_GRID)


def _interval(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(x) is not F for x in value) or value[0] > value[1]):
        raise ValueError('ordered exact rational interval required')
    return value


def _product(left, right):
    values = [a*b for a in left for b in right]
    return min(values), max(values)


def sampling_contract(statistics, *, samples=250, delta=F(1, 40)):
    """Deterministic rational radius, with an exact check of its final float.

    Directed log/root series choose the radius; platform libm sqrt/log do
    not choose receipt bits. The final binary64 conversion is independently
    certified at its own exact rational value by the concentration bound.
    """
    if (type(statistics) is not int or statistics < 1 or type(samples) is not int or samples < 1
            or type(delta) is not F or not 0 < delta < 1):
        raise ValueError('positive fixed statistic/shot counts and exact failure allowance required')
    logarithm = log_interval(F(2*statistics)/delta)[1]
    rational = root_interval(2*logarithm/samples)[1]+F(1, 10**15)
    width = F(float(rational))
    certificate = concentration_certificate(statistics, samples, float(width))
    if F(certificate['familywise_error_upper']) > delta:
        raise ArithmeticError('final sampling radius missed its exact failure allowance')
    return width, certificate


@lru_cache(maxsize=8192)
def _trig_point(value, cosine):
    """Degree-128 Taylor enclosure with derivative bounded by one."""
    if type(value) is not F or abs(value) > 32 or type(cosine) is not bool:
        raise ValueError('exact phase in [-32,32] and Boolean parity required')
    term = total = F(1) if cosine else value
    power = 0 if cosine else 1
    while power+2 <= 128:
        term *= -value*value/((power+1)*(power+2))
        total += term
        power += 2
    remainder = abs(value)**129/F(factorial(129))
    return _rounded((total-remainder, total+remainder))


@lru_cache(maxsize=8192)
def _trig_range(interval, cosine):
    lo, hi = _interval(interval)
    middle, radius = (lo+hi)/2, (hi-lo)/2
    lower, upper = _trig_point(middle, cosine)
    # Both functions are 1-Lipschitz; no extrema or branch are missed.
    return max(F(-1), lower-radius), min(F(1), upper+radius)


@lru_cache(maxsize=32768)
def _amplitude(interval, radius_squared):
    lo, hi = _interval(interval)
    if lo < 0 or type(radius_squared) is not int or not 0 <= radius_squared <= 200:
        raise ValueError('nonnegative damping interval and declared grid radius required')
    return (exponential_interval(-hi*radius_squared, bits=128)[0],
            exponential_interval(-lo*radius_squared, bits=128)[1])


def characteristic_box(cell, key):
    """Enclose exp(-a*(i*i+j*j))*exp(i*theta*j) on one shared cell."""
    if (type(cell) is not tuple or len(cell) != 2
            or type(key) is not tuple or len(key) != 2
            or any(type(i) is not int or not 0 <= i <= 10 for i in key)):
        raise ValueError('one damping/phase cell and declared grid key required')
    damping, theta = map(_interval, cell)
    if damping[0] < 0 or theta[0] < -PI_UPPER or theta[1] > PI_UPPER:
        raise ValueError('nonnegative damping and principal phase cover required')
    i, j = key
    amplitude = _amplitude(damping, i*i+j*j)
    phase = (theta[0]*j, theta[1]*j)
    return {'real': _product(amplitude, _trig_range(phase, True)),
            'imag': _product(amplitude, _trig_range(phase, False))}


def observation_boxes(means, width, *, keys):
    if (type(means) is not dict or set(means) != keys
            or type(width) is not F or not 0 <= width <= 2):
        raise ValueError('exactly the declared observations and an exact width in [0,2] required')
    _axes(means) if keys == AXIS_KEYS else None
    rows = {}
    for key, value in means.items():
        if type(value) not in (float, complex, int) or type(value) is bool:
            raise ValueError('finite original numeric means required')
        value = complex(value)
        if not all(-1 <= x <= 1 for x in (value.real, value.imag)):
            raise ValueError('finite spin means in [-1,1] required')
        rows[key] = {name: (max(F(-1), F(x)-width), min(F(1), F(x)+width))
                     for name, x in (('real', value.real), ('imag', value.imag))}
    return rows


def _compatible(cell, boxes):
    for key, observed in boxes.items():
        predicted = characteristic_box(cell, key)
        if any(predicted[part][1] < interval[0] or interval[1] < predicted[part][0]
               for part, interval in observed.items()):
            return False
    return True


def axis_cover(axis_means, *, extent=4, width=F(1, 4), depth=10):
    """A finite outer cover, calibrated without accepting validation values.

    a=(nbar+1/2)*(extent/10)^2, theta=2*alpha*extent/10.
    Grid phases alias modulo 2*pi. [-22/7,22/7] covers every alias class.
    Real-axis amplitudes bound a above; no thermal cutoff is guessed.
    """
    if type(extent) is not int or extent not in (4, 5) or type(depth) is not int or not 0 <= depth <= 14:
        raise ValueError('published extent 4 or 5 and subdivision depth 0..14 required')
    boxes = observation_boxes(axis_means, width, keys=AXIS_KEYS)
    lower, upper = F(extent*extent, 200), None
    for i in range(1, 11):
        lo, hi = boxes[i, 0]['real']
        if hi <= 0:
            return ()
        lower = max(lower, log_interval(1/hi)[0]/(i*i))
        if lo > 0:
            bound = log_interval(1/lo)[1]/(i*i)
            upper = bound if upper is None else min(upper, bound)
    if upper is None:
        raise ValueError('axis boxes do not bound the thermal domain; conclusion is unresolved')
    if lower > upper:
        return ()
    lower, upper = _rounded((lower, upper))
    initial = ((lower, upper), (-PI_UPPER, PI_UPPER))
    active = [initial]
    for level in range(depth+1):
        survivors = [cell for cell in active if _compatible(cell, boxes)]
        if level == depth:
            return tuple(survivors)
        active = []
        axis = level % 2
        for cell in survivors:
            lo, hi = cell[axis]
            middle = (lo+hi)/2
            for half in ((lo, middle), (middle, hi)):
                child = list(cell)
                child[axis] = half
                active.append(tuple(child))
    raise AssertionError('unreachable subdivision')


def forward_envelopes(cells, keys=VALIDATION_KEYS):
    if not cells:
        return {}
    envelopes = {}
    for key in sorted(keys):
        predictions = [characteristic_box(cell, key) for cell in cells]
        envelopes[key] = {part: (min(row[part][0] for row in predictions),
                                max(row[part][1] for row in predictions))
                          for part in ('real', 'imag')}
    return envelopes


def common_witness(parameters, axis_means, validation_means, *, extent, axis_width, validation_width):
    """Check one axis-selected state on all observations; never optimize here."""
    if (type(parameters) is not dict or set(parameters) != {'alpha_real', 'thermal_occupation'}
            or any(type(x) is not float for x in parameters.values())):
        raise ValueError('original binary64 source preparation parameters required')
    alpha, occupation = (F(parameters[key]) for key in ('alpha_real', 'thermal_occupation'))
    if occupation < 0 or type(extent) is not int or extent not in (4, 5):
        raise ValueError('physical thermal occupation and published extent required')
    step = F(extent, 10)
    a, theta = (occupation+F(1, 2))*step*step, 2*alpha*step
    cell = ((a, a), (theta, theta))
    records, compatible = [], True
    for role, means, width, keys in (
            ('calibration', axis_means, axis_width, AXIS_KEYS),
            ('validation', validation_means, validation_width, VALIDATION_KEYS)):
        observed = observation_boxes(means, width, keys=keys)
        margins = []
        for key, boxes in observed.items():
            predicted = characteristic_box(cell, key)
            for part in ('real', 'imag'):
                lo, hi = predicted[part]
                margin = min(lo-boxes[part][0], boxes[part][1]-hi)
                margins.append(margin)
                compatible &= margin >= 0
        records.append({'role': role, 'real_statistics': 2*len(means),
                        'minimum_exact_box_margin': str(min(margins))})
    return {'one_state_used_for_every_point_and_quadrature': True,
            'parameters': parameters, 'source_parameters_treated_as_exact_binary64': True,
            'all_exact_prediction_intervals_inside_simultaneous_boxes': bool(compatible),
            'checks': records}


def source_witness_link(run, parameters, validation_means, width):
    """Check the already executed 64-level source uses this very same state.

    This binds an existing receipt, not a fresh source replay or a portable
    authentication of its backend residuals. Its normal owner replays it.
    """
    source = run['source_execution']
    if (run['calibration']['parameters'] != parameters or source['dimension'] != 64
            or source['analytic_characteristic_used_for_prediction'] is not False
            or source['preparation']['alpha_real'] != parameters['alpha_real']
            or source['preparation']['thermal']['mean_thermal_occupation_parameter'] != parameters['thermal_occupation']
            or len(source['records']) != 200
            or any(row['both_detector_outcomes_retained'] is not True or row['postselection'] is not False
                   for row in source['records'])
            or len(run['predictions']) != 100):
        raise ValueError('retained native source and shared state dictionary differ')
    boxes = observation_boxes(validation_means, width, keys=VALIDATION_KEYS)
    indices = [tuple(row['index']) for row in run['predictions']]
    if len(set(indices)) != 100 or set(indices) != VALIDATION_KEYS:
        raise ValueError('complete source readout menu required')
    records = {(row['u'], row['v'], row['theta']): row for row in source['records']}
    if len(records) != 200:
        raise ValueError('unique complete source detector settings required')
    margins = []
    for row in run['predictions']:
        i, j = row['index']
        u, v = float(F(run['coordinate_extent']*i, 10)), float(F(run['coordinate_extent']*j, 10))
        if (row['u'], row['v']) != (u, v):
            raise ValueError('source coordinates differ from the fixed grid dictionary')
        for part, phase in (('real', 0.), ('imag', -pi/2)):
            executed = records.get((u, v, phase))
            if executed is None or executed['mean'] != row['source_'+part]:
                raise ValueError('source quadratures differ from the fixed detector dictionary')
            mean = F(row['source_'+part])
            lo, hi = boxes[tuple(row['index'])][part]
            margins.append(min(mean-lo, hi-mean))
    return {'same_binary_preparation_and_fixed_detector_dictionary': True,
            'both_complete_detector_quadratures_retained': True,
            'all_retained_source_means_inside_validation_boxes': min(margins) >= 0,
            'minimum_exact_validation_margin_of_retained_source_means': str(min(margins)),
            'source_receipt': 'model/scalar_observations/experimental_readout/executable_results.json',
            'source_receipt_reexecuted_here': False,
            'source_numerical_residuals_authenticated_here': False,
            'physical_service_error_certified': False}


def build(*, depth=10):
    observations = load_observations()
    axis, validation = observations['calibration'], observations['validation']
    sampling = {role: sampling_contract(count)
                for role, count in (('calibration', 42), ('validation', 200))}
    widths = {role: row[0]+REPORTED_MEAN_ROUNDING for role, row in sampling.items()}
    simultaneous_failure = sum((F(row[1]['familywise_error_upper']) for row in sampling.values()), F())
    if simultaneous_failure > F(1, 20):
        raise ArithmeticError('combined simultaneous sampling allowance exceeds one twentieth')
    source_receipt = strict_json.loads((FOLDER/'executable_results.json').read_text(encoding='utf-8'))
    source_runs = {row['coordinate_extent']: row for row in source_receipt['runs']}
    if set(source_runs) != {4, 5} or len(source_receipt['runs']) != 2:
        raise ValueError('both published source coordinate conventions required')
    runs = []
    # Both predictions and chosen witnesses are sealed before scoring off-axis
    # observations. Validation never changes the state or the confidence cover.
    for extent in (4, 5):
        calibration = calibrate_state(axis, extent=extent)
        cells = axis_cover(axis, extent=extent, width=widths['calibration'], depth=depth)
        envelopes = forward_envelopes(cells)
        prediction = {'extent': extent, 'parameters': calibration['parameters'],
                      'cells': [[[str(x) for x in interval] for interval in cell] for cell in cells],
                      'envelopes': [{'index': list(key), **{part: [str(x) for x in interval]
                                     for part, interval in row.items()}}
                                    for key, row in envelopes.items()]}
        seal = _hash(prediction)
        heldout_boxes = observation_boxes(validation, widths['validation'], keys=VALIDATION_KEYS)
        remaining = tuple(cell for cell in cells if _compatible(cell, heldout_boxes))
        witness = common_witness(calibration['parameters'], axis, validation, extent=extent,
                                 axis_width=widths['calibration'], validation_width=widths['validation'])
        source_link = source_witness_link(source_runs[extent], calibration['parameters'],
                                          validation, widths['validation'])
        if _hash(prediction) != seal:
            raise ArithmeticError('axis-only prediction changed after opening validation')
        feasible = witness['all_exact_prediction_intervals_inside_simultaneous_boxes']
        if feasible and not remaining:
            raise ArithmeticError('outer cover lost its independently checked common witness')
        runs.append({'prediction': prediction, 'prediction_sha256_before_validation': seal,
                     'axis_cover_cells': len(cells), 'cells_not_excluded_by_all_validation_boxes': len(remaining),
                     'shared_state_witness': witness,
                     'same_state_native_source_link': source_link,
                     'verdict': ('one_axis_selected_state_is_simultaneously_compatible' if feasible else
                                 'shared_family_excluded_under_declared_contract' if not remaining else
                                 'outer_cover_intersection_is_inconclusive')})
    return {'schema': 'oph.scalar_observations.shared_calibration.v1',
            'classification': 'retrospective conditional shared-preparation comparison',
            'dataset': {'doi': observations['provenance']['doi'], 'version': 'v1',
                        'paper': observations['provenance']['paper'],
                        'simulated_duration_ms': 0, 'calibration_complex_points': 21,
                        'validation_complex_points': 100, 'mirror_copies_counted_as_trials': False},
            'fixed_dictionary': {'family': 'one real-displaced thermal oscillator',
                'mean': 'exp(-(nbar+1/2)*(u*u+v*v)+2i*alpha*v)',
                'grid': 'u=extent*i/10, v=extent*j/10',
                'phase_alias': 'theta=2*alpha*extent/10 modulo 2*pi on the integer grid',
                'same_alpha_and_nbar_for_both_quadratures_and_every_point': True,
                'fitted_gain_offset_or_per_point_parameters': False,
                'coordinate_extent_discrepancy_resolved': False},
            'sampling': {role: row[1] for role, row in sampling.items()},
            'subdivision_depth': depth, 'runs': runs,
            'scope': {'raw_shots_or_accepted_counts_available': False,
                      '250_independent_accepted_shots_per_component_assumed': True,
                      'additional_physical_bias_assumed_zero': True,
                      'reported_mean_rounding_bound_assumed': str(REPORTED_MEAN_ROUNDING),
                      'simultaneous_failure_probability_upper_exact': str(simultaneous_failure),
                      'simultaneous_failure_probability_admission_ceiling': '1/20',
                      'authenticated_experimental_95_percent_agreement_established': False,
                      'independence_between_different_means_required': False,
                      'parameter_cover_is_outer_and_not_fitted_parameter_standard_error': True,
                      'phase_alias_representatives_only_and_no_off_grid_extrapolation': True,
                      'witness_checks_one_unrestricted_semantic_state': True,
                      'confidence_cover_cells_source_executed': False,
                      'floating_source_error_or_physical_service_error_certified': False,
                      'physical_clock_or_interacting_action_tested': False,
                      'preparation_or_action_derived_from_OPH': False,
                      'prospective_prediction': False, 'empirical_confirmation': False},
            'input_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                             for path in SOURCE_PATHS})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', type=Path)
    action.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = build()
    if args.check:
        if strict_json.loads(args.check.read_text(encoding='utf-8')) != result:
            raise SystemExit('shared-preparation comparison drift')
    else:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(_canonical(result))
    print(_canonical({str(run['prediction']['extent']): run['verdict'] for run in result['runs']}), end='')


if __name__ == '__main__':
    main()
