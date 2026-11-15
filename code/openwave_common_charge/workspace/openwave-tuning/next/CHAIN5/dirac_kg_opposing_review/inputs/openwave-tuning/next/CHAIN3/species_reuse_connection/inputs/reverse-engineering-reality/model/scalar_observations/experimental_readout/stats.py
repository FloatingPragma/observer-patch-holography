"""Frozen-axis predictions and assumption-conditioned characteristic diagnostics.

Each reported real or imaginary characteristic mean is treated as the mean
of ``samples`` independent observations in [-1, 1]. Processed or postselected
source means do not establish that premise. The resulting coverage statements
are conditional diagnostics, not authenticated experimental confidence claims.
Integer coordinate keys identify grid points without inferring physical units.
"""
from __future__ import annotations

import cmath
from fractions import Fraction as F
import math
from numbers import Real


def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number')
    return float(value)


def _key(value):
    if (not isinstance(value, tuple) or len(value) != 2
            or any(type(component) is not int for component in value)):
        raise ValueError('coordinate keys must be pairs of integer grid indices')
    return value


def _means(values):
    if not isinstance(values, dict) or not values:
        raise ValueError('a nonempty dictionary of characteristic means is required')
    result = {}
    for key, value in values.items():
        _key(key)
        if isinstance(value, bool) or not isinstance(value, (Real, complex)):
            raise ValueError('each characteristic mean must be a finite complex number')
        value = complex(value)
        if any(not math.isfinite(part) or not -1 <= part <= 1
               for part in (value.real, value.imag)):
            raise ValueError('each real and imaginary mean must lie in [-1, 1]')
        result[key] = value
    return result


def _sample_contract(statistics, samples):
    if type(statistics) is not int or statistics < 1:
        raise ValueError('a positive integer number of real statistics is required')
    if type(samples) is not int or samples < 1:
        raise ValueError('a positive integer sample count is required')


def concentration_certificate(statistics, samples, width):
    """Rational upper bound on ``2*m*exp(-n*width**2/2)``.

    The exact rational value of the supplied floating width is used. A
    positive Taylor partial sum bounds the exponential from below. Flooring
    that exact sum to 64 fractional bits makes the retained rational numbers
    shorter and only enlarges the failure bound. No floating logarithm,
    exponential or square root participates in this certificate.
    """
    _sample_contract(statistics, samples)
    width = _real(width, 'interval half-width')
    if width < 0:
        raise ValueError('interval half-width cannot be negative')
    exact_width = F(width)
    argument = samples*exact_width*exact_width/2
    term = total = F(1)
    degree = 128
    for order in range(1, degree+1):
        term *= argument/order
        total += term
    denominator = 1 << 64
    lower = F((total.numerator*denominator)//total.denominator, denominator)
    failure = min(F(1), 2*statistics/lower)
    return {'simultaneous_real_statistics': statistics,
            'samples_per_real_mean_assumed': samples,
            'width_exact': str(exact_width), 'exponential_argument': str(argument),
            'taylor_degree': degree, 'taylor_lower_rounding_fractional_bits': 64,
            'exp_lower_partial_sum_rounded_down': str(lower),
            'familywise_error_upper': str(failure),
            'coverage': 'conditional_on_independent_bounded_raw_observations'}


def hoeffding_width(statistics, samples=250, delta=.025):
    """Certified simultaneous half-width for bounded [-1, 1] means.

    Independence is required within each mean's repeated observations. The
    union bound requires no independence between different scored means.
    ``statistics`` counts real components, not complex grid points.
    """
    _sample_contract(statistics, samples)
    delta = _real(delta, 'failure probability')
    if not 0 < delta < 1:
        raise ValueError('failure probability must lie strictly between zero and one')
    # A sum of logs avoids overflowing 2*statistics/delta unnecessarily.
    width = math.sqrt(2*(math.log(2)+math.log(statistics)-math.log(delta))/samples)
    required = F(str(delta))
    attempts = 0
    while F(concentration_certificate(statistics, samples, width)['familywise_error_upper']) > required:
        # Usually at most a few adjacent floats are needed. Extreme failure
        # budgets can require a materially larger exponent when using the
        # fixed Taylor degree, so avoid an unbounded adjacent-float search.
        attempts += 1
        width = math.nextafter(width if attempts <= 8 else 2*width, math.inf)
    return width


def _interval(value):
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError('an interval requires two ordered finite endpoints')
    lower, upper = (_real(endpoint, 'interval endpoint') for endpoint in value)
    if lower > upper:
        raise ValueError('interval endpoints are reversed')
    return lower, upper


def _box(value):
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        raise ValueError('a complex box requires real and imag intervals')
    return {part: _interval(value[part]) for part in ('real', 'imag')}


def mean_box(mean, width):
    """Componentwise interval, intersected with the known [-1, 1] range."""
    mean = _means({(0, 0): mean})[(0, 0)]
    width = _real(width, 'interval half-width')
    if width < 0:
        raise ValueError('interval half-width cannot be negative')
    if width == 0:
        return {'real': (mean.real, mean.real), 'imag': (mean.imag, mean.imag)}
    return {name: (max(-1., math.nextafter(part-width, -math.inf)),
                   min(1., math.nextafter(part+width, math.inf)))
            for name, part in (('real', mean.real), ('imag', mean.imag))}


def _product(left, right):
    products = [a*b for a in left for b in right]
    if any(not math.isfinite(value) for value in products):
        raise ValueError('interval product overflowed finite arithmetic')
    return math.nextafter(min(products), -math.inf), math.nextafter(max(products), math.inf)


def box_product(left, right):
    """Enclose every complex product, including all signs of its components.

    Coordinate extrema of a bilinear product occur at rectangle vertices.
    Each ordinary arithmetic step is rounded outward with ``nextafter``.
    """
    left, right = _box(left), _box(right)
    ac, bd = _product(left['real'], right['real']), _product(left['imag'], right['imag'])
    ad, bc = _product(left['real'], right['imag']), _product(left['imag'], right['real'])
    return {'real': (math.nextafter(ac[0]-bd[1], -math.inf),
                     math.nextafter(ac[1]-bd[0], math.inf)),
            'imag': (math.nextafter(ad[0]+bc[0], -math.inf),
                     math.nextafter(ad[1]+bc[1], math.inf))}


def conjugate_box(box):
    box = _box(box)
    return {'real': box['real'], 'imag': (-box['imag'][1], -box['imag'][0])}


def calibrate_axes(axis_means, samples=250, delta=.025):
    """Freeze only measured nonnegative axes; reject any withheld off-axis data.

    The origin can be included as an explicit calibration observation. It
    contributes to the conservative simultaneous count. No measured value is
    divided by the origin and no off-axis response selects a parameter.
    """
    means = _means(axis_means)
    if any(i < 0 or j < 0 or (i != 0 and j != 0) for i, j in means):
        raise ValueError('calibration accepts only nonnegative coordinate axes')
    width = hoeffding_width(2*len(means), samples, delta)
    return {'axis_means': means,
            'axis_boxes': {key: mean_box(value, width) for key, value in means.items()},
            'real_statistics': 2*len(means), 'samples_per_real_mean_assumed': samples,
            'failure_probability_assumed': float(delta), 'simultaneous_half_width': width,
            'concentration_certificate': concentration_certificate(2*len(means), samples, width),
            'coverage': 'conditional_on_independent_bounded_raw_observations',
            'processed_means_establish_sampling_premise': False,
            'empirical_confirmation': False,
            'model_assumption': 'aligned Gaussian axis factorization',
            'calibration_uses_withheld_values': False}


def predict(calibration, keys):
    """Predict off-axis means by the product of two frozen axis means.

    Negative axes use characteristic conjugation, not new observations.
    At zero a characteristic function equals one; the recorded origin is
    preserved for diagnostics and is not used to renormalize any response.
    """
    if not isinstance(calibration, dict) or not {'axis_means', 'axis_boxes'} <= set(calibration):
        raise ValueError('a frozen axis calibration is required')
    means = _means(calibration['axis_means'])
    if any(i < 0 or j < 0 or (i != 0 and j != 0) for i, j in means):
        raise ValueError('calibration contains non-axis or negative-axis observations')
    if not isinstance(calibration['axis_boxes'], dict) or set(calibration['axis_boxes']) != set(means):
        raise ValueError('every frozen axis mean requires its corresponding box')
    boxes = {key: _box(value) for key, value in calibration['axis_boxes'].items()}

    def axis(index, coordinate):
        if index == 0:
            return 1+0j, {'real': (1., 1.), 'imag': (0., 0.)}
        key = (abs(index), 0) if coordinate == 0 else (0, abs(index))
        if key not in means:
            raise ValueError(f'missing calibration axis {key}')
        mean, box = means[key], boxes[key]
        return (mean, box) if index > 0 else (mean.conjugate(), conjugate_box(box))

    result = {}
    for key in keys:
        i, j = _key(key)
        if key in result:
            raise ValueError('prediction coordinates must not be duplicated')
        horizontal, vertical = axis(i, 0), axis(j, 1)
        result[key] = {'mean': horizontal[0]*vertical[0],
                       'box': box_product(horizontal[1], vertical[1])}
    if not result:
        raise ValueError('at least one prediction coordinate is required')
    return result


def score(predictions, observed_means, samples=250, delta=.025):
    """Compare withheld complex responses with calibration-propagated boxes.

    Rejection is a disjoint component interval, using simultaneous sampling
    boxes for withheld observations. Calibration failure probability is a
    separate budget and must be added for a joint conditional guarantee.
    The numeric diagnostics alone establish neither premise nor provenance.
    """
    observed = _means(observed_means)
    if not isinstance(predictions, dict) or set(predictions) != set(observed):
        raise ValueError('predictions and withheld observations must have identical coordinates')
    width = hoeffding_width(2*len(observed), samples, delta)
    rows = {}
    squared_errors, errors, violations = [], [], []
    for key, value in observed.items():
        prediction = predictions[key]
        if not isinstance(prediction, dict) or set(prediction) != {'mean', 'box'}:
            raise ValueError('each prediction requires its complex mean and box')
        mean = prediction['mean']
        if isinstance(mean, bool) or not isinstance(mean, (Real, complex)):
            raise ValueError('predicted means must be finite complex numbers')
        mean = complex(mean)
        if not math.isfinite(mean.real) or not math.isfinite(mean.imag):
            raise ValueError('predicted means must be finite complex numbers')
        target, measured = _box(prediction['box']), mean_box(value, width)
        parts = [part for part in ('real', 'imag')
                 if measured[part][1] < target[part][0] or measured[part][0] > target[part][1]]
        error = abs(value-mean)
        errors.append(error)
        squared_errors.append(error*error)
        rows[key] = {'prediction': mean, 'observation': value, 'prediction_box': target,
                     'observation_box': measured, 'disjoint_components': parts,
                     'absolute_complex_error': error}
        if parts:
            violations.append(key)
    return {'complex_rmse': math.sqrt(math.fsum(squared_errors)/len(observed)),
            'complex_mae': math.fsum(errors)/len(observed),
            'withheld_points': len(observed), 'real_statistics': 2*len(observed),
            'samples_per_real_mean_assumed': samples, 'failure_probability_assumed': float(delta),
            'simultaneous_half_width': width, 'box_violations': violations,
            'concentration_certificate': concentration_certificate(2*len(observed), samples, width),
            'box_violation_count': len(violations), 'rows': rows,
            'coverage': 'conditional_on_independent_bounded_raw_observations',
            'processed_means_establish_sampling_premise': False,
            'calibration_failure_budget_included': False,
            'physical_service_error_budget_included': False,
            'empirical_confirmation': False, 'provenance_verified': False,
            'verdict': 'conditional_box_violation' if violations else 'not_rejected_under_assumed_sampling'}


def nominal_gaussian(u, v, thermal=.04, amplitude=1.5):
    """Supplied displaced thermal characteristic, without fitting data.

    Default convention: exp[-.54*(u*u+v*v)-i*sqrt(2)*1.5*v].
    This function cannot identify the source dataset's coordinate units.
    """
    u, v = _real(u, 'u'), _real(v, 'v')
    thermal, amplitude = _real(thermal, 'thermal occupation'), _real(amplitude, 'amplitude')
    if thermal < 0:
        raise ValueError('thermal occupation cannot be negative')
    return cmath.exp(-(.5+thermal)*(u*u+v*v)-1j*math.sqrt(2)*amplitude*v)
