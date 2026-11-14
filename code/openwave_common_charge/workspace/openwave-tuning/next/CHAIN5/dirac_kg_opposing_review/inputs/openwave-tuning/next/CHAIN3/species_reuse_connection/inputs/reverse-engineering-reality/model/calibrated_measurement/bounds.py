"""Counts-to-moment certificates using exact scalar arithmetic.

Interface bounds are trusted results of fresh operator certification, not
authenticated caller-authored claims. Physical service bounds remain supplied.
Every count is retained; no fitted quantum state or simulated shot is used.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math


FIRST = ('null_plus', 'null_minus', 'h_real', 'h_imag')
SECOND = ('2h_real', '2h_imag')
SETTINGS = FIRST + SECOND
GRID = 10**18


def _fraction(value, name, *, positive=False):
    if type(value) is not F or value < 0 or (positive and value == 0):
        raise ValueError(f'exact nonnegative Fraction {name} required')
    return value


def root_interval(value):
    """Fixed rational enclosure, including exact squares and zero."""
    _fraction(value, 'radicand')
    numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if numerator*numerator == value.numerator and denominator*denominator == value.denominator:
        root = F(numerator, denominator)
        return root, root
    lower = math.isqrt(value.numerator*GRID*GRID//value.denominator)
    return F(lower, GRID), F(lower+1, GRID)


def log_interval(value, terms=32):
    """Enclose log(value)>=0 by range reduction and the atanh series.

    For 0<=t<=1/3 the omitted positive tail is at most
    2*t^(2m+1)/((2m+1)*(1-t^2)). No floating logarithm enters the result.
    """
    _fraction(value, 'logarithm input', positive=True)
    if value < 1 or type(terms) is not int or not 1 <= terms <= 128:
        raise ValueError('logarithm input >=1 and 1..128 terms required')
    power = 0
    while value >= 2:
        value /= 2
        power += 1
    def reduced(x):
        t = (x-1)/(x+1)
        lower = 2*sum((t**(2*j+1)/(2*j+1) for j in range(terms)), F(0))
        tail = 2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
        return lower, lower+tail
    a, b = reduced(value), reduced(F(2))
    return a[0]+power*b[0], a[1]+power*b[1]


def sampling_radius(trials, failure_probability=F(1, 100)):
    """Simultaneous guarantee for the fixed six-setting family.

    Independent +/-1 trials within each setting; fixed trial budgets and
    stable means. Independence between settings is not needed by the union
    bound. Four-setting prediction prefixes retain the same six allocations.
    """
    if type(trials) is not int or not 1 <= trials <= 10**12:
        raise ValueError('1..10^12 predetermined trials required')
    _fraction(failure_probability, 'failure probability', positive=True)
    if failure_probability >= 1:
        raise ValueError('failure probability must be less than one')
    log_upper = log_interval(2*len(SETTINGS)/failure_probability)[1]
    return root_interval(2*log_upper/trials)[1]


@dataclass(frozen=True)
class InterfaceBounds:
    mass: tuple[F, ...]
    null: tuple[F, F]
    moment: F

    def __post_init__(self):
        if type(self.mass) is not tuple or len(self.mass) != 6 or type(self.null) is not tuple or len(self.null) != 2:
            raise ValueError('six mass and two null bounds required')
        for value in (*self.mass, *self.null, self.moment):
            _fraction(value, 'interface error')
        if any(value >= 1 for value in self.mass):
            raise ValueError('effect-sum defect must be below one')

    def completion(self, setting):
        c = self.mass[SETTINGS.index(setting)]
        # T=sum K* K, L=E0-E1. The complete comparison instrument K T^-1/2
        # has observable T^-1/2 L T^-1/2, within c(1+c)/(1-c) of L.
        return c*(1+c)/(1-c)


def _counts(counts):
    if type(counts) is not dict or set(counts) not in (set(FIRST), set(SETTINGS)):
        raise ValueError('the four prediction settings or all six settings are required')
    means, sizes = {}, {}
    for setting, row in counts.items():
        if type(row) not in (tuple, list) or len(row) != 4 or any(type(n) is not int or n < 0 for n in row):
            raise ValueError('four nonnegative integer spin/bit counts per setting required')
        size = sum(row)
        if not 1 <= size <= 10**12:
            raise ValueError('positive bounded count total required')
        sizes[setting] = size
        means[setting] = F(row[0]-row[1]+row[2]-row[3], size)
    return means, sizes


def _product_interval(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def analyze_counts(counts, interface, *, service_mean_errors, failure_probability=F(1, 100)):
    """A conditional decision from complete counts and justified mean errors.

    A None service entry leaves the physical conclusion unresolved. Explicit
    zeros describe the complete mathematical comparison instrument, unless
    independent evidence justifies those zeros for an apparatus.
    """
    if type(interface) is not InterfaceBounds:
        raise ValueError('fresh trusted InterfaceBounds required')
    means, sizes = _counts(counts)
    if type(service_mean_errors) is not dict or set(service_mean_errors) != set(SETTINGS):
        raise ValueError('all six physical service mean-error obligations required')
    for value in service_mean_errors.values():
        if value is not None:
            _fraction(value, 'service mean error')
    radii = {name: sampling_radius(size, failure_probability) for name, size in sizes.items()}
    result = {'raw_means': means, 'trials': sizes, 'sampling_radii': radii,
              'family_failure_probability_upper': failure_probability,
              'status': 'unbounded_services', 'rejects_joint_premises': False}
    if any(value is None for value in service_mean_errors.values()):
        return result
    service = service_mean_errors
    gain, offset = (means['null_plus']-means['null_minus'])/2, (means['null_plus']+means['null_minus'])/2
    calibration = sum((radii[name]+service[name]+interface.null[i]+interface.completion(name)
                       for i, name in enumerate(FIRST[:2])), F(0))/2
    result.update({'estimated_gain': gain, 'estimated_offset': offset,
                   'gain_and_offset_error_upper': calibration,
                   'true_gain_interval': (gain-calibration, min(F(1), gain+calibration))})
    if gain <= calibration:
        result['status'] = 'inconclusive_gain'
        return result
    corrected = {name: (means[name]-offset)/gain for name in means if name not in FIRST[:2]}
    errors = {name: (radii[name]+service[name]+interface.completion(name)
                    +(2+interface.mass[SETTINGS.index(name)])*calibration)/gain
              for name in corrected}
    z = corrected['h_real'], corrected['h_imag']
    z_squared = z[0]*z[0]-z[1]*z[1], 2*z[0]*z[1]
    z_norm_squared = z[0]*z[0]+z[1]*z[1]
    z_norm_upper = root_interval(z_norm_squared)[1]
    e1 = root_interval(errors['h_real']**2+errors['h_imag']**2)[1]
    # Enclose the possible actual B moment using only the first complex read.
    radius = 1-z_norm_squared+interface.moment+4*z_norm_upper*e1+e1*e1
    first_impossible = z_norm_squared > (root_interval(1+interface.moment)[1]+e1)**2
    raw_predictions = {}
    for index, name in enumerate(SECOND):
        completion = interface.completion(name)
        low = max(F(-1), z_squared[index]-max(F(0), radius)-completion)
        high = min(F(1), z_squared[index]+max(F(0), radius)+completion)
        if low > high:
            first_impossible = True
            low, high = F(-1), F(1)
        product = _product_interval((gain-calibration, min(F(1), gain+calibration)), (low, high))
        raw_predictions[name] = (max(F(-1), offset-calibration+product[0]-service[name]),
                                 min(F(1), offset+calibration+product[1]+service[name]))
    result.update({'corrected_means': corrected, 'quadrature_error_upper': errors,
                   'first_complex_error_upper': e1, 'first_moment_incompatible': first_impossible,
                   'second_observed_mean_predictions': raw_predictions,
                   'prediction_uses_second_setting_counts': False,
                   'status': 'prediction_only', 'rejects_joint_premises': first_impossible})
    if set(counts) == set(SETTINGS):
        w = corrected['2h_real'], corrected['2h_imag']
        e2 = root_interval(errors['2h_real']**2+errors['2h_imag']**2)[1]
        distance_lower = root_interval(sum(((w[i]-z_squared[i])**2 for i in range(2)), F(0)))[0]
        excess = distance_lower+z_norm_squared-1
        allowance = interface.moment+e2+4*z_norm_upper*e1+2*e1*e1
        reject = first_impossible or excess > allowance
        result.update({'second_complex_error_upper': e2, 'moment_excess_lower': excess,
                       'moment_allowance_upper': allowance, 'rejection_margin_lower': excess-allowance,
                       'status': 'reject' if reject else 'not_rejected', 'rejects_joint_premises': reject})
    return result


def json_safe(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: json_safe(item) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [json_safe(item) for item in value]
    return value
