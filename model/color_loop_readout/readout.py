"""Complete rational probability enclosures and bare coefficient inversion."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

import sympy as sp

from .source import Controls, Sector, add, fraction, sqrt_upper


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        object.__setattr__(self, 'lo', F(self.lo))
        object.__setattr__(self, 'hi', F(self.hi))
        if self.lo > self.hi:
            raise ValueError('ordered interval required')

    def __add__(self, other):
        if not isinstance(other, Interval):
            other = Interval(other, other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    def scale(self, value):
        x, y = self.lo*value, self.hi*value
        return Interval(min(x, y), max(x, y))

    def widen(self, radius):
        if radius < 0:
            raise ValueError('nonnegative error required')
        return Interval(self.lo-radius, self.hi+radius)

    def divide_positive(self, other):
        if other.lo <= 0:
            raise ValueError('strictly positive denominator enclosure required')
        values = [a/b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    def positive_square(self):
        if self.lo < 0:
            raise ValueError('nonnegative interval required')
        return Interval(self.lo*self.lo, self.hi*self.hi)

    def positive_sqrt(self, bits=48):
        if self.lo < 0:
            raise ValueError('nonnegative interval required')
        upper = sqrt_upper(self.hi, bits)
        lower_upper = sqrt_upper(self.lo, bits)
        lower = lower_upper if lower_upper*lower_upper == self.lo else lower_upper-F(1, 1 << bits)
        return Interval(lower, upper)

    def contains(self, value):
        return self.lo <= value <= self.hi

    def packet(self, bits=128):
        """Compact outward dyadic endpoints; serialization never narrows a bound."""
        scale = 1 << bits
        a, b = self.lo*scale, self.hi*scale
        return [str(F(a.numerator//a.denominator, scale)),
                str(F(-(-b.numerator//b.denominator), scale))]


def complete_probabilities(sector, t, order=28, powers=None):
    """All joint (Higgs occupation, link block pair) records, no conditioning."""
    t = fraction(t)
    if type(order) is not int or order < 1:
        raise ValueError('positive Taylor order required')
    powers = sector.powers(order) if powers is None else powers
    if len(powers) != order+1:
        raise ValueError('complete Taylor vector catalog required')
    real, imag = {}, {}
    for k, v in enumerate(powers):
        factor = t**k/F(factorial(k))
        target = real if k % 2 == 0 else imag
        factor *= (1, -1, -1, 1)[k % 4]
        for i, value in v.items():
            add(target, i, factor*value)
    radius = (abs(t)*sector.norm_bound)**(order+1)/factorial(order+1)
    # For the complete instrument, the same isometry bound gives a trace-norm
    # error at most 2r+r^2 on arbitrary inputs, including untouched references.
    error = 2*radius+radius*radius
    values = {}
    for key, indices in sector.records.items():
        keep = set(indices)
        p = sector.inner(real, real, keep)+sector.inner(imag, imag, keep)
        values[key] = Interval(max(F(0), p-error), min(F(1), p+error))
    return values, radius


def coarse_probability(records, block):
    out = Interval(0, 0)
    for (nh, left, right), value in records.items():
        if left == right == block:
            out = out+value
    return out


@lru_cache(None)
def interpolation_weights(count):
    """Exact first two coefficient rows of the even-power Vandermonde inverse."""
    if type(count) is not int or not 2 <= count <= 8:
        raise ValueError('between two and eight complete time settings required')
    inverse = sp.Matrix([[j**(2*k) for k in range(1, count+1)]
                         for j in range(1, count+1)]).inv()
    return tuple(tuple(F(int(inverse[k, j].p), int(inverse[k, j].q))
                       for j in range(count)) for k in (0, 1))


def fit_coefficients(probabilities, step, norm_bound, *, observation_error=F(0)):
    """Enclose c2 and c4 from positive-time probabilities with rigorous tails."""
    step, observation_error = fraction(step), fraction(observation_error)
    if step <= 0 or norm_bound < 0 or observation_error < 0:
        raise ValueError('positive step and nonnegative bounds required')
    count = len(probabilities)
    weights = interpolation_weights(count)
    result = []
    for k, row in enumerate(weights, start=1):
        interval = Interval(0, 0)
        for j, (weight, probability) in enumerate(zip(row, probabilities, strict=True), start=1):
            # Real H, seed and record projectors make p(t) even. Taylor's
            # integral remainder uses ||ad_H^(2K+2)(P)|| <= (2B)^(2K+2).
            remainder = (2*norm_bound*j*step)**(2*count+2)/factorial(2*count+2)
            interval = interval+probability.widen(remainder+observation_error).scale(weight/step**(2*k))
        result.append(interval)
    return tuple(result)


def invert_color_coefficients(coefficients):
    """c3=4 and 3 distinguish the electric coefficient independently of Higgs."""
    m2 = coefficients[0][0].scale(F(1, 9))
    # The two selected representations have d=3 and c^2 difference 16-9=7.
    contrast = coefficients[0][1].scale(F(1, 9))+coefficients[3][1].scale(F(-1, 9))
    product = contrast.scale(F(-3, 7))  # e^2 m^2
    e2 = product.divide_positive(m2)
    ratio_squared = product.divide_positive(m2.positive_square())
    return {'magnetic_squared':m2, 'electric_squared':e2,
            'electric_over_abs_magnetic_squared':ratio_squared,
            'electric_over_abs_magnetic':ratio_squared.positive_sqrt()}


def sample_budget(outcomes, settings, tolerance, failure_probability):
    """Integer Hoeffding/union bound; log(x) <= ceil(log2(x)) for x >= 1."""
    if type(outcomes) is not int or type(settings) is not int or min(outcomes, settings) < 1:
        raise ValueError('positive finite observation catalog required')
    tolerance, failure_probability = fraction(tolerance), fraction(failure_probability)
    if not 0 < tolerance < 1 or not 0 < failure_probability < 1:
        raise ValueError('probability and tolerance in (0,1) required')
    target = F(2*outcomes*settings, 1)/failure_probability
    k = 0
    while 1 << k < target:
        k += 1
    bound = F(k, 2)/tolerance**2
    return (bound.numerator+bound.denominator-1)//bound.denominator


def witness(q, *, controls=Controls(), step=F(1, 1024), settings=6, order=28):
    if controls.electric <= 0 or controls.magnetic == 0:
        raise ValueError('positive electric coefficient and nonzero loop coupling required')
    sector = Sector(q, controls)
    powers = sector.powers(order)
    probability_tables = []
    selected = {r:[] for r in (0, 2, 3)}
    radii = []
    for j in range(1, settings+1):
        records, radius = complete_probabilities(sector, j*step, order, powers)
        for r in selected:
            selected[r].append(coarse_probability(records, r))
        if not sum(x.lo for x in records.values()) <= 1 <= sum(x.hi for x in records.values()):
            raise ArithmeticError('complete probability enclosure lost unit mass')
        radii.append(radius)
        probability_tables.append({'time':str(j*step), 'amplitude_error':str(radius),
            'records':[{'N_H':nh, 'link_blocks':[left, right], 'probability':p.packet()}
                       for (nh, left, right), p in records.items()]})
    coefficients = {r:fit_coefficients(p, step, sector.norm_bound) for r, p in selected.items()}
    values = invert_color_coefficients(coefficients)
    expected = {'magnetic_squared':controls.magnetic**2,
                'electric_squared':controls.electric**2,
                'electric_over_abs_magnetic_squared':(controls.electric/controls.magnetic)**2,
                'electric_over_abs_magnetic':controls.electric/abs(controls.magnetic)}
    if not all(values[k].contains(v) for k, v in expected.items()):
        raise ArithmeticError('full source coefficient readout misses specified controls')
    observed = {r:fit_coefficients(p, step, sector.norm_bound, observation_error=F(1, 2**80))
                for r, p in selected.items()}
    sampled = invert_color_coefficients(observed)
    return {'Qstar':q, 'dimension':sector.dimension, 'seed_basis_index':sector.seed,
            'controls':controls.packet(), 'Hamiltonian_norm_upper_bound':str(sector.norm_bound),
            'complete_record_count':len(sector.records),
            'complete_record_rank_sum':sum(map(len, sector.records.values())),
            'sparse_generator_entries':len(sector.numerators),
            'clock_step':str(step), 'settings':settings, 'Taylor_order':order,
            'all_record_probabilities':probability_tables,
            'coefficient_intervals':{str(r):{'t2':v[0].packet(), 't4':v[1].packet()} for r, v in coefficients.items()},
            'bare_coefficient_readout':{k:v.packet() for k, v in values.items()},
            'sampled_probability_tolerance':'1/'+str(2**80),
            'sampled_coefficient_readout':{k:v.packet() for k, v in sampled.items()},
            'samples_per_setting_for_two_sectors_failure_at_most_1e_minus6':str(
                sample_budget(3, 2*settings, F(1, 2**80), F(1, 10**6))),
            'complete_reference_instrument_trace_norm_error_upper_bound':str(2*max(radii)+max(radii)**2)}
