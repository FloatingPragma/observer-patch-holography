"""Certified occupation cutoff for the complete declared scalar instrument.

The reference is the infinite *two-oscillator* Gaussian action, at model
time at most one and native current -1, 0 or 1. A number-moment bound controls
the projected Duhamel residual along the impulse, evolution and both pointer
arms. This includes terminal field states, without conditional normalization.
All certificate calculations are rational, including directed root/exp bounds.
Finite splitting, numerical execution and physical service errors are separate.
"""
from fractions import Fraction as F
from functools import lru_cache
import math

from .bounds import root_upper


def _root_lower(value, scale=10**30):
    return F(math.isqrt(value.numerator*scale*scale//value.denominator), scale)


@lru_cache(maxsize=128)
def exponential_upper(value, digits=30):
    """Rational upper bound on exp(value) by a Taylor/geometric remainder.

    After term m, subsequent term ratios are at most value/(m+2).
    The geometric bound on the tail is positive and requires that ratio < 1.
    Rounding the final bound upward prevents unwieldy returned denominators.
    """
    if not isinstance(value, F) or value < 0 or type(digits) is not int or digits < 1:
        raise ValueError("nonnegative exact rational and positive integer digits required")
    term, total, order = F(1), F(1), 0
    resolution = F(1, 10**digits)
    while True:
        following = term*value/(order+1)
        ratio = value/(order+2)
        if ratio < 1:
            remainder = following/(1-ratio)
            if remainder <= resolution:
                bound = total+remainder
                scale = 10**digits
                return F(-(-bound.numerator*scale//bound.denominator), scale)
        order += 1
        term = following
        total += term


def number_moment_bound(q, displacement_squared):
    """Upper bound on <q**N> for a displaced coupled vacuum.

    The common mode has frequency 1, the difference mode frequency 5/4 and
    squeezing parameter 1/9. The squared coherent displacement is bounded
    in the *original* unit-frequency number convention. Real and imaginary
    normal-mode components are covered by their largest exact coefficient.
    """
    if (not isinstance(q, F) or not F(1) < q < 9
            or not isinstance(displacement_squared, F) or displacement_squared < 0):
        raise ValueError("exact rational 1 < q < 9 and nonnegative displacement required")
    prefactor = root_upper(F(80)/(81-q*q), scale=10**30)
    exponent = (q-1)*8*displacement_squared/(9-q)
    return prefactor*exponential_upper(exponent)


def _sine_interval(value, digits=30):
    """Directed rational enclosure from a decreasing alternating series."""
    if not isinstance(value, F) or not 0 <= value <= F(5, 4):
        raise ValueError("exact rational sine argument between zero and 5/4 required")
    term, total, order = value, value, 0
    scale = 10**digits
    while True:
        following = -term*value*value/((2*order+2)*(2*order+3))
        if abs(following) <= F(1, scale):
            lower, upper = sorted((total, total+following))
            return (F(lower.numerator*scale//lower.denominator, scale),
                    F(-(-upper.numerator*scale//upper.denominator), scale))
        order += 1
        term = following
        total += term


def reference_probability_interval(current=1):
    """Exact rational probability enclosure at model time one.

    This is the specified infinite two-oscillator Gaussian instrument,
    with no measured probability, fitted parameter or floating trigonometry.
    Negative current uses the exact parity relation; zero current is fair.
    """
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("native current -1, 0 or 1 required")
    sine_one = _sine_interval(F(1))
    sine_five_fourths = _sine_interval(F(5, 4))
    mean_lower = (sine_one[0]-F(4, 5)*sine_five_fourths[1])/2
    mean_upper = (sine_one[1]-F(4, 5)*sine_five_fourths[0])/2
    if not 0 <= mean_lower <= mean_upper <= 1:
        raise ArithmeticError("the mean enclosure must certify sine monotonicity")
    signal_lower = _sine_interval(mean_lower)[0]
    signal_upper = _sine_interval(mean_upper)[1]
    exponential = exponential_upper(F(9, 40))
    # The Taylor tail contributes at most 1e-30 and final upward rounding
    # contributes at most 1e-30. Subtracting their sum bounds exp from below.
    exponential_lower = exponential-F(2, 10**30)
    attenuation_lower, attenuation_upper = 1/exponential, 1/exponential_lower
    lower = (1+attenuation_lower*signal_lower)/2
    upper = (1+attenuation_upper*signal_upper)/2
    if current == -1:
        lower, upper = 1-upper, 1-lower
        mean_lower, mean_upper = -mean_upper, -mean_lower
    elif current == 0:
        lower = upper = F(1, 2)
        mean_lower = mean_upper = F(0)
    return {
        "model_time": "1", "native_current": current,
        "mean_lower": str(mean_lower), "mean_upper": str(mean_upper),
        "attenuation_lower": str(attenuation_lower),
        "attenuation_upper": str(attenuation_upper),
        "plus_probability_lower": str(lower), "plus_probability_upper": str(upper),
        "scope": "rational Gaussian probability enclosure for the supplied infinite two-oscillator action and pointer; not empirical data",
    }


def _number_tail(d, displacement_squared):
    # q >= 5/2 ensures n**k/q**n decreases for n>=d>=2 and k=0,1,2.
    # Choosing a minimum of independently valid rational bounds is valid.
    candidates = []
    for numerator in range(5, 17):
        q = F(numerator, 2)
        moment = number_moment_bound(q, displacement_squared)
        candidates.append((moment/q**d, q, moment))
    tail, q, moment = min(candidates)
    return {"q": str(q), "number_exponential_moment_upper": str(moment),
            "occupation_tail_upper": str(tail)}, tail


def _preparation_tail(d):
    # The actual preparation retains total occupations 0,2,...,d-2.
    total, term = F(0), F(1)
    for pair in range(d//2):
        total += term
        term *= F(2*pair+1, 81*(2*pair+2))
    # A certified lower sqrt(5) gives an upper omitted probability.
    return 1-F(4, 9)*_root_lower(F(5), scale=10**(2*d+30))*total


def occupation_cutoff_bound(d, time=F(1), current=1):
    """Ideal finite continuous instrument versus its infinite reference.

    The comparison starts with the existing normalized total-occupation
    vacuum preparation, applies exp(i*j*Q0), evolves for model time t, then
    applies both pointer arms exp(+-i*Q1/2) and reads Pauli Y. It bounds the
    joint classical outcome/terminal quantum field state, after embedding
    the finite registers in the infinite number basis. It never divides by
    a Born weight. This is not a spatial continuum or laboratory certificate.
    """
    if type(d) is not int or d < 2 or d % 2:
        raise ValueError("an even onsite occupation dimension at least two is required")
    if not isinstance(time, F) or not 0 <= time <= 1:
        raise ValueError("exact rational model time between zero and one required")
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("native current -1, 0 or 1 required")

    # Conservation of the displacement's classical oscillator energy gives
    # |alpha|^2 <= j^2/2 during the kick and subsequent H evolution.
    radius = F(current*current, 2)
    evolution_certificate, evolution_tail = _number_tail(d, radius)
    # At time t, |p1|=|j*(cos(t)-cos(5*t/4))/2| <= |j|*t/8.
    # A pointer arm adds momentum +-s/2, 0<=s<=1, giving this uniform radius.
    pointer_radius = radius+abs(current)*time/16+F(1, 8)
    pointer_certificate, pointer_tail = _number_tail(d, pointer_radius)

    # Increase directed-rounding precision with d, so it creates no positive
    # floor in the occupation-refinement limit.
    root_scale = 10**(d+15)
    preparation_tail = _preparation_tail(d)
    preparation = root_upper(2*preparation_tail, scale=root_scale)
    impulse = abs(current)*root_upper(F(d, 2)*evolution_tail, scale=root_scale)
    evolution = time*F(9, 32)*d*root_upper(evolution_tail, scale=root_scale)
    pointer = root_upper(F(d, 8)*pointer_tail, scale=root_scale)
    projected = preparation+impulse+evolution+pointer
    # The projected error and final outside component are orthogonal.
    stinespring = root_upper(projected*projected+pointer_tail, scale=root_scale)
    trace_distance = min(F(1), stinespring)
    return {
        "dimension": d,
        "model_time": str(time),
        "native_current": current,
        "impulse_and_evolution_displacement_squared_upper": str(radius),
        "pointer_displacement_squared_upper": str(pointer_radius),
        "impulse_and_evolution_number_tail": evolution_certificate,
        "pointer_number_tail": pointer_certificate,
        "preparation_omitted_probability_upper": str(preparation_tail),
        "preparation_vector_error_upper": str(preparation),
        "impulse_projected_residual_integral_upper": str(impulse),
        "evolution_projected_residual_integral_upper": str(evolution),
        "pointer_projected_residual_integral_upper": str(pointer),
        "projected_path_vector_error_upper": str(projected),
        "final_outside_probability_upper": str(min(F(1), pointer_tail)),
        "full_stinespring_vector_error_upper": str(stinespring),
        "instrument_trace_distance_upper": str(trace_distance),
        "detector_probability_error_upper": str(trace_distance),
        "scope": "complete ideal continuous finite instrument versus the infinite two-oscillator instrument; includes terminal field branches; excludes finite splitting, numerical execution, physical service and spatial continuum errors",
    }
