"""Certified rational majorants for the centered composite's energy escape.

The theorem's sharper analytic constants are retained by ``certificate``.
Evaluators intentionally use rational upper bounds for exponentials and roots;
the floating-point choice of delta is labelled as a diagnostic only.
"""
from fractions import Fraction as F
from math import ceil, isfinite, isqrt, log, log1p, sqrt


def rational(value, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('a nonnegative exact rational is required')
    return F(value)


def sqrt_upper(value, bits=80):
    """Relative-precision rational upper square root, including tiny inputs."""
    x = rational(value)
    if type(bits) is not int or not 1 <= bits <= 256:
        raise ValueError('bits must be an integer between 1 and 256')
    if not x:
        return F(0)
    shift = bits + max(0, (x.denominator.bit_length()-x.numerator.bit_length()+1)//2)
    scale = 1 << shift
    n = isqrt(x.numerator*scale*scale//x.denominator)
    if F(n*n, scale*scale) < x:
        n += 1
    return F(n, scale)


def exp_upper(value):
    """exp(x) <= 2**ceil(3x/2), since log(2) >= 2/3 by Jensen."""
    x = rational(value)
    exponent = ceil(3*x/2)
    if exponent > 1_000_000:
        raise ValueError('exponential exceeds the exact calculator budget')
    return F(1 << exponent)


def occupation_defect_upper(coupling_ratio):
    u = rational(coupling_ratio)
    return min(F(2), 11*u*u/8)


def contraction_constants():
    return {'cubic_cubic': 4*9*8, 'cubic_linear': 2*3*2,
            'linear_cubic': 2*3*2, 'total': 312,
            'centered_degree_three_derivation': F(3, 2),
            'free_absolute_kernel_growth': 3}


def microtime_upper(coupling_ratio, time, norm=1):
    """Majorize [3u|s|+312sqrt(Delta)exp(9|s|)]||f||^2."""
    u, t, n = map(rational, (coupling_ratio, time, norm))
    if not u or not t or not n:
        return F(0)
    # sqrt(Delta)<=sqrt(11/8)u<3u/2, with no absolute rounding floor.
    return min(2*n*n, (3*u*t+468*u*exp_upper(9*t))*n*n)


def free_low_energy_upper(epsilon, norm=1):
    e, n = map(rational, (epsilon, norm))
    return min(F(1), 32*e**6)*n*n


def free_gaussian_upper(delta, norm=1):
    d = rational(delta, positive=True)
    n = rational(norm)
    return min(F(1), 1536*d**6)*n*n


def gaussian_low_energy_upper(kinetic_scale, coupling, cutoff, delta, norm=1):
    """Exact upper bound on combined d addition/removal mass in [-R,R].

    For nonzero U this evaluates a conservative rational version of (C):
    sqrt(2/pi)<=1 and sqrt(Delta)<=3u/2. U=0 uses the sharper free bound.
    """
    c = rational(kinetic_scale, positive=True)
    U, R, n = map(rational, (coupling, cutoff, norm))
    d = rational(delta, positive=True)
    if not U:
        return free_low_energy_upper(R/c, n)
    if not n:
        return F(0)
    u = U/c
    terms = 1536*d**6+3*u/d+936*u*exp_upper(F(81, 2)/(d*d))
    return min(F(1), exp_upper(R*R/(2*c*c*d*d))*terms)*n*n


def asymptotic_delta_diagnostic(kinetic_scale, coupling, cutoff):
    """Floating diagnostic only; None at U=0, where no smoothing is needed."""
    values = (kinetic_scale, coupling, cutoff)
    if any(type(x) not in (int, float, F) for x in values):
        raise ValueError('finite real diagnostic parameters required')
    try:
        c, U, R = map(float, values)
    except (OverflowError, ValueError) as exc:
        raise ValueError('parameters exceed floating diagnostic range') from exc
    if not all(map(isfinite, (c, U, R))) or c <= 0 or U < 0 or R < 0:
        raise ValueError('c>0, U>=0 and R>=0 must be finite')
    if U == 0:
        return None
    if U >= c:
        raise ValueError('asymptotic choice requires 0<U/c<1')
    # Widely separated values need a difference of logs to avoid U/c
    # underflow; adjacent large floats need log1p to avoid cancellation.
    logarithm = log1p((c-U)/U) if U >= c/2 else log(c)-log(U)
    result = max(R/c, 9/sqrt(logarithm))
    if not isfinite(result):
        raise ValueError('delta exceeds floating diagnostic range')
    return result


def vacuum_duhamel_upper(coupling_bound, time, f_norm, hf_norm, cutoff, mass_upper):
    """Uniform vacuum-vector error, for both field and adjoint.

    mass_upper bounds the combined d spectral mass for EVERY normalized input.
    It is not a bound inferred from just the particular packet f.
    """
    U, T, n, hn, M = map(rational, (coupling_bound, time, f_norm, hf_norm, mass_upper))
    R = rational(cutoff, positive=True)
    if M > 1:
        raise ValueError('normalized spectral mass upper bound exceeds one')
    if not n and hn:
        raise ValueError('a zero packet cannot have nonzero kinetic norm')
    if not U or not T or not n:
        return F(0)
    refined = U/2*(T*sqrt_upper(M)*n+(2*n+T*hn)/R)
    return min(U*T*n/2, refined)


def certificate():
    examples = []
    # These very small lattice spacings exhibit convergence of deliberately
    # conservative constants. Report only exact rational comparison targets,
    # not enormous expanded denominators or floating underflows.
    for denominator, target in ((4, F(751, 1000)), (8, F(3, 250)), (16, F(1, 5000))):
        d = F(1, denominator)
        exponent = ceil(F(243, 4)/(d*d))
        power = 2*exponent
        value = gaussian_low_energy_upper(1 << power, 1, 1, d)
        assert value < target
        examples.append({'delta': str(d), 'kinetic_scale_power_of_two': power,
                         'coupling': '1', 'physical_cutoff': '1',
                         'mass_upper_less_than': str(target), 'exact_comparison': True})
    return {
        'scope': 'centered full ground; c->infinity, U/c->0; all normalized lattice smears',
        'occupation_defect_upper': 'min(2,11(U/c)^2/8)',
        'microtime_characteristic_difference': '3u|s|+312sqrt(Delta)exp(9|s|)',
        'microtime_constants': {k: str(v) for k, v in contraction_constants().items()},
        'free_composite_low_energy_mass': 'min(1,32epsilon^6)',
        'free_composite_has_no_linear_wick_component': True,
        'all_eight_weyl_nodes_retained': True,
        'gaussian_bound': 'exp[R^2/(2c^2delta^2)](1536delta^6+3u sqrt(2/pi)/delta+624sqrt(Delta)exp[81/(2delta^2)])',
        'asymptotic_delta': 'max(R/c,9/sqrt(log(1/u))) for 0<u<1',
        'zero_coupling_endpoint': 'min(1,32(R/c)^6)',
        'uniform_finite_energy_composite_escape': True,
        'vacuum_duhamel_bound': '(Umax/2)[T sqrt(M_a(R))||f||+(2||f||+T||h_a f||)/R]',
        'vacuum_duhamel_requires_uniform_input_mass_bound': True,
        'vacuum_duhamel_applies_to_both_adjoint_branches': True,
        'exact_evaluator_majorants': ['exp(x)<=2^ceil(3x/2)', 'sqrt(Delta)<=3u/2', 'sqrt(2/pi)<=1'],
        'floating_delta_is_diagnostic_only': True,
        'exact_refinement_examples': examples,
    }
