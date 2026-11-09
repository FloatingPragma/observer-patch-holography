"""Finite mass-coefficient stability and an actual-source Koide control.

Physical self-energies, a selected gauge vacuum and measured masses are inputs
to a matching theorem; this module supplies none of those inputs.
"""
from fractions import Fraction as F
from math import isqrt

import sympy as s

from model.standard_model_source_map import koide_source


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, F)):
        raise ValueError('an exact integer, rational string or Fraction is required')
    q = F(value)
    return F(q.numerator, q.denominator)


def sqrt_enclosure(value, bits=80):
    """Exact dyadic enclosure; no floating square root enters an admission."""
    q = _rational(value)
    if q < 0 or isinstance(bits, bool) or not isinstance(bits, int) or bits < 1:
        raise ValueError('nonnegative radicand and positive integer precision required')
    scale = 1 << bits
    k = isqrt((q.numerator * scale * scale) // q.denominator)
    low = F(k, scale)
    high = low if low * low == q else F(k + 1, scale)
    return low, high


def coefficient_bound(kinetic_left, kinetic_right, relative_mass_correction,
                      *, matching_error=F(0), bits=80):
    """Conditional relative bound for each ordered positive singular value.

    eps_L/R bound ||Z_L/R-I||; eta bounds ||E||/sigma_min(M0).
    matching_error is a separately supplied relative observable-matching bound,
    and its use does not establish that the bound holds for a physical source.
    """
    left, right, eta, match = map(_rational, (
        kinetic_left, kinetic_right, relative_mass_correction, matching_error))
    if any(not 0 <= x < 1 for x in (left, right, eta, match)):
        raise ValueError('all four relative bounds must belong to [0,1)')
    upper_den = sqrt_enclosure((1 + left) * (1 + right), bits)[1]
    lower_den = sqrt_enclosure((1 - left) * (1 - right), bits)[0]
    if lower_den <= 0:
        raise ValueError('increase precision to resolve the positive kinetic lower bound')
    a = (1 - eta) / upper_den
    b = (1 + eta) / lower_den
    a_match, b_match = (1 - match) * a, (1 + match) * b
    return {'kinetic_left_bound': left, 'kinetic_right_bound': right,
            'relative_mass_correction_bound': eta,
            'effective_lower_multiplier': a, 'effective_upper_multiplier': b,
            'separate_matching_error_hypothesis': match,
            'matched_lower_multiplier': a_match, 'matched_upper_multiplier': b_match,
            'ratio_distortion_upper': b_match / a_match,
            'actual_physical_matching_bound_established': False}


def koide_interval(initial_lower, initial_upper, lower_multiplier, upper_multiplier):
    """Scale-free interval for three positive masses under componentwise bounds."""
    ql, qu, a, b = map(_rational, (
        initial_lower, initial_upper, lower_multiplier, upper_multiplier))
    if not F(1, 3) <= ql <= qu <= 1 or not 0 < a <= b:
        raise ValueError('three-mass Koide interval and positive multipliers required')
    return {'lower': max(F(1, 3), a * ql / b),
            'upper': min(F(1), b * qu / a),
            'initial_lower': ql, 'initial_upper': qu,
            'common_scale_cancels': True}


def absolute_singular_intervals(intervals, kinetic_left, kinetic_right,
                                correction_norm, *, bits=80):
    """Ordered coefficient intervals when only an absolute correction is given."""
    left, right, error = map(_rational, (kinetic_left, kinetic_right, correction_norm))
    if not 0 <= left < 1 or not 0 <= right < 1 or error < 0:
        raise ValueError('positive kinetic lower bounds and nonnegative error required')
    den_u = sqrt_enclosure((1 + left) * (1 + right), bits)[1]
    den_l = sqrt_enclosure((1 - left) * (1 - right), bits)[0]
    if den_l <= 0:
        raise ValueError('increase precision for the kinetic lower bound')
    rows = []
    for interval in intervals:
        if len(interval) != 2:
            raise ValueError('each ordered singular value needs two endpoints')
        lo, hi = map(_rational, interval)
        if not 0 <= lo <= hi:
            raise ValueError('invalid singular-value interval')
        rows.append((max(F(0), lo - error) / den_u, (hi + error) / den_l))
    if not rows:
        raise ValueError('at least one singular-value interval is required')
    return rows


def koide_log_gradient(masses):
    """Exact first variation: d log Q = sum gradient_i d log m_i."""
    m = [s.sympify(x) for x in masses]
    if len(m) != 3 or any(x.is_positive is not True for x in m):
        raise ValueError('three certified positive masses are required')
    roots = [s.sqrt(x) for x in m]
    return tuple(s.simplify(x / sum(m) - r / sum(roots)) for x, r in zip(m, roots))


def source_nonpreservation(root_increment=F(1, 10000)):
    """Exact countercontrol from the native coherent coefficient at phase zero.

    The perturbation is a permitted coefficient change, not a computed loop.
    """
    e = _rational(root_increment)
    if not 0 < e < 1:
        raise ValueError('a positive root increment below one is required')
    ep = s.Rational(e.numerator, e.denominator)
    source = koide_source.exact_source(1)
    y = source['action_coefficient']
    p = s.ImmutableMatrix((s.eye(3) + koide_source.P + koide_source.P**2) / 3)
    if p * p != p or p.H != p or p.rank() != 1:
        raise ArithmeticError('the actual cyclic heavy projector changed')
    z = 4 + 3 * s.sqrt(2)
    light = s.simplify((s.sqrt(s.Rational(10, 72)) - s.sqrt(s.Rational(5, 72)))**2 / 3)
    heavy = s.simplify(light * z**2)
    if (y - light * (s.eye(3) - p) - heavy * p).applyfunc(s.simplify) != s.zeros(3):
        raise ArithmeticError('the source mass eigenspaces changed')
    corrected = (y + ((1 + ep)**2 - 1) * y * p).applyfunc(s.simplify)
    x = (1 + ep) * z
    ratio = s.simplify((2 + x*x) / (2 + x)**2)
    difference = s.simplify(ep * ((2 + ep) * z*z - 8*z) / (3 * (2 + x)**2))
    if s.simplify(ratio - s.Rational(2, 3) - difference) != 0 or difference.is_positive is not True:
        raise ArithmeticError('strict source Koide violation was not established')
    eta = s.simplify(((1 + ep)**2 - 1) * z**2)
    gradient = koide_log_gradient([light, light, heavy])
    return {'coefficient': s.ImmutableMatrix(y), 'heavy_projector': p,
            'corrected_coefficient': s.ImmutableMatrix(corrected),
            'source_singular_values': (light, light, heavy),
            'corrected_singular_values': (light, light, s.simplify((1 + ep)**2 * heavy)),
            'source_Koide': s.Rational(2, 3), 'corrected_Koide': ratio,
            'strict_Koide_increase': difference,
            'relative_mass_correction_norm': eta, 'root_increment': e,
            'log_mass_gradient': gradient,
            'perturbation_is_an_actual_computed_self_energy': False,
            'empirical_mass_inputs': []}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, s.MatrixBase):
        return [[str(value[i, j]) for j in range(value.cols)] for i in range(value.rows)]
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def certificate():
    """JSON-safe finite theorem controls, with no inferred physical inputs."""
    conditional = coefficient_bound(F(1, 1000), F(1, 2000), F(1, 10000))
    a, b = conditional['effective_lower_multiplier'], conditional['effective_upper_multiplier']
    return _json({'theorem': 'two-sided kinetic whitening and coefficient perturbation',
                  'conditional_control': conditional,
                  'conditional_Koide_interval': koide_interval(F(2, 3), F(2, 3), a, b),
                  'actual_source_nonpreservation_control': source_nonpreservation(),
                  'scope': {'actual_source_loop_correction_computed': False,
                            'selected_interacting_gauge_vacuum_supplied': False,
                            'physical_threshold_or_resonance_matching_proved': False,
                            'physical_particle_assignments': [],
                            'empirical_mass_inputs': []}})
