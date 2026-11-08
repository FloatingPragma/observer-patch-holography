"""Exact bounds for the same ordered-spin continuum and complete readouts."""
from fractions import Fraction as F


def rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(f'{name} must be an exact integer or Fraction')
    return F(value)


def continuum_response(channel, m0, radius, response=F(1), kappa_upper=F(1, 2)):
    m, r, b, k = [rational(x, n) for x, n in
                   zip((m0, radius, response, kappa_upper), ('m0', 'radius', 'response', 'kappa_upper'))]
    if not 0 < m <= F(1, 2) or r <= 0 or not m*m/2 <= k <= F(1, 2):
        raise ValueError('require positive radius, 0<m0<=1/2 and m0^2/2<=kappa_upper<=1/2')
    if channel == 'P':
        n, h, upper = m*m*r/24, k*r*r/2, F(1, 12)
    elif channel == 'Q':
        n, h, upper = m*m/r, 6*k, 4/(r*r)
    else:
        raise ValueError('channel must be P or Q')
    lower = n*n/h
    if not 0 < lower <= upper:
        raise ArithmeticError('inconsistent source susceptibility bounds')
    return {'channel': channel, 'm0': m, 'momentum_inner': r/2,
            'momentum_outer': r, 'unit_test_L2_norm': True,
            'excitation_norm_squared_lower': n, 'first_energy_upper': h,
            'inverse_energy_lower': lower, 'inverse_energy_upper': upper,
            'tangent_response': b,
            'minimum_tangent_energy_lower': b*b/(4*upper),
            'minimum_tangent_energy_upper': b*b/(4*lower),
            'minimum_is_infimum_without_inverse_vector_domain': True,
            'sharp_shell_or_speed_inferred': False}


def regularized_response(mass, first, eta):
    n, h, e = [rational(x, name) for x, name in
               zip((mass, first, eta), ('mass', 'first', 'eta'))]
    if n <= 0 or h <= 0 or e <= 0:
        raise ValueError('positive spectral mass, first moment and regularizer required')
    lower = n*n/(h+e*n)
    return {'regularized_inverse_lower': lower,
            'regularized_inverse_upper': n/e,
            'inverse_square_lower': lower*lower/n,
            'filtered_vector_norm_squared_lower': e*e*lower*lower/n}


def residual_certificate(dual_value, residual_squared, eta):
    d, r, e = [rational(x, name) for x, name in
               zip((dual_value, residual_squared, eta), ('dual_value', 'residual_squared', 'eta'))]
    if r < 0 or e <= 0:
        raise ValueError('nonnegative squared residual and positive regularizer required')
    if d+r/e < 0:
        raise ValueError('residual certificate contradicts positivity of the inverse response')
    return {'regularized_inverse_lower': max(F(), d),
            'regularized_inverse_upper': d+r/e,
            'absolute_dual_defect_upper': r/e,
            'requires_ground_gap': False}


def filter_budget(eta, cutoff, step, energy_band, coefficient_error=F(0)):
    e, t, h, band, coefficients = [rational(x, name) for x, name in
        zip((eta, cutoff, step, energy_band, coefficient_error),
            ('eta', 'cutoff', 'step', 'energy_band', 'coefficient_error'))]
    if min(e, t, h) <= 0 or band < 0 or coefficients < 0:
        raise ValueError('positive regularizer/cutoff/step and nonnegative errors required')
    # pi>2; the exact Fourier omitted mass is at most 2/(pi*eta*T).
    return {'omitted_time_mass_upper': min(F(1), 1/(e*t)),
            'midpoint_error_upper': band*h/2,
            'coefficient_error_upper': coefficients,
            'successful_operator_error_upper': min(F(1), 1/(e*t))+band*h/2+coefficients,
            'absolute_ground_origin_used_by_filter': False,
            'native_compilation_error_included': False}


def finite_readout_error(eta, energy_interval_radius, complete_l1_error,
                        input_trace_error=F(0), far_label_probability=F(0)):
    e, width, detector, state, leakage = [rational(x, n) for x, n in
        zip((eta, energy_interval_radius, complete_l1_error, input_trace_error, far_label_probability),
            ('eta', 'energy_interval_radius', 'complete_l1_error', 'input_trace_error',
             'far_label_probability'))]
    if e <= 0 or min(width, detector, state, leakage) < 0 or leakage > 1:
        raise ValueError('positive regularizer and nonnegative certified errors required')
    return {'bounded_score_error_upper': width/e+leakage+detector+state,
            'score_range': (F(0), F(1)),
            'far_label_probability_upper': leakage,
            'spectral_label_coupling_certificate_required': True,
            'all_reject_scores_zero': True,
            'renormalization_by_success_used': False}


def bounded_probe_error(eta, amplitude, complete_score_error):
    e, s, error = [rational(x, n) for x, n in
        zip((eta, amplitude, complete_score_error), ('eta', 'amplitude', 'complete_score_error'))]
    if e <= 0 or s == 0 or error < 0:
        raise ValueError('positive regularizer, nonzero amplitude and nonnegative error required')
    return {'regularized_inverse_error_upper': error/(e*s*s),
            'known_score_factor': e*s*s,
            'divides_by_observed_success': False,
            'continuum_limit_order': ['fixed_amplitude_and_regularizer_source_subsequence',
                                     'probe_amplitude_to_zero', 'regularizer_to_zero'],
            'effective_continuum_modulus_proved': False}


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) in (list, tuple):
        return [encode(x) for x in value]
    if type(value) is dict and all(type(k) is str for k in value):
        return {k: encode(v) for k, v in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('exact typed response receipt required')
