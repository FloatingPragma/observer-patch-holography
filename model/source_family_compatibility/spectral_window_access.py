"""Complete source Ramsey records and certified original Lorentzian windows."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import sympy as s

from .local_matter import boundary_upper
from .local_spectral import dimer_source


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('an exact rational is required')
    value = F(value)
    if value < 0 or (positive and not value):
        raise ValueError('a nonnegative exact rational is required')
    return value


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('an integer in the admitted domain is required')
    return value


def interval(value, *, probability=False):
    if type(value) not in (tuple, list) or len(value) != 2:
        raise ValueError('a two-endpoint exact interval is required')
    if any(type(x) not in (int, F) for x in value):
        raise ValueError('exact interval endpoints are required')
    lo, hi = map(F, value)
    if lo > hi or (probability and not 0 <= lo <= hi <= 1):
        raise ValueError('an ordered interval in the admitted range is required')
    return lo, hi


def geometric_times(eta, step, count):
    eta = rational(eta, positive=True); step = rational(step, positive=True)
    count = integer(count, 1)
    r = 1/(1+eta*step)
    return {'ratio': r, 'weights': tuple((1-r)*r**j for j in range(count)),
            'dummy_weight': r**count, 'duration': (count-1)*step,
            'complete_records': 4*count+2}


def time_error(eta, step, count, second_moment_upper, first_moment_upper):
    second = rational(second_moment_upper)
    first = rational(first_moment_upper)
    if first*first < second:
        raise ValueError('the first-moment allowance must dominate the second-moment square root')
    times = geometric_times(eta, step, count)
    return {'Euler_Lorentzian_error': rational(step)*first/2,
            'discarded_time_score_error': times['dummy_weight'],
            'original_score_error': min(F(2), rational(step)*first/2+times['dummy_weight']),
            'dummy_records_retained': True, 'UV_alias_cut_required': False}


def cayley_service_bound(duration, step, hamiltonian_norm, subdivisions):
    t = rational(duration); dt = rational(step, positive=True)
    h = rational(hamiltonian_norm); n = integer(subdivisions, 1)
    unitary = t*dt*dt*h**3/(12*n*n)
    return {'unitary_operator_error': min(F(2), unitary),
            'complete_output_halftrace_error': min(F(1), 2*unitary)}


def source_second_moment(coupling, contact, degree=6, *, local_mode=False):
    c = rational(coupling); u = rational(contact); d = integer(degree)
    # The local identity is valid on every ground; the packet identity consumes
    # the inherited invariant full-ground CAR reference.
    return (u*u+d*c*c)/4 if local_mode else (u*u+(d*c)**2)/4


def complete_score_interval(measured, eta, step, count, second_moment_upper,
                            first_moment_upper, radius, support_radius,
                            support_sites, coupling, *, packet_error=0,
                            marginal_error=0, preparation_error=0,
                            service_error=0, repair_error=0):
    lo, hi = interval(measured)
    if lo < -1 or hi > 1:
        raise ValueError('the complete signed score must lie in [-1,1]')
    packet = rational(packet_error)
    errors = tuple(rational(x) for x in (marginal_error, preparation_error,
                                         service_error, repair_error))
    if any(x > 1 for x in errors) or packet > 2:
        raise ValueError('complete half-trace and unit-packet errors outside their domains')
    clock = geometric_times(eta, step, count)
    te = time_error(eta, step, count, second_moment_upper, first_moment_upper)
    boundary = boundary_upper(radius, support_radius, support_sites,
                              clock['duration'], rational(coupling))
    complete = min(F(1), boundary+sum(errors, F()))
    error = 2*packet+te['original_score_error']+2*complete
    lower, upper = max(F(), lo-error), min(F(1), hi+error)
    if lower > upper:
        raise ArithmeticError('the supplied score and error bounds are inconsistent')
    return {'original_Lorentzian_interval': (lower, upper),
            'original_score_error': error,
            'complete_output_halftrace_error': complete,
            'boundary_operator_error': boundary,
            'time_error': te, 'complete_records': clock['complete_records'],
            'conditional_normalization_applied': False,
            'same_ground_marginal_required': True,
            'packet_is_a_finite_cube_momentum_eigenvector': False}


def window_from_score(score, eta, energy_upper):
    lo, hi = interval(score, probability=True)
    eta = rational(eta, positive=True); w = rational(energy_upper, positive=True)
    outside = eta*eta/(eta*eta+w*w)
    return {'source_energy_window': (F(), w),
            'original_total_charge_weight_lower': max(F(), (lo-outside)/(1-outside)),
            'outside_Lorentzian_score_upper': outside,
            'zero_energy_component_allowed': True,
            'separate_charge_equality_assumed': False,
            'physical_mass_assigned': None}


def pressure_score_interval(side, coupling, field, energy_zero, energy_plus,
                            energy_minus, spans, *, stride=2):
    """A single algebraic even-score diagnostic on actual translation cells."""
    side = integer(side, 2); stride = integer(stride, 1)
    c = rational(coupling); h = rational(field, positive=True)
    if side % stride or side % 2 or type(spans) is not tuple or len(spans) != 3:
        raise ValueError('an even cube tiled by the specified translation cells is required')
    if any(type(x) is not int or not 0 <= x < side for x in spans):
        raise ValueError('the local observable spans must fit in the diagnostic cube')
    z = interval(energy_zero); plus = interval(energy_plus); minus = interval(energy_minus)
    cell_count = (side//stride)**3
    contained = 1
    for span in spans:
        contained *= (side-1-span)//stride+1
    omitted = 1-F(contained, cell_count)
    density_inverse = stride**3
    boundary = F(3*density_inverse, side)*c/h
    lo = max(F(-1), F(density_inverse, side**3)*(plus[0]-z[1])/h-boundary-omitted)
    hi = min(F(1), F(density_inverse, side**3)*(z[1]-minus[0])/h+boundary+omitted)
    if lo > hi:
        raise ArithmeticError('the pressure enclosures contradict the ground variational bound')
    return {'score_interval': (lo, hi),
            'translation_stride': stride, 'perturbation_density': F(1, stride**3),
            'omitted_translate_fraction': omitted,
            'surface_secant_error': boundary,
            'complete_observable_norm_upper': F(1),
            'isolated_cut_ground_substituted': False,
            'infinite_ground_prepared': False}


def pressure_window_certificate(side, coupling, field, energy_zero, energy_plus,
                                energy_minus, spans, *, stride, eta, step,
                                count, second_moment_upper, first_moment_upper,
                                radius, support_radius, support_sites,
                                clock_service_error, energy_upper,
                                packet_error=0):
    """Join one certified ground pressure diagnostic to an original window.

    The three finite energy intervals must enclose the exact diagnostic
    matrices for this same algebraic score observable. No isolated finite
    ground or measured score is accepted as a replacement for those bounds.
    """
    pressure = pressure_score_interval(side, coupling, field, energy_zero,
        energy_plus, energy_minus, spans, stride=stride)
    score = complete_score_interval(pressure['score_interval'], eta, step,
        count, second_moment_upper, first_moment_upper, radius,
        support_radius, support_sites, coupling, packet_error=packet_error,
        service_error=clock_service_error)
    score['same_ground_marginal_required'] = False
    score['ground_score_certified_by_pressure'] = True
    score['a_finite_prepared_state_marginal_has_been_identified'] = False
    window = window_from_score(score['original_Lorentzian_interval'], eta, energy_upper)
    return {'pressure': pressure, 'score': score, 'window': window,
            'strict_positive_original_weight_certified':
                window['original_total_charge_weight_lower'] > 0,
            'thermodynamic_ground_preparation_claimed': False}


def default_terminating_instance():
    """An actual all-ground positive fixed-source target for the pressure search."""
    c, u, eta, dt, m = F(1, 8), F(1), F(1), F(1, 256), 2048
    second = source_second_moment(c, u, local_mode=True)
    if second != F(35, 128):
        raise ArithmeticError('actual local CAR second moment differs')
    te = time_error(eta, dt, m, second, F(3, 4))
    if te['original_score_error'] >= F(1, 256):
        raise ArithmeticError('exact geometric time budget failed')
    times = geometric_times(eta, dt, m)
    radius = max(1, (96*c*times['duration']).__ceil__())
    while boundary_upper(radius, 0, 1, times['duration'], c) > F(1, 512):
        radius += 1
    sites = (2*radius+1)**3
    bonds = 3*(2*radius)*(2*radius+1)**2
    hnorm = c*bonds+u*sites/4
    need = times['duration']*dt*dt*hnorm**3*F(512, 6)
    n = isqrt(need.numerator//need.denominator)+1
    service = cayley_service_bound(times['duration'], dt, hnorm, n)
    if service['complete_output_halftrace_error'] > F(1, 512):
        raise ArithmeticError('algebraic source clock subdivision budget failed')
    # Keep the exact 2048th-power admission, but publish a short outward
    # rational envelope. Its exact numerator is thousands of digits long.
    grid = 10**12
    time_upper = F((te['original_score_error']*grid).__ceil__(), grid)
    if not te['original_score_error'] <= time_upper < F(1, 256):
        raise ArithmeticError('compact outward time budget failed')
    error = time_upper+F(1, 128)
    jensen = 1/(1+second)
    if not error < F(3, 256) < F(1, 64) or jensen-error <= F(3, 4):
        raise ArithmeticError('uniform pressure-stopping margin failed')
    window = window_from_score((F(47, 64), F(1)), 1, 1)
    return {'source_c': c, 'source_U': u, 'local_second_moment': second,
            'exact_Jensen_Lorentzian_lower': jensen,
            'analytic_Markov_original_window_lower': 1-second,
            'operational_floor_improves_analytic_moment_bound': False,
            'eta': eta, 'time_step': dt, 'time_count': m,
            'complete_records': times['complete_records'],
            'time_error_upper': time_upper,
            'cut_radius': radius, 'cut_sites': sites,
            'cut_Hamiltonian_norm_upper': hnorm,
            'Cayley_subdivisions_per_time_step': n,
            'algebraic_complete_clock_error': service['complete_output_halftrace_error'],
            'uniform_score_lower': jensen-error,
            'pressure_search_strict_target': F(3, 4),
            'pressure_search_termination_proved': True,
            'pressure_search_executed': False,
            'published_strict_original_window_floor': window['original_total_charge_weight_lower'],
            'window': window['source_energy_window'],
            'native_large_clock_tape_executed': False,
            'shrinking_window_uniform_weight_proved': False}


def _simple(matrix):
    return s.ImmutableMatrix(matrix.applyfunc(s.simplify))


def _radical_interval(value, scale=10**10):
    """Outward exact interval for the source field Q(sqrt(17))."""
    value = s.radsimp(s.expand(value))
    root = s.sqrt(17)
    a = s.expand(value).coeff(root, 0)
    b = s.expand(value).coeff(root, 1)
    if not a.is_Rational or not b.is_Rational or s.simplify(value-a-b*root) != 0:
        raise ArithmeticError('source expectation left its exact radical field')
    a, b = F(int(a.p), int(a.q)), F(int(b.p), int(b.q))
    grid = 10**32
    lower = F(isqrt(17*grid*grid), grid); upper = lower+F(1, grid)
    if not lower*lower <= 17 <= upper*upper:
        raise ArithmeticError('exact square-root enclosure failed')
    lo, hi = (a+b*lower, a+b*upper) if b >= 0 else (a+b*upper, a+b*lower)
    return F((lo*scale).__floor__(), scale), F((hi*scale).__ceil__(), scale)


@lru_cache(maxsize=1)
def _finite_control():
    h, rho, a, projections = dimer_source()
    shift = h-s.Rational(1, 2)*s.eye(16)
    eta, dt, count = F(1), F(1, 4), 16
    times = geometric_times(eta, dt, count)
    v = _simple((s.eye(16)-s.I*s.Rational(1, 8)*shift)
                *(s.eye(16)+s.I*s.Rational(1, 8)*shift).inv())
    if _simple(v.H*v) != s.eye(16):
        raise ArithmeticError('finite source Cayley clock is not unitary')
    power = s.eye(16); effect = s.zeros(16); score = s.zeros(16); rows = []
    for j, weight in enumerate(times['weights']):
        w = s.Rational(weight.numerator, weight.denominator)
        for label, op in (('addition', a.H), ('removal', a)):
            evolved = power*op*power.H
            for sign in (1, -1):
                k = (op+sign*evolved)/2
                e = _simple(w*k.H*k)
                effect += e; score += sign*e
                probability = s.simplify(s.trace(rho*e))
                rows.append({'time_index': j, 'charge': label, 'sign': sign,
                             'original_probability_interval': _radical_interval(probability)})
        power = _simple(v*power)
    w = s.Rational(times['dummy_weight'].numerator, times['dummy_weight'].denominator)
    for label, op in (('addition', a.H), ('removal', a)):
        e = _simple(w*op.H*op); effect += e
        rows.append({'time_index': None, 'charge': label, 'sign': 0,
                     'original_probability_interval': _radical_interval(s.trace(rho*e))})
    if _simple(effect) != s.eye(16) or score != score.H:
        raise ArithmeticError('complete source instrument failed')
    score_value = s.simplify(s.trace(rho*score))
    exact = sum(s.simplify(s.trace(rho*op.H*p*op)/(1+energy**2))
                for op in (a.H, a) for energy, p in projections)
    exact = s.simplify(exact)
    score_box = _radical_interval(score_value); exact_box = _radical_interval(exact)
    second = source_second_moment(F(1, 8), 1, degree=1, local_mode=True)
    te = time_error(eta, dt, count, second, F(9, 16))
    service = cayley_service_bound(times['duration'], dt, F(5, 8), 1)
    error = te['original_score_error']+2*service['complete_output_halftrace_error']
    if max(abs(score_box[0]-exact_box[1]), abs(score_box[1]-exact_box[0])) > error:
        raise ArithmeticError('actual full-source score exceeds its certified error')
    measured_q = (max(F(), score_box[0]-error), min(F(1), score_box[1]+error))
    window = window_from_score(measured_q, 1, 1)
    if window['original_total_charge_weight_lower'] <= F(1, 3):
        raise ArithmeticError('finite source original-window floor failed')
    return {'source_sites': 2, 'source_Fock_rows': 16,
            'source_c': F(1, 8), 'source_U': F(1), 'eta': eta,
            'time_step': dt, 'time_count': count,
            'all_records': rows, 'complete_record_count': len(rows),
            'effect_sum_equals_identity_exactly': True,
            'algebraic_clock_unitarity_exact': True,
            'source_score_interval': score_box,
            'exact_source_Lorentzian_interval': exact_box,
            'certified_score_error': error,
            'original_Lorentzian_interval_from_records': measured_q,
            'source_window': window,
            'all_source_occupation_sectors_retained': True,
            'thermodynamic_marginal_substitution_made': False,
            'native_pulse_execution_claimed': False,
            'ideal_algebraic_complete_instrument_evaluated': True}


def _detach(value):
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is dict:
        return {key: _detach(x) for key, x in value.items()}
    if type(value) is tuple:
        return tuple(_detach(x) for x in value)
    if type(value) is list:
        return [_detach(x) for x in value]
    return value


def finite_control():
    return _detach(_finite_control())


def certificate():
    return {'complete_Ramsey_Lorentzian_readout': finite_control(),
            'actual_fixed_source_terminating_pressure_instance': default_terminating_instance(),
            'physical_species_assigned': False,
            'physical_mass_assigned': None,
            'uniform_shrinking_window_weight_proved': False,
            'effective_growing_ground_marginal_modulus_proved': False}
