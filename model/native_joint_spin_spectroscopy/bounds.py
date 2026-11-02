"""Exact finite-resource and small-success-weight accounting."""
from fractions import Fraction as F


def conditional_bound(trace_error, success_lower):
    if type(trace_error) is not F or type(success_lower) is not F or not 0 <= trace_error < success_lower <= 1:
        raise ValueError('certified positive success and error below it are required')
    return {'approximate_success_lower': success_lower - trace_error / 2,
            'conditional_unhalved_trace_error': min(F(2), 2 * trace_error / success_lower)}


def leakage_upper(clock_size, bins, phase_margin):
    if type(clock_size) is not int or clock_size < 2 or type(bins) is not int or not 0 <= bins <= clock_size:
        raise ValueError('finite clock and counted event required')
    if type(phase_margin) is not F or not 0 < phase_margin <= 3:
        raise ValueError('certified circular separation in (0,3] required')
    # sin(delta/2)>=delta/pi>=delta/4 for 0<delta<=3<pi.
    return min(F(1), F(16 * bins, clock_size**2) / phase_margin**2)


def finite_program_error(clock_size, tau, coupling, field, edges, overlaps,
                         max_step, origin_error, pulses, precision_bits):
    if type(clock_size) is not int or clock_size < 2:
        raise ValueError('finite clock required')
    if any(type(x) is not F or x < 0 for x in (tau, coupling, field, max_step, origin_error)) or min(tau, coupling, max_step) <= 0:
        raise ValueError('exact nonnegative source data, positive time/coupling/step required')
    if any(type(x) is not int or x < 0 for x in (edges, overlaps, pulses, precision_bits)) or min(pulses, precision_bits) == 0:
        raise ValueError('nonnegative finite graph counts and positive pulse precision required')
    duration = (clock_size - 1) * tau
    commutator_half_sum = overlaps * coupling**2 + edges * coupling * field
    return {'clock_duration_sum': duration,
            'operator_splitting_error': commutator_half_sum * duration * max_step,
            'operator_origin_error': duration * origin_error,
            'operator_native_error': F(29 * pulses, 1 << precision_bits)}


def certificate():
    return {'conditional_example': {k: str(v) for k, v in conditional_bound(F(1, 10000), F(1, 100)).items()},
            'separated_clock_leakage_example': str(leakage_upper(64, 1, F(1))),
            'arithmetic_example_not_a_compiled_cube': {k: str(v) for k, v in finite_program_error(4, F(1, 2), F(1), F(1, 4), 12, 24, F(1, 8192), F(1, 10**12), 100000, 64).items()},
            'conditional_state_defined_without_positive_success': False,
            'packet_coefficient_error_included_in_complete_dilation': True}
