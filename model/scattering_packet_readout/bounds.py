"""Exact cutoffs and complete-output budgets for the source packet experiment."""
from fractions import Fraction as F
from math import comb, factorial

from ..contact_continuum.continuum import rational
from ..moving_contact_scattering.bounds import scattering_certificate


MAX_SPACING = F(1, 65536)


def ceiling(value):
    value = F(value)
    return -(-value.numerator//value.denominator)


def log2_ceiling(value):
    value = rational(value, 'positive logarithm argument', positive=True)
    if value <= 1:
        return 0
    power = max(0, value.numerator.bit_length()-value.denominator.bit_length())
    if value.denominator*(1 << power) < value.numerator:
        power += 1
    return power


def profile_certificate():
    integral = sum((F(2*(-1)**j*comb(6, j), 2*j+1) for j in range(7)), F())
    if integral != F(2048, 3003):
        raise ArithmeticError('packet normalization changed')
    return {'energy_profile': 'sqrt(3003/512)*(1-16*(s-1)^2)^3 on [3/4,5/4], zero elsewhere',
            'momentum_profile': 'product_j sqrt(3003/2048)*(1-P_j^2)^3 on [-1,1]^3, zero elsewhere',
            'energy_squared_norm': F(1), 'momentum_squared_norm': F(1),
            'profile_L1_bounds': [2, 48, 192],
            'normalized_total_momentum_derivative_L2_upper': 20,
            'all_eight_edge_valleys_retained': True,
            'preparation_uses_interaction_phase': False}


def analytic_certificate():
    a2 = MAX_SPACING**2
    if F(1, 9)-32*a2 <= F(1, 10):
        raise ArithmeticError('complex radial ray separation failed')
    if min((1-32*a2)/4-2*a2, (1-13*a2)/4-2*a2) <= F(1, 5):
        raise ArithmeticError('arc or exterior denominator separation failed')
    if F(9, 640) <= F(1, 100):
        raise ArithmeticError('positive source density floor failed')
    slope, curvature = F(512*44, 7), F(8192*44, 7)+4*256**2*F(44, 7)**2
    if not slope < 3220 or not curvature < 11000000:
        raise ArithmeticError('scattering derivative majorant failed')
    return {'complex_energy_disk': '|s-1|<=1/2',
            'complex_momentum_disks': '|K_j-K_j(real)|<=a, |K_j(real)|<=a',
            'ray_denominator_lower': '(x^2+a^2)/10',
            'arc_and_exterior_denominator_lower': F(1, 5),
            'analytic_Green_upper': '32*a^2',
            'Green_energy_derivative_upper': ['256*a^2', '4096*a^2'],
            'density_bounds': ['a^3/100', 'a^3/2'],
            'density_energy_derivative_upper': ['86*a^2', '1366*a^2'],
            'density_lattice_momentum_derivative_upper': '11*a',
            'scattering_energy_derivative_upper': ['3220/a', '11000000/a^2'],
            'unrounded_scattering_derivative_coefficients': [slope, curvature],
            'incoming_plus_outgoing_Cook_numerator_upper': '40000000/sqrt(a)',
            'Cook_operator_on_packet_error_upper': '160000000/(a^3*T)',
            'normalized_spatial_truncation_error_upper': '40000000/(a^5*R)',
            'finite_box_error_upper': '2^(2-L) for L>=288*T/a^2'}


def plan(spacing=MAX_SPACING, component_error=F(1, 100)):
    a = rational(spacing, 'spacing', positive=True)
    error = rational(component_error, 'component error', positive=True)
    if a > MAX_SPACING or error > F(1, 100):
        raise ValueError('requires a<=1/65536 and component error<=1/100')
    duration = ceiling(160000000/(a**3*error))
    support = ceiling(40000000/(a**5*error))
    padding = ceiling(288*duration/a**2)+log2_ceiling(4/error)
    side = 2*(support+padding)+1
    sites = side**3
    bits = 2*sites
    term_norm_sum = 10*sites/a**2
    steps = ceiling(6*duration**2*term_norm_sum**2/error)
    # The amplitude tree count is exponential in bits. Its logarithmic bound
    # is computed without allocating either that integer or the native tape.
    preparation_log = bits+(11*bits).bit_length()
    evolution_pulses = 768*steps*sites*(bits+1)+10
    pulse_log_upper = max(preparation_log, evolution_pulses.bit_length())+1
    port_bits = pulse_log_upper+log2_ceiling(58/error)
    phase = scattering_certificate(a, 1, 1, F(5, 4))['uniform_scattering_operator_error_upper']
    plus_upper = (1-F(49, 50)+phase)/2+7*error
    free_lower = F(1, 2)-7*error
    if plus_upper >= F(1, 8) or free_lower-plus_upper <= F(3, 10):
        raise ArithmeticError('complete recorded phase contrast failed')
    return {'spacing': a, 'gamma': 1, 'component_error': error,
            'absolute_echo_segment_duration': duration,
            'normalized_packet_support_radius': support,
            'box_padding': padding, 'box_side': side, 'box_sites': sites,
            'field_bits': bits, 'interference_bits': 1,
            'first_order_steps_per_echo_segment': steps,
            'native_pulse_count_formula': '11*M*2^M+768*n*V*(M+1)+10',
            'native_pulse_log2_upper': pulse_log_upper,
            'dyadic_integrated_port_fraction_bits': port_bits,
            'angle_approximation_error_allocation': error/2,
            'dyadic_port_execution_error_allocation': error/2,
            'complete_dilation_error_upper': 7*error,
            'continuum_phase_error_upper': phase,
            'valid_Y_plus_probability_upper': plus_upper,
            'free_valid_Y_plus_probability_lower': free_lower,
            'complete_law_event_separation_lower': free_lower-plus_upper,
            'giant_tape_or_exponential_pulse_integer_allocated': False,
            'fine_spacing_native_experiment_executed': False}


def sampling_certificate():
    positive_sum = sum((F(10**n, factorial(n)) for n in range(9)), F())
    if positive_sum <= 4000:
        raise ArithmeticError('positive exponential certificate failed')
    return {'independent_complete_histories_per_law': 2000,
            'number_of_laws': 2, 'simultaneous_event_probability_error': F(1, 20),
            'failure_upper': F(4)/positive_sum,
            'requested_failure_upper': F(1, 1000),
            'rejected_histories_counted': True,
            'rejected_states_recycled_as_fresh_inputs': False,
            'random_shot_records_fabricated': False}
