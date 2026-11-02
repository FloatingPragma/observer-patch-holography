"""Success, conditioning, quadrature and explicit universal native budgets."""
from fractions import Fraction as F
from math import isqrt
from .filters import exact, fourier_bounds, goldstone_filter


def ceil(x):
    x = exact(x)
    return -(-x.numerator//x.denominator)


def ceil_log2(x):
    x = exact(x, positive=True)
    n = x.numerator.bit_length()-x.denominator.bit_length()
    power = F(1 << n) if n >= 0 else F(1, 1 << (-n))
    return n+(power < x)


def locality_radius(support_radius, coupling, duration, error):
    """Certify the inherited LR bound without constructing 2**(huge time)."""
    if type(support_radius) is not int or support_radius < 0:
        raise ValueError('nonnegative integer support radius required')
    j, t, eps = (exact(x, positive=True) for x in (coupling, duration, error))
    base = ceil(96*j*t)+ceil_log2(6*j*t/eps)
    radius = max(16, support_radius+max(base, 0)+1)
    while radius-support_radius < base+2*ceil_log2(F(2*radius+1)):
        radius = support_radius+base+2*ceil_log2(F(2*radius+1))+1
    return {'radius': radius, 'support_radius': support_radius,
            'time_exponent_upper': ceil(96*j*t),
            'log2_prefactor_over_error_upper': ceil_log2(6*j*t/eps),
            'certified_Pauli_operator_cutoff_error': eps,
            'applies_to_each_controlled_site_block': True}


def sqrt_upper(x):
    x = exact(x)
    if x < 0:
        raise ValueError('nonnegative square root required')
    scale = 1 << 64
    n = (x.numerator*scale*scale+x.denominator-1)//x.denominator
    k = isqrt(n)
    return F(k+(k*k < n), scale)


def packet_success(p, envelope_bound, m=F(1, 1000)):
    p, lam, m = (exact(x, positive=True) for x in (p, envelope_bound, m))
    if p > 1:
        raise ValueError('the certified low-momentum regime requires p<=1')
    filt = goldstone_filter(m)
    normalization = filt['kernel_L1_upper']+1
    q = m*m*p**4/(96*lam*lam*normalization**2)
    if q > 1:
        raise ValueError('the supplied envelope bound is incompatible with a probability')
    return {'p': p, 'envelope_L1_sampling_bound': lam,
            'time_LCU_normalization': normalization,
            'ideal_success_lower': q,
            'independent_fresh_input_attempts_mean_upper': ceil(1/q),
            'packet_weight_lower': m*m*p/96}


def conditioned_errors(success_lower, success_amplitude_error, input_trace_error=F()):
    q = exact(success_lower, positive=True)
    eta = exact(success_amplitude_error)
    delta = exact(input_trace_error)
    if q > 1 or min(eta, delta) < 0:
        raise ValueError('valid probability and nonnegative errors required')
    # A rational lower square-root bound avoids inflating the denominator.
    root_upper = sqrt_upper(q)
    root_lower = q/root_upper
    if eta >= root_lower or delta >= (root_lower-eta)**2:
        raise ValueError('positive actual success is not certified')
    pure_lower = (root_lower-eta)**2
    return {'actual_success_lower': pure_lower-delta,
            'conditional_trace_error_upper': min(F(2), 4*eta/root_lower+2*delta/pure_lower),
            'convention': 'unhalved trace norm; fresh supplied input state'}


def quadrature_plan(a, b, d, p, error, field_over_coupling=F()):
    p, eps = exact(p, positive=True), exact(error, positive=True)
    ratio = exact(field_over_coupling)
    if ratio < 0:
        raise ValueError('the selection-field ratio is nonnegative')
    if p > 1 or eps >= 1:
        raise ValueError('require p<=1 and an error below one')
    f = fourier_bounds(a, b, d)
    # Error of the UNNORMALIZED time-filter operator on a normalized site block.
    cutoff = ceil(sqrt_upper(3*f['tail_numerator']/eps))
    derivative = f['kernel_derivative_upper']
    c = f['kernel_L1_upper']
    max_step = eps/(3*(derivative*cutoff/2+(6+ratio/2)*c/p))
    cells = max(1, ceil(2*cutoff/max_step))
    step = F(2*cutoff, cells)
    tail = f['tail_numerator']/cutoff**2
    coefficient = derivative*cutoff*step/2
    evolution = (6+ratio/2)*c*step/p
    # Each complex cell coefficient is enclosed/rounded with total error eps/3.
    rounding = eps/3
    return {'time_cutoff': cutoff, 'cells': cells, 'step': step,
            'maximum_physical_time_per_inverse_J': F(cutoff)/p,
            'kernel_tail_error': tail, 'midpoint_coefficient_error': coefficient,
            'Heisenberg_midpoint_error': evolution,
            'coefficient_rounding_allowance': rounding,
            'total_time_filter_error_upper': tail+coefficient+evolution+rounding,
            'coefficient_L1_upper': c+coefficient+rounding,
            'per_coefficient_complex_error_allowance': rounding/cells}


def native_plan(sites, edges, overlap_pairs, field_terms, labels, duration,
                coupling, field, product_error, native_error):
    integers = (sites, edges, overlap_pairs, field_terms, labels)
    if any(type(x) is not int or x < 0 for x in integers) or sites < 1 or labels < 1:
        raise ValueError('finite nonnegative graph counts and positive carriers required')
    t, j, pe, ne = (exact(x, positive=True) for x in (duration, coupling, product_error, native_error))
    b = exact(field)
    if b < 0 or pe >= 1 or ne >= 1:
        raise ValueError('nonnegative field and errors below one required')
    q = max(1, (labels-1).bit_length())
    dimension = 1 << q
    commutators = overlap_pairs*j*j+edges*j*b
    steps = max(1, ceil(2*commutators*t*t/pe))
    pauli_terms = 2*steps*(4*edges+field_terms)+2
    # Explicit controlled-pattern expansion; this deliberately does not hide
    # the exponential-in-ancilla, polynomial-in-labels SELECT overhead.
    prep = 4*dimension*dimension*(8*q+1)
    select = labels*dimension*pauli_terms*(8*(2*sites+q)+1)
    pulses = prep+select
    bits = 0
    while 30*pulses*F(1, 1 << bits) > ne:
        bits += 1
    return {'label_qubits': q, 'padded_label_dimension': dimension,
            'Trotter_steps_per_evolution': steps,
            'maximum_conjugation_product_error': 2*commutators*t*t/steps,
            'PREP_and_inverse_native_pulses_upper': prep,
            'SELECT_native_pulses_upper': select,
            'filter_native_pulses_upper': pulses, 'port_fraction_bits': bits,
            'angle_synthesis_error_upper': pulses*F(1, 1 << bits),
            'native_port_isometry_error_upper': 29*pulses*F(1, 1 << bits),
            'native_complete_isometry_error_upper': 30*pulses*F(1, 1 << bits),
            'input_ground_preparation_cost_included': False,
            'literal_large_tape_executed': False}
