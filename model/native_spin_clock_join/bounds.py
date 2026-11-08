"""Exact budgets; trace errors use one half of the trace norm throughout."""
from fractions import Fraction as F


def rational(value, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('nonnegative exact rational required')
    return F(value)


def count(value, *, positive=False):
    if type(value) is not int or value < 0 or (positive and not value):
        raise ValueError('nonnegative integer required')
    return value


def ceil(value):
    return -(-value.numerator // value.denominator)


def ceil_log2(value):
    value = rational(value, positive=True)
    n = value.numerator.bit_length() - value.denominator.bit_length()
    power = F(1 << n) if n >= 0 else F(1, 1 << (-n))
    return n + (power < value)


def finite_history(field, sites, absolute_times, *, input_error=F(),
                   native_error=F(), product_error=F()):
    """Interleaved complete channels; controlled times use their block maximum.

    Each time is the nonnegative absolute duration for one evolution slot.
    A conjugation's forward and backward evolutions are separate slots.
    Field/reference inputs agree apart from the declared input halftrace error.
    """
    b = rational(field)
    n = count(sites, positive=True)
    times = tuple(rational(t) for t in absolute_times)
    errors = tuple(rational(e) for e in (input_error, native_error, product_error))
    if any(e > 1 for e in errors):
        raise ValueError('individual halftrace errors must be at most one')
    total_time = sum(times, F())
    clock_error = b * n * total_time / 2
    total = clock_error + sum(errors, F())
    return {'sites': n, 'field': b, 'absolute_duration': total_time,
            'evolution_slots': len(times), 'clock_halftrace_upper': clock_error,
            'input_halftrace_upper': errors[0], 'native_halftrace_upper': errors[1],
            'product_halftrace_upper': errors[2],
            'complete_halftrace_upper': min(F(1), total),
            'unclipped_halftrace_upper': total}


def ground_energy_bound(field, sites, bond_count, coupling=F(1), *,
                        prepared_halftrace_error=F()):
    """H_J energy above its own ground, for an approximate H_B ground.

    H_J is a sum of negative singlet projectors, so its spectral range is
    at most J times the number of bonds. No H_B-to-H_J vacuum swap occurs.
    """
    b, j = rational(field), rational(coupling, positive=True)
    n, edges = count(sites, positive=True), count(bond_count)
    eps = rational(prepared_halftrace_error)
    if eps > 1:
        raise ValueError('halftrace error must be at most one')
    return {'exact_HB_ground_HJ_energy_excess_upper': b*n,
            'prepared_HJ_energy_excess_upper': b*n+j*edges*eps,
            'prepared_HJ_energy_excess_lower': F(),
            'HJ_variance_for_exact_HB_ground_upper': (b*n/2)**2,
            'HJ_vacuum_selected': False}


def locality_radius(support_radius, support_sites, coupling, duration,
                    operator_error):
    """Certify the inherited spin LR cutoff using logarithmic integer tests.

    B_R <= 6 J s (2R+1)^2 T 2^{-(R-r)} exp(48JT).
    Onsite selection fields do not change this interaction-picture bound.
    """
    r = count(support_radius)
    s = count(support_sites, positive=True)
    if s > (2*r+1)**3:
        raise ValueError('support exceeds its cube')
    j, t = rational(coupling, positive=True), rational(duration)
    eps = rational(operator_error, positive=True)
    if not t:
        return {'radius': r+1, 'operator_cutoff_upper': F(),
                'integer_exponent_upper': 0, 'log_prefactor_upper': 0}
    exponent = ceil(96*j*t)
    prefactor = ceil_log2(6*j*s*t/eps)
    base = exponent+prefactor
    radius = max(16, r+max(base, 0)+1)
    while radius-r < base+2*ceil_log2(F(2*radius+1)):
        radius = r+base+2*ceil_log2(F(2*radius+1))+1
    return {'radius': radius, 'operator_cutoff_upper': eps,
            'integer_exponent_upper': exponent,
            'log_prefactor_upper': prefactor}


def local_clock_plan(*, support_radius, support_sites, coupling, duration,
                     interventions, error):
    """Choose a positive rational field for two local infinite histories.

    The initial state is identical for both clocks. A finite field-selected
    ground's distance from an infinite ordered marginal is a separate input.
    Global controlled-energy probes must use finite_history on their full box.
    """
    m = count(interventions)
    eps = rational(error, positive=True)
    if eps >= 1:
        raise ValueError('error must lie strictly between zero and one')
    t = rational(duration)
    op_error = eps/(2*(2*m+1))
    local = locality_radius(support_radius, support_sites, coupling, t, op_error)
    n = (2*local['radius']+1)**3
    field = eps/(n*t) if t else F(1)
    boundary = (2*m+1)*local['operator_cutoff_upper']
    clock = field*n*t/2
    if boundary+clock > eps:
        raise ArithmeticError('local clock budget exceeded')
    return dict(local, sites=n, field_upper=field, duration=t,
                interventions=m, two_cutoffs_halftrace_upper=boundary,
                finite_clock_halftrace_upper=clock,
                total_halftrace_upper=boundary+clock,
                infinite_ground_marginal_error_supplied=False,
                global_controlled_spectrometer_covered=False)


def spin_port_error(pulses, fraction_bits):
    """The declared dyadic integrated-port implementation, complete failures."""
    n, bits = count(pulses), count(fraction_bits)
    return min(F(1), F(29*n, 1 << bits))


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: encode(item) for key, item in value.items()}
    if type(value) in (list, tuple):
        return [encode(item) for item in value]
    return value
