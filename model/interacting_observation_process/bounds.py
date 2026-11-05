"""Exact budgets for adaptive instruments interleaved with native memory steps."""
from fractions import Fraction as F


COMPACT_NORM = F(7869, 32)
MEMORY_CONTRACTION = F(1, 12)


def rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(name+' must be an exact rational')
    return F(value)


def segment(norm, duration, contraction, steps, *, stationary_product=False):
    m, t, p = (rational(x, name) for x, name in
               ((norm, 'norm'), (duration, 'duration'), (contraction, 'contraction')))
    if m < 0 or t < 0 or not 0 <= p < 1 or type(steps) is not int or steps < 1:
        raise ValueError('nonnegative bounds, contraction below one and positive steps required')
    if type(stationary_product) is not bool:
        raise ValueError('stationary product preparation must be explicit')
    transient = F() if stationary_product else 4*m*t/(steps*(1-p))
    quadratic = (8/(1-p)+4)*m*m*t*t/steps
    return {'unread_trace_norm_upper': transient+quadratic,
            'complete_record_trace_norm_upper_squared': min(F(4), 2*(transient+quadratic)),
            'transient': transient, 'quadratic': quadratic}


def schedule(durations, error=F(1, 100), *, norm=COMPACT_NORM,
             contraction=MEMORY_CONTRACTION, first_stationary_product=False):
    """Equal error allocations; only an explicitly prepared first interval may improve."""
    if type(durations) not in (tuple, list) or not durations or len(durations) > 1024:
        raise ValueError('one through 1024 interval duration bounds are required')
    if type(first_stationary_product) is not bool:
        raise ValueError('the first preparation flag must be Boolean')
    epsilon = rational(error, 'error')
    if not 0 < epsilon <= 2:
        raise ValueError('the complete-state error must lie in (0,2]')
    allocated = epsilon/len(durations)
    result = []
    for index, duration in enumerate(durations):
        stationary = first_stationary_product and index == 0
        initial = segment(norm, duration, contraction, 1, stationary_product=stationary)
        required = 2*initial['unread_trace_norm_upper']/(allocated*allocated)
        steps = max(1, -(-required.numerator//required.denominator))
        bound = segment(norm, duration, contraction, steps, stationary_product=stationary)
        if bound['complete_record_trace_norm_upper_squared'] > allocated*allocated:
            raise ArithmeticError('the exact per-interval error allocation failed')
        result.append({'duration_upper': str(F(duration)), 'native_steps': steps,
                       'allocated_complete_trace_error': str(allocated),
                       'complete_record_error_squared_upper': str(bound['complete_record_trace_norm_upper_squared']),
                       'stationary_product_used': stationary})
    return {'intervals': result, 'core_norm_upper': str(F(norm)),
            'memory_contraction_upper': str(F(contraction)),
            'complete_adaptive_output_trace_error_upper': str(epsilon),
            'total_native_steps': sum(row['native_steps'] for row in result),
            'intermediate_memory_reset': False,
            'arbitrary_untouched_reference': True,
            'all_native_and_instrument_records_retained': True,
            'uniform_in_exterior_volume': True,
            'uniform_in_lattice_spacing': False}


def equal_grid_squared(norm, duration, contraction, intervals, total_steps):
    """Cauchy-free squared sum for equal intervals and equal integer step counts."""
    if (type(intervals) is not int or intervals < 1 or type(total_steps) is not int
            or total_steps < intervals or total_steps % intervals):
        raise ValueError('equal positive integral step allocations are required')
    one = segment(norm, rational(duration, 'duration')/intervals, contraction,
                  total_steps//intervals)
    return min(F(4), intervals*intervals*one['complete_record_trace_norm_upper_squared'])


def certificate():
    m = 90*F(7, 240)+90*(1+1)+8100*F(1, 128)
    if m != COMPACT_NORM:
        raise ArithmeticError('the actual compact source norm bound changed')
    plan = schedule((F(1, 3),)*3)
    # A stationary product can be used at the first boundary only. The
    # source can correlate memory and field before every later intervention.
    prepared = schedule((F(1, 3),)*3, first_stationary_product=True)
    if [r['stationary_product_used'] for r in prepared['intervals']] != [True, False, False]:
        raise ArithmeticError('stationary-product preparation was restarted')
    return {'compact_annihilators': 180, 'compact_r_upper': '7/240',
            'feedback_gain_upper': '1', 'quartic_coupling': '1/128',
            'core_norm_upper': str(m), 'memory_contraction_upper': '1/12',
            'three_interval_one_percent_schedule': plan,
            'initial_stationary_product_control': prepared,
            'equal_grid_complete_error_squared_before_cap':
                '8*M*T*m^2/(N*(1-p))+2*(8/(1-p)+4)*M^2*T^2*m/N',
            'sufficient_growing_intervention_condition_at_fixed_T_M_p': 'm=o(sqrt(N))',
            'unlimited_microstep_adaptation_claimed': False}
