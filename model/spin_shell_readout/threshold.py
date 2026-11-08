"""Original-weight spectral edges, without assuming a one-particle pole."""
from fractions import Fraction as F

from model.spin_dynamical_memory.closure import native_dimer_control
from model.spin_multipacket_comparison.bounds import local_band, filter_bounds
from .bounds import ceil


def rational(value, name, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(name+' must be an exact nonnegative rational')
    return F(value)


def interval(value):
    if type(value) not in (tuple, list) or len(value) != 2:
        raise ValueError('a probability interval is required')
    lo, hi = (rational(x, 'probability') for x in value)
    if not lo <= hi <= 1:
        raise ValueError('ordered probabilities in [0,1] required')
    return lo, hi


def bracket(x_low, cdf_low, x_high, cdf_high, epsilon, amplitude, power):
    """Conditional on nu([theta,theta+e]) >= amplitude*e**power.

    The caller must certify that premise for every 0<e<=epsilon. Counts
    alone cannot supply it. Coordinates need not be nonnegative.
    """
    if type(x_low) not in (int, F) or type(x_high) not in (int, F) or x_low >= x_high:
        raise ValueError('ordered exact coordinates required')
    if type(power) is not int or power < 1:
        raise ValueError('positive integer tail power required')
    eps = rational(epsilon, 'tail radius', True)
    a = rational(amplitude, 'tail amplitude', True)
    low, high = interval(cdf_low), interval(cdf_high)
    if low[0] > high[1]:
        raise ValueError('CDF intervals contradict monotonicity')
    if a*eps**power > 1:
        raise ValueError('tail premise exceeds total original probability')
    if low[1] >= a*eps**power or high[0] <= 0:
        raise ValueError('strict original-probability evidence is insufficient')
    return {'threshold_lower_strict':F(x_low)-eps,
            'threshold_upper':F(x_high),
            'tail_radius':eps, 'tail_probability_floor':a*eps**power,
            'tail_modulus_is_required_source_premise':True,
            'pole_or_nonzero_atom_required':False,
            'conditional_success_renormalization_used':False}


def sampled_cdf(trials, success, counts, epsilon, failure_probability,
                complete_output_error=F()):
    """Simultaneous Chebyshev intervals for a predeclared ordered scan.

    All trials (including preparation failure and overflow) are counted.
    Trials are independent complete experiments. Systematic complete-output
    half-trace error is added to each event probability, not its success ratio.
    """
    if type(trials) is not int or trials < 1 or type(success) is not int:
        raise ValueError('positive trial count and integer successes required')
    if type(counts) not in (tuple, list) or not counts:
        raise ValueError('a nonempty predeclared CDF count list is required')
    if any(type(c) is not int for c in counts) or list(counts) != sorted(counts):
        raise ValueError('ordered integer CDF counts required')
    if not 0 <= counts[0] <= counts[-1] <= success <= trials:
        raise ValueError('complete CDF counts cannot exceed original successes')
    eps = rational(epsilon, 'sampling error', True)
    beta = rational(failure_probability, 'failure probability', True)
    error = rational(complete_output_error, 'complete-output error')
    if eps >= 1 or beta >= 1 or error > 1:
        raise ValueError('probability budgets outside their allowed range')
    required = ceil(F(len(counts), 1)/(4*beta*eps**2))
    if trials < required:
        raise ValueError('too few trials for the entire declared scan')
    return {'required_trials':required, 'total_trials':trials,
        'preparation_failure_count':trials-success,
        'success_above_last_threshold_count':success-counts[-1],
        'original_CDF_intervals':[(max(F(),F(c,trials)-eps-error),
                                  min(F(1),F(c,trials)+eps+error)) for c in counts],
        'simultaneous_failure_probability_upper':beta}


def non_atomic_control(level):
    if type(level) is not int or level < 1:
        raise ValueError('positive integer resolution level required')
    theta, weight, delta = F(9,16), F(1,5), F(1,2**level)
    eta = delta**2/100
    # Analytic subprobability CDF w*clamp(x-theta,0,1)^2, with no atom.
    report = bracket(theta-delta,(0,eta),theta+delta,
                     (weight*delta**2-eta,weight*delta**2+eta),delta,weight,2)
    return {'kind':'analytic_nonatomic_control_not_an_OPH_source_state',
        'threshold':theta, 'original_success_weight':weight,
        'complete_failure_weight':1-weight, 'atom_weight':F(),
        'original_weight_below_upper_probe':weight*delta**2,
        'CDF_error_budget':eta, 'bracket':report,
        'independent_trials_sufficient_for_two_CDFs':ceil(2/(4*F(1,100)*eta**2)),
        'large_trial_plan_executed':False, 'physical_mass':None}


def native_control():
    source = native_dimer_control()
    energies = [F(x) for x in source['original_success_energies']]
    weights = [F(x) for x in source['original_success_spectral_weights']]
    edge = min(energies)
    floor = sum(w for e,w in zip(energies,weights,strict=True) if e == edge)
    # For 0<epsilon<=1, the exact first atom proves nu([edge,edge+e])>=floor*e.
    low, high, eta, eps = F(281,250), F(563,500), F(1,10000), F(1,1000)
    cdf = lambda x: sum((w for e,w in zip(energies,weights,strict=True) if e <= x),F())
    measured = bracket(low,(max(F(),cdf(low)-eta),cdf(low)+eta),high,
                       (cdf(high)-eta,cdf(high)+eta),eps,floor,1)
    return {'actual_source_control':source,
        'source_prepared_response_edge':edge,
        'tail_amplitude':floor, 'tail_power':1, 'tail_radius_max':1,
        'complete_original_probability_bracket':measured,
        'additional_native_tape_executed':False,
        'continuum_threshold_modulus_proved':False,
        'spatial_momentum_or_physical_mass_assigned':False}


def gapless_charge_obstruction(target_gap=F(1), energy_conversion=F(1)):
    """Conditional whole-sector comparison, using the actual ordered-spin bound."""
    gap = rational(target_gap,'target same-charge energy gap',True)
    conversion = rational(energy_conversion,'fixed source-to-target energy factor',True)
    source = local_band()
    a,b = source['dimensionless_energy_window']
    filt = filter_bounds(a,b)
    mass = filt['time_filter_L1_upper']
    original_weight = source['each_charge_original_joint_band_probability_lower']/mass**2
    radius = min(F(1), gap/(4*conversion*b))
    energy_upper = conversion*2*b*radius
    if not energy_upper < gap:
        raise ArithmeticError('chosen source window did not lie below target gap')
    return {'source_family':'ordered-zero-field spin, compact positive-transfer creator',
        'source_global_internal_charges':[1,-1],
        'source_momentum_radius':radius,
        'target_energy_upper_of_source_event':energy_upper,
        'assumed_target_charge_gap':gap,
        'fixed_energy_conversion':conversion,
        'complete_dictionary_half_trace_error_lower':original_weight,
        'target_preserves_full_source_charge_sector_and_vacuum_energy':True,
        'electromagnetic_charge_identification_proved':False,
        'massive_subbranches_or_different_sources_excluded':False}


def certificate():
    return {'non_atomic_controls':[non_atomic_control(k) for k in (2,4,8)],
        'actual_native_finite_threshold':native_control(),
        'gapless_charge_target_comparison':gapless_charge_obstruction(),
        'infrared_threshold_route_is_a_conditional_readout_theorem':True,
        'continuum_infraparticle_or_measured_mass_identified':False}
