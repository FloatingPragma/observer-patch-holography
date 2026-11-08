"""Rational source-gap, energy and concentration bounds for one joint family."""
from fractions import Fraction as F
from math import ceil

from model.spin_canonical_field.bounds import sqrt_upper, energy_after, transport_error_squared
from model.spin_fixed_field_locality import bounds as locality


def rational(value, name, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(name+' must be an exact '+('positive' if positive else 'nonnegative')+' rational')
    return F(value)


def family(n, coupling=F(1)):
    if type(n) is not int or n < 2:
        raise ValueError('integer scale n>=2 required')
    j = rational(coupling, 'coupling', True)
    p, b = F(1,n*n), j/F(n**5)
    return {'n':n, 'p':p, 'coupling':j, 'field':b, 'gap':b/8,
            'field_over_J_p_squared':F(1,n), 'J_p_cubed_over_field':F(1,n),
            'physical_clock_per_scaled_unit':1/(j*p)}


def order_concentration(n, support_radius, supremum, coupling=F(1)):
    row = family(n,coupling)
    radius = rational(support_radius,'support radius',True)
    if radius < 1:
        raise ValueError('support radius must be at least one')
    sup = rational(supremum,'test supremum')
    volume = (2*radius+1)**3
    energy = 6*row['coupling']*volume*sup**2*row['p']**3
    variance = energy/row['gap']
    return {**row,'support_count_coefficient':volume,'test_supremum':sup,
            'ground_energy_form_upper':energy,'variance_upper':variance,
            'center_is_the_actual_positive_field_expectation':True,
            'magnetization_convergence_rate_supplied':False}


def bond_concentration(n, coefficient, l1_volume, coupling=F(1)):
    row = family(n,coupling)
    a = rational(coefficient,'coefficient'); v = rational(l1_volume,'l1 volume')
    j,b,p = row['coupling'],row['field'],row['p']
    energy = (99*j+3*b)*a*a*v*p**3
    return {**row,'coefficient':a,'l1_volume':v,
            'ground_energy_form_upper':energy,'variance_upper':energy/row['gap'],
            'simple_variance_upper':816*a*a*v/n,
            'selection_double_commutator_included':True,
            'mean_is_the_actual_source_bond_expectation':True}


def bond_transport(n, coefficient, l1_volume, scaled_time, scaled_energy,
                   phase=F(1), coupling=F(1)):
    row = bond_concentration(n,coefficient,l1_volume,coupling)
    t = rational(scaled_time,'absolute scaled time')
    e = rational(scaled_energy,'scaled energy')
    s = rational(phase,'absolute phase')
    bound = transport_error_squared(t/(row['coupling']*row['p']),
        e*row['coupling']*row['p'],s*s*row['ground_energy_form_upper'])
    return {**row,'distance_squared_upper':bound,'scaled_time':t,'scaled_energy':e}


def field_energy(n, kind, exchange_scaled_cost, sampled_square_integral,
                 coupling=F(1)):
    """Cost is delta^2/(Jp), including the actual selection Hamiltonian."""
    row = family(n,coupling)
    if kind not in ('Q','P'):
        raise ValueError('Q or P field required')
    exchange = rational(exchange_scaled_cost,'exchange scaled cost')
    l2 = rational(sampled_square_integral,'p^3 sum of test squares')
    factor = row['field']/row['coupling']
    if kind == 'P':
        factor /= row['p']**2
    selection = factor*l2/4
    return {**row,'kind':kind,'exchange_scaled_cost':exchange,
            'sampled_square_integral':l2,'selection_scaled_cost':selection,
            'total_scaled_cost':exchange+selection}


def history_energy(costs):
    value = F()
    for cost in costs:
        value = energy_after(value,rational(cost,'scaled kick'))
    return value


def cube_counts(radius):
    if type(radius) is not int or radius < 1:
        raise ValueError('positive integer source radius required')
    side = 2*radius+1
    return {'radius':radius,'sites':side**3,'edges':3*(side-1)*side**2,
            'odd_sites':(side**3-1)//2 if radius%2 == 0 else (side**3+1)//2}


def ground_preparation_counts(sites, edges, odd_sites, coupling, field, accuracy):
    """The literal field_gap ramp budget, evaluated without allocating a graph."""
    if (type(sites) is not int or sites<1 or type(edges) is not int or edges<0
        or type(odd_sites) is not int or not 0<=odd_sites<=sites):
        raise ValueError('valid finite graph counts required')
    j = rational(coupling,'coupling',True); b = rational(field,'field',True)
    eps = rational(accuracy,'accuracy',True)
    if eps>=1:
        raise ValueError('accuracy must be less than one')
    d,g = j*edges,b/8
    coefficient = 2*d/g**2+7*d*d/g**3
    time = 3*coefficient/eps if d else F()
    norm_sum = d+b*sites/2
    numerator = (d*time+norm_sum**2*time**2)/2
    steps = ceil(3*numerator/eps) if d else 0
    pulses = 16*sites+2*odd_sites+steps*(62*edges+11*sites)
    bits = max(1,(87*pulses*eps.denominator//eps.numerator).bit_length())
    while 29*pulses*F(1,1<<bits)>eps/3:
        bits += 1
    native = 29*pulses*F(1,1<<bits)
    adiabatic = coefficient/time if time else F()
    schedule = numerator/steps if steps else F()
    total = adiabatic+schedule+native
    if total>eps:
        raise ArithmeticError('native preparation bound exceeds target')
    return {'sites':sites,'edges':edges,'odd_sites':odd_sites,'gap_lower':g,
            'source_clock_duration':time,'slices':steps,'native_pulses_upper':pulses,
            'dyadic_port_bits':bits,'adiabatic_half_trace_upper':adiabatic,
            'schedule_operator_error_upper':schedule,
            'native_complete_isometry_error_upper':native,
            'total_half_trace_upper':total,'requested_accuracy':eps,
            'ground_postselection':False,'all_native_failure_records_retained':True,
            'large_circuit_executed':False}


def local_preparation(n, observation_radius, observed_sites, accuracy, coupling=F(1)):
    row = family(n,coupling)
    eps = rational(accuracy,'accuracy',True)
    if eps>=1:
        raise ValueError('accuracy must be less than one')
    selected = locality.choose_box(row['coupling'],row['field'],observation_radius,
                                   observed_sites,eps/2)
    # Even radius aligns the literal graph's root bipartition with epsilon_x.
    radius = selected['source_radius']
    if radius%2:
        radius *= 2
        selected = locality.box(row['coupling'],row['field'],observation_radius,
                                radius,observed_sites)
        while selected['half_trace_upper']>eps/2:
            radius *= 2
            selected = locality.box(row['coupling'],row['field'],observation_radius,
                                    radius,observed_sites)
    counts = cube_counts(radius)
    prep = ground_preparation_counts(counts['sites'],counts['edges'],counts['odd_sites'],
                                      row['coupling'],row['field'],eps/2)
    total = selected['half_trace_upper']+prep['total_half_trace_upper']
    if total>eps:
        raise ArithmeticError('prepared source marginal exceeds target')
    return {**row,'observation_radius':observation_radius,'observed_sites':observed_sites,
            'finite_source':counts,'marginal':selected,'preparation':prep,
            'prepared_marginal_half_trace_upper':total,'requested_accuracy':eps,
            'finite_source_matrix_allocated':False,
            'target':'the unique infinite-volume H_B ground at this positive B',
            'ordered_zero_field_state_accuracy_claimed':False,
            'initially_independent_reference_allowed':True}
