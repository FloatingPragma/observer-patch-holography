"""Volume-uniform local collision calibration of the actual cubic CAR action."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import factorial

from model.field_symmetry.symmetry import link_sign
from model.source_family_compatibility.local_matter import boundary_upper
from model.source_strong_coupling.readout import Interval


def neighbors(x):
    for axis in range(3):
        for step in (-1, 1):
            y = list(x); y[axis] += step
            yield tuple(y)


def sign(x, y):
    changed = [j for j in range(3) if x[j] != y[j]]
    if len(changed) != 1 or abs(x[changed[0]]-y[changed[0]]) != 1:
        raise ValueError('one actual cubic edge required')
    j = changed[0]
    return link_sign(x if x[j] < y[j] else y, j)


def paths(origin, axis):
    if (type(origin) is not tuple or len(origin) != 3 or
            any(type(k) is not int or k < 2 for k in origin) or
            type(axis) is not int or axis not in range(3)):
        raise ValueError('interior integer source coordinate and axis required')
    endpoint = list(origin); endpoint[axis] += 1; endpoint = tuple(endpoint)
    collar = {origin, endpoint, *neighbors(origin), *neighbors(endpoint)}
    rows = []
    for z in neighbors(origin):
        for w in neighbors(z):
            if endpoint in neighbors(w):
                path = (origin, z, w, endpoint)
                if not set(path) <= collar:
                    raise ArithmeticError('three-step source return escaped the one-collar')
                product_sign = sign(origin,z)*sign(z,w)*sign(w,endpoint)
                rows.append((path, product_sign))
    normalized = [v*sign(origin,endpoint) for _,v in rows]
    if (len(rows), normalized.count(1), normalized.count(-1), sum(normalized)) != (15,11,4,7):
        raise ArithmeticError('actual pi-flux path identity failed')
    # h/c has edge coefficient1/2; (h/c)^2_xx=3/2.
    a = F(1,2)
    b = 2*(F(7,16)+3*F(3,2)*F(1,4))
    return {'origin': origin, 'axis': axis, 'collar_sites': len(collar),
            'length_three_paths': tuple(rows), 'signed_sum_in_edge_gauge': 7,
            'positive_paths': 11, 'negative_paths': 4,
            'quadratic_coefficient_over_c_squared': a,
            'hopping_quartic_inner_product_over_c_fourth': b,
            'contact_subtraction_coefficient': 8*b}


def local_path_certificate():
    rows = [paths(tuple(2+k for k in parity),axis)
            for parity in product((0,1),repeat=3) for axis in range(3)]
    return {'parity_and_direction_cases': len(rows), 'cases': rows,
            'all_cases_contact_subtraction_coefficient': F(25),
            'one_graph_step_collar_suffices_for_coefficients': True,
            'all_large_box_statement_uses_local_path_exhaustion': True}


def collision_effect():
    """Actual complete raw8 effect on all sixteen local occupation states."""
    import sympy as s
    from model.defect_field_bridge.full_domain import finite_branches
    from model.field_symmetry import symmetry
    from model.full_fock_field.field import FieldGeometry
    expected = [int((mask & 3) == 3 or (mask >> 2) == 3) for mask in range(16)]
    for law in ('C', 'E'):
        rows = finite_branches(sites=2, flavors=2, law=law, seam=(0, 1))
        total, collision = [0]*16, [0]*16
        for (coin, raw), columns in rows.items():
            for origin, output, phase in columns:
                total[origin] += 1
                if raw == 8:
                    collision[origin] += 1
                    if output != origin or phase != 1:
                        raise ArithmeticError('actual collision branch is not a hold')
        if total != [2]*16 or collision != [2*x for x in expected]:
            raise ArithmeticError('complete source collision effect changed')
    effect = s.diag(*expected)
    for axis in range(3):
        sides = tuple(2 if j == axis else 1 for j in range(3))
        source = symmetry.operators(FieldGeometry(sides))
        frame = s.Matrix(source.frame.real.toarray())+s.I*s.Matrix(source.frame.imag.toarray())
        if frame*effect*frame.H != effect:
            raise ArithmeticError('collision effect is not invariant under the actual source frame')
    return {'full_local_Fock_dimension': 16, 'raw_record': 8,
            'effect': 'D_x+D_y-D_x D_y', 'two_particle_restriction': 'D_x+D_y',
            'uniform_law_coin_raw_records_per_seam': 64,
            'all_padded_records_retained': True,
            'all_three_actual_frame_effect_checks': True,
            'full_fine_instrument_covariance_asserted': False}


def _positive(value):
    if type(value) not in (int,F) or value <= 0:
        raise ValueError('positive exact rational required')
    return F(value.numerator,value.denominator)


def local_probability_error(radius, time, coupling_upper, marginal_error=F()):
    """Collision event, arbitrary exterior; not a source/reference channel norm."""
    t,c = _positive(time),_positive(coupling_upper)
    if type(marginal_error) not in (int,F) or not 0 <= marginal_error <= 1:
        raise ValueError('local marginal half-trace error in[0,1] required')
    # The tested edge is in the radius-one support cube, support size two.
    return min(F(1),F(marginal_error)+boundary_upper(radius,1,2,t,c))


def identify(time, p1, p2, coupling_upper, contact_upper):
    """Complete local-event probability enclosures with a volume-free tail."""
    t,c,u = map(_positive,(time,coupling_upper,contact_upper))
    if any(not isinstance(p,Interval) or p.low < 0 or p.high > 1 for p in (p1,p2)):
        raise ValueError('complete collision probabilities in[0,1] required')
    norm=6*c+u
    error_c=F(32,27)*norm**6*t**4
    error_q=F(544,45)*norm**6*t**2
    first=Interval((8*p1.low-p2.high/2)/(3*t**2),
                   (8*p1.high-p2.low/2)/(3*t**2)).expanded(error_c)
    second=Interval((8*p1.low-2*p2.high)/t**4,
                    (8*p1.high-2*p2.low)/t**4).expanded(error_q)
    if first.low <= 0 or second.low < 0:
        raise ValueError('finite data do not certify positive source coefficients')
    lo=second.low/first.high**2-25
    hi=second.high/first.low**2-25
    if hi < 0:
        raise ValueError('data exclude the specified cubic action/preparation')
    return {'coupling_squared': first,
            'coupling_squared_times_U_squared_plus_25c_squared': second,
            'contact_over_coupling_squared': Interval(max(F(),lo),hi),
            'finite_volume_independent_C_error':error_c,
            'finite_volume_independent_Q_error':error_q}


def collar_plan(time, coupling_upper, probability_error):
    t,c,epsilon=map(_positive,(time,coupling_upper,probability_error))
    if epsilon > 1:
        raise ValueError('probability error must not exceed one')
    radius=2
    while local_probability_error(radius,2*t,c)>epsilon:
        radius+=1
    return {'collar_radius': radius, 'collar_sites':(2*radius+1)**3,
            'blank_local_modes_for_supplied_reset':2*(2*radius+1)**3,
            'both_time_settings_local_probability_error':local_probability_error(radius,2*t,c),
            'arbitrary_exterior_state_allowed':True,
            'global_vacuum_preparation_required':False,
            'pure_two_CAR_collar_input_required':True,
            'reset_pulse_tape_executed':False}


def uniform_attempt_budget(edge_count, collision_error, failure_probability):
    """Full-attempt sampling, including all seams without postselection."""
    from model.source_family_compatibility.source_action_identification import sampling_budget
    if type(edge_count) is not int or edge_count < 1:
        raise ValueError('positive exact source edge count required')
    epsilon, delta = map(_positive, (collision_error, failure_probability))
    if epsilon >= 1 or delta >= 1:
        raise ValueError('probability and confidence errors must be below one')
    budget = sampling_budget(epsilon/edge_count, delta)
    return {**budget, 'source_edges': edge_count,
            'actual_indicator': 'chosen seam and raw8',
            'indicator_mean': 'p(t)/source_edges',
            'score_multiplier': edge_count,
            'collision_probability_error': epsilon,
            'uniform_seam_selector_required': True,
            'other_seam_records_discarded_or_postselected': False,
            'repair_attempts_between_preparation_and_readout': 0}


def pair_hop(state, sites, coupling, contact):
    """Exact invariant one-up/one-down sector in the spin-first CAR basis."""
    c,u=map(_positive,(coupling,contact)); sites=frozenset(sites)
    out=defaultdict(F)
    for (x,y),value in state.items():
        if x not in sites or y not in sites:
            raise ValueError('input outside the finite source collar')
        if x==y:out[x,y]+=u*value
        for z in neighbors(x):
            if z in sites:out[z,y]+=c*sign(x,z)*value/2
        for z in neighbors(y):
            if z in sites:out[x,z]+=c*sign(y,z)*value/2
    return {key:value for key,value in out.items() if value}


def probability_enclosure(radius, coupling, contact, time, order=8):
    """Sparse rational source Taylor control; no full-Fock allocation."""
    if type(radius) is not int or radius<2 or type(order) is not int or not 1<=order<=12:
        raise ValueError('finite radius>=2 and Taylor order1..12 required')
    c,u,t=map(_positive,(coupling,contact,time))
    center=(radius+2,)*3; x=center; y=(center[0]+1,center[1],center[2])
    sites=set(product(*(range(k-radius,k+radius+1) for k in center)))
    power={(x,y):F(1)}; real=defaultdict(F); imag=defaultdict(F)
    for n in range(order+1):
        part=real if n%2==0 else imag
        factor=(-1)**((n+1)//2)*t**n/factorial(n)
        for key,value in power.items():part[key]+=factor*value
        if n<order:power=pair_hop(power,sites,c,u)
    value=sum(real[z,z]**2+imag[z,z]**2 for z in (x,y))
    remainder=((6*c+u)*t)**(order+1)/factorial(order+1)
    error=2*remainder+remainder**2
    return Interval(max(F(),value-error),min(F(1),value+error))


def certificate():
    witnesses = []
    for exponent, lo, hi in ((10, F(62), F(66)),
                              (12, F(639, 10), F(641, 10)),
                              (14, F(6399, 100), F(6401, 100))):
        t = F(1, 2**exponent)
        first = probability_enclosure(2, F(1, 8), F(1), t, 6)
        second = probability_enclosure(2, F(1, 8), F(1), 2*t, 6)
        derived = identify(t, first, second, F(1, 8), F(1))
        interval = derived['contact_over_coupling_squared']
        if not lo <= interval.low <= F(64) <= interval.high <= hi:
            raise ArithmeticError('source ratio does not satisfy the admitted interval')
        witnesses.append({'time': t, 'ratio_squared_enclosure': (lo, hi),
                          'actual_sparse_rational_source_admission': True,
                          'finite_source_radius': 2, 'finite_source_sites': 125,
                          'is_probability_of_the_larger_collar_plan': False})
    t = F(1, 2**12)
    return {'actual_cubic_path_exhaustion': local_path_certificate(),
            'complete_collision_detector': collision_effect(),
            'configured_action': {'c': F(1, 8), 'U': F(1)},
            'finite_time_ratio_witnesses': witnesses,
            'separate_arbitrary_exterior_collar_plan': collar_plan(t, F(1, 8), t**6),
            'uniform_global_attempt_example': uniform_attempt_budget(300, F(1, 100), F(1, 100)),
            'known_local_calibration_preparation_consumed': True,
            'complete_history_or_reference_error_from_marginal_inferred': False,
            'interacting_ground_preparation_inferred': False,
            'capacity_selects_action_or_clock': False,
            'large_native_experiment_executed': False,
            'physical_particles': [], 'physical_masses': []}
