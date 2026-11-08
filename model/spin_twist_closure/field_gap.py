"""Fixed-selection-field source gap and complete adiabatic preparation budgets."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.certified_reference import certificate as cert
from model.native_ordered_ground import source
from model.native_spin_observer import circuits, precision


GAP_REFERENCE = 'https://arxiv.org/abs/2607.10765v1'
ADIABATIC_REFERENCE = 'https://arxiv.org/abs/quant-ph/0603175v3'


def gap(field):
    field = F(field)
    if field <= 0:
        raise ValueError('a strictly positive selection field is required')
    return field/8


def _odd_sites(sites, edges):
    colors = {}
    adjacency = [[] for _ in range(sites)]
    for x,y in edges:
        adjacency[x].append(y); adjacency[y].append(x)
    for root in range(sites):
        if root in colors:
            continue
        colors[root] = 0
        queue = [root]
        while queue:
            x = queue.pop()
            for y in adjacency[x]:
                if y not in colors:
                    colors[y] = 1-colors[x]; queue.append(y)
                elif colors[y] == colors[x]:
                    raise ValueError('the actual source requires a bipartite graph')
    return tuple(x for x in range(sites) if colors[x])


@lru_cache(maxsize=1)
def source_certificate():
    rows = []
    for sites,edges in ((1,()),(2,((0,1),)),(3,((0,1),(1,2))),
                        (4,((0,1),(1,2),(2,3),(0,3))),
                        (8,source.cube()[0])):
        h,scale = source.matrix(sites,edges,F(1),F(1,4))
        dimension = 1 << sites
        walsh = np.array([[(-1)**((i&j).bit_count()) for j in range(dimension)]
                          for i in range(dimension)],dtype=np.int64)
        target = np.zeros_like(h)
        coupling,field = scale//2,scale//8
        for mask in range(dimension):
            target[mask,mask] -= field*(sites-2*mask.bit_count())
            for x,y in edges:
                if not (((mask>>x)^(mask>>y))&1):
                    target[mask,mask] -= coupling
                    target[mask^(1<<x)^(1<<y),mask] += coupling
        if not np.array_equal(walsh@h@walsh,dimension*target):
            raise ArithmeticError('actual source does not equal the Lee-Yang substitution')
        rows.append({'sites':sites,'edges':len(edges),'all_matrix_rows':dimension,
                     'integer_Hadamard_conjugation_exact':True})
    return {'primary_gap_theorem':GAP_REFERENCE,
            'theorem':'Definition 1 and Theorem 11, Eq.41',
            'gap_theorem_is_an_explicit_external_mathematical_input':True,
            'source_couplings_after_global_Hadamard':
                {'w_z':'sJ/4','w_x':'sJ/4','w_y':'-sJ/4','w_xy':'0','w_yx':'0','h':'B/2'},
            'Lee_Yang_inequality':'sJ/4=(abs(sJ/2)+abs(0))/2',
            'constant_removed':'-sJ*edge_count/4',
            'uniform_ramp_parameter_interval':['0','1'],
            'full_finite_ground_unique':True,'full_finite_gap_lower':'B/8',
            'graph_connectivity_required':False,
            'finite_source_matrix_controls':rows,
            'zero_field_gap_or_ground_marginal_modulus_claimed':False}


def resource_plan(sites, edges, coupling=F(1), field=F(1,4), error=F(1,100)):
    source.parameters(sites,edges,coupling,field)
    error = F(error)
    if not 0 < error < 1:
        raise ValueError('an error strictly between zero and one is required')
    odd = _odd_sites(sites,edges)
    derivative = coupling*len(edges)
    g = gap(field)
    adiabatic_coefficient = 2*derivative/g**2+7*derivative**2/g**3
    time = 3*adiabatic_coefficient/error if derivative else F()
    norm_sum = derivative+field*sites/2
    # Ramp freezing plus first-order unitary product error over all slices.
    numerator = (derivative*time+norm_sum**2*time**2)/2
    steps = math.ceil(3*numerator/error) if derivative else 0
    schedule_error = numerator/steps if steps else F()
    pulses = 16*sites+2*len(odd)+steps*(62*len(edges)+11*sites)
    bits = 1
    while 29*pulses*F(1,1<<bits) > error/3:
        bits += 1
    native_error = 29*pulses*F(1,1<<bits)
    adiabatic_error = adiabatic_coefficient/time if time else F()
    total = adiabatic_error+schedule_error+native_error
    if total > error:
        raise ArithmeticError('complete finite ground-preparation budget exceeded')
    return {'sites':sites,'edges':len(edges),'coupling':coupling,'selection_field':field,
            'physical_odd_sites':list(odd),'gap_lower':g,
            'primary_adiabatic_theorem':ADIABATIC_REFERENCE,
            'adiabatic_theorem':'Theorem 3, Eq.6, rank-one band',
            'derivative_norm_upper':derivative,'second_derivative_exact':F(),
            'adiabatic_coefficient':adiabatic_coefficient,'source_clock_duration':time,
            'piecewise_constant_slices':steps,'native_pulses_upper':pulses,
            'dyadic_port_bits':bits,'adiabatic_half_trace_error_upper':adiabatic_error,
            'ramp_and_product_operator_error_upper':schedule_error,
            'complete_native_isometry_error_upper':native_error,
            'total_prepared_output_half_trace_error_upper':total,
            'requested_error':error,'all_first_failure_records_retained':True,
            'ground_postselection_used':False,'large_ramp_executed':False,
            'finite_volume_target':'the actual physical H_B ground at the specified B',
            'infinite_volume_or_zero_field_accuracy_included':False,
            'admitted_coherent_controls_and_blank_registers_required':True,
            'refined_rational_port_interface_required':True,
            'fixed_binary64_error_sent_to_zero':False}


def preparation_word(sites, edges, coupling, field, duration, steps):
    """Lazy literal word; budgets do not materialize the large certified ramp."""
    source.parameters(sites,edges,coupling,field)
    duration = F(duration)
    if type(steps) is not int or steps < 0 or duration < 0 or bool(duration) != bool(steps):
        raise ValueError('positive duration and slices, or the zero-duration product, required')
    odd = set(_odd_sites(sites,edges))
    yield from circuits.encode(sites)
    for site in range(sites):
        yield from circuits.spin_rotation('Y',cert._pi(-1,4),site,sites)
        if site in odd:
            yield ('I',(2*site,),cert._pi(-1,2))
            yield ('Z',(2*site,),cert._pi(1,2))
    if not steps:
        return
    delta = duration/steps
    for index in range(steps):
        strength = coupling*F(index,steps)
        if strength:
            for edge in edges:
                yield from circuits.exchange(delta,strength,sites,edge)
        for site in range(sites):
            sign = -1 if site in odd else 1
            yield from circuits.spin_rotation('X',delta*field*sign/2,site,sites)


def dyadic_preparation_word(sites, edges, coupling, field, duration, steps, bits):
    for kind,targets,angle in preparation_word(sites,edges,coupling,field,duration,steps):
        yield {'targets':list(targets),'port_command':precision.dyadic_command(kind,angle,bits)}


@lru_cache(maxsize=1)
def small_native_execution():
    """Execute a short source ramp; it audits compilation, not adiabatic accuracy."""
    from collections import Counter
    from model.centered_native_preparation import tree
    from model.native_spin_observer import algebra, execute
    from scipy.linalg import expm

    sites,edges,j,b,time,steps = 2,((0,1),),F(1),F(1,4),F(1,16),2
    template = list(preparation_word(sites,edges,j,b,time,steps))
    word = list(circuits.emit(template))
    consumed,angle_error,variables = cert.consume(word,0,template)
    if consumed != len(word) or variables:
        raise ArithmeticError('the source ramp has an unresolved angle')
    port_error = tree.port_bound(Counter((event['kind'],event['angle']) for event in word))
    error = angle_error+port_error
    vector = np.array([1,1,-1,-1],complex)/2  # +X on site0, -X on site1
    bond = (sum(algebra.tensor(axis*2) for axis in 'XYZ')-np.eye(4))/4
    delta = float(time/steps)
    for index in range(steps):
        vector = expm(-1j*delta*float(j)*index/steps*bond)@vector
        for site,sign in ((0,1),(1,-1)):
            label = 'XI' if site == 0 else 'IX'
            vector = expm(1j*delta*float(b)*sign/2*algebra.tensor(label))@vector
    target = algebra.inclusion(sites)@vector[:,None]
    initial = np.zeros((16,1),complex); initial[0] = 1
    actual = execute.complete(word,4,initial)
    distance = float(np.sqrt(np.linalg.norm(actual.success-target)**2+actual.failure_probability))
    if distance > float(error)+execute.ADMISSION:
        raise ArithmeticError('complete native ramp disagrees with its literal source schedule')
    if max(distance,actual.failure_probability,abs(actual.normalization_residual)) >= 1e-10:
        raise ArithmeticError('short native ramp exceeds its diagnostic enclosure')
    return {'sites':sites,'slices':steps,'duration':str(time),'pulses':len(word),
            'word_sha256':cert.word_digest(word),'all_work_rows':16,
            'angle_error_bound':str(angle_error),'port_error_bound':str(port_error),
            'complete_isometry_error_bound':str(error),
            'complete_distance_diagnostic_interval':['0','1/10000000000'],
            'complete_first_failure_probability_interval':['0','1/10000000000'],
            'normalization_residual_absolute_upper':'1/10000000000',
            'failure_slots':actual.native_failure_slots,
            'short_run_is_a_certified_ground_preparation':False,
            'long_certified_adiabatic_budget_executed':False,
            'all_failure_records_retained':True}
