"""Exact finite source energy moments and local clock comparison bounds.

Research attachment to the existing fair binary source. The operator
comparison concerns controlled unitary continuations, not unread repair.
"""
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from math import comb, isqrt

import numpy as np

from model.observer import native


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('an exact integer in range is required')
    return x


def graph_data(ports, edges):
    integer(ports, 4)
    if (type(edges) is not tuple or edges != tuple(sorted(set(edges)))
            or any(type(e) is not tuple or len(e) != 2
                   or any(type(v) is not int for v in e)
                   or not 0 <= e[0] < e[1] < ports for e in edges)):
        raise ValueError('canonical simple graph edges are required')
    degree = Counter()
    for a,b in edges:
        degree[a] += 1; degree[b] += 1
    if not edges:
        raise ValueError('a nonempty graph is required')
    return len(edges), max(degree.values()), sum(d*(d-1) for d in degree.values())


def moments(ports, seams, adjacent_ordered, total):
    """Normalized trace moments on the complete fixed-number sector.

    Only graph counts are used. No global carrier is allocated.
    """
    integer(ports, 4); integer(seams, 1); integer(adjacent_ordered); integer(total)
    if total > ports or adjacent_ordered > seams*(seams-1):
        raise ValueError('occupation or edge-pair count outside range')
    q, v, m = total, ports, seams
    u = F(q*(v-q),v*(v-1))
    w = F(q*(q-1)*(v-q)*(v-q-1),v*(v-1)*(v-2)*(v-3))
    disjoint_ordered = m*(m-1)-adjacent_ordered
    first = 2*m*u
    second = 4*m*u+adjacent_ordered*u+4*disjoint_ordered*w
    variance = second-first*first
    c = 1-u
    mu_rho = first-variance/(2*m*c)
    if variance < 0 or mu_rho < 0:
        raise ValueError('graph counts do not give valid source moments')
    return {'trace_G': first, 'trace_G_squared': second,
            'canonical_variance': variance, 'trace_P': c,
            'rho_mean_G': mu_rho, 'rho_variance_upper': variance/c,
            'local_clock_multiplier': 1/(1-mu_rho/(2*m)),
            'density_limit_multiplier': 1/(1-F(q,v)*(1-F(q,v)))}


def sqrt_upper(value, bits=40):
    if type(value) not in (int,F) or value < 0:
        raise ValueError('a nonnegative exact rational is required')
    integer(bits, 1)
    if bits > 128:
        raise ValueError('at most 128 square-root bits are supported')
    value = F(value); scale = 1 << bits
    n = isqrt((value.numerator*scale*scale)//value.denominator)
    if F(n*n,scale*scale) < value:
        n += 1
    result = F(n,scale)
    if result*result < value:
        raise ArithmeticError('square-root enclosure failed')
    return result


def local_clock_bound(ports,seams,adjacent_ordered,total,degree,support,time,
                      floor=F(10,21)):
    """Rational upper bound in ||.||_(2,rho) for norm-one local A.

    The supplied floor is an explicit premise on P, authenticated by the
    owning federation theorem for the default 10/21 value.
    """
    integer(degree, 1); integer(support, 1)
    if (support > ports or type(time) not in (int,F)
            or type(floor) not in (int,F) or not 0 < floor <= 1):
        raise ValueError('finite local support, exact time and positive floor required')
    row = moments(ports,seams,adjacent_ordered,total)
    b1 = 4*support*degree
    b2 = 16*support*(support+1)*degree*degree
    sigma = sqrt_upper(row['rho_variance_upper'])
    rate = (b1*sigma+F(b2,2))/(2*seams*floor*floor)
    difference = abs(row['local_clock_multiplier']-row['density_limit_multiplier'])
    return dict(row, local_support=support, maximum_degree=degree,
        commutator_norm_upper=b1, double_commutator_norm_upper=b2,
        rho_standard_deviation_upper=sigma, generator_comparison_upper=rate,
        controlled_local_dynamics_upper=min(F(2),abs(time)*rate),
        density_clock_dynamics_upper=min(F(2),abs(time)*(rate+difference*b1)),
        dominates_reference_expectation_factor_squared=1/floor)


def matrices(ports,edges,total):
    """Bounded exact integer source matrices, independently of sparse transfer."""
    integer(ports, 4); integer(total)
    if ports>24 or total>ports or comb(ports,total)>256:
        raise ValueError('at most 24 ports and 256 complete sector states are supported')
    graph_data(ports,edges)
    states=tuple(tuple(int(i in occupied) for i in range(ports))
                 for occupied in combinations(range(ports),total))
    if not states or len(states)>256:
        raise ValueError('at most 256 nonempty sector basis states are supported')
    index={x:i for i,x in enumerate(states)}; d=len(states)
    identity=np.eye(d,dtype=np.int64); local=[]
    for a,b in edges:
        swap=np.zeros((d,d),dtype=np.int64)
        for j,x in enumerate(states):
            y=list(x); y[a],y[b]=y[b],y[a]
            swap[index[tuple(y)],j]=1
        local.append(identity-swap)
    return states,tuple(local),sum(local,np.zeros((d,d),dtype=np.int64))


def source_numerator(matrix,states,edges):
    """4M times unread fair C/E Schrödinger map, retaining raw partitions.

    This is an exact contraction of the complete instrument. It does not
    claim an expanded execution or omit zero-probability raw slots.
    """
    d=len(states); index={x:i for i,x in enumerate(states)}
    if matrix.shape!=(d,d) or matrix.dtype!=np.dtype('int64'):
        raise ValueError('bounded exact int64 source matrix is required')
    if d*max(1,int(abs(matrix).max()))*4*len(edges)>=2**62:
        raise OverflowError('source contraction exceeds the exact integer budget')
    out=np.zeros_like(matrix)
    for edge in edges:
        for coin in (0,1):
            target=[]; current=[]; parity=[]
            for state in states:
                image,(_,j)=native.native(state,edge,coin)
                target.append(index[image]); current.append(j)
                parity.append((state[edge[0]]+state[edge[1]])%2)
            for law in ('C','E'):
                for i in range(d):
                    for j in range(d):
                        if current[i]==current[j] and (law=='C' or parity[i]==parity[j]):
                            out[target[i],target[j]]+=matrix[i,j]
    return out


def exact_energy_certificate(ports,edges,total):
    m,degree,adjacent=graph_data(ports,edges)
    states,local,g=matrices(ports,edges,total); d=len(states)
    identity=np.eye(d,dtype=np.int64)
    if not np.array_equal(source_numerator(identity,states,edges),4*m*identity):
        raise ArithmeticError('complete source channel is not unital')
    phi_numerator=source_numerator(g,states,edges)
    dissipation=4*m*g-phi_numerator
    trace_p=F(d)-F(int(np.trace(g)),2*m)
    heating=F(int(np.trace(g@dissipation)),8*m*m*trace_p)
    expected=moments(ports,m,adjacent,total)
    if F(int(np.trace(g)),d)!=expected['trace_G'] or F(int(np.trace(g@g)),d)!=expected['trace_G_squared']:
        raise ArithmeticError('independent complete sector moments disagree')
    currents=[]
    for e,a in enumerate(local):
        divergence=np.zeros_like(g)
        for f,b in enumerate(local):
            # The actual Hermitian current is i times this real antisymmetric matrix.
            j=a@b-b@a
            if not set(edges[e])&set(edges[f]) and np.any(j):
                raise ArithmeticError('energy current exceeds its native seam support')
            if not np.array_equal(j,-j.T):
                raise ArithmeticError('energy current is not Hermitian after multiplication by i')
            divergence+=j
            if e<f and np.any(j): currents.append((e,f))
        if np.any(g@a-a@g+divergence):
            raise ArithmeticError('exact local energy continuity failed')
    # The positive transfer reference exists only when the independently
    # established transfer floor applies. This exact coefficient is retained
    # even for controls lacking such a floor; no log is evaluated there.
    return {'ports':ports,'seams':m,'sector':total,'dimension':d,
            'adjacent_ordered_edge_pairs':adjacent,
            'source_moments':expected,'trace_P_unnormalized':trace_p,
            'repair_exchange_energy_change_on_P_reference':heating,
            'nonzero_unordered_energy_currents':len(currents),
            'all_seam_continuity_identities_exact':True,
            'repair_energy_conservation':not np.any(dissipation)}


def h66_detector_heating_bound(label_error=F(1,10**30), preparation_error=F(0)):
    """Conditional complete two-setting energy-score resolution.

    Ideal rho=P/56 is the input. Nonzero preparation error must be a
    separately supplied complete half-trace bound, not a numerical tolerance.
    The existing 66-word measurement retains every failure with score zero.
    """
    for value in (label_error,preparation_error):
        if type(value) not in (int,F) or value<0:
            raise ValueError('exact nonnegative error bounds required')
    maximum=F(30); instrument=F(1,12000); spectral=F(27,10**8)
    deduction=4*maximum*(instrument+preparation_error)+2*(spectral+label_error)
    return {'exact_transfer_energy_increase_lower':F(37,840),
            'complete_measurement_score_upper':maximum,
            'two_setting_systematic_error_upper':deduction,
            'complete_measured_score_increase_lower':F(37,840)-deduction,
            'ideal_prepared_reference_assumed':preparation_error==0}


def complete_rational_instrument(expected_labels, branches):
    """Check a complete finite real-rational Kraus menu, including zero slots.

    This bounded matrix checker authenticates the supplied finite instrument.
    Its locality on the declared source ports is a separate attachment.
    """
    if (type(expected_labels) is not tuple or not expected_labels
            or len(expected_labels)>64
            or any(type(s) is not str or not s for s in expected_labels)
            or len(set(expected_labels))!=len(expected_labels)
            or type(branches) is not tuple
            or tuple(x[0] for x in branches if type(x) is tuple and len(x)==2)!=expected_labels):
        raise ValueError('every declared distinct outcome, including zero slots, is required')
    matrices=[];shape=None
    for _,value in branches:
        if (type(value) is not tuple or not value or len(value)>64
                or any(type(row) is not tuple or len(row)!=len(value) for row in value)
                or any(type(x) not in (int,F) for row in value for x in row)):
            raise ValueError('bounded square exact rational Kraus matrices are required')
        if shape is None:shape=len(value)
        if shape!=len(value):raise ValueError('all branches must share one carrier')
        matrices.append(np.array(value,dtype=object))
    total=sum((k.T@k for k in matrices),np.zeros((shape,shape),dtype=object))
    if not np.array_equal(total,np.eye(shape,dtype=object)):
        raise ValueError('the complete Kraus menu is not trace preserving')
    return tuple(matrices)


def finite_history_bound(ports,seams,adjacent_ordered,total,degree,supports,
                         outcome_labels,times,floor=F(10,21)):
    """Complete fixed classical history TV bound, conditional on local menus.

    One exact local Kraus instrument follows each controlled time interval.
    Every listed outcome is retained, including zero outcomes and failures.
    The fixed menu must be attached to a complete instrument independently.
    """
    integer(degree,1)
    if (type(supports) is not tuple or not 1<=len(supports)<=16
            or type(outcome_labels) is not tuple or type(times) is not tuple
            or len(supports)!=len(outcome_labels) or len(times)!=len(supports)
            or any(type(s) is not int or not 1<=s<=ports for s in supports)
            or any(type(t) not in (int,F) for t in times)
            or type(floor) not in (int,F) or not 0<floor<=1):
        raise ValueError('fixed local supports, full menus, exact times and transfer floor required')
    count=1
    for labels in outcome_labels:
        if (type(labels) is not tuple or not 1<=len(labels)<=64
                or any(type(s) is not str or not s for s in labels)
                or len(set(labels))!=len(labels)):
            raise ValueError('distinct complete declared outcome labels required')
        count*=len(labels)
    if count>2**32:raise ValueError('declared full outcome menu exceeds its arithmetic budget')
    row=moments(ports,seams,adjacent_ordered,total)
    sigma=sqrt_upper(row['rho_variance_upper'])
    dominance=sqrt_upper(1/F(floor))
    b=tuple(4*s*degree for s in supports)
    c=tuple(16*s*(s+1)*degree*degree for s in supports)
    terms=[]
    for j,t in enumerate(times):
        B=sum(b[j:]);C=sum(c[j:])
        rate=(2*B*sigma+C+2*B*B)/(2*seams*floor*floor)
        terms.append(abs(t)*dominance*rate)
    per_history=sum(terms,F())
    return {'full_declared_history_count':count,
            'complete_outcome_labels':outcome_labels,
            'per_history_probability_error_upper':min(F(1),per_history),
            'complete_history_total_variation_upper':min(F(1),count*per_history/2),
            'interval_probability_error_uppers':tuple(terms),
            'clock_multiplier':row['local_clock_multiplier'],
            'initial_reference':'rho=P/Tr(P)',
            'scope':'fixed finite classical output menu; terminal source/reference states are not compared'}
