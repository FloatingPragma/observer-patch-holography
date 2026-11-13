"""Outward integer arithmetic for the exact native event process.

Uniformization is an identity for the original event channel. Its virtual
P steps are a calculation device; they do not replace the native event
clock or its full record instrument.
"""
from fractions import Fraction as F
from hashlib import sha256
import json

from model.openwave_native_process import source as native
from model.openwave_native_response.source import coherence_multiplier

BITS=128
SCALE=1<<BITS
ORDER=128
TIMES=(F(1,4),F(1),F(4))
WIDTHS=(0,1,2)


def ceil_div(a,b):return -(-a//b)


def rational_interval(n,d,scale=SCALE):
    if n<0 or d<=0:raise ValueError('nonnegative rational required')
    return n*scale//d,ceil_div(n*scale,d)


def multiply(a,b,scale=SCALE):
    if min(*a,*b)<0:raise ValueError('positive intervals required')
    return a[0]*b[0]//scale,ceil_div(a[1]*b[1],scale)


def power_interval(n,d,exponent,scale=SCALE):
    if type(exponent) is not int or exponent<0:raise ValueError('nonnegative integer power required')
    base=rational_interval(n,d,scale);out=(scale,scale)
    while exponent:
        if exponent&1:out=multiply(out,base,scale)
        exponent//=2
        if exponent:base=multiply(base,base,scale)
    return out


def binomial_intervals(attempts,denominator,order=ORDER,scale=SCALE):
    """Binomial(attempts,6/denominator), including certified omitted tail."""
    if not (type(attempts) is int and attempts>=0 and denominator>6 and order>=0):
        raise ValueError('native finite attempt/degree domain required')
    weight=power_interval(denominator-6,denominator,attempts,scale)
    weights=[weight]
    for k in range(min(order,attempts)):
        n=(attempts-k)*6;d=(k+1)*(denominator-6)
        weight=(weight[0]*n//d,ceil_div(weight[1]*n,d))
        weights.append(weight)
    tail=max(0,scale-sum(w[0] for w in weights))
    assert sum(w[0] for w in weights)<=scale
    return weights,tail


def distances(g,origin=0):
    if type(origin) is not int or not 0<=origin<g.ports:raise ValueError('source port required')
    d=[-1]*g.ports;d[origin]=0;queue=[origin]
    for i in queue:
        for j in g.neighbours[i]:
            if d[j]<0:d[j]=d[i]+1;queue.append(j)
    assert min(d)>=0
    return d


def local_numerator_step(g,vector,allowed=None):
    """Six times P=I-L/6, or its exact unarrived block, without truncation."""
    if len(vector)!=g.ports:raise ValueError('full port vector required')
    output=[(6-len(g.neighbours[i]))*v for i,v in enumerate(vector)]
    for a,b in g.seams:
        output[a]+=vector[b];output[b]+=vector[a]
    if allowed is not None:
        output=[v if allowed[i] else 0 for i,v in enumerate(output)]
    assert min(output)>=0
    return output


def outward_pair(lower,upper):
    assert lower<=upper
    return [str(F(lower,SCALE)),str(F(upper,SCALE))]


def summarize_packet(level,width):
    if width not in WIDTHS:raise ValueError('frozen packet width required')
    g=native.graph(level);M=len(g.seams);d=distances(g)
    selected=[i for i,x in enumerate(d) if x<=width];size=len(selected)
    allowed=[x<4 for x in d]
    vector=[int(i in selected) for i in range(g.ports)]
    unarrived=list(vector);den=size
    schedules=[]
    for tau in TIMES:
        attempts=2*M*tau;assert attempts.denominator==1
        weights,tail=binomial_intervals(int(attempts),2*M)
        schedules.append({'tau':tau,'attempts':int(attempts),'weights':weights,'tail':tail,
                          'lower':[0]*g.ports,'upper':[0]*g.ports,
                          'nohit':[0,0]})
    # Every intermediate P^k vector is retained over the full finite graph.
    # No absorbing physical boundary is substituted: unarrived is a second
    # observer-flag block; the arrived block is full minus unarrived.
    for k in range(ORDER+1):
        assert sum(vector)==den and all(a<=b for a,b in zip(unarrived,vector,strict=True))
        nohit=sum(unarrived)
        for row in schedules:
            if k>=len(row['weights']):continue
            lo,hi=row['weights'][k]
            row['nohit'][0]+=lo*nohit//den
            row['nohit'][1]+=ceil_div(hi*nohit,den)
            for i,v in enumerate(vector):
                row['lower'][i]+=lo*v//den
                row['upper'][i]+=ceil_div(hi*v,den)
        if k<ORDER:
            vector=local_numerator_step(g,vector)
            unarrived=local_numerator_step(g,unarrived,allowed)
            den*=6
    results=[]
    for row in schedules:
        tail=row['tail'];low=row['lower'];high=[min(SCALE,v+tail) for v in row['upper']]
        assert sum(low)<=SCALE<=sum(high)
        nohitlo=row['nohit'][0];nohithi=min(SCALE,row['nohit'][1]+tail)
        arrival=(max(0,SCALE-nohithi),SCALE-nohitlo)
        radial=(sum(v*r for v,r in zip(low,d,strict=True)),sum(v*r for v,r in zip(high,d,strict=True)))
        target=(sum(low[i] for i in range(g.ports) if not allowed[i]),
                min(SCALE,sum(high[i] for i in range(g.ports) if not allowed[i])))
        support=(sum(low[i] for i in selected),min(SCALE,sum(high[i] for i in selected)))
        # The complete coherent packet has diagonal1/size and all selected
        # offdiagonals1/size. Native coherences never acquire another address.
        cohlo=cohhi=0
        for i in selected:
            for j in selected:
                if i==j:continue
                lam=coherence_multiplier(g,i,j)
                l,h=power_interval(lam.numerator,lam.denominator,row['attempts'])
                cohlo+=l;cohhi+=h
        fidelity=( (support[0]*size+cohlo)//(size*size),
                    ceil_div(support[1]*size+cohhi,size*size))
        # A comparison process agreeing inside radius8 can differ only after
        # at least9-width virtual moves. The original full graph is not cut.
        radius_tail=max(0,SCALE-sum(x[0] for x in row['weights'][:9-width]))
        payload={'lower':low,'upper':high}
        results.append({'intensity_time_exact':str(row['tau']),'original_attempt_count':row['attempts'],
            'full_source_probability_interval_exact':outward_pair(sum(low),sum(high)),
            'full_port_interval_sha256':sha256(json.dumps(payload,separators=(',',':')).encode()).hexdigest(),
            'omitted_binomial_tail_upper_exact':str(F(tail,SCALE)),
            'mean_graph_distance_interval_exact':outward_pair(*radial),
            'terminal_distance_at_least_four_interval_exact':outward_pair(*target),
            'first_current_record_arrival_distance_four_interval_exact':outward_pair(*arrival),
            'return_to_preparation_support_interval_exact':outward_pair(*support),
            'coherent_preparation_fidelity_interval_exact':outward_pair(*fidelity),
            'sum_absolute_offdiagonals_interval_exact':outward_pair(cohlo//size,ceil_div(cohhi,size)),
            'comparison_exit_radius_eight_probability_upper_exact':str(F(min(SCALE,radius_tail),SCALE)),
            'intervals_decimal_for_display':{'first_arrival':[float(F(x,SCALE)) for x in arrival],
                 'mean_graph_distance':[float(F(x,SCALE)) for x in radial],
                 'coherent_fidelity':[float(F(x,SCALE)) for x in fidelity]}})
    return {'level':level,'ports':g.ports,'seams':M,'packet_graph_radius':width,
            'preparation_support_size':size,'held_out_width':width==2,
            'all_ports_retained':g.ports,'graph_eccentricity_from_source':max(d),
            'coherence_support_never_leaves_preparation_support':True,'times':results}


def run():
    return {'schema':'oph.native_packet_propagation.v1',
        'verdict':'PASS: complete native packet propagation and recorded arrivals with outward bounds',
        'contract':{'process':'model/openwave_native_process; independent uniform existing seams, total load1',
            'packet':'uniform coherent superposition on the graph ball about port0; source preparation is supplied',
            'tested_widths':[0,1],'predeclared_held_out_width':2,
            'Nature_target_exposure':'No physical target enters the process or width choice; this is an internal source holdout, not a blind Nature prediction.',
            'common_clock':'tau=attempts/(2*number_of_seams), an explicit event-intensity coordinate, not physical seconds',
            'times_exact':list(map(str,TIMES)),
            'arrival':'first raw crossing-current record into graph-distance>=4; the source continues after arrival',
            'record_complement':'all no-crossing, zero-current, law/coin/padding and hold outcomes remain; no successful branch is normalized'},
        'reduction':{'identity':'T^t=sum_{k=0}^t Binomial(t,k)(6/(2M))^k(1-6/(2M))^(t-k)P^k; P=I-L/6',
            'exact_local_power_count':ORDER,'outward_fixed_point_bits':BITS,
            'precision_contract':'All intervals use floor/ceiling integer arithmetic. Decimal fields are displays, not acceptance comparisons.',
            'source_vectors':'Complete finite graph vectors are propagated; hashes bind all outward port intervals.',
            'first_arrival_flag':'unarrived evolves under Q T Q; arrived is the complementary full-source block, not an absorbed or discarded source',
            'cost':'Original attempt count and blank raw record resources are unchanged; only marginal computation is shortened.'},
        'packets':[summarize_packet(level,width) for level in (0,1,2) for width in WIDTHS],
        'physics':{'persistent_coherent_moving_packet':False,
            'reason':'The exact native offdiagonal matrix units decay without changing their two source addresses. Occupation spreads, but no coherent translation of this packet occurs.',
            'common_physical_speed_derived':False,
            'finite_volume_error':'The calculation is exact up to its stated outward arithmetic/tail errors on each complete finite graph. A comparison P kernel agreeing inside radius8 and using the same binomial scheduling weights has the displayed exit-probability bound. Different M changes those weights even at the same intensity time; no cross-size, infinite-volume or continuum identification follows.',
            'energy_work':'No source Hamiltonian or physical event-clock scale is provided by these repair records. Input preparation and blank recorded memory are supplied resources; no zero-energy or free-reset assertion is made.',
            'new_Nature_validation':False}}
