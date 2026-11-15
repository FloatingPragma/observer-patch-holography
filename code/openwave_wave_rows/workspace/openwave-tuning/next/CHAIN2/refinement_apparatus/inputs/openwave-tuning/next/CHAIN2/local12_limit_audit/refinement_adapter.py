"""Additive finite-member native LOCAL12 adapter with scalable rational controls.

The frozen size-two compiler is never modified. This adapter retains its
native event/recoder and M6 proper-code services, but declares an exact-
rational control serialization instead of its binary64-only pulse API.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PROVIDER = ROOT / 'openwave-tuning/auditor/native_walk_bridge'
sys.path.insert(0, str(ROOT / 'reverse-engineering-reality'))
sys.path.insert(0, str(PROVIDER))
import common_walk as common
import compile_walk as frozen
from model.native_spatial_refinement import source as spatial
from model.openwave_coherent_native_architecture import source as architecture
from model.observer.source_scalar import compiler
from model.certified_reference import certificate as templates
from model.certified_reference import source as certified_source

ADAPTER_ID = 'native-local12-rational-refinement-adapter-v1'


def layout(size):
    if type(size) is not int or size < 2:
        raise ValueError('integer cubic member size at least two required')
    data = 8 + (size**3 - 1).bit_length()
    # Routing scratch uses data-3 bits; endpoint extraction uses data-3.
    return {'size': size, 'data_bits': data, 'source': tuple(range(data, data+3)),
            'raw': tuple(range(data+3, data+6)), 'law': data+6,
            'scratch': tuple(range(data+7, 2*data+4)), 'memory_bits': 2*data+4}


def member_geometry(size):
    """Actual pinned graph, with all padded addresses retained."""
    g = spatial.graph(size)
    d = layout(size)['data_bits']
    target = list(range(1 << d))
    edges = []
    colour = {tuple(edge): c for c,row in enumerate(common.INTERNAL_MATCHINGS) for edge in row}
    for a,b in g.seams:
        direction = colour[(a%12,b%12)] if a//12 == b//12 else 5
        x,y = 256*(a//12)+16*(a%12)+2*direction+1, 256*(b//12)+16*(b%12)+2*direction+1
        target[x],target[y] = target[y],target[x]
        edges.append((a,b,direction))
    return g, target, edges


def extraction_gates(size):
    info = layout(size)
    _,_,edges = member_geometry(size)
    gates = []
    for a,b,direction in edges:
        for vertex, target in ((a,info['source'][1]), (b,info['source'][0])):
            position = 16*(vertex//12)+vertex%12
            values = {j+1:(direction >> j)&1 for j in range(3)}
            values.update({j+4:(position >> j)&1 for j in range(info['data_bits']-4)})
            gates += frozen.multi_control(values, target, info['scratch'])
    return gates


@lru_cache(None)
def reversible_program(size):
    info = layout(size)
    _, target, _ = member_geometry(size)
    extraction = extraction_gates(size)
    route = architecture.compile_permutation(target)
    remap = {q:q for q in range(info['data_bits'])}
    remap.update({info['data_bits']+j:q for j,q in enumerate(info['scratch'])})
    route_gates = frozen.remap_gates(route['gates'],remap)
    local_map = {0:info['source'][0],1:info['source'][1],2:info['source'][2],
                 3:info['raw'][0],4:info['raw'][1],5:info['raw'][2],6:0,7:info['law']}
    local_map.update({8+j:info['scratch'][j] for j in range(5)})
    native = frozen.remap_gates(architecture.compile_permutation(architecture.event_permutation())['gates'], local_map)
    recoder = frozen.remap_gates(architecture.compile_permutation(architecture.record_recoder())['gates'], local_map)
    gates = extraction + [('cx',(info['source'][1],0))] + native + recoder + route_gates + list(reversed(extraction))
    readable = extraction + list(reversed(recoder)) + list(reversed(extraction))
    return {'layout':info,'gates':gates,'inverse_record_recoder':readable,
            'routing_gates':route_gates,'routing_permutation':target,
            'component_gate_counts':{'extraction':len(extraction),'native_event':len(native),
                'recoder':len(recoder),'routing':len(route_gates)},
            'all_padding_and_dirty_scratch_defined':True}


def execute_bitplanes(gates,inputs,memory_bits):
    """Exact integer execution, packing columns as bits instead of float arrays."""
    values=list(map(int,inputs))
    if len(set(values))!=len(values) or any(x<0 or x>=1<<memory_bits for x in values):
        raise ValueError('distinct complete-register columns required')
    planes=[sum(1<<j for j,x in enumerate(values) if x>>q&1) for q in range(memory_bits)]
    mask=(1<<len(values))-1
    for kind,targets in gates:
        if kind=='flip': planes[targets[0]]^=mask
        elif kind=='cx': planes[targets[1]]^=planes[targets[0]]
        elif kind=='ccx': planes[targets[2]]^=planes[targets[0]]&planes[targets[1]]
        else: raise ValueError('unknown reversible primitive')
    dtype=np.int64 if memory_bits<63 else object
    return np.array([sum(((plane>>j)&1)<<q for q,plane in enumerate(planes)) for j in range(len(values))],dtype=dtype)


# Exact rational outward intervals. Precision is an API argument, not a fixed
# 160-bit ceiling or a host float conversion.
def round_interval(v, bits):
    scale=1 << bits
    return (F(v[0].numerator*scale//v[0].denominator,scale),
            F(-((-v[1].numerator*scale)//v[1].denominator),scale))


def add(a,b,bits): return round_interval((a[0]+b[0],a[1]+b[1]),bits)
def neg(a): return -a[1],-a[0]
def mul(a,b,bits):
    products=[x*y for x in a for y in b]
    return round_interval((min(products),max(products)),bits)
def scale(a,q,bits): return mul(a,(F(q),F(q)),bits)
def div(a,b,bits):
    if b[0] <= 0: raise ValueError('positive interval divisor required')
    return mul(a,(1/b[1],1/b[0]),bits)


def sqrt_interval(value,bits):
    value=F(value)
    if value<0: raise ValueError('nonnegative radicand required')
    scale=1 << bits
    below=isqrt(value.numerator*scale*scale//value.denominator)
    exact=below*below*value.denominator == value.numerator*scale*scale
    return F(below,scale),F(below if exact else below+1,scale)


def atan_interval(x,bits):
    if not 0 <= x[0] <= x[1] < 1: raise ValueError('atan series requires [0,1)')
    def endpoint(q):
        total=(F(0),F(0)); power=(q,q); q2=(q*q,q*q); j=0
        while True:
            total=add(total,scale(power,F((-1)**j,2*j+1),bits),bits)
            next_power=mul(power,q2,bits)
            tail=next_power[1]/(2*j+3)
            if tail < F(1,1 << (bits-4)):
                return round_interval((total[0]-tail,total[1]+tail),bits)
            power=next_power; j+=1
    return endpoint(x[0])[0],endpoint(x[1])[1]


@lru_cache(None)
def constants(bits):
    if type(bits) is not int or bits<32: raise ValueError('at least32 precision bits required')
    pi=add(scale(atan_interval((F(1,5),F(1,5)),bits),16,bits),
           scale(atan_interval((F(1,239),F(1,239)),bits),-4,bits),bits)
    root2=sqrt_interval(2,bits)
    ratio=div(add(root2,(F(-1),F(-1)),bits),add(root2,(F(1),F(1)),bits),bits)
    alpha=add(scale(pi,F(1,4),bits),atan_interval(ratio,bits),bits)
    return {'pi':pi,'alpha':alpha,'root2':root2,'root5':sqrt_interval(5,bits)}


# Angles are rational linear combinations q_pi*pi + q_alpha*atan(sqrt2).
def angle(qpi=0,qalpha=0,rational=0): return F(qpi),F(qalpha),F(rational)
def negate_template(word): return [(k,t,tuple(-v for v in a)) for k,t,a in reversed(word)]
def convert_template(word):
    return [(k,t,angle(a.value) if isinstance(a,templates.PiMultiple) else a) for k,t,a in word]
def pauli(label,a): return convert_template(templates.pauli(label,a))
def cnot(a,b): return convert_template(templates.cnot(a,b))
def diagonal(phases):
    bits=(len(phases)-1).bit_length(); result=[]
    for mask in range(len(phases)):
        q=sum((F(x)*(-1)**((i&mask).bit_count()) for i,x in enumerate(phases)),F())/len(phases)
        result += pauli(''.join('Z' if mask>>j&1 else 'I' for j in range(bits)),angle(q))
    return result


@lru_cache(None)
def primitive_templates():
    return {'flip':[('I',(0,),angle(F(-1,2))),('X',(0,),angle(F(1,2)))],
            'hadamard':[('I',(0,),angle(F(-1,2))),('H',(0,),angle(F(1,2)))],
            'cx':cnot(0,1),
            'ccx':[('H',(2,),angle(F(1,2)))]+diagonal([0]*7+[1])+[('H',(2,),angle(F(-1,2)))]}


@lru_cache(None)
def uniform_template():
    result=[]
    for a,b,theta in ((0,1,angle(0,-1)),(1,3,angle(F(-1,4))), (3,2,angle(F(-1,2)))):
        bit=(a^b).bit_length()-1; orient=(-1)**((a>>bit)&1)
        for mask in range(4):
            if mask & (a^b): continue
            q=F(orient*(-1)**((a&mask).bit_count()),2)
            coefficient=tuple(q*v for v in theta)
            label=''.join('Y' if j==bit else 'Z' if mask>>j&1 else 'I' for j in range(2))
            result += [(k,tuple(j+2 for j in t),v) for k,t,v in pauli(label,coefficient)]
    for bit in (0,1): result += [(k,(bit,),a) for k,t,a in primitive_templates()['hadamard']]
    return result


def reflection_template():
    u=uniform_template()
    return negate_template(u)+diagonal([0]+[1]*15)+u


@lru_cache(None)
def source_constants(bits):
    certified_source.specification()  # exact all-ten control identities
    roots=constants(bits)
    coefficients={}
    for k,rows in compiler.source_summary()['exact_base_coefficients'].items():
        coefficients[k]=[add((F(a),F(a)),scale(roots['root5'],F(b),bits),bits) for a,b in rows]
    coefficients['H']=[div(add(a,b,bits),roots['root2'],bits) for a,b in zip(coefficients['X'],coefficients['Z'])]
    # Full-six-level anti-Hermitian source norm upper bounds, hence Duhamel.
    norms=[]
    for generator in certified_source.specification():
        rows=[]
        for row in generator:
            bound=F()
            for re,im in row:
                for a,b in (re,im):
                    lo,hi=add((a,a),scale(roots['root5'],b,bits),bits)
                    bound+=max(abs(lo),abs(hi))
            rows.append(bound)
        norms.append(max(rows))
    return coefficients,norms


@lru_cache(None)
def certified_pulse(kind,targets,desired,bits):
    """Emit twelve exact rational controls and a full-response norm bound."""
    if kind not in ('X','Y','Z','I','H','J') or len(targets)!=(2 if kind=='J' else 1) or len(set(targets))!=len(targets) or min(targets)<0:
        raise ValueError('native kind and distinct nonnegative targets required')
    c=constants(bits)
    theta=add(scale(c['pi'],desired[0],bits),scale(c['alpha'],desired[1],bits),bits)
    constant=desired[2] if len(desired)==3 else F()
    theta=add(theta,(constant,constant),bits)
    coefficients,norms=source_constants(bits)
    actual=[]; error=F()
    for base,norm in zip(coefficients[kind],norms):
        lo,hi=mul(theta,base,bits)
        center=(lo+hi)/2
        actual.append(center)
        error+=max(center-lo,hi-center)*norm
    return {'kind':kind,'targets':list(targets),'controls_exact':list(map(str,actual)),
            'desired_angle_exact':list(map(str,desired)),
            'angle_interval_exact':list(map(str,theta)), 'full_response_operator_error_exact':str(error),
            'precision_bits':bits, 'export':'coherent proper-code export with every complement and first-failure record retained'}


def verify_pulse(event):
    """Reject coefficient changes before relying on the emitted certificate."""
    kind=event['kind']; bits=event['precision_bits']; desired=tuple(map(F,event['desired_angle_exact']))
    c=constants(bits)
    theta=add(scale(c['pi'],desired[0],bits),scale(c['alpha'],desired[1],bits),bits)
    constant=desired[2] if len(desired)==3 else F()
    theta=add(theta,(constant,constant),bits)
    coefficients,norms=source_constants(bits); expected=[]; bound=F()
    for base,norm in zip(coefficients[kind],norms):
        lo,hi=mul(theta,base,bits); midpoint=(lo+hi)/2
        expected.append(midpoint); bound+=(hi-lo)*norm/2
    if (list(map(F,event['controls_exact']))!=expected
            or tuple(map(F,event['angle_interval_exact']))!=theta
            or F(event['full_response_operator_error_exact'])!=bound):
        raise ArithmeticError('native controls differ from their exact emitted certificate')
    return True


def template_bound(word,bits):
    return sum((F(certified_pulse(k,tuple(range(2 if k=='J' else 1)),a,bits)['full_response_operator_error_exact']) for k,t,a in word),F())


def certificate(size,step_tolerance):
    """Compile precision for a finite member; delta_W can be set to h^2."""
    tol=F(step_tolerance)
    if tol<=0: raise ValueError('positive per-word tolerance required')
    program=reversible_program(size)
    counts=Counter(k for k,t in program['gates'])
    primitive=primitive_templates(); reflection=reflection_template()
    pulses=len(reflection)+sum(n*len(primitive[k]) for k,n in counts.items())
    bits=max(48,(pulses*tol.denominator//tol.numerator+1).bit_length()+32)
    while True:
        errors={k:template_bound(word,bits) for k,word in primitive.items()}
        bound=template_bound(reflection,bits)+sum(n*errors[k] for k,n in counts.items())
        if bound<=tol: break
        bits*=2
    return {'adapter_id':ADAPTER_ID,'size':size,'layout':program['layout'],'precision_bits':bits,
            'native_pulses':pulses,'load_export_services':2*pulses,
            'reversible_gate_counts':dict(counts),'component_gate_counts':program['component_gate_counts'],
            'whole_W_operator_error_exact':str(bound),'requested_step_tolerance_exact':str(tol),
            'preparation_uniform12_error_exact':str(template_bound(uniform_template(),bits)),
            'complete_stinespring_error_bound':str(bound), 'complete_channel_TV_and_half_diamond_upper':str(min(F(1),bound)),
            'physical_duration_claimed':False,'serialization':'twelve exact rational integrated native controls per pulse; not the frozen binary64 API'}


def pulse_stream(size,precision_bits):
    """Actual chronological W generator; no ideal matrix is an execution input."""
    for k,t,a in reflection_template(): yield certified_pulse(k,t,a,precision_bits)
    primitive=primitive_templates()
    for name,targets in reversible_program(size)['gates']:
        for k,local,a in primitive[name]:
            yield certified_pulse(k,tuple(targets[j] for j in local),a,precision_bits)


def execute_full_response(event):
    """Numerical development evaluator of emitted controls, all six levels.

    Certification uses exact intervals, not this float diagnostic. Return
    every coherent export Kraus branch on the admitted proper-code input.
    """
    from scipy.linalg import expm
    from math import sqrt
    controls=tuple(F(x) for x in event['controls_exact'])
    if len(controls)!=12: raise ValueError('all twelve controls required')
    spec=certified_source.specification()
    full=np.zeros((6,6),complex)
    for c,g in zip(controls,spec):
        full+=float(c)*np.array([[complex(float(r[0])+sqrt(5)*float(r[1]),float(i[0])+sqrt(5)*float(i[1])) for r,i in row] for row in g])
    U=expm(full); d=4 if event['kind']=='J' else 2
    return [K@U[:,:d] for K in compiler.transfer_kraus(d,True)]


def reset_repetition_plan(registers,total_timeout):
    """Uniform all78-input reset certificate, repeated without postselection."""
    if type(registers) is not int or registers<1 or not 0<F(total_timeout)<1: raise ValueError('register count and positive failure budget below one required')
    receipt=json.loads((PROVIDER/'COMMON_WALK_RECEIPT.json').read_text())
    base=receipt['preparation_reset']['certificate']
    q=F(base['failure_upper'])
    repeats=1
    while registers*q**repeats>F(total_timeout): repeats+=1
    return {'registers':registers,'blocks_per_register':repeats,'base_uniform_failure_exact':str(q),
            'complete_timeout_union_bound_exact':str(registers*q**repeats),
            'protocol':'Each register executes up to r base reset blocks, stopping at its first heralded success; retain every fine letter and final timeout. A failed block remains in the native78-state domain.',
            'conditioning_applied':False}


def resource_upper_bound(size):
    """Symbolic family estimate; does not construct or execute a graph."""
    d=layout(size)['data_bits']
    local=architecture.compile_permutation(architecture.event_permutation())['gates']
    recode=architecture.compile_permutation(architecture.record_recoder())['gates']
    # There are36L^3 disjoint moving pairs and72L^3 extraction endpoints.
    gates=144*size**3*(4*d-7)+36*size**3*(2*d-1)*(4*d-7)+len(local)+len(recode)+1
    per_gate=max(map(len,primitive_templates().values()))
    return {'size':size,'data_bits':d,'memory_bits':2*d+4,
            'reversible_gates_upper':gates,'native_pulses_upper':len(reflection_template())+per_gate*gates,
            'growth':'O(L^3 log^2 L) native pulses; O(log L) memory qubits',
            'enumerated':False}


def clock_tick_pulse(size,cell_spacing,precision_bits,clock_wire=None):
    """Fresh common-clock pulse after the entire completed member word."""
    h=F(cell_spacing); memory=layout(size)['memory_bits']
    wire=memory if clock_wire is None else clock_wire
    if h<=0 or type(wire) is not int or wire<memory:
        raise ValueError('positive h and fresh clock address outside full source required')
    return certified_pulse('Z',(wire,),angle(0,0,-h/16),precision_bits)
