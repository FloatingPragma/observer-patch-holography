"""Full signed amplitude tree using Gray-code uniformly controlled Y words."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math
import mpmath
import numpy as np

from ..observer.source_scalar import compiler
from ..certified_reference import certificate as cert, source, interval as iv
from . import target

BITS=16


def walsh(values):
    result=list(values);step=1
    if not result or len(result)&(len(result)-1):raise ValueError('power of two required')
    while step<len(result):
        for start in range(0,len(result),2*step):
            for j in range(start,start+step):
                a,b=result[j],result[j+step];result[j],result[j+step]=a+b,a-b
        step*=2
    return result


def exact_floats(values):
    ratios=[float(v).as_integer_ratio() for v in values]
    denominator=max(d for _,d in ratios)
    return [n*(denominator//d) for n,d in ratios],denominator


def design(values,bits):
    if type(bits) is not int or not 1<=bits<=16 or len(values)!=1<<bits or any(type(v) is not int for v in values) or not any(values):
        raise ValueError('complete nonzero integer amplitude vector and 1..16 bits required')
    prefix=[0]
    for v in values:prefix.append(prefix[-1]+v*v)
    context=mpmath.mp.clone();context.prec=192
    levels=[]
    for bit in reversed(range(bits)):
        raw=[]
        for label in range(1<<(bits-bit-1)):
            start=label*(1<<(bit+1));half=1<<bit
            if bit==0:
                a,b=values[start:start+2]
                x,y=context.mpf(a),context.mpf(b)
            else:
                a=prefix[start+half]-prefix[start];b=prefix[start+2*half]-prefix[start+half]
                x,y=context.sqrt(a),context.sqrt(b)
            raw.append(float(-context.atan2(y,x)) if x or y else 0.)
        numerators,denominator=exact_floats(raw)
        coefficients=np.array([float(F(v,denominator*len(raw))) for v in walsh(numerators)])
        raw=np.array(raw);raw.flags.writeable=False;coefficients.flags.writeable=False
        levels.append((bit,raw,coefficients))
    return tuple(levels)


@lru_cache(maxsize=1)
def ground_design():return design(target.ground()[0],BITS)


def level_word(bit,coefficients,bits=BITS):
    controls=bits-bit-1;n=1<<controls
    if len(coefficients)!=n:raise ValueError('complete uniformly controlled level required')
    for j in range(n):
        gray=j^(j>>1)
        yield compiler.pulse('Y',float(coefficients[gray]),[bit])
        if controls:
            following=((j+1)%n)^(((j+1)%n)>>1)
            toggle=gray^following
            if not toggle or toggle&(toggle-1):raise ArithmeticError('Gray cycle lost one-bit transition')
            control=bit+1+toggle.bit_length()-1
            yield from compiler.cnot(control,bit)


def word():
    for bit,_,coefficients in ground_design():yield from level_word(bit,coefficients)


def blocks(events,max_pulses=32):
    """Partition the actual tape; never replace its rotations by ideal gates."""
    from ..protected_particle_readout import execute
    pending=[];support=set();offset=0
    def pack():
        targets=tuple(sorted(support));index={v:j for j,v in enumerate(targets)}
        tape=tuple((e['kind'],e['angle'],tuple(index[t] for t in e['targets'])) for e in pending)
        return execute.compile_block(tape,len(targets)),targets
    for event in events:
        new=support|set(event['targets'])
        if pending and (len(new)>3 or len(pending)>=max_pulses):
            instrument,targets=pack();yield offset,instrument,targets
            offset+=len(pending);pending=[];support=set();new=set(event['targets'])
        pending.append(event);support=new
    if pending:
        instrument,targets=pack();yield offset,instrument,targets


def level_error(values,bit,raw,coefficients,bits):
    """Raw angles are checked, not trusted; rounded Gray coefficients are inverted exactly."""
    if len(raw)!=1<<(bits-bit-1) or len(coefficients)!=len(raw):raise ValueError('level dimensions changed')
    numerators,denominator=exact_floats(coefficients)
    recovered=walsh(numerators)
    prefix=[0]
    for v in values:prefix.append(prefix[-1]+v*v)
    maximum=F();coefficient_error=F()
    for label,angle in enumerate(raw):
        start=label*(1<<(bit+1));half=1<<bit
        a,b=(values[start],values[start+1]) if bit==0 else (prefix[start+half]-prefix[start],prefix[start+2*half]-prefix[start+half])
        maximum=max(maximum,cert.two_vector_error(F(float(angle)),a,b,bit==0))
        coefficient_error=max(coefficient_error,abs(F(recovered[label],denominator)-F(float(angle))))
    return maximum,coefficient_error


@lru_cache(maxsize=1)
def generator_norms():
    return tuple(max(sum(source.field_interval(part).absolute_upper() for cell in row for part in cell) for row in matrix)
                 for matrix in source.specification())


def port_bound(counts):
    """Exact scalar reduction of full-generator errors for many different angles.

    G_actual(theta)-theta*T = theta*(G_base-T)+sum(delta_p D_p).
    delta_p uses actual binary64 control outputs, with no rounding-model assumption.
    """
    norms=generator_norms();bases={};total=F()
    for (kind,angle),count in counts.items():
        if kind not in bases:
            base=source.literal_controls(kind,1.)
            bases[kind]=(base,source.generator_error(kind,F(1),base))
        base,error=bases[kind];theta=F(angle)
        actual=source.literal_controls(kind,angle)
        bound=abs(theta)*error+sum(abs(a-theta*b)*norm for a,b,norm in zip(actual,base,norms,strict=True))
        total+=count*bound
    return total


@lru_cache(maxsize=1)
def certificate():
    values,_,_=target.ground();rows=[];count=0;counts=Counter();digest=__import__('hashlib').sha256()
    for bit,raw,coefficients in ground_design():
        angle,rounding=level_error(values,bit,raw,coefficients,BITS)
        pulses=0
        for event in level_word(bit,coefficients):
            counts[event['kind'],event['angle']]+=1;pulses+=1
            digest.update((event['kind']+'|'+event['angle'].hex()+'|'+','.join(map(str,event['targets']))+';').encode())
        count+=pulses
        rows.append({'target':bit,'conditional_nodes':len(raw),'pulses':pulses,
                     'signed_node_error':str(angle),'Gray_coefficient_error':str(rounding)})
    cnots=(1<<BITS)-2
    fixed=cnots*(2*iv.pi_angle_error(F(math.pi/2),F(1,2))+iv.pi_angle_error(F(-math.pi),F(-1)))
    native=port_bound(counts)
    error=fixed+native+sum((F(r['signed_node_error'])+F(r['Gray_coefficient_error']) for r in rows),F())
    if count!=4*((1<<BITS)-1)-3:raise ArithmeticError('complete Gray tree pulse count changed')
    if error>=F(1,10**7):raise ArithmeticError('complete tree misses declared precision')
    return {'bits':BITS,'pulses':count,'CNOTs':cnots,'conditional_nodes':(1<<BITS)-1,
            'word_sha256':digest.hexdigest(),'levels':rows,'fixed_CNOT_error':str(fixed),
            'actual_port_generator_error':str(native),'complete_tree_error':str(error),
            'unique_kind_angles':len(counts),'all_zero_rotation_requests_retained':True,
            'Gray_cycle_final_CNOT_retained':True,'full_six_level_generators_used':True}
