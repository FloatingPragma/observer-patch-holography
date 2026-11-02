"""Refinable dyadic port commands from the exact readable source identities.

This command format is deliberately not the old binary64 pulse interface.
All approximation bounds here concern exact rational integrated port inputs.
The separately executed legacy words retain their actual binary64 error.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from ..observer.source_scalar import compiler
from ..certified_reference import source, certificate as cert


def add(a,b): return a[0]+b[0],a[1]+b[1]
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)


def sqrt_interval(n,bits):
    if type(n) is not int or n<0 or type(bits) is not int or bits<1:
        raise ValueError('nonnegative integer radicand and positive bit precision required')
    scale=1<<bits;lo=isqrt(n*scale*scale)
    return F(lo,scale),F(lo if lo*lo==n*scale*scale else lo+1,scale)


def atan_interval(inverse,bits):
    if type(inverse) is not int or inverse<2 or type(bits) is not int or bits<1:
        raise ValueError('integer inverse >=2 and positive bit precision required')
    total=F();k=0;target=F(1,1<<bits)
    while True:
        term=F((-1)**k,(2*k+1)*inverse**(2*k+1));total+=term;k+=1
        following=F((-1)**k,(2*k+1)*inverse**(2*k+1))
        if abs(following)<=target:return min(total,total+following),max(total,total+following)


@lru_cache(maxsize=128)
def pi_interval(bits):
    # Machin: pi=16 atan(1/5)-4 atan(1/239), verified by tangent addition
    # and 3<pi<4; the alternating remainders have their displayed signs.
    if type(bits) is not int or bits<1:raise ValueError('positive integer precision required')
    return add(mul((F(16),F(16)),atan_interval(5,bits+5)),
               mul((F(-4),F(-4)),atan_interval(239,bits+5)))


def angle_interval(angle,bits):
    if type(angle) is F:return angle,angle
    if type(angle) is cert.PiMultiple and type(angle.value) is F:return mul((angle.value,angle.value),pi_interval(bits))
    raise ValueError('exact rational angle or rational multiple of pi required')


@lru_cache(maxsize=1)
def generator_norms():
    # sqrt(5)<3 gives an entirely rational absolute-row majorant.
    return tuple(max(sum(abs(p[0])+3*abs(p[1]) for cell in row for p in cell)
                     for row in matrix) for matrix in source.specification())


def coefficient_intervals(kind,bits):
    if kind not in compiler.KINDS:raise ValueError('declared native kind required')
    source.specification()  # independently recheck all exact source identities
    exact=compiler._source()[2];s5=sqrt_interval(5,bits)
    if kind=='H':
        coefficients=[(x[0]+z[0],x[1]+z[1]) for x,z in zip(exact['X'],exact['Z'],strict=True)]
    else:coefficients=exact[kind]
    values=[add((a,a),mul((b,b),s5)) for a,b in coefficients]
    if kind=='H':
        a,b=sqrt_interval(2,bits);values=[mul(v,(1/b,1/a)) for v in values]
    return values


def dyadic_command(kind,angle,bits):
    """Return twelve exact dyadics, each within 2^-bits of its ideal port.

    Interval precision is increased until the requested property is proved;
    there is no floating-point stopping test or fixed-precision floor.
    """
    if type(bits) is not int or bits<1:raise ValueError('positive integer output bits required')
    target=F(1,1<<bits);working=bits+8
    while True:
        a=angle_interval(angle,working)
        intervals=[mul(a,c) for c in coefficient_intervals(kind,working)]
        scale=1<<bits;commands=[];errors=[]
        for lo,hi in intervals:
            midpoint=(lo+hi)/2;rounded=(midpoint*scale+F(1,2)).numerator//(midpoint*scale+F(1,2)).denominator
            q=F(rounded,scale);commands.append(q);errors.append(max(abs(q-lo),abs(q-hi)))
        if max(errors)<=target:break
        working*=2
    bound=sum(e*n for e,n in zip(errors,generator_norms(),strict=True))
    return {'kind':kind,'fraction_bits':bits,'working_interval_bits':working,
            'integrated_port_coefficients':[str(q) for q in commands],
            'port_error_bounds':[str(e) for e in errors],
            'generator_and_complete_output_error':str(bound),
            'uniform_error_ceiling':str(sum(generator_norms())*target)}


def bits_for(pulses,trace_budget):
    if type(pulses) is not int or pulses<1 or type(trace_budget) is not F or not 0<trace_budget<2:
        raise ValueError('positive pulse count and exact trace budget in (0,2) required')
    bits=1;scale=2*pulses*sum(generator_norms())
    while scale/F(1<<bits)>trace_budget:bits+=1
    return bits


def certificate():
    examples=[dyadic_command(k,cert._pi(1,4),p) for k in ('X','H','J') for p in (24,80,160)]
    return {'exact_source':'twelve anti-Hermitian Q(sqrt(5)) generators and checked target identities',
            'generator_norm_majorants':[str(v) for v in generator_norms()],
            'sum_generator_norm_majorants':str(sum(generator_norms())),
            'coefficient_algorithm':'integer-square-root enclosures, Machin alternating series, interval refinement and dyadic midpoint rounding',
            'per_pulse_operator_error':'B*2^(-p), B=sum port-generator norm majorants',
            'N_pulse_complete_channel_trace_error':'2*N*B*2^(-p)',
            'examples':examples,'fixed_binary64_error_sent_to_zero':False,
            'refined_rational_commands_executed_by_legacy_binary64_API':False,
            'refined_source_interface':'exact rational integrated port coefficients in the same exponential source law; physical analogue precision remains an apparatus requirement'}
