"""Exact source-derived diagonal services; an explicit admission extension."""
from fractions import Fraction as F
import json
import numpy as np
from scipy.linalg import expm
from model.certified_reference import interval as ci
from . import derived_control as beam

SCHEMA='oph.derived_occupation_diagonal.v1'
KINDS=('conditional_phase','rail_phase')


def exact_identity(kind):
    if kind not in KINDS:raise ValueError('declared diagonal source direction required')
    coefficients=(F(1),)*12 if kind=='conditional_phase' else (F(1),)*4+(F(),)*8
    spec=beam.source.specification();diagonal=[]
    wanted=((F(1),F()),)*3+((F(),F()),)*3 if kind=='conditional_phase' else (
        (F(-4,3),F(-1,3)),(F(2,3),F(-1,3)),(F(5,3),F(2,3)),
        (F(),F()),(F(),F()),(F(),F()))
    for i in range(6):
        for j in range(6):
            real=tuple(sum(c*g[i][j][0][k]for c,g in zip(coefficients,spec,strict=True))for k in (0,1))
            imag=tuple(sum(c*g[i][j][1][k]for c,g in zip(coefficients,spec,strict=True))for k in (0,1))
            if real!=(F(),F())or imag!=(wanted[i]if i==j else (F(),F())):
                raise ArithmeticError('complete exact diagonal source identity failed')
        diagonal.append(wanted[i])
    return coefficients,tuple(diagonal)


def declare(kind,angle,targets=(0,1),cell=0,level=0):
    coefficients,_=exact_identity(kind)
    value=beam.declare(angle,targets,cell,level)
    value.update(schema=SCHEMA,kind=kind,integrated_controls=[str(F(angle)*x)for x in coefficients])
    for _ in range(8):
        size=len(beam.native.canonical(value).encode())
        if size==value['resources']['canonical_declaration_bytes']:return value
        value['resources']['canonical_declaration_bytes']=size
    raise ArithmeticError('canonical declaration length did not stabilize')


def authenticate(value):
    if type(value)is not dict:raise ValueError('complete diagonal declaration required')
    try:expected=declare(value['kind'],F(value['angle']),tuple(value['local_targets']),value['cell'],value['federation_level'])
    except (KeyError,TypeError,ZeroDivisionError)as error:
        raise ValueError('typed canonical diagonal service required')from error
    if beam.native.canonical(value)!=beam.native.canonical(expected):
        raise ValueError('diagonal declaration or source evidence changed')
    return expected


def phase_enclosures(pair,theta):
    x=theta*(ci.Interval.point(pair[0])+pair[1]*ci.sqrt_fraction(F(5)))
    mid=F(x.lo+x.hi,2*ci.SCALE);radius=F(x.hi-x.lo,2*ci.SCALE)
    sine,cosine=ci.sine_cosine(mid)
    # Both derivatives have absolute value at most one.
    expansion=ci.ceiling(radius.numerator*ci.SCALE,radius.denominator)
    return (ci.Interval(cosine.lo-expansion,cosine.hi+expansion),
            ci.Interval(sine.lo-expansion,sine.hi+expansion))


def execute(value):
    value=authenticate(value);theta=F(value['angle']);_,identity=exact_identity(value['kind'])
    # Full source generator; lower complement entries are included. Its
    # binary64 approximation is checked against the exact Q(sqrt5) law below.
    full=expm(np.diag([1j*float(theta)*(float(a)+float(b)*np.sqrt(5))for a,b in identity]))
    if not np.isfinite(full).all()or any(full[i,j]!=0 for i in range(6)for j in range(6)if i!=j):
        raise ArithmeticError('complete stored diagonal source response required')
    error=F();entries=[]
    for i,pair in enumerate(identity):
        real,imag=phase_enclosures(pair,theta)
        aa,bb=F(float(full[i,i].real)),F(float(full[i,i].imag))
        error=max(error,(ci.Interval.point(aa)-real).absolute_upper()+(ci.Interval.point(bb)-imag).absolute_upper())
        entries.append((aa,bb))
    if error>beam.GATE:raise ArithmeticError('exact full-source diagonal error gate failed')
    return {'declaration':value,'stored_diagonal':tuple(entries),'complement_rows':((F(),)*4,(F(),)*4),
            'full_complement_diagonal':tuple(entries[4:]),'operator_error_upper':beam.GATE,
            'complete_halfdiamond_upper':beam.GATE*(2+beam.GATE)/2,
            'complete_map_gain_upper':(1+beam.GATE)**2,
            'exact_stored_number_preservation':True,'literal_exponential_error_diagnostic_upper':error}


def lift(value,parts,carrier):
    """Apply a complete diagonal service to exact real/imaginary operator parts.

    The native interpreter defines complex operators through two exact real
    calls. No real-only projection is applied here; all imaginary entries are
    retained explicitly in the output pair.
    """
    if type(parts)is not tuple or len(parts)!=2:raise ValueError('exact real and imaginary operator parts required')
    result=execute(value);c=beam.native.carrier(carrier)
    declaration=result['declaration']
    if c['kind']!='federation'or c['level']!=declaration['federation_level']:
        raise ValueError('authenticated finite source geometry required')
    real,imag=(beam.native._validate_density(x,c)for x in parts)
    a,b=(12*declaration['cell']+v for v in declaration['local_targets'])
    index=lambda x:2*((x>>a)&1)+((x>>b)&1)
    rr={};ii={}
    for key in set(real)|set(imag):
        x,_,y,_=key
        ax,bx=result['stored_diagonal'][index(x)];ay,by=result['stored_diagonal'][index(y)]
        factor_r=ax*ay+bx*by;factor_i=bx*ay-ax*by
        out_r=real.get(key,F())*factor_r-imag.get(key,F())*factor_i
        out_i=real.get(key,F())*factor_i+imag.get(key,F())*factor_r
        if out_r:rr[key]=out_r
        if out_i:ii[key]=out_i
    return {'success':(rr,ii),'complements':({'label':4,'parts':({}, {})},{'label':5,'parts':({}, {})}),
            'absolute_targets':(a,b),'complete_halfdiamond_upper':result['complete_halfdiamond_upper'],
            'complete_map_gain_upper':result['complete_map_gain_upper']}
