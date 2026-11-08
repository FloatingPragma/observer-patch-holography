"""Rational residual certificates for the source-derived spectral logarithm."""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from tools import strict_json
from model.native_record_response import response, transfer

DIM=66
PROPOSAL=Path(__file__).with_name('proposal.json')
GRAM_GATE=F(1,10**9)
RESIDUAL_GATE=F(1,10**9)


def fraction(x):
    if type(x) in (int,F):return F(x)
    if hasattr(x,'p') and hasattr(x,'q'):return F(int(x.p),int(x.q))
    raise ValueError('an exact rational number is required')


def load(path=PROPOSAL):
    p=strict_json.load(path)
    if set(p)!={'schema','denominator','columns','eigenvalues'} or p['schema']!='oph.native_transfer_control.proposal.v1':
        raise ValueError('the typed spectral proposal schema is required')
    d=p['denominator']; cols=p['columns']; values=p['eigenvalues']
    if type(d) is not int or d!=2**44 or type(cols) is not list or len(cols)!=DIM or type(values) is not list or len(values)!=DIM:
        raise ValueError('the fixed rational proposal dimensions are required')
    if any(type(c) is not list or len(c)!=DIM or any(type(x) is not int or abs(x)>d for x in c) for c in cols):
        raise ValueError('bounded exact integer proposal columns are required')
    if any(type(x) is not int or not transfer.LOWER*d<=x<=d for x in values):
        raise ValueError('positive rational proposed transfer values are required')
    return tuple(tuple(F(x,d) for x in c) for c in cols),tuple(F(x,d) for x in values)


def source_intervals(data=None):
    data=response.data() if data is None else data
    P=data['P']; pi=tuple(map(fraction,data['pi']))
    if (P.shape!=(DIM,DIM) or len(pi)!=DIM or sum(pi)!=1
            or any(x<=0 for x in pi) or min(P[i,i] for i in range(DIM))<F(791,990)):
        raise ArithmeticError('the native stationary source contract failed')
    if any(sum((fraction(P[i,j]) for i in range(DIM)),F())!=1 for j in range(DIM)):
        raise ArithmeticError('the source population transfer is not stochastic')
    if any(fraction(P[i,j])*pi[j]!=fraction(P[j,i])*pi[i] for i in range(DIM) for j in range(DIM)):
        raise ArithmeticError('the actual stationary source is not reversible')
    # Every sqrt(7505*pi_i) is 25, 10 sqrt(2), or 8.
    q=2**100; k=isqrt(2*q*q); root=(F(k,q),F(k+1,q))
    scale=[]
    for p in pi:
        z=7505*p
        if z==625:scale.append((F(25),F(25)))
        elif z==200:scale.append((10*root[0],10*root[1]))
        elif z==64:scale.append((F(8),F(8)))
        else:raise ArithmeticError('the authoritative weighted reference changed')
    rows=[]
    for i in range(DIM):
        row=[]
        for j in range(DIM):
            p=fraction(P[i,j])
            if p<0:raise ArithmeticError('the native transfer lost positivity')
            if p:row.append((j,p*scale[j][0]/scale[i][1],p*scale[j][1]/scale[i][0]))
        rows.append(tuple(row))
    return tuple(rows),pi


def certify(columns=None,values=None,data=None):
    if columns is None and values is None:columns,values=load()
    elif columns is None or values is None:raise ValueError('both proposal components are required')
    if (len(columns)!=DIM or len(values)!=DIM or any(len(c)!=DIM for c in columns)
            or any(type(x) is not F for c in columns for x in c)
            or any(type(x) is not F or not transfer.LOWER<=x<=1 for x in values)):
        raise ValueError('exact full source proposal components are required')
    rows,pi=source_intervals(data)
    gram=[[sum((a*b for a,b in zip(columns[i],columns[j],strict=True)),F())-int(i==j)
           for j in range(DIM)] for i in range(DIM)]
    g=max(sum(map(abs,row),F()) for row in gram)
    residual=[]
    for i,row in enumerate(rows):
        bound=[]
        for c,lam in zip(columns,values,strict=True):
            lo=hi=-lam*c[i]
            for j,a,b in row:
                x,y=(a*c[j],b*c[j]) if c[j]>=0 else (b*c[j],a*c[j])
                lo+=x;hi+=y
            bound.append(max(abs(lo),abs(hi)))
        residual.append(bound)
    r=max(max(sum(row,F()) for row in residual),max(sum(col,F()) for col in zip(*residual,strict=True)))
    if g>GRAM_GATE or r>RESIDUAL_GATE:
        raise ArithmeticError('the proposal does not certify the native transfer')
    return {'gram_operator_upper':str(GRAM_GATE),'residual_operator_upper':str(RESIDUAL_GATE),
            'positive_floor':str(transfer.LOWER),'source_dimension':DIM,
            'spectral_eigenvalue_gaps_required':False,'polar_basis_is_executed':False}


def time_parameter(value):
    if type(value) is not F or abs(value)>1 or max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64:
        raise ValueError('a 64-bit rational transfer parameter in [-1,1] is required')
    return value


def log_interval(value,degree=96):
    if type(value) is not F or not transfer.LOWER<=value<=1 or type(degree) is not int or degree<1:
        raise ValueError('a positive bounded rational transfer value and degree are required')
    x=1-value; total=sum((x**n/n for n in range(1,degree+1)),F())
    tail=x**(degree+1)/((degree+1)*(1-x))
    return total,total+tail


def phase_grid(value,time):
    from model.native_effect_synthesis.intervals import AQ
    time=time_parameter(time);a,b=log_interval(value)
    lo,hi=sorted((-time*a,-time*b));n=round((lo+hi)*AQ/2)
    if max(abs(F(n,AQ)-lo),abs(F(n,AQ)-hi))>F(1,AQ):
        raise ArithmeticError('the rational logarithm phase exceeds its grid allowance')
    return n


def ideal_basis_error(time):
    time=time_parameter(time)
    # Polar comparison is a proof device; actual words use normalized rational rays.
    return abs(time)*(RESIDUAL_GATE+2*GRAM_GATE)/transfer.LOWER+8*DIM*GRAM_GATE
