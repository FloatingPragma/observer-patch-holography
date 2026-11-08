"""Rational source-mesh controls for a translated compact C3 preparation."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as s


def bump(x):return (1-x*x)**4 if abs(x)<1 else F(0)


@lru_cache(maxsize=1)
def integrals():
    x=s.symbols('x');b=(1-x*x)**4
    i0=F(s.integrate(b*b,(x,-1,1)))
    i1=F(s.integrate(s.diff(b,x)**2,(x,-1,1)))
    return i0,i1


def control(channel, divisions, center=(F(1),F(0),F(0)), axis=0):
    if channel not in ('Q','P')or type(divisions)is not int or divisions<2:
        raise ValueError('Q/P and integer mesh divisions at least two required')
    if type(axis)is not int or axis not in range(3)or len(center)!=3:
        raise ValueError('three-dimensional source direction required')
    if any(type(c)is not F or (c*divisions).denominator!=1 for c in center):
        raise ValueError('exact lattice-aligned rational center required')
    if sum(int(c*divisions)for c in center)%2:
        raise ValueError('selected-state even source translation required')
    p=F(1,divisions);tables=[]
    for c in center:
        start=int(c*divisions)-divisions-2;stop=int(c*divisions)+divisions+3
        values=[(p*j,bump(p*j-c),bump(p*(j+1)-c))for j in range(start,stop)]
        table={}
        for power in (2,4):
            for edge in (False,True):
                entries=[]
                for x,a,b in values:
                    value=(a+b if channel=='Q'else (b-a)/p)if edge else a
                    position=x+p/2 if edge else x
                    entries.append((position,value**power))
                table[power,edge,'sum']=sum(v for _,v in entries)
                table[power,edge,'weighted']=sum(x*v for x,v in entries)
                table[power,edge,'absolute']=sum(abs(x)*v for x,v in entries)
        tables.append(table)
    quadratic=energy=absolute_fourth=F(0)
    for edge_axis in range(3):
        q=e=fourth=F(1)
        for j,table in enumerate(tables):
            edge=j==edge_axis
            q*=table[2,edge,'weighted'if j==axis else'sum']
            e*=table[2,edge,'sum']
            fourth*=table[4,edge,'absolute'if j==axis else'sum']
        quadratic+=p**3*q/2;energy+=p**3*e/2
        # Both source relative bond angles equal p² times the edge value
        # in these tables, including the P difference quotient.
        absolute_fourth+=p**7*fourth/24
    i0,i1=integrals();limit_energy=6*i0**3 if channel=='Q'else F(3,2)*i1*i0*i0
    limit=center[axis]*limit_energy
    assert quadratic==center[axis]*energy,'exact source translation/work identity failed'
    return {'channel':channel,'divisions':divisions,'axis':axis,'center':list(center),
            'quadratic_polarization_per_kappa':quadratic,
            'quadratic_energy_per_kappa':energy,
            'exact_source_polarization_lower_per_kappa':quadratic-absolute_fourth,
            'exact_source_polarization_upper_per_kappa':quadratic+absolute_fourth,
            'cosine_remainder_per_kappa':absolute_fourth,
            'continuum_polarization_per_kappa':limit,
            'continuum_energy_per_kappa':limit_energy,
            'exact_affine_translation_identity':True,
            'actual_kappa_numerical_value_assigned':False}
