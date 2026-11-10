"""Exact algebraic evaluation of the complete charged pair denominator."""
from functools import lru_cache
import sympy as s
from ..ground_composite_response.model import source
from ..atomic_reference_control.ground_exchange import actual_spin_coefficients
from ..family_cell.quotient import evaluate,outward
from .model import moment_data


@lru_cache(maxsize=1)
def field_data():
    e=s.Symbol('e')
    f=s.Poly((e*(e-1))**3-313*(e*(e-1))**2+31708*e*(e-1)-1041468,e)
    assert f.is_irreducible
    root=s.CRootOf(f.as_expr(),0)
    return e,f,s.QQ.algebraic_field(root)


def field_value(polynomial):
    _,_,field=field_data();value=field.zero
    for coefficient in polynomial.all_coeffs():value=value*field.unit+field.convert(coefficient)
    return value


def polynomial(value):
    e,_,field=field_data()
    coefficients=[s.Rational(c.numerator,c.denominator) for c in value.to_list()]
    return s.Poly.from_list(coefficients,e)


def encoded(value):
    return [str(c) for c in polynomial(value).all_coeffs()]


@lru_cache(maxsize=None)
def response(family):
    e,f,field=field_data();p,q,trace_rows,blocks=moment_data(family)
    x=p.gen
    coefficients=[]
    for i in range(p.degree()):
        c=s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1,p.degree()+1)),x)
        coefficients.append(field_value(c)/field_value(p.diff()))
    # The normalized selected source projector has rank one in this magnetic block.
    h=source(1)['H'];power=s.eye(h.rows);trace=field.zero
    for coefficient in coefficients:
        trace+=coefficient*field.convert(s.trace(power));power=power*h
    assert trace==field.one
    moments=tuple(sum((c*field.convert(t) for c,t in zip(coefficients,row,strict=True)),field.zero)
                  for row in trace_rows)
    actual=actual_spin_coefficients()
    beta=field_value(s.Poly(actual['beta_ec_numerator'],e))/field_value(s.Poly(actual['beta_ec_denominator'],e))
    assert moments[0]==(beta if family=='ec' else (field.one-beta)/2)
    z=s.Symbol('z')
    denominator=s.Poly.from_dict({(j,):field.convert(q.nth(j)) for j in range(q.degree()+1)},z,domain=field)
    entries={}
    for i,moment in enumerate(moments):
        for j in range(i+1,q.degree()+1):
            power=j-i-1
            entries[(power,)]=entries.get((power,),field.zero)+field.convert(q.nth(j))*moment
    numerator=s.Poly.from_dict(entries,z,domain=field)
    common=s.gcd(numerator,denominator)
    reduced_n=numerator.exquo(common);reduced_d=denominator.exquo(common)
    return {'p':p,'q':q,'moments':moments,'numerator':numerator,'denominator':denominator,
            'reduced_numerator':reduced_n,'reduced_denominator':reduced_d,'blocks':blocks}


@lru_cache(maxsize=None)
def pair_integral(family):
    """Apply the signed response functional to its reflected rational resolvent."""
    data=response(family);e,f,field=field_data()
    q=data['denominator'];n=data['numerator']
    result,proof=pair_functional(q,n,data['moments'],field.one-2*field.unit)
    return result,{'response':data,**proof}


def pair_functional(q,n,moments,offset):
    """Evaluate sum_ij w_i w_j/(offset+epsilon_i+epsilon_j) exactly."""
    if q.gen!=n.gen or q.domain!=n.domain or len(moments)!=q.degree():
        raise ValueError('incompatible rational response and complete moments')
    if n.degree()>=q.degree():raise ValueError('a proper response resolvent is required')
    field=q.domain;z=q.gen
    shift=s.Poly.from_dict({(1,):-field.one,(0,):-offset},z,domain=field)
    reflected_q=q.compose(shift);reflected_n=n.compose(shift)
    inverse=s.invert(reflected_q,q)
    assert (reflected_q*inverse).rem(q)==s.Poly(1,z,domain=field)
    evaluated=(-reflected_n*inverse).rem(q)
    assert (reflected_q*evaluated+reflected_n).rem(q).is_zero
    result=sum((evaluated.rep.nth(i)*m for i,m in enumerate(moments)),field.zero)
    return result,{'reflected_denominator':reflected_q,
                   'inverse':inverse,'evaluated_resolvent':evaluated}


@lru_cache(maxsize=1)
def exchange():
    e,f,field=field_data();parts={family:pair_integral(family)[0] for family in ('ec','L')}
    total=(parts['ec']+2*parts['L'])/4
    root=tuple(map(s.Rational,actual_spin_coefficients()['ground_energy_isolating_interval']))
    interval=evaluate(polynomial(total),root)
    assert interval[0]>0
    return {'exact_J_coefficients_descending_in_E0':encoded(total),
            'exact_J_interval':outward(interval,48),
            'ec_integral_coefficients':encoded(parts['ec']),
            'L_integral_coefficients':encoded(parts['L']),
            'native_internal_multiplicities':{'ec':1,'L':2},
            'actual_four_hopping_coefficient_squared':'1/4',
            'positive_pair_denominator_keeps_all_charged_states':True}
