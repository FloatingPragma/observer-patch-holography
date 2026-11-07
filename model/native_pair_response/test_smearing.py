from fractions import Fraction as F
from itertools import combinations
from math import comb
import pytest
from . import smearing as w,response as r


@pytest.mark.parametrize('V',(4,5,6,8))
def test_primitive_subtraction_on_every_configuration_and_number(V):
    f=w.localized_profile(V,{0:2,1:-1});g=w.localized_profile(V,{0:1,2:-2})
    h=w.projection(f,g);assert not r.s.lower(h,V,2)
    for q in range(2,V-1):
        chi=F(q*(V-q),V*(V-1));a=F(V-2*q,V-2)
        states=tuple(combinations(range(V),q));ys=[];corrections=[]
        for x in states:
            actual=w.field(f,g,q,x)
            included=sum(v for pair,v in h.items() if set(pair)<=set(x))/(V*chi)
            assert actual==included;ys.append(actual)
            correction=a*sum(f[i]*g[i]*(int(i in x)-F(q,V)) for i in range(V))/(V*chi)
            corrections.append(correction)
        assert sum(ys)==0
        for i in range(V):assert sum(y for x,y in zip(states,ys) if i in x)==0
        assert sum(z*z for z in corrections)/len(states)==w.correction_bound(f,g,q)['correction_mean_square']


def test_finite_projection_needs_one_body_correction_off_half_filling():
    f=w.localized_profile(6,{0:2,1:-1});g=w.localized_profile(6,{0:1,2:-2})
    assert w.correction_bound(f,g,2)['correction_mean_square']>0
    assert w.correction_bound(f,g,3)['correction_mean_square']==0
    for V in (12,24,48,96,192):
        f=tuple(F(1 if i<V//2 else -1) for i in range(V))
        g=tuple(F(1 if i%3==0 else -1,1 if i%3==0 else 2) for i in range(V))
        error=w.correction_bound(f,g,V//3)
        assert error['correction_mean_square']<=error['correction_mean_square_upper']<=F(5,V)


@pytest.mark.parametrize('V',(4,6,8))
@pytest.mark.parametrize('m',(1,2,4))
def test_uniform_even_moment_bound_against_complete_canonical_enumeration(V,m):
    f=tuple(F(1 if i<V//2 else -1) for i in range(V))
    for q in range(1,V):
        chi=F(q*(V-q),V*(V-1))
        exact=sum(sum(f[i] for i in x)**(2*m) for x in combinations(range(V),q))/(comb(V,q)*(V*chi)**m)
        assert exact<=w.even_moment_bound(V,q,F(1),m)


def test_eighth_moment_constant_and_nonGaussian_finite_snapshot():
    assert w.even_moment_bound(240,80,F(1),4)<3888
    f=(3,-1,-1,-1);q=2;chi=F(1,3)
    fourth=sum(sum(f[i] for i in x)**4 for x in combinations(range(4),q))/(6*(4*chi)**2)
    assert fourth==9 and fourth!=3*(sum(x*x for x in f)/4)**2


@pytest.mark.parametrize('call',(
    lambda:w.projection((1,-1,True,-1),(1,-1,1,-1)),
    lambda:w.projection((1.,-1,1,-1),(1,-1,1,-1)),
    lambda:w.projection((1,-1,1,0),(1,-1,1,-1)),
    lambda:w.field((1,-1,1,-1),(1,-1,1,-1),2,(True,2)),
    lambda:w.field((1,-1,1,-1),(1,-1,1,-1),1,(0,)),
    lambda:w.even_moment_bound(6,2,1.0,4),
    lambda:w.even_moment_bound(6,0,F(1),4),
    lambda:w.even_moment_bound(6,2,F(1),True),
    lambda:w.localized_profile(6,{False:1}),
))
def test_quadratic_domains_reject_type_aliases_and_unsupported_profiles(call):
    with pytest.raises(ValueError):call()
