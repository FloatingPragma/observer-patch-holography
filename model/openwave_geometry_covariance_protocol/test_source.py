from fractions import Fraction as F
from math import comb
import pytest
from . import source as m


def test_complete_pair_records_and_exact_fourth_moments():
    for n in (1,2,3,5):
        r=m.pair_records(n);assert len(r['complete_count_pair_records'])==(n+1)**2
        assert sum(F(x['probability_exact'])for x in r['complete_count_pair_records'])==1
        assert F(r['variance_mean_exact'])==m.Q*m.Q/n
        assert F(r['variance_second_moment_exact'])==(m.Q*m.Q/n)**2*(3-F(1,n))


def test_every_complete_operator_and_memory_outcome():
    for n in (1,2,3):
        r=m.complete_instrument_control(n)
        assert r['complete_operator_work_identity_error_exact']=='0'
        assert len(r['all_count_probabilities'])==n+1


def test_nonasymptotic_sampling_probability_and_failure_controls():
    r=m.sampling_plan();k=r['independent_blocks']
    exact=sum((F(comb(k,j))*F(1,4)**j*F(3,4)**(k-j)for j in range((k+1)//2,k+1)),F(0))
    assert F(r['median_failure_upper_exact'])==exact< F(1,1000)
    for e,k in ((0,3),(1,3),(F(1,10),4)):
        with pytest.raises(ValueError):m.sampling_plan(e,k)


def test_common_unit_interval_includes_all_endpoints():
    c=(F(1,100),F(11,1000));s=(F(1),F(11,10));p=(F(1),F(6,5));u=(F(1,2),F(3,5));v=(F(9,10),F(1))
    lo,hi=m.propagate_strength(c,s,p,u,v)
    for cc in c:
        for ss in s:
            for pp in p:
                for uu in u:
                    for vv in v:
                        xlo,xhi=m.propagate_strength((cc,cc),(ss,ss),(pp,pp),(uu,uu),(vv,vv))
                        assert lo<=xlo<=xhi<=hi
    with pytest.raises(ValueError):m.propagate_strength((0,1),s,p,u,v)


def test_count_freedom_and_resolution_are_exposed():
    for n in (1,2,4,8):
        lo,hi=m.normalized_strength_interval(n);a,b=m.normalized_strength_interval(2*n)
        assert 0<a<=b<lo<=hi
        assert F(m.quantization_bias(n,0)['covariance_bias_upper_exact'])==0
        assert F(m.quantization_bias(n,F(1,10000))['covariance_bias_upper_exact'])>0
    for n in (0,-1,True,F(1)):
        with pytest.raises(ValueError):m.normalized_strength_interval(n)
