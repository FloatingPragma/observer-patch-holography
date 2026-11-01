from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest

from model.interacting_vacuum import ground
from model.retarded_matter_response import response, spectrum
from . import matrix as m


def test_full_source_columns_share_exact_polynomial_and_matrix_algebra():
    data=m.exact_data()
    assert m.annihilates_columns(data['polynomial'])
    damaged=list(data['polynomial'])
    a,b=damaged[-1];damaged[-1]=(a+1,b)
    assert not m.annihilates_columns(damaged)
    assert data['moments'][0]==((1296,0),(0,0))
    assert data['moments'][1][1]==(1296,0)
    assert data['identity_numerator']==spectrum.exact_data()['numerator']


def test_rank_polynomial_identity_requires_the_actual_edge_numerator():
    data=m.exact_data()
    square=m._convolution(data['identity_numerator'],data['identity_numerator'])
    damaged=list(data['adjacency_numerator'])
    a,b=damaged[-1];damaged[-1]=(a+1,b)
    other=m._convolution(damaged,damaged)
    residual=m._remainder([m._add(x,m._scale(y,-3)) for x,y in zip(square,other)],data['polynomial'])
    assert any(a or b for a,b in residual)


def test_matrix_rational_response_matches_independent_complete_source_resolvent():
    data=m.exact_data();source=response.dilation('ground','down')
    adjacency=ground.band_data()[0]
    v=(source.vectors[0]+np.sqrt(3)*source.vectors[1])/36
    k=source.h16.toarray()/16+np.diag(source.sqrt3_shift*np.sqrt(3)/16)
    with mp.workdps(120):
        for z in (complex(-.2,.3),complex(.7,.5)):
            x=16*mp.mpc(z.real,z.imag)-4*mp.sqrt(3)
            def evaluate(values):
                result=mp.mpc(0)
                for a,b in values:result=result*x+a+b*mp.sqrt(3)
                return result
            factor=mp.mpf(16)/(1296*evaluate(data['polynomial']))
            proposed=(complex(factor*evaluate(data['identity_numerator']))*np.eye(8)
                      +complex(factor*evaluate(data['adjacency_numerator']))*adjacency)
            actual=v.T@np.linalg.solve(z*np.eye(560)-k,v)
            assert np.linalg.norm(actual-proposed)<1e-11
            inverse=np.linalg.inv(actual)
            mask=(adjacency==0)&~np.eye(8,dtype=bool)
            assert np.linalg.norm(inverse[mask])<1e-11


def test_exact_static_inverse_recovers_both_identity_and_adjacency_coefficients():
    result=m.static_certificate()
    c0,c1,y0,y1=(tuple(map(F,result[key])) for key in ('c0','c1','inverse_c0','inverse_c1'))
    assert m._add(m._mul(c0,y0),m._scale(m._mul(c1,y1),3))==(1,0)
    assert m._add(m._mul(c0,y1),m._mul(c1,y0))==(0,0)
    assert m._field_interval(c0)[0]>12
    assert m._field_interval(c1)[1]<-4
    assert result['floating_linear_solver_used'] is False


def test_cyclic_bound_remainders_have_the_stated_domain():
    assert m.CYCLIC_GAP>response.GAP
    assert m.memory_remainder(F(1,100),8)<m.memory_remainder(F(1,100),4)
    assert m.memory_remainder(0,0)==0
    for radius in (m.CYCLIC_GAP,-1,.01,True):
        with pytest.raises(ValueError):m.memory_remainder(radius,4)
    for order in (-1,33,True):
        with pytest.raises(ValueError):m.memory_remainder(F(1,100),order)


def test_all_residue_bands_are_positive_rank_four_and_exhaust_observed_space():
    result=m.certificate()
    rows=result['matrix_residue_rows']
    assert len(rows)==78
    assert result['observed_cyclic_dimension']==312
    assert result['orthogonal_unobserved_dimension']==248
    assert result['band_line_counts']=={'plus':40,'minus':38}
    for band in (-1,1):
        selected=[row for row in rows if row['band']==band]
        weights=[tuple(map(F,row['local_weight_interval'])) for row in selected]
        assert 2*sum(x[0] for x in weights)<=1<=2*sum(x[1] for x in weights)
        assert all(row['rank']==4 for row in selected)
    assert F(rows[0]['gap_above_ground_interval'][0])>m.CYCLIC_GAP
    assert F(rows[-1]['gap_above_ground_interval'][1])<m.CYCLIC_UPPER
    assert result['gap_bound_applied_to_unobserved_sector'] is False
    assert result['common_Lorentz_metric_or_physical_masses_identified'] is False
