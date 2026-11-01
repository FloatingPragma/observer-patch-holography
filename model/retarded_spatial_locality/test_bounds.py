from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import bounds as b
from ..retarded_matter_response import response


def test_exact_geometric_moments_against_known_half_ratio_sums():
    expected=(1,2,6,26,150,1082)
    for n,total in enumerate(expected):
        assert b.geometric_power_sum(F(1,2),n)==total
        assert sum(F(d**n,2**d) for d in range(1,100))<total


def test_three_dimensional_shell_count_and_rational_tail_majorant():
    shells={d:0 for d in range(5)}
    for point in product(range(-4,5),repeat=3):
        shells[max(map(abs,point))]+=1
    assert shells[0]==1
    assert all(shells[d]==24*d*d+2 for d in range(1,5))
    rate=F(1,3);q=1/(1+rate)
    for n in range(4):
        partial=sum(F(count)*d**n*q**d for d,count in shells.items())
        assert partial<b.shell_moment_upper(rate,n)


def test_every_source_constant_retains_the_neumann_margin():
    for eta in (F(1,100),F(1,4),F(1),F(8)):
        value=b.constants(eta,max(eta,F(2)))
        assert 0<value['inverse_decay']<=value['G_half_decay']<=F(1,8)
        assert value['inverse_decay']*value['G_prefactor']*value['weighted_sum']<=1/(2*value['inverse_norm'])
        assert value['memory_operator_norm']==1/(4*eta)


def test_real_axis_can_destroy_invertibility_even_for_an_exact_atomic_response():
    # Atomic Hubbard response with occupation r: (1-r)/z+r/(z-U).
    # At the real z=U(1-r), G vanishes between its poles.
    r,z,u=F(1,4),F(3,4),F(1)
    assert (1-r)/z+r/(z-u)==0
    with pytest.raises(ValueError):b.constants(0,2)


def test_same_source_finite_resolvent_obeys_moment_and_memory_bounds():
    # Independent numerical diagnostic, not an admission criterion in packet.
    data=response.dilation('ground','down')
    v=(data.vectors[0]+np.sqrt(3)*data.vectors[1])/36
    h=data.h16.toarray()/16+np.diag(data.sqrt3_shift*np.sqrt(3)/16)
    m1=response.first_numerator('ground','down')/48
    for eta in (F(1,2),F(1)):
        z=0.25+1j*float(eta)
        value=b.constants(eta,2)
        g=v.T@np.linalg.solve(z*np.eye(len(h))-h,v)
        inverse=np.linalg.inv(g)
        memory=z*np.eye(8)-m1-inverse
        assert np.linalg.eigvalsh(-(g-g.conj().T)/(2j)).min()>1/float(value['inverse_norm'])
        assert np.linalg.norm(inverse,2)<float(value['inverse_norm'])
        assert np.linalg.norm(memory,2)<=float(value['memory_operator_norm'])*(1+1e-12)


def test_taylor_remainders_vanish_with_increment_and_keep_frequency_boundary():
    assert b.spatial_taylor_remainder(1,2,0,1)==0
    x=b.spatial_taylor_remainder(1,2,F(1,10**25),1)
    y=b.spatial_taylor_remainder(1,2,F(1,2*10**25),1)
    assert y==x/4
    assert b.frequency_taylor_remainder(1,2,0,2)==0
    assert b.frequency_taylor_remainder(1,2,F(1,8),2)<b.frequency_taylor_remainder(1,2,F(1,4),2)
    with pytest.raises(ValueError):b.frequency_taylor_remainder(1,2,F(1,2),2)


@pytest.mark.parametrize('call',[
    lambda:b.constants(1.0,2),lambda:b.constants(True,2),
    lambda:b.constants(1,F(1,2)),lambda:b.constants(1,2,0),
    lambda:b.constants(1,2,F(1,8),-1),lambda:b.geometric_power_sum(1,3),
    lambda:b.geometric_power_sum(F(1,2),-1),lambda:b.shell_moment_upper(0,1),
    lambda:b.spatial_taylor_remainder(1,2,-1,2),
    lambda:b.frequency_taylor_remainder(1,2.0,F(1,8),2),
])
def test_inexact_or_out_of_domain_inputs_fail_closed(call):
    with pytest.raises(ValueError):call()
