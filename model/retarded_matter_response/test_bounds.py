from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import solve

from model.retarded_matter_response import bounds as b, response as r


def test_exact_frequency_error_domains_and_geometric_decay():
    assert b.memory_taylor_remainder(0,0)==0
    for p in range(5):
        assert b.memory_taylor_remainder(r.GAP/4,p+1)==b.memory_taylor_remainder(r.GAP/4,p)/4
    assert b.large_frequency_remainder(8,4)==F(1,6144)
    assert b.certificate()['gap_extended_to_excited_reference_or_infinite_volume'] is False


def test_computed_exact_moments_enclose_direct_complex_resolvent():
    data=r.dilation();a,c=data.vectors
    v=(a+np.sqrt(3)*c)/36
    k=data.h16.toarray()/16+np.diag(data.sqrt3_shift*np.sqrt(3)/16)
    z=8j
    approximation=np.zeros((8,8),dtype=complex)
    for n in range(5):
        encoded=data.moment(n)
        matrix=np.array([[float(F(x))+np.sqrt(3)*float(F(y)) for x,y in row] for row in encoded])
        approximation+=matrix/z**(n+1)
    direct=v.T@solve(z*np.eye(560)-k,v)
    assert np.linalg.norm(direct-approximation,2)<float(b.large_frequency_remainder(8,4))


@pytest.mark.parametrize('call',[
    lambda:b.memory_taylor_remainder(r.GAP,2),lambda:b.memory_taylor_remainder(2*r.GAP,2),
    lambda:b.memory_taylor_remainder(-1,2),lambda:b.memory_taylor_remainder(0.0,2),
    lambda:b.memory_taylor_remainder(0,True),lambda:b.memory_taylor_remainder(0,-1),
    lambda:b.large_frequency_remainder(2,4),lambda:b.large_frequency_remainder(1,4),
])
def test_frequency_boundaries_reject_unproved_domains(call):
    with pytest.raises(ValueError):call()
