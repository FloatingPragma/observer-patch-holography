from fractions import Fraction as F
from itertools import combinations
import math
import numpy as np
import pytest
from scipy.linalg import expm
from . import collar as c
from . import regional as r

@pytest.mark.parametrize('x,radius',[(F(0),0),(F(0),3),(F(1,3),0),(F(1,3),5),(F(3),2),(F(12),30)])
def test_rational_tail_bound(x,radius):
 value=c.tail_upper(x,radius)
 # Independent direct positive-series sum, avoiding subtractive cancellation.
 truth=sum(float(x)**n/math.factorial(n)for n in range(radius,150))
 assert float(value)>=truth*(1-1e-14)
 assert float(value)<=truth*(1+1e-10)+1e-14

@pytest.mark.parametrize('p',(F(1,2),F(1,3)))
def test_noncommuting_finite_boundary_flow(p):
 V=8;edges=tuple((i,i+1)for i in range(V-1))+((0,V-1),);A=tuple(range(5));q=4
 basis,G=r.laplacian(V,edges,q);_,GA=r.laplacian(V,tuple(e for e in edges if set(e)<=set(A)),q)
 N=np.diag([((x>>0)&1)+((x>>4)&1)for x in basis]);D=GA+(1-2*p)*N;J=1/(1-p*(1-p));O=np.diag([(x>>2)&1 for x in basis]);time=F(1,100)
 def flow(H):
  U=expm(1j*float(time*J)*np.array(H,float));return U@O@U.conj().T
 error=np.linalg.norm(flow(D)-flow(G),2)
 budget=c.locality_bound(1,2,2,2,time,p)
 assert 0<error<=float(budget)

def test_uniform_region_bound_dominates_every_local_sector():
 for level in(1,2,6):
  V=240*4**level;M=630*4**level
  bound=c.uniform_modular_error(V,M,V//2,12,30,3,F(1,2))
  for k in range(13):assert bound>=r.modulus_bound(V,M,V//2,12,30,3,k,F(1,2))['generator_error']

def test_collar_tail_beats_exponential_boundary_count():
 rows=[c.locality_bound(1,6,6**radius,radius,F(1,100),F(1,2))for radius in(10,20,40,80)]
 assert all(a>b for a,b in zip(rows,rows[1:]))and rows[-1]<F(1,10**50)

def test_invalid_inputs_are_not_silent_bounds():
 with pytest.raises(ValueError):c.tail_upper(F(-1),3)
 with pytest.raises(ValueError):c.tail_upper(F(1),True)
 with pytest.raises(ValueError):c.uniform_modular_error(20,10,10,12,8,2,F(1,2))
