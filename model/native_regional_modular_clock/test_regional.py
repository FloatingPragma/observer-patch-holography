from fractions import Fraction as F
from itertools import combinations
import numpy as np
import pytest
from . import regional as r

@pytest.mark.parametrize('V,edges,q',[
 (4,((0,1),(1,2),(2,3)),2),
 (5,((0,1),(0,4),(1,2),(2,3),(3,4)),1),
 (6,tuple(combinations(range(6),2)),3),
 (7,((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(0,6)),4),
])
def test_all_supported_and_zero_blocks_are_literal_partial_trace(V,edges,q):
 basis,rho,Z=r.reference(V,edges,q)
 for A in tuple(a for a in ((0,1),(0,1,2)) if len(a)<=V-2):
  direct=r.partial_trace(V,basis,rho,A)
  for k in range(len(A)+1):
   local,block,_=r.regional_block(V,edges,q,A,k)
   assert local==direct[k][0] and np.array_equal(block,direct[k][1])
  assert sum(sum(block[i,i]for i in range(len(block)))for _,block in direct.values())==1

@pytest.mark.parametrize('q',range(7))
def test_exact_swap_moment_identity(q):
 V=6;edges=((0,1),(0,2),(0,3),(1,2),(2,3),(3,4),(4,5))
 basis,G=r.laplacian(V,edges,q);mom=r.canonical_variance(V,edges,q);d=len(basis)
 assert F(sum(G[i,i]for i in range(d)),d)==mom['mean']
 GG=G@G;assert F(sum(GG[i,i]for i in range(d)),d)==mom['second']
 _,rho,_=r.reference(V,edges,q)
 mean=sum((rho@G)[i,i]for i in range(d));assert mean==mom['rho_mean']
 var=sum((rho@(G-mean*np.eye(d,dtype=object))@(G-mean*np.eye(d,dtype=object)))[i,i]for i in range(d))
 assert 0<=var<=mom['rho_variance_upper']

def test_fixed_region_overlap_does_not_disappear():
 assert r.overlap_witness()['half_density_rescaled_modular_overlap_commutator_norm']==F(4,3)

def test_local_modular_bound_converges_without_full_carrier():
 errors=[]
 for ell in range(1,8):
  V=240*4**ell;M=630*4**ell
  errors.append(max(r.modulus_bound(V,M,V//2,12,30,3,k,F(1,2))['generator_error']for k in range(13)))
 assert all(a>b for a,b in zip(errors,errors[1:])) and errors[-1]<F(1,1000)

def test_invalid_gate_and_boolean_parameters_fail():
 with pytest.raises(ValueError):r.modulus_bound(12,30,2,10,30,0,1,F(1,2))
 with pytest.raises(ValueError):r.modulus_bound(True,30,2,2,1,8,1,F(1,2))

def test_exact_commutator_count():
 V=6;edges=((0,1),(0,2),(0,3),(1,2),(2,3),(3,4),(4,5));A=(0,1);q=3
 basis,rho,_=r.reference(V,edges,q);N=np.diag([sum((x>>i)&1 for i in A)for x in basis]);C=rho@N-N@rho
 assert sum(x*x for x in C.flat)==r.cross_sector_hs_squared(V,len(edges),q,len(r.split(V,edges,A)[2]))

def test_numeric_logarithm_within_exact_rational_bound():
 # A finite control; the theorem uses the exact convergent logarithm series.
 V=9;edges=tuple((i,i+1)for i in range(8))+((0,8),);A=(0,1);q=4
 for k in range(3):
  _,R,data=r.regional_block(V,edges,q,A,k);D=np.array(data['D'],float)
  values,basis=np.linalg.eigh(D);a=float(data['a']);M=len(edges)
  H=(basis*(-2*M*np.log1p(-values/(2*M*a))))@basis.T
  # Compare the logarithmic term to its exact first order term.
  bound=r.modulus_bound(V,M,q,2,1,2,k,F(1,2))
  assert np.linalg.norm(H-D/a,2)<=float(bound['logarithm_remainder'])+1e-12


@pytest.mark.parametrize('bad',(True,1.5,'1'))
def test_region_sector_types_fail_closed(bad):
 with pytest.raises(ValueError):r.regional_block(5,((0,1),(1,2),(2,3),(3,4)),bad,(0,1),1)
 with pytest.raises(ValueError):r.regional_block(5,((0,1),(1,2),(2,3),(3,4)),2,(0,1),bad)

def test_matrix_allocation_and_incomplete_source_fail_closed():
 with pytest.raises(ValueError):r.masks(64,2)
 with pytest.raises(ValueError):r.masks(20,10)
 basis,rho,_=r.reference(4,((0,1),(1,2),(2,3)),2)
 with pytest.raises(ValueError):r.partial_trace(4,basis[:-1],rho[:-1,:-1],(0,1))
 with pytest.raises(ValueError):r.partial_trace(4,basis,rho.astype(float),(0,1))


def test_local_binomial_ratio_uses_only_the_retained_region():
 for V in (4,9,12):
  for q in range(V+1):
   for v in range(min(4,V)+1):
    for k in range(v+1):
     assert r.local_tracial_eigenvalue(V,q,v,k)==F(r.choose(V-v,q-k),r.choose(V,q))
 V=2**400
 rows=[r.local_tracial_eigenvalue(V,V//2,12,k)for k in range(13)]
 assert sum(F(r.choose(12,k))*rows[k]for k in range(13))==1
 assert all(x>0 for x in rows)


@pytest.mark.parametrize('q',range(6))
def test_positive_reference_sampler_keeps_all_rejection(q):
 V=5;edges=((0,1),(0,4),(1,2),(2,3),(3,4))
 row=r.coherent_reference_mixture(V,edges,q)
 assert row['accepted_probability']+row['rejected_probability']==1
 assert row['accepted_probability']>=F(1,2)
 assert row['coherent_mixer_is_supplied']
