from fractions import Fraction as F
import numpy as np
import sympy as s
from .ring import construction,reduction,finite_volume_bound,memory_kernel,fine_records,evolve

def test_complete_parent_spectrum_and_symmetry_across_lapack_drivers():
 from .ring import hadron,hadron_eigensystem
 d=hadron();labels=d['hadron_labels'];energies=[]
 for driver in('evr','evd','ev','evx'):
  rows,e,U=hadron_eigensystem(d,1,1,driver=driver)
  H=d['Hh'].toarray()[np.ix_(rows,rows)]
  assert len(e)==16 and np.linalg.norm(H@U-U*e)<2e-12
  assert np.max(abs(e-np.linalg.eigvalsh(H)))<2e-12
  assert abs(e[0]-(d['mass'][1]-1+(181-np.sqrt(5473))/48))<2e-12
  assert e[1]-e[0]>1/50
  columns=[]
  for xa in range(2):
   v=np.zeros(16)
   for xq,aq in((0,xa),(1,1-xa)):
    for sq,sa in((0,1),(1,0)):
     v[next(j for j,i in enumerate(rows)if labels[i]==(1,xq,sq,1,aq,sa))]=(1 if sq==0 else -1)/2
   columns.append(v)
  B=np.array(columns).T
  assert np.linalg.norm(U[:,0]-B@(B.T@U[:,0]))<2e-12
  energies.append(e[0])
 assert max(energies)-min(energies)<2e-12

def test_native_two_site_source_and_record_transfer():
 from model.openwave_meson_weak_current.decay import construction as previous,complete_records
 a=construction(2);b=previous()
 assert np.array_equal(a['keep'],b['keep'])
 for key in('H','H0','V'):assert np.max(abs((a[key]-b[key]).toarray()))<3e-11
 rows,_=evolve(a,(720,));r=fine_records(a,rows[0]['position_state']);old=complete_records(b,720)
 assert np.max(abs(np.array(r['conditional_position_probabilities'])*a['success']/3-np.array(old['complete_position_colour_triplet_probability_each'])))<1e-11
 assert np.max(abs(np.array(r['conditional_momentum_probabilities'])*a['success']/3-np.array(old['complete_momentum_colour_triplet_probability_each'])))<1e-11
 assert abs(r['unconditional_position_readout_work']-old['total_position_readout_mean_work'])<1e-9
 assert abs(r['unconditional_momentum_readout_work']-old['total_momentum_readout_mean_work'])<1e-9
 assert abs(sum(r['conditional_position_probabilities'])-1)<1e-12
 assert max(abs(x-1)for x in r['both_complete_original_probability_sums'])<1e-12

def test_exact_native_spin_injection_identity():
 pauli=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1));bars=(s.eye(2),*(-p for p in pauli));eta=(1,-1,-1,-1)
 T=s.zeros(8,2)
 for t in range(2):
  for b in range(2):
   for a in range(2):
    for q in range(2):T[4*t+2*b+a,q]=sum(eta[k]*bars[k][t,q]*bars[k][b,a]for k in range(4))
 assert T.H*T==8*s.eye(2)
 for L in(2,4):
  d=construction(L);B,h,error=reduction(d);m=memory_kernel(d,B,h)
  assert abs(m['initial_quadratic_pair_probability_coefficient']-3/512)<1e-15
  assert m['Schur_resolvent_numeric_residual']<1e-10

def test_rational_finite_volume_tail_and_boundaries():
 r=finite_volume_bound(32,F(1,100));bound=F(r['prepared_domain_vector_error_upper_exact'])
 assert 0<bound<F(2546754,10**18)
 assert r['first_omitted_hop_order']==15
 assert not r['floating_evolution_roundoff_included']
 assert finite_volume_bound(1024,F(1))['first_omitted_hop_order']==511
