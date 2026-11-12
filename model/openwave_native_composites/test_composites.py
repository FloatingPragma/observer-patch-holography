import numpy as np
import scipy.sparse as sp
from .spectrum import sector,selected,native,U,J,T,MASS,norm

def test_native_casimirs_and_intertwiners():
 q,a,joins,_,_=native()
 for family in q+a:
  assert np.allclose(sum(g@g for g in family),4*np.eye(3)/3,atol=1e-13)
 for f in range(2):
  for k in range(8):assert np.max(abs(q[0][k]@joins[f]+joins[f]@a[f][k].T))<1e-13

def test_full_sector_embedding_and_complete_colour_projector():
 for kind in ('baryon','meson'):
  a=sector(kind);B=a['B'];h=a['h']
  assert norm(a['H']@B-B@sp.csr_matrix(h))<1e-12
  assert norm(B.conj().T@B-sp.eye(B.shape[1]))<1e-12
  assert np.linalg.eigvalsh(h)[0]>=MASS*(3 if kind=='baryon' else 2)

def test_binding_against_same_source_fragment_thresholds():
 for charge in (0,1):
  r=selected('baryon',charge,.5)
  assert r['ground_energy']<3*MASS+3*J/4+3*T
  assert 3*MASS+8*U/3+3*J/4-r['ground_energy']>8*U/3-3*T
 for charge in (-1,0,1):
  r=selected('meson',charge,0)
  assert r['ground_energy']<2*MASS+2*T
  assert 2*MASS+8*U/3+3*J/2-r['ground_energy']>8*U/3+3*J/2-2*T

def test_flavour_degeneracy_and_complete_preparation():
 b=[selected('baryon',q,.5) for q in (0,1)]
 assert b[0]['ground_energy']==b[1]['ground_energy']
 for row in b+[selected('meson',q,0) for q in (-1,0,1)]:
  p=row['spectral_preparation'];assert 0<p['success_probability']<1
  assert abs(p['success_probability']+p['failure_probability']-1)<1e-12
  assert abs(sum(row['successful_full_occupation_record_probabilities'])-1)<1e-9
  assert abs(sum(row['cell_charge_expectations'])-row['charge'])<1e-10
  assert row['eigen_residual']<1e-11
  complete=row['complete_unconditional_records']
  assert abs(sum(complete['preparation_success_and_every_full_sector_occupation'])+complete['preparation_failure']-1)<1e-9
  for conditional,unconditional in zip(row['successful_full_occupation_record_probabilities'],complete['preparation_success_and_every_full_sector_occupation']):
   assert abs(unconditional-p['success_probability']*conditional)<2e-12
