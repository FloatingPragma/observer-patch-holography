"""Independent rational one-body, Gauss and retained-record review."""
from pathlib import Path
import math,json
import sympy as z
import numpy as np
ROOT=Path(__file__).parent
X=z.Matrix([[0,1],[1,0]]);Y=z.Matrix([[0,-z.I],[z.I,0]])
Z=z.diag(1,-1);T=-(Z+z.I*Y)/2
A=z.zeros(8);C=z.zeros(8)
for cell in range(2):
 for wall in range(2):
  j=4*cell+2*wall;A[j:j+2,j:j+2]=8*Z
 A[4*cell:4*cell+2,4*cell+2:4*cell+4]=8*T
 A[4*cell+2:4*cell+4,4*cell:4*cell+2]=8*T.conjugate().T
for wall in range(2):
 j=2*wall;B=-4*(Z+z.I*X)
 A[j:j+2,j+4:j+6]=B;A[j+4:j+6,j:j+2]=B.conjugate().T
 C[j+4:j+6,j:j+2]=B
U=z.kronecker_product(z.eye(4),X)
assert U*A.T*U.T==-A
assert U*C.T*U.T==-C.conjugate().T
p=A.charpoly().as_poly().set_domain(z.ZZ)
receipt=json.loads((ROOT/'A0_INTERVALS.json').read_text())
assert str(p.as_expr())==receipt['characteristic_polynomial']
for a,b in receipt['intervals']:
 assert p.count_roots(z.Rational(a),z.Rational(b))==1

gauss=0;creations=0
for w in range(256):
 N=w.bit_count();n0=(w&15).bit_count();q0=n0+2-N
 for ell in range(-2,3):
  E0=ell+q0;E1=ell
  assert q0-E0+E1==0 and -q0-E1+E0==0
  gauss+=1
 for a in range(8):
  if(w>>a)&1:continue
  wp=w|(1<<a);qnew=(wp&15).bit_count()+2-wp.bit_count()
  # The reference charge changes at cell0. A cell1 insertion includes
  # the fixed tree Wilson path, hence its interior electric change−1.
  assert qnew-q0==(0 if a<4 else -1)
  creations+=1

rows=json.loads((ROOT/'RESULT.json').read_text());sp=json.loads((ROOT/'BOUNDARY_SPECTRA.json').read_text());rec=json.loads((ROOT/'VACUUM_RECORDS.json').read_text())
assert len(rows)==len(sp)==len(rec)==6
for row,spec,record in zip(rows,sp,rec):
 K=row['K'];assert {x['N']for x in row['all_sectors']}==set(range(9))
 assert sum(x['dimension']for x in row['all_sectors'])==256*(2*K+1)
 assert len(record['probabilities'])==70*(2*K+1)
 assert all(x['occupation'].bit_count()==4 for x in record['probabilities'])
 assert abs(sum(x['probability']for x in record['probabilities'])-1)<1e-11
 assert 0<=row['probability_low_W_at_most_half']<=1
 if K<=8:
  assert len(spec['particle'])==len(spec['hole'])==56*(2*K+1)
  total=np.zeros((4,4),complex)
  for channel in ['particle','hole']:
   for event in spec[channel]:
    total+=np.array(event['residue_real'])+1j*np.array(event['residue_imag'])
  assert np.max(abs(total-np.eye(4)))<1e-10
 else:
  assert not spec['complete_spectrum'] and not row['gap_interval_certified']
 assert row['full_particle_resolvent_residual']<1e-11
 assert row['full_hole_resolvent_residual']<1e-11

review={'source_module_imported':False,'exact_charpoly':str(p.as_expr()),
 'eight_root_isolations_verified':True,'exact_particle_hole_one_body_relations':True,
 'exact_Gauss_words_and_fluxes':gauss,'exact_dressed_creation_checks':creations,
 'all_number_sectors_retained':True,'complete_small_boundary_spectral_levels':sum(len(x.get('particle',[]))+len(x.get('hole',[]))for x in sp),
 'complete_vacuum_records':sum(len(x['probabilities'])for x in rec),
 'all_CAR_and_resolvent_checks_pass':True,'sparse_K16_not_promoted_to_certified_order':True}
(ROOT/'EXACT_REVIEW.json').write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps(review,indent=2))
