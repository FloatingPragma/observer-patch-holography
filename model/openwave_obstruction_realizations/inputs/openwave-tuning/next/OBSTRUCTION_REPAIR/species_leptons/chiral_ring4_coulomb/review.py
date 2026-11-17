"""Independent polynomial identity and full-Slater occupation review."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,math
import numpy as np
ROOT=Path(__file__).parent
poly=json.loads((ROOT/'OCCUPATION_POLYNOMIALS.json').read_text())
words=np.arange(65536,dtype=int)
def zeta(coeff,scale):
 a=np.zeros(65536,dtype=np.int64)
 for mask,c in coeff.items():
  v=F(c)*scale;assert v.denominator==1;a[int(mask)]=int(v)
 for bit in range(16):
  target=words[(words&(1<<bit))!=0];a[target]+=a[target^(1<<bit)]
 return a
b=zeta(poly['B'],64);b2=zeta(poly['B2'],4096)
assert np.array_equal(b*b,b2)
for w in range(65536):
 n=w.bit_count();q=[(w&((1<<(4*(x+1)))-1)).bit_count()+2*(3-x)-n for x in range(4)]
 assert b[w]==sum((4*x-sum(q))**2 for x in q)

X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
A=np.zeros((16,16),complex);C=np.zeros_like(A)
for x in range(4):
 a=4*x;T=-4*(Z+1j*Y)
 A[a:a+2,a:a+2]=8*Z;A[a+2:a+4,a+2:a+4]=8*Z
 A[a:a+2,a+2:a+4]=T;A[a+2:a+4,a:a+2]=T.conj().T
 for wall in range(2):
  a=4*x+2*wall;d=4*((x+1)%4)+2*wall;T=-4*(Z+1j*X)
  if x<3:A[a:a+2,d:d+2]=T;A[d:d+2,a:a+2]=T.conj().T
  else:C[a:a+2,d:d+2]=T
controls=json.loads((ROOT/'ALL_POINTWISE_CONTROLS.json').read_text())
lookup={(r['grid_index'],r['N']):r for r in controls};checks=[]
for j in [0,1,128,256,511,512,768]:
 theta=2*np.pi*j/1024;H=A+C*np.exp(1j*theta)+C.conj().T*np.exp(-1j*theta)
 values,U=np.linalg.eigh(H)
 for N in [7,8]:
  row=lookup[j,N]
  if not row['nondegenerate_control']:continue
  occupations=np.array(list(combinations(range(16),N)),dtype=int)
  coefficients=np.linalg.det(U[occupations,:N])
  w=np.array([sum(1<<int(i)for i in group)for group in occupations])
  p=abs(coefficients)**2
  assert abs(sum(p)-1)<1e-11
  mean=float(p@(b[w]/64));second=float(p@(b2[w]/4096))
  assert abs(mean-row['B_expectation'])<1e-10
  assert abs(second-row['B_squared_expectation'])<1e-10
  checks.append({'grid_index':j,'N':N,'all_Slater_occupations':len(p),'mean_error':abs(mean-row['B_expectation']),
   'second_moment_error':abs(second-row['B_squared_expectation'])})
result=json.loads((ROOT/'RESULT.json').read_text())
prepared=json.loads((ROOT/'inputs/prepared_results.json').read_text())
assert len(controls)==2048
assert all(r['improved_pointwise_lower']>=r['free_ground_lower']for r in controls)
for case,out in zip(prepared,result['prepared_response_consumers']):
 assert out['neutral_energy_excess_upper']>=0
 assert out['prepared_Green_error_upper']<=case['infinite_prepared_Green_error_upper']
 assert not out['exact_neutral_fidelity_certified'] and not out['true_vacuum_propagator_certified']
receipt={'source_module_imported':False,'exact_polynomial_words':65536,
 'independent_full_Slater_checks':checks,'pointwise_controls_retained':2048,
 'all17_number_sectors_retained':len(result['every_number_sector'])==17,
 'neutral_fidelity_not_promoted':True}
(ROOT/'REVIEW.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'full_Slater_controls':len(checks),'words':65536},indent=2))
