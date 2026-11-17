"""Rational exterior bounds and independent complete-record checks."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json,math
ROOT=Path(__file__).parent

def atan_bounds(q,n=40):
 a=sum(((-1)**j*F(1,(2*j+1)*q**(2*j+1))for j in range(n)),F(0))
 assert n%2==0
 return a,a+F(1,(2*n+1)*q**(2*n+1))

@lru_cache(None)
def pi_lower():
 a,b=atan_bounds(5);c,d=atan_bounds(239)
 return 16*a-4*d

@lru_cache(None)
def cosine_upper(S,j):
 x=pi_lower()*F(j,2*S+1)
 return sum(((-1)**k*x**(2*k)/math.factorial(2*k)for k in range(11)),F(0))

def sqrt_upper(x,digits=30):
 x=F(x);scale=10**digits
 n=math.isqrt(x.numerator*scale*scale//x.denominator)
 return F(n+1,scale)

def exact_exterior(S,N,K,m):
 n=min(N,2*S-N)
 amin=-16*sum((cosine_upper(S,j)for j in range(1,n+1)),F(0))-abs(F.from_float(m))
 c=sqrt_upper(64*((n+1)**2//4))
 return amin,c,F((K+1)**2,4)+amin-2*c

r=json.loads((ROOT/'RESULT.json').read_text());spectra=json.loads((ROOT/'ALL_BOUNDARY_SPECTRA.json').read_text())
rows=[]
for case,sp in zip(r,spectra):
 assert (case['S'],case['mass_label'],case['K'])==(sp['S'],sp['mass'],sp['K'])
 S=case['S'];K=case['K'];m=case['mass']
 assert {x['N']for x in case['all_number_sectors']}==set(range(2*S+1))
 for row in case['all_number_sectors']:
  N=row['N'];amin,c,d=exact_exterior(S,N,K,m)
  b=row['infinite_tail']
  assert F.from_float(b['coupling_norm_upper'])>=c
  assert F.from_float(b['discarded_block_lower'])<=d
  assert row['dimension']==math.comb(2*S,N)*(2*K+1)
 assert len(sp['particle'])==math.comb(2*S,S+1)*(2*K+1)
 assert len(sp['hole'])==math.comb(2*S,S-1)*(2*K+1)
 total=[[0.,0.],[0.,0.]]
 for channel in ['particle','hole']:
  for event in sp[channel]:
   assert event['gap']>0
   for a in range(2):
    for b in range(2):total[a][b]+=event['residue_real'][a][b]
 assert max(abs(total[a][b]-(a==b))for a in range(2)for b in range(2))<1e-10
 rows.append({'S':S,'mass':case['mass_label'],'K':K,'rational_exterior_bounds_pass':True,
  'all_number_sectors_and_particle_hole_levels_retained':True,'CAR_sum_pass':True,
  'source_module_imported':False})
save={'cases':rows,'case_count':len(rows),'complete_particle_hole_levels':sum(len(x['particle'])+len(x['hole'])for x in spectra),
 'finite_eigensystem_bounds':'conditional on the stated standard-IEEE arithmetic model; no external interval BLAS',
 'physical_lepton_pole_claim':False}
(ROOT/'EXACT_REVIEW.json').write_text(json.dumps(save,indent=2)+'\n')
print(json.dumps({'cases':len(rows),'levels':save['complete_particle_hole_levels']},indent=2))
