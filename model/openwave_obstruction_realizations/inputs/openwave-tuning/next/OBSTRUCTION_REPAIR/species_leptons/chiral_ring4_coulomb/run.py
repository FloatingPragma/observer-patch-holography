"""Keep the exact Coulomb variance in the complete-source scalar bound."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import importlib.util,json,math,time
import numpy as np
import mpmath as mp
ROOT=Path(__file__).parent;GRID=1024;AMPS=[0,.25,-.25,.5,-.5,1.,-1.,2.,-2.]
mp.iv.dps=40
spec=importlib.util.spec_from_file_location('frozen_source',ROOT/'inputs/source.py')
src=importlib.util.module_from_spec(spec);spec.loader.exec_module(src)

def add(poly,mask,c):poly[mask]=poly.get(mask,F(0))+c
B={}
for x in range(4):
 constant=3-2*x
 f=[F(int(i//4<=x))-F(4-i//4,4)for i in range(16)]
 add(B,0,F(constant**2,4))
 for i in range(16):add(B,1<<i,F(constant,2)*f[i]+f[i]**2/4)
 for i,j in combinations(range(16),2):add(B,(1<<i)|(1<<j),f[i]*f[j]/2)
B={m:c for m,c in B.items()if c}
B2={}
for a,c in B.items():
 for b,d in B.items():add(B2,a|b,c*d)
B2={m:c for m,c in B2.items()if c}
words=np.arange(1<<16,dtype=np.int64)
counts=np.array([[int(w&((1<<(4*(x+1)))-1)).bit_count()for x in range(4)]for w in words],dtype=int)
N=np.array([int(w).bit_count()for w in words]);q=counts+np.array([6,4,2,0])-N[:,None]
direct=np.sum((4*q-q.sum(axis=1)[:,None])**2,axis=1) #64B
def evaluate_integer(poly,scale):
 out=np.zeros(len(words),dtype=np.int64)
 for mask,c in poly.items():
  val=c*scale;assert val.denominator==1
  out+=int(val)*((words&mask)==mask)
 return out
assert np.array_equal(evaluate_integer(B,64),direct)
assert np.array_equal(evaluate_integer(B2,4096),direct**2)
groups={}
for degree in range(5):
 masks=sorted({m for m in B}|{m for m in B2})
 masks=[m for m in masks if m.bit_count()==degree]
 ix=np.array([[i for i in range(16)if(m>>i)&1]for m in masks],dtype=int)
 groups[degree]=(ix,np.array([float(B.get(m,0))for m in masks]),np.array([float(B2.get(m,0))for m in masks]))

def moments(P):
 b=0j;s2=0j
 for degree,(ix,c,d)in groups.items():
  if not len(c):continue
  determinants=np.ones(len(c))if degree==0 else np.linalg.det(P[ix[:,:,None],ix[:,None,:]])
  b+=c@determinants;s2+=d@determinants
 assert abs(b.imag)<1e-10 and abs(s2.imag)<1e-10
 return float(b.real),float(s2.real)

A,C,_=src.one_particle(0)
DELTA=1e-10;MOMENT_ROUND=1e-6
rows=[];bounds={n:np.full(len(AMPS),np.inf)for n in [7,8]}
started=time.time()
for j in range(GRID):
 theta=2*np.pi*j/GRID
 H=A+C*np.exp(1j*theta)+C.conj().T*np.exp(-1j*theta)
 values,U=np.linalg.eigh(H)
 residual=float(np.linalg.norm(H@U-U*values));orth=float(np.linalg.norm(U.conj().T@U-np.eye(16)))
 assert residual<1e-11 and orth<1e-12
 theta_iv=2*mp.iv.pi*j/GRID
 trial_terms=[]
 for a in AMPS:
  c=mp.iv.mpf(a)*mp.iv.cos(theta_iv)-mp.iv.mpf(a)**2*mp.iv.sin(theta_iv)**2
  trial_terms.append(np.nextafter(float(c.a),-np.inf))
 for n in [7,8]:
  e0=float(sum(values[:n]));gap=float(values[n]-values[n-1])
  base=e0-n*DELTA;maxB=float(max(direct[N==n]))/64
  entry={'grid_index':j,'N':n,'free_ground_lower':base,'free_Fock_gap_diagnostic':gap,
   'B_operator_norm_exact':str(F(int(max(direct[N==n])),64)),
   'one_particle_residual':residual,'one_particle_Gram_error':orth,
   'nondegenerate_control':gap>4*DELTA,'improved_pointwise_lower':base}
  if gap>4*DELTA:
   P=U[:,:n]@U[:,:n].conj().T
   b,s2=moments(P)
   projector_error=2*DELTA/(gap-2*DELTA)
   b_error=4*maxB*math.sqrt(n)*projector_error+MOMENT_ROUND
   s2_error=4*maxB**2*math.sqrt(n)*projector_error+MOMENT_ROUND
   b_lower=max(0,b-b_error);gap_lower=gap-2*DELTA
   off_squared=max(0,s2+s2_error-b_lower**2)
   increase=max(0,(b_lower+gap_lower-math.sqrt((b_lower-gap_lower)**2+4*off_squared))/2)
   entry.update({'B_expectation':b,'B_squared_expectation':s2,'projector_error_bound':projector_error,
    'B_expectation_error':b_error,'B_squared_error':s2_error,
    'off_block_norm_squared_upper':off_squared,'two_block_increase':increase,
    'improved_pointwise_lower':base+increase})
  for ia,a in enumerate(AMPS):
   bounds[n][ia]=min(bounds[n][ia],np.nextafter(entry['improved_pointwise_lower']+trial_terms[ia],-np.inf))
  rows.append(entry)
 if j%128==0:print('grid',j,'seconds',time.time()-started,flush=True)

sector_results=[]
for n in [7,8]:
 tests=[]
 for ia,a in enumerate(AMPS):
  correction=(16*mp.iv.sqrt(2)+abs(a)+a*a)*mp.iv.pi/GRID
  correction_upper=np.nextafter(float(correction.b),np.inf)
  tests.append({'a':a,'grid_lower':float(bounds[n][ia]),'global_correction_upper':float(correction_upper),
   'full_rotor_lower':float(np.nextafter(bounds[n][ia]-correction_upper,-np.inf))})
 sector_results.append({'N':n,'all_Barta_controls':tests,'best_lower':max(t['full_rotor_lower']for t in tests)})
sector_results.append({'N':9,'particle_hole_of':7,'best_lower':sector_results[0]['best_lower']})
base=json.loads((ROOT/'inputs/base_bounds.json').read_text())
all_sectors=[]
for r in base['all_number_sectors']:
 new=next((x['best_lower']for x in sector_results if x['N']==r['N']),r['best_Barta_infinite_rotor_lower'])
 all_sectors.append({'N':r['N'],'base_lower':r['best_Barta_infinite_rotor_lower'],
  'complete_source_lower':max(new,r['best_Barta_infinite_rotor_lower'])})
prepared=json.loads((ROOT/'inputs/prepared_results.json').read_text());consumers=[]
for case in prepared:
 m=case['mass'];lo8=all_sectors[8]['complete_source_lower']-4*abs(m)
 charged=min(all_sectors[7]['complete_source_lower'],all_sectors[9]['complete_source_lower'])-4*abs(m)
 gap=charged-case['variational_energy_upper']
 ratio=case['charged_infinite_resolvent_gap_lower']/gap
 # Remove/add the stored fixed product margin; each solution error is
 # residual/g and therefore rescales by this exact ratio of admitted gaps.
 err=(case['infinite_prepared_Green_error_upper']-1e-8)*ratio+1e-8
 consumers.append({'case':case['case'],'neutral_energy_excess_upper':case['variational_energy_upper']-lo8,
  'charged_resolvent_gap_lower':gap,'prepared_Green_error_upper':err,
  'exact_neutral_fidelity_certified':False,'true_vacuum_propagator_certified':False})
result={'unchanged_source_manifest':'9cbfad0966962e59ad6424a59925cae61ca53b7f4dbb553e62d37228d507d0ff',
 'occupation_words_verified':65536,'B_polynomial_terms':len(B),'B_squared_polynomial_terms':len(B2),
 'polynomial_max_degrees':[max(x.bit_count()for x in B),max(x.bit_count()for x in B2)],
 'moment_projector_arithmetic':{'one_body_comparison_envelope':DELTA,'moment_roundoff_envelope':MOMENT_ROUND},
 'grid':GRID,'refined_sectors':sector_results,'every_number_sector':all_sectors,
 'prepared_response_consumers':consumers,'elapsed_seconds':time.time()-started,
 'new_Hamiltonian_parameter':False,'empirical_mass_used':False}
(ROOT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'ALL_POINTWISE_CONTROLS.json').write_text(json.dumps(rows,indent=2)+'\n')
(ROOT/'OCCUPATION_POLYNOMIALS.json').write_text(json.dumps({'B':{str(m):str(c)for m,c in B.items()},'B2':{str(m):str(c)for m,c in B2.items()}},indent=2)+'\n')
print(json.dumps(result,indent=2))
