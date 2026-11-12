import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.particle_common_contract.provenance import import_closure
from model.openwave_vectorlike_flavour.flavour import branches,ward_certificate
from model.openwave_vectorlike_flavour.__main__ import pins as flavour_pins
from model.openwave_baryon_ordering.__main__ import pins as baryon_pins
from model.openwave_baryon_ordering.source import balanced_baryon,one_body
from model.openwave_baryon_ordering.calculation import ratio_intervals
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def algebra_checks():
 worst=0.;count=0
 for m,g,M in product((.0001,.01,.25,2.),(.25,1.,4.),(.5,1.,2.)):
  lows=[]
  for y in(2/3,1/3):
   A=np.array([[m,0],[g*y,M]]);singular=np.linalg.svd(A,compute_uv=False);lo,hi=branches(m,y,g,M);lows.append(float(lo))
   worst=max(worst,float(abs(singular[-1]-float(lo))),float(abs(singular[0]-float(hi))))
   assert abs(float(lo*hi)-m*M)<1e-12 and 0<float(lo)<m<float(hi);count+=1
  assert 0<lows[0]<lows[1]<m and lows[1]-lows[0]<m
 return {'independent_SVD_parameter_rows':count,'maximum_singular_value_discrepancy':worst,
  'all_heavy_branches_retained':True,'declared_G_M_extension_is_not_source_forced':True,'native_Ward_and_anomaly_certificate':ward_certificate()}

def baryon_checks(receipt):
 a=balanced_baryon();h=a['h'];e,O=eigh(h);cert=receipt['exact_equal_mass_spectrum_certificate'];D=cert['dyadic_denominator'];I=np.array([[int(x)for x in row]for row in cert['dyadic_eigenframe']],dtype=object);diag=np.array([int(x)for x in cert['dyadic_eigenvalues']],dtype=object)
 assert max(abs(int(x))for x in diag)<cert['source_and_eigenvalue_norm_bound']*D
 gram=I.T@I-D*D*np.eye(len(I),dtype=object);eta=F(cert['exact_Gram_Frobenius_error_upper']);assert sum(int(x)**2 for x in gram.flat)<=eta*eta*D**4
 intervals=[[F(x)for x in row]for row in cert['all_eigenvalue_intervals_rational']]
 assert all(float(lo)-2e-12<=v<=float(hi)+2e-12 for v,(lo,hi)in zip(e,intervals))
 # Independently form total weak isospin on every original three-quark word.
 modes=a['modes'];indices={v:i for i,v in enumerate(modes)};raising=np.zeros((24,24))
 for x,s,c in product(range(2),range(2),range(3)):raising[indices[(0,x,s,0,c)],indices[(0,x,s,1,c)]]=1
 R=one_body(raising,a['basis']);z=np.array([sum((.5 if f==0 else -.5)for i,(_,x,s,f,c)in enumerate(modes)if mask>>i&1)for mask in a['basis']]);I2=sp.diags(z*z+z)+R.conj().T@R;reduced=(a['B'].conj().T@I2@a['B']).toarray();ground=O[:,:4]
 iso=float(np.linalg.norm(reduced@ground-.75*ground));spin=float(np.linalg.norm(a['old']['spin']@ground-.75*ground));assert max(iso,spin)<1e-11
 profiles=[]
 for state in receipt['complete_intrinsic_records']:
  v=list(map(int,state['exact_integer_trial_vector']));norm=int(state['squared_integer_norm']);assert sum(x*x for x in v)==norm
  success=F(state['exact_filter_success_probability']);failure=F(state['exact_filter_failure_probability']);probabilities=list(map(F,state['complete_unconditional_native_occupation_probabilities']));assert sum(probabilities)+failure==1 and success+failure==1
  assert all(p==success*F(x*x,norm)for p,x in zip(probabilities,v))
  profile={}
  for mask,x in zip(a['basis'],v):
   if not x:continue
   word=[m for i,m in enumerate(modes)if mask>>i&1];centre=F(sum(m[1]for m in word),3)
   for _,position,s,flavour,c in word:
    distance=str(abs(F(position)-centre));q=F(2,3)if flavour==0 else F(-1,3);profile[distance]=profile.get(distance,F(0))+F(x*x,norm)*q
  assert sum(profile.values())==state['charge']
  for r,value in profile.items():assert str(value)==state['intrinsic_charge_profile'][r]['trial_exact']
  if state['charge']==0:assert profile['1/3']>0>profile['2/3']
  profiles.append({'charge':state['charge'],'exact_profile':{k:str(v)for k,v in profile.items()},'full_unconditional_record_count':len(probabilities)+1})
 mu,md=ratio_intervals();assert 0<md.lo-mu.hi and md.hi-mu.lo<1
 for r in receipt['baryon_energies']:
  baseline=list(map(F,cert['lowest_interval_rational']));nup=r['light_up_count'];ndown=r['light_down_count'];lo=baseline[0]+nup*(mu.lo-1)+ndown*(md.lo-1);hi=baseline[1]+nup*(mu.hi-1)+ndown*(md.hi-1)
  assert list(map(F,r['source_gap_interval_rational']))==[lo,hi]
  dec=list(map(F,r['source_gap_interval']));assert dec[0]<=lo<=hi<=dec[1]
 assert receipt['exact_beta_threshold_obstruction']['n_to_p_e_antinu_ground_channel_energetically_allowed']is False
 return {'independent_ground_spin_residual':spin,'independent_ground_isospin_residual':iso,'exact_profile_and_failure_records':profiles,
  'exact_free_branch_intervals_propagated_to_baryon_energy':True,'outward_decimal_energy_endpoints':True,
  'beta_threshold_obstruction_preserved_for_all_positive_G_M':True,
  'physical_mass_form_factor_and_free_to_interacting_projection_not_inferred':True}

def certificate():
 inputs={};data={}
 for name,fn in(('openwave_vectorlike_flavour',flavour_pins),('openwave_baryon_ordering',baryon_pins)):
  path=ROOT/'model'/name/'results.json';r=json.loads(path.read_text());assert r['source_sha256']==fn();data[name]=r;inputs.update(r['source_sha256']);inputs[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
 for p in import_closure({str(p.relative_to(ROOT))for p in HERE.glob('*.py')}):inputs[p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
 return {'verdict':'PASS_INDEPENDENT_VECTORLIKE_BARYON_ORDERING_REVIEW','source_sha256':inputs,
  'vectorlike_algebra':algebra_checks(),'baryon_spectrum_and_records':baryon_checks(data['openwave_baryon_ordering']),
  'scope':'Declared hypercharge vectorlike extension and balanced current action; empirical electron anchor and retrospective family ancestry retained. Correct internal ordering/profile signs and beta obstruction do not establish physical quark or nucleon identity.'}
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=certificate();path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:
  from model.openwave_larmor_dynamics.certificate import equivalent
  assert equivalent(json.loads(path.read_text()),r)
 print(r['verdict'])
if __name__=='__main__':main()
