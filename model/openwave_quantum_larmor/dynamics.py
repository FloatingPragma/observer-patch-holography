"""Full field evolution and exact-angular-sector interval branch witnesses."""
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
import sympy as s
from mpmath import iv
from model.openwave_quantum_rotational.source import source,matter_exact,array
from model.openwave_quantum_rotational.enclosure import exact_sparse,interval_entries,multiply,inner,norm_upper,endpoints,number,normalized_dyadic,scalar

CAP=2;MASS=16

def prior():
 d=json.loads((Path(__file__).resolve().parents[1]/'openwave_quantum_rotational/results.json').read_text())
 return next(x for x in d['interval_certificates']if x['cap']==CAP)

@lru_cache(maxsize=8)
def perturbed(field):
 d=source(CAP);b=float(F(field));external=np.kron(array(matter_exact()['external'][2]),np.eye(len(d['photon']['basis'])));B=d['B']+b*external
 values,vec=np.linalg.eigh(B);kin=(values*values)/(np.sqrt(MASS*MASS+values*values)+MASS)
 number=d['Nm']+d['Nf']
 h=(vec*kin)@vec.conj().T+number;energies,vectors=np.linalg.eigh(h)
 branches={}
 for k in range(2):
  v=vectors[:,k];j=float(np.vdot(v,d['J'][2]@v).real);target=1 if j>0 else -1
  assert abs(j-target/2)<1e-7;branches[target]=(float(MASS+energies[k]),v)
 return d,h,energies,vectors,branches


def branch_certificate(field,sign):
 iv.dps=50;scalar.cache_clear();d,h,e,v,branches=perturbed(field);b=F(field)
 B0,E0,J0,diag=exact_sparse(CAP);Bdict=dict(B0)
 for key,z in E0.items():Bdict[key]=Bdict.get(key,0)+s.Rational(b.numerator,b.denominator)*z
 B=interval_entries(Bdict);J=interval_entries(J0)
 phi,integers,D=normalized_dyadic(branches[sign][1]);target=F(sign,2)
 # Spectral polynomial is the exact Jz projector; rounded witnesses need not be eigenvectors.
 for twice in range(-9,10,2):
  mu=F(twice,2)
  if mu==target:continue
  jphi=multiply(J,phi);phi=[(x-number(mu)*y)/number(target-mu)for x,y in zip(jphi,phi)]
 projection_norm=norm_upper(phi);assert float(projection_norm.a)>.9
 # norm_upper returns an interval containing the exact norm, not just its endpoint.
 phi=[x/projection_norm for x in phi]
 powers=[phi]
 for k in range(6):powers.append(multiply(B,multiply(B,powers[-1])))
 hv=[(MASS+n)*x for n,x in zip(diag,phi)]
 for k in range(1,7):
  c=F(s.binomial(s.Rational(1,2),k))/MASS**(2*k-1)
  hv=[x+number(c)*y for x,y in zip(hv,powers[k])]
 ray=inner(phi,hv).real;residual=norm_upper([x-ray*y for x,y in zip(hv,phi)])
 a=F(5,2)+abs(b)*F(4,5);delta=number(F(33,2048)*a**14/MASS**13)
 high=ray.b+delta;gap=17-high;assert float(gap.a)>0
 lo=ray.a-(residual.b+delta)**2/gap
 return {'field_exact':str(b),'Jz_exact':str(target),'trial_integer_coefficients':integers,'trial_squared_integer_norm':str(D),
  'exact_Jz_projector_eigenvalues':[str(F(k,2))for k in range(-9,10,2)],
  'sixth_polynomial_residual_upper':endpoints(residual)[1],'positive_square_root_remainder_upper':endpoints(delta)[1],
  'same_source_branch_energy_interval':[endpoints(lo)[0],endpoints(high)[1]]}


def time_bound(field,time):
 p=prior();b=abs(F(field));a=F(5,2);e=F(4,5);m=F(16)
 v=a*e*b/m+e*e*b*b/(2*m);L2=(e*e+((a+b*e)*e/m)**2)/m;r=L2*b*b/2
 gap=F(17)-F.from_float(p['actual_ground_energy_interval'][1]);sep=gap-2*v;assert sep>0
 D=r+v*v/sep;eta=2*np.sqrt(2)*float(v/sep)+abs(float(time))*float(D)
 return {'field_exact':str(F(field)),'time':float(time),'Hamiltonian_change_norm_upper':float(v),
  'second_order_effective_energy_error_upper':float(D),'full_source_reference_trace_distance_upper':min(2.,2*eta),
  'analytic_precession_static_g_difference_upper':4*p['actual_ground_energy_interval'][1]*float(D/b)}


def joint_outputs(d,state):
 photon=d['photon'];basis=photon['basis'] if isinstance(photon,dict)else photon
 # Basis order is matter slow, photon fast throughout the declared action.
 numbers=np.array([sum(x)for x in basis]);nfield=len(basis);jx=d['J'][0];out=[]
 for n in range(CAP+1):
  idx=np.array([k for k in range(len(state))if numbers[k%nfield]==n]);values,vectors=np.linalg.eigh(jx[np.ix_(idx,idx)]);amp=vectors.conj().T@state[idx]
  for twice in sorted(set(int(round(2*x))for x in values)):
   selected=np.array([int(round(2*x))==twice for x in values]);prob=float(np.sum(abs(amp[selected])**2))
   out.append({'photon_number':n,'total_Jx':str(F(twice,2)),'probability':prob})
 assert abs(sum(x['probability']for x in out)-np.vdot(state,state).real)<1e-11
 return out

def complex_rows(v):return [[float(z.real),float(z.imag)]for z in v]


def records():
 d=source(CAP);up=d['psi'];down=(d['J'][0]-1j*d['J'][1])@up;down/=np.linalg.norm(down)
 assert abs(np.vdot(up,down))<1e-12
 mu=float(np.vdot(up,d['current_derivative']@up).real);initial=(up+down)/np.sqrt(2);seed=d['seed'];a=np.vdot(up,seed);success=abs(a)**2;failure=seed-a*up
 rows=[]
 for field in('1/1000','1/10000'):
  _,h,energies,vectors,branches=perturbed(field);b=float(F(field));projected=vectors.conj().T@initial;fp=vectors.conj().T@failure
  for tau in(1,10,100):
   t=tau/b;phase=np.exp(-1j*(energies-(d['E']-MASS))*t);actual=vectors@(phase*projected)
   ideal=(np.exp(-1j*b*mu*t)*up+np.exp(1j*b*mu*t)*down)/np.sqrt(2)
   overlap=abs(np.vdot(actual,ideal));distance=2*np.sqrt(max(0.,1-min(1.,overlap)**2));bound=time_bound(field,t)
   assert distance<=bound['full_source_reference_trace_distance_upper']+2e-7
   fail=vectors@(phase*fp);outputs=joint_outputs(d,np.sqrt(success)*actual)
   assert abs(sum(x['probability']for x in outputs)+np.vdot(fail,fail).real-1)<1e-11
   rows.append({'field_exact':field,'scaled_time_b_t':tau,'source_time':t,'bound':bound,
    'numerical_full_pure_state_trace_distance':float(distance),'complete_success_joint_Jx_photon_number_outputs':outputs,
    'failure_probability':float(np.vdot(fail,fail).real),'full_joint_success_state_coefficients':complex_rows(np.sqrt(success)*actual),
    'full_joint_failure_state_coefficients':complex_rows(fail),'idle_compensating_reference_charge':1,
    'full_other_native_Fock_sectors_probability':0.,'numerical_precession_g':2*d['E']*abs(branches[1][0]-branches[-1][0])/b})
 return rows


def certificate():
 p=prior();rows=[]
 for field in('1/100','1/1000'):
  branches=[branch_certificate(field,sign)for sign in(1,-1)]
  plus=iv.mpf(branches[0]['same_source_branch_energy_interval']);minus=iv.mpf(branches[1]['same_source_branch_energy_interval'])
  gap=abs(plus-minus);g=2*iv.mpf(p['actual_ground_energy_interval'])*gap/number(F(field))
  rows.append({'field_exact':field,'branches':branches,'finite_field_precession_g_interval':endpoints(g)})
 return {'schema':'oph.quantum_larmor_same_source.v1','cap':CAP,'unperturbed_ground_energy_interval':p['actual_ground_energy_interval'],
  'static_current_g_interval':p['rest_energy_normalized_static_g_interval'],'certified_finite_field_splittings':rows,
  'complete_finite_time_records':records(),'all_joint_Maxwell_matter_reference_outputs_retained':True,
  'same_rest_unit_current_state_field_and_clock':True,'field_strength_and_source_couplings_physically_calibrated':False,
  'scope':'Full finite quantum source and all field/reference outputs, exact angular sectors and analytic finite-time bounds. The photon cap and confined spatial regulator are retained; no observed-electron g or field calibration is inferred.',
  'verdict':'PASS_QUANTUM_LARMOR_FULL_OUTPUT_CERTIFICATE'}
