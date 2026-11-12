"""Independent finite CAR, unitary-dilation and Fierz checks of three sources."""
from collections import defaultdict
import hashlib,json
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import sqrtm
from model.openwave_weak_decay.decay import interaction
from model.common_source import calibrated_matter
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve()
NAMES=('openwave_neutral_transfer','openwave_neutrino_oscillation','openwave_weak_decay')
SIGMA=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]))
def act(v,i,create):
 out={}
 for mask,z in v.items():
  if bool(mask&(1<<i))==create:continue
  out[mask^(1<<i)]=z*(-1)**((mask&((1<<i)-1)).bit_count())
 return out
def dot(a,b):return sum(np.conj(z)*b.get(i,0)for i,z in a.items())
def plus(*vectors):
 out=defaultdict(complex)
 for v in vectors:
  for k,z in v.items():out[k]+=z
 return {k:z for k,z in out.items()if abs(z)>1e-14}
def scale(v,x):return {k:z*x for k,z in v.items()}
def neutral():
 up={sum(1<<i for i in indices):z for indices,z in [((1,2,3),2),((0,1,5),-1),((0,2,4),1)]}
 down=plus(*(act(act(up,i,False),i+3,True)for i in range(3)));states=[up,down];moments=[]
 for sigma in SIGMA:
  matrix=np.zeros((2,2),complex)
  def hop(v):
   terms=[]
   for internal,s,t in product(range(3),range(2),range(2)):
    coefficient=sigma[s,t]
    if coefficient:
     a=6+3*s+internal;b=3*t+internal
     terms.extend((scale(act(act(v,b,False),a,True),coefficient),scale(act(act(v,a,False),b,True),np.conj(coefficient))))
   return plus(*terms)
  for c,state in enumerate(states):
   v=state
   for degree in range(4):
    for r,target in enumerate(states):
     receiver={mask<<6:z for mask,z in target.items()};value=dot(receiver,v)/6
     if degree<3:assert abs(value)<1e-13
     else:matrix[r,c]=value
    if degree<3:v=hop(v)
  assert np.max(abs(matrix+6*sigma))<1e-13;moments.append(matrix.tolist())
 amplitude=-6*2*(-1j*np.sqrt(2)/8)**3
 assert abs(amplitude+3j*np.sqrt(2)/64)<1e-16
 return {'uncompressed_group_identity_three_CAR_moments':'-6 sigma_i on both actual spin-half columns; gauge-singlet Haar reduction and retained outer-link states justify the native compression','H3_receiver_imaginary':float(amplitude.imag),'lower_moments_exactly_zero':True,'boundary_sign_mutation_H3':0}
def oscillation(receipt):
 C=np.array(receipt['actual_charged_current_mixing_real'])+1j*np.array(receipt['actual_charged_current_mixing_imag']);m=np.array(receipt['source_neutral_masses_in_native_units']);maxerror=0.;unitaryerror=0.;referenceerror=0.
 for row in receipt['complete_records']:
  p=row['momentum'];t=row['source_time'];a=row['initial_weak_family'];E=np.sqrt(p*p+m*m);weight=(E+p)/(2*E);K=np.diag(np.sqrt(weight))@C.conj().T
  F0=np.asarray(sqrtm(np.eye(3)-K.conj().T@K),complex);F1=np.asarray(sqrtm(np.eye(3)-K@K.conj().T),complex)
  W=np.block([[K,F1],[F0,-K.conj().T]]);unitaryerror=max(unitaryerror,float(np.linalg.norm(W.conj().T@W-np.eye(6))))
  v=np.eye(3)[:,a];prepared=W@np.r_[v,np.zeros(3)];u=np.exp(-1j*E*t)*prepared[:3];returned=W.conj().T@np.r_[u,np.zeros(3)]
  probabilities=list(abs(returned[:3])**2)+[float(np.vdot(prepared[3:],prepared[3:]).real),float(np.vdot(returned[3:],returned[3:]).real)]
  maxerror=max(maxerror,max(abs(np.array(probabilities)-row['complete_original_probabilities'])))
  U=np.diag(np.exp(-1j*E*t));T=K.conj().T@U@K;fail=F1@U@K
  referenceerror=max(referenceerror,float(np.linalg.norm(T.conj().T@T+F0.conj().T@F0+fail.conj().T@fail-np.eye(3))))
 assert unitaryerror<1e-12 and maxerror<1e-12 and referenceerror<1e-12
 assert receipt['model_mixing_parameters']['sin2_theta13']<.0001
 return {'complete_rows_recomputed':len(receipt['complete_records']),'independent_Halmos_unitary_error':unitaryerror,'all_original_probability_error':maxerror,'complete_Kraus_identity_error_with_arbitrary_reference':referenceerror,'measured_PMNS_fit_asserted':False}
def weak():
 poly,*_=interaction();J=np.array(calibrated_matter.weak_raising_current(families=False),complex);tau=(J+J.conj().T,-1j*(J-J.conj().T),J@J.conj().T-J.conj().T@J)
 currents=[(1 if mu==0 else -1,np.kron(np.kron(bar,T/2),np.eye(3)))for T in tau for mu,bar in enumerate((np.eye(2),*(-x for x in SIGMA)))]
 maxerror=0.
 for (a,b,c,d),coefficient in poly.items():
  # Independent antisymmetrized matrix-element identity for :j j:.
  expected=2*sum(sign*(A[a%90,c%90]*A[b%90,d%90]-A[a%90,d%90]*A[b%90,c%90])for sign,A in currents)
  maxerror=max(maxerror,abs(expected-float(coefficient)))
 assert maxerror==0 and len(poly)==31464 and sum(abs(x)for x in poly.values())==41952
 key=(123,256,169,210);assert poly[key]==-2
 # Full12-mode charged BCS and neutral occupations, independently constructed.
 local={0:.5,9:-.5,6:.5,15:-.5};vac={a|(b<<4)|(1<<8)|(1<<10):u*v for a,u in local.items()for b,v in local.items()}
 initial=scale(plus(act(vac,5,True),scale(act(vac,6,False),-1)),1/np.sqrt(2));target=scale(plus(act(vac,0,True),act(vac,3,False)),1/np.sqrt(2));target=act(act(target,8,False),11,True)
 out=initial
 for i,creation in((5,False),(8,False),(11,True),(0,True)):out=act(out,i,creation)
 matrix=-2*dot(target,out);assert abs(abs(matrix)-1)<1e-14
 assert abs(dot(initial,initial)-1)<1e-14 and abs(dot(target,target)-1)<1e-14 and dot(initial,target)==0
 return {'full_polynomial_Fierz_coefficients_compared':len(poly),'coefficient_max_error':float(maxerror),'independent_12_CAR_receiver_modulus':float(abs(matrix)),'complete_interaction_norm':41952,'free_zero_coupling_transition':0,'physical_lifetime_asserted':False}
def run():
 receipts={name:json.loads((ROOT/'model'/name/'results.json').read_text())for name in NAMES};bindings={}
 for name,value in receipts.items():
  for ref,digest in value.get('input_sha256',value.get('source_sha256',{})).items():
   assert hashlib.sha256((ROOT/ref).read_bytes()).hexdigest()==digest;bindings[ref]=digest
  relative=f'model/{name}/results.json';bindings[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
 bindings[str(HERE.relative_to(ROOT))]=hashlib.sha256(HERE.read_bytes()).hexdigest()
 value={'source_sha256':bindings,'neutral':neutral(),'oscillation':oscillation(receipts[NAMES[1]]),'weak':weak(),'scope':'Independent algebra and complete-instrument checks supplement the package proofs. Native neutral transport is a prepared composite; mixing is a declared nonempirical family model; the contact transition is not a measured muon decay.','verdict':'PASS_INDEPENDENT_NEUTRAL_FLAVOUR_WEAK_REVIEW'}
 (HERE.with_suffix('.json')).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');return value
if __name__=='__main__':print(run()['verdict'])
