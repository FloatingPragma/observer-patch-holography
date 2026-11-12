"""Independent complete CAR current lift and quantum-radiation accounting."""
from pathlib import Path
import hashlib,json
from itertools import product
import numpy as np
from scipy import sparse as sp
from scipy.linalg import expm
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3];TARGET=ROOT/'model/openwave_radiative_pair'

def annihilate(mode):
 columns=np.array([mask for mask in range(4096)if mask&(1<<mode)])
 rows=columns^(1<<mode);sign=np.array([(-1)**int(mask&((1<<mode)-1)).bit_count()for mask in columns])
 return sp.csr_matrix((sign,(rows,columns)),shape=(4096,4096))

def main():
 sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.]);sigma=(sx,sy,sz)
 frames=[]
 for p in (0.,.5,-.5):
  U=np.zeros((4,4),complex)
  for band,sign in enumerate((-1,1)):
   for spin in (0,1):
    h=1 if spin==0 else-1;E=sign*np.sqrt(1+p*p);v=np.array([1,E+h*p]);v/=np.linalg.norm(v)
    U[[spin,spin+2],2*band+spin]=v
  assert np.linalg.norm(U.conj().T@U-np.eye(4))<1e-14;frames.append(U)
 c=[annihilate(j)for j in range(12)]
 def lift(matrix):
  result=sp.csr_matrix((4096,4096),dtype=complex)
  for i,j in zip(*np.nonzero(abs(matrix)>1e-14)):result+=matrix[i,j]*(c[i].conj().T@c[j])
  return result
 ground=sum(1<<i for i in (0,1,4,5,8,9));vacuum=np.eye(1,4096,ground,dtype=complex).ravel()
 A=np.zeros((12,12),complex);A[2,0]=A[3,1]=1
 pair=lift(A)@vacuum/np.sqrt(2);assert abs(np.vdot(pair,pair)-1)<1e-14
 currents={}
 for q,pol in product((1,2),range(2)):
  J=np.zeros((12,12),complex);alpha=-np.kron(sz,sigma[pol])
  for p in range(3):J[4*((p-q)%3):4*((p-q)%3)+4,4*p:4*p+4]=frames[(p-q)%3].conj().T@alpha@frames[p]
  assert np.linalg.norm(J.conj().T@J-np.eye(12))<1e-14
  currents[q,pol]=lift(J)
 coefficients=[]
 for x,y in product(range(2),repeat=2):
  a,b=currents[1,x],currents[2,y]
  value=np.vdot(vacuum,(a@b+b@a)@pair)
  expected=0 if x==y else (1 if x==0 else-1)*2j*np.sqrt(2/5)
  assert abs(value-expected)<1e-14;coefficients.append([value.real,value.imag])
 spin=[]
 for sigma_axis in sigma:
  S=np.zeros((12,12),complex)
  for p in range(3):S[4*p:4*p+4,4*p:4*p+4]=frames[p].conj().T@np.kron(np.eye(2),sigma_axis/2)@frames[p]
  assert np.linalg.norm(S@A-A@S)<1e-14
  total=lift(S);pv=float(np.linalg.norm(total@pair)**2);gv=float(np.linalg.norm(total@vacuum)**2)
  spin.append({'pair_variance':pv,'vacuum_variance':gv})
 assert np.max(abs(np.array([r['pair_variance']for r in spin])-[.2,.2,0]))<1e-14
 # Independent finite oscillator and full photon-dephasing energy identity.
 rng=np.random.default_rng(139);dim=4;cap=2
 X=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim));matter=(X+X.conj().T)/2
 J=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim));a=np.diag(np.sqrt([1.,2.]),1)
 H0=np.kron(matter,np.eye(3))+np.kron(np.eye(dim),np.diag([0,1,2]));V=np.kron(J,a.conj().T)+np.kron(J.conj().T,a)
 H=H0+.03*V;v=rng.normal(size=12)+1j*rng.normal(size=12);v/=np.linalg.norm(v)
 v=expm(-1j*H*.7)@v;rho=np.outer(v,v.conj());dephased=np.zeros_like(rho)
 for n in range(3):
  P=np.kron(np.eye(dim),np.diag([int(j==n)for j in range(3)]));dephased+=P@rho@P
 work=np.trace(H@(dephased-rho)).real;interaction=.03*np.trace(V@rho).real
 assert abs(work+interaction)<1e-14
 # Exact rational finite-sector norm and Taylor certificate, without calling producer.
 eta=F(1,20);t=F(1,10**8);B=25
 assert F(4)+4*F(223607,100000)+8<21 # sqrt5 <2.23607
 assert 48*eta*F(141422,100000)<4 # sqrt2 <1.41422
 leading=F(3,5)*eta*eta*t*t;tail=F(B**3,6)*t**3;lower=2*(leading-tail)**2
 assert lower>F(4,10**38)
 print('INDEPENDENT_RADIATION_MATHEMATICS_PASS',flush=True)
 receipt_path=TARGET/'results.json'
 if not receipt_path.exists():raise RuntimeError('final producer receipt required for sealing review')
 receipt=json.loads(receipt_path.read_text())
 for p,h in receipt['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 assert F(receipt['exact_short_time_certificate']['original_two_photon_radiation_probability_lower'])==lower
 paths=[p for p in TARGET.iterdir()if p.suffix in ('.py','.md','.tex','.json')]
 report={'verdict':'PASS','independent_full_fermion_Fock_dimension':4096,
 'native_current_unitary_one_body_checks':4,'all_four_complete_two_photon_coefficients':coefficients,
 'spin_scalar_pair_generator_commutators':3,'retained_source_sea_spin_variances':spin,
 'independent_photon_number_dephasing_work_identity':True,'finite_total_sector_Taylor_probability_lower':str(lower),
 'scientific_scope':'Declared native-current coupling to finite transverse cochain modes converts a specified scalar particle-hole excitation to the free-source reference plus two quanta. The free reference is not the interacting ground; spatial scattering, physical particle identification, local Gauss QED and empirical coupling are not derived.',
 'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}}
 (Path(__file__).resolve().parents[1]/'radiation_review.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print('PASS_INDEPENDENT_NATIVE_CURRENT_TWO_PHOTON_AUDIT')
if __name__=='__main__':main()
