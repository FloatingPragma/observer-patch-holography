"""Independent oscillator quadrature, rotations, spectra and current checks."""
import hashlib,json
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import eigh,solve_sylvester
from scipy.special import eval_hermite,factorial,roots_hermite
from model.openwave_quantum_rotational.source import source
from model.openwave_quantum_rotational.__main__ import pins
HERE=Path(__file__).resolve();ROOT=HERE.parents[3];PACKAGE=HERE.parents[1]
def eps(i,j,k):
 if len({i,j,k})<3:return 0
 return 1 if(i,j,k)in((0,1,2),(1,2,0),(2,0,1))else -1
def osc(cap):
 basis=sorted((p for p in product(range(cap+1),repeat=3)if sum(p)<=cap),key=lambda p:(sum(p),p));where={p:i for i,p in enumerate(basis)};a=[]
 for k in range(3):
  v=np.zeros((len(basis),len(basis)),complex)
  for j,p in enumerate(basis):
   if p[k]:q=list(p);q[k]-=1;v[where[tuple(q)],j]=np.sqrt(p[k])
  a.append(v)
 J=[sum((-1j*eps(i,j,k)*a[j].conj().T@a[k]for j,k in product(range(3),repeat=2)),np.zeros_like(a[0]))for i in range(3)]
 return basis,a,J

def run():
 stored=json.loads((PACKAGE/'results.json').read_text());basis,a,L=osc(2);sigma=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])];I20=np.eye(20)
 p=[(v-v.conj().T)/1j/np.sqrt(2)for v in a];x=[(v+v.conj().T)/np.sqrt(2)for v in a]
 z,w=roots_hermite(8);z=z/np.sqrt(1.2);w=w/np.sqrt(1.2)
 def moment(i,j,power):return np.sum(w*eval_hermite(i,z)*eval_hermite(j,z)*z**power)/np.sqrt(np.pi*2.**(i+j)*factorial(i)*factorial(j))
 radial=[np.array([[np.prod([moment(u[d],v[d],int(d==k))for d in range(3)])for v in basis]for u in basis])for k in range(3)]
 T=sum((np.kron(sigma[k],p[k])for k in range(3)),np.zeros((20,20),complex));M=[];external=[]
 for k in range(3):
  M.append(sum((np.kron(sigma[i],sum((eps(i,k,d)*radial[d]/2 for d in range(3)),np.zeros((10,10))))for i in range(3)),np.zeros((20,20),complex)))
  external.append(sum((np.kron(sigma[i],sum((eps(i,k,d)*x[d]/2 for d in range(3)),np.zeros((10,10))))for i in range(3)),np.zeros((20,20),complex)))
 Jm=[np.kron(np.eye(2),L[k])+np.kron(sigma[k]/2,np.eye(10))for k in range(3)];rows=[]
 for cap in(1,2,3):
  fb,fa,FJ=osc(cap);f=len(fb);N=np.diag([sum(q)for q in fb]);b=np.kron(T,np.eye(f))
  for m,v in zip(M,fa):b+=np.kron(m,(v+v.conj().T)/8)
  eig,u=eigh(b);K=(u*np.sqrt(256+eig**2))@u.conj().T
  H=K+np.kron(np.kron(np.eye(2),np.diag([sum(q)for q in basis])),np.eye(f))+np.kron(I20,N)
  J=[np.kron(v,np.eye(f))+np.kron(I20,fv)for v,fv in zip(Jm,FJ)];J2=sum((j@j for j in J),np.zeros_like(H));e,V=eigh(H);low=V[:,:2];z,r=eigh(low.conj().T@J[2]@low);psi=low@r[:,1]
  source_error=float(np.max(abs(b-source(cap)['B'])));assert source_error<2e-14
  covariance=max(float(np.linalg.norm(H@j-j@H))for j in J);casimir=float(np.linalg.norm((J2-.75*np.eye(len(H)))@low));assert covariance<2e-11 and casimir<2e-12
  assert abs(e[1]-e[0])<1e-12 and e[2]>=17 and e[0]<17
  probe=np.kron(external[2],np.eye(f));X=solve_sylvester(K,K,b@probe+probe@b);g=float(4*e[0]*abs(np.vdot(psi,X@psi)));share=float(np.vdot(psi,np.kron(I20,FJ[2])@psi).real);assert share>1e-6
  # Full angular generator, not its matter restriction, commutes with H.
  missing_field=float(np.linalg.norm(H@np.kron(Jm[0],np.eye(f))-np.kron(Jm[0],np.eye(f))@H));assert missing_field>.01
  rows.append({'cap':cap,'ground_energy':float(e[0]),'next_distinct_energy':float(e[2]),'static_g':g,'Maxwell_Jz':share,'all_axis_covariance_error':covariance,'full_Casimir_error':casimir,'source_matrix_max_error':source_error,'omitted_field_torque_norm':missing_field})
 # Find interval records without relying on the writer's nesting layout.
 def nodes(value):
  if isinstance(value,dict):
   yield value
   for child in value.values():yield from nodes(child)
  elif isinstance(value,list):
   for child in value:yield from nodes(child)
 intervals={v['cap']:v for v in nodes(stored)if 'actual_ground_energy_interval'in v and 'cap'in v}
 assert set(intervals)=={1,2,3}
 for row in rows:
  c=intervals[row['cap']];lo,hi=c['actual_ground_energy_interval'];assert lo-5e-14<=row['ground_energy']<=hi+5e-14
  lo,hi=c['rest_energy_normalized_static_g_interval'];assert lo<=row['static_g']<=hi
 scientific=pins();scientific[str((PACKAGE/'results.json').relative_to(ROOT))]=hashlib.sha256((PACKAGE/'results.json').read_bytes()).hexdigest();scientific[str(HERE.relative_to(ROOT))]=hashlib.sha256(HERE.read_bytes()).hexdigest()
 return {'review':'Independent Gauss-Hermite matrix reconstruction, spectral calculation and Sylvester derivative, compared with source-bound interval receipts. Numerical checks supplement the exact equivariant-band and full-Fock positivity arguments; they are not new interval proofs.','source_sha256':scientific,'checks':rows,'scope':'Declared common charge -1, confined native Clifford/quantum Maxwell extension, one empirical electron anchor. Coupled matter+field spin is exact; Maxwell alone has integer spin, and no physical electron or three-family identity transfers.','verdict':'PASS_INDEPENDENT_QUANTUM_ROTATIONAL_REVIEW'}
if __name__=='__main__':
 value=run();(PACKAGE/'audits/independent_review.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(value['verdict'])
