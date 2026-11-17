"""Two spatial cells, complete wall CAR, and exact Gauss electric fields."""
from pathlib import Path
import importlib.util
import math
import numpy as np
from scipy import sparse
spec=importlib.util.spec_from_file_location('frozen_loop',Path(__file__).parent/'inputs/loop_source.py')
loop=importlib.util.module_from_spec(spec);spec.loader.exec_module(loop)
S=2;L=2;M=2*S*L

def matrices(m=0):
 Aone,Cone,Eone,Uone=loop.matrices(S,m)
 A=np.kron(np.eye(L),Aone)
 # Each oriented spatial link carries its own beta/2 endpoint term,
 # already summed in Aone. The interior and closing links are distinct.
 A[:2*S,2*S:]+=Cone;A[2*S:,:2*S]+=Cone.conj().T
 C=np.zeros_like(A);C[2*S:,:2*S]=Cone
 E=np.kron(np.eye(L),Eone)
 U=np.kron(np.eye(L*S),loop.X)
 for t in [0,.125,1.5]:
  H=A+C*np.exp(1j*t)+C.conj().T*np.exp(-1j*t)
  Hminus=A+C*np.exp(-1j*t)+C.conj().T*np.exp(1j*t)
  assert np.max(abs(U@Hminus.T@U.conj().T+H))<1e-12
 return A,C,E

def block(N,K,m=0):
 A,C,E=matrices(m);AN=loop.lift(A,N);CN=loop.lift(C,N)
 words=loop.words(M,N);n=len(words);flux=np.arange(-K,K+1)
 n0=np.array([(w&((1<<(2*S))-1)).bit_count()for w in words])
 q0=n0+S-N
 electric=((flux[:,None]+q0[None,:])**2+flux[:,None]**2)/4
 shift=sparse.diags(np.ones(2*K),-1,format='csr')
 H=sparse.kron(sparse.eye(2*K+1),sparse.csr_matrix(AN))+sparse.kron(shift,sparse.csr_matrix(CN))+sparse.kron(shift.T,sparse.csr_matrix(CN.conj().T))
 H+=sparse.diags(electric.ravel())
 return H.tocsr(),AN,CN,E,q0

def particle_hole(N,K):
 low=loop.words(M,N);high=loop.words(M,M-N);index={w:i for i,w in enumerate(high)}
 rows=[];cols=[];values=[]
 for j,w in enumerate(low):
  full=(1<<M)-1;sign=1
  for a in reversed([a for a in range(M)if(w>>a)&1]):
   b=a^1;sign*=(-1)**((full&((1<<b)-1)).bit_count());full^=1<<b
  rows.append(index[full]);cols.append(j);values.append(sign)
 C=sparse.csr_matrix((values,(rows,cols)),shape=(len(high),len(low)))
 return sparse.kron(sparse.csr_matrix(np.eye(2*K+1)[::-1]),C,format='csr')

def creation(N,coefficient):
 return loop.creation(S*L,N,coefficient)

def electric_outside_lower(N,K):
 # Enumerate every allowed occupation charge. The convex quadratic has
 # its exterior minimum at a boundary or a nearest exterior vertex.
 words=loop.words(M,N);qs={((w&15).bit_count()+S-N)for w in words}
 candidates=[]
 for q in qs:
  center=-q/2
  ells={-(K+1),K+1,math.floor(center),math.ceil(center)}
  candidates += [((ell+q)**2+ell**2)/4 for ell in ells if abs(ell)>=K+1]
 return min(candidates)

def coupling_norm_bound(N):
 n=min(S,N,M-N)
 return 8*math.sqrt(((n+1)**2)//4)
