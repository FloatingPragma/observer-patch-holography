"""Full fixed-number CAR blocks of the declared compact wall-loop source."""
from itertools import combinations
import math
import numpy as np
from scipy import sparse

X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex);PP=(np.eye(2)+X)/2;PM=(np.eye(2)-X)/2
h=1/8

def words(M,N):return tuple(sum(1<<j for j in c)for c in combinations(range(M),N))

def lift(one,N):
 M=len(one);basis=words(M,N);index={w:i for i,w in enumerate(basis)}
 out=np.zeros((len(basis),len(basis)),complex)
 for j,w in enumerate(basis):
  for a,b in zip(*np.nonzero(one)):
   a,b=int(a),int(b)
   if not(w>>b)&1:continue
   v=w^(1<<b);sign=(-1)**((w&((1<<b)-1)).bit_count())
   if(v>>a)&1:continue
   sign*=(-1)**((v&((1<<a)-1)).bit_count())
   out[index[v|(1<<a)],j]+=sign*one[a,b]
 return out

def matrices(S,m):
 M=2*S;A=np.kron(np.eye(S),Z);T=-(Z+1j*Y)/2
 for j in range(S-1):
  A[2*j:2*j+2,2*j+2:2*j+4]+=T
  A[2*j+2:2*j+4,2*j:2*j+2]+=T.conj().T
 E=np.zeros((M,2),complex);E[:2]=PP;E[-2:]=PM
 A+=h*m*(E@Z@E.conj().T)
 C=-np.kron(np.eye(S),Z+1j*X)/2
 U=np.kron(np.eye(S)[::-1],Y)
 assert np.max(abs(U@A.T@U.conj().T+A))<1e-12
 assert np.max(abs(U@C.T@U.conj().T+C))<1e-12
 return A/h,C/h,E,U

def block(S,N,K,m):
 A,C,E,U=matrices(S,m);An=lift(A,N);Cn=lift(C,N)
 n=len(An);flux=np.arange(-K,K+1)
 shift=sparse.diags(np.ones(2*K),-1,format='csr')
 H=sparse.kron(sparse.eye(2*K+1),sparse.csr_matrix(An))+sparse.kron(shift,sparse.csr_matrix(Cn))+sparse.kron(shift.T,sparse.csr_matrix(Cn.conj().T))
 H+=sparse.kron(sparse.diags(flux**2/4),sparse.eye(n))
 return H.tocsr(),An,Cn,E

def creation(S,N,coefficient):
 low=words(2*S,N);high=words(2*S,N+1);idx={w:i for i,w in enumerate(high)}
 out=np.zeros((len(high),len(low)),complex)
 for j,w in enumerate(low):
  for a,c in enumerate(coefficient):
   if(w>>a)&1:continue
   sign=(-1)**((w&((1<<a)-1)).bit_count())
   out[idx[w|(1<<a)],j]+=sign*c
 return sparse.csr_matrix(out)

def exterior_bound(S,N,K,m):
 # A at m=0 is a length2S open chain with eigenvalues2cos(j pi/(2S+1))/h.
 # Boundary input changes its fixed-N lift by norm at most |m|.
 n=min(N,2*S-N)
 amin=-16*sum(math.cos(j*math.pi/(2*S+1))for j in range(1,n+1))-abs(m)
 # The charged shift is a sum of n singly occupied two-state raising operators.
 c=8*math.sqrt(((n+1)**2)//4)
 d=(K+1)**2/4+amin-2*c
 return amin,c,d
