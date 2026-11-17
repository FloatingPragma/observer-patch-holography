"""Matrix-free exact-Gauss actions for all sixteen CAR modes."""
from itertools import combinations
from functools import lru_cache
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
L=4;S=2;M=16;H_CLOCK=1/8;MAX_K=4;MAX_SECTOR_DIM=250000
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex);PP=(np.eye(2)+X)/2;PM=(np.eye(2)-X)/2
PARITY=np.array([(-1)**j.bit_count()for j in range(1<<M)],dtype=np.int8)

@lru_cache(None)
def basis(N):
 words=np.array([sum(1<<j for j in c)for c in combinations(range(M),N)],dtype=np.int32)
 index=np.full(1<<M,-1,dtype=np.int32);index[words]=np.arange(len(words))
 return words,index

def lift(one,N):
 words,index=basis(N);rows=[];cols=[];data=[]
 for a,b in zip(*np.nonzero(one)):
  a,b=int(a),int(b)
  if a==b:
   sel=np.flatnonzero((words>>a)&1);target=sel;sign=np.ones(len(sel),dtype=int)
  else:
   sel=np.flatnonzero(((words>>b)&1)&(1-((words>>a)&1)))
   target=index[words[sel]^(1<<a)^(1<<b)]
   lo,hi=sorted((a,b));between=((1<<hi)-1)^((1<<(lo+1))-1)
   sign=PARITY[words[sel]&between]
  rows.append(target);cols.append(sel);data.append(one[a,b]*sign)
 n=len(words)
 return sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n)).tocsr()

def one_particle(m=0):
 A=np.zeros((M,M),complex);C=np.zeros_like(A)
 T=-(Z+1j*Y)/2;hop=-(Z+1j*X)/2
 for cell in range(L):
  for wall in range(S):
   j=4*cell+2*wall;A[j:j+2,j:j+2]+=Z
  j=4*cell;A[j:j+2,j+2:j+4]+=T;A[j+2:j+4,j:j+2]+=T.conj().T
  E=np.zeros((M,2),complex);E[j:j+2]=PP;E[j+2:j+4]=PM
  A+=H_CLOCK*m*E@Z@E.conj().T
 for cell in range(L):
  for wall in range(S):
   j=4*cell+2*wall;k=4*((cell+1)%L)+2*wall
   if cell<L-1:
    A[j:j+2,k:k+2]+=hop;A[k:k+2,j:j+2]+=hop.conj().T
   else:C[j:j+2,k:k+2]+=hop
 U=np.kron(np.eye(L*S),X)
 assert np.max(abs(U@A.T@U.T+A))<1e-12
 assert np.max(abs(U@C.T@U.T+C.conj().T))<1e-12
 E=np.zeros((M,2),complex);E[:2]=PP;E[2:4]=PM
 return A/H_CLOCK,C/H_CLOCK,E

class Sector:
 def __init__(self,N,K,m=0):
  assert 0<=N<=M and 0<=K<=MAX_K
  self.N=N;self.K=K;self.words,self.index=basis(N)
  self.n=len(self.words);self.f=2*K+1;self.dim=self.n*self.f
  assert self.dim<=MAX_SECTOR_DIM
  A,C,self.E=one_particle(m);self.A=lift(A,N);self.C=lift(C,N);self.Ct=self.C.conj().T
  self.flux=np.arange(-K,K+1)
  q=[]
  for x in range(L):
   mask=(1<<(4*(x+1)))-1
   counts=np.array([int(w&mask).bit_count()for w in self.words])
   q.append(counts+S*(L-x-1)-N)
  self.q=np.array(q,dtype=int).T
  self.electric=np.sum((self.flux[:,None,None]+self.q[None,:,:])**2,axis=2)/4
  self.diagonal=self.electric+self.A.diagonal().real[None,:]
  self.operator=LinearOperator((self.dim,self.dim),matvec=self.matvec,rmatvec=self.matvec,matmat=self.matmat,dtype=np.complex128)
 def matmat(self,V):
  k=V.shape[1];v=V.reshape(self.f,self.n,k)
  y=(self.A@v.transpose(1,0,2).reshape(self.n,-1)).reshape(self.n,self.f,k).transpose(1,0,2)
  y=y+self.electric[:,:,None]*v
  if self.f>1:
   y[1:]+=(self.C@v[:-1].transpose(1,0,2).reshape(self.n,-1)).reshape(self.n,self.f-1,k).transpose(1,0,2)
   y[:-1]+=(self.Ct@v[1:].transpose(1,0,2).reshape(self.n,-1)).reshape(self.n,self.f-1,k).transpose(1,0,2)
  return y.reshape(self.dim,k)
 def matvec(self,v):return self.matmat(np.asarray(v).reshape(self.dim,1)).ravel()
 def shifted(self,e):
  return LinearOperator((self.dim,self.dim),matvec=lambda x:self.matvec(x)-e*x,rmatvec=lambda x:self.matvec(x)-e*x,dtype=np.complex128)
 def exterior(self,v):
  a=np.asarray(v).reshape(self.f,self.n,-1)
  return np.vstack([self.Ct@a[0],self.C@a[-1]])

def create(N,coeff):
 low,_=basis(N);high,idx=basis(N+1);rows=[];cols=[];vals=[]
 for a,c in enumerate(coeff):
  if not c:continue
  selected=np.flatnonzero(1-((low>>a)&1))
  rows.append(idx[low[selected]|(1<<a)]);cols.append(selected)
  vals.append(c*PARITY[low[selected]&((1<<a)-1)])
 return sparse.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(high),len(low))).tocsr()

@lru_cache(None)
def particle_hole_fermion(N):
 low,_=basis(N);high,idx=basis(M-N);rows=[];vals=[]
 for w in low:
  w=int(w);full=(1<<M)-1;sign=1
  for a in reversed([a for a in range(M)if(w>>a)&1]):
   b=a^1;sign*=int(PARITY[full&((1<<b)-1)]);full^=1<<b
  rows.append(idx[full]);vals.append(sign)
 return sparse.csr_matrix((vals,(rows,np.arange(len(low)))),shape=(len(high),len(low)))

def between(op,state,n,f):
 return (op@state.reshape(f,n).T).T.ravel()
