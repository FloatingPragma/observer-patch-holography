"""One autonomous source for native matter, harmonic geometry and scalar waves."""
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as sy
from scipy.linalg import eigh,expm
from model.openwave_quantum_geometry.source import energies
KAPPA=100

@lru_cache(None)
def shape(n=7):
 if n not in(5,7,9):raise ValueError('declared odd planar grid required')
 m=n*n;L=sy.zeros(m)
 for x,y in product(range(n),repeat=2):
  i=n*x+y;L[i,i]=4
  for a,b in((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
   if 0<=a<n and 0<=b<n:L[i,n*a+b]=-1
 c=n//2;sites=[n*(c-1)+c,n*(c+1)+c];B=sy.zeros(m,2)
 for j,i in enumerate(sites):B[i,j]=1
 G=L.inv()*B;F=G*G[sites,:].inv();S=F.T*L*F;M=F.T*F
 assert F[sites,:]==sy.eye(2) and S[0,0]==S[1,1] and M[0,0]==M[1,1]
 assert all((L*F)[i,j]==0 for i in range(m)if i not in sites for j in range(2))
 return {'n':n,'L_exact':L,'F_exact':F,'S_exact':S,'M_exact':M,'L':np.array(L,float),'F':np.array(F,float),'S':np.array(S,float),'M':np.array(M,float),'source_sites':sites}

@lru_cache(None)
def geometry(n=7,cap=1):
 if cap not in(1,2):raise ValueError('complete product cap1 or2 required')
 a=np.diag(np.sqrt(np.arange(1,cap+1)),1);one=np.eye(cap+1);operators=[np.kron(a,one),np.kron(one,a)];sh=shape(n)
 M=sh['M'];S=sh['S'];omega=np.sqrt([(S[0,0]+S[0,1])/(M[0,0]+M[0,1]),(S[0,0]-S[0,1])/(M[0,0]-M[0,1])]);mu=np.array([M[0,0]+M[0,1],M[0,0]-M[0,1]])
 c=1/np.sqrt(4*KAPPA*omega*mu);normal=[operators[i]+operators[i].T for i in range(2)]
 Q=(c[0]*normal[0]+c[1]*normal[1],c[0]*normal[0]-c[1]*normal[1]);states=list(product(range(cap+1),repeat=2));Hg=np.diag([omega@v for v in states]);d=len(states)
 # The one-mode coordinate diagonalization gives the entire commuting algebra.
 ev,W=eigh(a+a.T);coordinate_W=np.kron(W,W);qvalues=np.stack([c[0]*ev[i]+c[1]*ev[j] for i,j in states]);rvalues=np.stack([c[0]*ev[i]-c[1]*ev[j]for i,j in states])
 return {'cap':cap,'dim':d,'Q':Q,'Hg':Hg,'omega':omega,'mu':mu,'normal_coefficients':c,'a':operators,'coordinate_W':coordinate_W,'coordinate_values':np.column_stack([qvalues,rvalues]),'vacuum':np.eye(d)[:,0]}

@lru_cache(None)
def field(n=7,cap=1):
 sh=shape(n);g=geometry(n,cap);m=n*n;d=g['dim'];W=g['coordinate_W'];blocks=[];Kblocks=[];qprofiles=[]
 for z in g['coordinate_values']:
  q=sh['F']@z;s=np.exp(2*q);K=s[:,None]*sh['L']*s[None,:];ev,V=eigh(K);assert ev[0]>0
  blocks.append((V*np.sqrt(ev))@V.T);Kblocks.append(K);qprofiles.append(q)
 V=np.kron(W,np.eye(m));diagonal=np.zeros((d*m,d*m))
 for a,B in enumerate(blocks):diagonal[a*m:(a+1)*m,a*m:(a+1)*m]=B
 return {'H':V@diagonal@V.T,'coordinate_blocks':blocks,'K_blocks':Kblocks,'qprofiles':qprofiles}

@lru_cache(None)
def block(n=7,cap=1,left=1,right=0):
 g=geometry(n,cap);E=energies();vac=g['Hg']+E[left]*expm(g['Q'][0])+E[right]*expm(g['Q'][1]);wave=np.kron(vac,np.eye(n*n))+field(n,cap)['H'];return {'vacuum':vac,'wave':wave}

def packet_integers(n=7,reflect=False):
 x,y=np.meshgrid(np.arange(n),np.arange(n),indexing='ij');c=n//2;offset=-1 if reflect else 1
 psi=np.exp(-((x-1)**2+(y-(c+offset))**2)/(2*.8**2)+1j*1.2*x).ravel()
 # These rounded integers define the exact preparation; their Gaussian
 # proposer is not an additional physical input to the accepted state.
 return [(int(round(z.real*2**30)),int(round(z.imag*2**30)))for z in psi]
def packet(n=7,reflect=False):
 integer=np.array([complex(a,b)for a,b in packet_integers(n,reflect)]);return integer/np.linalg.norm(integer)

def diagonalize(n=7,cap=1,left=1,right=0):return eigh(block(n,cap,left,right)['wave'])
def evolution(e,V,psi,t):return V@(np.exp(-1j*t*e)*(V.conj().T@psi))

def probe(n=7,cap=1,left=1,right=0):
 g=geometry(n,cap);psi=np.kron(g['vacuum'],packet(n));e,V=diagonalize(n,cap,left,right);m=n*n;Y=np.tile(np.arange(n),n);rows=[]
 for t in(1,2,3,4,5,6,8):
  state=evolution(e,V,psi,t).reshape(g['dim'],m);p=abs(state)**2;rho=state@state.conj().T
  rows.append((t,float(p.sum(axis=0)@Y),[float(np.trace(rho@Q).real)for Q in g['Q']],float(np.trace(rho@g['Hg']).real),float(np.trace(rho@rho).real)))
 # Repeated complete first-arrival instrument, including retained no-arrival.
 U=(V*np.exp(-.5j*e))@V.conj().T;v=psi.copy();receiver=np.tile(np.repeat(np.arange(n)>=n-2,n),g['dim']);arr=[]
 for step in range(1,17):
  v=U@v;p=float(np.vdot(v[receiver],v[receiver]).real);arr.append(p);v[receiver]=0
 return rows,arr,float(np.vdot(v,v).real)
