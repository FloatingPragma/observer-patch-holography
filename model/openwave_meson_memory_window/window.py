"""One emitted block for every native family, with exact contact injection."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
from model.openwave_meson_ring.ring import hadron,hadron_eigensystem,singlets,G,BAR,ETA,construction,reduction,finite_volume_bound

def operator(factor,A):
 out=np.ones((1,1),complex)
 for j in range(4):out=np.kron(out,A if j==factor else np.eye(2))
 return out

@lru_cache(None)
def local():
 d=hadron();labels=d['hadron_labels'];S=singlets().T;target=[i for i,x in enumerate(labels)if x[0]==0 and x[3]==0]
 parent,e,U=hadron_eigensystem(d,1,1);p=np.zeros(96,complex);p[parent]=U[:,0]
 blocks={};H=d['Hh'].toarray()
 for xq,xa,yq,ya in product(range(2),repeat=4):
  rows=[next(i for i in target if labels[i]==(0,xq,sq,0,xa,sa))for sq,sa in product(range(2),repeat=2)];cols=[next(i for i in target if labels[i]==(0,yq,sq,0,ya,sa))for sq,sa in product(range(2),repeat=2)]
  blocks[xq,xa,yq,ya]=S.conj().T@np.kron(H[np.ix_(rows,cols)],np.eye(4))@S
 Ve=np.zeros((2,2),complex);Vh=Ve.copy()
 for mu in range(4):
  Ve-=.5*ETA[mu]*(S.conj().T@operator(0,BAR[mu])@operator(2,BAR[mu])@S)
  Vh-=.5*ETA[mu]*(S.conj().T@operator(0,BAR[mu])@operator(3,BAR[mu].T)@S)
 injection={}
 j=d['raising'][np.ix_([10,11],[10,11])]
 for x,xa in product(range(2),repeat=2):
  spin=np.zeros(16,complex)
  for mu in range(4):
   had=d['qops'][x][mu][0]@p/d['C'][0,1]
   for sq,sa,se,sn in product(range(2),repeat=4):
    hi=next(i for i in target if labels[i]==(0,x,sq,0,xa,sa));spin[((2*sq+sa)*2+se)*2+sn]+=ETA[mu]*had[hi]*BAR[mu][se,sn]*j[0,1].conjugate()
  injection[x,xa]=S.conj().T@spin
 return dict(d,spin_basis=S,blocks=blocks,Ve=Ve,Vh=Vh,injection=injection,parent_contact_energy=float(e[0]-3*G))

def universal(L):
 assert L>=2 and L%2==0;d=local();q=d['q_two'];n=8*L*L
 def ix(xq,xa,xe,xh,k):return (((2*xq+xa)*L+xe)*L+xh)*2+k
 rows=[];cols=[];data=[];v=np.zeros(n,complex)
 for xq,xa,xe,xh in product(range(2),range(2),range(L),range(L)):
  ids=[ix(xq,xa,xe,xh,k)for k in range(2)]
  localV=3*np.eye(2)+(d['Ve']if xe==xq else 0)+(d['Vh']if xh==xq else 0)
  for i,k in product(range(2),repeat=2):
   z=(1+q if i==k else 0)+G*localV[i,k]
   if z!=0:rows.append(ids[i]);cols.append(ids[k]);data.append(z)
  for k in range(2):
   for coordinate in(0,1):
    for step in(-1,1):
     pe=(xe+step)%L if coordinate==0 else xe;ph=(xh+step)%L if coordinate==1 else xh;rows.append(ids[k]);cols.append(ix(xq,xa,pe,ph,k));data.append(-q/4)
  for yq,ya in product(range(2),repeat=2):
   M=d['blocks'][xq,xa,yq,ya]
   for i,k in product(range(2),repeat=2):
    if abs(M[i,k])>1e-14:rows.append(ids[i]);cols.append(ix(yq,ya,xe,xh,k));data.append(M[i,k])
  if xe==xh==xq:v[ids]=d['injection'][xq,xa]
 H=sp.coo_matrix((data,(rows,cols)),shape=(n,n)).tocsr();assert np.max(abs((H-H.conj().T).data),initial=0)<1e-12
 # Exactly the same even-reflection ordering as the complete native reduction.
 br=[];bc=[];bv=[];count=0
 for xa,xe,xh,k in product(range(2),range(L),range(L),range(2)):
  br.extend((ix(0,xa,xe,xh,k),ix(1,1-xa,(1-xe)%L,(1-xh)%L,k)));bc.extend((count,count));bv.extend((1/np.sqrt(2),1/np.sqrt(2)));count+=1
 B=sp.coo_matrix((bv,(br,bc)),shape=(n,count)).tocsr();D=(B.T@H@B).tocsr();psi=B.T@v
 assert abs(np.linalg.norm(psi)**2-8)<2e-12
 return dict(d,L=L,D=D,psi=np.asarray(psi).ravel(),unreduced_dimension=n,dimension=count)

def compare_native(L):
 d=construction(L);B,H,error=reduction(d);u=universal(L);n=u['dimension'];D=H[6:6+n,6:6+n];a=B.T@d['initial'];C=(H[6:6+n,:6]@a[:6])/G
 derror=float(np.max(abs((D-u['D']).data),initial=0));cerror=float(np.linalg.norm(C-d['C'][0,1]*u['psi']));assert derror<3e-10 and cerror<2e-12
 for f,ell in product(range(3),repeat=2):
  start=6+(3*f+ell)*n;actual=H[start:start+n,start:start+n];shift=d['mass'][f]+d['mass'][ell]-2
  assert np.max(abs((actual-u['D']-shift*sp.eye(n)).data),initial=0)<4e-10
 return {'ring_sites':L,'complete_native_prepared_dimension':H.shape[0],'universal_emitted_block_dimension':n,'full_nine_family_block_numeric_residual':derror,'actual_native_injection_numeric_residual':cerror,'all_nine_blocks_retained':True}

def first_light(u,t):
 M=u['D']-u['parent_contact_energy']*sp.eye(u['dimension']);psi=u['psi'];n=len(psi)
 augmented=sp.bmat([[M,sp.csr_matrix((G*psi)[:,None])],[None,sp.csr_matrix((1,1))]],format='csr');v=np.zeros(n+1,complex);v[-1]=1
 out=expm_multiply(-1j*float(t)*augmented,v)[:n]
 return out


def exact_constants():
 """Algebraic coefficients before dyadic numerical evaluation."""
 import sympy as sy
 from model.openwave_weak_decay.decay import mass_data
 boxes=mass_data()[2];mid=[sum(map(F,b),F(0))/2 for b in boxes];mu=mid[2]/mid[1]
 mu=sy.Rational(mu.numerator,mu.denominator);g=sy.Rational(1,64);q=(mu-2)/2-11*g/4
 root=sy.sqrt(5473);c=sy.sqrt((1+73/root)/2);d=sy.sqrt((1-73/root)/2)
 return {'q':q,'parent_energy':mu-1+(181-root)/48-3*g,'sqrt3':sy.sqrt(3),'c':c,'d':d,'sqrt2c':sy.sqrt(2)*c,'sqrt6c':sy.sqrt(6)*c,'sqrt2d':sy.sqrt(2)*d,'sqrt6d':sy.sqrt(6)*d}

def certified_universal(L):
 """Exact dyadic surrogate, with an interval bound to the full source."""
 from model.openwave_quantum_rotational.enclosure import scalar
 from mpmath import iv
 import sympy as sy
 iv.dps=60;scalar.cache_clear();D=2**40;constants=exact_constants();quant={};errors={}
 for name,z in constants.items():
  center=F(round(float(z.evalf(18))*D),D);interval=scalar(z)
  lo=F.from_float(float(np.nextafter(float(interval.a),-np.inf)));hi=F.from_float(float(np.nextafter(float(interval.b),np.inf)))
  quant[name]=center;errors[name]=max(abs(center-lo),abs(center-hi))
 assert max(errors.values())<F(1,10**11)
 # Integer reflection orbits allow exact assembly without floating sqrt(2).
 n=8*L*L
 def ix(xq,xa,xe,xh,k):return (((2*xq+xa)*L+xe)*L+xh)*2+k
 rows=[];cols=[];data=[];q=float(quant['q']);rt=float(quant['sqrt3']);ep=float(quant['parent_energy'])
 Ve=np.array([[-.5,-rt/2],[-rt/2,-1.5]]);Vh=np.array([[-.5,-rt/2],[-rt/2,.5]])
 for xq,xa,xe,xh in product(range(2),range(2),range(L),range(L)):
  ids=[ix(xq,xa,xe,xh,k)for k in range(2)];localV=3*np.eye(2)+(Ve if xe==xq else 0)+(Vh if xh==xq else 0)
  for i,k in product(range(2),repeat=2):
   hdiag=(2.25+.5*k)if xq==xa else 127/24
   # 127/24 is quantized separately below; all remaining scalar arithmetic is dyadic.
   hdiag=float(F(round(hdiag*D),D))
   z=(1+q-ep+hdiag if i==k else 0)+G*localV[i,k]
   if z:rows.append(ids[i]);cols.append(ids[k]);data.append(z)
  for k in range(2):
   for coordinate in(0,1):
    for step in(-1,1):
     pe=(xe+step)%L if coordinate==0 else xe;ph=(xh+step)%L if coordinate==1 else xh
     rows.append(ids[k]);cols.append(ix(xq,xa,pe,ph,k));data.append(-q/4)
   for yq,ya in((1-xq,xa),(xq,1-xa)):
    rows.append(ids[k]);cols.append(ix(yq,ya,xe,xh,k));data.append(-.125)
 H=sp.coo_matrix((data,(rows,cols)),shape=(n,n)).tocsr();br=[];bc=[];count=0;psi=[]
 for xa,xe,xh,k in product(range(2),range(L),range(L),range(2)):
  br.extend((ix(0,xa,xe,xh,k),ix(1,1-xa,(1-xe)%L,(1-xh)%L,k)));bc.extend((count,count));count+=1
  key=('sqrt2'if k==0 else 'sqrt6')+('c'if xa==0 else 'd');psi.append(float(quant[key])if xe==xh==0 else 0.)
 B=sp.coo_matrix((np.ones(len(br)),(br,bc)),shape=(n,count)).tocsr();M=(B.T@H@B)/2;psi=np.array(psi)
 assert max(np.diff(M.indptr))<=16
 row=float(abs(M).sum(axis=1).max());col=float(abs(M).sum(axis=0).max());assert max(row,col)<206
 # q multiplies the two-particle Laplacian whose norm is2; the energy is scalar.
 error=2*errors['q']+errors['parent_energy']+F(1,2*D)+F(1,64)*errors['sqrt3']
 source_error=4*max(errors[k]for k in('sqrt2c','sqrt6c','sqrt2d','sqrt6d'))
 return {'L':L,'M':M.tocsr(),'psi':psi,'dimension':count,'dyadic_denominator':D,'coefficient_intervals_error':{k:str(v)for k,v in errors.items()},'operator_error':error,'injection_error':source_error,'row_sum_upper':206}

def taylor_first(u,t):
 """Segmented inhomogeneous evolution with a proved binary64 bound."""
 from math import factorial
 t=F(t);steps=20 if t>F(1,100)else 2;dt=t/steps;degree=24;n=u['dimension'];v=np.zeros(n+1,complex);v[-1]=1
 K=sp.bmat([[u['M'],sp.csr_matrix((G*u['psi'])[:,None])],[None,sp.csr_matrix((1,1))]],format='csr');A=(-1j*float(dt)*K).tocsr()
 assert max(np.diff(A.indptr))<=17 and max(float(abs(A).sum(axis=1).max()),float(abs(A).sum(axis=0).max()))<2
 for _ in range(steps):
  term=v.copy();total=v.copy()
  for k in range(1,degree+1):term=(A@term)/k;total+=term
  v=total
 # ||exp(-i K s)||<=exp(G||psi||s)<2; polynomial intermediate sums <exp(2)<8.
 # gamma<=2^-44 covers each complex sparse row, division and addition.
 # At most24 stages give per-step roundoff <= 24*256*2^-44; amplification<2.
 roundoff=F(2*steps*degree*256,2**44)
 x=F(207)*dt;assert x<2
 truncation=2*steps*x**(degree+1)/F(factorial(degree+1))/(1-x/F(degree+2))
 source=t*(u['operator_error']+F(1,64)*u['injection_error'])*2
 error=roundoff+truncation+source
 assert error<F(1,10**7)
 return v[:n],{'segmented_Taylor_degree':degree,'segments':steps,'matrix_norm_upper':207,'binary64_gamma_upper':'1/17592186044416','roundoff_vector_error_upper_exact':str(roundoff),'Taylor_tail_vector_error_upper_exact':str(truncation),'source_surrogate_vector_error_upper_exact':str(source),'total_numeric_vector_error_upper_exact':str(error),'total_numeric_vector_error_upper':float(error)}
