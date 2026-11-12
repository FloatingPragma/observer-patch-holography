"""Integer current Casimirs and exact27-dimensional invariant reduction."""
from itertools import product,combinations,permutations
from collections import Counter
from functools import lru_cache
from math import gcd
import numpy as np
import sympy as sy
import scipy.sparse as sp
from scipy.linalg import qr
from model.coherent_defect_wave.physical_attachment.chiral import representation

MODES=list(product(range(2),range(2),range(2),range(3)))
INDEX={v:i for i,v in enumerate(MODES)}


@lru_cache(None)
def native_Fierz():
 subsets,_,g,_=representation();ids=[[subsets.index((c,w))for c in range(3)]for w in(3,4)]
 generators=list(g[:7])+[(g[6]+2*g[7])/sy.sqrt(3)]
 T=[[a.extract(ii,ii)for a in generators]for ii in ids]
 assert T[0]==T[1]
 for i,j,k,l in product(range(3),repeat=4):
  value=sum(a[i,j]*a[k,l]for a in T[0]);expected=sy.Rational(int(i==l and j==k),2)-sy.Rational(int(i==j and k==l),6)
  assert sy.simplify(value-expected)==0
 charges=[[(g[10]+g[11])[i,i]for i in ii]for ii in ids];assert charges==[[sy.Rational(2,3)]*3,[-sy.Rational(1,3)]*3]
 return {'exact_actual_native_Fierz_entries':81,'same_colour_frame_on_both_Q_components':True,'native_Q_charges':['2/3','-1/3']}

@lru_cache(None)
def basis(nu=3,nup=4):
 groups=[[i for i,v in enumerate(MODES)if v[1]==s and v[2]==f]for s,f in product(range(2),repeat=2)];out=[]
 for a in range(7):
  counts=[a,nup-a,nu-a,6-nu-nup+a]
  if any(n<0 or n>6 for n in counts):continue
  masks=[[sum(1<<i for i in ids)for ids in combinations(g,n)]for g,n in zip(groups,counts)]
  out.extend(x|y|z|w for x,y,z,w in product(*masks))
 return tuple(sorted(out))

def op(terms,src,dst):
 lookup={b:i for i,b in enumerate(dst)};rr=[];cc=[];vv=[]
 for col,b in enumerate(src):
  for i,j,v in terms:
   if not(b>>j)&1:continue
   b1=b^(1<<j)
   if(b1>>i)&1:continue
   target=b1|(1<<i);row=lookup.get(target)
   if row is None:raise ValueError('declared occupation sector is not invariant')
   sign=(-1)**((b&((1<<j)-1)).bit_count()+(b1&((1<<i)-1)).bit_count())
   rr.append(row);cc.append(col);vv.append(sign*v)
 return sp.csr_matrix((vv,(rr,cc)),shape=(len(dst),len(src)),dtype=np.int64)

def local_number(src,x,axis=None):
 return np.array([sum((1 if axis is None or v[axis]==0 else -1)for i,v in enumerate(MODES)if(x is None or v[0]==x)and(mask>>i)&1)for mask in src],dtype=np.int64)

def colour48(src,cells):
 n=len(src);C=sp.csr_matrix((n,n),dtype=np.int64)
 for a,b in((0,1),(0,2),(1,2)):
  X=[];Y=[]
  for x,ss,f in product(cells,range(2),range(2)):
   i=INDEX[x,ss,f,a];j=INDEX[x,ss,f,b];X.extend(((i,j,1),(j,i,1)));Y.extend(((i,j,1),(j,i,-1)))
  xx=op(X,src,src);yy=op(Y,src,src);C+=12*(xx@xx-yy@yy)
 D1=[];D8=[]
 for mask in src:
  cs=[sum(1 for i,v in enumerate(MODES)if v[0]in cells and v[3]==c and(mask>>i)&1)for c in range(3)];D1.append(cs[0]-cs[1]);D8.append(cs[0]+cs[1]-2*cs[2])
 return C+sp.diags(12*np.array(D1,dtype=np.int64)**2+4*np.array(D8,dtype=np.int64)**2,dtype=np.int64)

def raising(src,cells,spin,nu=3,nup=4):
 terms=[]
 for x,other,c in product(cells,range(2),range(3)):
  a=(x,0,other,c)if spin else(x,other,0,c);b=(x,1,other,c)if spin else(x,other,1,c);terms.append((INDEX[a],INDEX[b],1))
 return op(terms,src,basis(nu,nup+1)if spin else basis(nu+1,nup))

@lru_cache(None)
def full_source():
 native_Fierz();src=basis();n=len(src);C=colour48(src,(0,))+colour48(src,(1,));S=sp.csr_matrix((n,n),dtype=np.int64);I=S.copy()
 for x in range(2):
  for spin in(True,False):
   z=local_number(src,x,1 if spin else 2);plus=raising(src,(x,),spin);term=sp.diags(3*(z*z+2*z),dtype=np.int64)+12*(plus.T@plus)
   if spin:S+=term
   else:I+=term
 hop=[]
 for i,(x,ss,f,c)in enumerate(MODES):hop.extend(((i,i,1),(i,INDEX[1-x,ss,f,c],-1)))
 H=288*sp.eye(n,dtype=np.int64)+6*op(hop,src,src)+C+S+I
 Cglobal=colour48(src,(0,1)).tocsr();assert np.all(Cglobal.data%16==0);C3=Cglobal.copy();C3.data//=16
 ps=raising(src,(0,1),True);Sg=2*sp.eye(n,dtype=np.int64)+ps.T@ps;pi=raising(src,(0,1),False);Ig=pi.T@pi
 colour_masks=[sum(1<<i for i,v in enumerate(MODES)if v[3]==c)for c in range(3)]
 keep=np.array([i for i,b in enumerate(src)if all((b&m).bit_count()==2 for m in colour_masks)],dtype=int)
 def restrict(A):
  sub=A[keep,:].tocsr();outside=np.ones(n,bool);outside[keep]=False;assert sub[:,outside].nnz==0
  return sub[:,keep].tocsr()
 return {'basis':src,'keep':keep,'Hnum':restrict(H),'unbalanced_Hnum':restrict(H-I),'C3':restrict(C3),'Sg':restrict(Sg),'Ig':restrict(Ig),
  'full_dimension':n,'zero_Cartan_dimension':len(keep),'full_H_trace_numerator':int(H.diagonal().sum())}

def exact_product(A,B):
 A=A.tocsr();B=B.tocsr()
 bound=int(np.max(abs(A.data),initial=0))*int(np.max(abs(B.data),initial=0))*max(1,int(np.max(np.diff(A.indptr))))
 assert bound<2**62,'integer sparse product overflow contract'
 out=(A@B).tocsr();out.eliminate_zeros();return out

@lru_cache(None)
def reduction():
 d=full_source();n=len(d['keep']);eye=sp.eye(n,dtype=np.int64,format='csr');P=eye.copy();den=1;steps=[]
 factors=[(c*eye-d['C3'],c)for c in(9,18,24,36,54)]+[(c*eye-d['Sg'],c-2)for c in(6,12)]+[(c*eye-d['Ig'],c)for c in(2,6,12)]
 for A,b in factors:
  P=exact_product(A.tocsr(),P);den*=b;common=gcd(den,*map(int,P.data));P.data//=common;den//=common;steps.append({'denominator':den,'nonzero_entries':P.nnz})
 assert (P-P.T).nnz==0
 assert (exact_product(P,P)-den*P).nnz==0
 assert exact_product(d['C3'],P).nnz==0 and(exact_product(d['Sg'],P)-2*P).nnz==0 and exact_product(d['Ig'],P).nnz==0
 assert (exact_product(d['Hnum'],P)-exact_product(P,d['Hnum'])).nnz==0
 assert int(P.diagonal().sum())==27*den
 rng=np.random.default_rng(7);Y=P@rng.normal(size=(n,27));_,_,pivot=qr(Y.T,pivoting=True,mode='economic');B=P[:,pivot[:27]].tocsr();G=exact_product(B.T,B).toarray();Gram=sy.Matrix(G.tolist());L,D=Gram.LDLdecomposition(hermitian=True);assert all(x>0 for x in D.diagonal())
 A=exact_product(B.T,exact_product(d['Hnum'],B)).toarray();A0=exact_product(B.T,exact_product(d['unbalanced_Hnum'],B)).toarray();Li=L.inv();T=Li*sy.Matrix(A.tolist())*Li.T;T0=Li*sy.Matrix(A0.tolist())*Li.T;diag=list(D.diagonal())
 H=sy.Matrix(27,27,lambda i,j:T[i,j]/(48*sy.sqrt(diag[i]*diag[j])));H0=sy.Matrix(27,27,lambda i,j:T0[i,j]/(48*sy.sqrt(diag[i]*diag[j])))
 transform=np.asarray(Li.T.evalf(18),float)@np.diag([1/np.sqrt(float(x))for x in diag]);h=np.asarray(H.evalf(18),float);h0=np.asarray(H0.evalf(18),float)
 return {**d,'Pnum':P,'Pden':den,'projector_steps':steps,'Binteger':B,'basis_columns':pivot[:27].tolist(),'Gram':Gram,'LDL_diagonal':diag,'H_exact':H,'H0_exact':H0,'h':h,'h0':h0,'coefficient_transform':transform}
