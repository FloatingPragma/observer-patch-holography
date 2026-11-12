"""Charge sectors of the same native contact on two stated strong actions."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations,product
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_axial_state import source as axial
from model.openwave_calibrated_hadrons import source as mass_only
from model.openwave_baryon_weak_formfactors.native import currents,quartet,PAULI
from model.openwave_native_composites.spectrum import one_body
from model.common_source import calibrated_matter as internal
G=F(1,64)
def lepton_index(w,x,s):return 4*w+2*x+s

@lru_cache(None)
def common():
 c=currents();masks=[sum(1<<i for i in v)for v in combinations(range(8),4)];reference=masks.index(15)
 ne=np.array([(m>>4).bit_count()for m in masks]);lz=np.array([sum(1-2*s for w,x,s in product(range(2),repeat=3)if m>>lepton_index(w,x,s)&1)for m in masks])
 q=(mass_only.parameters()['delta']-1)/2;K=np.array([[1,-1],[-1,1]],dtype=int);a=np.zeros((8,8),int)
 for w,s in product(range(2),repeat=2):
  ids=[lepton_index(w,x,s)for x in range(2)];a[np.ix_(ids,ids)]=(1 if w else -1)*K
 K2=one_body(a,masks)+4*sp.eye(70);Hl=sp.diags(ne,dtype=np.int64)+float(q)/2*K2
 raising=np.array(internal.weak_raising_current(families=False)).astype(complex)[np.ix_([10,11],[10,11])];third=(raising@raising.conj().T-raising.conj().T@raising)/2
 V=sp.csr_matrix((8400,8400),dtype=complex)
 for x,mu in product(range(2),range(4)):
  lp,lzop=[] ,[]
  operators=[]
  for weak in(raising,third):
   A=np.zeros((8,8),complex)
   for w,v,s,z in product(range(2),repeat=4):A[lepton_index(w,x,s),lepton_index(v,x,z)]=PAULI[mu][s,z]*weak[w,v]
   operators.append(one_body(A,masks))
  lp,lzop=operators;qp,qz=map(sp.csr_matrix,c['local'][x][mu]);V+=(1 if mu==0 else -1)*(sp.kron(qp,lp.conj().T)+sp.kron(qp.conj().T,lp)+2*sp.kron(qz,lzop))
 # Exact site exchange includes all fermionic exterior-power signs.
 lookup={m:i for i,m in enumerate(masks)};P=np.zeros((70,70),int)
 for j,m in enumerate(masks):
  mapped=[lepton_index(w,1-x,s)for w,x,s in product(range(2),repeat=3)if m>>lepton_index(w,x,s)&1]
  P[lookup[sum(1<<i for i in mapped)],j]=(-1)**sum(mapped[a]>mapped[b]for a in range(4)for b in range(a+1,4))
 return dict(currents=c,masks=masks,reference=reference,ne=ne,lz=lz,q=q,K2=K2,Hl=Hl,V=V.tocsr(),parity=sp.kron(sp.csr_matrix(c['parity']),sp.csr_matrix(P)).tocsr())

@lru_cache(None)
def construction(member='axial',charge=0,coupling=G):
 provider=axial if member=='axial'else mass_only;native=provider.sector('baryon');c=common();ne=c['ne'];qh=np.rint(np.diag(native['q']).real).astype(int);zh=np.rint(2*np.diag(native['z']).real).astype(int)
 charges=(qh[:,None]-ne).ravel();spins=(zh[:,None]+c['lz']).ravel();keep=np.flatnonzero((charges==charge)&(spins==1));outside=np.flatnonzero((charges!=charge)|(spins!=1))
 full0=sp.kron(sp.csr_matrix(native['h']),sp.eye(70))+sp.kron(sp.eye(120),c['Hl']);H0=full0.tocsr()[keep][:,keep];V=c['V'][keep][:,keep];H=H0+float(coupling)*V
 assert np.max(abs(full0.tocsr()[keep][:,outside].data),initial=0)<3e-10 and np.max(abs(c['V'][keep][:,outside].data),initial=0)<3e-12
 parity=c['parity'][keep][:,keep];assert np.max(abs((H@parity-parity@H).data),initial=0)<3e-10
 W,_=quartet(native['h']);input_full=np.kron(W[:,0 if charge==0 else 2],np.eye(70)[:,c['reference']]);initial=input_full[keep]
 assert np.linalg.norm(parity@initial-initial)<1e-10 and abs(np.linalg.norm(initial)-1)<1e-10
 return dict(c,member=member,provider=provider,native=native,coupling=F(coupling),charge=charge,keep=keep,H0=H0,V=V,H=H,parity=parity,initial=initial,charge_dimension=int(np.sum(charges==charge)),spin_dimension=len(keep))

def parity_basis(d,sign=1):
 P=d['parity'].tocsr();n=P.shape[0];rows=[];cols=[];data=[];count=0
 for i in range(n):
  row=P.getrow(i);assert row.nnz==1;j=int(row.indices[0]);s=int(round(row.data[0].real));assert abs(s)==1
  if j<i:continue
  if j==i:
   if s==sign:rows.append(i);cols.append(count);data.append(1.);count+=1
  else:
   rows.extend((i,j));cols.extend((count,count));data.extend((1/np.sqrt(2),sign*s/np.sqrt(2)));count+=1
 B=sp.coo_matrix((data,(rows,cols)),shape=(n,count)).tocsr();assert np.max(abs((B.T@B-sp.eye(count)).data),initial=0)<1e-14
 return B

def exact_source(d):
 native=d['native'];parts=native['integer_reduced_parts'];weights=d['provider'].weights(d['provider'].parameters());g=native['Gram'];n=len(g);hadron=sy.zeros(n)
 for name,w in weights.items():hadron+=sy.Rational(w.numerator,w.denominator)*sy.Matrix(parts[name].tolist())
 keep=d['keep'];gram=np.repeat(g,70)[keep];scaled=d['V'].multiply(np.sqrt(gram)[:,None]).multiply(np.sqrt(gram)[None,:]).tocoo();vdict={}
 for i,j,z in zip(scaled.row,scaled.col,scaled.data):
  assert abs(z.imag)<1e-10 and abs(2*z.real-round(2*z.real))<1e-9;vdict[int(i),int(j)]=sy.Rational(int(round(2*z.real)),2)
 K2=d['K2'].toarray();assert np.max(abs(K2.imag))<1e-12 and np.max(abs(K2.real-np.rint(K2.real)))<1e-12
 q=sy.Rational(d['q'].numerator,d['q'].denominator);Gq=sy.Rational(d['coupling'].numerator,d['coupling'].denominator);H=sy.MutableSparseMatrix(len(keep),len(keep),{})
 for i,original_i in enumerate(keep):
  hi,li=divmod(int(original_i),70)
  for j,original_j in enumerate(keep):
   hj,lj=divmod(int(original_j),70);value=hadron[hi,hj]/sy.sqrt(int(g[hi])*int(g[hj]))if li==lj else sy.Integer(0)
   if hi==hj:value+=(int(d['ne'][li])if li==lj else 0)+q*int(round(K2[li,lj].real))/2
   if(i,j)in vdict:value+=Gq*vdict[i,j]/sy.sqrt(int(gram[i])*int(gram[j]))
   if value!=0:H[i,j]=value
 assert np.max(abs(np.asarray(H.evalf(17),float)-d['H'].toarray().real))<3e-10
 return H

def exact_parity_matrix(d,sign=1):
 H=exact_source(d);B=parity_basis(d,sign).tocoo();S=sy.MutableSparseMatrix(B.shape[0],B.shape[1],{})
 for i,j,v in zip(B.row,B.col,B.data):S[int(i),int(j)]=sy.sign(v)*(1 if abs(v)>.9 else 1/sy.sqrt(2))
 return S.T*H*S,S
