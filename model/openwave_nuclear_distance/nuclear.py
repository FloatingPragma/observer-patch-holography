"""Exact native Casimir/parity census with full six-quark spectra."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import gcd
import numpy as np
import scipy.sparse as sp
import sympy as sy
from scipy.linalg import eigh,qr,cholesky,solve_triangular
from scipy.optimize import brentq
from model.openwave_six_quark import source as n
from model.openwave_axial_state.source import parameters,weights,fit_at_colour,matrix,INPUTS
from model.openwave_baryon_weak_formfactors.native import formfactors

@lru_cache(None)
def carrier():
 src=n.basis(3,3);N=len(src);counts=[sum(1<<i for i,v in enumerate(n.MODES)if v[3]==c)for c in range(3)];keep=np.array([i for i,b in enumerate(src)if all((b&m).bit_count()==2 for m in counts)])
 def restrict(A):return A.tocsr()[keep][:,keep].tocsr()
 Cg=n.colour48(src,(0,1));assert np.all(Cg.data%16==0);C3=restrict(Cg);C3.data//=16
 Sp=n.raising(src,(0,1),True,3,3);Ip=n.raising(src,(0,1),False,3,3);Sg=restrict(Sp.T@Sp);Ig=restrict(Ip.T@Ip)
 C=n.colour48(src,(0,))+n.colour48(src,(1,));S=sp.csr_matrix((N,N),dtype=np.int64);I=S.copy()
 for cell in(0,1):
  for spin in(True,False):
   z=n.local_number(src,cell,1 if spin else 2);R=n.raising(src,(cell,),spin,3,3);term=sp.diags(z*z+2*z,dtype=np.int64)+4*(R.T@R)
   if spin:S+=term
   else:I+=term
 hop=[t for i,(x,s,f,c)in enumerate(n.MODES)for t in((i,i,1),(i,n.INDEX[1-x,s,f,c],-1))];K=n.op(hop,src,src);n0=n.local_number(src,0)
 parts={'N':6*sp.eye(N,dtype=np.int64,format='csr'),'D':sp.csr_matrix((N,N),dtype=np.int64),'K':K,'C48':C,'Cg48':Cg,'S4':S,'I4':I,'B9':sp.diags(n0*n0+(6-n0)**2,dtype=np.int64),'Ig4':4*(Ip.T@Ip)}
 lookup={b:i for i,b in enumerate(src)};rr=[];vv=[]
 for b in src:
  occupied=[i for i in range(24)if b>>i&1];mapped=[n.INDEX[1-n.MODES[i][0],*n.MODES[i][1:]]for i in occupied];target=sum(1<<i for i in mapped);sign=(-1)**sum(mapped[i]>mapped[j]for i in range(6)for j in range(i+1,6));rr.append(lookup[target]);vv.append(sign)
 parity=restrict(sp.csr_matrix((vv,(rr,range(N))),shape=(N,N),dtype=np.int64));eye=sp.eye(len(keep),dtype=np.int64,format='csr');P=eye.copy();den=1
 for c in(9,18,24,36,54):
  P=n.exact_product(c*eye-C3,P);den*=c;g=gcd(den,*map(int,P.data));P.data//=g;den//=g
 assert(n.exact_product(P,P)-den*P).nnz==0 and n.exact_product(C3,P).nnz==0
 p={k:restrict(a)for k,a in parts.items()}
 for a in p.values():assert(n.exact_product(a,P)-n.exact_product(P,a)).nnz==0
 return {'src':src,'keep':keep,'P':P,'den':den,'Sg':Sg,'Ig':Ig,'parity':parity,'parts':p,'full_dimension':N,'colour_weight_dimension':int(P.diagonal().sum())//den,'n0':n0[keep]}

@lru_cache(None)
def channel(S,I,parity):
 d=carrier();eye=sp.eye(len(d['keep']),dtype=np.int64,format='csr');P=d['P'].copy();den=d['den']
 for op,j in((d['Sg'],S),(d['Ig'],I)):
  for c in(0,2,6,12):
   if c==j*(j+1):continue
   P=n.exact_product(op-c*eye,P);den*=j*(j+1)-c;g=gcd(abs(den),*map(int,P.data));P.data//=g;den//=g
 P=n.exact_product(eye+parity*d['parity'],P);den*=2;g=gcd(abs(den),*map(int,P.data));P.data//=g;den//=g
 if den<0:P=-P;den=-den
 assert(P-P.T).nnz==0 and(n.exact_product(P,P)-den*P).nnz==0
 rank=int(P.diagonal().sum())//den
 if not rank:return {'rank':0}
 assert (n.exact_product(d['Sg'],P)-S*(S+1)*P).nnz==0
 assert(n.exact_product(d['Ig'],P)-I*(I+1)*P).nnz==0
 rng=np.random.default_rng(100*S+10*I+parity+2);Y=P@rng.normal(size=(len(d['keep']),rank));_,_,piv=qr(Y.T,pivoting=True,mode='economic');B=P[:,piv[:rank]].tocsr();G=n.exact_product(B.T,B).toarray();Gram=sy.Matrix(G.tolist());_,D=Gram.LDLdecomposition(hermitian=True);assert all(x>0 for x in D.diagonal())
 # Exact rank and invariance certify the subspace; floating orthonormalization is only used for reported spectral diagnostics.
 Q,_=np.linalg.qr(B.toarray().astype(float));parts={k:Q.T@(a@Q)for k,a in d['parts'].items()}
 assert np.linalg.norm(Q.T@Q-np.eye(rank))<2e-10
 return {'rank':rank,'projector_denominator':den,'parts':parts,'Q':Q,'Binteger':B,'Gram':Gram,'spin':S,'isospin':I,'parity':parity}

def spectrum(p,S,I,parity):
 d=channel(S,I,parity)
 if not d['rank']:return []
 H=sum(float(v)*d['parts'][k]for k,v in weights(p).items());return eigh(H,eigvals_only=True).tolist()

@lru_cache(None)
def fitted_eta(eta):
 p0=parameters();D=float(8*p0['U']/3-4*p0['A']/9);constant=float(3*p0['mbar']+p0['A']+3*p0['eta']/4)
 from model.openwave_hadron_spin_constraints.constraints import energy
 def par(A):return {**p0,'eta':eta,'A':A,'U':3*(D+4*A/9)/8,'mbar':(constant-A-.75*eta)/3}
 target=float(INPUTS['mass_calibration_targets_GeV']['charged_pion'])/float(INPUTS['electron_unit_GeV'])
 A=brentq(lambda x:energy('meson',par(x),1,0)-target,1000.,2000.,xtol=1e-10)
 p=par(A);assert abs(formfactors(matrix('baryon',p))['gA']-1.2753)<2e-8
 return p

@lru_cache(None)
def exact_channel(S,I,parity):
 from model.openwave_six_quark.certificate import spectral_enclosure
 d=carrier();c=channel(S,I,parity);B=c['Binteger'];Gram=c['Gram'];L,D=Gram.LDLdecomposition(hermitian=True);Li=L.inv();diag=list(D.diagonal());A=sy.zeros(c['rank'])
 for k,v in weights(parameters()).items():A+=sy.Rational(v.numerator,v.denominator)*sy.Matrix(n.exact_product(B.T,n.exact_product(d['parts'][k],B)).toarray().tolist())
 T=Li*A*Li.T;H=sy.Matrix(c['rank'],c['rank'],lambda i,j:T[i,j]/sy.sqrt(diag[i]*diag[j]));return spectral_enclosure(H)

@lru_cache(None)
def nuclear_certificate():
 d=carrier();p=parameters();rows=[];total=0;magnetic_total=0
 for S,I,par in product(range(4),range(4),(-1,1)):
  c=channel(S,I,par);rank=c['rank'];total+=rank;magnetic_total+=rank*(2*S+1)*(2*I+1)
  if rank:rows.append({'spin':S,'isospin':I,'cell_parity':par,'fixed_Ms0_Mi0_dimension':rank,'all_eigenvalues':spectrum(p,S,I,par),'projector_denominator':c['projector_denominator']})
 assert total==d['colour_weight_dimension'];assert magnetic_total==2520
 def bindings(par):
  e=eigh(matrix('baryon',par),eigvals_only=True);# fixed p+n uses charge-resolved minima, supplied through exact delta symmetry.
  from model.openwave_hadron_spin_constraints.constraints import energy
  threshold=energy('baryon',par,0,1)+energy('baryon',par,1,1);out=[]
  for S,I in((1,0),(0,0),(0,1)):
   for parity in(-1,1):
    es=spectrum(par,S,I,parity)
    if es:out.append({'spin':S,'isospin':I,'cell_parity':parity,'sector_ground':es[0],'two_baryon_threshold':threshold,'binding':threshold-es[0],'binding_GeV':(threshold-es[0])*float(INPUTS['electron_unit_GeV'])})
  return out
 controls=[]
 for eta in(1.,2.,4.,6.):
  pp=p if eta==1 else fitted_eta(eta);controls.append({'eta':eta,'parameters':{k:float(v)for k,v in pp.items()},'ancestry':'Four empirical inputs refitted for each eta; eta is a disclosed control choice.','channels':bindings(pp)})
 target=float(INPUTS['inherited_mass_input_ancestry']['held_out_comparisons_GeV']['deuteron_binding'])
 def residual(eta):return max(x['binding_GeV']for x in bindings(fitted_eta(eta))if(x['spin'],x['isospin'])==(1,0))-target
 eta_fit=brentq(residual,4.,6.,xtol=1e-8);fitted=fitted_eta(eta_fit);fit={'eta':eta_fit,'parameters':{k:float(v)for k,v in fitted.items()},'fifth_empirical_input_GeV':target,'channels':bindings(fitted)}
 return {'full_N6_Ms0_Mi0_occupation_dimension':len(d['src']),'zero_Cartan_dimension':len(d['keep']),'exact_complete_colour_singlet_weight_dimension':total,'complete_colour_singlet_dimension_including_all_magnetic_components':magnetic_total,'exact_ground_channel_certificates':{'S1_I0_even':exact_channel(1,0,1),'S0_I0_odd':exact_channel(0,0,-1)},'all_spin_isospin_parity_channels':rows,'fixed_four_input_source_binding_comparisons':bindings(p),'eta_controls':controls,'deuteron_fit_control':fit,
  'parity_scope':'Exact cell-swap parity of all six native fermions. It is not identified with physical continuum orbital parity. No channel is discarded.',
  'calibration_scope':'The eta fit consumes deuteron binding as a fifth empirical input and is a separate source. It does not change the fixed eta1 source used by the common beta join.',
  'arithmetic':'Exact integer Casimir/parity ranks and source intertwiners, plus exact outward source certificates for S1I0-even and S0I0-odd. Other full level values and eta binding controls use floating diagonalization.'}
