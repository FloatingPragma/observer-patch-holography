"""Positive native colour, spin and weak-current source; exact rational sectors."""
from functools import lru_cache
from itertools import product
from collections import Counter
from math import factorial,isqrt,lcm
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from mpmath import iv
from model.openwave_native_composites.spectrum import sector,one_body

def integral(a):
 a=a.tocsr();assert np.max(abs(a.data.imag),initial=0)==0
 assert np.max(abs(a.data.real-np.rint(a.data.real)),initial=0)==0
 return sp.csr_matrix((np.rint(a.data.real).astype(np.int64),a.indices,a.indptr),shape=a.shape)

@lru_cache(None)
def balanced_baryon():
 old=sector('baryon');basis=old['basis'];modes=old['modes'];index={v:i for i,v in enumerate(modes)};n=len(basis);d=len(modes)
 zero=sp.csr_matrix((n,n),dtype=np.int64);Cnum=zero.copy();Snum=zero.copy();Inum=zero.copy();Iglobal=[]
 # Native Q is an actual fundamental colour module in the ordered exterior basis.
 for x in range(2):
  for a,b in ((0,1),(0,2),(1,2)):
   X=np.zeros((d,d));Y=X.copy()
   for ss,f in product(range(2),repeat=2):
    i=index[(0,x,ss,f,a)];j=index[(0,x,ss,f,b)];X[i,j]=X[j,i]=1;Y[i,j]=1;Y[j,i]=-1
   xx=integral(one_body(X,basis));yy=integral(one_body(Y,basis));Cnum+=12*(xx@xx-yy@yy)
  D1=[];D8=[]
  for mask in basis:
   cs=[sum(1 for i,(_,xx,ss,f,cc) in enumerate(modes) if xx==x and cc==c and mask>>i&1) for c in range(3)]
   D1.append(cs[0]-cs[1]);D8.append(cs[0]+cs[1]-2*cs[2])
  Cnum+=sp.diags(12*np.array(D1)**2+4*np.array(D8)**2,dtype=np.int64)
  for spin in (True,False):
   Z=[];A=np.zeros((d,d))
   for mask in basis:
    Z.append(sum((1 if (ss if spin else f)==0 else -1) for i,(_,xx,ss,f,c) in enumerate(modes) if xx==x and mask>>i&1))
   for other,c in product(range(2),range(3)):
    if spin:i=index[(0,x,0,other,c)];j=index[(0,x,1,other,c)]
    else:i=index[(0,x,other,0,c)];j=index[(0,x,other,1,c)]
    A[i,j]=1
   plus=integral(one_body(A,basis));Z=np.array(Z)
   term=sp.diags(3*(Z*Z+2*Z),dtype=np.int64)+12*(plus.T@plus)
   if spin:Snum+=term
   else:Inum+=term
 hop=np.zeros((d,d))
 for i,(_,x,ss,f,c) in enumerate(modes):hop[i,i]=1;hop[i,index[(0,1-x,ss,f,c)]]=-1
 Knum=6*integral(one_body(hop,basis))
 Hnum=144*sp.eye(n,dtype=np.int64)+Knum+Cnum+Snum+Inum
 gram=np.array([6*np.prod([factorial(v) for v in Counter(label).values()]) for label in old['labels']],dtype=np.int64)
 R=old['B']@sp.diags(np.sqrt(gram));R=sp.csr_matrix(R)
 assert np.max(abs(R.data.imag),initial=0)==0 and np.max(abs(R.data.real-np.rint(R.data.real)),initial=0)<1e-13
 R=sp.csr_matrix((np.rint(R.data.real).astype(np.int64),R.indices,R.indptr),shape=R.shape)
 assert (R.T@R-sp.diags(gram,dtype=float)).nnz==0
 A=(R.T@Hnum@R).toarray();L=lcm(*map(int,gram))
 assert (L*(Hnum@R)-R@sp.diags(L//gram,dtype=float)@sp.csr_matrix(A)).nnz==0
 h=A/(48*np.sqrt(np.outer(gram,gram)))
 B=R@sp.diags(1/np.sqrt(gram));assert np.max(abs((Hnum/48)@B-B@sp.csr_matrix(h)))<1e-12
 assert np.max(abs((old['H']*48+Inum-Hnum).data),initial=0)<1e-12
 return {'Hnum':Hnum,'denominator':48,'basis':basis,'modes':modes,'Binteger':R,'Gram':gram,'reduced_numerator':A,'h':h,'B':B,'old':old,'weak_term_numerator':Inum}

def ceil_sqrt(x):
 n=isqrt(int(x));return n if n*n==x else n+1

def frobenius_ceil(a):return ceil_sqrt(sum(int(x)*int(x) for x in np.asarray(a,dtype=object).ravel()))

@lru_cache(None)
def spectral_certificate():
 src=balanced_baryon();h=src['h'];A=src['reduced_numerator'];g=src['Gram'];n=len(g);D=2**40
 iv.dps=50
 H=np.rint(h*D).astype(np.int64);source_error=F(0)
 for i,j in zip(*np.nonzero(A)):
  z=iv.mpf(int(A[i,j]))/(48*iv.sqrt(int(g[i])*int(g[j])))
  lo=F.from_float(float(np.nextafter(float(z.a),-np.inf)));hi=F.from_float(float(np.nextafter(float(z.b),np.inf)))
  mid=F(int(H[i,j]),D);source_error=max(source_error,abs(lo-mid),abs(hi-mid))
 delta=n*source_error
 e,o=eigh(h);O=np.rint(o*D).astype(object);O=np.array([[int(v) for v in row] for row in O],dtype=object);E=np.array([int(round(x*D)) for x in e],dtype=object)
 gram_error=O.T@O-D*D*np.eye(n,dtype=object);eta=F(frobenius_ceil(gram_error),D*D);assert eta<F(1,2)
 residual=np.asarray(H,dtype=object)@O-O*E
 rho=F(frobenius_ceil(residual),D*D)+delta*(1+eta)
 bound=max(sum(abs(int(x)) for x in row) for row in H)
 K=(F(bound,D)+delta).__ceil__()+1
 assert max(abs(int(x)) for x in E)<K*D
 epsilon=2*rho+4*K*eta
 intervals=[[str(F(int(x),D)-epsilon),str(F(int(x),D)+epsilon)] for x in E]
 return {'source_dimension':n,'dyadic_denominator':D,'source_entry_to_dyadic_operator_error_upper':str(delta),
  'exact_Gram_Frobenius_error_upper':str(eta),'exact_residual_Frobenius_error_upper':str(rho),
  'all_eigenvalue_error_upper':str(epsilon),'source_and_eigenvalue_norm_bound':K,
  'lowest_interval_rational':intervals[0],'lowest_interval':[float(F(v)) for v in intervals[0]],
  'all_eigenvalue_intervals_rational':intervals,'dyadic_eigenvalues':[str(int(x)) for x in E],
  'dyadic_eigenframe':[[str(int(x)) for x in row] for row in O],
  'method':'Exact integer Gram/residual products plus outward scalar algebraic entry intervals; polar-unitary spectral comparison.'}

def selected(charge):
 a=balanced_baryon();old=a['old'];h=a['h'];q=old['q'];s=old['spin'];z=old['z'];n=len(h)
 filtermat=(q-charge*np.eye(n))@(q-charge*np.eye(n))+(s-.75*np.eye(n))@(s-.75*np.eye(n))+(z-.5*np.eye(n))@(z-.5*np.eye(n))
 e,v=eigh(filtermat);P=v[:,e<1e-8];w,u=eigh(P.T.conj()@h@P)
 gg=P@u[:,abs(w-w[0])<1e-9];proj=gg@gg.T.conj();diagonal=np.diag(proj).real;seed=int(np.flatnonzero(diagonal>max(diagonal)-1e-9)[0]);success=diagonal[seed];psi=proj[:,seed]/np.sqrt(success)
 full=a['B']@psi;probs=abs(full)**2;profile={};quadrupole=0.;coloc=0.
 for mask,p in zip(a['basis'],probs):
  occ=[v for i,v in enumerate(a['modes']) if mask>>i&1];centre=sum(v[1] for v in occ)/3
  if len({v[1] for v in occ})==1:coloc+=float(p)
  for _,x,ss,f,c in occ:
   charge_i=F(2,3) if f==0 else -F(1,3);distance=abs(F(x)-F(sum(v[1] for v in occ),3));key=str(distance)
   profile[key]=profile.get(key,0.)+float(p)*float(charge_i)
   quadrupole+=float(p)*2*float(charge_i)*float(distance)**2
 return {'charge':charge,'baseline_equal_mass_energy':round(float(w[0]),12),'internal_spin':.5,'co_location_probability':round(coloc,12),
  'intrinsic_charge_profile_by_number_centre_distance':{k:round(v,12) for k,v in sorted(profile.items())},'intrinsic_axial_quadrupole':round(quadrupole,12),
  'complete_unconditional_occupation_probabilities':[round(float(success*p),14) for p in probs],
  'preparation_failure_probability':round(float(1-success),14),'spectral_preparation_success_probability':round(float(success),14),
  'record_count':len(probs)+1,'preparation_seed':seed,'full_sector_eigen_residual':float(np.linalg.norm(h@psi-w[0]*psi))}

def integer_matvec(A,v):
 A=A.tocsr();return np.array([sum(int(A.data[k])*int(v[A.indices[k]]) for k in range(A.indptr[i],A.indptr[i+1])) for i in range(A.shape[0])],dtype=object)

def rational_sqrt_upper(q):
 if q==0:return F(0)
 return F(ceil_sqrt(q.numerator),max(1,isqrt(q.denominator)))

@lru_cache(None)
def exact_record_witness(charge):
 a=balanced_baryon();old=a['old'];h=a['h'];n=len(h);q=old['q'];s=old['spin'];z=old['z']
 f=(q-charge*np.eye(n))@(q-charge*np.eye(n))+(s-.75*np.eye(n))@(s-.75*np.eye(n))+(z-.5*np.eye(n))@(z-.5*np.eye(n))
 e,u=eigh(f);P=u[:,e<1e-8];ev,v=eigh(P.T.conj()@h@P);psi=P@v[:,0]
 assert np.max(abs(psi.imag))<1e-10;psi=psi.real
 # A rational colour-singlet trial in the original integer epsilon basis.
 coefficients=np.array([int(round(x/y**.5*2**45)) for x,y in zip(psi,a['Gram'])],dtype=object)
 for i,label in enumerate(old['labels']):
  nu=sum(t[2]==0 for t in label);nup=sum(t[1]==0 for t in label)
  if nu-1!=charge or nup!=2:coefficients[i]=0
 vfull=integer_matvec(a['Binteger'],coefficients)
 # Exact spin-one-half projector, with its harmless scalar denominator omitted.
 modes=a['modes'];index={x:i for i,x in enumerate(modes)};raising=np.zeros((24,24))
 for x,fl,c in product(range(2),range(2),range(3)):raising[index[(0,x,0,fl,c)],index[(0,x,1,fl,c)]]=1
 Rspin=integral(one_body(raising,a['basis']));S4=3*sp.eye(len(vfull),dtype=np.int64)+4*(Rspin.T@Rspin)
 vfull=15*vfull-integer_matvec(S4,vfull)
 assert any(vfull)
 norm=sum(int(x)**2 for x in vfull);Hv=integer_matvec(a['Hnum'],vfull);energy=F(sum(int(x)*int(y) for x,y in zip(vfull,Hv)),48*norm)
 residual=[energy.denominator*int(y)-48*energy.numerator*int(x) for x,y in zip(vfull,Hv)]
 rho2=F(sum(x*x for x in residual),48**2*energy.denominator**2*norm);rho=rational_sqrt_upper(rho2)
 spec=spectral_certificate();intervals=[[F(z) for z in row] for row in spec['all_eigenvalue_intervals_rational']]
 assert max(x[1] for x in intervals[:4])<intervals[4][0]
 gamma=intervals[4][0]-energy;assert gamma>0
 trace=2*rho/gamma
 # The isolated four-dimensional band is one spin-half/isospin-half irrep.
 profile={};quad=F(0)
 for mask,amp in zip(a['basis'],vfull):
  if amp==0:continue
  prob=F(int(amp)**2,norm);occ=[v for i,v in enumerate(modes) if mask>>i&1];centre=F(sum(v[1] for v in occ),3)
  assert sum(1 for v in occ if v[3]==0)-1==charge
  for _,x,ss,fl,c in occ:
   qi=F(2,3) if fl==0 else -F(1,3);r=abs(F(x)-centre);profile[r]=profile.get(r,F(0))+prob*qi;quad+=2*prob*qi*r*r
 # Exact rank-one filter on an explicit original colour-singlet epsilon seed.
 overlaps=integer_matvec(a['Binteger'].T,vfull)
 candidates=[F(int(x)**2,norm*int(g)) for x,g in zip(overlaps,a['Gram'])];seed=max(range(len(candidates)),key=lambda i:candidates[i]);success=candidates[seed]
 assert 0<success<=1
 rows={str(r):{'trial_exact':str(p),'ground_interval':[float(p-2*trace),float(p+2*trace)]} for r,p in sorted(profile.items())}
 if charge==0:
  assert F(1,3) in profile and profile[F(1,3)]>2*trace
  assert profile[F(2,3)]< -2*trace
 return {'charge':charge,'exact_integer_trial_vector':[str(int(x)) for x in vfull],'squared_integer_norm':str(norm),
  'trial_Rayleigh_exact':str(energy),'trial_Rayleigh':float(energy),'source_residual_norm_upper':str(rho),'next_band_distance_lower':str(gamma),
  'trace_distance_to_exact_charge_spin_ground_line_upper':str(trace),'ground_band_dimension':4,
  'intrinsic_charge_profile':rows,'intrinsic_quadrupole_trial_exact':str(quad),
  'exact_preparation_seed_epsilon_column':seed,'exact_filter_success_probability':str(success),'exact_filter_failure_probability':str(1-success),
  'complete_unconditional_native_occupation_probabilities':[str(success*F(int(x)**2,norm)) for x in vfull],
  'record_count':len(vfull)+1,'clock':'t=0 complete records; at arbitrary t the full-record trace distance from these records is at most min(2,2*abs(t)*rho)',
  'failure_quantum_output':'(I-|v><v|/norm) times the normalized epsilon seed, retained with its probability; idle opposite-charge reference unchanged'}
