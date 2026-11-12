"""Padded exact oscillator moments and enclosed infinite-field Ritz witnesses."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
from math import fsum,sqrt,comb
from pathlib import Path
import json
import numpy as np
import sympy as s
from scipy.sparse import coo_matrix,bmat,csr_matrix
from mpmath import iv
from model.openwave_quantum_rotational.source import matter_exact,occupations
from model.openwave_quantum_rotational.enclosure import scalar

MAT_EPS=2.**-40
DOT_EPS=2.**-48
TINY=1e-300

def upper(x):return float(np.nextafter(float(x),np.inf))
def lower(x):return float(np.nextafter(float(x),-np.inf))
def isum(values):return fsum(float(x)for x in values)
def dot_enclosure(v,r,w,z):
 products=v*w;mid=isum(products)
 scale=isum(abs(products));error=isum(abs(v)*z+abs(w)*r+r*z)
 return mid,upper((error+DOT_EPS*scale)*(1+16*DOT_EPS)+TINY)
def norm_upper(v,r):
 mid,error=dot_enclosure(v,r,v,r);return upper(sqrt(upper(mid+error)))

def plus(v,r,w,z,coefficient=1.):
 out=v+coefficient*w
 radius=(r+abs(coefficient)*z+DOT_EPS*(abs(v)+abs(coefficient)*abs(w)))*(1+16*DOT_EPS)+TINY
 return out,np.nextafter(radius,np.inf)

def multiply(M,R,v,r):
 # Real CSR arithmetic with a declared IEEE-754 guard larger than gamma_n.
 out=M@v
 major=abs(M)@abs(v)
 radius=(abs(M)@r+R@(abs(v)+r)+MAT_EPS*major)*(1+16*MAT_EPS)+TINY
 return out,np.nextafter(radius,np.inf)

@lru_cache(maxsize=8)
def padded_operator(cap):
 iv.dps=45;scalar.cache_clear();m=matter_exact();fb=occupations(cap);where={p:i for i,p in enumerate(fb)};d=len(fb);n=20*d
 assert 2*n<100000
 rows=[];cols=[];values=[];radii=[]
 cache={}
 def add(i,j,expr):
  expr=s.expand(expr)
  if expr not in cache:
   pair=[]
   for part in(s.re(expr),s.im(expr)):
    if part==0:
     pair.append((0.,0.));continue
    z=scalar(part);lo=lower(float(z.a));hi=upper(float(z.b));mid=float((lo+hi)/2);rad=upper(max(mid-lo,hi-mid))
    pair.append((mid,rad))
   cache[expr]=pair
  (re,er),(im,ei)=cache[expr]
  for ri,cj,value,radius in((i,j,re,er),(i+n,j+n,re,er),(i+n,j,im,ei),(i,j+n,-im,ei)):
   if value or radius:rows.append(ri);cols.append(cj);values.append(value);radii.append(radius)
 for (i,j),z in m['B'].todok().items():
  for k in range(d):add(i*d+k,j*d+k,z)
 for axis,mat in enumerate(m['derivative']):
  entries=list(mat.todok().items())
  for k,p in enumerate(fb):
   for direction in(-1,1):
    q=list(p);q[axis]+=direction
    if q[axis]<0 or tuple(q)not in where:continue
    ell=where[tuple(q)];coefficient=s.sqrt(max(p[axis],q[axis]))/8
    for (i,j),z in entries:add(i*d+ell,j*d+k,z*coefficient)
 M=coo_matrix((values,(rows,cols)),shape=(2*n,2*n)).tocsr();R=coo_matrix((radii,(rows,cols)),shape=(2*n,2*n)).tocsr()
 assert M.nnz==len(values), 'exact matrix assembly has duplicate entries'
 maxrow=int(np.max(np.diff(M.indptr)))
 assert maxrow<256 and np.max(abs((M-M.T).data),initial=0)<1e-14
 # The enclosure uses fewer than 256 terms per row, comfortably below its guard.
 assert (maxrow+4)*2.**-53/(1-(maxrow+4)*2.**-53)<MAT_EPS/8
 number=np.array([sum(m['basis'][(i//d)%10])+sum(fb[i%d])for i in range(n)],float)
 return M,R,number,fb,maxrow


def trial(cap,padded_cap):
 old=json.loads((Path(__file__).resolve().parents[1]/'openwave_quantum_rotational/results.json').read_text())
 row=next(x for x in old['interval_certificates']if x['cap']==cap)
 integers=row['trial_integer_coefficients'];D=sum(a*a+b*b for a,b in integers)
 assert str(D)==row['trial_squared_integer_norm']
 small=occupations(cap);large=occupations(padded_cap);where={p:i for i,p in enumerate(large)};n=20*len(large);values=np.zeros(2*n);radii=np.zeros(2*n);scale=iv.sqrt(iv.mpf(D))
 for i,(a,b)in enumerate(integers):
  j=(i//len(small))*len(large)+where[small[i%len(small)]]
  for address,k in((j,a),(j+n,b)):
   exact=iv.mpf(k)/scale;lo=lower(float(exact.a));hi=upper(float(exact.b));value=float((lo+hi)/2)
   values[address]=value;radii[address]=upper(max(value-lo,hi-value))
 return values,radii,integers,D


def certify(cap):
 padded=cap+14;M,R,number,basis,maxrow=padded_operator(padded);v,r,integers,D=trial(cap,padded)
 powers=[(v,r)]
 for k in range(14):powers.append(multiply(M,R,*powers[-1]))
 diag=np.r_[16+number,16+number];h=diag*v;hr=abs(diag)*r+DOT_EPS*abs(h)+TINY
 for k in range(1,7):
  exact=F(s.binomial(s.Rational(1,2),k))/16**(2*k-1);coefficient=float(exact);coefficient_error=float(abs(exact-F.from_float(coefficient)))
  w,z=powers[2*k];hr+=coefficient_error*(abs(w)+z)
  h,hr=plus(h,hr,w,z,coefficient)
 rayleigh,rayleigh_error=dot_enclosure(v,r,h,hr)
 residual,residual_radius=plus(h,hr,v,r,-rayleigh);residual_radius+=rayleigh_error*(abs(v)+r)
 residual_bound=norm_upper(residual,residual_radius)
 moment14=norm_upper(*powers[14]);coefficient=F(33,2048*16**13);delta=upper(float(coefficient)*moment14*(1+DOT_EPS))
 lowR=lower(rayleigh-rayleigh_error);highR=upper(rayleigh+rayleigh_error);separation=lower(17-highR-delta);assert separation>0
 error=upper(residual_bound+delta);lo=lower(lowR-error*error/separation);hi=upper(highR+delta)
 band_distance=upper(2*error/separation);tail=upper(2*(error/separation)**2)
 old=json.loads((Path(__file__).resolve().parents[1]/'openwave_quantum_rotational/results.json').read_text());anchor=F('0.00051099895069');oldcap=next(x for x in old['interval_certificates']if x['cap']==2)
 elo,ehi=map(F.from_float,oldcap['actual_ground_energy_interval']);unitlo=anchor/ehi;unithi=anchor/elo
 gevlo=lower(float(F.from_float(lo)*unitlo));gevhi=upper(float(F.from_float(hi)*unithi))
 return {'trial_photon_cap':cap,'padded_exact_moment_cap':padded,'real_matrix_rows':M.shape[0],'max_row_nonzeros':maxrow,
  'trial_integer_coefficients':integers,'trial_squared_integer_norm':str(D),
  'lower_polynomial_Rayleigh_interval':[lowR,highR],'lower_polynomial_residual_norm_upper':residual_bound,
  'infinite_B14_trial_norm_upper':moment14,'square_root_vector_remainder_upper':delta,
  'infinite_photon_ground_energy_interval':[lo,hi],'ground_band_trace_distance_upper':band_distance,
  'infinite_ground_photon_tail_above_trial_cap_upper':min(1.,tail),
  'GeV_interval_using_identical_cap_two_electron_unit':[gevlo,gevhi],
  'arithmetic_contract':'IEEE-754 binary64 round-to-nearest, no overflow; sparse rows below 256 terms; 2^-40 CSR guard, math.fsum and 2^-48 scalar guards, 1e-300 underflow allowance; exact algebraic entry intervals at 45 digits',
  'scope':'Full infinite-photon operator at fixed complete matter shells. Finite-cap eigensolvers propose exact normalized dyadic witnesses only; padded moments and a global scalar Taylor remainder certify the untruncated field.'}

def analytic_trial():
 m=matter_exact();B2=m['B']**2+s.Rational(1,64)*sum((x*x for x in m['derivative']),s.zeros(20));v0=s.simplify(B2[0,0]/32)
 assert B2[10,10]==B2[0,0] and B2[0,10]==0 and 0<v0<s.Rational(1,16)
 return {'vacuum_spin_trial_shift_exact':str(v0),'vacuum_spin_trial_shift_upper':'1/16','full_infinite_source_lower':'16*N_ferm+dGamma(N_matter)+N_photon','lowest_reduced_band_dimension':2,'next_reduced_energy_floor':'17','exact_dressed_Casimir':'3/4','full_native_Clifford_ground_multiplicity':16,'finite_spatial_regulator_retained':True}

def certificate():
 rows=[certify(c)for c in(1,2,3)];intersection=[max(r['infinite_photon_ground_energy_interval'][0]for r in rows),min(r['infinite_photon_ground_energy_interval'][1]for r in rows)];assert intersection[0]<=intersection[1]
 return {'schema':'oph.infinite_photon_quantum_particle.v1','analytic_infinite_field':analytic_trial(),'certified_trials':rows,
  'joint_infinite_ground_energy_interval':intersection,
  'strict_tail_improvement':rows[2]['infinite_ground_photon_tail_above_trial_cap_upper']<rows[0]['infinite_ground_photon_tail_above_trial_cap_upper'],
  'cap_compression_boundary':'The operator is sqrt(16^2+B_infinity^2)+N_matter+N_photon. Finite trial vectors from the truncated source are embedded and freshly evaluated; sqrt(PBP squared) is not identified with P sqrt(B squared) P.',
  'field_coupling_mass_trap_and_electron_unit_source_selected':False,
  'scope':'An infinite-photon theorem and finite energy/tail certificates on the fixed spatial oscillator regulator. No spatial continuum, trap removal, physical electron identity or infinite-field magnetic coefficient is asserted.',
  'verdict':'PASS_INFINITE_PHOTON_SPIN_ENERGY_AND_TAIL_CERTIFICATE'}
