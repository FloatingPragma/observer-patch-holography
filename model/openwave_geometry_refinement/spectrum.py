"""Finite spectral and Schur certificates for the fixed geometry source."""
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
from scipy.linalg import eigh
from .source import DEN,UNIT,gamma,upper,frob_integer,vacuum,block,compression,indices,exact_geometry,ends,native,iv,scalar

def sqrt_upper(x):
 x=F(x);return F(math.isqrt(x.numerator*x.denominator)+1,x.denominator)
def matrix_frobenius(A):return sqrt_upper(F(sum(int(x)**2 for x in A.flat),DEN**2))
def numeric_norm_upper(A):
 # A is a stored binary64 array. Squaring and summing at most1225^2
 # components is enclosed by gamma_(2size), with a square-root guard.
 raw=float(np.linalg.norm(A));assert np.isfinite(raw)
 return F(raw)*(1+gamma(2*A.size+2))+F(1,2**90)

def native_certificate(left,right):
 a=vacuum(left,right);HI=a['integer_H'].astype(object);e,V=eigh(a['H']);OI=np.rint(V*DEN).astype(np.int64).astype(object);DI=np.rint(e*DEN).astype(np.int64).astype(object);d=len(e)
 eta=frob_integer(OI.T@OI-np.eye(d,dtype=object)*DEN**2,DEN**2);R=frob_integer(HI@OI-OI*DI[None,:],DEN**2)+2*a['error'];epsilon=upper(2*R+128*eta)
 assert eta<F(1,2) and max(abs(int(x))for x in DI)<32*DEN
 return {'native_levels':[left,right],'dimension':d,'source_H_error_upper':str(upper(a['error'])),'exact_Gram_error_upper':str(eta),'exact_residual_upper':str(upper(R)),'spectral_error_upper':str(epsilon),
  'dyadic_denominator':DEN,'eigenframe_integer_rows':[[int(x)for x in row]for row in OI],'eigenvalue_integers':list(map(int,DI)),
  'all_eigenvalue_intervals':[[str(F(int(x),DEN)-epsilon),str(F(int(x),DEN)+epsilon)]for x in DI]}

def compression_spectrum(a):
 H=a['H'];e,O=eigh(H);d=len(e);OI=np.rint(O*DEN).astype(np.int64);DI=np.rint(e*DEN).astype(np.int64);V=np.array(OI,float)/DEN;D=np.array(DI,float)/DEN
 fV=matrix_frobenius(OI);fH=matrix_frobenius(a['integer_H']);maxD=max(abs(F(int(x),DEN))for x in DI);assert maxD<32
 gram=V.T@V-np.eye(d);eta=numeric_norm_upper(gram)+gamma(2*d)*fV*fV+gamma(2)*(fV*fV+d)
 residual=H@V-V*D;rounding=gamma(2*d)*fH*fV+gamma(3)*(fH*fV+fV*maxD);R=numeric_norm_upper(residual)+rounding+2*a['error'];eps=upper(2*R+128*eta)
 assert eta<F(1,2)
 return {'eigenvalues':DI,'eigenframe':OI,'error':eps,'orthogonality_upper':upper(eta),'residual_upper':upper(R),
  'proof_contract':'binary64 matrix products and Frobenius sums with gamma_n forward-error bounds; exact dyadic input matrices'}

def wave_schur(left=1,cap=3):
 full=block(left,0);small=compression(left,0,cap);spectral=compression_spectrum(small);ix=small['indices'];outside=np.setdiff1d(np.arange(1225),ix);B=full['integer_H'][np.ix_(outside,ix)]
 rows=max(sum(abs(int(v))for v in row)for row in B);cols=max(sum(abs(int(v))for v in B[:,j])for j in range(B.shape[1]));beta=upper(sqrt_upper(F(rows*cols,DEN**2))+full['error'])
 omega=ends(exact_geometry()['omega'][0])[0];minimum_native=ends(scalar(native()[2][left].lo)*iv.exp(-scalar(F(19,100))))[0];b=upper((cap+1)*omega+minimum_native+F(1,5),80)-F(1,2**80)
 a2=F(int(spectral['eigenvalues'][1]),DEN)-spectral['error'];second_lower=(a2+b-sqrt_upper((b-a2)**2+4*beta*beta))/2
 # A dyadic compressed ground vector is embedded before any observable
 # comparison. Its Rayleigh and squared residual are exact rationals.
 p=np.zeros(1225,dtype=object);p[ix]=spectral['eigenframe'][:,0].astype(object);norm=sum(int(v)**2 for v in p);w=full['integer_H'].astype(object)@p;A=sum(int(x)*int(y)for x,y in zip(p,w));rayleigh=F(A,DEN*norm);resnumer=sum((norm*int(x)-A*int(y))**2 for x,y in zip(w,p));res2=F(resnumer,DEN**2*norm**3);res=sqrt_upper(res2)+full['error'];low_ray=rayleigh-full['error'];high_ray=rayleigh+full['error'];assert second_lower>high_ray
 lower=low_ray-res*res/(second_lower-high_ray);upper_energy=high_ray
 return {'native_levels':[left,0],'nested_trial_cap':cap,'compression_dimension':len(ix),'complete_wave_dimension':1225,
  'complement_energy_lower':str(b),'offdiagonal_norm_upper':str(beta),'compressed_second_eigenvalue_lower':str(a2),'master_second_eigenvalue_lower':str(second_lower),
  'exact_trial_integer_vector':list(map(int,p)),'exact_trial_squared_norm':str(norm),'exact_dyadic_Rayleigh':str(rayleigh),'exact_dyadic_residual_squared':str(res2),
  'physical_source_residual_upper':str(upper(res)),'master_ground_interval_rational':[str(lower),str(upper_energy)],
  'master_ground_interval':[float(np.nextafter(float(lower),-np.inf)),float(np.nextafter(float(upper_energy),np.inf))],
  'compression_spectral_error_upper':str(spectral['error']),'compression_polar_contract':spectral['proof_contract'],
  'method':'finite Schur lower bound on the second master level, followed by a ground-state Temple enclosure using the exact embedded trial'}

def spectra():
 clocks=[native_certificate(i,j)for i,j in product(range(8),repeat=2)];c=exact_geometry()['c'];bound=native()[2][1].lo*sum(ends(x)[0]**4 for x in c)/12
 return {'all64_native_vacuum_spectra':clocks,'wave_Schur_ground_certificates':[wave_schur(i,c)for i in(0,1)for c in(1,2,3)],
  'strict_reexponentiation_counterexample':{'matrix_element':'native10 geometry00 diagonal; master cap4 compression minus cap1 re-exponentiation','difference_lower':str(bound),
   'proof':'positive closed four-step oscillator walks through occupation2 add(c_plus^4+c_minus^4)/12; all higher walk coefficients are nonnegative'},
  'infinite_geometry_or_spatial_continuum_limit_claimed':False}
