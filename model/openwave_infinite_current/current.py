"""Global differentiated remainder and padded current/energy witnesses."""
from fractions import Fraction as F
from pathlib import Path
from math import prod,fsum,sqrt
import json
import numpy as np
import sympy as s
from scipy.sparse import csr_matrix,kron,eye,bmat
from mpmath import iv
from model.openwave_quantum_rotational.source import matter_exact,source,occupations,ladders
from model.openwave_quantum_rotational.enclosure import scalar,interval_entries,multiply as ivmultiply,norm_upper as ivnorm,number as ivnumber
from model.openwave_infinite_photon.field import padded_operator,multiply,plus,dot_enclosure,norm_upper,upper,lower,DOT_EPS,TINY

TRIAL_CAP=4

def norm_controls():
 m=matter_exact();E=m['external'][2];assert s.trace(E*E)==5
 records=[]
 for label,M in [('external',E),*((str(i),x)for i,x in enumerate(m['derivative']))]:
  A=s.Rational(16,25)*s.eye(20)-M*M;L,D=A.LDLdecomposition(hermitian=True);assert L*D*L.H==A and all(x>0 for x in D.diagonal());records.append({'operator':label,'strict_norm_upper':'4/5','positive_LDL_diagonal':[str(x)for x in D.diagonal()]})
 assert m['B']**4==s.Rational(5,2)*m['B']**2
 return {'B0_norm_upper':'8/5','all_probe_and_field_matrix_controls':records,'probe_Hilbert_Schmidt_squared_exact':'5','uniform_infinite_current_norm_upper':'sqrt(5)','uniform_kernel_second_derivative_norm_upper':'141/400'}


def projected_trial(cap,padded):
 iv.dps=55;scalar.cache_clear();v=source(cap)['psi'];ints=[[int(round(z.real*2**40)),int(round(z.imag*2**40))]for z in v];D=sum(a*a+b*b for a,b in ints);phi=[iv.mpc(a,b)/iv.sqrt(D)for a,b in ints]
 m=matter_exact();fb,a=ladders(cap,True);d=len(fb);jf=-s.I*(a[0].H*a[1]-a[1].H*a[0]);J=s.kronecker_product(m['J'][2],s.eye(d))+s.kronecker_product(s.eye(20),jf);entries=interval_entries(J.todok());limit=2*cap+5
 for twice in range(-limit,limit+1,2):
  if twice==1:continue
  jphi=ivmultiply(entries,phi);phi=[(x-iv.mpf(twice)*y/2)/((1-twice)/2)for x,y in zip(jphi,phi)]
 nrm=ivnorm(phi);assert float(nrm.a)>.99;phi=[z/nrm for z in phi]
 large=occupations(padded);where={x:i for i,x in enumerate(large)};N=20*len(large);v=np.zeros(2*N);r=np.zeros(2*N)
 for i,z in enumerate(phi):
  j=(i//d)*len(large)+where[fb[i%d]]
  for index,part in((j,z.real),(j+N,z.imag)):
   lo=lower(float(part.a));hi=upper(float(part.b));mid=(lo+hi)/2;v[index]=mid;r[index]=upper(max(mid-lo,hi-mid))
 return v,r,ints,D


def probe(padded):
 E=matter_exact()['external'][2];mid=np.zeros((20,20),complex);re=np.zeros((20,20));im=np.zeros((20,20))
 for (i,j),z in E.todok().items():
  values=[]
  for part in(s.re(z),s.im(z)):
   if part==0:values.append((0.,0.));continue
   exact=scalar(part);lo=lower(float(exact.a));hi=upper(float(exact.b));v=(lo+hi)/2;values.append((v,upper(max(v-lo,hi-v))))
  mid[i,j]=values[0][0]+1j*values[1][0];re[i,j]=values[0][1];im[i,j]=values[1][1]
 d=len(occupations(padded));I=eye(d,format='csr');M=kron(csr_matrix(mid),I,format='csr');Rr=kron(csr_matrix(re),I,format='csr');Ri=kron(csr_matrix(im),I,format='csr')
 return bmat([[M.real,-M.imag],[M.imag,M.real]],format='csr'),bmat([[Rr,Ri],[Ri,Rr]],format='csr')


def certificate():
 controls=norm_controls();C=TRIAL_CAP;padded=C+14;M,R,number,basis,maxrow=padded_operator(padded);v,r,ints,D=projected_trial(C,padded);powers=[(v,r)]
 for _ in range(14):powers.append(multiply(M,R,*powers[-1]))
 diag=np.r_[16+number,16+number];h=diag*v;hr=abs(diag)*r+DOT_EPS*abs(h)+TINY
 for k in range(1,7):h,hr=plus(h,hr,*powers[2*k],float(F(s.binomial(s.Rational(1,2),k))/16**(2*k-1)))
 ray,rayerr=dot_enclosure(v,r,h,hr);res,resr=plus(h,hr,v,r,-ray);resr+=rayerr*(abs(v)+r);resnorm=norm_upper(res,resr);delta=upper(float(F(33,2048*16**13))*norm_upper(*powers[14])*(1+DOT_EPS));sep=lower(17-upper(ray+rayerr)-delta);error=upper(resnorm+delta);elo=lower(ray-rayerr-error*error/sep);ehi=upper(ray+rayerr+delta);trace=upper(2*error/sep)
 E,ER=probe(padded);ep=[multiply(E,ER,*p)for p in powers[:12]];values=[];errors=[]
 for k in range(1,7):
  coefficient=float(F(s.binomial(s.Rational(1,2),k))/16**(2*k-1))
  for j in range(2*k):
   value,err=dot_enclosure(*powers[j],*ep[2*k-1-j]);values.append(coefficient*value);errors.append(abs(coefficient)*err)
 current=fsum(values);arithmetic=upper(fsum(errors)+DOT_EPS*fsum(abs(x)for x in values))
 # Exact field-coordinate radial moment: <r^26> <= product_j(2(C+j)+3).
 moment=prod(2*(C+j)+3 for j in range(13));boundB=s.Rational(8,5)+s.sqrt(6)/10*s.Integer(moment)**s.Rational(1,26)
 derivative_remainder=s.Rational(14*33,2048*16**13)*s.sqrt(5)*boundB**13
 remainder=upper(float(scalar(derivative_remainder).b));state_error=upper(sqrt(5)*trace);muerror=upper(arithmetic+remainder+state_error);mulo=lower(current-muerror);muhi=upper(current+muerror);assert mulo>0
 glow=lower(4*elo*mulo);ghigh=upper(4*ehi*muhi)
 times=[]
 for b in(F(1,10**6),F(1,10**8)):
  change=iv.sqrt(5)*ivnumber(b);gap=17-iv.mpf(ehi)-2*change;Denergy=ivnumber(F(141,800)*b*b)+change*change/gap
  for tau in(1,10,100):
   time=tau/float(b);eta=2*iv.sqrt(2)*change/gap+ivnumber(F(tau)/b)*Denergy+tau*iv.mpf(muerror)
   times.append({'field_exact':str(b),'scaled_time_b_t':tau,'source_time':time,'complete_source_arbitrary_reference_trace_distance_upper':upper(min(2.,float((2*eta).b)))})
 return {'schema':'oph.infinite_photon_native_current_larmor.v1','exact_matrix_norm_controls':controls,
  'trial_photon_cap':C,'exact_moment_padding':padded,'exact_Jz_target':'1/2','trial_integer_coefficients':ints,'trial_squared_integer_norm':str(D),
  'infinite_ground_energy_interval':[elo,ehi],'infinite_ground_line_trace_distance_upper':trace,
  'lower_polynomial_current_expectation':current,'polynomial_current_arithmetic_error_upper':arithmetic,
  'field_radial_moment26_integer_upper':str(moment),'global_differentiated_remainder_upper':remainder,'ground_state_current_transfer_error_upper':state_error,
  'actual_infinite_current_derivative_interval':[mulo,muhi],'rest_energy_normalized_infinite_g_interval':[glow,ghigh],
  'complete_finite_time_Larmor_bounds':times,
  'all_photon_field_reference_and_failure_outputs_retained':True,
  'scope':'Full photon Fock source with fixed complete spatial shells, same mass/profile/coupling and source clock. B is unbounded; finite coordinate-matrix Hilbert-Schmidt and global remainder estimates replace a false global B-norm assumption. No physical electron anomaly or field calibration is inferred.',
  'verdict':'PASS_INFINITE_PHOTON_CURRENT_AND_FULL_REFERENCE_LARMOR'}
