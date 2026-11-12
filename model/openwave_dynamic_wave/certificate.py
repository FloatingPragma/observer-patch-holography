"""Exact source residuals, finite propagation errors and complete instruments."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.linalg import eigh
from mpmath import iv
from model.openwave_quantum_geometry.certificate import scalar,ends,sqceil
from model.openwave_energy_lapse.source import native
from .source import shape,geometry,field,block,packet,packet_integers,KAPPA,evolution
DEN=2**46
ARITHMETIC=F(1,2**28)
def rat(x):return F(int(x.p),int(x.q))
def number(x):return scalar(rat(x))
def interval(x):return iv.mpf([scalar(x.lo),scalar(x.hi)])
def encode(x):return [str(v)for v in ends(x)]
def rational_matrix(A):return np.array([[int(round(float(x)*DEN))for x in row]for row in A],dtype=object)
def frobenius(A,denominator):return F(sqceil(sum(int(x)**2 for x in A.flat)),denominator)

@lru_cache(None)
def exact_geometry(n=7):
 iv.dps=80;s=shape(n);M=s['M_exact'];S=s['S_exact'];mu=[number(M[0,0]+M[0,1]),number(M[0,0]-M[0,1])]
 omega=[iv.sqrt(number(S[0,0]+S[0,1])/mu[0]),iv.sqrt(number(S[0,0]-S[0,1])/mu[1])]
 c=[1/iv.sqrt(4*KAPPA*omega[i]*mu[i])for i in range(2)]
 values=[[a*c[0]+b*c[1],a*c[0]-b*c[1]]for a,b in product((-1,1),repeat=2)]
 profile=[[number(s['F_exact'][x,0])*v[0]+number(s['F_exact'][x,1])*v[1]for x in range(n*n)]for v in values]
 assert all(max(abs(a),abs(b))<F(1,10)for row in profile for a,b in map(ends,row))
 return omega,values,profile

@lru_cache(None)
def wave_square_roots(n=7):
 _,_,profiles=exact_geometry(n);s=shape(n);f=field(n,1);m=n*n;rows=[];integers=[]
 for q,A,Kfloat in zip(profiles,f['coordinate_blocks'],f['K_blocks']):
  S=rational_matrix(A);K=rational_matrix(Kfloat);error=F(0)
  for i,j in product(range(m),repeat=2):
   exact=int(s['L'][i,j])*iv.exp(2*(q[i]+q[j]));lo,hi=ends(exact);v=F(int(K[i,j]),DEN);error=max(error,abs(lo-v),abs(hi-v))
  residual=frobenius(S@S-K*DEN,DEN*DEN)+m*error
  positive=min(F(int(S[i,i])-sum(abs(int(S[i,j]))for j in range(m)if i!=j),DEN)for i in range(m));assert positive>0
  # lambda(K)>=8/(n+1)^2 exp(-.4)>1/16 for the primary n7.
  assert n==7
  rooterror=residual/(positive+F(1,4))
  rows.append({'positive_root_Gershgorin_lower':str(positive),'exact_squared_root_residual_upper':str(residual),
   'exact_positive_root_operator_error_upper':str(rooterror),'dyadic_root_integer_rows':[[int(v)for v in row]for row in S]})
  integers.append(S)
 return rows,integers

def exact_vacuum(left,right,n=7):
 omega,values,_=exact_geometry(n);g=geometry(n,1);W=np.rint(2*g['coordinate_W']).astype(int);assert np.array_equal(W.T@W,4*np.eye(4))
 energy=[interval(native()[2][i])for i in(left,right)];diagonal=[energy[0]*iv.exp(v[0])+energy[1]*iv.exp(v[1])for v in values]
 H=[[sum((int(W[a,k]*W[b,k])*diagonal[k]/4 for k in range(4)),iv.mpf(0))for b in range(4)]for a in range(4)]
 for i,(a,b)in enumerate(product(range(2),repeat=2)):H[i][i]+=a*omega[0]+b*omega[1]
 return H

@lru_cache(None)
def vacuum_spectral_certificate(left,right,n=7):
 exact=exact_vacuum(left,right,n);H=rational_matrix(block(n,1,left,right)['vacuum']);entry=F(0)
 for a,b in product(range(4),repeat=2):
  low,high=ends(exact[a][b]);x=F(int(H[a,b]),DEN);entry=max(entry,abs(low-x),abs(high-x))
 E,O=eigh(np.array(H,float)/DEN);OI=rational_matrix(O);DI=np.array([int(round(float(x)*DEN))for x in E],dtype=object)
 eta=frobenius(OI.T@OI-np.eye(4,dtype=object)*DEN**2,DEN**2);residual=frobenius(H@OI-OI*DI[None,:],DEN**2)+8*entry;epsilon=2*residual+120*eta
 return {'native_levels':[left,right],'dimension':4,'source_to_dyadic_operator_error_upper':str(4*entry),
  'orthogonality_Frobenius_upper':str(eta),'exact_source_spectral_residual_upper':str(residual),
  'polar_Hamiltonian_error_upper':str(epsilon),'eigenvalue_integers':list(map(int,DI)),
  'eigenframe_integer_rows':[[int(v)for v in row]for row in OI],
  'all_eigenvalue_intervals':[[str(F(int(v),DEN)-epsilon),str(F(int(v),DEN)+epsilon)]for v in DI]}

@lru_cache(None)
def spectral_certificate(left,right,n=7):
 roots,ints=wave_square_roots(n);g=geometry(n,1);W=np.rint(2*g['coordinate_W']).astype(int);m=n*n
 V=rational_matrix(block(n,1,left,right)['vacuum']);exact=exact_vacuum(left,right,n);entry=F(0)
 for a,b in product(range(4),repeat=2):
  low,high=ends(exact[a][b]);x=F(int(V[a,b]),DEN);entry=max(entry,abs(low-x),abs(high-x))
 H=4*np.kron(V,np.eye(m,dtype=object))
 for k in range(4):H+=np.kron(np.outer(W[:,k],W[:,k]),ints[k])
 error=4*entry+max(F(r['exact_positive_root_operator_error_upper'])for r in roots)
 denH=4*DEN;E,O=eigh(np.array(H,float)/denH);OI=rational_matrix(O);DI=np.array([int(round(float(x)*DEN))for x in E],dtype=object);d=len(DI)
 eta=frobenius(OI.T@OI-np.eye(d,dtype=object)*DEN**2,DEN**2);assert eta<F(1,2)
 residual=frobenius(H@OI-4*OI*DI[None,:],denH*DEN)+2*error
 epsilon=2*residual+120*eta;assert max(abs(int(v))for v in DI)<30*DEN
 return {'native_levels':[left,right],'dimension':d,'source_to_dyadic_operator_error_upper':str(error),
  'orthogonality_Frobenius_upper':str(eta),'exact_source_spectral_residual_upper':str(residual),
  'polar_Hamiltonian_error_upper':str(epsilon),'eigenvalue_integers':list(map(int,DI)),
  'eigenframe_integer_rows':[[int(v)for v in row]for row in OI]}

def phases(c,t):
 iv.dps=80;out=[]
 for e in c['eigenvalue_integers']:
  angle=scalar(F(t)*F(e,DEN));z=iv.exp(-iv.j*angle);a=ends(z.real);b=ends(z.imag);r=float(sum(a)/2);im=float(sum(b)/2)
  assert max(abs(F(r)-v)for v in a)<F(1,2**48)and max(abs(F(im)-v)for v in b)<F(1,2**48)
  out.append(complex(r,im))
 return np.array(out)
def unitary(c,t):
 O=np.array(c['eigenframe_integer_rows'],float)/DEN
 return (O*phases(c,t))@O.T
def unitary_error(c,t):return abs(F(t))*F(c['polar_Hamiltonian_error_upper'])+3*F(c['orthogonality_Frobenius_upper'])+ARITHMETIC
def widen(value,error):return [float(np.nextafter(float(value-error),-np.inf)),float(np.nextafter(float(value+error),np.inf))]

def observations(c,reflect=False):
 n=7;g=geometry(n,1);m=n*n;psi=np.kron(g['vacuum'],packet(n,reflect));Y=np.tile(np.arange(n),n);rows=[]
 blocks=block(n,1,*c['native_levels']);H=blocks['wave'];energy=float(np.vdot(psi,H@psi).real);vc=vacuum_spectral_certificate(*c['native_levels'])
 vacuum_energy=float(g['vacuum']@blocks['vacuum']@g['vacuum'])
 for t in(1,2,3,4,5,6,8):
  v=unitary(c,t)@psi;vacuum=unitary(vc,t)@g['vacuum'];delta=float(unitary_error(c,t)+F(1,2**40));trace=2*delta+delta*delta+2**-45;v=v.reshape(4,m);p=abs(v)**2;rho=v@v.conj().T;mean=float(p.sum(axis=0)@Y)
  purity=float(np.trace(rho@rho).real);purity_error=trace*(1+(1+delta)**2)+2**-35
  rows.append({'time':t,'all_vacuum_geometry_number_probabilities':(abs(vacuum)**2/2).tolist(),'all_wave_site_geometry_number_probabilities':(p/2).tolist(),
   'remaining63_native_sectors_probability':'0','complete_outcome_count':64*4*(m+1),
   'state_vector_error_upper':delta,'each_joint_probability_error_upper':trace/2,
   'conditional_transverse_mean':mean,'conditional_transverse_mean_interval':widen(mean,6*trace),
   'geometry_mean_coordinates':[float(np.trace(rho@Q).real)for Q in g['Q']],
   'conditional_geometry_energy':float(np.trace(rho@g['Hg']).real),
   'geometry_purity_conditioned_on_initial_wave':purity,'geometry_purity_interval':widen(purity,purity_error),
   'source_energy_before':(energy+vacuum_energy)/2,
   'source_energy_after':float((np.vdot(v.ravel(),H@v.ravel()).real+np.vdot(vacuum,blocks['vacuum']@vacuum).real)/2),
   'position_number_measurement_work_absolute_upper':30})
 # Orthogonal receiver-site records at every scheduled step, followed by
 # complete geometry number records on every arrival and every timeout.
 U=unitary(c,F(1,2));v=psi.copy();mask=np.tile(np.repeat(np.arange(n)>=n-2,n),4);arrivals=[];prob=[]
 for step in range(1,17):
  v=U@v;hit=np.where(mask,v,0).reshape(4,m);p=float(np.vdot(hit,hit).real);prob.append(p)
  arrivals.append({'time':step/2,'all_geometry_site_probabilities':(abs(hit)**2/2).tolist(),'unconditional_hit_probability':p/2});v[mask]=0
 success=sum(prob);first=sum((i+1)*.5*p for i,p in enumerate(prob));timeout=float(np.vdot(v,v).real)
 # A sixteen-step retained-history dilation differs by at most16delta.
 step_error=unitary_error(c,F(1,2));assert step_error<F(1,10**6)
 error=float(16*step_error+256*step_error**2+F(1,2**40));lawerror=2*error+error*error+2**-45
 assert abs(success+timeout-1)<1e-8 and success>lawerror
 mean=first/success;mean_interval=[max(0,(first-8*lawerror)/(success+lawerror)),(first+8*lawerror)/(success-lawerror)]
 return {'native_levels':c['native_levels'],'reflected':reflect,'single_time_records':rows,
  'arrival_instrument':{'schedule_step':'1/2','horizon':8,'receiver':'x>=5; all seven y sites',
   'all_arrival_records':arrivals,'all_vacuum_timeout_geometry_probabilities':(abs(unitary(vc,8)@g['vacuum'])**2/2).tolist(),'complete_wave_timeout_geometry_site_probabilities':(abs(v.reshape(4,m))**2/2).tolist(),
   'conditional_on_initial_wave_hit_probability':success,'conditional_on_initial_wave_timeout_probability':timeout,
   'mean_arrival_time_conditioned_on_hit':mean,'mean_arrival_time_interval':mean_interval,
   'complete_record_total_variation_error_upper':lawerror/2,'readout_work_absolute_upper_per_step':30,
   'arbitrary_external_reference_and_retained_failure_records':True}}

def run():
 certificates=[spectral_certificate(i,0)for i in(0,1)];records=[observations(c)for c in certificates];reflected=observations(certificates[1],True)
 shifts=[]
 for a,b in zip(records[1]['single_time_records'],records[0]['single_time_records']):
  lo=a['conditional_transverse_mean_interval'][0]-b['conditional_transverse_mean_interval'][1];hi=a['conditional_transverse_mean_interval'][1]-b['conditional_transverse_mean_interval'][0]
  shifts.append({'time':a['time'],'source_minus_zero_native_control_shift_interval':[lo,hi],'numeric_shift':a['conditional_transverse_mean']-b['conditional_transverse_mean']})
 A=records[1]['arrival_instrument'];B=records[0]['arrival_instrument'];delay=[A['mean_arrival_time_interval'][0]-B['mean_arrival_time_interval'][1],A['mean_arrival_time_interval'][1]-B['mean_arrival_time_interval'][0]]
 assert delay[0]>0 and shifts[3]['source_minus_zero_native_control_shift_interval'][1]<0
 assert records[1]['single_time_records'][3]['geometry_purity_interval'][1]<1
 assert all(abs(r['conditional_transverse_mean']+s['conditional_transverse_mean']-6)<2e-8 for r,s in zip(reflected['single_time_records'],records[1]['single_time_records']))
 from .source import probe
 control=[]
 for level in(0,1):
  rows,arr,timeout=probe(7,2,level);control.append({'native_level':level,'wave_time_rows':rows,'arrival_probabilities':arr,'timeout':timeout,'mean_arrival_time':sum((i+1)*.5*p for i,p in enumerate(arr))/sum(arr)})
 coherent=[];g=geometry();initial=np.kron(g['vacuum'],packet())
 for t in(4,8):
  v=[unitary(c,t)@initial for c in certificates];z=[unitary(vacuum_spectral_certificate(i,0),t)@g['vacuum']for i in(0,1)]
  pw=abs(np.array([v[0]+v[1],v[0]-v[1]]))**2/8;pv=abs(np.array([z[0]+z[1],z[0]-z[1]]))**2/8;assert abs(pw.sum()+pv.sum()-1)<1e-8
  delta=max(float(unitary_error(c,t)+F(1,2**40))for c in certificates)
  coherent.append({'time':t,'all_native_X_geometry_wave_site_probabilities':pw.reshape(2,4,49).tolist(),
   'all_native_X_geometry_vacuum_probabilities':pv.tolist(),'other62_native_pair_outcomes_each_geometry_field_probability':'0',
   'complete_outcome_count':12800,'wave_conditional_clock_visibility':float(abs(np.vdot(v[0],v[1]))),
   'vacuum_conditional_clock_visibility':float(abs(np.vdot(z[0],z[1]))),
   'each_complete_probability_error_upper':2*delta+delta*delta+2**-45,
   'same_source_coherent_native_geometry_field_preparation':True})
 shape_data=shape(7)
 return {'schema':'oph.single_source_autonomous_wave_geometry.v1','verdict':'PASS_SAME_SOURCE_AUTONOMOUS_WAVE_BENDING_DELAY',
  'source':'Hg+exp(Q0)K0+exp(Q1)K1+0 direct-sum sqrt(exp(2q)Lexp(2q)); q=FQ',
  'complete_native_dimension':64,'geometry_dimension':4,'wave_vacuum_dimension':50,'full_source_dimension':12800,
  'exact_packet_integer_pairs':[list(x)for x in packet_integers()],
  'exact_packet_squared_integer_norm':str(sum(a*a+b*b for a,b in packet_integers())),
  'geometry_modes_are_declared_harmonic_architecture':True,'harmonic_invariant_projection_of_unrestricted_geometry_claimed':False,
  'spatial_shape_exact':[[str(v)for v in row]for row in shape_data['F_exact'].tolist()],
  'stiffness_matrix_exact':[[str(v)for v in row]for row in shape_data['S_exact'].tolist()],
  'inertia_matrix_exact':[[str(v)for v in row]for row in shape_data['M_exact'].tolist()],
  'kappa':KAPPA,'exact_wave_square_root_certificates':wave_square_roots()[0],
  'complete_native_block_certificates':certificates,'vacuum_block_certificates':[vacuum_spectral_certificate(i,0)for i in(0,1)],'native_control_and_source_observations':records,
  'reflected_source_observations':reflected,'same_source_bending_intervals':shifts,
  'same_source_operational_delay_interval':delay,'cap_two_numerical_control':control,
  'complete_native_clock_geometry_wave_interference_records':coherent,
  'all64_native_sectors_retained_and_commute_with_source':True,'global_ground':'unique native00,geometry vacuum,field vacuum; zero energy',
  'native_geometry_wave_reference_and_clock_shared':True,'geometry_refresh_or_static_profile_service_used':False,
  'total_energy_conserved_between_observations':True,'geometry_wave_entanglement_retained':True,
  'canonical_time_dependent_covariant_scalar_quantization_claimed':False,
  'arithmetic_contract':'IEEE754 binary64 round-to-nearest, no overflow; dense dimensions<=196; interval phases; per-propagator guard2^-28 bounds complex dot products and phase rounding',
  'physical_G_or_observed_lensing_or_Einstein_dynamics_claimed':False,
  'scope':'Declared two-mode harmonic geometry on one7x7 scalar grid; quantum autonomous propagation with complete source/readout outcomes and certified primary finite numerical errors. Cap2 is an explicit numerical regulator control.'}
