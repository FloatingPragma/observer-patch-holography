"""Projected scalar stress, finite Ehrenfest corrections and dressed clocks."""
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.linalg import eigh,expm,solve_sylvester
from mpmath import iv
from model.openwave_energy_lapse.source import native
from model.openwave_dynamic_wave.source import shape,geometry,field,block,packet,energies,KAPPA
from model.openwave_dynamic_wave.certificate import spectral_certificate,vacuum_spectral_certificate,unitary,unitary_error,exact_geometry,scalar,ends,DEN

def operators(n=7,cap=1,left=1,right=0):
 sh=shape(n);g=geometry(n,cap);f=field(n,cap);d=g['dim'];m=n*n;I=np.eye(m)
 Q=[np.kron(q,I)for q in g['Q']];Hg=np.kron(g['Hg'],I);H=block(n,cap,left,right)['wave'];Hf=f['H'];Vm=H-Hg-Hf
 velocity=[1j*(g['Hg']@q-q@g['Hg'])for q in g['Q']]
 Psmall=[KAPPA*sum(sh['M'][a,b]*velocity[b]for b in range(2))for a in range(2)];P=[np.kron(p,I)for p in Psmall]
 W=np.kron(g['coordinate_W'],I);forces=[];sylvester=[]
 for a in range(2):
  D=np.diag(sh['F'][:,a]);out=np.zeros((d*m,d*m));res=0.
  for k,(root,K)in enumerate(zip(f['coordinate_blocks'],f['K_blocks'])):
   rhs=2*(D@K+K@D);G=solve_sylvester(root,root,rhs);res=max(res,float(np.linalg.norm(root@G+G@root-rhs)))
   out[k*m:(k+1)*m,k*m:(k+1)*m]=G
  forces.append(W@out@W.T);sylvester.append(res)
 native_force=[energies()[level]*np.kron(expm(g['Q'][a]),I)for a,level in enumerate((left,right))]
 actual=[1j*(H@p-p@H)for p in P];boundary=[1j*((Vm+Hf)@p-p@(Vm+Hf))+native_force[a]+forces[a]for a,p in enumerate(P)]
 spring=[-KAPPA*sum(sh['S'][a,b]*Q[b]for b in range(2))for a in range(2)]
 assert max(np.linalg.norm(actual[a]-spring[a]+native_force[a]+forces[a]-boundary[a])for a in range(2))<1e-10
 free_velocity_error=max(np.linalg.norm(1j*(H@Q[a]-Q[a]@H)-sum(np.linalg.inv(KAPPA*sh['M'])[a,b]*P[b]for b in range(2)))for a in range(2));assert free_velocity_error<1e-10
 return {'Q':Q,'P':P,'Psmall':Psmall,'H':H,'Hg':Hg,'Hf':Hf,'Vm':Vm,'field_force':forces,'native_force':native_force,'actual_momentum_rate':actual,'finite_boundary_correction':boundary,'spring':spring,'sylvester_residuals':sylvester,'velocity_residual':free_velocity_error}

def expval(v,A):return float(np.vdot(v,A@v).real)
def dynamics():
 c=spectral_certificate(1,0);g=geometry();o=operators();psi=np.kron(g['vacuum'],packet());vc=vacuum_spectral_certificate(1,0);vacH=block()['vacuum'];rows=[];instrument=[]
 for t in(1,2,3,4,5,6,8):
  v=unitary(c,t)@psi;vac=unitary(vc,t)@g['vacuum'];waveg=v.reshape(4,49);rho=waveg@waveg.conj().T
  rates=[expval(v,1j*(o['H']@o[k]-o[k]@o['H']))for k in('Hg','Vm','Hf')];assert abs(sum(rates))<1e-10
  rows.append({'time':t,'wave_induced_geometry_coordinate_difference':[expval(v,Q)-expval(vac,q)for Q,q in zip(o['Q'],g['Q'])],
   'wave_conditional_geometry_energy':expval(v,o['Hg']),'vacuum_conditional_geometry_energy':expval(vac,g['Hg']),
   'source_energy':expval(v,o['H']),'geometry_native_field_energy_rates':rates,
   'geometry_purity':float(np.trace(rho@rho).real),'native_force':[expval(v,A)for A in o['native_force']],
   'field_force':[expval(v,A)for A in o['field_force']],
   'spring_force':[expval(v,A)for A in o['spring']],
   'momentum_rate':[expval(v,A)for A in o['actual_momentum_rate']],
   'finite_boundary_correction':[expval(v,A)for A in o['finite_boundary_correction']],
   'state_vector_error_upper':float(unitary_error(c,t)+F(1,2**40)),
   'force_values_are_numerical_evaluations_not_interval_enclosures':True})
  if t in(4,8):
   for a in range(2):
    e,V=eigh(o['field_force'][a]);p=abs(V.conj().T@v)**2/2;pg,Z=eigh(o['Psmall'][a]);joint=abs(Z.conj().T@waveg)**2/2;pv=abs(Z.conj().T@vac)**2/2
    assert abs(p.sum()+.5-1)<1e-8 and abs(joint.sum()+pv.sum()-1)<1e-8
    instrument.append({'time':t,'source_coordinate':a,
     'complete_field_stress_eigenvalues':e.tolist(),'all196_field_stress_probabilities':p.tolist(),
     'all4_vacuum_geometry_stress_zero_probabilities':(abs(vac)**2/2).tolist(),
     'geometry_momentum_eigenvalues':pg.tolist(),'all_geometry_momentum_wave_site_probabilities':joint.tolist(),
     'all_geometry_momentum_vacuum_probabilities':pv.tolist(),'remaining63_native_sectors_probability':'0',
     'complete_outcome_count_per_instrument':12800,'source_work_absolute_upper_per_readout':30,
     'position_and_momentum_are_separate_complete_instruments':True})
 # The uniform ambient q variation gives dK/ds=4K and hence dHfield/ds=2Hfield.
 uniform=[]
 for H,K in zip(field()['coordinate_blocks'],field()['K_blocks']):
  Fwave=solve_sylvester(H,H,4*K);uniform.append(float(np.linalg.norm(Fwave-2*H)))
 assert max(uniform)<1e-10
 refined=operators(cap=2);gg=geometry(cap=2);ee,VV=eigh(refined['H']);initial=np.kron(gg['vacuum'],packet());vv=VV@(np.exp(-4j*ee)*(VV.T@initial))
 refinement={'cap':2,'time':4,'complete_wave_geometry_dimension':441,
  'field_force':[expval(vv,A)for A in refined['field_force']],
  'finite_boundary_correction':[expval(vv,A)for A in refined['finite_boundary_correction']],
  'cap_independent_limit_claimed':False}
 return {'response_rows':rows,'complete_stress_and_momentum_instruments':instrument,
  'field_derivative_Sylvester_residuals':o['sylvester_residuals'],'exact_free_geometry_velocity_residual':o['velocity_residual'],
  'uniform_field_scale_Ward_identity':'sum_x dHfield/dq_x = 2Hfield',
  'uniform_native_scale_Ward_identity':'sum_a dHnative/dQ_a = Hnative',
  'uniform_field_identity_numeric_residuals':uniform,
  'field_only_and_native_energy_rates_cancel_geometry_rate':True,
  'finite_CCR_boundary_correction_discarded':False,
  'cap_two_numerical_boundary_control':refinement,
  'coupling_off_control':'H=Hg+constant native energies+flat Hfield; product geometry vacuum is stationary and no geometry-field entanglement is generated'}

def clocks():
 certificates=[vacuum_spectral_certificate(i,j)for i,j in product(range(8),repeat=2)];rows=[];records=[]
 ground=lambda i,j:list(map(F,certificates[8*i+j]['all_eigenvalue_intervals'][0]))
 for j in range(1,8):
  gaps=[];vectors=[]
  for i in(0,1):
   Ej,E0=ground(i,j),ground(i,0);gap=[Ej[0]-E0[1],Ej[1]-E0[0]];assert gap[0]>0;gaps.append(gap)
   e,A=eigh(block(7,1,i,0)['vacuum']);f,B=eigh(block(7,1,i,j)['vacuum']);v=A[:,0];w=B[:,0]
   if v[0]<0:v=-v
   if w[0]<0:w=-w
   overlap=float(v@w);vectors.append((e[0],f[0],v,w,overlap))
   for t in(0,1,4,8):
    first=np.exp(-1j*t*e[0])*v;second=np.exp(-1j*t*f[0])*w;X=np.array([first+second,first-second])/2;p=abs(X)**2;assert abs(p.sum()-1)<1e-12
    records.append({'background_native_level':i,'clock_native_level':j,'time':t,
     'all_X_clock_geometry_number_probabilities':p.tolist(),'clock_visibility':overlap,
     'other_native_or_wave_outcomes_probability':'0','full_source_outcome_count':12800,
     'source_energy':float((e[0]+f[0])/2),'stationary_dressed_branch_preparation_is_declared':True})
  ratio=[gaps[1][0]/gaps[0][1],gaps[1][1]/gaps[0][0]]
  rows.append({'native_probe_level':j,'bare_native_energy':float(energies()[j]),
   'zero_background_dressed_gap_interval':list(map(str,gaps[0])),
   'native_background_dressed_gap_interval':list(map(str,gaps[1])),
   'same_source_clock_frequency_ratio_interval':list(map(str,ratio)),
   'same_source_clock_shift_ppm_interval':[float((x-1)*10**6)for x in ratio]})
 lo=max(F(r['same_source_clock_frequency_ratio_interval'][0])for r in rows)-min(F(r['same_source_clock_frequency_ratio_interval'][1])for r in rows)
 hi=max(F(r['same_source_clock_frequency_ratio_interval'][1])for r in rows)-min(F(r['same_source_clock_frequency_ratio_interval'][0])for r in rows);assert lo>0
 # Exact finite weak-clock domain: f(s)=ground(Hg+Esource expQ0+s expQ1).
 # qmax<.067, omega_min>.698 and total source coefficient<3.85 give
 # a uniform gap >.18. Radius(expQ1)<.068 gives the stated rational bound.
 omega,values,_=exact_geometry();assert ends(omega[0])[0]>F(349,500)
 qbound=max(max(abs(a),abs(b))for v in values for a,b in map(ends,v));assert qbound<F(67,1000)
 assert native()[2][1].hi+native()[2][7].hi<F(77,20)
 width=iv.exp(scalar(F(67,1000)))-iv.exp(-scalar(F(67,1000)));assert ends(width)[1]<F(671,5000)
 gap=F(349,500)-F(77,20)*F(671,5000);assert gap>F(9,50)
 radius=F(17,250);coefficient=radius*radius/F(9,50)
 return {'all64_native_vacuum_block_certificates':certificates,'all7_same_source_clock_ratios':rows,
  'strict_clock_ratio_spread_interval':list(map(str,(lo,hi))),'strict_clock_ratio_spread_ppm_interval':[float(lo*10**6),float(hi*10**6)],
  'complete_dressed_clock_records':records,
  'exact_test_clock_theorem':{'source_probe_coefficient_domain':'0 <= s <= E7; background E0 or E1',
   'uniform_geometry_gap_lower':'9/50','radius_of_exp_Q_about_scalar_upper':'17/250',
   'effective_lapse_error_upper_rule':'abs((f(s)-f(0))/s - <ground(0)|exp(Qprobe)|ground(0)>) <= coefficient*s',
   'coefficient_exact':str(coefficient),'at_probe_coefficient_one_thousandth_error_upper':str(coefficient/F(1000)),
   'proof':'finite variational concavity, Hellmann-Feynman and reduced-resolvent gap bound; no canonical-commutator or Einstein premise'},
  'exact_universal_finite_strength_clock_ratio':False,
  'physical_equivalence_principle_or_clock_calibration_claimed':False}

def run():
 return {'schema':'oph.autonomous_source_stress_clock_response.v1',
  'verdict':'PASS_PROJECTED_DYNAMIC_STRESS_AND_DRESSED_CLOCK_DOMAIN',
  'same_source_as':'model/openwave_dynamic_wave/source.py',
  'complete_native_geometry_scalar_vacuum_and_reference_domain':True,
  'stress_backreaction':dynamics(),'dressed_clocks':clocks(),
  'geometry_architecture':'two harmonic coordinates with exact induced stiffness F^T L F and inertia F^T F; product oscillator cap1',
  'physical_Einstein_stress_tensor_G_or_universal_observed_redshift_claimed':False,
  'scope':'Exact finite action derivatives and boundary-corrected Ehrenfest identity, numerical complete stress records, interval-certified clock gaps and quantitative weak-probe universality. No continuum or observed gravity identification.'}
