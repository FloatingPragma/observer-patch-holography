"""Native288-CAR realization of the calibrated neutral tensor."""
from itertools import product
from fractions import Fraction as F
import numpy as np
from model.openwave_dirac_neutral_current.neutral import ix,bog,complete_vertex,free_source_check,exact_spin_contraction,flight,ward
from model.openwave_weak_decay.decay import interaction
from .calibration import calibration,catalogue,rigid_obstruction
from .antiparticle import complete_antiparticle_vertex,exact_antiparticle_spin_contraction

def construction():
 c=calibration();charged=c['charged'];mu=c['neutral'];R=c['R'];delta=1-charged[0]
 q=np.sqrt((delta**2-(mu[0]+mu[1])**2)*(delta**2-(mu[0]-mu[1])**2))/(2*delta)
 assert abs(q*q-float(c['source_recoil_q_squared_exact']))<1e-15
 mapping={};energy=np.zeros(288)
 for p,s,k,f in product((-1,0,1),range(2),range(16),range(3)):
  mass=charged[f]if k in(9,11)else mu[f]if k in(10,15)else 0.
  energy[ix(p,s,k,f)]=np.hypot(p*q,mass)
 for p,s in product((-1,0,1),range(2)):
  sign=1-2*s
  for f in range(3):
   for component,partner,mass in((11,9,charged[f]),(10,15,mu[f])):
    E,u,v=bog(p*q,sign,mass);a=ix(p,s,component,f);b=ix(-p,1-s,partner,f);energy[a]=energy[b]=E
    if component==11:mapping[a]=[(a,False,u),(b,True,-v)]
    mapping[b]=[(b,False,u),(a,True,v)]
  for alpha in range(3):
   row=[]
   for f,mass in enumerate(mu):
    E,u,v=bog(p*q,sign,mass);a=ix(p,s,10,f);b=ix(-p,1-s,15,f)
    row.extend(((a,False,R[alpha,f]*u),(b,True,-R[alpha,f]*v)))
   mapping[ix(p,s,10,alpha)]=row
  for component in range(16):
   if component in(9,10,11,15):continue
   for f in range(3):a=ix(p,s,component,f);mapping[a]=[(a,sign*p*q>0,1.)]
 assert len(mapping)==288 and q>.49
 return {'charged':charged,'neutral':mu,'q':q,'R':R,'mapping':mapping,'energies':energy}

def vacuum_kernel(d,phase,antineutrino=False):
 # This comparison kernel is derived from the positive source dispersion.
 # Phase is Delta_m31^2 t/(2q), using the same specified source parameter.
 mu=d['neutral'];ratio=(mu*mu-mu[0]**2)/(mu[-1]**2-mu[0]**2)
 U=d['R'].conj()if antineutrino else d['R'];A=(U*np.exp(-1j*phase*ratio))@U.conj().T
 return abs(A)**2

def complete_flight(d,labels,W,source_flavours=(1,0)):
 A=np.zeros((3,3),complex);failure=0.
 for k,(o,i,j,se,sn,sbar)in enumerate(labels):
  if(o,se,sn,sbar)==(1,0,1,1):A[i,j]=W[k,1]
  else:failure+=float(abs(W[k,1])**2)
 mu=d['neutral'];E=np.hypot(d['q'],mu);gap=(mu*mu-mu[0]**2)/(E+E[0]);delta=mu[-1]**2-mu[0]**2
 nu=[d['R'][f:f+1,:]*np.array([bog(d['q'],-1,m)[1]for m in mu])for f in range(3)]+[np.diag([abs(bog(d['q'],-1,m)[2])for m in mu])]
 antinu=[d['R'][f:f+1,:].conj()*np.array([abs(bog(d['q'],1,m)[2])for m in mu])for f in range(3)]+[np.diag([bog(d['q'],1,m)[1]for m in mu])]
 for Ks in(nu,antinu):assert np.linalg.norm(sum(K.conj().T@K for K in Ks)-np.eye(3))<1e-14
 rows=[]
 for phase in (0.,np.pi/2,np.pi,3*np.pi/2,2*np.pi,64*np.pi):
  t=phase*2*d['q']/delta;v=np.exp(-1j*t*gap);state=v[:,None]*A*v[None,:]
  prob=np.array([[np.linalg.norm(K@state@J.T)**2 for J in antinu]for K in nu])
  ideal=np.zeros((4,4));ideal[:3,:3]=vacuum_kernel(d,phase)[:,source_flavours[0],None]*vacuum_kernel(d,phase,True)[None,:,source_flavours[1]]
  before=float(np.sum(abs(state)**2*(d['charged'][0]+E[:,None]+E[None,:])));after=0.
  for a,K in enumerate(nu):
   for b,J in enumerate(antinu):
    output=K@state@J.T;left=E if a==3 else np.zeros(1);right=E if b==3 else np.zeros(1)
    after+=float(np.sum(abs(output)**2*(d['charged'][0]+left[:,None]+right[None,:])))
  tv=(float(np.abs(prob-ideal).sum())+failure)/2;assert abs(prob.sum()+failure-1)<2e-13
  rows.append({'atmospheric_relativistic_phase':phase,'source_parameter':float(t),
   'complete_active_sterile_neutrino_by_antineutrino_records':prob.tolist(),'production_filter_failure':failure,
   'original_probability_total':float(prob.sum()+failure),'vacuum_oscillation_comparison_records':ideal.tolist(),
   'complete_probability_TV_difference_diagnostic':tv,'original_selected_free_energy_before':before,
   'original_selected_free_energy_after':after,'original_mean_detector_work':after-before,
   'neutrino_e_detection_probability':float(prob[0].sum()),'antineutrino_e_detection_probability':float(prob[:,0].sum()),
   'no_postselection_renormalization':True})
 # Uniform exact bound for all phases |phi|<=64pi. The mass/box bounds
 # are verified from rational calibration quantities in the proof.
 eta=F(1,4)*F(5,10**10)**2/F(49,100)**2
 from .calibration import sqrt_bounds
 phase_error=F(1,5*10**15);assert F(128*22,7)*eta<phase_error
 amplitude=8*sqrt_bounds(eta)[1]+F(1,10**14)
 assert amplitude<F(1,2*10**8)
 assert max(r['complete_probability_TV_difference_diagnostic']for r in rows)<float(amplitude)
 return {'source_flavours_neutrino_antineutrino':list(source_flavours),'declared_protocol':'native current instrument, selected charged recoil, free H0 flight with G switched off, then native active-coordinate absorption ports; quantum sterile and preparation complements retained',
  'phase_domain_absolute_upper':'64*pi','uniform_exact_complete_record_TV_bound':str(amplitude),
  'sterile_chirality_weight_upper':str(eta),'source_dispersion_phase_error_upper':str(phase_error),
  'full288_CAR_source_and_filled_reference_retained':True,'reference_port_charge_and_lepton_number_retained':True,
  'physical_baseline_clock_or_beam_flux_calibrated':False,'floating_point_diagnostic_error_included_in_analytic_bound':False,
  'complete_records':rows,
  'CP_comparison_at_pi':{'nu_mu_to_e':float(vacuum_kernel(d,np.pi)[0,1]),'antinu_mu_to_e':float(vacuum_kernel(d,np.pi,True)[0,1])}}

def certificate():
 d=construction();labels,W=complete_vertex(d);free=free_source_check(d);records=[]
 for spin in range(2):
  probabilities=abs(W[:,spin])**2;energies=np.array([d['charged'][0]+np.hypot(d['q'],d['neutral'][i])+np.hypot(d['q'],d['neutral'][j])for o,i,j,se,sn,sbar in labels])
  records.append({'parent_spin':spin,'all144_current_channel_probabilities':probabilities.tolist(),'total_probability':float(probabilities.sum()),'source_mean_current_apparatus_work':float(probabilities@energies-1)})
 labels_bar,Wbar=complete_antiparticle_vertex(d);normal=complete_flight(d,labels,W);anti=complete_flight(d,labels_bar,Wbar,(0,1))
 cp=normal['complete_records'][2]['neutrino_e_detection_probability']-anti['complete_records'][2]['antineutrino_e_detection_probability'];assert cp>.01
 G=F(1,64);t=F(1,10**8);norm=sum(abs(c)for c in interaction()[0].values());assert norm==41952
 error=(G+G*G*norm*norm)*t*t/2;lead=G*t;assert error<lead
 return {'schema':'oph.oscillation_calibrated_native_Dirac.v1','verdict':'PASS_CALIBRATED_NEUTRAL_TENSOR_NATIVE_CURRENT_AND_CONTROLLED_FLIGHT',
  'rigid_charged_ratio_obstruction':rigid_obstruction(),'calibrated_neutral_catalogue':catalogue(),
  'full_source':free,'exact_native_spin_contractions':exact_spin_contraction(),'native_gauge_identities':ward(),
  'native_current_words':len(interaction()[0]),'native_current_tensor_unchanged':True,
  'neutral_family_tensor_is_a_declared_calibration':True,'source_momentum_box':d['q'],
  'all144_channel_labels':[list(x)for x in labels],'actual_current_real':W.real.tolist(),'actual_current_imag':W.imag.tolist(),
  'complete_current_records':records,'complete_flight':normal,
  'exact_antiparticle_spin_contractions':exact_antiparticle_spin_contraction(),
  'actual_antiparticle_current_real':Wbar.real.tolist(),'actual_antiparticle_current_imag':Wbar.imag.tolist(),
  'complete_antiparticle_parent_records':[{'parent_spin':p,'all144_current_channel_probabilities':(abs(Wbar[:,p])**2).tolist(),'total_probability':float(np.sum(abs(Wbar[:,p])**2)),'source_mean_current_apparatus_work':float((abs(Wbar[:,p])**2)@energies-1)}for p in range(2)],
  'complete_antiparticle_flight':anti,
  'actual_native_mu_neutral_vs_antimu_antineutral_e_detection_contrast_at_pi':cp,
  'autonomous_unfiltered_short_time_control':{'G':str(G),'source_parameter':str(t),'full_interaction_norm_upper':str(norm),
   'state_amplitude_error_upper':str(error),'three_body_probability_interval':[str((lead-error)**2),str((lead+error)**2)],
   'full_source_complement_retained':True,'long_time_contact_propagation_claimed':False},
  'scope':'Six oscillation inputs calibrate one native Dirac family tensor. Exact Ward/current algebra and the full source dispersion generate complete finite observer records. The comparison with those six inputs is a calibration; physical beam, matter effects, baseline apparatus and event-rate fits are outside this calculation.'}
