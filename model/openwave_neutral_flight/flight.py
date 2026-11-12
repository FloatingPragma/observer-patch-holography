"""Finite current contractions, recoil and complete same-clock observations."""
from functools import lru_cache
import numpy as np
from model.openwave_neutrino_oscillation.oscillation import construction

REFERENCE_COST=1/4
REFERENCE_STIFFNESS=1/8

def index(i,n,k,L):return (i*L+n)*L+k

def data(L=5,zero_phase=False):
 if L<3 or L%2!=1:raise ValueError('an odd periodic side at least three is required')
 native=construction();mass=native['masses'];C=np.eye(3,dtype=complex)if zero_phase else native['C']
 theta=2*np.pi*np.arange(L)/L;p=np.sin(theta);E=np.sqrt(mass[:,None]**2+p[None,:]**2);w=(E+p[None,:])/(2*E)
 er=REFERENCE_STIFFNESS*(1-np.cos(theta));A=np.zeros((3,3*L),complex)
 for i in range(3):A[i,i*L:(i+1)*L]=np.sqrt(w[i])/np.sqrt(L)
 assert np.max(abs(A@A.conj().T-.5*np.eye(3)))<1e-14
 K=np.zeros((3*L*L,3*L*L),complex)
 for i in range(3):
  for a in range(3):
   for n in range(L):
    for m in range(L):
     coefficient=np.conj(C[a,i])*np.sqrt(w[i,n]*w[a,m])/L
     for k in range(L):K[index(i,n,(k+m-n)%L,L),index(a,m,k,L)]=coefficient
 gram=K.conj().T@K;assert np.max(abs(gram@gram-gram/4))<1e-14
 failure=np.eye(len(K))+(2*np.sqrt(3)-4)*gram
 assert np.max(abs(gram+failure@failure-np.eye(len(K))))<2e-14
 # Total momentum is the exact cyclic-group character, without a branch cut.
 Z=np.diag([np.exp(2j*np.pi*(n+k)/L)for i in range(3)for n in range(L)for k in range(L)])
 assert np.max(abs(Z@K-K@Z))<2e-14
 energy=(E[:,:,None]+er[None,None,:]).reshape(-1)
 return dict(L=L,C=C,mass=mass,E=E,w=w,er=er,A=A,K=K,F=failure,energy=energy,
             native_current_frame_residual=native['frame_residual'],native_Ward_identities=native['ward'])

def recoil(L,kind):
 if kind=='localized':return np.ones(L,complex)/np.sqrt(L)
 if kind=='sharp':return np.eye(L,dtype=complex)[0]
 if kind=='partial':
  angle=2*np.pi*np.arange(L)/L;v=np.exp(2*np.cos(angle)).astype(complex);return v/np.linalg.norm(v)
 raise ValueError('unknown declared recoil preparation')

def receiver(d,x):
 L=d['L'];phase=np.exp(2j*np.pi*np.arange(L)*x/L);A=d['A']*np.tile(phase,3)[None,:]
 D=A.conj().T@d['C']@A
 gram=D.conj().T@D
 assert np.max(abs(gram@gram-gram/4))<2e-14
 F=np.eye(3*L)+(2*np.sqrt(3)-4)*gram
 return D,F

def mean(v,h):return float(np.vdot(v,h*v).real)

def experiment(d,initial_family,time,x,recoil_kind='localized',mass_tag=False):
 L=d['L'];r=recoil(L,recoil_kind);initial=np.zeros(3*L*L,complex)
 initial[initial_family*L*L:initial_family*L*L+L]=r
 seed=d['K']@initial;failed=d['F']@initial;w=float(np.vdot(seed,seed).real)
 assert abs(w-1/(4*L))<1e-14
 E=d['energy'];Hs=E+REFERENCE_COST
 original_energy=mean(initial,E)
 source_work=mean(seed,Hs)+mean(failed,E)-original_energy
 D,F=receiver(d,x);D=np.kron(D,np.eye(L));F=np.kron(F,np.eye(L));u=np.exp(-1j*time*Hs)
 seeds=[]
 if mass_tag:
  for i in range(3):v=seed.copy();v[:i*L*L]=0;v[(i+1)*L*L:]=0;seeds.append(v)
 else:seeds=[seed]
 full_probability=[];receive=np.zeros(3*L*L);fail=np.zeros(3*L*L);probe_work=0.;energy_error=0.
 for v in seeds:
  evolved=u*v;out=D@evolved;complement=F@evolved
  pr=abs(out)**2;pf=abs(complement)**2;receive+=pr;fail+=pf
  full_probability.append({'received':pr.tolist(),'detector_failure':pf.tolist()})
  # Detector recoil stays in its prepared position eigenstate. Its kinetic
  # mean is unchanged; a received charged output pays its positive charge cost.
  probe_work+=mean(out,Hs+REFERENCE_COST)+mean(complement,Hs)-mean(evolved,Hs)
  energy_error=max(energy_error,abs(mean(evolved,Hs)-mean(v,Hs)))
 source_failure=abs(failed)**2;total=float(source_failure.sum()+receive.sum()+fail.sum());assert abs(total-1)<3e-14
 byfamily=receive.reshape(3,L,L).sum(axis=(1,2))
 return {'initial_family':initial_family,'source_time':time,'detector_site':x,'source_recoil_preparation':recoil_kind,
  'mass_record_retained':mass_tag,'source_success_probability':w,
  'received_by_charged_family':byfamily.tolist(),'total_received_probability':float(receive.sum()),
  'complete_success_records_by_retained_mass_tag':full_probability,
  'complete_source_failure_probabilities':source_failure.tolist(),
  'record_axes':'charged/neutral family, source momentum, source-recoil momentum; mass tag is additional when supplied',
  'complete_probability_sum':total,'source_preparation_mean_energy_work':source_work,
  'detector_mean_energy_work':probe_work,'free_source_energy_error':energy_error,
  'free_initial_source_plus_recoil_energy':original_energy,
  'source_and_detector_charge_reference_cost_each':REFERENCE_COST}

def certificate():
 d=data(5);rows=[];profiles={}
 for kind in('localized','sharp','partial'):
  for t in(0.,3.,9.):
   group=[experiment(d,1,t,x,kind)for x in range(5)];rows+=group
   profile=[v['total_received_probability']for v in group];profiles[f'{kind}:{t}']=profile
   if kind=='sharp':assert np.ptp(profile)<2e-15
 assert np.ptp(profiles['localized:0.0'])>1e-3
 # Complete records for all three input families and a retained mass-record control.
 rows += [experiment(d,a,7.,2,'localized',tag)for a in range(3)for tag in(False,True)]
 tagged=[sum(experiment(d,1,t,x,'localized',True)['received_by_charged_family'][b]for x in range(5))for t in(0.,9.)for b in range(3)]
 assert np.max(abs(np.array(tagged[:3])-tagged[3:]))<2e-14
 coherent=[sum(experiment(d,1,t,x,'localized')['received_by_charged_family'][b]for x in range(5))for t in(0.,9.)for b in range(3)]
 contrast=float(np.max(abs(np.array(coherent[:3])-coherent[3:])));assert contrast>1e-8
 zero=data(5,True);control=[experiment(zero,a,9.,2)for a in range(3)]
 assert all(max(v['received_by_charged_family'][b]for b in range(3)if b!=a)<1e-25 for a,v in enumerate(control))
 d7=data(7);volume=[experiment(d7,1,3.,x)for x in range(7)]
 return {'schema':'oph.complete_native_current_neutral_flight.v1',
  'native_internal_carrier':'15 internal times3families times2spin on each site; complete quadratic Fock space, filled negative bands and retained zero-mode convention',
  'source_action':'same native charged-Y and neutral V*YV^dagger Higgs kernel, homogeneous periodic Weyl kinetic symbol sin(p); exact one-quasiparticle positive-band current contractions',
  'sites':5,'physical_CAR_modes':450,'active_one_quasiparticle_recoil_dimension':75,
  'source_family_masses':d['mass'].tolist(),'actual_native_weak_current_frame_real':d['C'].real.tolist(),'actual_native_weak_current_frame_imag':d['C'].imag.tolist(),
  'native_current_frame_residual':d['native_current_frame_residual'],'native_LH_Ward_identities':d['native_Ward_identities'],
  'exact_aperture_identity':'K^dagger K=(1/4)P; F=sqrt(I-K^dagger K)=I+(2sqrt(3)-4)K^dagger K',
  'exact_source_success':'1/(4L) for every rest-family input and every normalized recoil register',
  'exact_total_momentum_conservation':'source momentum plus charged-reference momentum modulo the reciprocal period',
  'complete_experiments':rows,'arrival_profiles':profiles,
  'profile_normalization':'Each x is a complete separately chosen detector aperture; uniform x selection divides every row probability by L and retains x on success and both failures. The yes profile alone is not a normalized position law.',
  'coherent_mass_phase_contrast_in_spatially_summed_family_records':contrast,
  'mass_tag_control_spatially_summed_flavour_difference':float(np.max(abs(np.array(tagged[:3])-tagged[3:]))),
  'zero_port_phase_control_has_zero_other_family_probability':True,
  'side_seven_complete_records':volume,
  'arbitrary_reference_domain':'the exact Stinespring instruments on the admitted one-quasiparticle/reference sector, extended by the identity no-operation outcome on the full Fock complement',
  'scope':'Current/spectral apparatus, free Higgs background, finite spatial regulator and charged recoil reservoirs are declared. The same source produces complete flight/mixing records. No autonomous weak decay, observed PMNS fit, neutrino GeV scale, or transfer of the separate PW119 composite is asserted.',
  'verdict':'PASS_COMMON_NATIVE_CURRENT_PRODUCTION_FLIGHT_AND_RECOIL_RECORDS'}
