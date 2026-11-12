"""Dirac neutral completion and all native active/sterile current channels."""
from functools import lru_cache
from itertools import product
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from model.openwave_weak_decay.decay import interaction,mass_data,labels
from model.openwave_current_averaging.averaging import act
from model.native_flavour_transport.transport import source,word_unitary,numeric
from model.common_source import calibrated_matter as internal

EPS=F(1,2**32)
def ix(p,spin,component,family):return (p+1)*96+spin*48+component*3+family

def bog(momentum,sign,mass):
 E=np.hypot(momentum,mass);z=-sign*momentum
 if z>=0:u=np.sqrt((E+z)/(2*E));v=sign*mass/np.sqrt(2*E*(E+z))
 else:u=mass/np.sqrt(2*E*(E-z));v=sign*np.sqrt((E-z)/(2*E))
 return E,u,v

def native_frame():
 e,U=eigh(numeric(source()['y']));R=U.conj().T@numeric(word_unitary((2,3,10)))@U
 assert np.linalg.norm(R.conj().T@R-np.eye(3))<1e-13
 return R

def construction(epsilon=float(EPS)):
 m,M,_=mass_data();charged=np.array([m,1.,M]);mu=epsilon*charged;delta=1-m
 assert mu[0]+mu[1]<delta
 q=np.sqrt((delta**2-(mu[0]+mu[1])**2)*(delta**2-(mu[0]-mu[1])**2))/(2*delta)
 R=native_frame();mapping={};energy=np.zeros(288)
 for p,s,c,f in product((-1,0,1),range(2),range(16),range(3)):
  mass=charged[f]if c in(9,11)else mu[f]if c in(10,15)else 0.
  energy[ix(p,s,c,f)]=np.hypot(p*q,mass)
 for p,s in product((-1,0,1),range(2)):
  sign=1-2*s;momentum=p*q
  for f in range(3):
   for component,partner,mass in((11,9,charged[f]),(10,15,mu[f])):
    E,u,v=bog(momentum,sign,mass)
    a=ix(p,s,component,f);b=ix(-p,1-s,partner,f)
    energy[a]=energy[b]=E
    if component==11:
     mapping[a]=[(a,False,u),(b,True,-v)];mapping[b]=[(b,False,u),(a,True,v)]
    else:
     mapping[b]=[(b,False,u),(a,True,v)]
  # Native active coordinates are the full unitary family rotation.
  for nativefamily in range(3):
   row=[]
   for f,mass in enumerate(mu):
    E,u,v=bog(momentum,sign,mass);a=ix(p,s,10,f);b=ix(-p,1-s,15,f)
    row.extend(((a,False,R[nativefamily,f]*u),(b,True,-R[nativefamily,f]*v)))
   mapping[ix(p,s,10,nativefamily)]=row
  for component in range(16):
   if component in(9,10,11,15):continue
   for f in range(3):a=ix(p,s,component,f);mapping[a]=[(a,sign*momentum>0,1.)]
 assert len(mapping)==288
 return dict(epsilon=epsilon,charged=charged,neutral=mu,q=q,R=R,mapping=mapping,energies=energy)

@lru_cache(None)
def leptonic_words():
 out=[]
 for key,c in interaction()[0].items():
  if all(labels(i)[2]in(10,11)for i in key):out.append((tuple(ix(*labels(i))for i in key),float(c)))
 return out

def channel_amplitude(d,i,j,se,sn,sbar,parent_spin=1,orientation=1):
 parent=1<<ix(0,parent_spin,11,1);target=(1<<ix(0,se,11,0))|(1<<ix(orientation,sn,10,i))|(1<<ix(-orientation,sbar,15,j));changed={ix(0,parent_spin,11,1),ix(0,se,11,0),ix(orientation,sn,10,i),ix(-orientation,sbar,15,j)}
 out=0j;count=0
 for(a,b,c,e),coefficient in leptonic_words():
  choices=[]
  for k,dagger in((a,True),(b,True),(e,False),(c,False)):
   row=[(mode,create!=dagger,np.conj(z)if dagger else z)for mode,create,z in d['mapping'][k]if mode in changed]
   if not row:break
   choices.append(row)
  if len(choices)!=4:continue
  for term in product(*choices):
   k,sign=act(parent,tuple((mode,create)for mode,create,z in term))
   if k==target and sign:
    count+=1;out+=coefficient*sign*np.prod([z for mode,create,z in term])
 return out,count

def complete_vertex(d):
 labels=list(product((-1,1),range(3),range(3),range(2),range(2),range(2)));W=np.zeros((len(labels),2),complex)
 for k,(o,i,j,se,sn,sbar)in enumerate(labels):
  for parent_spin in range(2):W[k,parent_spin]=channel_amplitude(d,i,j,se,sn,sbar,parent_spin,o)[0]
 gram=W.conj().T@W;assert np.linalg.norm(gram-np.eye(2))<3e-13
 assert all(abs(W[k,p])<1e-14 for k,(o,i,j,se,sn,sbar)in enumerate(labels)for p in range(2)if sbar!=p or se==sn)
 expected=np.zeros_like(W)
 for k,(o,i,j,se,sn,sbar)in enumerate(labels):
  for p in range(2):
   if sbar!=p or se==sn:continue
   u=bog(o*d['q'],1-2*sn,d['neutral'][i])[1]
   v=bog(o*d['q'],1-2*(1-p),d['neutral'][j])[2]
   expected[k,p]=-o*(-1)**(p+se)*u*v*d['R'][1,i].conj()*d['R'][0,j]
 assert np.linalg.norm(W-expected)<3e-13
 return labels,W

@lru_cache(None)
def exact_spin_contraction():
 # Treat u,v and the two family coefficients as independent scalar factors.
 # The charged rest coefficients contribute their exact product1/2.
 rows=[]
 for o,p,se,sn,sbar in product((-1,1),range(2),range(2),range(2),range(2)):
  parent=1<<ix(0,p,11,1);target=(1<<ix(0,se,11,0))|(1<<ix(o,sn,10,1))|(1<<ix(-o,sbar,15,0))
  mapping={ix(0,p,11,1):(ix(0,p,11,1),False,1),ix(0,se,11,0):(ix(0,se,11,0),False,1),
    ix(o,sn,10,1):(ix(o,sn,10,1),False,1),ix(o,1-sbar,10,0):(ix(-o,sbar,15,0),True,-1)}
  coefficient=F(0)
  for (a,b,c,e),z in leptonic_words():
   if any(k not in mapping for k in(a,b,c,e)):continue
   ops=[];factor=1
   for k,dagger in((a,True),(b,True),(e,False),(c,False)):
    mode,create,sgn=mapping[k];ops.append((mode,create!=dagger));factor*=sgn
   output,sign=act(parent,tuple(ops))
   if output==target:coefficient+=F(z)*sign*factor/2
  expected=F(-o*(-1)**(p+se))if sbar==p and se!=sn else F(0)
  assert coefficient==expected
  rows.append({'orientation':o,'parent_spin':p,'charged_spin':se,'neutral_spin':sn,'antineutral_spin':sbar,'exact_native_CAR_coefficient':str(coefficient)})
 return rows

def free_source_check(d):
 A=np.zeros((288,288),complex);B=A.copy()
 for original,terms in d['mapping'].items():
  for mode,create,z in terms:(B if create else A)[original,mode]+=z
 car=max(np.linalg.norm(A@A.conj().T+B@B.conj().T-np.eye(288)),np.linalg.norm(A@B.T+B@A.T))
 W=np.block([[A,B],[B.conj(),A.conj()]]);H=(W*np.r_[d['energies'],-d['energies']])@W.conj().T
 h=np.diag([-float(1-2*s)*p*d['q']for p,s,c,f in product((-1,0,1),range(2),range(16),range(3))]);pair=np.zeros_like(A)
 for p,s,f in product((-1,0,1),range(2),range(3)):
  sign=1-2*s;pair[ix(p,s,11,f),ix(-p,1-s,9,f)]=sign*d['charged'][f]
  pair[ix(-p,1-s,9,f),ix(p,s,11,f)]=-sign*d['charged'][f]
  for j in range(3):
   pair[ix(p,s,10,j),ix(-p,1-s,15,f)]=sign*d['R'][j,f]*d['neutral'][f]
   pair[ix(-p,1-s,15,f),ix(p,s,10,j)]=-sign*d['R'][j,f]*d['neutral'][f]
 discrepancy=np.linalg.norm(H-np.block([[h,pair],[-pair.conj(),-h.T]]));assert car<2e-12 and discrepancy<2e-12
 return {'full_CAR_modes':288,'complete_Nambu_dimension':576,'CAR_numeric_residual':car,'declared_Weyl_Dirac_BdG_numeric_residual':discrepancy,'massless_zero_momentum_CAR_coordinates':int(np.sum(d['energies']==0)),'all_three_momentum_spectator_coordinates_retained':True}

def ward():
 import sympy as sy
 from model.common_source.chiral_candidate import action
 L=internal._LH_tensor();J=internal.weak_raising_current(families=False);T=[(J+J.H)/2,(J-J.H)/(2*sy.I),(J*J.H-J.H*J)/2];higgs=[]
 for t in T:
  S=-(L.H*L).inv()*L.H*t.T*L;assert t.T*L+L*S==sy.zeros(15,2);higgs.append(S)
 Y=sy.diag(*[sy.Rational(r['q6Y'])/6 for r in action.ordering()['all_left_components']]);assert Y.T*L+L*sy.eye(2)/2==sy.zeros(15,2)
 q=action.electric_charge().extract(range(0,45,3),range(0,45,3));assert q.T*(L*sy.Matrix([0,1]))==sy.zeros(15,1)
 return {'actual_native_weak_Higgs_intertwiner_identities':3,'actual_native_hypercharge_Higgs_identity':True,'added_Weyl_hypercharge_colour_weak_generators':'0','gauge_anomaly_polynomials_unchanged':'the added neutral module contributes0 to every gauge/mixed gauge anomaly','exact_conserved_Dirac_lepton_number':'active L:+1, charged/neutral singlet conjugates:-1; every Yukawa pair and native current has total0','Majorana_term':'absent by declared conserved Dirac lepton number'}

def flight(d,labels,W):
 # Parent spin down, positive neutral momentum and the spin-up charged recoil.
 # The coherent selected block retains all neutral and antineutral masses.
 A=np.zeros((3,3),complex);production_failure=0.
 for k,(o,i,j,se,sn,sbar)in enumerate(labels):
  if(o,se,sn,sbar)==(1,0,1,1):A[i,j]=W[k,1]
  else:production_failure+=float(abs(W[k,1])**2)
 production=float(np.linalg.norm(A)**2);mu=d['neutral'];E=np.hypot(d['q'],mu);differences=(mu*mu-mu[0]**2)/(E+E[0]);scale=differences[-1];assert scale>0
 u=np.array([bog(d['q'],-1,x)[1]for x in mu]);v=np.array([abs(bog(d['q'],1,x)[2])for x in mu]);
 neutrino=[d['R'][f:f+1,:]*u for f in range(3)]+[np.diag([abs(bog(d['q'],-1,x)[2])for x in mu])]
 antinu=[d['R'][f:f+1,:].conj()*v for f in range(3)]+[np.diag([bog(d['q'],1,x)[1]for x in mu])]
 assert np.linalg.norm(sum(K.conj().T@K for K in neutrino)-np.eye(3))<1e-13
 assert np.linalg.norm(sum(K.conj().T@K for K in antinu)-np.eye(3))<1e-13
 rows=[]
 for phase in(0.,np.pi/2,np.pi,2*np.pi):
  phases=np.exp(-1j*phase*differences/scale);state=phases[:,None]*A*phases[None,:]
  probabilities=np.array([[np.linalg.norm(K@state@J.T)**2 for J in antinu]for K in neutrino])
  assert abs(probabilities.sum()+production_failure-1)<2e-13
  before=float(np.sum(abs(state)**2*(d['charged'][0]+E[:,None]+E[None,:])));after=0.
  for r,K in enumerate(neutrino):
   for s,J in enumerate(antinu):
    output=K@state@J.T;left=E if r==3 else np.zeros(1);right=E if s==3 else np.zeros(1)
    after+=float(np.sum(abs(output)**2*(d['charged'][0]+left[:,None]+right[None,:])))
  rows.append({'heavy_light_phase':phase,'source_time':float(phase/scale),'complete3active_plus_sterile_neutrino_by_antineutrino_probabilities':probabilities.tolist(),'production_filter_failure_probability':production_failure,'original_record_probability_sum':float(probabilities.sum()+production_failure),'original_selected_branch_free_energy_before':before,'original_selected_branch_free_energy_after':after,'original_mean_detector_work':after-before,'no_postselection_renormalization':True})
 return {'controlled_protocol':'same288-CAR H0, same neutral tensor and native current; source/current filter then free propagation with G switched off and declared native active-coordinate annihilation-port instruments','retained_charged_recoil':'light rest spin up; other recoil/spin/orientation records are the production complement','production_selected_probability':production,'three_active_detector_rows_and_quantum_sterile_failure_each_particle':True,'detector_reference':'An active annihilation success leaves the source vacuum in that one-particle factor and transfers its conserved lepton number to a retained port reference of declared zero energy; the sterile failure keeps the mass-state quantum factor. References retain+1 for a neutrino absorption and-1 for an antineutrino absorption. No physical charged-current apparatus or energy calibration is inferred.','complete_records':rows,'physical_PMNS_fit_claimed':False,'long_time_unfiltered_contact_propagation_claimed':False}

def conditional_masses():
 from decimal import Decimal,localcontext
 from model.lepton_mass_compiler.compiler import electron_calibration
 with localcontext()as ctx:
  ctx.prec=40;data=electron_calibration();eps=Decimal(1)/Decimal(2**32);rows=[]
  for f,r in enumerate(data['rows']):
   mass=eps*Decimal(r['electron_anchored_GeV']);rows.append({'family':f+1,'conditional_neutrino_and_antineutrino_GeV':str(mass),'conditional_eV':str(mass*Decimal(10**9)),'helicities_per_particle':2,'independent_sterile_mass_branch':False})
  return {'sole_empirical_energy_anchor_GeV':data['rows'][0]['electron_anchored_GeV'],'declared_neutral_to_charged_mass_scale':str(EPS),'rows':rows,'active_and_sterile_chiralities_share_each_Dirac_mass':True,'source_scale_or_mixing_fitted_to_neutrino_data':False,'physical_neutrino_mass_identity_or_oscillation_fit_claimed':False,'inherited_lepton_phase_and_comparison_ancestry':'phase2/9 and existing electron calibration; retrospective charged-lepton comparison is retained'}

def certificate():
 spin_coefficients=exact_spin_contraction()
 results=[]
 for epsilon in(float(EPS),.25):
  d=construction(epsilon);labels,W=complete_vertex(d);free=free_source_check(d);records=[]
  for parent_spin in range(2):
   probabilities=abs(W[:,parent_spin])**2;energies=np.array([d['charged'][0]+np.hypot(d['q'],d['neutral'][i])+np.hypot(d['q'],d['neutral'][j])for o,i,j,se,sn,sbar in labels])
   records.append({'prepared_parent_spin':parent_spin,'all144_native_channel_probabilities':probabilities.tolist(),'outside_parent_preparation_failure_probability':0.,'complete_record_sum':float(probabilities.sum()),'free_source_mean_apparatus_work':float(probabilities@energies-1)})
  results.append({'neutral_scale':epsilon,'full_free_source':free,'source_box_q':d['q'],'massive_neutral_source_energies':d['neutral'].tolist(),'channel_labels':[list(x)for x in labels],'actual_native_vertex_real':W.real.tolist(),'actual_native_vertex_imag':W.imag.tolist(),'native_current_intertwiner_Gram_real':(W.conj().T@W).real.tolist(),'complete_prepared_parent_records':records,'all_free_kinematic_channels_open':[[bool(d['neutral'][i]+d['neutral'][j]<1-d['charged'][0])for j in range(3)]for i in range(3)],'closed_free_channels_removed_from_current_instrument':False,'complete_joint_flight_records':flight(d,labels,W)if epsilon==float(EPS)else None})
 G=F(1,64);t=F(1,10**8);L=sum(abs(c)for c in interaction()[0].values());assert L==41952
 error=(G+G*G*L*L)*t*t/2;lead=G*t;assert error<lead
 return {'schema':'oph.native_Dirac_neutral_current_join.v1','verdict':'PASS_COMPLETE_DIRAC_NEUTRAL_SOURCE_CURRENT_AND_RECORDS','source':'complete native contact embedded in288 CAR coordinates, with one gauge-neutral Weyl field per family and a declared Dirac Higgs Yukawa tensor','native_quartic_word_count':len(interaction()[0]),'weak_active_CAR_coordinates':144,'weak_singlet_CAR_coordinates':144,'same_native_current_added_sterile_entries_zero':True,'gauge_and_number_identities':ward(),'source_conventions':{'charged_masses':'same three analyzer ratios, middle charged mass1','neutral_mass_matrix':'epsilon R*diag(m,1,M), R=F^dagger V_(2,3,10)F in the ordered real native Y frame','Higgs_tensor':'native LH tensor times H in the annihilation bilinear, paired with the sterile module and spin epsilon; H=(0,1); the declared creation pair has Delta_(nu_j,N_i)=sigma_s R_ji mu_i','reference':'filled negative charged and neutral Dirac bands, massless zero bands vacant; epsilon-positive neutral zero-band reference is distinct from the massless zero-vacancy convention','current':'exact31464 native weak-contact quartic words, with all original spectators retained','neutral_module_or_scale_uniqueness_required':False},'exact_native_spin_contractions':spin_coefficients,'complete_massive_members':results,'conditional_GeV_catalogue':conditional_masses(),'exact_current_isometry':'K=P_all144out V P_parent; K^dagger K=P_parent from actual spin contractions, family unitarity and u^2+v^2=1. K_failure=I-P_parent. This is a declared current apparatus, not installed source decay dynamics.','autonomous_unfiltered_short_time_control':{'G':str(G),'source_time':str(t),'primary_all_target_free_detuning_bound':'1','full_native_interaction_norm_upper':'41952','state_amplitude_error_upper':str(error),'total_three_body_probability_interval':[str((lead-error)**2),str((lead+error)**2)],'complete_three_body_and_entire_source_complement_retained':True,'finite_time_dyson_proof':'Interaction-picture Duhamel remainder <=G^2 L^2 t^2/2; detuning of the exact isometry contributes <=G t^2/2. No positive-frequency vertex replacement.'},'no_physical_weak_rate_PMNS_fit_or_full_particle_dictionary_claimed':True,'scope':'Complete Dirac spectrum/current join, exact native current isometry and controlled observer protocol with conditional mass values. Added neutral module, epsilon, port word, fixed Higgs, preparation/detection and switches are allowed declared inputs. The full unfiltered source is controlled only in the stated short-time comparison; no physical lifetime or sterile discovery is inferred.'}
