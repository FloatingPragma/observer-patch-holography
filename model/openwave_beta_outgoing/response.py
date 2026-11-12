from functools import lru_cache
from pathlib import Path
import json
from itertools import product
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from model.openwave_axial_state.source import sector
from model.openwave_baryon_weak_formfactors.native import currents,quartet,PAULI
from model.openwave_neutrino_calibration.calibration import calibration
from model.openwave_dirac_neutral_current.neutral import bog
from model.openwave_common_mass_source.source import local_norm_certificate

G=F(1,64)
@lru_cache(None)
def source():
 d=sector('baryon');W,_=quartet(d['h']);n=W[:,0];p=W[:,2:];E,U=eigh(d['h']);En=float(np.vdot(n,d['h']@n).real);Ep=float(np.vdot(p[:,0],d['h']@p[:,0]).real);cur=currents()['local'];T={};P={}
 for kind in('plus','minus','third'):
  T[kind]=np.zeros((2,4,120),complex);P[kind]=np.zeros((2,4,2),complex)
  for x,mu in product(range(2),range(4)):
   J=cur[x][mu][0]if kind=='plus'else cur[x][mu][0].conj().T if kind=='minus'else cur[x][mu][1]
   T[kind][x,mu]=U.conj().T@J@n;P[kind][x,mu]=p.conj().T@J@n
 c=calibration();charged=np.asarray(c['charged'])/c['charged'][0];neutral=np.asarray(c['neutral'])/c['charged'][0]
 assert max(c['source_mass_squared_exact'])/c['charged_exact'][0]**2<F(1,10**7)**2
 axial=json.loads((Path(__file__).resolve().parents[1]/'openwave_axial_state/results.json').read_text())
 alo,ahi=map(F,axial['current_ground_certificate']['native_ground_axial_interval_rational']);assert F('1.2752')<alo<ahi<F('1.2754')
 ga=float(abs(sum(P['plus'][x,3,0]for x in range(2))))
 assert 1.2752<ga<1.2754
 # Every actual local proton current is half the unchanged global current.
 for x in range(2):
  assert np.linalg.norm(P['plus'][x,0]-np.array([.5,0]))<2e-10
  for mu in range(1,4):assert np.linalg.norm(P['plus'][x,mu]-.5*ga*PAULI[mu][:,0])<2e-10
 return dict(native=d,W=W,initial=n,proton=p,E=E,U=U,En=En,Ep=Ep,T=T,P=P,charged=charged,neutral=neutral,R=c['R'],ga=ga)

def grid(L):
 assert L>=8 and L%8==0
 theta=2*np.pi*np.arange(L)/L;p=2*np.sin(theta);p[0]=p[L//2]=0
 return theta,p

def pair_response(L,t):
 d=source();theta,mom=grid(L);phase=np.exp(1j*(theta[None,:]-theta[:,None]));rows=[];total=0.;normV=0.;selected=np.zeros((L,L));selected_by_neutral=[];energy=0.;by_type={};spectral=np.zeros(120)
 channels=[]
 for f,i in product(range(3),repeat=2):
  channels.append(('electron_antineutral',f,i,'plus',d['charged'][f],d['neutral'][i],d['R'][f,i]))
  channels.append(('neutral_positron',i,f,'minus',d['neutral'][i],d['charged'][f],d['R'][f,i].conjugate()))
 for f in range(3):
  channels.append(('electron_positron',f,f,'third',d['charged'][f],d['charged'][f],-1.))
  channels.append(('neutral_antineutral',f,f,'third',d['neutral'][f],d['neutral'][f],1.))
 for name,f,i,kind,ma,mb,family in channels:
  probability=0.;zero_column=0.;protonmatrix=np.zeros((L,L));Es=np.hypot(mom,ma);Et=np.hypot(mom,mb)
  delta=d['E'][:,None,None]-d['En']+Es[None,:,None]+Et[None,None,:];timefactor=float(t)*np.sinc(delta*float(t)/(2*np.pi))*np.exp(-.5j*delta*float(t))
  for s,r in product(range(2),repeat=2):
   ua=np.array([bog(p,1-2*s,ma)[1]for p in mom]);vb=np.array([bog(p,1-2*r,mb)[2]for p in mom]);weights=family*ua[:,None]*vb[None,:]/L
   J=np.array([sum((1 if mu==0 else -1)*d['T'][kind][x,mu]*PAULI[mu][s,r]for mu in range(4))for x in range(2)])
   column=(J[0,:,None,None]+J[1,:,None,None]*phase[None,:,:])*weights[None,:,:]
   probability_array=abs(float(G)*column*timefactor)**2;probability+=float(probability_array.sum());zero_column+=float(np.linalg.norm(column)**2);spectral+=probability_array.sum(axis=(1,2));energy+=float(np.sum(probability_array*(d['E'][:,None,None]+Es[None,:,None]+Et[None,None,:])))
   if name=='electron_antineutral'and f==0:
    jp=np.array([sum((1 if mu==0 else -1)*d['P'][kind][x,mu]*PAULI[mu][s,r]for mu in range(4))for x in range(2)])
    gap=d['Ep']-d['En']+Es[:,None]+Et[None,:];tf=float(t)*np.sinc(gap*float(t)/(2*np.pi))*np.exp(-.5j*gap*float(t));amp=float(G)*(jp[0,:,None,None]+jp[1,:,None,None]*phase[None,:,:])*weights[None,:,:]*tf
    protonmatrix+=np.sum(abs(amp)**2,axis=0)
  if name=='electron_antineutral'and f==0:selected+=protonmatrix;selected_by_neutral.append(protonmatrix.tolist())
  rows.append({'pair_type':name,'particle_family':f,'antiparticle_family':i,'first_order_probability':probability,'exact_first_column_squared_norm_floating':zero_column});total+=probability;normV+=zero_column;by_type[name]=by_type.get(name,0)+probability
 forward=float(selected[mom>0].sum());backward=float(selected[mom<0].sum());zero=float(selected[mom==0].sum());pselected=forward+backward+zero
 assert abs(sum(spectral)-total)<1e-15 and len(rows)==24
 return {'ring_sites':L,'source_time':str(t),'native_CAR_coordinate_count':192+24*(L-2),'active_three_quark_CAR_coordinates':24,'full_Dirac_lepton_CAR_coordinates':24*L,'complete_three_quark_singlet_coordinates':120,'all_first_pair_channels':rows,'first_order_pair_probability':total,'first_order_pair_probability_by_type':by_type,'complete_hadron_spectral_energies':d['E'].tolist(),'complete_hadron_spectral_pair_weights':spectral.tolist(),'first_pair_energy_weighted_sum':energy,'full_native_first_column_squared_norm_floating':normV,'selected_actual_proton_light_electron_all_neutrals':{'positive_momentum':forward,'negative_momentum':backward,'zero_momentum':zero,'odd_difference':forward-backward,'conditional_first_order_asymmetry':(forward-backward)/pselected,'complete_momentum_probability_by_neutral_mass':selected_by_neutral,'particle_momenta':mom.tolist(),'antineutral_physical_momentum':'minus the canonical hole momentum indexed by the second grid coordinate'},'multiple_pair_and_no_pair_and_failed_source_branches_retained_in_exact_complement':True}

def parity_bound(L,t):
 t=F(t);glo=F(1594,1250);ghi=F(6377,5000);mu=F(1,10**7);delta=(G*54*t)**2/2
 # All momenta in [pi/4,3pi/4] have p/E>4/5, at least L/4 modes.
 # The massive neutral correction uses p_min>=8/L for L>=8: sin(2pi/L)>=4/L.
 pmin=F(8,L);mass_correction=24*(G*t)**2*(mu*mu/(4*pmin*pmin)+t*mu)
 # |energy mismatch|<7 in the selected proton/electron/neutral source.
 sinc_lower=1-(7*t/2)**2/3
 leading_lower=(G*t)**2*glo*(glo-1)/5*sinc_lower-mass_correction
 probability_error=(8*G*t+delta)*delta
 lower=leading_lower-probability_error;assert lower>0
 return {'source_time':str(t),'ring_sites':L,'native_full_contact_norm_upper':'54','full_unfiltered_state_remainder_upper':str(delta),'selected_proton_emission_first_order_amplitude_upper':str(4*G*t),'massive_neutral_vs_zero_mass_limit_probability_error_upper':str(mass_correction),'negative_first_order_parity_magnitude_lower':str(leading_lower),'full_unfiltered_odd_probability_error_upper':str(probability_error),'strict_negative_exact_parity_difference_upper':str(-lower),'strict_negative_exact_parity_difference_upper_decimal':-float(lower),'original_input_parity_magnitude_lower':str(lower/120),'all_L_multiple8_in_domain':'8<=L<=64; same proof bounds each finite member; no infinite-volume angular limit is inferred'}

def certificate():
 d=source();norm=local_norm_certificate();rows=[pair_response(L,F(1,10000))for L in(8,16,32,64)];bounds=[parity_bound(L,F(1,10000))for L in(8,16,32,64)]
 for r,b in zip(rows,bounds):
  assert r['selected_actual_proton_light_electron_all_neutrals']['odd_difference']<float(F(b['strict_negative_exact_parity_difference_upper']))
  wick=wick_column_norm(r['ring_sites']);assert abs(wick-r['full_native_first_column_squared_norm_floating'])<2e-9
  r['independent_full_Gaussian_Wick_column_squared_norm']=wick
  err=float(F(b['full_unfiltered_odd_probability_error_upper']));sel=r['selected_actual_proton_light_electron_all_neutrals'];r['complete_original_terminal_record_intervals']={'preparation_failure':[119/120,119/120],**{k:[max(0,sel[k]-err)/120,min(1,sel[k]+err)/120]for k in('positive_momentum','negative_momentum','zero_momentum')},'prepared_quantum_complement':[(1-sum(sel[k]for k in('positive_momentum','negative_momentum','zero_momentum'))-err)/120,(1-sum(sel[k]for k in('positive_momentum','negative_momentum','zero_momentum'))+err)/120]}
 return {'schema':'oph.complete_native_beta_ring_parity.v1','verdict':'PASS_COMPLETE_NATIVE_POLARIZED_PARITY_RESPONSE','source':'full axial-calibrated native baryon action and unfiltered semileptonic cross-current; all charged and calibrated neutral Dirac bands on a periodic ring','declared_dispersion':'p_k=2 sin(2pi k/L); all lattice valleys and both helicities retained','interaction_coefficient':'1/64','active_contact_cells':[0,1],'neutral_calibration_inputs':'the six disclosed neutral/oscillation inputs and chosen lightest mass from openwave_neutrino_calibration; the current is not fitted separately','strong_calibration_inputs':'p,n,charged pion and |gA|; no effective axial current factor','reference':'the full free massive Dirac Gaussian vacuum on the ring; opposite momenta give local active covariance I/2 and zero local weak-current mean; nonlocal covariances are recomputed by the complete pair sums','same_source_mean_initial_energy_GeV':d['En']*.00051099895069,'native_axial_diagnostic':d['ga'],'native_local_norm_certificate':norm,'all_first_pair_types_retained':['electron+antineutral','neutral+positron','electron+positron','neutral+antineutral'],'finite_ring_executions':rows,'full_unfiltered_parity_certificates':bounds,'scalar_kinetic_control':'The filled-neutral scalar kinetic source has exact spatial reflection; its even polarized preparation gives zero positive/negative momentum contrast at every time. The chiral Dirac ring is an explicitly different kinetic/reference member.','complete_records':{'preparation_success':'1/120 on the uniform complete native baryon input and declared Gaussian lepton input','preparation_failure_quantum_branch':'119/120','terminal_instrument':'actual strong proton band with light electron/any neutral pair, resolved by all momenta and helicities; the entire orthogonal full-Fock source is one retained quantum outcome','full_source_multi_pair_channels':'included in the exact unitary and bounded remainder, never truncated in the action','source_energy':'conserved under the autonomous full Hamiltonian','nonselective_preparation_mean_work':'0','terminal_readout_absolute_work_upper':'2 G ||V|| <=27/16 because this spectral/momentum instrument commutes with H0; the full-current work is retained'},'physical_neutron_lifetime_or_three_dimensional_angular_coefficient_claimed':False,'scope':'An actual parity-odd finite-time native observer statistic, with complete-source remainder and finite-volume numerical comparison. The source time and G are declared units, not a calibrated physical weak lifetime; autonomous finite-time response is not an irreversible released-particle theorem.'}


def wick_column_norm(L):
 """Independent full Gaussian current covariance, including off-site entries."""
 d=source();theta,mom=grid(L);rho=np.zeros((24,24),complex)
 def ix(x,s,w,f):return ((2*x+s)*2+w)*3+f
 for x,y,s,w in product(range(2),range(2),range(2),range(2)):
  masses=d['neutral']if w==0 else d['charged'];diag=[]
  for mass in masses:diag.append(sum(np.exp(1j*(x-y)*theta[k])*bog(mom[k],1-2*s,mass)[2]**2 for k in range(L))/L)
  block=d['R']@np.diag(diag)@d['R'].conj().T if w==0 else np.diag(diag)
  for f,j in product(range(3),repeat=2):rho[ix(x,s,w,f),ix(y,s,w,j)]=block[f,j]
 J=[];Q=[];factors=[];cur=currents()['local'];weak=(np.array([[0,1],[1,0]])/2,np.array([[0,-1j],[1j,0]])/2,np.diag([1.,-1.])/2)
 for x,mu,a in product(range(2),range(4),range(3)):
  matrix=np.zeros((24,24),complex)
  for s,r,w,v,f in product(range(2),range(2),range(2),range(2),range(3)):matrix[ix(x,s,w,f),ix(x,r,v,f)]=PAULI[mu][s,r]*weak[a][w,v]
  qp,qz=cur[x][mu];q=((qp+qp.conj().T)/2,(qp-qp.conj().T)/(2j),qz)[a];J.append(matrix);Q.append(q@d['initial']);factors.append(2*(1 if mu==0 else -1))
 assert np.max(abs(np.array([np.trace(j@rho)for j in J])))<2e-14
 total=0j
 for a,b in product(range(24),repeat=2):total+=factors[a]*factors[b]*np.vdot(Q[a],Q[b])*np.trace(J[a]@(np.eye(24)-rho)@J[b]@rho)
 assert abs(total.imag)<1e-10
 return float(total.real)
