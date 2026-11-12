"""Finite spectral measures with exact dyadic frame and complete records."""
from fractions import Fraction as F
import numpy as np
from .source import construction,exact_parity_matrix,parity_basis,G,mass_only
from model.openwave_six_quark.certificate import spectral_enclosure

def measure(member,charge):
 d=construction(member,charge);H,S=exact_parity_matrix(d);cert=spectral_enclosure(H);D=cert['dyadic_denominator'];O=np.array([[int(x)for x in row]for row in cert['dyadic_eigenframe']],dtype=object);E=[F(int(x),D)for x in cert['dyadic_eigenvalues']];epsilon=F(cert['all_eigenvalue_error_upper']);eta=F(cert['exact_Gram_error_upper'])
 B=parity_basis(d);a=np.asarray(B.T@d['initial']).ravel();v=np.array([int(round(float(x.real)*2**48))for x in a],dtype=object);N=sum(x*x for x in v);assert N>0
 vector=np.array(v,dtype=float)/np.sqrt(float(N));exact_state=B@vector;deviation=float(np.linalg.norm(exact_state-d['initial']));assert deviation<1e-11
 overlap=O.T@v;weights=[F(int(x*x),int(D*D*N))for x in overlap];assert abs(sum(weights)-1)<=eta
 groups=[];start=0
 for i in range(1,len(E)+1):
  if i==len(E)or E[i]-E[i-1]>8*epsilon:groups.append(list(range(start,i)));start=i
 rows=[];unit=F(mass_only.INPUTS['electron_unit_GeV'])
 for indices in groups:
  lo,hi=indices[0],indices[-1];gaps=[]
  if lo:gaps.append(E[lo]-E[lo-1])
  if hi+1<len(E):gaps.append(E[hi+1]-E[hi])
  gap=min(gaps)if gaps else F(10**10);delta=4*eta+4*epsilon/gap;w=sum(weights[i]for i in indices);p=[max(F(0),w-delta),min(F(1),w+delta)]
  energy=[E[lo]-epsilon,E[hi]+epsilon]
  rows.append({'eigenvalue_indices':indices,'energy_source_interval':list(map(str,energy)),'conditional_energy_GeV_interval':[str(x*unit)for x in energy],'dyadic_frame_probability':str(w),'conditional_probability_interval':list(map(str,p)),'conditional_probability_diagnostic':float(w),'projector_probability_error_upper':str(delta),'outside_cluster_dyadic_gap':str(gap)})
 complete={'uniform_input_dimension':d['spin_dimension'],'preparation_success_probability':str(F(1,d['spin_dimension'])),'preparation_failure_probability':str(1-F(1,d['spin_dimension'])),'prepared_parity_complement_dimension':d['spin_dimension']-H.rows,'prepared_parity_complement_probability':'0','all_other_source_charge_spin_and_Fock_outcomes':'retained with zero probability on the admitted input','actual_preparation_even_basis_integer_coefficients':[str(x)for x in v],'actual_preparation_squared_integer_norm':str(N),'nominal_strong_ground_vector_deviation_diagnostic':deviation,'original_spectral_probability_rule':'exact prepared spectral weight / uniform_input_dimension; every complementary flag retained','spectral_and_nonselective_preparation_mean_work':'0 exactly on uniform input; spectral readout commutes with full H','opposite_charge_idle_reference':-charge,'reference_energy':'0 declared'}
 mean=float(np.vdot(exact_state,d['H']@exact_state).real);mean0=float(np.vdot(exact_state,d['H0']@exact_state).real);contact=float(np.vdot(exact_state,d['V']@exact_state).real);assert abs(contact-(2*charge-1))<1e-11
 return {'source_member':member,'electric_charge':charge,'complete_charge_sector_dimension':d['charge_dimension'],'fixed_twice_spin_projection_one_dimension':d['spin_dimension'],'exact_even_parity_dimension':H.rows,'actual_strong_parameters':{k:str(x)for k,x in d['provider'].parameters().items()},'contact_coefficient':str(G),'outgoing_kinetic_coefficient':str(d['q']),'complete_prepared_spectrum_certificate':cert,'complete_spectral_clusters':rows,'complete_records':complete,'source_expectation_GeV':mean*float(unit),'strong_expectation_GeV':mean0*float(unit),'exact_reference_contact_expectation':str(2*charge-1),'reference_expectation_shift_GeV':str((2*charge-1)*G*unit),'a_single_interacting_baryon_mass_inferred_from_the_expectation':False}

def calculate():
 members=[]
 for member in('mass_only','axial'):
  n=measure(member,0);p=measure(member,1);pl=max(p['complete_spectral_clusters'],key=lambda r:r['conditional_probability_diagnostic']);intervalp=list(map(F,pl['energy_source_interval']));unit=F(mass_only.INPUTS['electron_unit_GeV']);budgets=[]
  for line in n['complete_spectral_clusters']:
   if line['conditional_probability_diagnostic']<.01:continue
   E=list(map(F,line['energy_source_interval']));budget=[(E[0]-intervalp[1]-1)*unit,(E[1]-intervalp[0]-1)*unit]
   budgets.append({'neutron_line_indices':line['eigenvalue_indices'],'neutron_conditional_weight_interval':line['conditional_probability_interval'],'proton_receiver_line_indices':pl['eigenvalue_indices'],'finite_line_electron_energy_budget_GeV_interval':list(map(str,budget)),'positive_energy_budget':budget[0]>0})
  assert len(budgets)>=3 and all(r['positive_energy_budget']for r in budgets)
  members.append({'member':member,'neutron':n,'proton':p,'strong_calibration_targets_retained':True,'proton_dominant_line_probability_interval':pl['conditional_probability_interval'],'complete_significant_neutron_line_energy_budgets':budgets,'omitted_small_weight_lines':'all retained in the complete_spectral_clusters arrays; this display threshold changes no source or record'})
 return {'schema':'oph.complete_native_contact_baryon_spectral_lines.v1','verdict':'PASS_COMPLETE_DRESSED_LINES_AND_REFERENCE_BOUND_ENERGY_BUDGETS','source':'complete native two-site strong baryon action plus the unchanged full semileptonic current and filled neutral reference; mass-only and axial-calibrated members labelled separately','members':members,'same_source_charged_rest_gap':{'electron_source_energy':'1','electron_GeV':mass_only.INPUTS['electron_unit_GeV'],'proof':'The quark-vacuum sector is exactly reducing and the cross-current vanishes. Adding the charged zero-momentum mode to the retained free reference costs exactly1.','neutral_hole_rest_energy':'0 in this filled-reference member'},'energy_budget_scope':'Differences of discrete finite spectral lines measured with common source units and retained charge references, subtracting the exact charged rest gap. No asymptotic proton-electron factorization, continuum decay endpoint or relativistic physical pole is established.','calibration_ancestry':'mass-only fits p,n,charged pion; axial member additionally fits |gA|. Weak G and outgoing box are declared. Contact spectral shifts are computed without refitting these source parameters.','complete_records_and_orthogonal_Fock_complements_retained':True,'physical_interacting_nucleon_poles_or_neutron_lifetime_claimed':False}
