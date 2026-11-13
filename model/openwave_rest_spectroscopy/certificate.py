"""Certified finite energy lines and complete Fourier-register records."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from .source import construction,exact_parity_matrix,parity_basis,G,gradient
from model.openwave_six_quark.certificate import spectral_enclosure
@lru_cache(None)
def measure(charge):
 d=construction(charge);H,S=exact_parity_matrix(d);cert=spectral_enclosure(H);D=cert['dyadic_denominator'];O=np.array([[int(x)for x in row]for row in cert['dyadic_eigenframe']],dtype=object);E=[F(int(x),D)for x in cert['dyadic_eigenvalues']];epsilon=F(cert['all_eigenvalue_error_upper']);eta=F(cert['exact_Gram_error_upper'])
 B=parity_basis(d);a=np.asarray(B.T@d['initial']).ravel();v=np.array([int(round(float(x.real)*2**48))for x in a],dtype=object);N=sum(x*x for x in v);assert N>0
 vector=np.array(v,dtype=float)/np.sqrt(float(N));exact_state=B@vector;deviation=float(np.linalg.norm(exact_state-d['initial']));assert deviation<1e-11
 overlap=O.T@v;weights=[F(int(x*x),int(D*D*N))for x in overlap];assert abs(sum(weights)-1)<=eta
 groups=[];start=0
 for i in range(1,len(E)+1):
  if i==len(E)or E[i]-E[i-1]>8*epsilon:groups.append(list(range(start,i)));start=i
 rows=[];unit=F(gradient.INPUTS['electron_unit_GeV'])
 for indices in groups:
  lo,hi=indices[0],indices[-1];gaps=[]
  if lo:gaps.append(E[lo]-E[lo-1])
  if hi+1<len(E):gaps.append(E[hi+1]-E[hi])
  gap=min(gaps)if gaps else F(10**10);delta=4*eta+4*epsilon/gap;w=sum(weights[i]for i in indices);p=[max(F(0),w-delta),min(F(1),w+delta)]
  energy=[E[lo]-epsilon,E[hi]+epsilon]
  rows.append({'eigenvalue_indices':indices,'energy_source_interval':list(map(str,energy)),'conditional_energy_GeV_interval':[str(x*unit)for x in energy],'dyadic_frame_probability':str(w),'conditional_probability_interval':list(map(str,p)),'conditional_probability_diagnostic':float(w),'projector_probability_error_upper':str(delta),'outside_cluster_dyadic_gap':str(gap)})
 complete={'uniform_input_dimension':d['spin_dimension'],'preparation_success_probability':str(F(1,d['spin_dimension'])),'preparation_failure_probability':str(1-F(1,d['spin_dimension'])),'prepared_parity_complement_dimension':d['spin_dimension']-H.rows,'prepared_parity_complement_probability':'0','all_other_source_charge_spin_and_Fock_outcomes':'retained with zero probability on the admitted input','actual_preparation_even_basis_integer_coefficients':[str(x)for x in v],'actual_preparation_squared_integer_norm':str(N),'nominal_strong_ground_vector_deviation_diagnostic':deviation,'original_spectral_probability_rule':'exact prepared spectral weight / uniform_input_dimension; every complementary flag retained','spectral_and_nonselective_preparation_mean_work':'0 exactly on uniform input; spectral readout commutes with full H','opposite_charge_idle_reference':-charge,'reference_energy':'0 declared'}
 mean=float(np.vdot(exact_state,d['H']@exact_state).real);mean0=float(np.vdot(exact_state,d['H0']@exact_state).real);contact=float(np.vdot(exact_state,d['V']@exact_state).real);assert abs(contact-(2*charge-1))<1e-11
 return {'source_member':'selected_spin_gradient','electric_charge':charge,'complete_charge_sector_dimension':d['charge_dimension'],'fixed_twice_spin_projection_one_dimension':d['spin_dimension'],'exact_even_parity_dimension':H.rows,'actual_strong_parameters':{k:str(x)for k,x in d['provider'].parameters().items()},'contact_coefficient':str(G),'outgoing_kinetic_coefficient':str(d['q']),'complete_prepared_spectrum_certificate':cert,'complete_spectral_clusters':rows,'complete_records':complete,'source_expectation_GeV':mean*float(unit),'strong_expectation_GeV':mean0*float(unit),'exact_reference_contact_expectation':str(2*charge-1),'reference_expectation_shift_GeV':str((2*charge-1)*G*unit),'a_single_interacting_baryon_mass_inferred_from_the_expectation':False}


def pointer(member,M):
    cert=member['complete_prepared_spectrum_certificate'];D=cert['dyadic_denominator']
    energies=np.array([int(x)/D for x in cert['dyadic_eigenvalues']])
    O=np.array(cert['dyadic_eigenframe'],float)/D
    v=np.array(member['complete_records']['actual_preparation_even_basis_integer_coefficients'],float)
    v/=np.linalg.norm(v);weights=abs(O.T@v)**2
    # Rational source window is certified against every eigenvalue, not selected peaks.
    offset=F(1792);width=F(256)
    intervals=[[F(x)for x in row]for row in cert['all_eigenvalue_intervals']]
    assert min(x[0]for x in intervals)>offset and max(x[1]for x in intervals)<offset+width
    k=np.arange(M);delta=(energies[:,None]-float(offset))/float(width)-k[None,:]/M
    delta=(delta+.5)%1-.5
    kernel=(np.sinc(M*delta)/np.sinc(delta))**2
    probabilities=weights@kernel
    eta=F(cert['exact_Gram_error_upper']);error=F(cert['all_eigenvalue_error_upper'])
    # pi<22/7; finite Duhamel plus polar-frame perturbation, with no gap hypothesis.
    tmax=F(44,7)*F(M-1)/width;tv=4*eta+4*tmax*error
    assert np.all(probabilities>=0)and np.all(np.isfinite(probabilities))
    assert abs(probabilities.sum()-1)<float(eta)+1e-10
    order=np.argsort(probabilities)[-8:][::-1];unit=F(gradient.INPUTS['electron_unit_GeV'])
    return {'register_dimension':M,'source_energy_offset':str(offset),'source_energy_width':str(width),
      'all_conditional_pointer_probabilities':probabilities.tolist(),
      'all_original_input_pointer_probabilities':(probabilities/member['fixed_twice_spin_projection_one_dimension']).tolist(),
      'preparation_failure_probability':member['complete_records']['preparation_failure_probability'],
      'unrounded_pointer_energy_rule':'offset + width*k/M, k=0,...,M-1; multiply by electron_GeV',
      'energy_bin_width_GeV':str(width*unit/M),
      'largest_pointer_records':[{'k':int(i),'energy_GeV':float((offset+width*F(int(i),M))*unit),'conditional_probability':float(probabilities[i])}for i in order],
      'complete_probability_sum_diagnostic':float(probabilities.sum()),
      'exact_instrument_total_variation_error_upper':str(tv),
      'error_scope':'Exact mathematical Fourier instrument vs recorded dyadic-frame approximation; floating sinc summation additionally checked at 1e-10 absolute total probability. No normalization or omitted bins.',
      'conditional_terminal_state_rule':'K_k psi/sqrt(p_k), K_k=M^-1 sum_n exp(2pi*i*k*n/M) exp(-2pi*i*n*(H-offset)/width); H and preparation are completely specified by the pinned exact source and integer vector.',
      'nonselective_instrument_mean_work_exact':'0: every K_k commutes with H',
      'outside_admitted_charge_spin_parity_or_Fock_space':'identity complementary flag, retained with zero conditional weight',
      'time_interpretation':'Finite controlled source action with a register; no external physical time or universal clock is inferred.'}

def calculate():
    members=[];unit=F(gradient.INPUTS['electron_unit_GeV'])
    for charge in(0,1):
        m=measure(charge);m['operational_pointer_records']=[pointer(m,M)for M in(1024,4096,16384)]
        d=construction(charge,coupling=0);a=d['initial'];variance=float(np.linalg.norm(d['H']@a-np.vdot(a,d['H']@a)*a)**2)
        assert variance<1e-17
        m['zero_contact_ground_variance_diagnostic']=variance;members.append(m)
    assert members[0]['complete_charge_sector_dimension']==1400 and members[1]['complete_charge_sector_dimension']==360
    dominant=[]
    for m in members:
        rows=sorted(m['complete_spectral_clusters'],key=lambda r:r['conditional_probability_diagnostic'],reverse=True)
        dominant.append({'electric_charge':m['electric_charge'],'lines':[r for r in rows if r['conditional_probability_diagnostic']>.01],
        'all_small_weight_lines_retained':True})
    return {'schema':'oph.same_source_operational_baryon_rest_spectroscopy.v1',
      'verdict':'PASS_COMPLETE_CONTACT_REST_LINES_AND_FINITE_FOURIER_INSTRUMENT',
      'source':'Selected five-input spin-gradient strong action, full native semileptonic contact G=1/64, one light charged family and filled massless-neutral reference in a two-cell box.',
      'members':members,'significant_native_lines':dominant,
      'same_source_unit':{'electron_source_rest_energy':'1','electron_GeV':str(unit),'reference':'quark-vacuum charged zero-mode addition to the free filled reference; cross-current vanishes exactly'},
      'calibration_ancestry':'Strong p,n,charged-pion,gA and even nuclear binding are five inputs. The full weak-contact source is not refitted. Spectral splitting and pointer records follow from that fixed source.',
      'all_source_charge_spin_parity_and_Fock_complements_retained':True,
      'scope':'Discrete finite-box rest-energy spectral records in electron-calibrated GeV. Pointer resolution is finite apparatus resolution, not decay width. No infinite-volume pole, physical lifetime or source transfer to the fitted three-family massive Dirac reference is inferred.'}
