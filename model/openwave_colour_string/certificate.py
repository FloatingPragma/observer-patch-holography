"""Complete chain records, exact colour-string laws and controlled spectra."""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from .source import chain,native,native_transfer
from .fixed_witness import spectral_enclosure

WITNESSES=json.loads((Path(__file__).parent/'spectral_witnesses.json').read_text())

@lru_cache(None)
def execute(L):
    d=chain(L);cert=spectral_enclosure(d['exact'],WITNESSES[str(L)]);e,v=eigh(np.array(d['exact'].evalf(18),float));eps=F(cert['all_eigenvalue_error_upper']);rank=1
    while rank<len(e) and e[rank]-e[rank-1]<=max(1e-4,8*float(eps)):rank+=1
    state=np.zeros((len(d['labels']),rank));state[d['ids'],:]=v[:,:rank]
    prob=np.sum(abs(state)**2,axis=1)/rank;angular=d['angular'];S=angular(tuple(range(L)),True);I=angular(tuple(range(L)),False)
    spin=float(np.vdot(state,S@state)/rank);iso=float(np.vdot(state,I@state)/rank);assert abs(spin-.75)<1e-8 and abs(iso-.75)<1e-8
    weights=np.bincount(d['spans'],weights=prob,minlength=L);p=d['parameters'];unit=F(native.INPUTS['electron_unit_GeV'])
    energy=[F(cert['all_eigenvalue_intervals'][0][0]),F(cert['all_eigenvalue_intervals'][rank-1][1])];slope=F(4,3)*p['kappa_colour'];tails=[]
    for R in range(1,L):
        bound=min(F(1),(energy[1]-d['floor'])/(slope*R));actual=float(weights[R:].sum());assert actual<=float(bound)+1e-9
        tails.append({'minimum_span_cells':R,'actual_probability_diagnostic':actual,'exact_uniform_Markov_upper':str(bound)})
    # Dyadic eigenvector line has no degenerate spin partners inside the selected projection.
    gap=F(cert['all_eigenvalue_intervals'][rank][0])-energy[1];assert gap>0
    eta=F(cert['exact_Gram_error_upper']);ground_record_error=min(F(1),8*eta+8*eps/gap)
    Ip=0*d['H'];Jp=0*d['H']
    for x in range(L):
        for ss in range(2):
            term=d['E'](d['index'][x,ss,0],d['index'][x,ss,1]);Ip+=term;Jp+=(1 if ss==0 else -1)*term
    gA=float(np.vdot(Ip@state,Jp@state)/rank);gV=float(np.linalg.norm(Ip@state)/np.sqrt(rank))
    raw=np.asarray(d['B']@state);N=len(raw);diagonal=np.sum(abs(raw)**2,axis=1);success=diagonal/N;failure=(1-diagonal)/N
    assert abs(success.sum()+failure.sum()-1)<1e-12
    return {'sites':L,'preparation':'Uniform mixture on the complete isolated lowest spectral cluster, including any unresolved near-degeneracy.','full_native_N3_dimension':N,'all_colour_singlet_dimension':len(state),
      'fixed_twice_spin_projection_one_and_charge_zero_dimension':len(d['ids']),
      'complete_fixed_projection_spectrum_certificate':cert,
      'all_colour_singlet_eigenvalues_diagnostic':eigh(d['H'].toarray(),eigvals_only=True).tolist(),
      'lowest_cluster_rank':rank,'prepared_lowest_band_GeV_interval':[str(x*unit)for x in energy],
      'selected_ground_GeV_interval':[str(F(x)*unit)for x in cert['all_eigenvalue_intervals'][0]],
      'prepared_lowest_band_mean_GeV_diagnostic':float(e[:rank].mean())*float(unit),
      'band_grouping_resolution_source':'1/10000',
      'selected_ground_GeV_diagnostic':e[0]*float(unit),
      'selected_spin_squared_diagnostic':spin,'selected_isospin_squared_diagnostic':iso,
      'native_gV_diagnostic':gV,'native_gA_diagnostic':gA,
      'all_conditional_span_probabilities_diagnostic':weights.tolist(),
      'complete_prepared_band_record_error_upper':str(ground_record_error),
      'uniform_positive_source_floor':str(d['floor']),'span_tail_bounds':tails,
      'native_isometry_generator_errors':d['generator_errors'],'global_colour_annihilation_error':d['colour_error'],
      'complete_native_records':{'input_domain':'Uniform full native N3 carrier before colour preparation. Successful singlet records admit the zero-flux Gauss realization; coloured failures belong to the retained kinematic extension.',
       'occupation_probe':'Compress each native occupation effect to the complete singlet, Gauss-dress it and take its positive square root as Kraus operator. This physical conditional probe has the listed probabilities; only preparation work is asserted zero.',
       'uniform_input_dimension':N,'preparation_success_exact':str(F(rank,N)),
       'preparation_failure_exact':str(1-F(rank,N)),
       'all_native_occupation_success_probabilities_diagnostic':success.tolist(),
       'all_native_occupation_failure_probabilities_diagnostic':failure.tolist(),
       'all_colour_spin_charge_and_Fock_complements':'retained; zero success and unchanged failure outside the selected preparation',
       'idle_reference_charge':'0 for the selected neutral three-quark preparation',
       'nonselective_preparation_mean_work':'0 exactly on the uniform admitted input'},
      'parameters_fixed_without_graph_refitting':{k:str(v)for k,v in p.items()},
      'adjusted_local_colour_coefficient':str(d['adjusted_local_colour'])}

def calculate():
    rows=[execute(L)for L in(2,3,4,5)];p=native.parameters();unit=F(native.INPUTS['electron_unit_GeV']);slope=F(4,3)*p['kappa_colour']
    target=F(native.INPUTS['mass_calibration_targets_GeV']['neutron']);assert F(rows[0]['selected_ground_GeV_interval'][0])<=target<=F(rows[0]['selected_ground_GeV_interval'][1])
    return {'schema':'oph.native_Gauss_colour_string_chains.v1','verdict':'PASS_EXACT_COLOUR_STRING_AND_FIXED_SOURCE_CHAIN_RECORDS',
      'exact_native_transfer_execution':native_transfer(),'unit_elimination_control':unit_elimination(),'uniform_relative_confinement':uniform_confinement(),'source_parameters':{k:str(v)for k,v in p.items()},'same_two_cell_source_transfer':{
       'local_colour_subtraction':str(p['kappa_colour']/2),'adjusted_U':str(p['U']-p['kappa_colour']/2),
       'exact_identity':'On every two-cell global singlet, C0^2=C1^2 and (U-kappa/2)(C0^2+C1^2)+kappa*C0^2=U(C0^2+C1^2). The global-colour term vanishes. All colour-preserving contacts preserve this domain.',
       'domain':'Every two-cell global colour singlet at any particle number, including the three-quark, six-quark and colour-neutral contact sectors. No equality of coloured complements or growing-chain spectra is asserted.'},
      'Gauss_law_source':{'open_boundary_left_flux':'0','link_electric_generator':'E_x=sum_{y<=x} C_y','right_boundary_flux':'0 on the global singlet',
       'electric_action':'kappa sum_x E_x^2, defined by the open-chain axial-gauge Gauss elimination',
       'exact_native_3Q_prefix_identity':'C_X^2=(2/3)n_X(3-n_X)',
       'exact_string_energy':'(4*kappa/3)*(max_position-min_position)',
       'slope_source_units_per_cell':str(slope),'slope_GeV_per_cell':str(slope*unit),
       'static_quark_diquark_records':[{'separation_cells':r,'excess_electric_energy_source':str(slope*r),'excess_electric_energy_GeV':str(slope*r*unit)}for r in range(0,17)],
       'gauge_scope':'Finite Gauss-eliminated one-dimensional matter/electric construction with fixed N3. No independent magnetic link dynamics or 3+1-dimensional Yang-Mills source is claimed.'},
      'fixed_source_chains':rows,
      'calibration_ancestry':'Exactly the selected five-input two-cell source is matched at L2 by the colour identity. Coefficients stay fixed as chain size changes; no additional masses or lengths are fitted.',
      'scope':'Finite colour-singlet mass bands and native position-span records with an exact growing string cost. No physical length calibration, infinite-volume particle pole, pair creation or string-breaking law is established.'}


def uniform_confinement():
    d=chain(2);p=d['parameters'];trial=3*p['mbar']+p['delta']/2+p['A']+F(3,4)*(p['Js']+p['Jw']+p['eta']+p['spin_gradient'])+3*p['t']
    C=(trial-d['floor'])/(F(4,3)*p['kappa_colour']);assert C>0
    return {'all_chain_lengths_at_least_two':True,'exact_endpoint_colocated_spin_half_trial_energy':str(trial),
      'exact_uniform_ground_mean_span_upper':str(C),
      'all_ground_state_span_tail_rule':'P(span>=R)<=min(1,C/R), C the exact reported bound, independent of chain length.',
      'domain':'Fixed N3 total colour singlet and the admitted charge/spin projection; relative-coordinate confinement only, without a limit or fixed absolute centre.'}

def unit_elimination():
    import json
    from pathlib import Path
    inputs=json.loads((Path(__file__).parent/'inputs.json').read_text());parent=json.loads((Path(__file__).parents[1]/'openwave_spin_gradient_nuclei/results.json').read_text())
    q=list(map(F,parent['selected_nuclear_record']['quadrupole_graph_squared_interval']));p=native.parameters();unit=F(native.INPUTS['electron_unit_GeV']);Q=F(inputs['quadrupole_benchmark_fm_squared']);sigma=F(inputs['string_tension_benchmark_GeV_per_fm'])
    slope=F(4,3)*p['kappa_colour']*unit;invariant=[slope*slope*x for x in q];target=sigma*sigma*Q
    maxsigma2=(F(8,3)*p['U']*unit)**2*q[1]/Q
    minkappa2=(3*sigma/(4*unit))**2*Q/q[1]
    assert invariant[1]<target and maxsigma2<sigma*sigma and minkappa2>(2*p['U'])**2
    return {'comparison_inputs':inputs,'unit_free_identity':'sigma^2 Q_physical = (4*kappa*electron_GeV/3)^2 Q_graph',
      'fixed_kappa_unit_free_product_GeV_squared_interval':list(map(str,invariant)),
      'benchmark_product_GeV_squared':str(target),
      'quadrupole_matching_cell_length_squared_fm_interval':[str(Q/q[1]),str(Q/q[0])],
      'fixed_kappa_string_tension_at_that_length_GeV_per_fm_interval_diagnostic':[(float(slope*slope*x/Q))**.5 for x in q],
      'positivity_preserving_kappa_upper':str(2*p['U']),
      'any_positivity_preserving_kappa_string_tension_squared_upper':str(maxsigma2),
      'any_positivity_preserving_kappa_string_tension_upper_diagnostic':float(maxsigma2)**.5,
      'required_kappa_squared_lower':str(minkappa2),
      'simultaneous_benchmark_match_possible_with_nonnegative_adjusted_U':False,
      'scope':'Conditional comparison identifying the finite graph quadrupole with the supplied physical benchmark and using one shared length. A source-specific obstruction for this transfer family; no general OPH or confinement no-go theorem.'}
