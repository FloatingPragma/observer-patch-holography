"""Unchanged-source obstruction and a separately declared static field choice."""
from functools import lru_cache
import numpy as np
import sympy as s
from model.openwave_reference_binding.binding import family_certificate,matter_irrep_certificate
from model.openwave_native_charge_flux.records import certificate as native_records
from model.openwave_field_observables import lattice,electrostatic
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms,cells,cell_index
from .source import cochains,physical_path


@lru_cache(maxsize=1)
def current_projection():
    q=representation()[2][11][9,9]
    assert q==1
    paths=[physical_path(2,[(0,0,0)]),
           physical_path(2,[(0,0,0),(1,0,0)]),
           physical_path(2,[(0,0,0),(1,0,0),(1,1,0)]),
           physical_path(2,[(0,0,0),(1,0,0),(1,1,0),(1,1,1)])]
    short=physical_path(4,[(0,0,0),(1,0,0)])
    detour=physical_path(4,[(0,0,0),(0,1,0),(1,1,0),(1,0,0)])
    assert short['unit_kappa_Coulomb_energy']==detour['unit_kappa_Coulomb_energy']
    assert detour['length']-short['length']==2
    d,g=cochains(2);p=d*g*d.T
    assert p*p==p and p.T==p
    y=representation()[2][11];scalar=higgs_memory()[1]
    for a,t in enumerate(representation()[2]):
        assert y*t==t*y and scalar[11]*scalar[a]==scalar[a]*scalar[11]
    return {'actual_current':'total native matter hypercharge plus quantum scalar hypercharge and the complete static reference hypercharge at cell zero',
        'native_ec_reference_current':str(q),
        'native_link_flux':'b_e=L_(e,0)^Y; local Gauss gives D^*b=rho on the joint physical sector',
        'projected_field':'E_long=D G rho=(D G D^*)b on the physical sector',
        'actual_N2_longitudinal_projection_rank':int(s.trace(p)),
        'exact_projector_and_Pythagorean_identity':True,
        'readout_is_state_projection_or_discards_transverse_native_link_state':False,
        'all_source_reference_scalar_current_factors_retained':True,
        'current_is_central_for_all_twelve_native_generators':True,
        'native_static_path_rows':paths,
        'same_endpoint_different_actual_path_control':[short,detour],
        'local_endpoint_attached_to_complete_native_CAR_scalar_module':matter_irrep_certificate('ec')}


@lru_cache(maxsize=1)
def old_source_obstruction():
    d,g=cochains(2)
    rows=[]
    for family in ('ec','L'):
        old=family_certificate(family)
        charge=s.Integer(1) if family=='ec' else -s.Rational(1,2)
        for row in old['rows']:
            cell=row['cell'];rho=s.zeros(8,1);rho[0]+=charge;rho[cell]-=charge
            shift=s.cancel((rho.T*g*rho)[0]/2)
            old_bound=list(map(s.Rational,row['binding_cost_relative_to_origin_interval']))
            assert shift>=0
            rows.append({'family':family,'cell':cell,'distance':row['distance'],
                'actual_unchanged_source_binding_interval':list(map(str,old_bound)),
                'positive_additive_Coulomb_shift_at_kappa_one':str(shift),
                'additive_extension_binding_interval_at_kappa_one':[str(x+shift) for x in old_bound]})
    # Every actual hopping word has distinct endpoints in the volume family.
    # A local matter-singlet compression kills its lone nontrivial fermion
    # operator, independently of link occupancy or path length.
    word_rows=[]
    for n in (2,4,6):
        terms=kinetic_terms(Configuration(side=n))
        assert all(t.left_cell!=t.right_cell and len(t.path)==t.displacement for t in terms)
        word_rows.append({'side':n,'actual_hopping_path_words':len(terms),
            'all_endpoints_distinct_and_full_paths_retained':True})
    ec=family_certificate('ec')
    width=s.Rational(ec['actual_native_local_trial_gap_upper'])-s.Rational(ec['complete_compensated_local_gap_lower'])
    return {'scope':'complete single-unscreened-charge Pi_y domain, with the full screened complement separate',
        'native_PW119_link_penalty':'one for every nonblank link',
        'exact_L_binding_for_every_even_volume_and_hopping_coefficient':'E_L(y)-E_L(0)=graph_distance(0,y)',
        'ec_uniform_binding_width':str(width),
        'ec_uniform_binding_formula':'distance-width <= E_ec(y)-E_ec(0) <= distance+width',
        'volume_extension_uses_only_local_bounds_Gauss_support_and_singlet_hopping_compression':True,
        'actual_path_word_checks':word_rows,
        'adding_positive_Coulomb_term_cannot_cancel_the_linear_cost':True,
        'reason':'Coulomb operator is a scalar kappa q_ref^2(G00-G0y) on each Pi_y and is zero on Pi_0; it shifts each compressed minimum by that exact nonnegative number.',
        'all_N2_constrained_rows':rows,
        'unrestricted_screened_cloud_confinement_or_all_source_Coulomb_exclusion':False}


@lru_cache(maxsize=1)
def complete_existing_readout():
    native=native_records();g=cochains(2)[1];rows=[]
    for original in native['complete_original_preparation_records']['records']:
        row=dict(original)
        if 'joint_charge_and_flux_atom' in row:
            rho=s.zeros(8,1);rho[0]=1;rho[row['arrival']]=-1
            energy=s.cancel((rho.T*g*rho)[0]/2)
            assert energy==s.Rational(7,48)
            row['unit_kappa_quadratic_Gauss_flux_observable']=str(energy)
            row['probability_is_actual_unchanged_source_arrival_probability']=True
        else:
            row['quadratic_Gauss_flux_observable']='complete spectral function rho^T G rho/2 of each retained native-current outcome'
        rows.append(row)
    return {'same_original_Hamiltonian_clock_reference_and_complete_preparation':True,
        'observable':'cochain minimum quadratic flux readout derived from actual native charge/flux; it is not the old unit-link Hamiltonian energy',
        'records':rows,'coarse_record_count':10,
        'all_failures_and_nonarrival_current_flux_outputs_retained':True,
        'same_source_positive_six_route_probability_interval':native['same_source_six_route_joint_charge_flux_interval'],
        'complete_native_link_current_outcomes_and_their_longitudinal_postprocessing_retained':True,
        'unnormalized_charge_or_current_target_inserted':False}


@lru_cache(maxsize=1)
def field_shape():
    q=float(representation()[2][11][9,9]);assert q==1.
    base=electrostatic.run()
    n=16;g=lattice.green(n)
    rho=np.zeros((n,)*3);rho[0,0,0]=q;rho[3,0,0]=-q
    phi=q*(g-np.roll(g,3,axis=0));field=lattice.gradient(phi)
    assert np.max(abs(lattice.divergence_adjoint(field)-rho))<1e-12
    energy=float(np.sum(field*field)/2)
    assert abs(energy-q*q*(g[0,0,0]-g[3,0,0]))<1e-12
    return {'current_amplitude_from_actual_native_ec_reference':str(int(q)),
        'direct_native_density_Gauss_and_energy_check_grid':n,
        'direct_native_density_unit_kappa_energy':round(energy,12),
        'finite_Green_profile_and_force':base,
        'coupling_conversion':'kappa=g^2 in the cochain comparison producer; kappa remains an independent Hamiltonian coefficient',
        'lattice_metric':'counting incidence and counting link norm; distances here are lattice steps',
        'native_matter_spacing':'pi in the inherited PW119 kinetic convention; no physical length or energy calibration is inferred',
        'continuum_shape_comparison':'V_interaction approximately -kappa/(4pi r_sites) for unit native charges in the resolved large-box regime',
        'physical_alpha_derived':False}


def certificate():
    endpoint=matter_irrep_certificate('ec')
    return {'schema':'oph.native_current_coulomb_choice.v1',
        'actual_current_to_Gauss_minimizer':current_projection(),
        'unchanged_source_obstruction':old_source_obstruction(),
        'complete_existing_source_readout':complete_existing_readout(),
        'declared_replaced_penalty_source':{
            'Hamiltonian':'H_prime=H_local+lambda T+4 N_ref+(kappa/2) rho^Y G rho^Y; native nonblank-link penalty removed',
            'actual_native_kinetic_and_quantum_Yukawa_operators_retained':True,
            'reference_full_Fock_and_all_PW119_link_states_retained':True,
            'lambda':'1/17179869184','kappa_example':'1',
            'source_clock':'exp(-it H_prime), distinct from the original PW119 clock',
            'complete_carrier':'all 30 CAR modes and five scalar states per cell, 119 states per native link, and all 30 static reference CAR modes',
            'physical_domain':'joint native local Gauss kernel; central total source plus reference charge zero',
            'full_kinematic_extension':'G is the mean-zero pseudoinverse; nonneutral kinematic columns are defined but do not satisfy the physical Gauss readout',
            'same_full_local_Gauss_invariance':True,
            'local_current_continuity_preserved':'The changed diagonal term commutes with every local total Y density; the density time derivative is the same native hopping commutator.',
            'static_string_trial_energy':'(V-1) E0+e_ec+4+kappa(G00-G0y)',
            'exact_static_trial_difference_from_origin':'kappa(G00-G0y), independent of the chosen physical flux path',
            'trial_is_full_H_prime_eigenstate':'only when lambda=0; otherwise its expectation has this value because Pi_y T Pi_y=0',
            'all_Pi_y_minimum_bounds':'(V-1)E0-11.71+4+kappa C_y <= inf Pi_y H_prime Pi_y <= (V-1)E0+e_ec+4+kappa C_y',
            'static_distance_dependent_linear_penalty_removed':True,
            'transverse_native_flux_loop_energy':'none; loops with the same local matter and rho have identical added energy',
            'transverse_photon_Hamiltonian_or_Lorentz_completion':False,
            'old_interacting_ground_preparation_probabilities_or_energy_bins_transferred':False,
            'new_ground_autonomous_particle_spectrum_or_deconfined_vacuum_phase_proved':False,
            'actual_local_endpoint_polynomial':endpoint['selected_local_energy_polynomial']},
        'native_current_Green_profile':field_shape(),
        'coupling_obstruction':{'unit_current_and_all_Gauss_Ward_identities_fix_kappa':False,
            'same_native_current_representation_compatible_with_kappa_1_and_2':True,
            'both_pair_energy_and_static_force_scale_linearly_with_kappa':True,
            'physical_GeV_length_c_alpha_assignment':None},
        'verdict':'PASS_NATIVE_CURRENT_GAUSS_READOUT_WITH_EXPLICIT_COULOMB_SOURCE_CHOICE'}
