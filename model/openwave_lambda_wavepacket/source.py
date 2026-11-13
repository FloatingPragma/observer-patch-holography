"""Normalizable strange-particle trial, full pole weight and native CAR records."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as sy
from model.openwave_cluster_bandwidth import lambda_source as lam
from model.openwave_strange_shell_current.source import lambda_pole,trial
from model.openwave_strange_shell_current.native import contraction
from model.openwave_scattering_current_support.native import channel
from model.openwave_baryon_dispersion.source import repair_matrix
from model.openwave_baryon_poles.source import trial_moments,real_fraction
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_normalizable_particle_records.source import moment,cosine_mean
from model.openwave_capacity_stress.source import enclosure
KAPPA=F(1,512)
POSITION_RADIUS=8192


@lru_cache(None)
def packet():
    rest=lambda_pole();blocks,G,data=lam.compact()
    A=sum(blocks.values(),sy.zeros(7));v=trial(rest['actual_trial_coordinates'])
    n,mean,variance=trial_moments(A,G,v)
    zeta=repair_matrix()[2]
    J=sy.Rational(zeta)*G*A.inv()*G
    assert J==J.T and A*G.inv()*J==sy.Rational(zeta)*G
    e0lo,e0hi=map(F,rest['all_seven_compact_source_intervals'][0])
    norm=zeta/e0lo;gap=F(rest['certificate']['complete_ground_gap_lower_source'])
    perturb=6*data['p']['t']*KAPPA+norm*KAPPA*KAPPA
    assert 2*perturb<gap
    distance=perturb/(gap-perturb)
    sine=sqrt_upper(1-F(rest['certificate']['prepared_trial_pole_weight_lower']))+distance
    band_weight=1-sine*sine;assert 0<band_weight<1
    C=F(rest['complete_relative_Q_floor_source'])
    q=(6*data['p']['t'])**2/(C-e0hi-perturb)**2
    assert 0<q<1 and e0hi+perturb<C
    # Even |f|² removes the odd imaginary Fourier block. The exact remaining
    # source expectation uses the same common inverse-rest propagation.
    Api=blocks[0]-blocks[1]-blocks[-1]+4*J
    coefficient=real_fraction((v.T*(Api-A)*v)[0])/(2*n)
    assert coefficient>0
    clo,chi=cosine_mean(KAPPA)
    energy=((mean+coefficient*(1-chi))*data['unit'],(mean+coefficient*(1-clo))*data['unit'])
    full_band=(F(rest['full_ground_GeV_interval'][0])-perturb*data['unit'],
               F(rest['full_ground_GeV_interval'][1])+perturb*data['unit'])
    c=channel();M=c['projected_Gram'];T=c['total_Gram'];assert c['source_Gram']==G
    cab=contraction()['Cabibbo_squared']
    projected=cab*real_fraction((v.T*M*v)[0])/n/1296
    total=cab*real_fraction((v.T*T*v)[0])/n/1296
    assert 0<projected<total<1
    root=sqrt_upper(F(54409));a=(F(277)-root)/200;b=(F(277)+root)/200
    pole_projected=(cab*a*(1-q)/1296,cab*b/1296)
    pole_other=(cab*(6-b)*(1-q)/1296,cab*(7-a)/1296)
    pole_failure=(1-7*cab/1296-cab*q,1-6*cab*(1-q)/1296)
    position_var=F(5)/(2*KAPPA*KAPPA)
    tail=position_var/POSITION_RADIUS**2
    conditioned_tail=tail+sqrt_upper(1-band_weight)
    assert 0<conditioned_tail<1
    return {'rest':rest,'G':G,'J':J,'trial':v,'trial_norm':n,'trial_mean':mean,
        'trial_variance':variance,'coefficient':coefficient,'kinetic_norm':norm,'gap':gap,
        'perturbation':perturb,'projector_distance':distance,'band_weight':band_weight,
        'full_Q_upper':q,'energy':energy,'full_band_energy':full_band,'unit':data['unit'],
        'trial_channel_probability':projected,'trial_total_success':total,
        'pole_channel_probability':pole_projected,'pole_other_probability':pole_other,
        'pole_failure_probability':pole_failure,'Cabibbo_squared':cab,
        'position_variance':position_var,'position_tail':tail,'conditioned_position_tail':conditioned_tail,
        'velocity':6*data['p']['t']+2*norm,
        'position_work_bound':2*(6*data['p']['t']+4*norm)*data['unit']}


def run():
    d=packet();k=sy.Symbol('k',real=True);kap=sy.Rational(KAPPA)
    density=sy.Rational(15,16)/kap*(1-k*k/kap**2)**2
    assert sy.integrate(density,(k,-kap,kap))==1
    assert sy.integrate(k*k*density,(k,-kap,kap))==sy.Rational(moment(2,KAPPA))
    # All literal words remove s and recreate u+d+anti-u at the same cell:
    # the leftmost occupied cell is unchanged, so the finite compact current
    # and channel act as I_translation tensor the exact CAR internal maps.
    from model.openwave_strange_shell_current.native import vertex_words
    rows=[]
    for coefficient,ops in vertex_words():
        removed=[i//24 for create,i in ops if not create]
        created=[(i%48)//24 for create,i in ops if create]
        assert len(removed)==1 and all(x==removed[0] for x in created)
        rows.append(removed[0])
    assert len(rows)==72
    trial_outcomes={'native_proton_pion_fields':d['trial_channel_probability'],
                    'other_compact_N5':d['trial_total_success']-d['trial_channel_probability'],
                    'noncompact_N5':F(0),'probe_failure':1-d['trial_total_success']}
    assert sum(trial_outcomes.values())==1
    original={'preparation_failure':F(6,7),**{name:value/7 for name,value in trial_outcomes.items()}}
    assert sum(original.values())==1
    rest=d['rest']
    return {'schema':'oph.normalizable_native_lambda_packet.v1',
        'verdict':'PASS_NORMALIZABLE_LAMBDA_TRIAL_NATIVE_CHANNEL_RECORDS_AND_UNIFORM_FULL_POLE_WEIGHT',
        'normalizable_state':{'momentum_half_width_radians_exact':str(KAPPA),
            'momentum_amplitude':'sqrt(15/(16*kappa))*(1-(k/kappa)^2) on |k|<=kappa, zero outside; Lebesgue dk normalization.',
            'exact_L2_norm_squared':'1','momentum_variance_exact':str(moment(2,KAPPA)),
            'position_variance_cells_squared_exact':str(d['position_variance']),
            'source_state':'f(k) times the inherited exact seven-coordinate native Lambda trial, with exact Gram normalization; infinite-chain L2, not a finite ring or delta momentum.',
            'native_trial_coordinates':rest['actual_trial_coordinates'],
            'trial_Gram_norm_squared':str(d['trial_norm']),
            'translation_carrier':'Leftmost occupied cell in Z; every relative configuration and all Fock complements retained.',
            'actual_normalizable_Lambda_trial_constructed':True},
        'calibrated_energy':{'source_rest_pole_GeV_interval':rest['full_ground_GeV_interval'],
            'source_rest_pole_GeV_display':rest['full_ground_GeV_display'],
            'trial_packet_mean_GeV_interval':list(map(str,d['energy'])),
            'trial_packet_mean_GeV_display':enclosure(d['energy']),
            'mean_energy_is_rest_mass':False,
            'all_ground_band_energies_on_packet_support_GeV_interval':list(map(str,d['full_band_energy'])),
            'all_ground_band_energies_on_packet_support_GeV_display':enclosure(d['full_band_energy']),
            'inverse_rest_rule_shared_with_all_compact_clusters':True,
            'inverse_rest_operator_norm_upper_source':str(d['kinetic_norm']),
            'new_mass_or_length_fit':False},
        'complete_source_pole':{'k_zero_ground_gap_lower_source':str(d['gap']),
            'uniform_H_k_minus_H_zero_norm_upper':str(d['perturbation']),
            'uniform_ground_gap_lower_source':str(d['gap']-2*d['perturbation']),
            'uniform_projector_distance_from_rest_upper':str(d['projector_distance']),
            'unfiltered_trial_ground_band_probability_lower':str(d['band_weight']),
            'unfiltered_trial_ground_band_probability_display':enclosure((d['band_weight'],F(1))),
            'unfiltered_trial_other_spectrum_probability_upper':str(1-d['band_weight']),
            'full_relative_complement_probability_upper_on_ground_band':str(d['full_Q_upper']),
            'infinite_relative_complement_removed':False},
        'actual_trial_current_records':{'literal_native_vertex_words_checked':72,
            'translation_factor_exact':'I_translation tensor W_compact; the weak vertex preserves every occupied cell and its leftmost cell.',
            'native_component_probe_success_operator':'W/36, with sqrt(I-W*W/1296) on the complete source and Fock complement.',
            'conditional_on_prepared_trial_exact_probabilities':{k:str(v) for k,v in trial_outcomes.items()},
            'original_P7_over_7_input_exact_probabilities':{k:str(v) for k,v in original.items()},
            'native_proton_pion_event_probability_display':enclosure((d['trial_channel_probability'],d['trial_channel_probability'])),
            'all_probabilities_sum_exactly_to_one':True,
            'exact_ground_band_filter_required_for_these_trial_records':False,
            'current_channel_conditioning_preserves_spatial_packet_exactly':True,
            'source_trial_is_exact_Lambda_eigenstate':False},
        'position_records':{'radius_cells':POSITION_RADIUS,
            'unfiltered_trial_and_each_nonzero_internal_current_branch_outside_probability_upper':str(d['position_tail']),
            'inside_probability_lower':str(1-d['position_tail']),
            'complete_terminal_outcomes':['left tail','inside interval','right tail'],
            'branch_conditioning_changes_spatial_distribution':False,
            'joint_sharp_position_and_momentum_readout_claimed':False,
            'complete_unfiltered_time_t_second_moment_bound':'(sqrt(5/(2*kappa^2))+abs(t)*velocity)^2',
            'velocity_upper_cells_per_source_time':str(d['velocity']),
            'complete_position_readout_mean_work_absolute_GeV_bound':str(d['position_work_bound'])},
        'optional_ground_band_conditioning':{'bounded_direct_integral_projector':True,
            'success_probability_lower':str(d['band_weight']),'failure_probability_upper':str(1-d['band_weight']),
            'native_proton_pion_probability_interval':list(map(str,d['pole_channel_probability'])),
            'other_compact_N5_probability_interval':list(map(str,d['pole_other_probability'])),
            'noncompact_N5_probability_interval':['0',str(d['Cabibbo_squared']*d['full_Q_upper'])],
            'probe_failure_probability_interval':list(map(str,d['pole_failure_probability'])),
            'original_mixed_input_trial_band_current_channel_probability_lower':str(d['band_weight']*d['pole_channel_probability'][0]/7),
            'conditioned_position_outside_probability_upper':str(d['conditioned_position_tail']),
            'band_filter_nonselective_energy_work_exact':'0 because it commutes with the source Hamiltonian.',
            'finite_local_hardware_compilation_of_exact_band_filter_proved':False,
            'normalizability_preserved':True},
        'complete_protocol_work':{'trial_preparation_success_exact':'1/7','trial_preparation_failure_exact':'6/7',
            'nonselective_trial_preparation_work_exact':'0 on the specified spatial packet tensor P7/7 input.',
            'branch_conditioned_preparation_work_claimed_zero':False,
            'current_and_channel_work':'Tr[H(Phi(rho)-rho)] for the complete instrument; neither probe nor channel readout is asserted energy conserving.',
            'finite_local_compilation_of_entire_protocol_asserted':False,
            'arbitrary_idle_external_reference_retained':True},
        'finite_theorem_shortcut':'Seven-coordinate rest-pole residue and Schur bounds plus bounded momentum perturbation and Parseval; exact CAR Gram contraction factors through translation, so a normalizable trial/current experiment needs no exact ground-state filter.',
        'inherited_inputs':'Electron GeV anchor, common strong and axial calibrations, strange coefficients, Cabibbo angle, supplied physical speed and strong length. Packet width and position radius are declared readout choices.',
        'new_independent_match_with_nature':False,'physical_3plus1_Lambda_identification':False,
        'autonomous_on_shell_decay_or_lifetime_proved':False,
        'scope':'A normalizable native Lambda-dominated trial, exact complete current/channel records, and an optional conditioned full-source ground-band packet in the infinite one-dimensional model. Interacting outgoing weak transfer and physical continuum identification remain open.'}
