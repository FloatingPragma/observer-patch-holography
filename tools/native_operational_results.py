"""Compact source-qualified projections; no diagnostic is promoted to a physical fit."""
NAMES=('lambda_spectral_measure','finite_phase_clock','localized_quark_record',
       'nuclear_em_transfer','fock_field_transport','energy_schur',
       'collective_finite_strength','gravity_noise_normalization','capacity_normalization_transfer')

def project(d):
    def require(name,path,expected):
        v=d[name]
        for key in path:v=v[key]
        if v is not expected:raise ValueError('operational scope changed: '+name+'.'+'.'.join(path))
    for name,path in [
        ('lambda_spectral_measure',('new_match_with_nature',)),
        ('lambda_spectral_measure',('on_shell_native_weak_density_proved',)),
        ('lambda_spectral_measure',('actual_successful_current_energy','conditioned_energy_change_is_nonselective_probe_work')),
        ('finite_phase_clock',('scope','literal_openwave_row3_validated')),
        ('localized_quark_record',('physical_quark_rest_mass_GeV_derived',)),
        ('localized_quark_record',('self_sustained_fractional_texture_winding_proved',)),
        ('nuclear_em_transfer',('new_mass_fit',)),
        ('nuclear_em_transfer',('physical_nuclear_spectrum_or_common_alpha_attachment_proved',)),
        ('fock_field_transport',('interacting_continuum_or_physical_particle_dictionary_proved',)),
        ('fock_field_transport',('inherited_GeV_energies_changed',)),
        ('energy_schur',('nontrivial_continuum_or_Clay_result',)),
        ('energy_schur',('uniform_global_offblock_bound_claimed',)),
        ('collective_finite_strength',('scope','independent_G_prediction')),
        ('gravity_noise_normalization',('scope','independent_Newton_prediction')),
        ('gravity_noise_normalization',('scope','physical_metric_noise_identified')),
        ('capacity_normalization_transfer',('mode_and_observation_scope','new_Nature_match')),
        ('capacity_normalization_transfer',('mode_and_observation_scope','OpenWave_row28_never_fitted_amplitude_obligation_met')),
        ('capacity_normalization_transfer',('mode_and_observation_scope','metric_ratio_is_parent_finite_k_matter_growth_ratio'))]:require(name,path,False)
    for name,path in [
        ('localized_quark_record',('gauge_invariant_localized_order_parameter_record_constructed',)),
        ('fock_field_transport',('full_native_H_current_charge_reference_work_transport',)),
        ('collective_finite_strength',('scope','N_inherited_from_observed_G_calibration')),
        ('gravity_noise_normalization',('scope','observed_G_selects_inherited_N'))]:require(name,path,True)
    l=d['lambda_spectral_measure'];c=d['finite_phase_clock']['autonomous_finite_clock_controller']
    n=d['nuclear_em_transfer'];f=d['fock_field_transport'];g=d['collective_finite_strength']
    q=d['gravity_noise_normalization']['conditional_calibrated_member']
    nuclei=[{k:r[k]for k in ('name','full_dynamic_ground_GeV_interval','target_GeV','target_role','target_in_certified_interval')}for r in n['complete_ground_bounds']]
    if any(r['target_in_certified_interval']for r in nuclei if r['name']in ('helion','alpha')):raise ValueError('held-out nuclear failure omitted')
    return {
      'lambda':{'successful_current_mean_GeV_interval':l['actual_successful_current_energy']['mean_GeV_interval'],
                'energy_is_rest_mass':False,'on_shell_density_proved':False,'scope':l['scope']},
      'clock':{'reducing_carrier':c['executed_reducing_carrier_dimension'],
               'complete_control_cases':len(c['native_charged_mass_controls']),
               'maximum_demodulated_drift_diagnostic':max(r['native_frequency_demodulated_drift_diagnostic']for r in c['native_charged_mass_controls']),
               'native_complement':c['native_complement'],'numerical_interval_certificate':False,'self_starting_sustained_clock':False},
      'string':{'localized_gauge_invariant_record':True,'physical_quark_mass_derived':False,'scope':d['localized_quark_record']['scope']},
      'nuclei':nuclei,
      'fock':{'modes':f['native_transport']['native_CAR_modes'],
              'pair_implementers':f['native_transport']['commuting_even_pair_implementers'],
              'additional_poles':f['additional_common_source_poles_identified'],
              'onto_native_rotor_unitary':f['source_equivalence_obstruction']['onto_unitary_between_these_complete_carriers']},
      'gauge':[{'source':r['source'],'gap_exact':r['full_gap_interval_exact'],
                'width_improvement_diagnostic':r['gap_interval_width_improvement_factor'],
                'complement_floor_exact':r['full_Q_floor_exact']}for r in d['energy_schur']['records']],
      'gravity':{'pairs':g['complete_native_force']['positive_native_source_probe_pairs'],
                 'force_relative_error_upper_exact':g['complete_native_force']['central_force_relative_error_upper_exact'],
                 'mutual_energy_relative_error_upper_exact':g['native_finite_energy_certificate']['mutual_energy_relative_error_upper_exact'],
                 'quench_record_relative_error_upper_exact':g['two_time_coarse_record']['relative_error_upper_exact'],
                 'rms_Q_interval_exact':q['rms_Q_interval_exact'],'coordinate_work_GeV_exact':q['coordinate_work_GeV_exact'],
                 'observed_G_calibration':True,'physical_noise_match':False},
      'capacity':{'present_metric': d['capacity_normalization_transfer']['validated_metric_quadrature']['records'][0],
                  'normalization':d['capacity_normalization_transfer']['exact_normalization_theorem'],
                  'scope':d['capacity_normalization_transfer']['mode_and_observation_scope']},
      'new_independent_nature_validations':0,
    }

def render(v,interval,decimal):
    lines=['### Native clocks, field transport and spectral resolution','',
      'A [finite autonomous controller](../model/openwave_finite_phase_clock/PROOF.md) executes '+str(v['clock']['complete_control_cases'])+
      ' mass/preparation controls on a '+str(v['clock']['reducing_carrier'])+'-dimensional reducing carrier with the native complement retained. '
      'The free coherence frequency is the source mass. Phase references, initial clock and complete terminal measurement work are explicit; '
      'the executed trajectories are numerical controls, not interval certificates or self-starting sustained particle clocks.','',
      'The [actual Lambda current](../model/openwave_lambda_spectral_measure/PROOF.md) has successful-branch mean energy '+
      interval(v['lambda']['successful_current_mean_GeV_interval'])+' GeV. This probe-conditioned five-body energy is not its rest mass, '
      'decay energy or nonselective work. Broad-window positivity does not establish a nonzero shrinking-window density or decay lifetime.','',
      'A [localized native string](../model/openwave_localized_quark_record/PROOF.md) gives a reference-supported fractional-charge coherence '
      'with complete flux and work records. The prepared endpoint is not an autonomous physical quark. '
      'The [full CAR transport](../model/openwave_fock_field_transport/PROOF.md) acts on all '+str(v['fock']['modes'])+
      ' native modes by local even implementers, retaining Hamiltonian, currents, state and work. '
      'Its separate integer-rotor E+A construction requires a different carrier; it supplies no additional particle pole.','',
      'The [dynamical nuclear Gauss extension](../model/openwave_nuclear_em_transfer/PROOF.md) retains every integer cycle flux:', '',
      '| Nucleus | Full source interval / GeV | Target role; admitted |','| --- | ---: | --- |']
    for r in v['nuclei']:
        lines.append('| '+r['name']+' | '+interval(r['full_dynamic_ground_GeV_interval'])+' | '+r['target_role']+'; '+str(r['target_in_certified_interval'])+' |')
    lines+=['','The helion and alpha comparisons fail. A positive field addition with only a common vacuum subtraction cannot lower '
      'the inherited alpha excess. No mass target is refitted.','',
      'The [shifted full-Haar Schur certificate](../model/openwave_energy_schur/PROOF.md) sharpens complete-source gaps:', '',
      '| Source | Full gap / source units | Width improvement, diagnostic |','| --- | ---: | ---: |']
    for r in v['gauge']:lines.append('| '+r['source']+' | '+interval(r['gap_exact'])+' | '+f"{r['width_improvement_diagnostic']:.3f}"+' |')
    lines+=['','All omitted representations enter the remainder. The fixed global energy projector has an extensive off-block norm; '
      'local observable estimates are required for uniform-volume refinement. No physical GeV glueball or Clay result follows.','',
      'The [finite-strength collective theorem](../model/openwave_collective_finite_strength/PROOF.md) controls all '+str(v['gravity']['pairs'])+
      ' positive native source/probe pairs, using two spin-sector Jacobi blocks and a complete complement bound. '
      'At the inherited Newton-calibrated register count, the mutual-energy relative remainder is below 1e-18. '
      'The [response/noise/work identities](../model/openwave_gravity_noise_normalization/PROOF.md) give '
      'g_N=S Cov(Q,exp(Q))/(2π ω) and W Var(Q)=q²ωu/2. An independent covariance would fix the weak strength; '
      'the present count uses observed G, so its fluctuations are conditional and are not observed metric noise.','',
      'The [validated regular long-wavelength metric coefficient](../model/openwave_capacity_normalization_transfer/PROOF.md) '
      'is '+interval(v['capacity']['present_metric']['capacity']['metric_transfer_interval'])+
      ', versus '+interval(v['capacity']['present_metric']['lambda']['metric_transfer_interval'])+
      ' for the same-density Lambda control. The strict zero-wavenumber mode is a residual adiabatic coordinate mode with zero comoving density. '
      'The normalized metric ratio is not the finite-wavenumber matter growth ratio or a CMB observation. '
      'Absolute Hubble/density and primordial amplitude freedoms remain; the scalar and Einstein attachment are supplied.','',
      'These operational and spectral results add zero independent Nature-validation rows.','']
    return lines
