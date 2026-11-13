"""Authenticated projections of source-qualified spectra and field records."""
from fractions import Fraction
from pathlib import Path
import hashlib
import importlib
import json

PACKAGES = (
    'hadron_refinement_flow', 'strange_composites', 'light_nuclear_curve',
    'native_centre_string', 'dressed_core_response', 'dirac_maxwell_refinement',
    'autonomous_clock_channel', 'transverse_gravity', 'su2_volume_gap',
    'capacity_stress',
    'clock_phase_resource', 'quotient_homotopy', 'resolved_core',
    'dirac_field_connection', 'gravity_unit_constraint', 'horizon_capacity',
    'su2_propagation',
    'baryon_dispersion', 'cluster_bandwidth', 'nuclear_isospin',
    'baryon_poles', 'pole_weak_current', 'strange_shell_current',
    'coupled_haar_source', 'coupled_haar_spectrum', 'euclidean_refinement',
    'source_ward', 'einstein_error', 'radiation_capacity',
    'uniform_pole_bands', 'normalizable_particle_records', 'scattering_current_support',
    'haar_coupling_flow', 'cubic_haar_source', 'gauge_refinement_clock',
    'staggered_ward', 'replicated_geometry', 'capacity_perturbations',
    'full_strong_current', 'lambda_wavepacket', 'cluster_weyl_sequence',
    'refined_haar_source', 'refined_haar_spectrum', 'haar_memory_closure',
    'tensor_record_transport', 'collective_gravity_strength', 'capacity_transfer',
)


def decimal(value, upper=False, places=12):
    q = Fraction(value); scale = 10**places
    n = -((-q.numerator*scale)//q.denominator) if upper else q.numerator*scale//q.denominator
    sign = '-' if n < 0 else ''
    whole, part = divmod(abs(n),scale)
    return f'{sign}{whole}.{part:0{places}d}'


def interval(pair):
    return f'[{decimal(pair[0])}, {decimal(pair[1],upper=True)}]'


def difference_with_error(value, error):
    """Keep the printed centre-plus/minus interval outside the source interval."""
    value, error = Fraction(value), Fraction(error)
    if error < 0:
        raise ValueError('an error bound must be nonnegative')
    centre = f'{float(value):.10g}'
    radius = error + abs(value-Fraction(centre))
    return centre, decimal(radius,upper=True)


def collect(root):
    root=Path(root).resolve(); data={}; receipts={}
    for name in PACKAGES:
        relative=f'model/openwave_{name}/results.json'; path=root/relative
        value=json.loads(path.read_text())
        owner=importlib.import_module(f'model.openwave_{name}.__main__')
        expected=owner.pins()
        if not expected or value.get('source_sha256') != expected:
            raise ValueError('physical result inventory changed: '+name)
        if not str(value.get('verdict','')).startswith('PASS'):
            raise ValueError('physical result producer failed: '+name)
        for source,digest in expected.items():
            target=(root/source).resolve()
            if not target.is_relative_to(root) or hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
                raise ValueError('physical result source drift: '+source)
        data[name]=value; receipts[relative]=hashlib.sha256(path.read_bytes()).hexdigest()
    flow=data['hadron_refinement_flow']; strange=data['strange_composites']
    nuclear=data['light_nuclear_curve']; gravity=data['transverse_gravity']
    gap=data['su2_volume_gap']; capacity=data['capacity_stress']
    if capacity['openwave_row28_promotion'] or gravity['independent_energy_length_attachment']['physical_G_strength_match']:
        raise ValueError('unsupported cosmological or Newton-strength promotion')
    if not all(nuclear['mass_comparisons'][n]['matches'] is False for n in ('helion','alpha')):
        raise ValueError('independent helium failures omitted')
    if strange['native_nonleptonic_vertex_and_obstruction']['strong_channel_is_kinematically_open']:
        raise ValueError('closed strange-decay channel promoted')
    return {
        'receipts_sha256':receipts,
        'bulk_baryons':{'neutron_GeV_interval':flow['selected_bulk_neutron']['full_source_GeV_interval'],
                        'proton_GeV_interval':flow['selected_bulk_proton_GeV_interval'],
                        'axial_interval':flow['selected_bulk_neutron']['full_source_axial_interval'],
                        'ancestry':flow['calibration_ancestry'],'scope':flow['scope']},
        'strange_masses':[{k:r[k] for k in ('name','computed_GeV','comparison_GeV','comparison_role','difference_GeV')}
                          for r in strange['mass_rows']],
        'strange_ancestry':strange['calibration_ancestry'],
        'strange_decay_obstruction':{k:strange['native_nonleptonic_vertex_and_obstruction'][k] for k in
            ('strong_N4q_1antiquark_channel_floor_GeV','strong_Lambda_to_channel_gap_GeV',
             'gap_after_two_full_weak_norm_shifts_GeV','rigorous_channel_gap_above_Lambda_ground_GeV','interpretation')},
        'nuclear_masses':nuclear['mass_comparisons'],
        'nuclear_ancestry':nuclear['additional_calibration'],'nuclear_scope':nuclear['scope'],
        'transverse_photon_and_clocks':gravity['same_source_photon_comparison'],
        'independent_unit_attachment':gravity['independent_energy_length_attachment'],
        'gravity_scope':gravity['scope'],
        'pure_gauge_volume_limit':{'source':gap['source'],'conclusions':gap['proved_conclusions'],
                                  'physical_nontriviality':gap['physical_nontriviality'],
                                  'limits':gap['scope_limits']},
        'capacity_source':capacity['source'],
        'conditional_expansion_comparisons':[{k:r[k] for k in ('table','rows','anchors','sps_column',
            'fitted_chronometer_parameters','chi_square','chi_square_anchor_and_radiation_box',
            'max_radiation_shift_km_s_Mpc_upper','systematics_clamped_at_upper_endpoint_rows')}
            for r in capacity['comparisons']],
        'cosmology_calibration_boundary':capacity['calibration_boundary'],
        'native_scale_and_structure': scale_and_structure(data),
        'propagation_and_nuclei': propagation_and_nuclei(data),
        'pole_currents_and_full_field_records': pole_currents_and_full_fields(data),
        'normalizable_particles_and_spatial_fields': normalizable_particles_and_spatial_fields(data),
        'complete_dynamics_and_constraints': complete_dynamics_and_constraints(data),
        'all_sources_identified_as_one_physical_theory':False,
    }


def propagation_and_nuclei(data):
    baryon=data['baryon_dispersion']; strange=data['cluster_bandwidth']; nuclear=data['nuclear_isospin']
    if strange['native_weak_scope']['complete_outgoing_weak_decay_or_lifetime_proved']:
        raise ValueError('energy shell promoted to an executed strange decay')
    if nuclear['physical_nuclear_mass_spectrum_closed']:
        raise ValueError('failed nuclear comparisons promoted to a closed physical spectrum')
    return {
        'baryon': {'unit_dictionary':baryon['unit_dictionary'],
            'unmodified_relativistic_comparison':baryon['relativistic_comparison'],
            'unmodified_zone_edge_energy_gain_GeV':baryon['zone_edge_energy_gain_GeV'],
            'positive_correlated_propagation':{k:v for k,v in baryon['positive_correlated_propagation'].items()
                if k not in ('phase_records','proton_repaired_phase_records','complete_native_correlated_rule')}},
        'strange':{'source':strange['finite_range_source'],'bulk_energies_role':strange['bulk_energies_role'],
            'Lambda_rest_GeV_interval':strange['Lambda_rest']['full_Lambda_ground_GeV_interval'],
            'pion_rest':strange['pion_rest'],'kaon_rest':strange['kaon_rest'],
            'energy_shell':{k:v for k,v in strange['energy_shell_bracket'].items() if k not in ('below','above')},
            'locality_only_control':strange['locality_only_control'],'native_weak_scope':strange['native_weak_scope']},
        'nuclei':{k:nuclear[k] for k in ('exact_isospin_obstruction','declared_electrostatic_member',
            'complete_finite_spectra','infinite_volume_breakup_or_halo_proved','physical_nuclear_mass_spectrum_closed')},
    }


def scale_and_structure(data):
    """Project qualified results without duplicating large complete record tapes."""
    core=data['resolved_core']; gravity=data['gravity_unit_constraint']
    horizon=data['horizon_capacity']; ym=data['su2_propagation']
    if core['scope']['reference_bound_energies_are_observed_particle_poles']:
        raise ValueError('reference-bound core promoted to physical particle pole')
    if gravity['scope']['physical_G_and_c_match'] or horizon['physical_horizon_or_Lambda_amplitude_derived']:
        raise ValueError('unsupported physical gravity/cosmology identification')
    if ym['scope']['Clay_problem_solved'] or ym['scope']['physical_GeV_gap']:
        raise ValueError('fixed-spacing gauge result promoted to continuum or GeV gap')
    return {
        'phase': {k:data['clock_phase_resource'][k] for k in ('general_theorem','scope')},
        'quotient': {k:data['quotient_homotopy'][k] for k in ('continuous_native_lift','scope')},
        'core': {'architecture':core['architecture'],'selection':core['regulator_selection'],
            'scope':core['scope'],'members':[
                {'family':m['family'],'energy_reference':m['energy_reference'],
                 'reference_bound_energy_GeV_interval':m['reference_bound_energy_GeV_interval_exact'],
                 'radius_cells_interval':m['complete_core_probe']['radius_x_cells_interval_exact'],
                 'maximum_suppression_interval':m['complete_core_probe']['largest_first_Brillouin_zone_suppression_interval_exact'],
                 'half_suppression_wavevector_over_pi_interval':m['complete_core_probe']['half_suppression_wavevector_over_pi_interval_exact'],
                 'static_scalar_probe':{k:v for k,v in m['same_source_scalar_probe_coupling'].items()
                                        if k!='complete_interventions'}} for m in core['members']]},
        'dirac_field': {k:data['dirac_field_connection'][k] for k in
            ('source','exact_connection','correct_transformed_field','all_flux_blocks_retained_analytically','scope')},
        'unit_constraint': {k:gravity[k] for k in
            ('formula','light_matched_G_interval','required_stiffness_ratio_interval_exact',
             'energy_unit_required_if_electron_anchor_abandoned_GeV','diagnostic_float',
             'external_reference_inputs','scope')},
        'capacity_histories': {'anchors':horizon['anchors'],
            'branch_consequences':horizon['branch_consequences'],
            'comparisons':[{k:v for k,v in row.items() if k!='records'} for row in horizon['comparisons']],
            'physical_horizon_or_Lambda_amplitude_derived':horizon['physical_horizon_or_Lambda_amplitude_derived'],
            'same_data_alternatives_are_independent':horizon['same_data_alternatives_are_independent']},
        'full_Haar_dynamics':{k:ym[k] for k in ('source','commutator_bound','infinite_volume_error',
            'full_flux_propagation_bounds','SU_N_full_Haar_maps','proved','scope')},
    }


def render(value):
    b=value['bulk_baryons']; lines=['## Spatial spectra and common field records','',
        'The [complete path refinement](../model/openwave_hadron_refinement_flow/PROOF.md) '
        'uses two degree-dependent coefficients selected with the inherited neutron and axial targets. '
        'They vanish on the complete two-cell action. Exact compact reduction and the full infinite '
        'relative-position complement give these one-dimensional calibrated rest bands:', '',
        '| Native assignment | Energy / GeV |', '| --- | ---: |',
        f"| neutron band | {interval(b['neutron_GeV_interval'])} |",
        f"| proton band | {interval(b['proton_GeV_interval'])} |", '',
        f"The native axial interval is {interval(b['axial_interval'])}; the vector matrix element is exactly one. "
        'These source bands are not identified three-dimensional scattering poles.', '',
        '[Strange-family spectra](../model/openwave_strange_composites/PROOF.md) retain the same electron unit '
        'and light two-cell calibration domains. Two additional masses select the strange coefficients; '
        'a separate supplied Cabibbo magnitude enters current interventions.', '',
        '| Native assignment | Computed / GeV | Compared / GeV | Role |', '| --- | ---: | ---: | --- |']
    for r in value['strange_masses']:
        lines.append(f"| {r['name']} | {r['computed_GeV']:.12g} | {r['comparison_GeV']} | {r['comparison_role']} |")
    o=value['strange_decay_obstruction']
    lines += ['', 'The full two-cell four-quark/one-antiquark channel has energy at least '
              f"{decimal(o['strong_N4q_1antiquark_channel_floor_GeV'])} GeV, above the Lambda input by "
              f"{decimal(Fraction(o['strong_Lambda_to_channel_gap_GeV'])*1000)} MeV. "
              f"The actual channel gap above the Lambda ground is at least {decimal(Fraction(o['rigorous_channel_gap_above_Lambda_ground_GeV'])*1000)} MeV; "
              'it survives the full weak norm bound. The actual nonzero '
              'native nonleptonic vertex is a virtual response, without an on-shell strange decay.', '',
              'The [finite nuclear source](../model/openwave_light_nuclear_curve/PROOF.md) adds one '
              'occupancy-gated triton calibration and preserves all operators with at most six quarks. '
              'Both unused helium comparisons fail.', '',
              '| Native assignment | Source interval / GeV | Compared / GeV | Role |',
              '| --- | ---: | ---: | --- |']
    for name,r in value['nuclear_masses'].items():
        lines.append(f"| {name} | {interval(r['source_GeV_interval'])} | {r['target_GeV']} | {r['role']} |")
    g=value['transverse_photon_and_clocks']; units=value['independent_unit_attachment']
    gravity_rows=[]
    for label,centre_key,error_key in (
        ('mean arrival parameter','conditional_arrival_delay','arrival_delay_error_upper'),
        ('transverse position / cell','transverse_position_difference','position_difference_error_upper'),
        ('Fourier direction / rad','separate_Fourier_direction_deflection','direction_error_upper'),
    ):
        centre,error=difference_with_error(g[centre_key],g[error_key])
        gravity_rows.append(f'| {label} | {centre} | {error} |')
    lines += ['', value['nuclear_scope'], '',
        'One [transverse Maxwell, clock and quantum-geometry source](../model/openwave_transverse_gravity/PROOF.md) '
        'retains all 540 divergence-free field coordinates, both polarization controls, the photon vacuum, '
        'and every detector failure and work record. The declared constitutive source gives:', '',
        '| Observable | Difference | Certified numerical error bound |', '| --- | ---: | ---: |',
        *gravity_rows, '',
        'These are finite source observables, without a physical lensing-angle or seconds attachment. '
        f"The separate electron-energy/strong-length map gives G={units['conditional_G_natural_GeV_minus2']:.9g} GeV^-2 "
        f"({units['conditional_to_observed_strength_ratio']:.5g} times observed) and photon speed "
        f"{units['conditional_photon_speed_over_observed_c']:.9g}c. Both unit comparisons fail. "
        'No common strong/geometry Hamiltonian transfer is claimed.', '',
        'The [full-Haar pure SU(2) source](../model/openwave_su2_volume_gap/PROOF.md) has a gap uniform '
        'in periodic spatial volume for sufficiently strong coupling at fixed spacing, using Yarotsky\'s '
        'stability theorem and local norm 64/(g²)². The threshold is existential; this is not a numerical '
        'admission of g²=8 or a spatial-continuum Clay result.', '',
        'The [capacity-volume source](../model/openwave_capacity_stress/PROOF.md) derives mechanical '
        'w=-1 at fixed active capacity and accounts for every capacity-change work stroke. Its expansion '
        f"comparison supplies H0={value['conditional_expansion_comparisons'][0]['anchors']['H0_km_s_Mpc']} "
        f"and Omega_m={value['conditional_expansion_comparisons'][0]['anchors']['Omega_m']} plus the flat Friedmann branch. It does not "
        'derive Lambda\'s amplitude. The complete correlated fifteen-point BC03 statistic is '
        f"{interval(value['conditional_expansion_comparisons'][0]['chi_square'])}; "
        'the alternative M11 reduction gives '
        f"{interval(value['conditional_expansion_comparisons'][1]['chi_square'])}. "
        'These shared-data alternatives are not independent likelihoods.', '',
        'The [native centre-string](../model/openwave_native_centre_string/PROOF.md) identifies exact '
        'dual-electric phase jumps, with a separate continuous-texture interpretation boundary. '
        '[Complete pump/loss clocks](../model/openwave_autonomous_clock_channel/PROOF.md) retain exact '
        'native frequency, finite linewidth, batteries and the probe-created-phase control. The '
        '[Pauli/Maxwell refinement](../model/openwave_dirac_maxwell_refinement/PROOF.md) retains the '
        'magnetic commutator and all integer fluxes. Its [dressed current/core response](../model/openwave_dressed_core_response/PROOF.md) '
        'at Wilson coefficient 1/16 varies with scale, while that member has an unresolved spatial core. '
        'These sources are not pooled into one identified physical theory.', '']
    lines += render_propagation_and_nuclei(value['propagation_and_nuclei'])
    lines += render_scale_and_structure(value['native_scale_and_structure'])
    lines += render_complete_dynamics_and_constraints(value['complete_dynamics_and_constraints'])
    return lines


def render_propagation_and_nuclei(value):
    b=value['baryon']['positive_correlated_propagation']; s=value['strange']; n=value['nuclei']
    lines=['### Propagating rest bands and a local strange energy shell','',
        'The [positive correlated kinetic source](../model/openwave_baryon_dispersion/PROOF.md) uses '
        'every native compact colour-singlet triquark channel, with one inverse-rest-operator rule '
        'J = zeta M^-1. The physical speed, electron unit and strong length are inputs. Its complete '
        'three-quark rest fibre is unchanged; the neutron inertial mass is enclosed by '
        f"{interval(b['repaired_inertial_mass_GeV_interval'])} GeV and its ground-band width by "
        f"{interval(b['available_ground_bandwidth_GeV_interval'])} GeV. "
        'Full-Fock positivity does not transfer the higher-number spectra or prove exact Lorentz dispersion.', '',
        'The [local strange extension](../model/openwave_cluster_bandwidth/PROOF.md) preserves the '
        'complete two-cell action and makes separated colour-singlet clusters additive. With the same '
        'positive kinetic rule, its complete bulk Lambda energy is '
        f"{interval(s['Lambda_rest_GeV_interval'])} GeV. The proton and pion bands cross that energy "
        'at opposite momenta within '
        f"{interval(s['energy_shell']['opposite_momenta_GeV_interval'])} GeV. "
        'This establishes a strong-sector energy shell. The native weak current has no proved transfer '
        'to complete outgoing scattering states here; a strange decay rate is not supplied. The '
        'bulk rest shifts and inherited mass calibration are retained.', '',
        '| Local bulk assignment | Complete energy interval / GeV |',
        '| --- | ---: |',
        f"| Lambda | {interval(s['Lambda_rest_GeV_interval'])} |",
        f"| charged pion | {interval(s['pion_rest']['full_ground_GeV_interval'])} |",
        f"| charged kaon | {interval(s['kaon_rest']['full_ground_GeV_interval'])} |", '',
        '### Nuclear charge breaking','',
        'The [exact isospin identity](../model/openwave_nuclear_isospin/PROOF.md) proves that another '
        'isospin-scalar interaction cannot repair both mirror masses. The required charge-breaking '
        'gap is '
        f"{decimal(Fraction(n['exact_isospin_obstruction']['additional_charge_breaking_gap_needed_GeV'])*1000)} MeV. "
        'A declared collinear Coulomb pair member uses inherited electromagnetic and length inputs, '
        'with no additional nuclear mass fit. Complete spin and high-complement bounds give:', '',
        '| Assignment | Coulomb-member energy / GeV | Compared / GeV | Target contained |',
        '| --- | ---: | ---: | --- |']
    for row in n['complete_finite_spectra']:
        lines.append(f"| {row['name']} | {interval(row['ground_GeV_interval'])} | {row['target_GeV']} | {row['target_in_interval']} |")
    lines += ['', 'The helion and alpha comparisons fail. This finite electrostatic member has no '
        'asserted dynamical Gauss-field realization, infinite-volume breakup theorem or transfer of '
        'the correlated kinetic source to its nuclear spectra.', '']
    return lines


def render_scale_and_structure(value):
    core=value['core']; units=value['unit_constraint']; histories=value['capacity_histories']
    lines=['### Resolved charge profiles and their energies','',
        'The [full-rotor charge source](../model/openwave_resolved_core/PROOF.md) uses Wilson coefficient '
        f"{core['architecture']['Wilson_coefficient_exact']}, selected from the complete displayed dyadic scan "
        'by a half-suppression resolution test. No measured radius or running-charge value selects it. '
        'The same static charge intervention and independent position instrument establish the scale join. '
        'All integer fluxes, opposite charges, spin partners and preparation failures are retained.', '',
        '| Native family | Vacuum-subtracted energy / GeV | Radius / cell | Maximum charge suppression |',
        '| --- | ---: | ---: | ---: |']
    for row in core['members']:
        lines.append(f"| {row['family']} | {interval(row['reference_bound_energy_GeV_interval'])} | "
                     f"{interval(row['radius_cells_interval'])} | {interval(row['maximum_suppression_interval'])} |")
    lines += ['',
        'These reference-bound energies subtract the true full-rotor vacuum and retain the compensating '
        'quantum reference. They are not observed particle poles. The resolved profile and static effective '
        'charge do not establish a universal renormalized QED or QCD charge. The changed regulator does '
        'not inherit the mass or response fits of a different action.', '',
        'The [conditional clock theorem](../model/openwave_clock_phase_resource/PROOF.md) proves that '
        'an energy-invariant source and battery, an energy-conserving joint unitary and energy-compatible '
        'recorded effects leave every conditional source state energy-invariant. Degenerate relational '
        'coherence is allowed; it has zero free native frequency. A separately declared interaction '
        'frequency is not thereby a rest mass.', '',
        'The [continuous quotient lift](../model/openwave_quotient_homotopy/PROOF.md) identifies '
        '(SU(3) × U(1))/Z3 with U(3), proves fundamental group Z and closes the native fractional-charge '
        'generator with its colour phase. The electric eigenstates have uniform conjugate angle density; '
        'this is not a measured localized classical texture or an additional mass fit.', '',
        'The [dynamical Dirac field member](../model/openwave_dirac_field_connection/PROOF.md) retains '
        'both signed bands and every integer rotor flux. Diagonalizing its matter matrix gives the '
        'connection A = −i U†∂U, and the field energy becomes α(E+A)². Omitting A changes the spectrum. '
        'This analytical comparison member has no asserted equivalence to the full native Fock action.', '',
        '### Common units and expansion histories','',
        'The [exact unit elimination](../model/openwave_gravity_unit_constraint/PROOF.md) gives '
        'G/(c_model/c) = g/u². With unchanged electron energy unit u, length alone cannot match '
        'both light speed and Newton strength. Matching light gives G in '
        f"{interval(units['light_matched_G_interval'])} GeV^-2. "
        'The required stiffness change is about '
        f"{units['diagnostic_float']['required_stiffness_ratio']:.6g}; no such gravity calibration is made.", '',
        'The [capacity histories](../model/openwave_horizon_capacity/PROOF.md) compare the same fifteen '
        f"chronometer rows using supplied H0={histories['anchors']['H0_km_s_Mpc']} and "
        f"Omega_m={histories['anchors']['Omega_m']}. No parameter is fitted to these chronometers. "
        'The covariance and alternative stellar-population reductions are shared data, not independent likelihoods.', '',
        '| Capacity branch | Data reduction | Correlated chi-square interval |',
        '| --- | --- | ---: |']
    for row in histories['comparisons']:
        lines.append(f"| {row['branch']} | {row['table']} | {interval(row['chi_square_interval'])} |")
    lines += ['',
        'Fixed capacity gives w=−1. Comoving area gives w=−1/3. A Hubble-area branch with separately '
        'conserved dust gives w=0; retaining mechanical pressure instead requires energy exchange. '
        'For the future-event-area branch with gamma=1, a monotone rational equation replaces ODE '
        'integration and proves the future boundary condition. The physical horizon and cosmological '
        'amplitude are supplied choices, and the high-redshift dust control is not an early-Universe fit.', '',
        '### Complete gauge-field propagation','',
        'The [full-Haar propagation bound](../model/openwave_su2_propagation/PROOF.md) retains unbounded '
        'onsite electric energy and every representation. Bounded plaquette interactions give local '
        'finite-volume dynamics converging in norm, with explicit boundary errors. This propagation '
        'argument applies to compact gauge groups. For compact connected semisimple groups, the '
        'electric Casimir gap also supplies the fixed-spacing strong-coupling stability map; its '
        'coupling threshold is existential. At fixed coupling all energy gaps scale as 1/a, so a '
        'finite physical continuum gap requires a further scaling construction. Neither a physical '
        'GeV glueball mass nor the four-dimensional Yang-Mills Clay theorem follows.', '']
    return lines


def pole_currents_and_full_fields(data):
    """Typed source results; physical validation and calibration counts stay fixed."""
    poles=data['baryon_poles']; weak=data['pole_weak_current']; strange=data['strange_shell_current']
    source=data['coupled_haar_source']; spectrum=data['coupled_haar_spectrum']; euclidean=data['euclidean_refinement']
    ward=data['source_ward']; error=data['einstein_error']; radiation=data['radiation_capacity']
    forbidden=[
        poles['physical_3plus1_Lorentz_particle_pole_identification'],poles['mass_refit_performed'],
        weak['autonomous_weak_decay_or_neutron_lifetime_proved'],weak['physical_3plus1_current_or_particle_dictionary_proved'],
        weak['new_empirical_current_or_mass_input'],strange['physical_on_shell_decay_amplitude_proved'],
        strange['physical_Lambda_lifetime_proved'],strange['new_mass_or_current_calibration'],
        strange['complete_component_probe']['normalizable_wave_packet_transfer_executed'],
        source['scope']['four_dimensional_Yang_Mills_reconstruction'],source['scope']['physical_GeV_mass'],
        spectrum['scope']['Clay_problem_solved'],spectrum['scope']['spatial_continuum_Yang_Mills'],
        euclidean['scope']['four_dimensional_interacting_continuum_constructed'],euclidean['scope']['Clay_problem_solved'],
        ward['scope']['physical_metric_assignment_derived'],ward['scope']['nonlinear_einstein_equation_derived'],
        ward['scope']['physical_Newton_strength_matched'],error['scope']['continuum_or_nonlinear_einstein_theorem'],
        error['scope']['cosmological_constant_fixed_by_null_records'],error['scope']['physical_stress_identification_proved'],
        radiation['scope']['independent_observational_validation'],radiation['scope']['cosmological_amplitude_predicted'],
        radiation['scope']['early_universe_abundances_calculated'],radiation['scope']['standard_model_thermal_history_derived'],
    ]
    if any(x is not False for x in forbidden):
        raise ValueError('source result promoted to unproved physical validation')
    if strange['complete_component_probe']['probabilities_are_fixed_total_momentum_fibre_quantities'] is not True:
        raise ValueError('strange-current fixed-fibre domain omitted')
    return {
        'physical_comparison_count_added':0,
        'baryon_poles':{'rest_records':[
            {'assignment':r['assignment'],'source_pole_GeV_interval':r['source_pole_GeV_interval'],
             'prepared_trial_pole_weight_lower':r['pole']['prepared_trial_pole_weight_lower'],
             'original_input_pole_event_probability_lower':r['pole']['original_input_pole_event_probability_lower']}
            for r in poles['rest_source_poles']],
            'momentum_fibre_count':len(poles['momentum_fibre_records']),
            'matrix_spectral_measure':poles['matrix_spectral_measure'],
            'complete_observer_instrument':poles['complete_observer_instrument'],
            'input_ancestry':poles['input_ancestry'],'scope':poles['scope']},
        'pole_weak_current':{'momentum_records':[
            {'cosine':r['cosine'],'sine':r['sine'],'currents':{
                k:{name:v[name] for name in ('actual_source_pole_amplitude_magnitude_interval',
                  'current_success_then_proton_pole_probability_interval','certificate')}
                for k,v in r['currents'].items()}}for r in weak['momentum_records']],
            'complete_native_instrument':weak['complete_native_instrument'],
            'ancestry':weak['ancestry'],'scope':weak['scope']},
        'strange_current':{'Lambda_source_pole_GeV_interval':strange['Lambda_source_pole']['full_ground_GeV_interval'],
            'exact_native_vertex':strange['exact_native_vertex'],
            'full_pole_to_complete_compact_N5_response':strange['full_pole_to_complete_compact_N5_response'],
            'complete_component_probe':strange['complete_component_probe'],
            'same_source_energy_shell':strange['same_source_energy_shell'],
            'inherited_calibrations':strange['inherited_calibrations'],
            'remaining_premise':strange['remaining_premise'],'scope':strange['scope']},
        'coupled_Haar_source':{k:source[k]for k in ('graph','source','physical_basis','shared_electric_diagonal',
            'magnetic_operator','diagnostic_caps','complete_source_properties','scope')},
        'coupled_Haar_spectrum':{'source':spectrum['source'],'empirical_fit_inputs':spectrum['empirical_fit_inputs'],
            'records':[{k:r[k]for k in ('g_squared_exact','spacing_exact','shared_edge',
                'diagnostic_max_twice_loop_spin','retained_complete_dimension','full_Haar_gap_interval_exact',
                'connected_Wilson_covariance_interval_exact','Wilson_variance_interval_exact',
                'complete_joint_observer_instrument','scope')}for r in spectrum['records']],
            'physical_interpretation':spectrum['physical_interpretation'],'scope':spectrum['scope']},
        'Euclidean_response':{'representation_refinement_records':[
            {k:r[k]for k in ('proof_coordinate_key','source_Haar_representations_truncated',
                'vacuum_weight_retained','correlations','finite_positive_reflection_form_interval_exact')}
            for r in euclidean['representation_refinement_records']],
            'fixed_coupling_spacing_trajectory':euclidean['fixed_coupling_spacing_trajectory'],
            'continuum_acceptance':euclidean['continuum_acceptance'],'scope':euclidean['scope']},
        'native_lapse_and_Ward':{k:ward[k]for k in ('source_q_exact','source_stiffness_exact','source_sites',
            'boundary_convention','symbolic','native_clock_response','finite_geometry','scope')},
        'finite_Einstein_error':{k:error[k]for k in ('sharp_maximum_entry_error_constant',
            'operator_norm_error_upper_constant','example_uniform_charge_error_exact',
            'maximum_tensor_entry_error_bound_interval','ward_transport','scope')},
        'radiation_capacity':{k:radiation[k]for k in ('inputs','conservation_identity','equation_of_state',
            'integration_certificate','future_event_branch','fixed_capacity_branch','hubble_capacity_control',
            'early_scaling','scope')},
    }


def render_pole_currents_and_full_fields(value):
    poles=value['baryon_poles']; weak=value['pole_weak_current']; strange=value['strange_current']
    lines=['### Native source poles and their current responses','',
        'The [complete interacting baryon source](../model/openwave_baryon_poles/PROOF.md) '
        'has isolated momentum-fibre poles, with the electron unit and inherited strong/axial '
        'calibrations retained. The rational prepared trials have the following certified pole weights; '
        'these weights are spectral residues, not percentage agreements with Nature.', '',
        '| Source assignment | Pole energy / GeV | Prepared trial pole weight at least |',
        '| --- | ---: | ---: |']
    for row in poles['rest_records']:
        lines.append(f"| {row['assignment']} | {interval(row['source_pole_GeV_interval'])} | {decimal(row['prepared_trial_pole_weight_lower'])} |")
    rest=next(r for r in weak['momentum_records']if r['cosine']=='1'and r['sine']=='0')
    lines += ['',
        'The [same native CAR currents](../model/openwave_pole_weak_current/PROOF.md) give rest vector '
        f"amplitude {interval(rest['currents']['vector']['actual_source_pole_amplitude_magnitude_interval'])} "
        f"and axial amplitude {interval(rest['currents']['axial']['actual_source_pole_amplitude_magnitude_interval'])}. "
        'Every tested momentum has positive certified response. Preparation failure, current-produced '
        'channels, the charge reference and source work remain recorded. These fixed-spacing source '
        'poles are not identified four-dimensional relativistic particles or an autonomous beta-decay lifetime.', '',
        'The [native strange current](../model/openwave_strange_shell_current/PROOF.md) satisfies '
        '6 c² G <= W*W <= 7 c² G on every compact Lambda state, with the inherited Cabibbo factor '
        'c=sin(theta)cos(theta). Its complete compact five-body response norm squared is '
        f"{interval(strange['full_pole_to_complete_compact_N5_response']['amplitude_norm_squared_interval'])}; "
        'all noncompact outputs and probe failure are retained. The same strong source has an '
        'opposite proton–pion shell at '
        f"{interval(strange['same_source_energy_shell']['opposite_momentum_GeV_interval'])} GeV. "
        'The compact-current probabilities are total-momentum-zero relative-fibre quantities; transfer '
        'to a normalizable center-of-mass wave packet is not executed. A nonzero compact output is '
        'not the interacting on-shell scattering projector, and no Lambda decay lifetime follows.', '',
        '### Coupled full-Haar gauge records','',
        'The [shared-link two-plaquette source](../model/openwave_coupled_haar_source/PROOF.md) retains '
        'every admissible spin-network channel and the entire representation complement. Its '
        '[exact spectral bounds](../model/openwave_coupled_haar_spectrum/PROOF.md) use declared '
        'coupling and spacing, with no empirical fit:', '',
        '| g² | spacing | maximum twice-spin proof coordinate | Gap / source unit | Connected Wilson covariance |',
        '| ---: | ---: | ---: | ---: | ---: |']
    for row in value['coupled_Haar_spectrum']['records']:
        if row['shared_edge']:
            lines.append(f"| {row['g_squared_exact']} | {row['spacing_exact']} | {row['diagnostic_max_twice_loop_spin']} | {interval(row['full_Haar_gap_interval_exact'])} | {interval(row['connected_Wilson_covariance_interval_exact'])} |")
    e=value['Euclidean_response']; row=e['representation_refinement_records'][-1]
    quarter=next(r for r in row['correlations']if r['Euclidean_time_exact']=='1/4')
    lines += ['',
        'The independent-loop control retains zero covariance. Complete joint Wilson instruments '
        'retain all four outcomes and a nonzero apparatus-work budget. The '
        '[same-source Euclidean response](../model/openwave_euclidean_refinement/PROOF.md) at time 1/4 is '
        f"{interval(quarter['full_source_connected_correlation_interval_exact'])}. "
        'A fixed-coupling shrinking-spacing sequence has gap proportional to 1/a and collapses '
        'positive-time correlations in its unchanged clock. These are finite-graph full-Haar results, '
        'without a physical GeV glueball, a four-dimensional Yang–Mills continuum or a Clay-prize claim.', '',
        '### Native clock response, finite tensor error and radiation','',
        'The [actual weak-probe log lapse](../model/openwave_source_ward/PROOF.md) has source-energy '
        f"derivative {value['native_lapse_and_Ward']['native_clock_response']['formula']}, "
        'normalized by the nonunit vacuum lapse cosh(qF). Commuting centered differences give an exact '
        'linearized Bianchi identity. The metric/stress dictionary and physical Newton strength remain inputs or open identifications.', '',
        'The [nine-null tensor decoder](../model/openwave_einstein_error/PROOF.md) has sharp entrywise '
        f"error constant {value['finite_Einstein_error']['sharp_maximum_entry_error_constant']} "
        f"and operator-norm upper constant {value['finite_Einstein_error']['operator_norm_error_upper_constant']}. "
        'Ward errors propagate only on the actual difference graph; the retained centered-difference '
        'checkerboard prevents a nearest-neighbour constancy claim. Null data leave the metric multiple undetermined.', '']
    radiation=value['radiation_capacity']; last=radiation['future_event_branch'][-1]; inputs=radiation['inputs']
    lines += ['The [ideal-radiation capacity branch](../model/openwave_radiation_capacity/PROOF.md) uses supplied '
        f"H0={decimal(inputs['H0_km_s_Mpc'])}, Omega_m={decimal(inputs['Omega_m0'])} and "
        f"Omega_r={decimal(inputs['Omega_r0'])}. Exact outward Taylor bounds reach z={last['redshift']}, "
        f"where the event-area dark fraction lies in {interval(last['Omega_dark_interval'])}. "
        'The radiation-era fraction scales as a² for event-area capacity and a⁴ for fixed capacity. '
        'No cosmological amplitude, horizon selection, Standard Model thermal history, early abundance '
        'or independent observational validation is added. These source-qualified records add zero '
        'Nature-validation rows.', '']
    return lines



def normalizable_particles_and_spatial_fields(data):
    """Project complete new source records without adding empirical score rows."""
    bands=data['uniform_pole_bands']; packets=data['normalizable_particle_records']
    strange=data['scattering_current_support']; flow=data['haar_coupling_flow']
    cube=data['cubic_haar_source']; block=data['gauge_refinement_clock']
    ward=data['staggered_ward']; geometry=data['replicated_geometry']; capacity=data['capacity_perturbations']
    forbidden=(
        bands['physical_four_dimensional_particle_dictionary_proved'],
        packets['physical_four_dimensional_particle_identification'],
        packets['autonomous_beta_decay_or_lifetime_proved'],
        packets['complete_protocol']['finite_local_compilation_of_exact_band_filter_proved'],
        strange['physical_Lambda_decay_amplitude_proved'],strange['physical_Lambda_lifetime_proved'],
        strange['complete_Lambda_pole_transfer']['normalizable_Lambda_packet_executed'],
        strange['interacting_spectral_test']['Moller_wave_operator_constructed'],
        strange['interacting_spectral_test']['full_H5_on_shell_spectral_density_computed'],
        strange['actual_kinetic_substep']['complete_autonomous_strong_Hamiltonian_used'],
        flow['Yang_Mills_continuum_constructed'],cube['physical_gap_GeV_or_continuum_identification'],
        cube['uniform_spatial_volume_or_refinement_limit_proved'],
        block['scope']['autonomous_full_Hamiltonian_RG_derived'],block['scope']['continuum_noncollapse_proved'],
        ward['scope']['co_located_nine_null_observer_instrument_derived'],ward['scope']['physical_metric_or_stress_derived'],
        geometry['scope']['physical_Newton_strength_matched'],geometry['scope']['physical_Einstein_spacetime'],
        capacity['canonical_scalar_completion']['source_scalar_transfer_proved'],
        capacity['scope']['new_observational_match'],capacity['scope']['CMB_or_structure_likelihood_executed'],
    )
    if any(forbidden):raise ValueError('unproved spatial-field physical promotion')
    if (bands['mass_or_current_parameter_refitted']or packets['new_empirical_parameter_or_mass_fit']
            or strange['new_mass_or_current_calibration']or strange['new_independent_match_with_nature']
            or flow['input_nature_data_used']or cube['new_empirical_fit_inputs']
            or geometry['scope']['new_empirical_calibration']or capacity['scope']['new_empirical_fit']):
        raise ValueError('new calibration or validation changed the scoped source result')
    if (not packets['normalizable_packet']['carrier_is_infinite_chain']
            or packets['normalizable_packet']['finite_ring_or_delta_momentum_substitution']
            or not strange['actual_kinetic_substep']['full_infinite_N5_kinetic_carrier_retained']
            or not block['scope']['higher_coarse_operator_remainder_retained']
            or block['physical_clock']['same_source_Hamiltonian_time_rescaled']):
        raise ValueError('complete carrier or retained-clock domain omitted')
    pick=lambda row,keys:{k:row[k]for k in keys}
    return {
        'physical_comparison_count_added':0,
        'uniform_poles':{**pick(bands,('circle_cover','finite_theorem_shortcut','ancestry','scope')),
            'members':[{k:v for k,v in row.items()if k!='records'}for row in bands['members']]},
        'normalizable_packets':pick(packets,('normalizable_packet','members','complete_protocol',
            'uniform_native_weak_packet_records','ancestry','scope')),
        'strange_composite_current':{**pick(strange,('complete_Lambda_pole_transfer','complete_probe_records',
            'bare_separated_support_control','actual_kinetic_substep','same_source_open_energy_shell',
            'interacting_spectral_test','inherited_inputs','scope')),
            'channel':pick(strange['exact_CAR_channel'],('independent_composite_rank','raw_occupation_union_dimension',
                'generalized_eigenvalues_exact','generalized_eigenvalue_bounds_exact','Gram_order_exact',
                'complete_Pi_C_and_I_minus_Pi_C_instrument','fitted_mass_eigenvectors_used_to_define_channel'))},
        'Haar_coupling_response':pick(flow,('declared_couplings_exact','records','finite_coupling_steps',
            'exact_derivative_identity','physical_scope','spatial_refinement_performed')),
        'cubic_Haar_source':{**pick(cube,('graph','exact_source_checks','source','cube_symmetry',
            'representation_projection_is_source_cutoff','physical_gap_GeV_or_continuum_identification',
            'uniform_spatial_volume_or_refinement_limit_proved')),
            'records':[{k:v for k,v in row.items()if k!='complete_finite_spectrum'}for row in cube['records']]},
        'spatial_Haar_blocking':pick(block,('geometry','embedding','electric_compression','first_magnetic_compression',
            'exact_conditional_magnetic_expansion_through_order_four','uniform_remainder','declared_magnetic_steps',
            'physical_clock','scope')),
        'staggered_Ward':pick(ward,('construction','exact_identities','error_transport','complete_periodic_controls','scope')),
        'replicated_geometry':pick(geometry,('architecture','uniform_theorems','complete_register_records',
            'complete_ground_to_classical_limit','preserved_far_field_log_lapse_slope_interval','scope')),
        'capacity_scalar':pick(capacity,('inputs','uniform_past_certificate','required_completion',
            'canonical_scalar_completion','complete_finite_acoustic_control','scope')),
    }


def render_normalizable_particles_and_spatial_fields(value):
    bands=value['uniform_poles']; packets=value['normalizable_packets']; strange=value['strange_composite_current']
    lines=['### Normalizable source particles and composite currents','',
        'The [complete momentum-circle certificate](../model/openwave_uniform_pole_bands/PROOF.md) '
        'retains both high compact coordinates and every infinite relative configuration. '
        'Its three-field pole residue is positive at every momentum. The '
        '[normalizable infinite-chain packet](../model/openwave_normalizable_particle_records/PROOF.md) '
        'uses an explicit L² amplitude, with separate complete position and native-current records. '
        'Its energies below are packet mean energies including kinetic energy, not new rest masses. '
        'The electron unit and strong/axial calibration remain inherited.', '',
        '| Native assignment | All-momentum gap lower / source unit | Packet band-probability lower | Packet mean energy / GeV |',
        '| --- | ---: | ---: | ---: |']
    for b,p in zip(bands['members'],packets['members']):
        assert b['assignment']==p['assignment']
        lines.append(f"| {b['assignment']} | {decimal(b['all_momentum_gap_lower_source'])} | "
                     f"{decimal(p['uniform_band_projector_probability_lower'])} | "
                     f"{interval(p['complete_packet_trial_mean_energy_GeV_interval'])} |")
    lines += ['', 'The exact band filter is a declared nonlocal spectral service; finite local compilation '
        'and a physical four-dimensional particle dictionary remain open. Preparation, band, position '
        'and current complements, reference charge and work are retained.', '',
        'The [native proton–pion current](../model/openwave_scattering_current_support/PROOF.md) identifies '
        f"{strange['channel']['independent_composite_rank']} independent antisymmetrized composite fields. "
        'Its projected-current Gram eigenvalues are 6/5, 241/150 and (277±sqrt(54409))/200, '
        'before the inherited Cabibbo factor. The positive lower bound transfers to the complete '
        'k=0 Lambda source pole. A separate specified localized Lambda field has a rigorously '
        'positive disjoint-composite event after one full native kinetic substep; its original-input '
        'probability, including the probe success cost, is at least '
        f"{decimal(strange['actual_kinetic_substep']['original_input_probe_then_kinetic_separated_probability_lower'],places=16)}. "
        'This is not the complete autonomous strong evolution. A normalizable Lambda packet, the full '
        'interacting outgoing wave operator, on-shell overlap and physical lifetime remain open.', '',
        '### Complete cubic gauge records and actual spatial blocking','',
        'The [three-dimensional spatial cube](../model/openwave_cubic_haar_source/PROOF.md) retains all '
        'Haar representations and six plaquettes. Exact Gauss/Clebsch–Gordan contraction, a complete '
        'spectral frame and every outside trial successor give these source intervals:', '',
        '| Declared g² | Full source gap | Adjacent-face covariance envelope |', '| --- | ---: | ---: |']
    for row in value['cubic_Haar_source']['records']:
        adjacent=[tuple(map(Fraction,r['connected_covariance_interval_exact']))for r in row['all_fifteen_connected_covariances']if not r['opposite']]
        envelope=(min(r[0]for r in adjacent),max(r[1]for r in adjacent))
        lines.append(f"| {row['g_squared_exact']} | {interval(row['full_source_gap_interval_exact'])} | {interval(envelope)} |")
    block=value['spatial_Haar_blocking']
    lines += ['', 'All fifteen covariance pairs, 64 joint outcomes and the nonselective work budget remain. '
        'These are finite-graph source gaps, without physical GeV glueball or continuum identification. '
        'The [seven-coupling response](../model/openwave_haar_coupling_flow/PROOF.md) uses bounded magnetic '
        'expectations in the exact Feynman–Hellmann identity; its fixed-gap tangent is not an RG beta function.', '',
        'The [spatial Haar block](../model/openwave_gauge_refinement_clock/PROOF.md) adds '
        f"{block['geometry']['fine_links']} fine links and {block['geometry']['fine_faces']} faces in total. "
        'The conditional magnetic factor has exact expansion '
        '1+(9/2)κ²+(321/32)κ⁴+(κ⁴/64)Σ₆W+R, with '
        '‖R‖≤exp(36|κ|)(36|κ|)⁶/720. The induced Wilson coefficient is certified positive '
        'at the declared small steps. Electric matching requires g_f²=g_c²/4 at half spacing in '
        'the unchanged clock. This magnetic-factor blocking retains the full remainder and does '
        'not derive a full-Hamiltonian RG or four-dimensional continuum.', '',
        '### Finite tensor identities, a static classical register and scalar closure','',
        'The [compatible staggered tensor architecture](../model/openwave_staggered_ward/PROOF.md) '
        'removes the scalar checkerboard while retaining exact Bianchi/gauge identities and every halo. '
        'A co-located nine-null instrument and physical curvature/stress map remain open. The '
        '[replicated quantum geometry](../model/openwave_replicated_geometry/PROOF.md) has explicit '
        'O(1/N) static ground-energy errors and vanishing lapse noise, with its closing geometry gap '
        'retained. The failed G/c strength is unchanged; this is not an Einstein spacetime limit.', '',
        'The [capacity perturbation test](../model/openwave_capacity_perturbations/PROOF.md) rejects '
        'the adiabatic barotropic attachment. An established holographic-quintessence reconstruction '
        'uses the same supplied background and has sound speed squared '
        f"{value['capacity_scalar']['canonical_scalar_completion']['rest_frame_sound_speed_squared_exact']} "
        'and Lean-certified Vʺ/H² bounds '
        f"{interval(value['capacity_scalar']['canonical_scalar_completion']['V_second_derivative_over_H_squared_bounds_exact'])}. "
        'Its native-source scalar transfer, quantum ultraviolet completion and cosmological perturbation '
        'likelihood remain open. These source advances add zero Nature-validation rows.', '']
    return lines


def complete_dynamics_and_constraints(data):
    """Project complete source advances without promoting physical identification."""
    names=('full_strong_current','lambda_wavepacket','cluster_weyl_sequence',
           'refined_haar_source','refined_haar_spectrum','haar_memory_closure',
           'tensor_record_transport','collective_gravity_strength','capacity_transfer')
    flags=(
        ('full_strong_current',('physical_3plus1_decay_or_lifetime_proved',)),
        ('full_strong_current',('outgoing_Moller_wave_operator_or_on_shell_weak_spectral_density_proved',)),
        ('lambda_wavepacket',('new_independent_match_with_nature',)),
        ('lambda_wavepacket',('physical_3plus1_Lambda_identification',)),
        ('lambda_wavepacket',('calibrated_energy','mean_energy_is_rest_mass')),
        ('cluster_weyl_sequence',('outgoing_wave_operator_constructed',)),
        ('cluster_weyl_sequence',('nonzero_native_weak_on_shell_spectral_density_proved',)),
        ('refined_haar_source',('uniform_repeated_refinement_or_4D_limit_proved',)),
        ('refined_haar_spectrum',('physical_GeV_calibration_or_new_Nature_validation',)),
        ('refined_haar_spectrum',('derived_running_coupling_or_continuum_limit',)),
        ('haar_memory_closure',('retained_201_modes_exactly_autonomous',)),
        ('haar_memory_closure',('uniform_refinement_or_physical_4D_limit_proved',)),
        ('tensor_record_transport',('scope','nonlinear_or_continuum_Einstein_equation_proved')),
        ('tensor_record_transport',('scope','quantum_curvature_acquisition_from_native_source_proved')),
        ('capacity_transfer',('scope','new_Nature_validation')),
        ('capacity_transfer',('scope','CMB_large_scale_structure_likelihood')),
    )
    for name,path in flags:
        row=data[name]
        for key in path:row=row[key]
        if row:raise ValueError('unproved complete-dynamics physical promotion: '+name)
    gravity=data['collective_gravity_strength']
    if (gravity['calibrated_member']['electron_anchor_changed']
            or 'observed Newton' not in gravity['calibrated_member']['calibration_input']
            or gravity['all_N_theorem']['huge_carrier_enumerated_or_simulated']):
        raise ValueError('collective gravity calibration or carrier boundary omitted')
    packet=data['lambda_wavepacket'];weyl=data['cluster_weyl_sequence']
    if (not packet['normalizable_state']['actual_normalizable_Lambda_trial_constructed']
            or not weyl['constructive_Weyl_family']['full_infinite_chain_vectors_are_normalizable']
            or weyl['constructive_Weyl_family']['finite_local_observer_preparation_of_exact_spectral_vectors_compiled']
            or not data['refined_haar_spectrum']['fine_interior_modes_are_retained']
            or not data['haar_memory_closure']['actual_computed_self_energy']['all_Q_representations_and_multiplicities_in_remainder']):
        raise ValueError('complete source domain or preparation boundary omitted')
    return {'physical_comparison_count_added':0,
            'source_packages':{name:{k:v for k,v in data[name].items()if k!='source_sha256'}for name in names}}


def render_complete_dynamics_and_constraints(value):
    data=value['source_packages'];packet=data['lambda_wavepacket']
    energy=packet['calibrated_energy'];spectrum=data['refined_haar_spectrum']
    memory=data['haar_memory_closure'];gravity=data['collective_gravity_strength']['calibrated_member']
    lines=['### Normalizable Lambda records and the interacting breakup spectrum','',
        'The [normalizable Lambda trial](../model/openwave_lambda_wavepacket/PROOF.md) '
        'uses the infinite-chain source and a complete local current experiment without an exact band filter. '
        'Its ground-band probability exceeds '+decimal(packet['complete_source_pole']['unfiltered_trial_ground_band_probability_lower'])+
        '. These are source spectral weights, not percentage agreements with Nature.','',
        '| Quantity | Calibrated source interval / GeV |','| --- | ---: |',
        '| Lambda rest pole | '+interval(energy['source_rest_pole_GeV_interval'])+' |',
        '| Lambda packet mean, including kinetic energy | '+interval(energy['trial_packet_mean_GeV_interval'])+' |','',
        'The electron, strong and Lambda calibration inputs retain their ancestry; these values are conditional, '
        'not independent mass predictions. Optional exact pole conditioning is a nonlocal service.', '',
        'The [complete autonomous five-body current](../model/openwave_full_strong_current/PROOF.md) '
        'has a certified positive separated-field response with every strong term and probe failure retained. '
        'The [normalizable cluster Weyl construction](../model/openwave_cluster_weyl_sequence/PROOF.md) '
        'places the Lambda-assigned source energy in the full interacting proton–pion essential spectrum. '
        'Its exact spectral-vector coefficients define mathematical source states; their finite observer '
        'preparation is not supplied. The outgoing wave operator, on-shell weak spectral weight and physical '
        'decay lifetime are not established.', '',
        '### Interior gauge modes and computed memory','',
        'The [54-link, 36-face Haar source](../model/openwave_refined_haar_source/PROOF.md) '
        'has a complete 201-mode electric frame and retains every representation and intertwiner outside it. '
        'The [same-clock spectrum comparison](../model/openwave_refined_haar_spectrum/PROOF.md) gives:', '',
        '| Spatial source | Spacing | g² | Full source gap |','| --- | ---: | ---: | ---: |']
    for row in spectrum['records']:
        lines.append('| '+str(row['source'])+' | '+row['spacing_exact']+' | '+row['g_squared_exact']+
                     ' | '+interval(row['full_gap_interval_exact'])+' |')
    correction=memory['actual_computed_self_energy']
    lines += ['', 'The [computed self-energy](../model/openwave_haar_memory_closure/PROOF.md) has '
        +str(correction['nonzero_entries'])+' nonzero matrix entries. Its all-complement operator error is at most '
        +decimal(correction['full_Sigma0_minus_computed_matrix_norm_upper_exact'],upper=True)+
        ' source units. The 225-fold improvement at t=1/1024 concerns the specified Haar-vacuum defect bound. '
        'No exact autonomous compression, physical GeV glueball or four-dimensional continuum follows.', '',
        '### Finite tensor transport and conditional gravity strength','',
        'The [four-corner classical record map](../model/openwave_tensor_record_transport/PROOF.md) '
        'joins the staggered table to the existing nine-null decoder with the explicit error '
        'δ=(2+√3)ε+τ and adjacent scalar bound aw+8δ. Lean checks the transport inequality. '
        'This map does not acquire a quantum stress or curvature field.', '',
        'A [fixed-frequency collective geometry](../model/openwave_collective_gravity_strength/PROOF.md) '
        'uses the observed Newton coefficient to select approximately '
        +f"{gravity['diagnostic_float']['registers']:.6g}"+' registers. It admits conditional G≈'
        +f"{gravity['diagnostic_float']['conditional_G_GeV_minus2']:.7g}"+' GeV⁻² with the unchanged electron anchor. '
        'This is an explicit calibration, with complete readout work about '
        +f"{gravity['diagnostic_float']['readout_work_GeV']:.6g}"+' GeV. The huge carrier is evaluated analytically, '
        'not simulated; no shared strong/geometry continuum or independent G prediction is claimed.', '',
        '### Constrained scalar, matter and radiation transfer','',
        'The [linear Einstein/fluid/scalar calculation](../model/openwave_capacity_transfer/PROOF.md) '
        'proves C_N=−(3+2h)C and the spatial Einstein trace exactly, while evolving the metric and both fluids. '
        'For the inherited background and regular growing initial data:', '',
        '| k/H₀ | Present matter transfer / same-density Λ control |','| ---: | ---: |']
    for row in data['capacity_transfer']['numerical_transfer_comparisons']:
        lines.append(f"| {row['capacity']['k_over_H0']:g} | {row['present_comoving_matter_transfer_ratio_to_lambda']:.8f} |")
    lines += ['', 'These are numerical model differences with tolerance, initialization and Einstein–de Sitter '
        'controls, not observed matches or rigorous interval ODE certificates. Native scalar and primordial '
        'state, free-streaming species, recombination and a CMB likelihood are not supplied. '
        'These results add zero Nature-validation rows.', '']
    return lines
