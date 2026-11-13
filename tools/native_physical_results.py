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
