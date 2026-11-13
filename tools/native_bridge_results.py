"""Source-qualified field/record calculations and the actual M10 gate."""
from copy import deepcopy

NAMES=('lambda_narrow_weight','repeated_native_clock','quark_string_transport',
       'nuclear_field_attraction','matter_flux_source','local_gauge_memory',
       'stress_curvature_instruments','geometry_covariance_protocol',
       'capacity_finite_wavelength','m10_review_gate')

def project(data):
    selected={n:{k:deepcopy(v)for k,v in data[n].items()if k!='source_sha256'}for n in NAMES}
    for n,key in [('lambda_narrow_weight','physical_decay_rate_or_lifetime_proved'),
                  ('quark_string_transport','physical_quark_pole_or_GeV_rest_mass_proved')]:
        if selected[n][key] is not False:raise ValueError('unsupported physical promotion: '+n)
    counterterm=selected['matter_flux_source']['optional_charged_counterterm']
    if counterterm['physical_mass_prediction_or_native_process_closure'] or counterterm['unchanged_current_Ward_work_or_many_body_decay_transferred']:
        raise ValueError('optional rest counterterm promoted beyond its reducing sectors')
    clock=selected['repeated_native_clock']['scope']
    if clock['static_controller_exact_strength_independent_frequency'] or clock['sustained_replenished_clock'] or clock['literal_openwave_row3_validated']:
        raise ValueError('finite clock promoted beyond its complete controller')
    cosmology=selected['capacity_finite_wavelength']['scope']
    if cosmology['late_preparation_is_early_primordial_solution'] or cosmology['absolute_amplitude_predicted']:
        raise ValueError('late cosmological response promoted beyond its preparation')
    gate=selected['m10_review_gate']['official_status']
    if not gate['oph_admitted'] or gate['oph_model_id']!='M10' or gate['upstream_counts']['M10']!={'validated':0,'partial':0,'negative':0,'untested':31,'scored':0}:
        raise ValueError('actual pinned M10 status differs')
    return {'source_packages':selected,'physical_comparison_count_added':0,
        'native_twelve_port_implementation_established':False,
        'physical_premises_closed':[],
        'official_openwave_column':'M10','official_openwave_validated_rows':0,
        'official_openwave_untested_rows':31}

def render(value,interval,decimal):
    d=value['source_packages'];cosmo=d['capacity_finite_wavelength']['fine']
    lines=['### Complete source records and the native-process attachment','',
        'OPH is admitted to OpenWave as M10, with 31 untested physics rows. The '
        '[upstream arithmetic certification](../model/openwave_m10_review_gate/PROOF.md) '
        'reproduces seven exact dispersion quantities and awards no criterion cell. '
        'A criterion claim requires its observable to execute on the native twelve-port '
        'repair/memory process in-platform. Declared source controls retain their mathematical domains.','',
        'The [complete Lambda residual witness](../model/openwave_lambda_narrow_weight/PROOF.md) '
        'uses 32 full-source trial columns and exact translated-packet inner products. Its '
        'current spectral weight within 1 MeV of the declared Lambda centre is below 0.026622. '
        'The certified lower bound is zero; no shrinking-window density or lifetime follows.','',
        'The [repeated clock experiment](../model/openwave_repeated_native_clock/PROOF.md) '
        'retains two readouts, reset and finite phase/battery resources. Its complete '
        'static controller has a resolved frequency bias, distinct from the ideal '
        'isolated source frequency. The '
        '[autonomous quark-string source](../model/openwave_quark_string_transport/PROOF.md) '
        'has a proved reducing Jacobi chain inside its full Gauss carrier, with every '
        'endpoint, failure and work record retained. Neither calculation establishes '
        'a self-starting physical clock or a physical quark rest mass.','',
        'A [neutral dynamical field](../model/openwave_nuclear_field_attraction/PROOF.md) '
        'generates a signed four-site attraction by completing the oscillator square. '
        'It vanishes below 12 quarks and preserves the triton/electron input ancestry. '
        'The supplied coupling fails the held-out helion and alpha masses. The exact '
        'full-singlet Casimir trace replaces enormous identity-state enumeration for '
        'electromagnetic preparation work, independently checked on 2520 character weights.','',
        'The [full 192-mode matter and U(1) cycle](../model/openwave_matter_flux_source/PROOF.md) '
        'has a stable self-adjoint Hamiltonian with all integer fluxes. Exact Gauss '
        'dressing transports references and currents. Full-tail spectral bounds find '
        'a positive common charged-lepton energy shift; calibrated rest values cannot '
        'be transferred unchanged to this spatial extension. An optional source-defined '
        'subtraction restores all inherited charged one-lepton rest coefficients '
        'with one common number and preserves the global vacuum. It is a separately '
        'declared generator: physical masses, currents and many-body decays do not '
        'follow from this rest-sector identity.','',
        'The [local full-Haar memory bound](../model/openwave_local_gauge_memory/PROOF.md) '
        'uses no global electric cutoff. Fundamental and spin-one Wilson correlations '
        'exceed 0.97 on three refinements at the specified common time. Its error '
        'grows with refinement, so it supplies no fixed-time continuum limit. '
        'Complete temporal instruments retain backaction and uniform readout-work bounds.','',
        'The [stress/curvature instrument](../model/openwave_stress_curvature_instruments/PROOF.md) '
        'reduces an entire 34-site tensor stencil to ten same-state quantum means. '
        'All 60 source coordinates and 30 outcome probabilities are retained; sampling '
        'and backaction errors enter the existing sharp null-inverse and staggered '
        'transport theorems. The measured Einstein residual is nonzero.','',
        'The [geometry covariance protocol](../model/openwave_geometry_covariance_protocol/PROOF.md) '
        'retains complete paired measurements, finite-sample errors and their work. '
        'Independent register-count/stiffness information and a common rest/photon '
        'unit attachment are required for an unfitted Newton constant.','',
        'For the [specified late cosmological preparation](../model/openwave_capacity_finite_wavelength/PROOF.md), '
        'the k² metric coefficient lies in '+str(cosmo['metric_k_squared_coefficient_display'])+
        ' and the comoving dust coefficient in '+str(cosmo['dust_k_squared_coefficient_display'])+
        '. A complete k^4 remainder certifies two nonzero wavelengths. These responses '
        'are conditional on the supplied Einstein/scalar/fluid equations and their '
        'normalization; they are not a primordial-history or CMB fit.','',
        'These calculations close their stated finite source obligations. They add '
        'zero independent Nature-validation rows and discharge no physical '
        'identification by combining incompatible sources.','']
    return lines
