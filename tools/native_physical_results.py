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
        'all_sources_identified_as_one_physical_theory':False,
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
        'varies with scale, while the independent spatial core is unresolved below one cell. '
        'These sources are not pooled into one identified physical theory.', '']
    return lines
