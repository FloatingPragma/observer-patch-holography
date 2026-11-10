"""Same-source charged limits and declared collective field constructions."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from .canvas import AMBER, BLUE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
PATHS = {
    'charge': 'model/coherent_defect_wave/chiral_joint/interacting_charge_localization/results.json',
    'propagation': 'model/interacting_block_propagation/results.json',
    'field': 'model/interacting_block_quadrature/results.json',
    'collision': 'model/interacting_block_scattering/results.json',
    'charged_field': 'model/affine_gauss_full_spectrum/results.json',
    'charged_clock': 'model/affine_gauss_charged_clock/results.json',
    'observed_charge': 'model/affine_gauss_clock_observables/results.json',
    'native_energy': 'model/native_abelian_energy_records/results.json',
    'dirac': 'model/native_block_dirac/results.json',
    'higgs_join': 'model/native_higgs_spectral_join/results.json',
    'resolution': 'model/native_charge_resolution/results.json',
    'centroid': 'model/native_charge_centroid/results.json',
    'higgs_state': 'model/native_higgs_filtered_state/results.json',
    'weak_response': 'model/native_weak_current_response/results.json',
    'higgs_current': 'model/native_higgs_current_response/results.json',
    'coupled_field': 'model/native_wilson_maxwell/results.json',
    'photon_cutoff': 'model/native_wilson_maxwell_refinement/results.json',
    'quark_response': 'model/native_quark_response/results.json',
    'quark_higgs': 'model/native_quark_higgs_join/results.json',
    'coupled_ground': 'model/native_wilson_maxwell_ground/results.json',
    'ground_refinement': 'model/native_wilson_maxwell_ground_refinement/results.json',
    'local_current': 'model/native_local_current_readout/results.json',
    'vector': 'model/interacting_block_proca/results.json',
    'conjugation': 'model/coherent_defect_wave/chiral_joint/interacting_charge_conjugation/results.json',
}


def validate(data):
    field = data['field']
    if field['construction']['measured_inputs'] or not field['construction']['full_physical_complement_retained']:
        raise ValueError('Field construction inputs changed')
    theorem = field['theorems']
    if (not theorem['massive_free_KG_limit_in_increasing_capacity_subsector']
            or not theorem['three_components_are_internal_scalar_labels']
            or theorem['nontrivial_interacting_relativistic_limit_proved']
            or theorem['Lorentz_vector_or_Proca_identification_proved']):
        raise ValueError('Field interpretation changed')
    for row in data['charge']['full_source_charged_bins'].values():
        law = row['energy_filtered_localization']
        if (not law['every_normalized_state_in_bin_and_every_time']
                or law['finite_bin_complement_or_higher_states_truncated']
                or law['origin_probability_lower_inferred_by_discarding_screened_complement']):
            raise ValueError('Charged localization domain changed')
    collision = data['collision']['interval_certified_five_block_collision']
    if (collision['actual_interacting_ground_input_claim']
            or collision['three_dimensional_continuum_claim']
            or collision['exact_reducing_source_dimension'] != 1024):
        raise ValueError('Finite collision interpretation changed')
    charged = data['charged_field']
    if (charged['additional_state_truncation']
            or charged['physical_one_charge_dimension'] != 70928
            or charged['scope']['mass_GeV'] is not None
            or charged['scope']['Standard_Model_identity'] is not None
            or charged['scope']['translation_invariant_mass_shell_identified']):
        raise ValueError('Charged field interpretation changed')
    clock = data['charged_clock']
    if (not clock['phase_apparatus']['common_source_evolution_retained']
            or clock['phase_apparatus']['higher_source_charge_sectors_discarded']
            or clock['scope']['exact_ground_projection']
            or clock['scope']['GeV_mass'] is not None):
        raise ValueError('Charged detector domain changed')
    vector = data['vector']
    if (vector['construction']['Standard_Model_particle_assignment'] is not None
            or vector['construction']['measured_inputs']
            or vector['construction']['strict_continuum_interaction'] != 'free; fixed-capacity interactions are controlled errors'
            or not vector['complete_instruments']['all_failure_and_complement_records_retained']):
        raise ValueError('Vector readout interpretation changed')
    if data['conjugation']['reference_boundary']['same_reference_total_energy_equality_claimed']:
        raise ValueError('Conjugate reference energy was silently changed')
    observed=data['observed_charge']
    if (observed['exact_ground_projection_or_stationary_post_clock_state_claimed']
            or observed['positive_joint_record_lower_bounds']!=26):
        raise ValueError('Prepared position record changed')
    if data['native_energy']['physical_spin_Dirac_or_mass_identification']:
        raise ValueError('Native energy promoted to mass')
    spinor=data['dirac']['field_comparison']
    if (spinor['interacting_continuum_or_electromagnetic_Gauss_law_claimed']
            or spinor['SM_weak_chirality_or_observed_mass_identity_claimed']
            or spinor['mass_parameter_is_empirically_calibrated']
            or not spinor['same_source_state_reference_and_clock']):
        raise ValueError('Native Dirac scope changed')
    selector=data['higgs_join']['same_H_cross_residue_and_records']['positive_cross_selector']
    if (not selector['same_energy_elementary_composite_G1_doublet']
            or selector['maximum_elementary_E_star_common_composite_overlap_proved']
            or selector['physical_SM_weak_chirality_interactions_or_GeV_identity_proved']):
        raise ValueError('Higgs spectral join promoted')
    if (data['resolution']['original_maximum_residue_eigenvalue_relocalized']
            or data['resolution']['physical_spin_mass_or_SM_identity_claimed']):
        raise ValueError('Nested energy selection promoted')
    if (data['centroid']['actual_physical_preparation']['extra_P0_ground_projection_used']
            or data['centroid']['relativistic_mass_spin_or_SM_identification']
            or data['centroid']['clock']['exact_native_Eg_used_as_free_pointer_phase']):
        raise ValueError('Centroid preparation or energy origin changed')
    state=data['higgs_state']
    if (state['physical_preparation']['exact_interacting_ground_projection_assumed']
            or not state['complete_instruments']['every_fine_preparation_record_and_partner_retained']
            or state['typed_free_Dirac_corollary']['interacting_native_Lorentz_or_SM_chirality_transfer']):
        raise ValueError('Filtered-state domain changed')
    weak=data['weak_response']
    if (weak['actual_preparation']['exact_P0_added']
            or not weak['source']['same_source_clock_on_all_branches']
            or weak['scope']['weak_scattering_cross_section_or_SM_matching']
            or weak['scope']['W_boson_or_Lorentz_vector_identified']
            or Fraction(weak['finite_inelastic_experiment']['original_attempt_nonzero_phase_probability_lower'])<=0):
        raise ValueError('Weak response or physical interpretation changed')

    hc=data['higgs_current']
    if (hc['conjugation']['W_invariant']
            or hc['preparation']['exact_interacting_ground_projection_assumed']
            or hc['scope']['weak_boson_chiral_Lorentz_interaction_or_SM_identification']
            or not hc['apparatus']['full_fine_preparation_records_retained']
            or hc['apparatus']['complete_setting_pump_energy_current_alphabet_size']!=80):
        raise ValueError('Higgs current or complete detector changed')
    wm=data['coupled_field']
    if (wm['boundaries']['interacting_Lorentz_or_four_dimensional_QED_limit_proved']
            or wm['boundaries']['SM_particle_or_weak_chirality_identification']
            or wm['complete_interaction_records']['source_state_is_actual_coupled_ground']
            or not wm['actual_preparation']['no_initial_P0_or_code_measurement']
            or wm['count_instrument']['full_pointer_dimension']!=1161
            or Fraction(wm['complete_interaction_records']['exact']['actual_unconditional_target_probability_lower'])<=0):
        raise ValueError('Common field source or physical interpretation changed')
    cut=data['photon_cutoff']
    if (any(cut['boundaries'].values())
            or not cut['construction']['same_native_trial_state_reference_density_Gauss_and_clock']
            or cut['cap3_pair_photon_transfer']['finite_executed_pointer_labels']!=1161
            or cut['cap3_pair_photon_transfer']['common_overflow_padded_labels']!=1451
            or cut['cap3_pair_photon_transfer']['padding_changes_parent_executed_apparatus']
            or Fraction(cut['cap3_complete_state_and_record_bound']['complete_record_TV_upper'])<=0):
        raise ValueError('Photon cutoff source, records or interpretation changed')
    q=data['quark_response']
    if (not q['scope']['all_local_probe_roots_and_exact_zero_residues_retained']
            or any(v for k,v in q['scope'].items()if k!='all_local_probe_roots_and_exact_zero_residues_retained')
            or q['same_source']['additional_exact_ground_projection']
            or not q['computed_equal_coupling_uc_dc_resolvent_equality']):
        raise ValueError('Quark spectral source or interpretation changed')
    ground=data['coupled_ground'];completion=ground['declared_onsite_completion']
    if (ground['unchanged_source']['ground_in_code']
            or ground['unchanged_source']['active_ground_multiplicity']!=96
            or completion['global_scalar_energy_subtraction']!='none'
            or completion['old_source_receipt_or_H_mutated']
            or completion['actual_ground_radiation_or_SM_mass_claim']
            or not completion['complete_pattern_preserving_count_laws_exactly_invariant']
            or ground['complete_repaired_ground_preparation']['ground_energy_or_exact_ground_projector_used']
            or ground['exact_macrocell_and_Gauss_completion_control']['counterterm_is_forced_by_primitive_repair_law']):
        raise ValueError('Coupled ground, explicit completion or preparation changed')
    local=data['local_current']
    if (local['source']['full_Hext_cancellation_used_for_final_readout']
            or not local['source']['source_hamiltonian_continues_through_final_pulse']
            or local['control_resources']['full_matrix_preparation_and_exact_spectral_bin_premises_removed']
            or not local['control_resources']['coupling_is_declared_not_compiled_from_native_C_E']
            or local['complete_prepared_records']['complete_setting_pump_energy_current_alphabet_size']!=80
            or local['scope']['Standard_Model_species_or_GeV_mass_identified']):
        raise ValueError('Local current control or interpretation changed')
    qh=data['quark_higgs']
    if (qh['same_source']['additional_P0']
            or not qh['same_source']['all_complements_and_failure_records_retained']
            or qh['scope']['canonical_composite_CAR']
            or qh['scope']['physical_Lorentz_chirality_or_Standard_Model_matching']
            or qh['scope']['source_energies_are_calibrated_rest_masses']
            or qh['full_static_reference_controls']['full_joint_CAR_modes']!=60
            or any(Fraction(r['original_attempt_joint_contrast_lower'])<=0
                   for r in qh['full_interacting_interference_records'].values())):
        raise ValueError('Quark Higgs operator or complete record interpretation changed')
    qs=qh['normalized_admitted_states_and_histories']
    if (qs['additional_P0'] or qs['unconditional_complete_laws_claimed_close']
            or any(r['same_history_policy_may_read_pump_setting']
                   or r['original_pump_success_rates_equated']
                   or not r['bound_applies_to_source_and_common_retained_preparation_environment']
                   for r in qs['rows'].values())):
        raise ValueError('Quark admitted-state comparison domain changed')
    gr=data['ground_refinement'];gc=gr['fixed_complete_energy_clock']
    if (gr['same_source_family']['spatial_volume_or_spacing_refinement']
            or gr['boundaries']['ground_radiation_or_SM_mass_identification']
            or gr['boundaries']['uniform_spatial_or_relativistic_cone']
            or gc['imaginary_heat_time_is_physical_clock']
            or not gc['same_pointer_alphabet_for_all_K_and_mathematical_comparison']
            or not gc['all_M_basis_records_and_all_noop_failure_branches_retained']
            or gr['constructive_complete_preparation_history']['matrices_enumerated_or_resource_feasibility_claimed']):
        raise ValueError('Ground refinement clock, records or scope changed')


def load():
    data, pins = {}, {}
    for name, relative in PATHS.items():
        path = ROOT / relative
        data[name] = json.loads(path.read_text())
        local = data[name]['input_sha256']
        for key, expected in local.items():
            target = (ROOT / key).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != expected:
                raise ValueError('Field evidence changed: ' + key)
        pins.update(local)
        owned=dict(data[name].get('local_sha256', {}))
        owned.update(data[name].get('owned_source_sha256', {}))
        for key, expected in owned.items():
            target=(path.parent/key).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest()!=expected:
                raise ValueError('Field implementation changed: '+key)
            pins[target.relative_to(ROOT).as_posix()]=expected
        pins[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    validate(data)
    data['pins'] = pins
    return data


def draw(c, data, y):
    d = data['field_records']
    h = 415
    card(c, X0, y, CW, h, stroke=VIOLET)
    c.text(X0+28, y+43, 'A charged energy record is not automatically a moving particle', 27, TEXT, 700)
    c.para(X0+28, y+84, CW-56,
           'In the same eight-cell source, a static reference keeps the compensating charge. '
           'Moving one unscreened charge away requires excited links. A low-energy detector '
           'therefore sees a strong restriction on that motion.', 19, SOFT, 1.4)
    for i, (key, label) in enumerate((('ec', 'eᶜ removal'), ('L', 'L removal'))):
        bound = Fraction(d['charge']['full_source_charged_bins'][key]['energy_filtered_localization']
                         ['union_remote_single_unscreened_charge_probability_upper'])
        xx = X0+28+i*650
        c.text(xx, y+190, label, 20, VIOLET, 700)
        c.text(xx, y+226, f'remote single-charge probability ≤ {float(bound):.3e}', 22, TEAL, 600)
    c.para(X0+28, y+266, CW-56,
           'This bound holds after low-bin selection at every source time. The complete detector '
           'also keeps the origin and a multi-charge cloud; screening can live in that cloud. '
           'Positive higher-energy arrivals are a different setting. These facts do not prove '
           'that a free electron has been identified.', 19, SOFT, 1.4)
    c.link(X0+28, y+h-20, 'Exact Gauss, spectral-separation and screening proof',
           '../../'+PATHS['charge'].replace('results.json', 'PROOF.md'), 16)
    y += h+35
    y = charged_field_panel(c, d, y)+35
    y = native_energy_panel(c, d, y)+35
    y = higgs_current_panel(c, d, y)+35
    y = local_current_panel(c, d, y)+35
    y = quark_response_panel(c, d, y)+35
    y = dirac_panel(c, d, y)+35
    y = coupled_field_panel(c, d, y)+35
    y = photon_cutoff_panel(c, d, y)+35
    y = coupled_ground_panel(c, d, y)+35
    h = 640
    card(c, X0, y, CW, h, stroke=TEAL)
    c.text(X0+28, y+43, 'One computed neutral gap → a controlled relativistic field description', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'Take complete copies of the interacting eight-cell source. Its native current connects '
           'the ground to three neutral excitations. Declare paired exchange between those copies '
           'and retain every higher source state. The exact finite dynamics has a controlled '
           'large-capacity comparison.', 19, SOFT, 1.4)
    for j in range(3):
        xx = X0+95+j*190
        c.rect(xx, y+198, 125, 108, '#132636', TEAL, 1.5, 10)
        for a in range(2):
            for b in range(2):
                c.circle(xx+38+48*a, y+230+42*b, 10, [TEAL, BLUE, VIOLET, AMBER][2*a+b])
        if j<2:
            c.arrow(xx+132, y+251, xx+181, y+251, TEAL, 2, 7)
    c.text(X0+332, y+339, 'N full source blocks per spatial cell', 19, SOFT, 500, 'middle')
    c.arrow(X0+685, y+252, X0+775, y+252, VIOLET, 3, 9)
    c.text(X0+812, y+220, 'ω² → Δ² + (Δ²/6)|p|²', 29, TEAL, 700)
    c.para(X0+812, y+258, CW-850,
           'A massive Klein–Gordon description. Δ is the actual native gap, '
           'not a number chosen from a particle table.', 19, SOFT, 1.4)
    c.para(X0+28, y+391, CW-56,
           'Ground-state, finite-capacity, evolution and detector errors are bounded together. '
           'A retained finite orientation reference makes the field measurements relational. '
           'The unfiltered components are internal scalar fields. The next panel shows the '
           'additional readout and symmetry map that give a free vector description.', 19, SOFT, 1.4)
    c.para(X0+28, y+494, CW-56,
           'The paired exchange, replication and spacing-dependent clock are explicit design '
           'choices. The established refinement becomes free. A separate number-conserving '
           'exchange gives a Schrödinger limit. Neither construction supplies a Standard Model '
           'identity or a GeV calibration.', 19, MUTED, 1.4)
    c.link(X0+28, y+h-20, 'Complete collective field and reference proof',
           '../../'+PATHS['field'].replace('results.json', 'PROOF.md'), 16)
    c.link(X0+780, y+h-20, 'Distinct number-conserving propagation proof',
           '../../'+PATHS['propagation'].replace('results.json', 'PROOF.md'), 16)
    y += h+35
    y = vector_panel(c, d, y)+35
    h = 500
    card(c, X0, y, CW, h, stroke=AMBER)
    c.text(X0+28, y+43, 'Finite interactions: prepare two packets and read a complete detector', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'Capacity limits couple the neutral components. A large-capacity theorem identifies '
           'a nonzero quartic channel and bounds packet admission, preparation and higher-order '
           'errors below its amplitude. The selected hit and its entire complement are retained.', 19, SOFT, 1.4)
    frames = d['collision']['interval_certified_five_block_collision']['complete_time_series']
    for j, frame in enumerate(frames):
        xx = X0+100+j*250
        lo, hi = map(Fraction, frame['native_preparation_hit_probability_interval_exact'])
        c.text(xx, y+210, 'Δt = '+frame['clock_Delta_times_t'], 18, MUTED, 500)
        c.text(xx, y+246, f'{100*float(lo):.5f}% –', 21, AMBER, 700)
        c.text(xx, y+278, f'{100*float(hi):.5f}%', 21, AMBER, 700)
    c.para(X0+28, y+321, CW-56,
           'Shown: an independent five-block, 1,024-state diagnostic, including preparation '
           'error. Its finite projection changes the free comparison’s momentum supports, '
           'so this plot alone does not isolate a free-forbidden channel. Its reference-relative occupation '
           'records are quantum probabilities, not classical trajectories.', 19, SOFT, 1.4)
    contrast = Fraction(frames[-1]['native_minus_matched_free_hit_difference_interval_exact'][0])
    c.text(X0+28, y+449,
           f'Same prepared state, same detector: final source − free probability ≥ {100*float(contrast):.5f} percentage points.',
           19, TEAL, 600)
    c.link(X0+28, y+h-20, 'Finite collision theorem, comparison domains and exact intervals',
           '../../'+PATHS['collision'].replace('results.json', 'PROOF.md'), 16)
    return y+h


def charged_field_panel(c, d, y):
    h = 800
    card(c, X0, y, CW, h, stroke=BLUE)
    c.text(X0+28, y+43, 'A charge and its dynamical field: an actual finite energy', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'A declared electric/magnetic source lets charge move while its compensating reference '
           'stays fixed. Gauss constraints determine the longitudinal flux. Finite Peierls hopping '
           'couples that charge to every retained transverse field state.', 19, SOFT, 1.4)
    labels = [('13 matter modes', 'where the charge can be'),
              ('5,456 field states', 'the field can change with it'),
              ('70,928 states', 'the full one-charge calculation')]
    for i, (label, sub) in enumerate(labels):
        xx=X0+28+i*460
        c.rect(xx, y+186, 414, 86, '#132636', BLUE, 1.5, 10)
        c.text(xx+20, y+221, label, 23, BLUE, 700)
        c.text(xx+20, y+250, sub, 17, SOFT)
        if i<2:c.arrow(xx+420,y+230,xx+451,y+230,BLUE,2,7)
    lo, hi = map(Fraction, d['charged_field']['actual_charged_ground_energy_interval'])
    scale = 10**13
    lo_n, hi_n = (lo*scale).__floor__(), (hi*scale).__ceil__()
    lo_text = f'{lo_n//scale}.{lo_n%scale:013d}'
    hi_text = f'{hi_n//scale}.{hi_n%scale:013d}'
    c.text(X0+28, y+320, f'Total charged energy: [{lo_text}, {hi_text}] source units', 24, TEAL, 700)
    c.para(X0+28, y+351, CW-56,
           'The fixed reference contributes exactly 1. The field dressing lowers the zero-field '
           'compression by about 0.0001750058. This correction is computed with the same Hamiltonian; '
           'it is unrelated to a fitted lepton formula.', 18, SOFT, 1.4)
    c.text(X0+28, y+440, 'Where would a position detector click after the finite energy clock accepts?', 20, TEXT, 600)
    for i,row in enumerate(d['observed_charge']['post_clock_position_marginals']):
        lo,hi=map(Fraction,row['conditional_post_clock_probability_interval'])
        p=float((lo+hi)/2);xx=X0+40+i*104; bh=310*p
        color=AMBER if row['source_vertex']==1 else BLUE
        c.rect(xx,y+570-bh,62,bh,color,r=3)
        c.text(xx+31,y+595,str(row['source_vertex']),16,MUTED,500,'middle')
        c.text(xx+31,y+556-bh,f'{100*p:.1f}%',15,SOFT,500,'middle')
    c.text(X0+28,y+626,'All 26 position/field outcomes remain. Field excitation: 70.9795–73.5248 per million accepted attempts.',17,MUTED)
    clock=d['charged_clock']['finite_rational_gate_record_bounds']
    fidelity_n=(Fraction(clock['accepted_source_ground_fidelity_lower'])*10**11).__floor__()
    fidelity_text=f'{fidelity_n//10**11}.{fidelity_n%10**11:011d}'
    c.para(X0+28,y+658,CW-56,
           f'A finite clock detector retains all 14 × 2²⁰ records, including 75% preparation failure. '
           f'Its accepted source has ground fidelity above {fidelity_text}. '
           'The full energy and numerical errors are bounded. This finite boundary problem '
           'does not identify a relativistic particle, its spin or a GeV mass.',19,SOFT,1.4)
    c.link(X0+28,y+h-20,'Full charged Hamiltonian and independent residual certificate',
           '../../'+PATHS['charged_field'].replace('results.json','PROOF.md'),16)
    c.link(X0+810,y+h-20,'Complete finite clock and joint records',
           '../../'+PATHS['observed_charge'].replace('results.json','PROOF.md'),16)
    return y+h


def native_energy_panel(c,d,y):
    h=1160
    card(c,X0,y,CW,h,stroke=TEAL)
    c.text(X0+28,y+43,'Two probes can see the same native excitation',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'One probe adds an elementary charged fermion. Another combines the quantum Higgs '
           'operator with a lepton operator. The actual interacting equation of motion links them. '
           'Their positive spectral overlap proves that both can reach the same source energy.',19,SOFT,1.4)
    for i,(title,body,color) in enumerate([
            ('Elementary probe A','native eᶜ creator',BLUE),
            ('Same actual energy','positive common G₁ doublet',TEAL),
            ('Composite probe B','quantum Higgs × lepton',VIOLET)]):
        xx=X0+28+i*460
        c.rect(xx,y+188,414,94,'#132636',color,1.5,10)
        c.text(xx+20,y+224,title,22,color,700)
        c.text(xx+20,y+260,body,17,SOFT)
    c.arrow(X0+449,y+235,X0+479,y+235,TEAL,2,7)
    c.arrow(X0+942,y+235,X0+911,y+235,TEAL,2,7)
    c.text(X0+28,y+331,'[h, A] = B + λ[T, A]     ·     complete interference contrast > 0.07247007595',23,TEAL,700)
    c.para(X0+28,y+369,CW-56,
           'The scalar remains quantum. The experiment retains both interference settings, '
           'all spin/reference outcomes and failed preparations. A shared finite energy does '
           'not turn the composite into a canonical field on every occupation or establish weak scattering.',19,SOFT,1.4)
    c.text(X0+28,y+468,'Same prepared state: every common detector history differs by < 1.686 × 10⁻⁷',23,VIOLET,700)
    c.para(X0+28,y+507,CW-56,
           'A narrower energy window makes the elementary and composite preparations almost '
           'indistinguishable after success. Eight retained source trials and seven finite filters '
           'bound the preparation error. Their success rates differ and are kept in the records. '
           'The exact energy-bin operation is a declared apparatus service.',19,SOFT,1.4)
    c.text(X0+28,y+623,'A physical 20-bit clock: mean source energy 0.6916566792649 – 0.6916566826119',23,TEAL,700)
    c.para(X0+28,y+663,CW-56,
           'This predicted conditional average includes actual preparation and pointer errors, '
           'with no extra ground projection. It retains every clock label and costs the reference '
           'explicitly. A mean energy is not an identified mass or an individual spectral line.',19,SOFT,1.4)
    c.text(X0+28,y+768,'A dressed weak current responds on this same native source',23,VIOLET,700)
    c.para(X0+28,y+808,CW-56,
           'The complete current pump and finite energy-agreement detector give a nonzero-phase '
           'record with original probability above 3.79338%. A charged energy-filtered preparation '
           'also has a positive current insertion. Every failure is counted. This is not a '
           'W-boson line or a weak-scattering cross section.',19,SOFT,1.4)
    c.para(X0+28,y+929,CW-56,
           'On a distinct native charged-field attachment, position followed by the actual finite '
           'energy clock has original arrival-and-energy-acceptance probability above 4.50 × 10⁻⁶ '
           'at the apex and 1.13 × 10⁻⁶ at a boundary neighbor. Measurement backaction is included. '
           'This is an admitted-bin energy, not an identified relativistic rest mass.',19,SOFT,1.4)
    c.link(X0+28,y+h-52,'Quantum Higgs equation and complete interference detector',
           '../../'+PATHS['higgs_join'].replace('results.json','PROOF.md'),16)
    c.link(X0+810,y+h-52,'Same-H arrival then energy records',
           '../../'+PATHS['native_energy'].replace('results.json','PROOF.md'),16)
    c.link(X0+28,y+h-20,'Normalized state and complete-history bound',
           '../../'+PATHS['higgs_state'].replace('results.json','PROOF.md'),16)
    c.link(X0+810,y+h-20,'Finite prepared centroid and full phase clock',
           '../../'+PATHS['centroid'].replace('results.json','PROOF.md'),16)
    c.link(X0+810,y+h-84,'Actual weak-current moments and full inelastic detector',
           '../../'+PATHS['weak_response'].replace('results.json','PROOF.md'),16)
    return y+h


def dirac_panel(c,d,y):
    h=790
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'Native charged states → a complete free spin-½ field',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Four complete native constituents at each spatial site encode two spin states for '
           'each conjugate charge. Every occupation is present. A declared positive spatial '
           'coupling preserves the actual ground and moves these modes in the same source clock.',19,SOFT,1.4)
    for word in range(16):
        xx=X0+40+(word%8)*165; yy=y+191+(word//8)*61
        c.rect(xx,yy,139,48,'#132636',VIOLET,1,8)
        for bit in range(4):
            occupied=(word>>bit)&1
            c.circle(xx+25+bit*29,yy+24,7,VIOLET if occupied else '#253047',SOFT if occupied else MUTED,1)
    c.text(X0+28,y+340,'4 fermionic modes · all 16 occupations · native exchange signs · all complements retained',20,SOFT,600)
    c.text(X0+28,y+395,'Ωₐ² = E*² + (4/a²) Σ sin²(apᵢ/2)    →    ω² = E*² + |p|²',28,TEAL,700)
    c.para(X0+28,y+432,CW-56,
           'The rest frequency E* is an actual source-selected excitation, with no measured mass '
           'input. Spinor filters combine particle annihilation and antiparticle creation. '
           'Their opposite phases yield the free Dirac equation, spin-½ Lorentz action and causal '
           'anticommutator in the controlled record limit. Both source charges have positive energy.',19,SOFT,1.4)
    c.text(X0+28,y+560,'A two-time detector: second-only click 5% → 13.82355% at quarter-period pulse separation',22,VIOLET,700)
    c.para(X0+28,y+600,CW-56,
           'All 16 detector labels and preparation failures are retained; six complete attempts '
           'bound preparation-related record error below 10⁻⁹. The kinetic couplers have full '
           'spatial range and the apparatus uses declared programmable control. This vectorlike '
           'eᶜ code supplies no electron weak chirality, interacting charged limit or GeV identity.',19,MUTED,1.4)
    c.link(X0+28,y+h-20,'Complete native Fock code, same-clock spinor filters, Lorentz limit and detector errors',
           '../../'+PATHS['dirac'].replace('results.json','PROOF.md'),16)
    return y+h


def vector_panel(c,d,y):
    h=610
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'From three canonical fields to a massive vector readout',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Keep the same native blocks, actual ground, clock and orientation reference. '
           'Measure specific combinations of the fields and their conjugate momenta. '
           'Those filters supply three vector polarizations and the constrained time component.',19,SOFT,1.4)
    for i,(title,body,color) in enumerate([
            ('Canonical pairs','Φ₁, Φ₂, Φ₃ and Π₁, Π₂, Π₃',TEAL),
            ('Momentum-dependent filters','A₀ uses i p·Π / Ω',BLUE),
            ('Complete vector detector','16 records, including every flag',VIOLET)]):
        xx=X0+28+i*460
        c.rect(xx,y+194,414,104,'#132636',color,1.5,10)
        c.text(xx+20,y+231,title,21,color,700)
        c.text(xx+20,y+270,body,17,SOFT)
        if i<2:c.arrow(xx+420,y+248,xx+451,y+248,color,2,7)
    c.para(X0+28,y+342,CW-56,
           'The polarization sum, Proca equations and induced spin-one Lorentz action are proved '
           'for the free limiting description. The finite source’s Lorenz defect and complete '
           'measurement errors are bounded. The time component needs canonical momentum; '
           'three field labels alone do not establish this result.',19,SOFT,1.4)
    c.text(X0+28,y+454,'Rest frequency = the native gap Δ. Effective mass parameter = Δ / c*².',23,TEAL,700)
    c.para(X0+28,y+494,CW-56,
           'The explicit finite examples bound sixteen-outcome laws using a 19-mode Gaussian '
           'calculation and full native error estimates. Their large block counts are analytic '
           'capacity bounds. This readout does not identify W, Z or an interacting Standard Model sector.',19,MUTED,1.4)
    c.link(X0+28,y+h-20,'Same-source vector filters, boosts, full records and capacity bounds',
           '../../'+PATHS['vector'].replace('results.json','PROOF.md'),16)
    return y+h


def higgs_current_panel(c,d,y):
    h=600
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'Ask the same excitation the same current question',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Prepare the elementary excitation or its quantum-Higgs composite, admit the same '
           'energy window, then read the same dressed current. The detector has five signed '
           'answers. Negative current is a detector value, not negative energy.',19,SOFT,1.4)
    for i,(name,prob,col) in enumerate([
            ('Elementary A','at least 2.488191% of original attempts',BLUE),
            ('Higgs composite B/2','at least 0.297229% of original attempts',VIOLET)]):
        xx=X0+28+i*690
        c.rect(xx,y+185,650,150,'#132636',col,1.5,10)
        c.text(xx+22,y+224,name,23,col,700)
        c.text(xx+22,y+267,'Admitted current = −½',24,TEAL,700)
        c.text(xx+22,y+307,prob,19,SOFT)
    c.para(X0+28,y+373,CW-56,
           'Both conditional current means are strictly negative. Their successful readouts '
           'differ by less than 4 × 10⁻⁶ in total variation in this wider certified bin. '
           'Their different original success rates remain. Choosing a setting with probability '
           'one half halves its joint rate.',19,SOFT,1.4)
    c.para(X0+28,y+478,CW-56,
           'Conjugating the source also transforms the current operator: copying the opposite-charge '
           'value is incorrect. Four source intervals, all 80 apparatus labels and every preparation '
           'failure are retained. W names this current, not an identified W boson.',19,MUTED,1.4)
    c.link(X0+28,y+h-20,'Exact conjugation, current blocks, full-source bounds and complete apparatus',
           '../../'+PATHS['higgs_current'].replace('results.json','PROOF.md'),16)
    return y+h


def coupled_field_panel(c,d,y):
    h=850
    card(c,X0,y,CW,h,stroke=TEAL)
    c.text(X0+28,y+43,'A common finite spinor and quantum field: the interaction is calculated',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Eight spatial vertices contain 32 complete native constituents. Their fermion codes '
           'couple to five transverse modes. The same charge density enters the current, electric '
           'energy and Gauss constraints, including the boundary reference.',19,SOFT,1.4)
    steps=[('Prepare','deterministic native trial',BLUE),('Evolve','one coupled source clock',VIOLET),('Record','charge + conjugate + field quantum',TEAL)]
    for i,(title,desc,col) in enumerate(steps):
        xx=X0+28+i*460
        c.rect(xx,y+185,414,130,'#132636',col,1.5,10)
        c.text(xx+20,y+225,title,24,col,700)
        c.para(xx+20,y+253,374,desc,18,SOFT,1.35)
        if i<2:c.arrow(xx+420,y+250,xx+450,y+250,col,2,7)
    response=d['coupled_field']['complete_interaction_records']['exact']
    def lower(key):
        x=Fraction(response[key]);scale=10**13;return f'{(x*scale).__floor__()/scale:.5e}'
    c.text(X0+28,y+374,'Pair + field probability ≥ '+lower('actual_unconditional_target_probability_lower'),26,TEAL,700)
    c.text(X0+28,y+421,'Connected matter–field response ≥ '+lower('actual_unconditional_binary_connected_response_lower'),24,VIOLET,700)
    c.para(X0+28,y+461,CW-56,
           'These original-attempt bounds include the full nonlinear evolution remainder and the '
           'actual preparation weight. Two independent algebraic current calculations agree. '
           'With coupling zero the target record has probability zero. All 1,161 outcomes '
           'remain, including source complements and Gauss failure.',19,SOFT,1.4)
    c.text(X0+28,y+588,'g = 1/64 · mT = 1/64 · m = the shared native energy · displayed capacity < 10,479 nats',21,TEAL,700)
    c.para(X0+28,y+633,CW-56,
           'The input is a free-state preparation, not the coupled ground. This proves a finite '
           'interaction response, not vacuum radiation, an asymptotic scattering cross section '
           'or a GeV mass. Spacing a = 1/m, coupling, regulator, finite photon cap and full-matrix '
           'controls are declared choices; control synthesis workspace is counted separately.',19,SOFT,1.4)
    c.para(X0+28,y+744,CW-56,
           'A physical particle interpretation requires the coupled spectrum and a controlled '
           'spatial refinement with its local currents. The free Lorentz readout alone does '
           'not supply these interaction properties.',18,MUTED,1.4)
    c.link(X0+28,y+h-20,'Same-source Gauss, current, preparation and complete interaction proof',
           '../../'+PATHS['coupled_field'].replace('results.json','PROOF.md'),16)
    return y+h


def photon_cutoff_panel(c,d,y):
    h=565
    card(c,X0,y,CW,h,stroke=TEAL)
    c.text(X0+28,y+43,'More field memory: a controlled comparison on the same cube',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Matter, native preparation, gauge reference and clock stay fixed. Only the total '
           'photon cap grows. Each finite program has an error bound against five unbounded '
           'comparison oscillators; no infinite register is physically executed.',19,SOFT,1.4)
    rows=[r for r in d['photon_cutoff']['cutoff_table']if r['cap']in(1,2,3,4,8)]
    for i,row in enumerate(rows):
        xx=X0+28+i*275
        cap=row['cap'];value=Fraction(row['state_vector_error_upper'])
        # Fixed decimal exponent and integer ceiling keep every printed upper bound outward.
        exponent=-1
        while value<Fraction(10)**exponent:exponent-=1
        scale=Fraction(10)**(3-exponent)
        shown=f'{-(-(value*scale).numerator//(value*scale).denominator)/1000:.3f}e{exponent}'
        c.rect(xx,y+200,255,118,'#15383d'if cap==3 else '#132636',TEAL if cap==3 else BLUE,1.5,10)
        c.text(xx+16,y+235,f'cap {cap}'+(' · current'if cap==3 else ''),19,SOFT,700)
        c.text(xx+16,y+283,'error ≤ '+shown,20,TEAL,700)
    c.para(X0+28,y+357,CW-56,
           'At cap three, complete-record error is below 1.863e-9; the positive pair-plus-field '
           'response survives. The comparison keeps photon overflow: 1,161 executed labels '
           'are padded to 1,451 labels. Bounded electric/magnetic test histories also have '
           'complete-output bounds, with the actual finite commutator correction retained.',19,SOFT,1.4)
    c.para(X0+28,y+468,CW-56,
           'This closes an occupation-cutoff question on a fixed graph. Spatial refinement, '
           'a uniform local propagation cone and an interacting Lorentz limit remain distinct.',18,MUTED,1.4)
    c.link(X0+28,y+h-20,'Fixed-source photon capacity, field algebra and complete histories',
           '../../'+PATHS['photon_cutoff'].replace('results.json','PROOF.md'),16)
    return y+h


def decimal_bound(value, places=9, upper=False):
    """Print an outward rational bound without a binary64 rounding step."""
    v=Fraction(value);scale=10**places
    n=-(-(v*scale).numerator//(v*scale).denominator) if upper else (v*scale).__floor__()
    sign='−' if n<0 else '';a=abs(n)
    return f'{sign}{a//scale}.{a%scale:0{places}d}'


def quark_response_panel(c,d,y):
    h=830
    card(c,X0,y,CW,h,stroke=AMBER)
    c.text(X0+28,y+43,'Three color-current families: every local spectral line is retained',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Use the actual quark operators and their compensating quantum reference. Each family '
           'has a complete energy-response measure. A strong line means that this apparatus '
           'often prepares that energy; counting lines does not count particle species.',19,SOFT,1.4)
    for i,(key,label,nroots) in enumerate([('uc','uᶜ current',13),('dc','dᶜ current',13),('Q','Q doublet current',25)]):
        xx=X0+28+i*460
        rows=d['quark_response']['full_interacting_prepared_records'][key]['records']
        row=max(rows,key=lambda r:Fraction(r['original_attempt_probability_lower']))
        lo,hi=row['source_excitation_relative_actual_ground']
        c.rect(xx,y+200,420,185,'#2c251e',AMBER,1.5,10)
        c.text(xx+20,y+237,label,24,AMBER,700)
        c.text(xx+20,y+276,f'{nroots} local roots · strongest source bin',18,SOFT)
        c.text(xx+20,y+314,'['+decimal_bound(lo)+', '+decimal_bound(hi,upper=True)+']',20,TEXT,600)
        prob=decimal_bound(100*Fraction(row['original_attempt_probability_lower']),6)
        c.text(xx+20,y+358,'original attempts ≥ '+prob+'%',19,TEAL,700)
    c.para(X0+28,y+426,CW-56,
           'Equal chosen local Yukawa couplings give exactly equal uᶜ and dᶜ local measures. '
           'Q has 13 visible roots and 12 exact local zeros; the interacting detector still '
           'retains every bin, leakage, pump failure and preparation failure. Its static '
           'reference costs 4 source-energy units and has no spatial hopping.',19,SOFT,1.4)
    c.text(X0+28,y+556,'Actual Higgs dynamics: [h, A] = C + the bounded hopping term',23,VIOLET,700)
    c.para(X0+28,y+595,CW-56,
           'The elementary quark probe A and its quantum scalar/matter composite C share '
           'a positive contribution in the same bin. Complete interference contrasts exceed '
           '0.096796971 for uᶜ/dᶜ and 0.096907802 for Q, including actual preparation failures. '
           'Their successful states give common later histories within 5.193e-6; different '
           'initial setting labels and success rates remain. The scalar stays quantum.',19,SOFT,1.4)
    c.para(X0+28,y+699,CW-56,
           'These are source-energy windows with controlled native hopping errors. Their '
           'spin and gauge operator algebra is checked, but moving quark poles, observed '
           'up/down masses, generations and GeV conversion are not identified by this spectrum.',18,MUTED,1.4)
    c.link(X0+28,y+h-20,'Complete quark spectral records',
           '../../'+PATHS['quark_response'].replace('results.json','PROOF.md'),16)
    c.link(X0+650,y+h-20,'Actual quantum-Higgs operators and interference',
           '../../'+PATHS['quark_higgs'].replace('results.json','PROOF.md'),16)
    return y+h


def coupled_ground_panel(c,d,y):
    h=1030
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'Which state is really lowest? Audit the entire source',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'An interacting code can have a positive gap while another retained sector lies '
           'below it. The ground audit includes all native complements, Gauss constraints '
           'and reference costs before choosing a preparation.',19,SOFT,1.4)
    for i,(title,body,col) in enumerate([
            ('Original field source','96 active lowest states lie outside the intended code. A code-only gap cannot make its vacuum the ground.',AMBER),
            ('Explicit completed source','Add a positive onsite cost to each code complement. The active coupled ground is now unique; code energies stay unchanged.',TEAL)]):
        xx=X0+28+i*690
        c.rect(xx,y+185,650,165,'#132636',col,1.5,10)
        c.text(xx+22,y+224,title,23,col,700)
        c.para(xx+22,y+254,606,body,19,SOFT,1.4)
    c.text(X0+28,y+403,'H completed = H original + m c(g) Σₓ (I − Pₓ)',27,VIOLET,700)
    prep=d['coupled_ground']['complete_repaired_ground_preparation']
    success=decimal_bound(100*Fraction(prep['complete_tree_success_lower']),4)
    c.text(X0+28,y+452,'4 complete source trials → 3 finite filters → success ≥ '+success+'%',22,TEAL,700)
    c.para(X0+28,y+490,CW-56,
           'Conditional ground-state error is below 1.418e-7. All failed filters, partner sources '
           'and pointer labels remain. The complete physical gap is above 2.5207e-10 source units; '
           'it differs from the much larger gap inside the code.',19,SOFT,1.4)
    c.para(X0+28,y+590,CW-56,
           'The ground detector can record virtual matter–field occupation in its declared free '
           'basis. This is stationary dressing, not radiation from a vacuum. The penalty is an '
           'explicit architecture choice, not a primitive-repair consequence. Existing global '
           'kinetic guards remain; spatial locality and particle masses still need their proofs.',19,MUTED,1.4)
    c.text(X0+28,y+732,'The same interacting ground family survives increasing photon capacity',23,TEAL,700)
    for i,(title,body) in enumerate([
            ('K = 1024 · ideal ground','Two complete field-test histories: error < 3.094e-9.'),
            ('K = 2⁹⁶ · actual preparation','All preparation/count or field histories: error < 9.549e-9.')]):
        xx=X0+28+i*690
        c.rect(xx,y+766,650,100,'#132636',TEAL,1.2,10)
        c.text(xx+20,y+800,title,21,TEAL,700)
        c.text(xx+20,y+840,body,18,SOFT)
    c.para(X0+28,y+905,CW-56,
           'One fixed finite filter alphabet retains every failure at every cap. These huge '
           'matrices are not enumerated, and control/work costs grow with capacity. This is '
           'a fixed-cube field-occupation limit, not a spatial or relativistic continuum.',18,MUTED,1.4)
    c.link(X0+28,y+h-20,'Full-ground audit and actual preparation',
           '../../'+PATHS['coupled_ground'].replace('results.json','PROOF.md'),16)
    c.link(X0+650,y+h-20,'Uniform ground and complete-record capacity limit',
           '../../'+PATHS['ground_refinement'].replace('results.json','PROOF.md'),16)
    return y+h


def local_current_panel(c,d,y):
    h=550
    card(c,X0,y,CW,h,stroke=BLUE)
    c.text(X0+28,y+43,'A local apparatus reads a current while the source keeps running',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Couple the current at one carrier to a neutral five-state pointer for a short '
           'pulse. Different current values shift the pointer to different records. The '
           'source Hamiltonian continues throughout; the detector does not cancel it.',19,SOFT,1.4)
    for i,(title,body,col) in enumerate([
            ('One native carrier','all 80 current states',BLUE),
            ('Local coupling','f(W) × pointer generator',VIOLET),
            ('Five records','every current and failure kept',TEAL)]):
        xx=X0+28+i*460
        c.rect(xx,y+188,414,110,'#132636',col,1.5,10)
        c.text(xx+20,y+225,title,23,col,700)
        c.text(xx+20,y+266,body,18,SOFT)
        if i<2:c.arrow(xx+420,y+245,xx+450,y+245,col,2,7)
    c.para(X0+28,y+340,CW-56,
           'Complete-record error < 0.000268446640 at clock 3 + 2⁻²⁰. Both negative current '
           'means survive, including the quark terms in the drift bound. The source and '
           'apparatus preserve every Gauss constraint during the pulse.',19,SOFT,1.4)
    c.para(X0+28,y+438,CW-56,
           'This polynomial coupling is a declared local control. Its strong, short pulse '
           'has an explicit work bound. Earlier state preparation and spectral-bin services '
           'remain separate controls; efficient compilation from primitive repairs is open.',18,MUTED,1.4)
    c.link(X0+28,y+h-20,'Exact local pointer, uninterrupted source clock and complete-history error',
           '../../'+PATHS['local_current'].replace('results.json','PROOF.md'),16)
    return y+h
