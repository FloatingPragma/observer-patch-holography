"""Receipt-bound finite defect-wave interaction and its complete readouts."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as s

from .canvas import AMBER, BLUE, EDGE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / 'model/coherent_defect_wave/results.json'
EXTENSIONS = ('background', 'spectroscopy', 'spatial', 'species',
              'cubic_readout', 'energy', 'dispersion', 'physical_attachment', 'mobile_flux',
              'native_routing_join', 'many_body_route', 'mass_instrument',
              'field_identity_audit', 'propagating_native_join', 'chiral_source',
              'chiral_regulator', 'chiral_higgs', 'chiral_particle_map', 'chiral_joint',
              'chiral_joint/two_cell', 'chiral_joint/periodic_spectrum',
              'chiral_joint/periodic_family', 'chiral_joint/field_limit',
              'chiral_joint/higgs_composites', 'chiral_joint/spectral_edges',
              'chiral_joint/oriented_links', 'chiral_joint/neutral_mobility',
              'chiral_joint/family_cell', 'chiral_joint/neutral_local_spectrum',
              'chiral_joint/neutral_composite', 'chiral_joint/periodic_boson_records',
              'chiral_joint/peter_weyl_links', 'chiral_joint/family_transport',
              'chiral_joint/local_reference', 'chiral_joint/matter_induced_ring',
              'chiral_joint/matter_induced_plaquette', 'chiral_joint/finite_clock',
              'chiral_joint/atomic_reference_control', 'chiral_joint/atomic_even_remainder',
              'chiral_joint/charged_ground_exchange', 'chiral_joint/atomic_spin_geometry',
              'chiral_joint/atomic_momentum_response', 'chiral_joint/ground_composite_response',
              'chiral_joint/ground_scalar_response', 'chiral_joint/collective_spin_response',
              'chiral_joint/interacting_ground_preparation', 'chiral_joint/interacting_charged_response')
FORBIDDEN = {'physical_mass_identified', 'relativistic_mass_or_GeV_value',
             'GeV_mass', 'physical_particle_identity', 'SM_species_identified',
             'rest_masses_derived', 'interacting_relativistic_continuum',
             'GeV_prediction', 'physical_rest_mass_identified',
             'energies_are_physical_rest_masses', 'source_derived_GeV_masses',
             'physical_particle_or_Standard_Model_identification',
             'SM_particle_identification', 'GeV_masses', 'interacting_Lorentz_limit',
             'dynamical_chiral_gauge_theory', 'SM_identity_or_GeV_mass_derived',
             'interacting_relativistic_limit_proved', 'SM_particle_or_GeV_mass_identified',
             'positive_energy_many_body_vacuum_derived', 'physical_interacting_relativistic_limit',
             'scattering_or_SM_identity_or_GeV_mass', 'positive_energy_many_body_vacuum',
             'physical_statistics_or_chiral_matter', 'interacting_relativistic_limit_or_SM_or_GeV',
             'chiral_SM_particle_identification', 'interacting_relativistic_limit',
             'source_GeV_masses', 'physical_Standard_Model_identification',
             'interacting_relativistic_mass', 'GeV_mass_derived',
             'source_mass_renormalization_or_SM_identity_or_GeV_prediction',
             'interacting_relativistic_particle_identification',
             'background_is_a_quantum_Higgs_vacuum_theorem',
             'interacting_chiral_continuum', 'measured_species_identified',
             'onsite_gap_is_relativistic_particle_mass', 'full_reference_is_physical_Lorentz_vacuum',
             'physical_pole_mass_or_GeV_prediction', 'physical_rest_mass_or_GeV_identification',
             'observed_species_identified', 'finite_gap_is_Clay_mass_gap',
             'GeV_masses_or_observed_species_identified', 'interacting_continuum_limit_proved',
             'source_only_particle_momentum_or_relativistic_dispersion_proved',
             'named_particle_or_GeV_mass', 'interacting_relativistic_field_limit',
             'GeV_mass_prediction', 'interacting_Lorentz_continuum',
             'physical_spin_or_species_identified', 'physical_interacting_continuum',
             'observed_particles_identified', 'continuum_em_charge_or_GeV_mass_identified',
             'physical_particle_or_mass_identified', 'relativistic_limit_or_mass_claim',
             'source_operator_identified_with_measured_particle',
             'interacting_relativistic_refinement_proved', 'measured_particle_or_GeV_mass_identified',
             'source_energy_is_particle_rest_mass', 'spatial_propagation_or_GeV_mass',
             'moving_neutrino_or_relativistic_mass_identified', 'three_observed_neutrino_generations',
             'empirical_neutrino_identification', 'named_physical_bosons_identified',
             'Higgs_photon_W_Z_gluon_particle_identification', 'physical_Lorentz_spin_or_mass',
             'continuum_propagation_or_GeV_mass', 'physical_mass_or_relativistic_propagation_inferred',
             'representation_names_are_measured_particle_identifications',
             'prepared_states_identified_as_ground_or_free_particles',
             'physical_photon_gluon_W_or_Z_identified', 'physical_photon_or_W_Z_gluon_identified',
             'Coulomb_confining_or_relativistic_phase', 'physical_GeV_scale',
             'interacting_many_cell_vacuum_proved',
             'raw_charge_identified_with_measured_particle_number'}


def check_extension_scope(value, name):
    """Check nested claims too: a declared rest parameter is permitted, an SM mass is not."""
    if isinstance(value, dict):
        for key, item in value.items():
            if key=='measured_GeV_predictions' and item:
                raise ValueError('joined extension physical promotion: '+name+'/'+key)
            if key=='empirically_identified_species' and item is not None:
                raise ValueError('joined extension physical promotion: '+name+'/'+key)
            absent_mass_table=(key=='source_GeV_masses' and isinstance(item,(dict,list)) and not item)
            if key in FORBIDDEN and item is not False and not absent_mass_table:
                raise ValueError('joined extension physical promotion: ' + name + '/' + key)
            check_extension_scope(item, name + '/' + key)
    elif isinstance(value, list):
        for item in value:
            check_extension_scope(item, name)


def check_native_join(value):
    if not value['shared_tuple']['native_U_reads_the_same_link_register_as_routing']:
        raise ValueError('native routing join has no shared link register')
    finite=value['finite_join']
    if finite['Taylor_output_renormalized'] or not finite['complete_binary_instrument']:
        raise ValueError('native routing join substitutes an incomplete approximate readout')
    rows=finite['exact_unitary_probability_enclosures']
    if len(rows)!=2 or {row['changed_face'] for row in rows}!={False,True}:
        raise ValueError('native routing join omits a detector outcome')
    bounds={row['changed_face']:(Fraction(row['lower']),Fraction(row['upper'])) for row in rows}
    if (not all(0<=lo<=hi<=1 for lo,hi in bounds.values()) or bounds[True][0]<=0
            or bounds[False][0]+bounds[True][1]!=1 or bounds[False][1]+bounds[True][0]!=1):
        raise ValueError('native routing join probability bounds are not complete')



def check_mass_instrument(value, free):
    if value['source'] != free['program']:
        raise ValueError('mass apparatus and free field have different source clocks')
    if (not free['reference']['unique_positive_excitation_ground']
            or not free['scope']['CAR_statistics_are_declared_architecture']
            or not free['source']['all_occupations_retained']
            or free['algebra']['Wilson_clock_is_Floquet_logarithm']
            or value['continuum_rest_energy_in_source_units']!=free['reference']['full_Fock_gap']):
        raise ValueError('mass apparatus lacks its declared positive free reference')
    if value['unit_map']['source_energy_unit_in_GeV'] is not None:
        raise ValueError('mass apparatus substitutes an unowned GeV scale')
    for experiment in value['experiments']:
        records=experiment['complete_records']
        weights=[s.sympify(r['probability']) for r in records]
        if (len(records)!=36 or s.simplify(sum(weights,s.Integer(0)))!=1
                or any(p.is_nonnegative is not True for p in weights)):
            raise ValueError('mass apparatus omits complete reservoir records')
        if not experiment['filled_negative_band_reference']:
            raise ValueError('mass apparatus loses its positive reference')


def check_identity(value):
    dark=value['native_dark_current_module']
    if (not value['scope']['actual_native_R15_current_binding']
            or dark['module_prepared_or_bound_by_native_idle_law']
            or dark['native_K_on_module']!='0' or dark['internal_dimension']!=15
            or dark['empty_native_flux_code_actual_current_representation']!='trivial'):
        raise ValueError('native current module or empty flux identity changed')
    mass=value['mass_normalization_and_chirality']
    if (mass['gauge_invariant_Dirac_mass_parameter_dimension']!=5
            or not mass['Dirac_R15_contains_both_chiralities_in_same_internal_representation']):
        raise ValueError('charged Dirac identity loses its mass or chirality boundary')


def check_two_cell(value):
    if (value['schema']!='oph.coherent_defect_wave.chiral_joint.two_cell.v1'
            or value['source']['active_dimension']!=614400):
        raise ValueError('two-cell physical source changed')
    ground=next(row for row in value['full_spectra'] if row['boundary']==0)
    if (ground['full_physical_dimension']!=1053 or not ground['ground_is_certifiably_simple']
            or ground['ground_window_candidate_labels']!=[{'L':2,'twoSz':0,'eigen_index':0}]):
        raise ValueError('two-cell complete physical ground is missing')
    law=value['original_attempt_readout']
    rows=law['records']
    if (law['probabilities_renormalized'] or law['record_count_including_preparation_failure']!=len(rows)
            or len(rows)!=461 or rows[0]['outcome']!='preparation_failure'
            or rows[0]['probability_interval']!=['1052/1053']*2
            or law['preparation']['original_ground_success_probability']!='1/1053'):
        raise ValueError('two-cell original preparation failures are missing')
    if any(not 0<=s.Rational(r['probability_interval'][0])<=s.Rational(r['probability_interval'][1])<=1
           for r in rows):
        raise ValueError('two-cell invalid probability interval')
    low=sum(s.Rational(r['probability_interval'][0]) for r in rows)
    high=sum(s.Rational(r['probability_interval'][1]) for r in rows)
    if not low<=1<=high:raise ValueError('two-cell original probability law is incomplete')
    if {r['kind'] for r in value['conditional_readout']['records']}!={'add','add_fail','remove','remove_fail'}:
        raise ValueError('two-cell charged probe drops a branch')


def load():
    value = json.loads(RECEIPT.read_text())
    if value['schema'] != 'oph.coherent_defect_wave.v1':
        raise ValueError('joint defect-wave schema changed')
    if not value['input_sha256']:
        raise ValueError('joint source custody missing')
    for name, digest in value['input_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('stale joint defect-wave input: ' + name)
    for key in ('continuous_dynamical_gauge_theory', 'complete_boundary_gauss_constraints',
                'complete_joint_effect_tomography', 'interacting_3plus1_continuum',
                'physical_particle_species', 'physical_energy_or_GeV_mass'):
        if value['scope'][key]:
            raise ValueError('joint physical scope promotion: ' + key)
    if value['source']['dimension'] != 72000:
        raise ValueError('joint finite source domain changed')
    for key in ('wave_to_defect', 'defect_to_wave'):
        if not value['interaction'][key]['actual_full_tick_executed']:
            raise ValueError('joint interaction is not executed on the full tick')
    coherent = value['interaction']['constrained_coherence']
    if not (coherent['return_projector_commutes_with_Gamma'] and
            coherent['all_recorded_branches_preserve_positive_constraint']):
        raise ValueError('return apparatus leaves the admitted gauge sector')
    value['extensions'] = {}
    for name in EXTENSIONS:
        path = RECEIPT.parent / name / 'results.json'
        extension = json.loads(path.read_text())
        if not extension.get('input_sha256'):
            raise ValueError('missing joined extension custody: ' + name)
        for source, digest in extension['input_sha256'].items():
            target = (ROOT/source).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale joined extension input: ' + source)
        if (name in ('many_body_route','mass_instrument','field_identity_audit','propagating_native_join')
                and extension['schema']!='oph.coherent_defect_wave.'+name+'.v1'):
            raise ValueError('joined extension schema changed: '+name)
        check_extension_scope(extension, name)
        if name=='native_routing_join':check_native_join(extension)
        value['extensions'][name] = extension
    check_identity(value['extensions']['field_identity_audit'])
    check_two_cell(value['extensions']['chiral_joint/two_cell'])
    check_periodic(value['extensions']['chiral_joint/periodic_spectrum'],
                   value['extensions']['chiral_joint/higgs_composites'])
    check_neutral_motion(value['extensions']['chiral_joint/neutral_mobility'],
                         value['extensions']['chiral_joint/family_cell'])
    check_local_reference(value['extensions']['chiral_joint/neutral_local_spectrum'])
    check_neutral_detectors(value['extensions']['chiral_joint/neutral_composite'],
                            value['extensions']['chiral_joint/periodic_boson_records'])
    check_PW_transport(value['extensions'])
    check_interacting_reference(value['extensions'])
    check_charged_reference(value['extensions'])
    check_mass_instrument(value['extensions']['mass_instrument'],value['extensions']['many_body_route'])
    return value



def check_interacting_reference(extensions):
    ground=extensions['chiral_joint/atomic_spin_geometry']['full_interacting_transfer']
    records=extensions['chiral_joint/collective_spin_response']['same_source_interacting_records']
    line=records['actual_interacting_energy_line']
    if (not ground['unique_full_physical_interacting_ground_proved']
            or not ground['all_full_Fock_scalar_link_and_physical_Gauss_sectors_retained']
            or ground['unit_hopping_ground_statement_transferred']
            or ground['lowest_gap_is_a_zero_momentum_particle_rest_mass']
            or ground['named_species_or_GeV_masses']
            or ground['thermodynamic_or_relativistic_limit']
            or Fraction(ground['full_physical_ground_gap_interval_at_reported_coefficient'][0])<=0):
        raise ValueError('interacting ground domain or full complement changed')
    if (line['exact_full_source_cluster_rank']!=3
            or line['source_momenta_in_pi_units']!=[[1,0,0],[0,1,0],[0,0,1]]
            or not line['complete_outside_line_weight_strictly_positive']
            or Fraction(line['first_line_probability_conditioned_on_ground_preparation'][0])<=0
            or Fraction(line['outside_line_probability_conditioned_on_ground_preparation'][0])<=0
            or not line['ideal_finite_spectral_service_is_declared']):
        raise ValueError('complete collective line or declared apparatus changed')


def check_charged_reference(extensions):
    prep=extensions['chiral_joint/interacting_ground_preparation']['preparation']
    filt=prep['finite_energy_agreement_filter']
    if (not prep['scope']['constructive_finite_algebraic_service_schedule']
            or prep['scope']['physical_native_efficiency_or_mass_identification']
            or not 0<Fraction(filt['success_probability_lower'])<=1
            or not 0<Fraction(filt['conditional_one_source_trace_distance_upper'])<Fraction(1,1000000)):
        raise ValueError('physical preparation domain or error changed')
    charged=extensions['chiral_joint/interacting_charged_response']
    if charged['complete_probe_and_reference']['full_static_reference_modes']!=30:
        raise ValueError('charged reference occupations omitted')
    for row in charged['same_H_interacting_channels'].values():
        scope=row['scope']; law=row['full_interacting_energy_readout']
        if (scope['physical_rest_mass_or_free_charge_dispersion']
                or not scope['full_complement_error_control']
                or law['isolated_charged_eigenvalue_or_unchanged_multiplicity_claimed']
                or law['reference_energy_subtracted_silently']
                or Fraction(row['source']['reference_particle_cost'])!=4
                or not row['constructive_approximate_preparation_records']['all_preparation_success_and_failure_branches_obey_source_Gauss']
                or not row['reference_sector_comparison']['neutral_ground_tensor_reference_vacuum_is_unique_full_physical_ground_at_mu_four']):
            raise ValueError('charged source, reference or physical scope changed')
        for key in ('line','outside','failure'):
            lo,hi=map(Fraction,law['conditional_probabilities'][key])
            if not 0<lo<=hi<=1:
                raise ValueError('complete charged detector has no positive '+key)


def check_PW_transport(extensions):
    pw=extensions['chiral_joint/peter_weyl_links']['complete_neutral_transfer_experiment']
    family=extensions['chiral_joint/family_transport']
    dictionary=family['actual_native_family_dictionary']
    law=family['same_clock_full_H_probability_records']
    local=extensions['chiral_joint/local_reference']['global_physical_ground_and_gap']
    if (pw['extra_neutral_hopping_term'] or pw['discarded_source_states_or_dynamical_terms']
            or Fraction(pw['complete_arrival_probability_interval'][0])<=0):
        raise ValueError('Peter-Weyl source motion changed')
    if (dictionary['native_families']!=1
            or {r['six_Y'] for r in dictionary['native_irrep_rows']}!={-4,1,6,-3,2}
            or law['source_or_reference_states_discarded']
            or Fraction(law['strict_local_record_contrast_lower_bound'])<=0):
        raise ValueError('complete family transport changed')
    if (local['full_maximal_module_physical_dimension']!=104
            or local['global_physical_ground_rank']!=2 or local['dark_energy_one_rank']!=56
            or not local['full_physical_complement_excluded_by_strict_bound']
            or Fraction(local['global_physical_gap_interval'][0])<=0
            or local['interacting_many_cell_vacuum_proved']):
        raise ValueError('complete local reference or scope changed')


def check_periodic(value, higgs):
    neutral = next(row for row in value['complete_physical_spectra'] if row['boundary'] == [0, 0])
    law = value['compensated_charged_spectroscopy']
    if (neutral['physical_dimension'] != 3424 or not neutral['ground_certifiably_simple']
            or law['original_ground_preparation_success'] != '1/3424'
            or law['original_ground_preparation_failure'] != '3423/3424'
            or law['original_attempt_outcome_count_including_preparation_failure'] != 2627
            or law['probabilities_renormalized']
            or not law['joint_phase_energy_degeneracy_is_reference_position_structure']):
        raise ValueError('periodic ground or complete reference records changed')
    branches = law['charged_success_branches']
    if len(branches) != 2:
        raise ValueError('periodic conjugate branch omitted')
    first = []
    for branch in branches:
        if len(branch['complete_energy_windows']) != 328 or branch['exact_total_branch_probability'] != '1/4':
            raise ValueError('periodic energy alternatives omitted')
        row = branch['complete_energy_windows'][0]
        if Fraction(row['original_swap_probability_interval'][0]) <= 0:
            raise ValueError('periodic first line lacks certified weight')
        if row['joint_source_reference_phase_splitting_factors'] != ['0', '1/2', '0', '1/2']:
            raise ValueError('periodic phase apparatus changed')
        first.append(row['energy_transfer_interval'])
    if first[0] != first[1]:
        raise ValueError('periodic conjugate thresholds differ')
    bridge = higgs['same_periodic_H_spectral_bridge']
    if bridge['residues_are_instrument_probabilities'] or not bridge['quantum_Higgs_matrix_used']:
        raise ValueError('quantum Higgs residue interpretation changed')
    for branch in bridge['branches']:
        row = branch['first_residue_enclosure']
        if (not branch['exact_complete_Gauss_transfer_identity']
                or Fraction(row['quantum_composite_residue_interval'][0]) <= 0
                or row['residue_is_instrument_probability']):
            raise ValueError('quantum Higgs first-line connection lost')


def check_neutral_detectors(neutral, bosons):
    for value in (neutral, bosons):
        law = value['complete_readout']
        if (law['ground_preparation_success_probability'] != '1/3424'
                or law['retained_ground_preparation_failure_probability'] != '3423/3424'
                or not Fraction(law['complete_original_probability_enclosure'][0]) <= 1 <=
                Fraction(law['complete_original_probability_enclosure'][1])):
            raise ValueError('neutral detector loses complete reference records')
    if (neutral['complete_readout']['success_postselection_renormalization']
            or bosons['complete_readout']['success_probabilities_renormalized']
            or not bosons['complete_readout']['elastic_projectors_retained']):
        raise ValueError('neutral detector changes original response weights')
    rows = bosons['complete_readout']['settings']
    if ({r['setting'] for r in rows} != {'scalar_number', 'singlet_loop', 'doublet_loop'}
            or len({tuple(r['lowest_certified_positive_inelastic_record']['energy_transfer_interval'])
                    for r in rows}) != 1
            or any(Fraction(r['lowest_certified_positive_inelastic_record']
                            ['setting_success_energy_phase_probability_interval'][0]) <= 0 for r in rows)):
        raise ValueError('scalar and loop shared response is missing')


def check_local_reference(value):
    spectrum = value['complete_physical_spectrum']
    symmetry = value['physical_symmetries']
    records = spectrum['complete_energy_records']
    if (spectrum['full_seed_component_dimension'] != 113948
            or spectrum['full_physical_cyclic_dimension'] != 49
            or len(records) != 49
            or not spectrum['exact_full_graph_and_cyclic_intertwiners']
            or not spectrum['positive_Perron_overlap_of_physical_seed']
            or spectrum['global_ground_across_all_local_sectors_proved']
            or spectrum['reference_is_product_ground_or_particle_mass']
            or spectrum['full_component_excitation_gap_inferred_from_visible_line_spacing']):
        raise ValueError('local spectral reference scope changed')
    if (not all(0 < Fraction(r['probability_interval'][0]) <= Fraction(r['probability_interval'][1]) <= 1
                for r in records)
            or spectrum['outside_cyclic_sector_record'] != {'retained': True, 'probability_interval': ['0', '0']}
            or Fraction(spectrum['odd_minus_vacuum_component_ground_interval'][1]) >= 0):
        raise ValueError('local spectral reference comparison or complete weights changed')
    if (symmetry['centered_charge_L_minus_one'] != '0'
            or symmetry['raw_L_charge_conjugation'] != 'L -> 2-L'
            or not symmetry['charge_conjugation_pairs_ground_within_raw_L_one']
            or symmetry['raw_L_is_identified_with_measured_lepton_charge']
            or symmetry['spin_is_identified_with_a_relativistic_little_group']):
        raise ValueError('local reference charge or spin interpretation changed')


def check_neutral_motion(mobility, cell):
    law = mobility['complete_physical_arrival_experiment']
    if (not mobility['source']['all_base_interactions_retained']
            or law['old_ground_or_old_spectrum_consumed']
            or law['postselection_or_discarded_dynamical_sector']
            or law['complete_source_and_reference_readout'] != ['arrival', 'no_arrival']):
        raise ValueError('neutral motion changes its complete source or reference')
    lo, hi = map(Fraction, law['complete_arrival_probability_interval'])
    flo, fhi = map(Fraction, law['complete_no_arrival_probability_interval'])
    if not (0 < lo <= hi < 1 and lo + fhi == 1 and hi + flo == 1):
        raise ValueError('neutral motion loses its positive complete probability bound')
    spectrum = cell['complete_vacuum_energy_records']
    if (spectrum['full_physical_cyclic_dimension'] != 36
            or len(spectrum['complete_energy_records']) != 36
            or spectrum['numerical_eigensolver_used']
            or spectrum['global_other_sector_minimum_proved']
            or not spectrum['outside_cyclic_preparation_sector_record']['retained']):
        raise ValueError('full-family cyclic reference or complete records changed')


def input_pins(value):
    pins = {**value['input_sha256'], str(RECEIPT.relative_to(ROOT)):
            hashlib.sha256(RECEIPT.read_bytes()).hexdigest()}
    for name, extension in value['extensions'].items():
        path = RECEIPT.parent/name/'results.json'
        pins.update(extension['input_sha256'])
        pins[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def draw(c, data, y):
    value = data['joint']
    interaction, reference, gauss = value['interaction'], value['reference'], value['gauss']
    c.caps(X0, y + 10, 'A joint finite source: one protected region and one wave', TEAL, 14)
    c.text(X0, y + 62, 'Defects and waves change one another', 40, TEXT, 700, maxw=CW)
    y = c.para(X0, y + 105, CW,
        'The joint program retains the link memory in the wave’s quantum state. '
        'The links change the wave step; the wave’s port controls permitted defect moves. '
        'Both happen in one fixed reversible tick on ' + f"{value['source']['dimension']:,}" +
        ' basis states. No moved-or-held outcome is measured during idle propagation.',
        21, MUTED, 1.42) + 22

    h = 135
    card(c, X0, y, CW, h, stroke=TEAL)
    w = (CW - 120) / 3
    for i, (title, body, color) in enumerate((
            ('Retained link word g', 'Ten allowed words encode five defect positions.', AMBER),
            ('Wave step Wg', 'Coherent seam transport reads the actual link signs.', BLUE),
            ('Controlled link moves', 'The wave port selects local reversible repairs.', TEAL))):
        x = X0 + 26 + i * (w + 34)
        c.text(x, y + 43, title, 25, color, 700, maxw=w)
        c.para(x, y + 77, w - 14, body, 18, SOFT, 1.35)
        if i < 2:
            c.arrow(x + w - 5, y + 44, x + w + 20, y + 44, TEAL, 2.4, 9)
    y += h + 22

    h, gap, w = 364, 22, (CW - 44) / 3
    response = interaction['defect_to_wave']['yes_and_complement_probabilities']
    coherent = interaction['constrained_coherence']['coherent_return_and_complement'][0]
    recorded = interaction['constrained_coherence']['after_complete_face_readout_return_and_complement'][0]
    field = interaction['wave_to_defect']
    face_id = str(field['deterministic_final_face_indices'][0])
    face_probabilities = [row[face_id] for row in field['complete_gauge_invariant_face_distributions']]
    for i, (tag, title, body, number, color) in enumerate((
            ('Defect → wave', 'A different defect changes a click',
             'Keep the same two-port wave preparation. Changing the allowed link word '
             'changes the region’s defect position and the output of the same wave detector. '
             'Both detector complements are retained.',
             'Click probability: ' + response[0][0] + ' versus ' + response[1][0], BLUE),
            ('Wave → defect', 'A different wave changes the defect',
             'Prepare the inverse wave step of a port-5 or port-6 state. One complete tick '
             'then moves the defect to different faces. A detector of one regional face '
             'separates these two preparations.',
             'Face detector: ' + face_probabilities[0] + ' versus ' + face_probabilities[1], TEAL),
            ('Gauge-invariant interference', 'Recording a branch changes interference',
             'A wave superposition correlates with two defect positions. The apparatus uses '
             'the retained link frame to compare amplitudes. Undoing the controls restores '
             'the return. Reading defect position first halves it; both outcomes are retained.',
             'Return: ' + coherent + ' unread; ' + recorded + ' after readout', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 13)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 24, number, 17, color, 700, maxw=w - 48)
    y += h + 22

    h = 204
    card(c, X0, y, CW, h, stroke=BLUE)
    c.text(X0 + 28, y + 43, 'A reference and a gauge constraint for this same clock', 28, TEXT, 700)
    c.para(X0 + 28, y + 83, CW - 56,
        'The complete preparation filter uses this joint U, giving success probability ' +
        reference['success_probability_from_uniform_basis_ensemble'] +
        ' from its uniform basis ensemble. The admitted interior gauge constraint leaves ' +
        f"{gauss['physical_dimension']:,}" + ' basis states; imposing it first gives joint success ' +
        gauss['joint_success_from_uniform_ambient_input'] +
        '. Every failure is retained. The conditioned reference is stationary on the '
        'reduced joint source; the full record archive is not claimed stationary.', 19, SOFT, 1.4)
    c.link(X0 + 28, y + 180, 'Same-U preparation and complete outcome proof',
           '../../model/coherent_defect_wave/PREPARATION.md', 16)
    c.link(X0 + 805, y + 180, 'Interior gauge constraint', '../../model/coherent_defect_wave/GAUSS.md', 16)
    y += h + 22

    extensions = value['extensions']
    background = extensions['background']
    h, gap, w = 394, 22, (CW - 44) / 3
    for i, (tag, title, body, number, proof, color) in enumerate((
            ('Prepare → disturb → read', 'A fixed background with visible excitations',
             'Uniform link and wave amplitudes form a pure state fixed by the actual tick. '
             'A complete apparatus prepares it. Reading its defect face disturbs that state; '
             'a coherent filter then removes the unchanged background and retains every outcome.',
             'Fresh preparation success: ' + background['joint_success_per_fresh_ambient_trial'],
             'background', BLUE),
            ('Cell → neighbor → detector', 'Connected propagation with a finite cone',
             'Place cells on a declared cubic grid. A partial exchange carries wave information '
             'to a neighbor while keeping its link memory local and preserving its constraint. '
             'The arriving wave changes the neighbor’s defect motion. A finite cone bounds influence.',
             'Small-momentum spectator mode: quadratic', 'spatial', TEAL),
            ('One wave → several occupations', 'Several waves share the same links',
             'The occupation extension uses antisymmetric wave states and one shared defect register. '
             'Its one-wave sector is exactly the original program. Two-wave records are computed. '
             'Fermionic hopping and tensor-wave exchange are distinct spatial laws.',
             'Internal symmetry: U(3) × U(2)', 'species', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 13)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 58, number, 16, color, 700, maxw=w - 48)
        c.link(x + 24, y + h - 25, 'Exact construction and domain',
               '../../model/coherent_defect_wave/' + proof + '/PROOF.md', 16)
    y += h + 22

    cubic = extensions['cubic_readout']
    energies = extensions['dispersion']
    collision = extensions['energy']['finite_dark_collision']
    records = len(cubic['experiment']['complete_joint_records'])
    refined = len(cubic['higher_momentum_experiment']['complete_joint_records'])
    ratios = ' : '.join(f'{float(x):.6f}'.rstrip('0').rstrip('.')
                        for x in energies['three_nonzero_source_energy_ratios_approximate'])
    h, gap, w = 438, 22, (CW - 44) / 3
    for i, (tag, title, body, number, proof, color) in enumerate((
            ('Whole field → momentum record', 'Three spatial directions, one clock',
             'The cubic experiment translates each whole cell: its wave and retained link data. '
             'It records momentum and clock interference together, with every alternative kept. '
             'The larger cube adds complex translation phases; its complete probabilities are '
             'computed without materializing every branch vector. Preparation failures count too.',
             f'{records} records; {refined} on the larger cube', 'cubic_readout', BLUE),
            ('One explicit positive energy law', 'Exact energies; a failed relativistic test',
             'The declared positive clock has exact carrier energies and three-dimensional '
             'dispersion. Its three populated positive gaps have the ratios below. Their '
             'curvatures imply different limiting speeds; refining this additive branch '
             'does not give a common relativistic mass shell. These are source energies, not lepton masses.',
             ratios, 'dispersion', AMBER),
            ('Two waves → three detector outcomes', 'An exact finite collision',
             'Six actual source states form an invariant interacting sector of the same '
             'positive clock. From the specified pair preparation, two waves exchange '
             'occupation between cells. The detector '
             'keeps zero, one and two waves on the left. The source contact term changes '
             'this response compared with independent waves. The carrier flux stays in its cell.',
             'Maximum double occupation: ' +
             collision['maximum_double_occupation_probability'].replace('pi', 'π').replace('**2', '²'),
             'energy', TEAL))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 12)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 58, number, 16, color, 700, maxw=w - 48)
        c.link(x + 24, y + h - 25, 'Source, clock and complete readout',
               '../../model/coherent_defect_wave/' + proof + '/PROOF.md', 16)
    y += h + 22

    mobile = extensions['mobile_flux']
    paired, route = mobile['paired_flux'], mobile['protected_object_route']
    mass = route['positive_branch_rest_parameter']
    h, gap, w = 470, 22, (CW - 22) / 2
    for i, (tag, title, body, proof, color) in enumerate((
            ('Boundary data travel with the defect', 'The protected pair really changes cell',
             'The enlarged link domain includes a flat vacuum. An interface gate carries '
             'the interior flux and its exterior compensator together. Complete position '
             'records give equal left/right chances after one exchange and certain transfer '
             'after two. The moved link changes the receiver’s wave response. Five internal '
             'face locations describe one protected type, not five identified species.',
             'PROOF.md', AMBER),
            ('Same source object → free Dirac window', 'A rest frequency for that same flux pair',
             'Two routing lanes carry the literal paired-flux link code. Their local memory '
             'has ' + str(route['full_two_lane_local_dimension']) + ' states, including collisions. '
             'A chosen mass coin uses the actual count of ' + str(route['source_face_count']) +
             ' flipped faces. The same object’s position and Fourier readouts have a controlled '
             'free Dirac limit. Its spinor memory and coupling law are declared architecture.',
             'ROUTING_COMPLETION.md', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h, stroke=color)
        c.caps(x + 24, y + 35, tag, color, 12)
        yy = c.para(x + 24, y + 76, w - 48, title, 25, TEXT, 1.28, 700)
        c.para(x + 24, yy + 12, w - 48, body, 18, SOFT, 1.36)
        if i == 0:
            box_y = y + h - 143
            for box_x in (x + 44, x + w - 199):
                c.rect(box_x, box_y, 145, 62, 'none', EDGE, 1.7, 10)
                c.circle(box_x + 55, box_y + 24, 8, AMBER)
                c.circle(box_x + 90, box_y + 24, 8, 'none', AMBER, 2)
                c.line(box_x + 63, box_y + 24, box_x + 82, box_y + 24, AMBER, 1.4)
            c.arrow(x + 204, box_y + 25, x + w - 215, box_y + 25, color, 2.4, 9)
            c.text(x + 61, box_y + 52, 'inside + outside', 14, SOFT)
            c.text(x + w - 184, box_y + 52, 'inside + outside', 14, SOFT)
            c.text(x + 24, y + h - 55,
                   'P(left/right) = ' + paired['charge_position_records_after_one_exchange']['10'] +
                   ' / ' + paired['charge_position_records_after_one_exchange']['01'], 18, color, 700)
        else:
            c.text(x + 24, y + h - 115, 'ω² = |p|² + ' + str(mass*mass) + ';  |ω(0)| = ' +
                   str(mass) + ' (declared units)', 23, color, 700, maxw=w-48)
            c.para(x + 24, y + h - 79, w - 48,
                   'Bounded-window error: t·a·(3P² + ' + str(3*abs(mass)) +
                   'P). No GeV or Standard Model identity follows.',
                   17, MUTED, 1.3)
        c.link(x + 24, y + h - 25, 'Literal source encoding and proof',
               '../../model/coherent_defect_wave/mobile_flux/' + proof, 16)
    y += h + 22

    free=extensions['many_body_route']
    apparatus=extensions['mass_instrument']
    identity=extensions['field_identity_audit']
    rest=apparatus['continuum_rest_energy_in_source_units']
    gap=free['reference']['relative_charge_zero_gap']
    h,gap_x,w=434,22,(CW-44)/3
    for i,(tag,title,body,number,proof,color) in enumerate((
            ('A complete graded source', 'Positive particles and holes',
             'A separate declared local law retains all occupations of '
             +str(free['source']['literal_link_slots_per_cell'])+' graded flux slots per cell. '
             'Fill its negative-energy bands to obtain the lowest-energy reference. Adding an '
             'object or removing one costs positive energy; the missing object is a hole. '
             'Their free Dirac limit is controlled. The native gauge currents act trivially on these bare flux objects.',
             'Free rest energy '+rest+'; neutral-pair threshold '+gap,
             'many_body_route',BLUE),
            ('Actual source currents', 'Fifteen charged native states',
             'A two-address native wave combination is invariant under the idle law up to phase. '
             'Its even-occupation clouds form fifteen states with actual source currents and the '
             'Standard Model weight pattern. A declared gate transports the cloud with its '
             'graded flux object while conserving charge. The idle law gives zero interaction '
             'energy on this cloud and supplies no binding.',
             str(identity['native_dark_current_module']['actual_native_current_intertwiner_columns_checked'])+
             ' actual current identities; native K = 0', 'field_identity_audit',VIOLET),
            ('Two incoming objects, one link', 'A finite encounter with controlled error',
             'Added four-component memory puts both incoming objects on one clock with a '
             'common-speed free Dirac window. Their weak shared-link update changes the flux '
             'face with positive probability; complete detector bounds are retained. '
             'For one flux and one native occupation, bounded coupling gives vanishing '
             'contact response in a fixed-volume, fixed native-momentum window.',
             'Same source and clock; interacting limit open',
             'propagating_native_join',TEAL))):
        x=X0+i*(w+gap_x)
        card(c,x,y,w,h,stroke=color)
        c.caps(x+24,y+35,tag,color,12)
        yy=c.para(x+24,y+76,w-48,title,24,TEXT,1.28,700)
        c.para(x+24,yy+10,w-48,body,18,SOFT,1.36)
        c.text(x+24,y+h-58,number,16,color,700,maxw=w-48)
        c.link(x+24,y+h-25,'Exact source, reference and scope',
               '../../model/coherent_defect_wave/'+proof+'/PROOF.md',16)
    y+=h+22

    h=366
    card(c,X0,y,CW,h,stroke=BLUE)
    c.text(X0+28,y+43,'What does a source mass measurement look like?',28,TEXT,700)
    c.para(X0+28,y+83,CW-56,
        'A finite reservoir adds or removes an object from the filled reference. A battery '
        'supplies the work; the apparatus measures energy and relative pair count using the '
        'same Hamiltonian. Failed preparations are retained and can create neutral pairs. '
        'This gives an operational source mass, rather than reading a phase and calling it energy.',
        19,SOFT,1.4)
    diagram_y=y+180
    boxes=((X0+28,250,'Reservoir + battery','Prepare and supply work',TEAL),
           (X0+395,400,'Filled source reference','Add a particle / remove a hole',BLUE),
           (X0+910,350,'Energy + charge records',
            str(len(apparatus['experiments'][0]['complete_records']))+' outcomes, including failures',VIOLET))
    for xx,ww,title,subtitle,color in boxes:
        c.rect(xx,diagram_y,ww,72,'none',color,1.7,10)
        c.text(xx+ww/2,diagram_y+29,title,19,color,700,anchor='middle',maxw=ww-20)
        c.text(xx+ww/2,diagram_y+54,subtitle,15,SOFT,anchor='middle',maxw=ww-20)
    c.arrow(X0+290,diagram_y+36,X0+380,diagram_y+36,TEAL,2.4,9)
    c.arrow(X0+810,diagram_y+36,X0+895,diagram_y+36,BLUE,2.4,9)
    c.text(X0+28,y+294,'E² = |p|² + '+str(int(rest)**2)+
           '; rest energy = '+rest+' source units in this free limit',23,BLUE,700,maxw=CW-56)
    c.text(X0+28,y+326,'GeV value = '+rest+
           ' × energy unit in GeV. A physical clock reference fixes that unit; no measured particle is assigned.',
           17,MUTED,maxw=CW-56)
    c.link(X0+28,y+h-15,'Complete reservoir, battery and mass records',
           '../../model/coherent_defect_wave/mass_instrument/PROOF.md',16)
    y+=h+22

    field_map=extensions['chiral_particle_map']
    spectrum=field_map['background_spectrum']
    finite=extensions['chiral_joint/periodic_spectrum']
    first=finite['compensated_charged_spectroscopy']['charged_success_branches'][0]['complete_energy_windows'][0]
    interval=list(map(s.Rational,first['energy_transfer_interval']))
    threshold=f'{float(sum(interval)/2):.8f}'
    h=678
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'A chiral field dictionary from the actual source currents',28,TEXT,700)
    c.para(X0+28,y+83,CW-56,
        'A declared cell encoding retains all 30 spin/native fermion modes and every occupation. '
        'Its spectral left-Weyl propagation has a controlled free limit; its finite position '
        'kernel is nonlocal. In a specified neutral Higgs background, the source’s pair tensors '
        'join left fields to conjugate singlets. This is an exact change of field variables.',
        19,SOFT,1.4)
    table_y=y+193
    c.rect(X0+24,table_y,CW-48,182,'#142337',EDGE,1,8)
    for xx,label in ((X0+40,'Conditional field role'),(X0+540,'Electric charge'),(X0+800,'Rest parameter')):
        c.text(xx,table_y+28,label,18,VIOLET,700)
    roles=(('Up-type quark · 3 colors','2/3','up'),
           ('Down-type quark · 3 colors','−1/3','down'),
           ('Charged lepton · 1 field','−1','charged_lepton'),
           ('Neutral left lepton · 1 field','0','neutral'))
    for i,(role,charge,key) in enumerate(roles):
        yy=table_y+62+32*i
        c.text(X0+40,yy,role,20,TEXT,700)
        c.text(X0+540,yy,charge,20,SOFT)
        mass=spectrum['background_mass_energy'][key].replace('abs(', '|').replace(')', '|').replace('*','·')
        c.text(X0+800,yy,mass+' source units',20,SOFT)
    c.para(X0+28,y+408,CW-56,
        'In this background field dictionary, all fermion fields have spin ½. '
        'Weak currents act only on left components; color and electric '
        'charge act on both halves of the massive fields. Antiparticles have conjugate charges. '
        'The seven Dirac color/lepton components plus one neutral Weyl field form ONE family. '
        'No electron/muon/tau generation is selected, and v and the three Yukawa coefficients are inputs.',
        19,SOFT,1.4)
    c.para(X0+28,y+527,CW-56,
        'The periodic quantum lepton control has 3,424 neutral physical states and one certified ground. '
        'Its 2,627 complete records retain preparation failure (success 1/3424). Addition and removal '
        'both read '+threshold+' source units. An actual Higgs composite excites this same line, '
        'without a fitted scalar background. Joint source/reference phases supply only one squared momentum. '
        'Quarks are spectators here; relativistic masses, generations and GeV values are unproved.',
        18,MUTED,1.38)
    c.link(X0+28,y+h-25,'Source IDs, all currents and conditional mass formulas',
           '../../model/coherent_defect_wave/chiral_particle_map/PROOF.md',16)
    c.link(X0+805,y+h-25,'Quantum state and complete records',
           '../../model/coherent_defect_wave/chiral_joint/periodic_spectrum/PROOF.md',16)
    y += h + 35

    neutral = extensions['chiral_joint/neutral_composite']['complete_readout']
    bosons = extensions['chiral_joint/periodic_boson_records']['complete_readout']
    neutral_line = neutral['branches'][0]['complete_translation_branches'][0]['complete_energy_records'][0]
    boson_line = bosons['settings'][0]['lowest_certified_positive_inelastic_record']
    middle = lambda row: f"{float(sum(map(Fraction,row['energy_transfer_interval']))/2):.8f}"
    h=300
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'Several detectors read the same interacting source',28,TEXT,700)
    c.para(X0+28,y+84,CW-56,
        'The scalar, link loops and neutral composite are measured on the same 3,424-state periodic ground. '
        'Success and failure records are retained; scalar and loop probes keep elastic events. Preparation failure is retained too.',
        19,SOFT,1.4)
    c.text(X0+28,y+159,'Neutral scalar–lepton probe: '+middle(neutral_line)+' source energy',23,TEAL,700,maxw=CW-56)
    c.text(X0+28,y+197,'Scalar and two loop probes share a line at '+middle(boson_line),23,VIOLET,700,maxw=CW-56)
    c.para(X0+28,y+233,CW-56,
        'These finite energy transfers connect actual source operators. They do not assign a neutrino, Higgs, photon, W or Z mass.',
        18,MUTED,1.35)
    c.link(X0+28,y+h-15,'Complete neutral-composite records',
           '../../model/coherent_defect_wave/chiral_joint/neutral_composite/PROOF.md',16)
    c.link(X0+650,y+h-15,'Complete scalar and link-loop records',
           '../../model/coherent_defect_wave/chiral_joint/periodic_boson_records/PROOF.md',16)
    y+=h+22

    motion=extensions['chiral_joint/peter_weyl_links']['complete_neutral_transfer_experiment']
    family=extensions['chiral_joint/family_transport']['same_clock_full_H_probability_records']
    h=655
    card(c,X0,y,CW,h,stroke=TEAL)
    c.text(X0+28,y+43,'Quantum links let a neutral packet move',28,TEXT,700)
    c.para(X0+28,y+83,CW-56,
        'The full Peter–Weyl link rule includes 27 transitions between excited flux sectors. '
        'Those transitions release the conserved star tag of the minimal link rule. '
        'The source then transports a scalar–lepton packet using its actual Yukawa and hopping terms. '
        'No extra neutral hopping term is inserted.',19,SOFT,1.4)
    w=(CW-104)/3
    labels=(('Prepare a neutral packet','Scalar and lepton form a gauge singlet. '
             'Keep both spin branches and the fermionic reference.'),
            ('Let the source act','119 states per link. Five actual H factors '
             'give the first neighboring transfer; every matter block contributes.'),
            ('Read the next cell','Opposite-spin arrival has a strictly positive '
             'full-H bound. The complete no-arrival record is retained.'))
    for i,(title,body) in enumerate(labels):
        xx=X0+28+i*(w+24)
        c.rect(xx,y+181,w,144,'#142337',EDGE,1,8)
        c.text(xx+14,y+212,title,19,TEAL,700,maxw=w-28)
        c.para(xx+14,y+242,w-28,body,17,SOFT,1.32)
        if i<2:c.text(xx+w+5,y+259,'→',19,TEAL,700)
    c.text(X0+28,y+369,'All five matter representations also move with their compensating flux',23,TEAL,700,maxw=CW-56)
    c.para(X0+28,y+408,CW-56,
        'Up-conjugate quarks, the quark doublet, charged-lepton singlet, lepton doublet and '
        'down-conjugate quarks share one Hamiltonian and one dual-native reference. '
        'Every full preparation includes its gauge and spin partners. A local detector distinguishes '
        'prepared arrival from the preparation-off background at the same time.',19,SOFT,1.4)
    c.text(X0+28,y+517,'Common time '+family['common_source_time']+' · positive local contrast ≥ '+
        family['strict_local_record_contrast_lower_bound'],21,TEAL,700,maxw=CW-56)
    c.para(X0+28,y+553,CW-56,
        'These are source representation names and compensated finite motion. '
        'A stationary moving branch, physical mass shell and measured species require their own identification. '
        'The exact finite rational clock preserves the detector contrast.',18,MUTED,1.35)
    c.link(X0+28,y+h-18,'Actual link law and neutral propagation',
           '../../model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md',16)
    c.link(X0+660,y+h-18,'All-family transport and complete local detectors',
           '../../model/coherent_defect_wave/chiral_joint/family_transport/PROOF.md',16)
    y+=h+22

    local=extensions['chiral_joint/local_reference']
    gap=local['global_physical_ground_and_gap']['global_physical_gap_interval']
    h=420
    card(c,X0,y,CW,h,stroke=BLUE)
    c.text(X0+28,y+43,'Which interacting local state is the reference?',28,TEXT,700)
    c.text(X0+28,y+92,'All 5 × 2³⁰ local states  →  exact symmetry bounds  →  a ground doublet',25,BLUE,700,maxw=CW-56)
    c.para(X0+28,y+132,CW-56,
        'Exact collective algebra compares every zero-boundary-charge physical sector, including states '
        'invisible to the empty preparation. The decisive module has 104 physical states and 56 dark '
        'energy-one states. Its true physical local ground has spin ½ and energy about −11.23059 source units.',
        19,SOFT,1.4)
    c.text(X0+28,y+232,'First local excitation gap: '+f'{float(sum(map(Fraction,gap))/2):.11f}'+' source units',24,BLUE,700,maxw=CW-56)
    c.para(X0+28,y+270,CW-56,
        'The stationary ground density has mean native occupation fifteen. This filled interacting reference '
        'is not a claim of fifteen particles. The empty-register cyclic ground cannot be the universal vacuum; '
        'the lower reference energy is not a negative particle mass.',19,SOFT,1.4)
    c.para(X0+28,y+360,CW-56,
        'The controlled many-cell construction below includes link-compensated charges and actual hopping.',18,MUTED,1.35)
    c.link(X0+28,y+h-15,'Complete local ground, dark states, gap and occupations',
           '../../model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',16)
    y+=h+22

    h=355
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'A matter excitation leaves a closed gauge-loop record',28,TEXT,700)
    sx,sy=X0+75,y+111
    for dx,dy in ((0,0),(124,0),(124,124),(0,124)):
        c.circle(sx+dx,sy+dy,9,VIOLET)
    for dx1,dy1,dx2,dy2 in ((0,0,124,0),(124,0,124,124),(124,124,0,124),(0,124,0,0)):
        c.line(sx+dx1,sy+dy1,sx+dx2,sy+dy2,VIOLET,3)
    c.text(sx+42,sy-18,'1 →',20,VIOLET,700)
    c.text(sx+143,sy+69,'2 ↓',20,VIOLET,700)
    c.text(sx+43,sy+154,'← 3',20,VIOLET,700)
    c.text(sx-46,sy+69,'↑ 4',20,VIOLET,700)
    c.para(X0+300,y+91,CW-336,
        'Four actual hops bring the matter excitation back to its starting cell. '
        'The quantum links keep a closed, gauge-invariant trace. All spin signs, internal '
        'indices and the compensating reference are included. A contractible square is '
        'generated by the matter Hamiltonian; no magnetic loop term is inserted.',19,SOFT,1.4)
    c.text(X0+300,y+220,'Fourth moments: −1/256 (singlet), −1/768 (colored)',22,VIOLET,700,maxw=CW-336)
    c.para(X0+300,y+257,CW-336,
        'Complete loop and no-loop probabilities are bounded on the same source. '
        'A physical photon/gluon phase and a relativistic mass spectrum do not follow from this finite loop alone.',18,MUTED,1.35)
    c.link(X0+28,y+h-18,'Actual square, complete records and eliminated-state proof',
           '../../model/coherent_defect_wave/chiral_joint/matter_induced_plaquette/PROOF.md',16)
    c.link(X0+830,y+h-18,'Finite unitary clock and its error',
           '../../model/coherent_defect_wave/chiral_joint/finite_clock/PROOF.md',16)
    y+=h+22
    h=535
    card(c,X0,y,CW,h,stroke=TEAL)
    c.text(X0+28,y+43,'A finite interacting ground has a visible collective excitation',28,TEXT,700)
    c.para(X0+28,y+84,CW-56,
        'Eight cells retain every matter occupation, scalar state and quantum link. '
        'At a declared weak hopping strength, exact local bounds and a complete charged-state '
        'resolvent prove a unique physical ground. The first excited cluster has three cubic momenta. '
        'Its finite gap is computed with a bound for every complementary state.',19,SOFT,1.4)
    stages=[('PREPARE','Physical controller','bounded ground error'),('INTERVENE','Spin pump','all failures retained'),
            ('RESOLVE','Three X momenta','line and complement')]
    w=(CW-100)/3
    for i,(label,title,detail) in enumerate(stages):
        x=X0+32+i*(w+18)
        c.text(x,y+210,label,15,TEAL,700)
        c.text(x,y+243,title,23,TEXT,700)
        c.text(x,y+273,detail,18,SOFT)
        if i<2:c.text(x+w-12,y+243,'→',25,TEAL,700)
    c.text(X0+28,y+323,'Detected line: [8.20758, 8.21035] × 10⁻²² source energy units',25,TEAL,700,maxw=CW-56)
    c.para(X0+28,y+360,CW-56,
        'At hopping 2⁻³⁴, the line probability is 0.00204894–0.00211091 after successful ground '
        'preparation. A positive outside-line response is also bounded. Three reference axes, '
        'pump failure and preparation failure make a complete probability law.',19,SOFT,1.4)
    c.para(X0+28,y+447,CW-56,
        'This is an interacting collective energy record at nonzero momentum. '
        'A rest mass, measured species and GeV conversion require their own physical identification.',18,MUTED,1.35)
    c.link(X0+28,y+h-18,'Full interacting ground and momentum proof',
           '../../model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/PROOF.md',16)
    c.link(X0+730,y+h-18,'Complete source, apparatus and detector records',
           '../../model/coherent_defect_wave/chiral_joint/collective_spin_response/PROOF.md',16)
    y+=h+35
    charged=data['joint']['extensions']['chiral_joint/interacting_charged_response']
    h=675
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'Move a fermion into an apparatus; retain its charge',28,TEXT,700)
    c.para(X0+28,y+84,CW-56,
        'The same interacting ground admits a physical preparation from blank states. '
        'Retained spin references and a finite two-copy energy filter leave source error below '
        '2.36 × 10⁻⁷. The declared algebraic controller succeeds with probability above 0.99999994; '
        'its very large finite schedule is explicit. Every failure still satisfies Gauss constraints.',19,SOFT,1.4)
    c.circle(X0+100,y+252,42,TEAL)
    c.text(X0+100,y+259,'source',19,TEXT,700,'middle')
    c.arrow(X0+165,y+252,X0+405,y+252,VIOLET,3,12)
    c.text(X0+285,y+225,'one fermion transfer',18,VIOLET,600,'middle')
    c.rect(X0+425,y+208,228,89,'#151b2e',VIOLET,1.5,14)
    c.text(X0+539,y+246,'quantum reference',20,TEXT,700,'middle')
    c.text(X0+539,y+276,'30 modes, all occupations',16,SOFT,400,'middle')
    c.para(X0+700,y+224,CW-740,
        'Joint charge and total fermion parity stay conserved. The apparatus keeps the '
        'compensating fermion. Its energy cost is four source units, accounted separately.',18,SOFT,1.4)
    c.text(X0+28,y+351,'ACTUAL SOURCE REMOVAL',15,VIOLET,700)
    c.text(X0+400,y+351,'SOURCE ENERGY BIN',15,VIOLET,700)
    c.text(X0+835,y+351,'PROBABILITY IN THAT BIN',15,VIOLET,700)
    for i,(key,label) in enumerate((('ec','eᶜ · charged singlet'),('L','L · weak lepton doublet'))):
        row=charged['same_H_interacting_channels'][key]['full_interacting_energy_readout']
        energy=[float(Fraction(x)) for x in row['source_excitation_energy_bin_after_subtracting_declared_reference_cost']]
        probability=[float(Fraction(x)) for x in row['conditional_probabilities']['line']]
        yy=y+394+i*46
        c.text(X0+28,yy,label,21,TEXT,600)
        c.text(X0+400,yy,f'{energy[0]:.6f} – {energy[1]:.6f}',21,TEAL,600)
        c.text(X0+835,yy,f'{probability[0]:.6f} – {probability[1]:.6f}',21,TEAL,600)
    c.para(X0+28,y+488,CW-56,
        'These probabilities use the exact ground. The physical preparation adds its error bound. '
        'The remaining outcomes include the complete outside-bin signal and pump failure. '
        'No rare outcome is deleted. Total source-plus-apparatus energies add four to each bin.',19,SOFT,1.4)
    c.para(X0+28,y+575,CW-56,
        'This is a charged detector transition with a static compensator. It does not yet '
        'establish a freely propagating electron, neutrino or a rest mass in GeV.',18,MUTED,1.35)
    c.link(X0+28,y+h-20,'Physical preparation and finite error proof',
           '../../model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/PROOF.md',16)
    c.link(X0+740,y+h-20,'Full charged response and reference energy',
           '../../model/coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md',16)
    y+=h+35
    c.link(X0 + 24, y, 'Joint evolution and reciprocal detector proofs',
           '../../model/coherent_defect_wave/PROOF.md', 17)
    c.link(X0 + 805, y, 'Exact joint receipt', '../../model/coherent_defect_wave/results.json', 17)
    return y + 22
