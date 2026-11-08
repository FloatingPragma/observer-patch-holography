"""Catalogue protected objects with their exact completeness domains."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from model.completion import defect_algebra
from model.defect_protection import protection
from model.spin_joint_field_limit import packet as joint_field
from model.source_family_compatibility import hole_records, hole_dispersion, local_matter, charge_bands, composite
from . import pauli, charge, field, mechanism

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = tuple('model/'+name+'/results.json' for name in (
    'defect_protection', 'local_exchange_charge', 'common_spin_protection',
    'cube_eta_protection', 'compact_defect_pair', 'defect_localization_readout',
    'spin_canonical_field', 'spin_canonical_energy', 'spin_joint_field_limit',
    'source_family_compatibility'))


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def serialize(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {key:serialize(item)for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(item)for item in value]
    return value


def bindings():
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        pins = parent.get('source_bindings', parent.get('source_sha256'))
        if not isinstance(pins, dict) or not pins:
            raise ValueError('catalogue parent has no source custody: '+name)
        for path, digest in pins.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale catalogue parent: '+path)
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/protected_defect_catalog/PROOF.md',
        'source_receipts': list(PARENTS),
        'independent_checks': ['model/protected_defect_catalog/README.md']})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def rows():
    return [
        {'id':'bare_native_source', 'kind':'exhaustive_exclusion',
         'source':'complete all-seam C/E on H66 or H78; source-only charges',
         'protected_object':'only scalar branch-commuting charges',
         'completeness':'entire branch commutant, every quantum basis',
         'preparation_readout':'native complete occupation/current records',
         'propagation_energy':'no nontrivial charge supplied by this classification',
         'topology':'no protected source-only reducing charge in this domain',
         'boundary':'history, coherent memory, recovery, partial seam support and refinement are outside this domain',
         'evidence':'model/completion/defect_algebra.py'},
        {'id':'ordered_exchange_charge', 'kind':'protected_memory_charge',
         'source':'H66 with a coherent orientation bit; complete ordered-cover repair',
         'protected_object':'Q=I66 tensor X, two sharp eigenvalue sectors',
         'completeness':'entire all-seam commutant span{I,Q}',
         'preparation_readout':'coherent memory preparation; complete charge pointer; local qutrit realization retains the same Q',
         'propagation_energy':'exact pair-configuration transport; no relativistic rest-energy map',
         'topology':'Z2 character of the ordered exchange cover',
         'boundary':'sharp charge is not strictly localized with reference exterior; local nonsharp excess is available; memory dephasing destroys protection',
         'evidence':'model/defect_protection/README.md',
         'same_object_realizations':['model/local_exchange_charge/README.md','model/protected_reference/PROOF.md']},
        {'id':'selected_spin_charge', 'kind':'conserved_symmetry_charge',
         'source':'actual eight-site selected spin cube, J=1, B=1/4; inherited two-flavor repair',
         'protected_object':'integer U1 spin charge Q; global X parity is exp(i pi Q); local Z creates two sharp Q=+/-1 channels',
         'completeness':'all 65536 Pauli-product words plus exact full charge spectral resolution; not all conserved operators',
         'preparation_readout':'complete native local-Z followed by one X pointer prepares coherent Q=+/-1 branches; separate ideal 256-outcome X audit resolves charge weights and dephases within sectors',
         'propagation_energy':'positive finite same-H_B energy; contrast bound through ten repairs for T<=1/16 under the stated source scaling',
         'topology':'global internal U1 symmetry, not a winding number or identified electric charge',
         'boundary':'no nonidentity proper-support Pauli charge; coherent controls can change parity; no permanent localization or particle pole',
         'evidence':'model/defect_localization_readout/PROOF.md'},
        {'id':'single_occupancy_spin_memory', 'kind':'repair_protected_algebra',
         'source':'eight-site one-particle-per-site spin code with complete C/E repair',
         'protected_object':'full M256 operator algebra under repair alone',
         'completeness':'entire operator algebra of this code under the specified repair',
         'preparation_readout':'declared spin control and detector services',
         'propagation_energy':'exchange evolution acts nontrivially within the code',
         'topology':'no topological label supplied by identity repair',
         'boundary':'repair protection of every state does not make every state an energy eigenstate or a distinct species',
         'evidence':'model/defect_localization_readout/PROOF.md'},
        {'id':'joint_positive_field_canonical_sector', 'kind':'constructive_field_family',
         'source':'protected spin interaction H_B=J sum(S_x.S_y-1/4)-B sum epsilon_x S_x^1; supplied p=n^-2 and B=J n^-5; positive-B ground and retained H_B clock',
         'protected_object':'two regular canonical pairs (Q3,P2) and (Q2,P3), with common positive-energy time and full internal U1 on their bounded histories',
         'completeness':'specified field subnet along a common history subsequence; neither all continuum observables nor all protected defects are classified',
         'preparation_readout':'finite positive-field ground preparation, local marginal modulus, source spin-rotation observables and complete native local Pauli histories with all pointer and failure records',
         'propagation_energy':'continuous even-envelope translations, positive generator K_B/(Jp), and finite propagation cone; the volume-wide uniform mode has inverse-weighted centre and energy-square variance tending to zero, with a volume-dependent original probability bound; no nonzero-momentum shell or rest mass',
         'topology':'internal source spin rotation, no spatial winding or identified electromagnetic charge',
         'boundary':'B(n) is a supplied preparation prescription; unread repair does not select it; equality to the ordered zero-field history limit, physical Lorentz spin, sharp poles and scattering are unproved',
         'evidence':'model/spin_joint_field_limit/PROOF.md',
         'same_object_realizations':['model/spin_joint_field_limit/SECOND_PAIR.md',
                                     'model/spin_joint_field_limit/CLOCK_COMPARISON.md',
                                     'model/spin_joint_field_limit/UNIFORM_MODE.md']},
        {'id':'common_hidden_spin_factor', 'kind':'joint_noiseless_subsystem',
         'source':'same six-site Weyl/contact field and complete charged-block repair',
         'protected_object':'at N=4, S=1, a qutrit with 105 multiplicity states, total dimension315',
         'completeness':'exhibited common spin factor; maximal joint noiseless algebra not classified',
         'preparation_readout':'complete ideal joined histories with arbitrary spin/reference correlations',
         'propagation_energy':'multiplicity evolves and localized carrier leaks; repair need not conserve field energy',
         'topology':'internal spin factor, no topological boundary invariant',
         'boundary':'global factor preservation does not identify localized matter, Lorentz spin or a physical spin apparatus',
         'evidence':'model/common_spin_protection/README.md'},
        {'id':'cube_eta_range', 'kind':'joint_protected_sector',
         'source':'original eight-site Weyl/contact field with the stated four-dimer repair',
         'protected_object':'rank1820 image of eta+ from N4 to N6; eta>=2 sector',
         'completeness':'entire image of this injective map; all protected sectors not classified',
         'preparation_readout':'source-prepared eta reference and complete local response detector',
         'propagation_energy':'same interacting action intertwines with an energy shift U; one probe has eight inelastic lines plus elastic weight',
         'topology':'eta symmetry representation, no spatial winding invariant',
         'boundary':'eight response lines are not eight species; protection uses the complete outside outcome for this repair domain',
         'evidence':'model/cube_eta_protection/README.md'},
        {'id':'compact_opposite_texture_cores', 'kind':'specified_topological_background',
         'source':'prescribed two-core texture in the declared cubic two-spin field',
         'protected_object':'opposite shell degrees +1 and -1, constant exterior, four localized band modes',
         'completeness':'specified texture family and certified isolated band, not all native defects',
         'preparation_readout':'declared texture and transported band controls; autonomous source formation unproved',
         'propagation_energy':'uniform finite texture/vacuum-energy bounds and controlled band transport',
         'topology':'genuine shell winding of the prescribed order-parameter field',
         'boundary':'winding does not force zero modes, four generations or physical rest masses; stability needs the stated dynamics',
         'evidence':'model/compact_defect_pair/PROOF.md'},
        {'id':'one_hole_spin_readout', 'kind':'scoped_repair_sector',
         'source':'full two-flavour CAR repair on H_h: one empty site and one fermion at every other site; dimension V 2^(V-1)',
         'protected_object':'hole-position projectors and internal spin information carried by the recorded signed transport; conditional block dimension 2^(V-1)',
         'completeness':'all H_h inputs and spin/reference matrix elements under the stated C/E repair; no exhaustive classification of all source states',
         'preparation_readout':'supplied H_h input; complete E seam-cover records identify the initial hole; inverse signed transport gives sum_x P_x rho P_x and preserves each spin/reference block',
         'propagation_energy':'uniform global repair attempts have complete source/archive disturbance relative to the same propagation with identity repair at most 4/n on cubes of side n^2+1 with at most n^2 attempts and common H_h-preserving propagation',
         'topology':'one-unit occupation deficit, without a topological invariant or physical particle assignment',
         'boundary':'general Hubbard hopping leaves H_h through doublons; the fully aligned invariant restriction is free and has a separate vacancy-position, magnetic-translation and nonrelativistic continuum map; its reference gap is not a full-ground gap or physical mass; the matched repair chart, preparation and schedule are supplied',
         'evidence':'model/source_family_compatibility/REPAIR_TRAJECTORY.md',
         'same_object_realizations':['model/source_family_compatibility/ALIGNED_HOLE_DISPERSION.md']},
        {'id':'centered_full_CAR_charge_response', 'kind':'scoped_interacting_response_family',
         'source':'fixed-coupling original Weyl/contact action with bare full-domain repair; invariant centered full-ground subsequence; the signed Hubbard chart is a separate direct local instance',
         'protected_object':'odd addition/removal responses and even spin observables on the same full CAR algebra; total occupation is conserved under source propagation and repair',
         'completeness':'complete local observer records for the specified accessible algebra, not a census of all source excitations or the full repair archive',
         'preparation_readout':'local addition and removal each have original weight 1/2; odd transfers require a declared fermionic reference service; local finite-time response scores retain all local outcomes and finite errors',
         'propagation_energy':'same positive centered ground clock; both charges have a volume-independent positive original weight in an explicit finite-energy band; a nonzero connected four-point witness proves interacting local histories; source composite operators give a measured response block and a controlled spectral remainder',
         'topology':'source U1 charge grading, without an assigned winding invariant or Standard Model charge',
         'boundary':'the invariant ground subsequence has no effective preparation modulus; inaccessible remote raw values cannot drive feedback; a global sharp energy PVM is not supplied by the local bound; no selected relativistic scaling, particle pole or physical spin is established',
         'evidence':'model/source_family_compatibility/CHARGED_SPECTRAL_BANDS.md',
         'same_object_realizations':['model/source_family_compatibility/LOCAL_MATTER_READOUT.md',
                                     'model/source_family_compatibility/COMPOSITE_CAR.md']},
    ]


def joint_field_attachment():
    """Consume the owning source proof and complete finite record producer."""
    value = joint_field.build()
    if joint_field.canonical(value) != joint_field.canonical(strict_json.load(joint_field.OUTPUT)):
        raise ValueError('retained joint-field source receipt differs')
    scope = value['scope']
    required = ('sequence_is_a_chosen_source_preparation_prescription',
                'finite_positive_field_preparation_and_local_accuracy_computable',
                'second_canonical_pair_on_same_joint_history_subsequence',
                'full_internal_U1_on_both_canonical_pairs',
                'fixed_scale_complete_prepared_H_B_probe_history_constructed',
                'complete_native_failure_records_retained')
    if not all(scope[name] is True for name in required):
        raise ValueError('joint-field source attachment lacks an owning premise')
    exclusions = ('sequence_forced_by_unread_repair',
                  'equality_to_ordered_zero_field_state_clock_limit_proved',
                  'joint_continuum_subsequence_uniqueness_proved',
                  'sharp_particle_mass_shell_proved', 'Lorentz_covariance_proved')
    if any(scope[name] is not False for name in exclusions) or scope['physical_particle_assignment'] is not None:
        raise ValueError('joint-field attachment exceeds the owning physical scope')
    history = value['complete_prepared_local_history']
    uniform=value['uniform_mode_original_weight_and_concentration']
    spectral=value['complete_finite_uniform_spectral_readout']
    if (uniform['uniform_mode_centre_and_weighted_variance_vanish_on_joint_clock']is not True
            or uniform['uniform_original_probability_floor_proved']is not False
            or uniform['nonzero_momentum_or_physical_mass_shell_proved']is not False
            or spectral['fixed_box_table_substituted_for_growing_volume_field']is not False
            or spectral['physical_mass_prediction_count']!=0):
        raise ValueError('joint-field uniform spectral response exceeds its owning scope')
    return {'catalogue_entry':'joint_positive_field_canonical_sector',
            'owner_receipt':'model/spin_joint_field_limit/results.json',
            'owning_scientific_producer_replayed':True,
            'chosen_source_family':scope['source_field_sequence'],
            'state':scope['state'],'clock':scope['clock'],
            'source_observables':{'Q3':'p^2 sum epsilon_x F(px) S_x^3',
                                  'P2':'-p sum G(px) S_x^2',
                                  'Q2':'p^2 sum epsilon_x F(px) S_x^2',
                                  'P3':'+p sum G(px) S_x^3'},
            'canonical_pairs':[['Q3','P2'],['Q2','P3']],
            'symplectic_coefficient':'the same selected source order m>1/10',
            'internal_symmetry':'full U1 on both pairs and their same-clock bounded histories',
            'source_to_observable':'finite source spin rotations implement the field Weyl factors; native Pauli records read the prepared source',
            'repair_identity_generates_the_H_B_clock':False,
            'complete_finite_history_example':{
                'family_index':history['n'],
                'ideal_pointer_record_count':history['ideal_pointer_record_count'],
                'complete_output_half_trace_upper':history['complete_output_half_trace_upper'],
                'requested_accuracy':history['requested_accuracy'],
                'all_native_failure_and_pointer_records_retained':history['all_native_failure_and_pointer_records_retained'],
                'large_history_executed':history['large_history_executed']},
            'ordered_zero_field_history_identity_proved':False,
            'uniform_mode_response':{'weighted_centre_and_variance_tend_to_zero':True,
                'original_window_probability_lower':'m/(8N); finite order floor m>1/20',
                'uniform_original_continuum_residue_lower_proved':False,
                'nonzero_momentum_shell_proved':False,
                'finite_source_scale_indices':[row['n']for row in spectral['finite_cases']],
                'fixed_cube_spectrum_substituted_for_growing_field':False},
            'physical_species':[],'physical_masses':[],
            'full_source_field_classification_proved':False}


def hole_source_attachment():
    """Replay the owning full-record sector calculation, preserving its scope."""
    parent=strict_json.load(ROOT/'model/source_family_compatibility/results.json')
    value=hole_records.certificate()
    if canonical(value)!=canonical(parent['actual_complete_one_hole_repair']):
        raise ValueError('one-hole owning source receipt differs from the replay')
    if value['arbitrary_reference_preserved_in_each_internal_spin_block'] is not True:
        raise ValueError('one-hole attachment lacks complete spin/reference preservation')
    if (value['source']!='model.defect_field_bridge.full_domain.branch_action'
            or value['E_scan_corrected_complete_channel']!='sum_x P_x rho P_x'
            or value['Poisson_rate_supplied_not_selected_by_repair'] is not True):
        raise ValueError('one-hole source, position instrument or clock scope changed')
    excluded=('record_erasure_restores_hole_coherence',
              'one_hole_ground_or_relativistic_particle_assigned',
              'complete_word_simplex_and_counting_reference_premises_discharged',
              'correction_native_pulse_tape_executed')
    if any(value[key] is not False for key in excluded):
        raise ValueError('one-hole attachment exceeds the owning physical or execution scope')
    controls=value['global_fair_attempt_complete_record_workaround']
    if any(row['physical_clock_or_capacity_trajectory_selected'] is not False for row in controls):
        raise ValueError('one-hole attachment promotes a supplied source trajectory')
    propagation=value['actual_full_Hubbard_propagation_scope']
    if propagation['aligned_sector_identified_as_interacting_ground'] is not False:
        raise ValueError('one-hole aligned source is not an interacting ground attachment')
    return {'catalogue_entry':'one_hole_spin_readout',
            'owner_receipt':'model/source_family_compatibility/results.json',
            'owning_hole_scientific_producer_replayed':True,
            'full_owning_program_countercontrol_replayed':False,
            'source_domain':'H_h: one hole, one fermion on each other site; arbitrary internal spin and reference',
            'state':'supplied normalized input on H_h; no ground or vacuum selected',
            'clock':'literal global attempt count; any common H_h-preserving coherent propagation clock is supplied',
            'position_observable':'P_x projects onto the initially empty source site',
            'complete_corrected_position_instrument':value['E_scan_corrected_complete_channel'],
            'unconditional_hole_coherence_preserved_by_complete_scan':False,
            'arbitrary_reference_preserved_in_each_internal_spin_block':True,
            'signed_transport_correction_pulse_tape_executed':False,
            'allocated_zero_records':value['zero_history_representation'],
            'exact_source_cases':value['cases'],
            'global_attempt_complete_trace_bound':'min(2,2m sqrt((1+theta)d_max/|E|))',
            'counted_cube_trace_bound':'min(2,4/n); side n^2+1 and m<=n^2',
            'complete_source_archive_controls':controls,
            'actual_Hubbard_propagation_scope':propagation,
            'topological_invariant_assigned':False,
            'physical_species':[],'physical_masses':[],
            'all_source_states_classified':False}


def aligned_hole_attachment():
    parent=strict_json.load(ROOT/'model/source_family_compatibility/results.json')
    value=serialize(hole_dispersion.certificate())
    if canonical(value)!=canonical(parent['aligned_hole_position_dispersion_and_continuum']):
        raise ValueError('aligned-hole owning source receipt differs from the replay')
    if (value['full_interacting_ground_reference_selected']is not False
            or value['source_capacity_or_time_trajectory_selected']is not False
            or value['physical_particle_labels']or value['physical_masses']):
        raise ValueError('aligned-hole attachment exceeds its reference or physical scope')
    return {'catalogue_entry':'one_hole_spin_readout',
            'owner_receipt':'model/source_family_compatibility/results.json',
            'owning_scientific_producer_replayed':True,
            'source_chart':value['source_chart'],
            'reference':'filled hidden-up vector in the aligned Fock sector; not the full interacting ground',
            'clock':'centered signed Hubbard clock relative to that reference',
            'all_hole_numbers_have_exact_exterior_power_dynamics':True,
            'complete_vacancy_configuration_detector':value['complete_position_scans'],
            'position_translation':value['operational_position_translation'],
            'native_analyzer':value['native_complete_energy_analyzer'],
            'chosen_reference_gap':value['chosen_reference_gap'],
            'nonrelativistic_continuum_budget':value['continuum_budget'],
            'interacting_scattering_or_relativistic_mass_proved':False,
            'physical_species':[],'physical_masses':[]}


def centered_matter_attachment():
    parent=strict_json.load(ROOT/'model/source_family_compatibility/results.json')
    local=local_matter.certificate();bands=charge_bands.certificate()
    dressed=serialize(composite.certificate())
    if (canonical(local)!=canonical(parent['local_interacting_matter_history'])
            or canonical(bands)!=canonical(parent['centered_original_charged_bands'])
            or canonical(dressed)!=canonical(parent['source_generated_composite_CAR_response'])):
        raise ValueError('centered-matter owning source receipt differs from the replay')
    if (local['same_centered_thermodynamic_ground_subsequence']is not True
            or local['odd_matter_and_even_spin_share_full_CAR_action']is not True
            or local['original_Weyl_bare_repair_local_history_instantiated']is not True
            or bands['finite_and_infinite_ground_same_original_weight_bound']is not True):
        raise ValueError('centered-matter attachment lacks its common source premise')
    if any(local[name]is not False for name in(
            'effective_ground_subsequence_preparation_modulus','global_raw_archive_equality',
            'bare_repairs_identified_under_site_dependent_frame',
            'source_selected_capacity_or_repair_duration','relativistic_continuum_or_physical_particle_assignment')):
        raise ValueError('centered-matter attachment promotes its physical or archive scope')
    if bands['charged_pole_mass_or_Lorentz_covariance_derived']is not False:
        raise ValueError('centered-matter bands are not physical particle masses')
    if (dressed['independent_CAR_flavours_added']!=0 or dressed['physical_species_added']!=0
            or dressed['all_volume_first_composite_equation']is not True
            or dressed['all_volume_residual_discarded']is not False
            or dressed['dimer_or_aligned_reference_transferred_to_full_continuum']is not False):
        raise ValueError('centered-matter composite response exceeds its source or spectral scope')
    band=bands['cubic_readable_parameter_band']
    if (band['each_original_preparation_weight']!='1/2'
            or band['band_is_a_spectral_gap_or_particle_mass']is not False):
        raise ValueError('centered-matter original probability or band scope changed')
    return {'catalogue_entry':'centered_full_CAR_charge_response',
            'owner_receipt':'model/source_family_compatibility/results.json',
            'owning_local_and_charged_producers_replayed':True,
            'state':'particle-hole and hidden-spin invariant centered thermodynamic ground subsequence',
            'clock':'its positive centered GNS generator',
            'chart':local['chart'],
            'original_Weyl_bare_repair_instance':local['original_Weyl_bare_repair_local_instance'],
            'original_local_addition_and_removal_weights':['1/2','1/2'],
            'readable_source_charged_band':band,
            'complete_local_history_examples':local['local_history_plans'],
            'same_state_interaction_witness':local['interaction_example'],
            'source_generated_dimer_response':bands['actual_dimer'],
            'same_ground_composite_CAR_response':dressed,
            'retarded_response_error_is_arbitrary_input_physical_channel_bound':False,
            'global_raw_archive_equality':False,
            'global_sharp_energy_PVM_implemented_by_local_bound':False,
            'odd_transfer_native_service_derived':False,
            'effective_ground_subsequence_preparation_modulus':False,
            'topological_invariant_assigned':False,
            'ordered_spin_state_or_aligned_hole_reference_identity':False,
            'physical_species':[],'physical_masses':[]}


def build():
    pins = bindings()
    bare = defect_algebra.report()
    cover = {law:protection.algebra_certificate(law) for law in ('C','E')}
    bare_counts = {law:{domain:{key:bare['C_E'][law][domain][key] for key in (
        'dimension','generated_star_algebra_complex_dimension','commutant_complex_dimension')}
        for domain in ('full78','binary66')} for law in ('C','E')}
    cover_counts = {law:{key:record[key] for key in (
        'joint_dimension','generated_star_algebra_complex_dimension','commutant_complex_dimension',
        'complete_branch_X_commutation','ordered_cover_reachable_states','raw_maps_sha256')}
        for law,record in cover.items()}
    entries = rows()
    for entry in entries:
        entry['identified_Standard_Model_species'] = []
        entry['physical_mass'] = None
        entry['same_source_relativistic_particle_attachment_proved'] = False
    joint_attachment = joint_field_attachment()
    hole_attachment = hole_source_attachment()
    aligned_attachment = aligned_hole_attachment()
    centered_attachment = centered_matter_attachment()
    if bindings() != pins:
        raise RuntimeError('catalogue source changed during construction')
    return {'schema':'oph.protected_defect_catalog.v1', 'source_bindings':pins,
            'bare_complete_branch_classification':bare_counts,
            'ordered_cover_complete_branch_classification':cover_counts,
            'selected_spin_complete_Pauli_classification':pauli.census(),
            'same_source_charge_resolved_creation':charge.report(),
            'same_source_continuum_parity_attachment':field.preparation_bounds(),
            'same_source_continuum_configuration_charge':field.configuration_charge(),
            'observer_charge_and_creation_mechanism':mechanism.report(),
            'joint_positive_field_source_attachment':joint_attachment,
            'one_hole_source_attachment':hole_attachment,
            'aligned_hole_source_attachment':aligned_attachment,
            'centered_full_CAR_source_attachment':centered_attachment,
            'entries':entries,
            'catalogue_is_exhaustive_over_all_OPH_sources_or_defects':False,
            'full_parent_scientific_replays_performed':False,
            'parent_source_bindings_authenticated':True,
            'replayed_classifications':['bare all-seam branch commutant',
                                       'ordered-cover all-seam branch commutant',
                                       'actual selected-cube Pauli-product charges'],
            'physical_particle_map':[], 'physical_mass_prediction_count':0}


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('catalogue classification, scope or source custody changed')


def render(value):
    readouts = {
        'bare_native_source':'complete native records; no nontrivial source-only charge',
        'ordered_exchange_charge':'coherent charge pointer and pair-configuration transport',
        'selected_spin_charge':'coherent local charge creation; conjugate continuum configuration fields',
        'single_occupancy_spin_memory':'native spin controls and pointers; nontrivial exchange dynamics',
        'joint_positive_field_canonical_sector':'two regular canonical pairs and full internal U1; complete fixed-scale prepared H_B histories',
        'common_hidden_spin_factor':'complete ideal joint histories preserve the global spin factor',
        'cube_eta_range':'prepared eta reference and complete local spectral response',
        'compact_opposite_texture_cores':'declared transported four-mode band and current readout',
        'one_hole_spin_readout':'complete initial-hole position detector preserves internal spin/reference blocks; uniform global attempts give a 4/n archive bound relative to identical propagation with identity repair',
        'centered_full_CAR_charge_response':'complete local interacting records and positive original charged energy-band weight on one centered full-ground clock',
    }
    lines = ['# Protected sectors and defect candidates', '',
             'The catalogue records separate source architectures. It is not a list of physical species.',
             'Completeness applies only to the explicit domains in the table; there is no exhaustive census over all OPH realizations.', '',
             '| Object | Protecting source and object | Observer-facing connection | Physical species |',
             '| --- | --- | --- | --- |']
    for row in value['entries']:
        lines.append('| '+row['id']+' | '+row['protected_object']+' | '+readouts[row['id']]+' | Unidentified |')
    lines += ['', '## The source-to-observation chain', '',
              'The selected spin branch has the following proved connections:', '',
              '1. The signed native full-Fock charge pulls through the spin/Marshall inclusion to Q.',
              '2. Local Z followed by one X pointer prepares coherent Q=+1 or Q=-1 states, retaining both original probabilities.',
              '3. Complete [energy and clock readouts](../source_energy_effort/ENERGY_READOUT.md) measure their creation work.',
              '4. A [fixed charge detector](../charge_channel_lifetime/PROOF.md) retains nonzero contrast under charge-preserving evolution even when return survival decreases.',
              '5. On the prescribed ordered-state limit, the same U1 acts on Phi+=(Q3-iQ2) and Phi-=(Q3+iQ2). Source conjugation gives equal nonzero energy--momentum measures.',
              '6. [Bounded charged creators](../spin_charge_dispersion/PROOF.md) retain original positive joint-band probability at least1/1757600 at every momentum scale, using the same source state and clock.',
              '7. [Separated field preparations](../source_charge_scattering/CONTINUUM_PAIRS.md) realize every finite charge-sign pattern with positive original probability and an energy bound; all integer internal-charge sectors occur. These are sector labels, not a species count.', '',
              'The last connection is to collective field excitations. The finite local states are not identified with continuum poles.',
              'Dispersion, asymptotic scattering and physical quantum numbers are the required next connections to empirical species.',
              'The [Standard Model comparison table](../standard_model_source_map/PARTICLE_MAP.md) keeps that empirical assignment separate.', '']
    for row in value['entries']:
        lines += ['', '## '+row['id'].replace('_',' '), '',
                  '**Source:** '+row['source']+'.', '',
                  '**Classification domain:** '+row['completeness']+'.', '',
                  '**Preparation and readout:** '+row['preparation_readout']+'.', '',
                  '**Propagation and energy:** '+row['propagation_energy']+'.', '',
                  '**Topology:** '+row['topology']+'.', '',
                  '**Boundary:** '+row['boundary']+'.', '',
                  '[Owning proof or construction](../../'+row['evidence']+').']
    lines += ['', 'Source family, state, clock, physical spin/statistics, dispersion or resonance,',
              'detector coupling and experimental units must be joined before a catalogue',
              'entry can fill the [Standard Model correspondence](../standard_model_source_map/PARTICLE_MAP.md).', '']
    return '\n'.join(lines)
