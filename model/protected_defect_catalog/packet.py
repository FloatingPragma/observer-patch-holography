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
from model.high_field_particle import packet as high_field
from model.source_family_compatibility import hole_records, hole_dispersion, local_matter, charge_bands, composite, cube_response, charge_scaling, boost_current, cube_spectrum
from model.source_family_compatibility import momentum_response, charge_escape, uv_tail, local_spectral
from model.source_family_compatibility import kinetic_nodes, nodal_memory, ground_marginals, cube_memory
from model.source_family_compatibility import nodal_channels, spectral_residue, spectral_window_access
from model.source_family_compatibility import nodal_locality, nodal_moment_bounds, mass_scale_selection
from model.source_family_compatibility import source_action_identification, nodal_next_moments, native_lorentzian
from model.source_family_compatibility import bulk_action_identification
from . import pauli, charge, field, mechanism

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = tuple('model/'+name+'/results.json' for name in (
    'defect_protection', 'local_exchange_charge', 'common_spin_protection',
    'cube_eta_protection', 'compact_defect_pair', 'defect_localization_readout',
    'spin_canonical_field', 'spin_canonical_energy', 'spin_joint_field_limit',
    'source_family_compatibility', 'high_field_particle'))


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
         'evidence':'model/defect_localization_readout/PROOF.md',
         'spectral_attachments':['model/high_field_particle/PROOF.md']},
        {'id':'joint_positive_field_canonical_sector', 'kind':'constructive_field_family',
         'source':'protected spin interaction H_B=J sum(S_x.S_y-1/4)-B sum epsilon_x S_x^1; supplied p=n^-2 and B=J n^-5; positive-B ground and retained H_B clock',
         'protected_object':'two regular canonical pairs (Q3,P2) and (Q2,P3), with common positive-energy time and full internal U1 on their bounded histories',
         'completeness':'specified field subnet along a common history subsequence; neither all continuum observables nor all protected defects are classified',
         'preparation_readout':'finite positive-field ground preparation, local marginal modulus, source spin-rotation observables and complete native local Pauli histories with all pointer and failure records',
         'propagation_energy':'continuous even-envelope translations, positive generator K_B/(Jp), and finite propagation cone; bounded conjugate charge records each retain original probability at least 1/387200 in momentum [r/2,r] and energy [3r/387200,387200r], uniformly at small r; no sharp shell or rest mass',
         'topology':'internal source spin rotation, no spatial winding or identified electromagnetic charge',
         'boundary':'B(n) is a supplied preparation prescription; unread repair does not select it; equality to the ordered zero-field history limit, physical Lorentz spin, sharp poles and scattering are unproved',
         'evidence':'model/spin_joint_field_limit/PROOF.md',
         'same_object_realizations':['model/spin_joint_field_limit/SECOND_PAIR.md',
                                     'model/spin_joint_field_limit/CLOCK_COMPARISON.md',
                                     'model/spin_joint_field_limit/UNIFORM_MODE.md',
                                     'model/spin_joint_field_limit/ARCTAN_CHARGED_BANDS.md',
                                     'model/spin_joint_field_limit/DIRECTED_CHARGED_BANDS.md',
                                     'model/spin_joint_field_limit/QUANTUM_FORCE_SLACK.md']},
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
         'preparation_readout':'local addition and removal each have original weight 1/2; complete local phase and Lorentzian instruments retain all charge, time and complement records with boundary, marginal and service errors; one pressure diagnostic certifies a local score when all invariant grounds have a strict margin',
         'propagation_energy':'same positive centered ground clock, interacting local histories and positive original charged bands; source moments define every magnetic momentum; the actual preparation inversion and source symmetries reduce all 16 nodal components to one even scalar measure; composite memory and a positive fixed-source subcontact band are retained; locality makes positive nodal atom equivalent to bounded physical-window survival for packets supported within O(a) of a node on one fixed normalized reference, but actual atom positivity remains unproved',
         'topology':'source U1 charge grading, without an assigned winding invariant or Standard Model charge',
         'boundary':'the invariant ground subsequence has no effective preparation modulus; on its specified cubic average at U/c=8, the full-CAR paired-star bound gives liminf escaping weight at least 64/(95+sqrt(577)) above every fixed energy window for every normalized original packet under ballistic scaling; this excludes the entire unfiltered regular unit-mass limit but not filtered branches; any full fixed-reference original-field ballistic limit has zero weight on nonzero relativistic mass shells; convergent geometric refinement only enforces a discrete shell tower if a positive-mass shell exists; no actual field limit, physical spin or measured mass is established',
         'evidence':'model/source_family_compatibility/CHARGED_SPECTRAL_BANDS.md',
         'same_object_realizations':['model/source_family_compatibility/LOCAL_MATTER_READOUT.md',
                                     'model/source_family_compatibility/COMPOSITE_CAR.md',
                                     'model/source_family_compatibility/CUBE_RESPONSE.md',
                                     'model/source_family_compatibility/STRONG_CHARGE_SCALING.md',
                                     'model/source_family_compatibility/ENERGY_DIPOLE.md',
                                     'model/source_family_compatibility/CUBE_SPECTRUM.md',
                                     'model/source_family_compatibility/MOMENTUM_RESPONSE.md',
                                     'model/source_family_compatibility/LOCAL_CHARGE_ESCAPE.md',
                                     'model/source_family_compatibility/UV_TAIL.md',
                                     'model/source_family_compatibility/LOCAL_SPECTRAL_READOUT.md',
                                     'model/source_family_compatibility/NODAL_MEMORY.md',
                                     'model/source_family_compatibility/GROUND_MARGINALS.md',
                                     'model/source_family_compatibility/CUBE_MEMORY.md',
                                     'model/source_family_compatibility/NODAL_CHANNELS.md',
                                     'model/source_family_compatibility/SPECTRAL_RESIDUE.md',
                                     'model/source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md',
                                     'model/source_family_compatibility/NODAL_LOCALITY.md',
                                     'model/source_family_compatibility/NODAL_MOMENT_BOUNDS.md',
                                     'model/source_family_compatibility/MASS_SCALE_SELECTION.md',
                                     'model/source_family_compatibility/NODAL_NEXT_MOMENTS.md',
                                     'model/source_family_compatibility/NATIVE_LORENTZIAN_READOUT.md',
                                     'model/source_family_compatibility/BULK_ACTION_IDENTIFICATION.md',
                                     'model/source_family_compatibility/SOURCE_ACTION_IDENTIFICATION.md']},
        {'id':'kinetic_Weyl_nodes_and_free_fields', 'kind':'kinetic_topology_and_conditional_field_map',
         'source':'original cubic Weyl symbol -c sum sigma_j sin(k_j); centered full-ground family c=v/a and 0<=U=o(c^(2/7))',
         'protected_object':'eight isolated kinetic nodes, four of winding +1 and four of winding -1; conditional free vacuum map to eight Weyl fields',
         'completeness':'all zeros of the stated one-body symbol on its momentum torus; not all interacting excitations or repair-protected charges',
         'preparation_readout':'the same centered free-limit ground hierarchy, spatial wavepackets and positive clock; exact source spin-frame and Pauli-band projector dictionary',
         'propagation_energy':'free positive particle and antiparticle energies v|p|, zero rest mass and opposite helicities on the stated weak-contact family; can be paired into four massless Dirac fields',
         'topology':'momentum-space kinetic winding with explicit boundary-gap stability; period-two onsite mass channels couple opposite winding sectors',
         'boundary':'the default U/c=8 family is outside this free theorem; neither a selected mass coefficient nor an interacting pole nor a measured species is inferred from winding',
         'evidence':'model/source_family_compatibility/KINETIC_NODE_MAP.md'},
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
    quantum_force = value['source_quantum_correction_to_mean_spin_force']
    if (joint_field.canonical(quantum_force) != joint_field.canonical(
            joint_field.serialize(joint_field.quantum_force_slack.certificate()))
            or scope['strict_source_bond_minus_squared_order_margin_proved'] is not True
            or scope['unfiltered_microscopic_mean_spin_force_factorization_excluded'] is not True
            or scope['quantum_force_margin_excludes_corrected_infrared_shell'] is not False):
        raise ValueError('joint-field quantum force correction differs from its source scope')
    history = value['complete_prepared_local_history']
    uniform=value['uniform_mode_original_weight_and_concentration']
    spectral=value['complete_finite_uniform_spectral_readout']
    charged=value['bounded_joint_charged_energy_momentum_bands']
    directed=value['directed_joint_charged_bands']
    cap=directed['directed_joint_source_band']
    if (cap['each_charge_and_direction_original_probability_lower']!='1/774400'
            or cap['source_even_collective_internal_charge']is not True
            or cap['original_odd_CAR_identification']is not False
            or cap['energy_shell_or_plane_wave_atom_inferred']is not False
            or cap['effective_finite_source_to_continuum_error_claimed']is not False
            or any(directed['finite_periodic_source_control'][name]is not False for name in (
                'complete_native_pulse_tape_executed','finite_control_proves_continuum_floor',
                'finite_translation_labels_identified_with_continuum_momentum',
                'physical_particle_or_mass_assignment'))):
        raise ValueError('joint-field directed bands exceed their source scope')
    if (charged['continuum_band']['each_original_charge_probability_lower']!='1/387200'
            or charged['continuum_band']['joint_positive_field_clock_used']is not True
            or charged['continuum_band']['sharp_shell_or_physical_mass_inferred']is not False
            or charged['continuum_band']['old_zero_field_clock_substituted']is not False
            or charged['continuum_band']['original_odd_CAR_identification']is not False
            or charged['full_history_linear_force_closure_proved']is not False
            or charged['continuum_limit_error_rate_supplied']is not False
            or charged['physical_species']!=[] or charged['physical_mass_prediction']is not None):
        raise ValueError('joint-field charged response exceeds its owning scope')
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
            'bounded_charged_joint_bands':charged,
            'directed_joint_charged_bands':directed,
            'source_quantum_correction_to_mean_spin_force':quantum_force,
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
    cube=serialize(cube_response.certificate())
    scaling=serialize(charge_scaling.certificate())
    current=serialize(boost_current.certificate())
    spectrum=serialize(cube_spectrum.certificate())
    local_energy={
        'actual_cubic_momentum_response':serialize(momentum_response.certificate()),
        'default_original_charge_escape':serialize(charge_escape.certificate()),
        'uniform_original_energy_tails':serialize(uv_tail.certificate()),
        'complete_local_source_energy_readout':serialize(local_spectral.certificate()),
        'actual_nodal_response_memory':serialize(nodal_memory.certificate()),
        'ground_marginal_preparation_and_native_records':ground_marginals.certificate(),
        'full_cube_frequency_dependent_memory':serialize(cube_memory.certificate()),
        'actual_nodal_scalar_channels':serialize(nodal_channels.certificate()),
        'actual_spectral_residue_and_subcontact_band':serialize(spectral_residue.certificate()),
        'complete_Lorentzian_window_access':serialize(spectral_window_access.certificate()),
        'actual_nodal_locality_and_window_survival':serialize(nodal_locality.certificate()),
        'uniform_full_CAR_star_force_bound':serialize(nodal_moment_bounds.certificate()),
        'fixed_source_scaling_and_mass_shell_scope':serialize(mass_scale_selection.certificate()),
        'complete_source_action_identification':serialize(source_action_identification.certificate()),
        'actual_local_nodal_higher_moments':serialize(nodal_next_moments.certificate()),
        'executed_native_dimer_Lorentzian_readout':serialize(native_lorentzian.certificate()),
        'volume_uniform_cubic_action_identification':serialize(bulk_action_identification.certificate())}
    if any(canonical(value)!=canonical(parent[key]) for key,value in local_energy.items()):
        raise ValueError('centered-matter local energy or momentum receipt differs from the replay')
    if (local_energy['actual_cubic_momentum_response']['continuum_poles_or_physical_species_derived']is not False
            or local_energy['uniform_original_energy_tails']['positive_shrinking_window_weight_proved']is not False
            or local_energy['complete_local_source_energy_readout']['physical_mass_prediction']is not None
            or local_energy['default_original_charge_escape']['physical_particle_species']!=[]
            or local_energy['default_original_charge_escape']['physical_particle_masses']!=[]
            or local_energy['actual_nodal_response_memory']['physical_species_identified']!=[]
            or local_energy['actual_nodal_response_memory']['default_source_result']['zero_frequency_residue_inferred']is not False):
        raise ValueError('centered-matter energy readout promotes its physical scope')
    ground=local_energy['ground_marginal_preparation_and_native_records']
    if any(ground[name]is not False for name in (
            'thermodynamic_local_marginal_convergence_modulus_proved',
            'cut_ground_substituted_for_infinite_ground','physical_state_selection_proved')):
        raise ValueError('centered-matter finite preparation promotes its infinite-state scope')
    memory=local_energy['full_cube_frequency_dependent_memory']
    if (memory['physical_mass_assigned']is not None
            or memory['thermodynamic_charge_gap_proved']is not False
            or memory['native_large_filter_tape_executed']is not False):
        raise ValueError('centered-matter cube memory exceeds its finite-source scope')
    channels=local_energy['actual_nodal_scalar_channels']
    if any(channels[name]is not False for name in (
            'uniform_growing_packet_scaling_supplied',
            'physical_spin_or_species_identification','physical_mass_prediction')):
        raise ValueError('centered-matter nodal channels exceed their spectral scope')
    residue=local_energy['actual_spectral_residue_and_subcontact_band']
    if (residue['physical_mass_prediction']is not None or residue['physical_species']!=[]
            or residue['default_shrinking_window_weight_proved']is not False
            or residue['uniform_many_body_charge_gap_proved']is not False):
        raise ValueError('centered-matter spectral residue exceeds its proved source scope')
    access=local_energy['complete_Lorentzian_window_access']
    if (access['physical_species_assigned']is not False or access['physical_mass_assigned']is not None
            or access['uniform_shrinking_window_weight_proved']is not False
            or access['effective_growing_ground_marginal_modulus_proved']is not False):
        raise ValueError('centered-matter spectral access exceeds its preparation scope')
    locality=local_energy['actual_nodal_locality_and_window_survival']
    plan=locality['conditional_example']
    if (locality['physical_mass_assigned']is not None
            or any(locality[name]is not False for name in (
                'actual_nodal_atom_positivity_proved','sharp_dispersion_proved',
                'growing_native_instrument_executed'))
            or any(plan[name]is not False for name in (
                'actual_nodal_atom_lower_proved','native_growing_readout_error_met',
                'positive_energy_annulus_proved','dispersion_or_particle_identification_proved'))
            or plan['same_normalized_reference_at_every_spacing_required']is not True
            or plan['zero_energy_included']is not True):
        raise ValueError('centered-matter locality promotes its conditional window scope')
    star=local_energy['uniform_full_CAR_star_force_bound']
    budget=star['source_moment_budget']
    if (star['physical_species_identified']!=[] or star['physical_mass_prediction']is not None
            or budget['physical_species_identified']!=[] or budget['physical_mass_prediction']is not None
            or any(star[name]is not False for name in (
                'spin_only_ground_or_antiferromagnetic_order_assumed',
                'universal_local_bound_attained_in_the_actual_infinite_ground',
                'native_large_tape_executed'))
            or any(budget[name]is not False for name in (
                'actual_zero_atom_positive_proved','source_eighth_moment_computed'))):
        raise ValueError('centered-matter star bounds promote their one-sided source scope')
    dilation=local_energy['fixed_source_scaling_and_mass_shell_scope']
    if (dilation['physical_species_assigned']!=0 or dilation['physical_mass_assigned']is not None
            or any(dilation[name]is not True for name in (
                'fixed_reference_and_dimensionless_couplings_required',
                'full_joint_vague_limit_required'))
            or any(dilation[name]is not False for name in (
                'arbitrary_sparse_subsequence_excluded_by_theorem',
                'geometric_limit_has_continuous_dilation_covariance',
                'every_unstable_resonance_excluded','continuous_timelike_response_excluded',
                'actual_full_continuum_limit_proved','actual_nodal_atom_positivity_proved',
                'running_reference_or_observable_dictionary_excluded'))
            or dilation['actual_default_source_homogeneity']['native_large_detector_tape_executed']is not False
            or dilation['actual_default_source_homogeneity']['finite_dimer_is_continuum_momentum_mode']is not False
            or dilation['actual_default_source_homogeneity']['conditional_probability_renormalization_applied']is not False):
        raise ValueError('centered-matter dilation exceeds its fixed-source conditional scope')
    identification=local_energy['complete_source_action_identification']
    if (identification['source']['complete_record_count']!=64
            or identification['source']['postselection_used']is not False
            or identification['source_clock_and_known_initial_preparation_required']is not True
            or any(identification[name]is not False for name in (
                'coherent_clock_pulse_tape_executed_here',
                'physical_strong_coupling_or_QCD_alpha_selected',
                'capacity_selects_action_or_physical_clock',
                'thermodynamic_ground_selected_by_this_readout',
                'physical_particle_or_mass_identified'))):
        raise ValueError('source action identification exceeds its finite observation scope')
    moments=local_energy['actual_local_nodal_higher_moments']
    if (moments['physical_species_identified']!=[] or moments['physical_masses']!=[]
            or any(moments[name]is not False for name in (
                'source_selected_ground_moments_evaluated','actual_zero_atom_decided'))
            or moments['actual_source_memory_bound']['zero_frequency_atom_inferred']is not False
            or any(row['selected_ground_moment']is not None
                   for row in moments['actual_local_nodal_moment_readouts'])):
        raise ValueError('local moment observables exceed their selected-ground scope')
    native=local_energy['executed_native_dimer_Lorentzian_readout']
    if (native['source_sites']!=2 or native['full_CAR_rows']!=16
            or native['proper_terminal_record_count']!=8194
            or native['raw_pointer_terminal_record_count']!=36*native['time_count']+6
            or native['physical_first_failure_records']!=17449226
            or any(native[name]is not True for name in (
                'same_centered_source_excitation_clock',
                'all_native_pulses_on_potentially_nonzero_prefixes_executed',
                'all_source_rows_and_charge_columns_retained','dummy_clock_environment_retained',
                'zero_charge_prefixes_follow_same_declared_controller','selector_blanks_are_supplied'))
            or any(native[name]is not False for name in (
                'floating_contraction_diagnostics_used_as_interval_premise',
                'postselection_or_output_renormalization_used',
                'large_infinite_ground_collar_executed','hardware_shots_sampled',
                'physical_species_or_mass_assigned'))):
        raise ValueError('native energy experiment exceeds its complete finite-source scope')
    bulk=local_energy['volume_uniform_cubic_action_identification']
    if (bulk['known_local_calibration_preparation_consumed']is not True
            or bulk['physical_particles']!=[] or bulk['physical_masses']!=[]
            or any(bulk[name]is not False for name in (
                'complete_history_or_reference_error_from_marginal_inferred',
                'interacting_ground_preparation_inferred','capacity_selects_action_or_clock',
                'large_native_experiment_executed'))
            or any(row['is_probability_of_the_larger_collar_plan']is not False
                   for row in bulk['finite_time_ratio_witnesses'])
            or bulk['separate_arbitrary_exterior_collar_plan']['pure_two_CAR_collar_input_required']is not True
            or bulk['separate_arbitrary_exterior_collar_plan']['reset_pulse_tape_executed']is not False
            or bulk['complete_collision_detector']['full_fine_instrument_covariance_asserted']is not False
            or bulk['uniform_global_attempt_example']['uniform_seam_selector_required']is not True
            or bulk['uniform_global_attempt_example']['other_seam_records_discarded_or_postselected']is not False):
        raise ValueError('cubic action inverse exceeds its local calibration scope')
    if (canonical(local)!=canonical(parent['local_interacting_matter_history'])
            or canonical(bands)!=canonical(parent['centered_original_charged_bands'])
            or canonical(dressed)!=canonical(parent['source_generated_composite_CAR_response'])
            or canonical(cube)!=canonical(parent['default_cube_composite_response'])
            or canonical(scaling)!=canonical(parent['volume_uniform_charged_scaling'])
            or canonical(current)!=canonical(parent['full_source_energy_dipole'])
            or canonical(spectrum)!=canonical(parent['complete_cube_wavepacket_spectrum'])):
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
            'cubic_momentum_and_escape_state':'the explicit proper-cubic arithmetic average of that ground, used at every momentum without conditioning',
            'clock':'its positive centered GNS generator',
            'chart':local['chart'],
            'original_Weyl_bare_repair_instance':local['original_Weyl_bare_repair_local_instance'],
            'original_local_addition_and_removal_weights':['1/2','1/2'],
            'readable_source_charged_band':band,
            'complete_local_history_examples':local['local_history_plans'],
            'same_state_interaction_witness':local['interaction_example'],
            'source_generated_dimer_response':bands['actual_dimer'],
            'same_ground_composite_CAR_response':dressed,
            'finite_default_cube_composite_response':cube,
            'volume_uniform_charged_scaling':scaling,
            'full_source_energy_dipole':current,
            'finite_complete_cube_wavepacket_spectrum':spectrum,
            **local_energy,
            'retarded_response_error_is_arbitrary_input_physical_channel_bound':False,
            'global_raw_archive_equality':False,
            'global_sharp_energy_PVM_implemented_by_local_bound':False,
            'odd_transfer_native_service_derived':False,
            'effective_ground_subsequence_preparation_modulus':False,
            'topological_invariant_assigned':False,
            'ordered_spin_state_or_aligned_hole_reference_identity':False,
            'physical_species':[],'physical_masses':[]}


def kinetic_field_attachment():
    parent=strict_json.load(ROOT/'model/source_family_compatibility/results.json')
    value=serialize(kinetic_nodes.certificate())
    if canonical(value)!=canonical(parent['kinetic_node_to_free_field_map']):
        raise ValueError('kinetic-field owning source receipt differs')
    if (value['scope']['default_interacting_Green_function_winding_derived']is not False
            or value['scope']['physical_Standard_Model_assignment']is not False
            or value['scope']['all_native_repair_protected_charge']is not False
            or value['scope']['real_space_defect_classification']is not False
            or value['free_field_scope']['default_U_over_c_eight_in_this_class']is not False
            or any(row['physical_particle_assignment']is not None for row in value['nodes'])):
        raise ValueError('kinetic-field map exceeds its physical scope')
    return value


def high_field_attachment():
    """Attach a source band without introducing another protected object."""
    value = high_field.build()
    if canonical(value) != canonical(strict_json.load(high_field.OUTPUT)):
        raise ValueError('high-field owning source receipt differs')
    checks = {
        'exact_source_algebra': high_field.algebra.certificate(),
        'complete_local_charge_creation': high_field.readout.creation_certificate(),
        'finite_ground_preparation_plan': high_field.readout.finite_preparation_example(),
        'short_complete_native_preparation': high_field.preparation.execution_certificate(),
    }
    for key, expected in checks.items():
        if canonical(value[key]) != high_field.canonical(expected):
            raise ValueError('high-field independent source certificate differs: '+key)
    scope = value['scope']
    required = ('exact_finite_coefficients_through_second_order',
                'complete_local_charge_instrument_certified',
                'same_source_finite_ground_planner_reused',
                'native_preparation_tape_executed_in_this_packet',
                'bounded_J_over_B_finite_mass_local_relativistic_refinement_excluded',
                'exclusion_uses_same_ground_relative_clock_and_spatial_dictionary')
    excluded = ('finite_coupling_dispersion_remainder_certified',
                'all_source_regimes_excluded', 'measured_mass_inputs_used')
    contract = value['external_theorem_contract']
    if (any(scope[key] is not True for key in required)
            or any(scope[key] is not False for key in excluded)
            or scope['physical_particle_assignments'] != []
            or scope['physical_mass_predictions'] != []
            or scope['uniform_numeric_original_particle_residue_floor'] is not None
            or scope['executed_preparation_scope'] != 'finite eight-site J=1 B=400 product-residual certificate'
            or contract['band_existence'] != 'for sufficiently small positive J/B'
            or contract['explicit_high_field_admission_threshold'] is not None
            or contract['theorems_reproved_by_finite_matrix_tests'] is not False):
        raise ValueError('high-field spectral attachment exceeds its source scope')
    expected_speed = [
        {'spacing': str(a), 'J_over_B_ceiling': '1/100',
         'physical_rest_energy_ceiling': '1',
         'physical_speed_upper': str(high_field.readout.speed_bound(a, F(1, 100), F(1)))}
        for a in (F(1, 8), F(1, 64), F(1, 512), F(1, 4096))]
    if canonical(value['bounded_mass_speed_controls']) != canonical(expected_speed):
        raise ValueError('high-field same-clock speed bound differs')
    return {
        'catalogue_entry': 'single_occupancy_spin_memory',
        'attachment_kind': 'spectral band on a specified high-field source family',
        'additional_protected_objects_counted': 0,
        'additional_topological_classes_counted': 0,
        'owning_scientific_and_native_producer_replayed': True,
        'state_and_clock_regime': 'positive-field ground and its H_B clock, sufficiently small positive J/B',
        'identity_with_ordered_or_joint_small_field_limit_claimed': False,
        **{key: value[key] for key in checks},
        'external_theorem_contract': contract,
        'bounded_mass_speed_controls': value['bounded_mass_speed_controls'],
        'owning_scope': scope,
        'physical_species': [], 'physical_masses': [],
    }


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
    kinetic_attachment = kinetic_field_attachment()
    high_field_spectral = high_field_attachment()
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
            'kinetic_node_free_field_attachment':kinetic_attachment,
            'high_field_spin_spectral_attachment':high_field_spectral,
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
        'single_occupancy_spin_memory':'native spin controls and pointers; separate high-field family has a conditional isolated lattice band',
        'joint_positive_field_canonical_sector':'two regular canonical pairs and full internal U1; complete fixed-scale prepared H_B histories',
        'common_hidden_spin_factor':'complete ideal joint histories preserve the global spin factor',
        'cube_eta_range':'prepared eta reference and complete local spectral response',
        'compact_opposite_texture_cores':'declared transported four-mode band and current readout',
        'one_hole_spin_readout':'complete initial-hole position detector preserves internal spin/reference blocks; uniform global attempts give a 4/n archive bound relative to identical propagation with identity repair',
        'centered_full_CAR_charge_response':'complete local interacting records and positive original charged energy-band weight on one centered full-ground clock',
        'kinetic_Weyl_nodes_and_free_fields':'exact winding, chirality and particle/hole dictionary for the proved weak-contact free vacuum; no Standard Model assignment',
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
    lines += ['', '## High-field spin spectral attachment', '',
              'The [high-field source](../high_field_particle/PROOF.md) is a spectral attachment to the single-occupancy spin code, not an additional protected defect or topological class.',
              'For sufficiently small positive J/B, its positive-field ground has an isolated lattice band with conjugate internal charges -1 and +1. No numerical admission threshold is supplied.',
              'The complete local creator is CAR-even and retains failure; it does not identify an original odd physical fermion. Its finite-coupling success is not a pure band projection.',
              'The separate eight-site J=1,B=400 preparation executes 167 native pulses with 608 failure slots. Its exact-ground complete-output error is below 1/100; its ideal pointer success floor 155233/155236 is not an infinite-volume particle residue.',
              'On the same ground-relative clock and local spatial dictionary, bounded J/B and bounded physical rest energy force continuum propagation speed to zero. The ordered B/J->0 family is outside this result. No physical species or mass is assigned.', '',
              '## Exact kinetic-node to free-field map', '',
              'This table uses the centered weak-contact family c=v/a, 0<=U=o(c^(2/7)).',
              'Its free energies are v|p|. It assigns no measured particle; U/c=8 is outside this free theorem.', '',
              '| Source node / pi | Kinetic winding | Free particle helicity | Antiparticle helicity | Free rest mass |',
              '| --- | ---: | ---: | ---: | ---: |']
    for row in value['kinetic_node_free_field_attachment']['nodes']:
        lines.append('| '+str(tuple(row['source_node_pi_units']))+' | '+str(row['kinetic_winding'])
                     +' | '+row['positive_charge_free_helicity']+' | '+row['conjugate_free_helicity']+' | 0 |')
    lines += ['', 'Source family, state, clock, physical spin/statistics, dispersion or resonance,',
              'detector coupling and experimental units must be joined before a catalogue',
              'entry can fill the [Standard Model correspondence](../standard_model_source_map/PARTICLE_MAP.md).', '']
    return '\n'.join(lines)
