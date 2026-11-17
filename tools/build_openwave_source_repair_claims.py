#!/usr/bin/env python3
"""Bind finite native selection/response results to their portable evidence."""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_source_repairs'
OUTPUT = ROOT / 'claims/openwave_source_repairs.json'
sys.set_int_max_str_digits(0)  # Hash-bound exact rational certificates exceed 4300 digits.


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def generate():
    spec = importlib.util.spec_from_file_location('source_repair_custody', PACKAGE / 'replay.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    file_counts = module.custody()
    origins = read(PACKAGE / 'ORIGINS.json')
    charge = read(PACKAGE / 'charge/RECEIPT.json')
    rotation = read(PACKAGE / 'rotation/RECEIPT.json')
    coupling = read(PACKAGE / 'rotation/COUPLING_RECEIPT.json')
    bridge = read(PACKAGE / 'bridge/RECEIPT.json')
    contact = read(PACKAGE / 'contact/RECEIPT.json')
    scalar = read(PACKAGE / 'scalar/RECEIPT.json')
    clock = read(PACKAGE / 'clock/RECEIPT.json')
    apolar_join = read(PACKAGE / 'apolar_join/TRANSFER_RESULT.json')
    pair_loader = read(PACKAGE / 'pair_loader/RECEIPT.json')
    pair_route = read(PACKAGE / 'pair_route/RECEIPT.json')
    flow = read(PACKAGE / 'flow/workspace/openwave/shared/compact_logarithmic_flow/attempts/stage2_002/RECEIPT.json')['scientific']
    flow_replay = read(PACKAGE / 'reviews/flow/FULL_REPLAY_COMPARISON.json')
    energy_dir = PACKAGE / 'energy/workspace/openwave/row_02_electron_rest_energy/work/2026-11-15_codex_source_repair_v2'
    energy = read(energy_dir / 'RECEIPT.json')
    energy_routed = read(energy_dir / 'routed_join/RECEIPT.json')
    energy_routed = energy_routed.get('scientific', energy_routed)
    common = read(PACKAGE / 'common/RECEIPT.json')
    unknown = read(PACKAGE / 'unknown_total/RECEIPT.json')
    aligned = read(PACKAGE / 'aligned_pair/RECEIPT.json')
    local = read(PACKAGE / 'source_local_flow/workspace/source_local_flow/RECEIPT.json')
    common_v2 = read(PACKAGE / 'common_v2/RECEIPT.json')
    require(charge['all_gates_pass'] and all(charge['gates'].values()), 'charge gates')
    require(charge['literal_native_calls'] == 555001 and
            charge['local_endpoint_quadruples'] == 6561 and
            charge['primitive_source_branches'] == 7920, 'charge scope')
    require(not charge['physical_electric_identification'] and
            not charge['full_row_solved'] and not charge['g6_confirmed'], 'charge promotion')
    require(rotation['source_antipodal_identification'] and
            rotation['source_axis_preserved'] and
            rotation['signed_ordered_discovery_contexts'] == 120, 'rotation mechanism')
    require(coupling['conditional_on_each_coin'] and coupling['cases'] == 120 and
            coupling['future_source_processes_equal'] and
            coupling['historical_current_records_retained_distinct'], 'rotation coupling')
    require(not rotation['full_physical_spin_row_solved'] and
            not rotation['g6_other_group_confirmed'], 'rotation promotion')
    require(bridge['all_gates_pass'] and all(bridge['gates'].values()) and
            bridge['source240_actual_not_difference'] and
            bridge['literal_native_calls'] == 438693 and
            not bridge['coherent_state_selection_claimed'], 'actual-source bridge')
    require(contact['all_pass'] and all(contact['checks'].values()) and
            contact['contact']['uniform_total'] == '1/1260' and
            contact['outcome_totals']['native']['psi_click'] == '932/945' and
            not contact['scope']['radiation'], 'contact evidence')
    require(scalar['all_pass'] and not scalar['physical_spin_generated'] and
            scalar['unread_ray_fidelity'] == '11867/12144' and
            scalar['unread_escape_probability'] == '1/506', 'scalar scope')
    require(clock['all_gates_pass'] and all(clock['gates'].values()) and
            clock['native_single_generator_values'] == [6, 6, 6, 6] and
            clock['native_pair_generator_values'] == [12, 12, 12, 12] and
            clock['literal_native_calls'] == 242175 and
            clock['fixed_one_attempt_factor_gap'] == '1/44100', 'native clock source')
    require(not clock['physical_rest_energy_identified'] and
            not clock['coherent_UL_execution_in_this_packet'] and
            not clock['full_row_solved'], 'clock promotion')
    require(apolar_join['all_certified_stage_gates'] and
            apolar_join['native_implementation_transfer_complete'] and
            not apolar_join['failure_budget_deleted'] and
            not apolar_join['full_row_award_claimed'], 'same-object spatial join')
    require(pair_loader['all_gates_pass'] and all(pair_loader['gates'].values()) and
            pair_loader['literal_native_calls'] == 472146 and
            not pair_loader['complete_quantum_sector_AND_implemented_here'] and
            not pair_loader['physical_energy_or_spin_claimed'], 'actual signed-pair loading')
    require(pair_route['all_gates_pass'] and all(pair_route['gates'].values()) and
            pair_route['actual_output_pair_code'] == 14 and not pair_route['reset_used'] and
            pair_route['literal_native_calls'] == 150, 'actual common routed source')
    require(flow['native_processor_counts']['pulse_applications'] == 8725933 and
            flow['native_processor_counts']['native_matrix_executions'] == 196316 and
            flow['source_native_columns_from_registered_proposal'] == 302400 and
            not flow['physical_clock_calibration_claimed'] and not flow['row3_self_start_claimed'], 'full source-flow implementation')
    require(flow_replay['full_producer_exit'] == 0 and flow_replay['scientific_json_byte_equal'] and
            len(flow_replay['artifacts']) == 10 and
            all(a['byte_equal'] for a in flow_replay['artifacts'].values()), 'full isolated flow replay')
    require(energy['final_scientific_gates_pass'] and not energy['historical_raw_trace_gate_pass'] and
            not energy['other_group_G6'] and not energy['submission_ready'] and
            not energy_routed['full_row_solved'] and len(energy_routed['handoffs']) == 12,
            'selected source energy scope')
    require(Fraction(energy_routed['complete_join_error_upper']) < Fraction(1, 1000) and
            energy_routed['pair_to_single_interval'][0] < 2 < energy_routed['pair_to_single_interval'][1] and
            energy_routed['pair_keV_at511_anchor_interval'][1] < 1024.45,
            'complete source energy and clock discrimination')
    require(common['all_gates_pass'] and all(common['gates'].values()) and
            common['broad_loader_before_assay_choice'] and
            common['spatial_assay_before_dark_sector_diagnosis'] and
            common['fresh_copies_same_preparation'] and
            common['field_energy_density_max_difference'] < 1e-12,
            'same selected source and complete loading before both assays')
    require(Fraction(common['energy_complete_error_upper']) < Fraction(1, 1000) and
            Fraction(common['spatial_halfturn_TV_upper']) < Fraction(1, 1000) and
            Fraction(common['active_dephasing_TV_lower']) > Fraction(35, 1000) and
            common['pair_to_single_interval'][0] < 2 < common['pair_to_single_interval'][1] and
            common['pair_keV_at_single511_anchor_interval'][1] < 1024.45,
            'common instrument error and substantive controls')
    require(common['literal_new_native_calls'] == 0 and
            common['fresh_fullflow_executions'] == 0 and
            common['compressed_native_routine_applications'] == 68 and
            not common['other_group_G6'] and
            not common['upstream_award_or_pin_changed'] and
            not common['finite_exact_halfturn_claimed'], 'common execution and award scope')
    require(unknown['all_gates_pass'] and all(unknown['gates'].values()) and
            unknown['literal_native_calls'] == 5576 and
            unknown['complete_outcomes'] == dict(handoff=24, no_unit=1, selection_rejected=25) and
            unknown['controller_seed_independent_of_K'] and
            unknown['all_K_labels_hidden_from_controller'] and
            not unknown['total_one_promise_used'] and not unknown['reset_used'],
            'unknown-total current-only preparation')
    require(Fraction(unknown['preparation_inconclusive_upper']) == Fraction(8501, 2**58) and
            not unknown['independent_topology_established'] and
            not unknown['full_row_solved'] and not unknown['other_group_G6'],
            'unknown-total theorem and scope')
    require(aligned['all_gates_pass'] and all(aligned['gates'].values()) and
            aligned['full144_density_not_marginal_target'] and
            aligned['geometry_alignment_target_informed_and_supplied'] and
            aligned['terminal_native_bit_readers'] == 16 and
            aligned['joint_address_outcomes'] == 65536 and
            Fraction(aligned['unaligned_joint_ideal_TV']) == Fraction(1, 3),
            'whole-object spatial outcome and retained counterexample')
    require(Fraction(aligned['full_joint_halfturn_TV_upper']) < Fraction(1, 1000) and
            Fraction(aligned['active_dephasing_full_joint_TV_lower']) > Fraction(35, 1000) and
            aligned['literal_new_native_calls'] == 0 and
            not aligned['finite_exact_halfturn_claimed'] and not aligned['other_group_G6'],
            'whole-object complete bounds and execution scope')
    require(local['all_gates_pass'] and local['bounded_replay_read_only'] and
            local['gate_exit_tests'] == 23 and not local['compact_projector_or_energy_scalar_in_compiler'] and
            local['compiler_inputs'] == ['edge list', 'public time', 'refinement'] and
            local['shared_clock_slots_per_time'] == 154395440640 and
            local['clock_padding_source_channel_error'] == 0 and
            not local['other_group_G6'] and not local['submission_ready'],
            'source-local compiler, common clock and scope')
    require(local['final_successful_execution_counts']['pinned_endpoint_calls'] == 384324 and
            local['final_successful_execution_counts']['coupled_pair_truth_columns'] == 72273600 and
            local['final_successful_execution_counts']['M6_matrices'] == 22 and
            local['final_successful_execution_counts']['clock_zero_source_calls'] == 258,
            'local source literal versus represented counts')
    require(common_v2['all_gates_pass'] and all(common_v2['gates'].values()) and
            common_v2['source_join']['actual_handoffs'] == 24 and
            common_v2['single_calibration_join']['actual_positive_only_calibrations'] == 36 and
            common_v2['single_calibration_join']['independent_unit_chart_not_signed_pair_spectator'] and
            common_v2['full_selected_energy_spatial_density_max_difference'] < 1e-14 and
            common_v2['native_energy_pointer_max_difference'] < 1e-14,
            'actual common pair and separate physical calibration')
    require(all(Fraction(v) < Fraction(1, 1000) for v in common_v2['energy_complete_errors_upper'].values()) and
            Fraction(common_v2['full_joint_halfturn_TV_upper']) < Fraction(1, 1000) and
            Fraction(common_v2['active_dephasing_full_joint_TV_lower']) > Fraction(35, 1000) and
            common_v2['pair_keV_at511_anchor'][0] < 1022 < common_v2['pair_keV_at511_anchor'][1] < 1024.45 and
            common_v2['source_law_fringe_distance_centre'] > common_v2['source_law_two_arm_error_radius_upper'] and
            common_v2['literal_new_native_calls'] == common_v2['new_fullflow_executions'] == 0 and
            not common_v2['other_group_G6'] and not common_v2['finite_exact_halfturn_claimed'],
            'common complete error, source-sensitive observer and execution scope')
    claims = [
        dict(id='native_arbitrary_background_response', rows=[1], packet='charge',
             statement='Two integer sources under common native choices have exactly the native discrepancy/current law. Equal-total sources synchronize with finite mean time; compared stopped currents have mean B^T L^+ q. A selected relative unit drives this response.',
             domain='Finite connected graph; integer inputs; common independent uniform edges and fair coins; declared twin preparation and record-only probes. The theorem is classical.',
             falsifier='A mismatch in the complete conditional discrepancy/current law, a nonabsorbing neutral class, or failure of the exact current-reward identity.',
             measured=dict(native_calls=555001, local_quadruples=6561, exact_primitive_branches=7920, successful_native_pipelines=48)),
        dict(id='native_apolar_field_and_oriented_frame', rows=[7], packet='rotation',
             statement='A first-current-discovered bond couples opposite unit starts after one common coin and selects a stationary apolar position law. Actual body half-turns return that law and reverse its retained frame; nonmixing controls preserve polar distinction.',
             domain='Declared twelve-port body, fair native bond events, full oriented position reader and Euclidean chart; field means operational source distribution. No spinor or exchange theorem.',
             falsifier='Different native opposite-start endpoint laws, a polar control made identical by the reader, erased axis information, or failure of an exact body spatial half-turn.',
             measured=dict(contexts=120, exact_body_rotations=30, suffix_histories=184320,
                           native_calls=rotation['literal_native_calls'], coupling_native_calls=720)),
        dict(id='actual_old240_selected_unit_handoff', rows=[1, 7], packet='bridge',
             statement='Signed total-one multiload preparations become actual basis units on the full carrier. Current records identify the unit and the same live source continues into discovered-bond evolution without reset.',
             domain='First-cell scheduler in OLD240; total-one and exterior-zero preparation promises; native selector and census. Coherent address encoding is a declared basis isometry.',
             falsifier='Record reconstruction disagrees with the actual selected unit, a hidden reset changes the handoff, or a retained selector/census failure is counted as success.',
             measured=dict(native_calls=438693, native_handoffs=36, control_selections=0, selector_rejections=5061)),
        dict(id='full_graph_compact_native_contact', rows=[8], packet='contact',
             statement='A corrected24-coordinate compact antisymmetric pair has full-graph Dirichlet eigenvalue12 and native contact probability1/1260; HOLD and HOLD/SWAP give zero removal with the same reader.',
             domain='OLD240 signed positive instrument with complete current flags; prepared pair and ideal complete ray diagnostic. No radiation or autonomous energy-transfer claim.',
             falsifier='A full-graph eigen-residual, incorrect fine contact probability, branch collision or omitted complement.',
             measured=dict(contact_probability='1/1260', unread_ray_probability='932/945',
                           prepared_records=30240, literal_native_calls=30266)),
        dict(id='symmetric_scalar_compact_orbital', rows=[6], packet='scalar',
             statement='The symmetric74 graph has scalar compact eigenrays with Laplacian eigenvalue1 and60 fixed-graph rotations. The complete signed native instrument is checked.',
             domain='Declared74-port graph and positive memory lift. The spin factor is supplied; unread fidelity11867/12144 and escape1/506 are retained. No physical spin claim.',
             falsifier='A full-graph residual, non-scalar rotation action, incomplete source branch or wrong exact unread-channel probability.',
             measured=dict(vertices=74, seams=253, rotations=60, unread_fidelity='11867/12144', unread_escape='1/506')),
        dict(id='native_poisson_source_generator', rows=[2], packet='clock',
             statement='A native auxiliary unit supplies the event-count sampler. Poisson composition of the actual native transfer gives generator1260(I-T); separately enumerated single and signed-pair sectors have values6 and12. The fixed-count control remains nonadditive.',
             domain='Declared source-clock tick1/1260 and original OLD240 compact sectors; ideal uncapped theorem plus finite flagged count/address implementation. Coherent continuation and physical rest-energy identification are separate obligations.',
             falsifier='A nonzero full-coordinate native eigen-residual, incorrect count-bin instrument, omitted unresolved flag, or failure of the changed-clock control to exhibit its native discrepancy.',
             measured=dict(native_calls=242175, source_columns=100800, prefix_columns=131070,
                           actual_clocked_events=505, complete_histories=8,
                           single_generator_values=[6, 6, 6, 6], pair_generator_values=[12, 12, 12, 12],
                           retained_prefix_failure='1/8192', fixed_count_factor_gap='1/44100')),
        dict(id='same_compact_source_spatial_return', rows=[2, 7], packet='apolar_join',
             statement='Actual source-unit checkpoints continue into a native region-selection pump. Its compact energy-sector density has a directional apolar position field; committing coordinate records under the same active pump destroys the half-turn return beyond complete error bounds.',
             domain='First OLD240 cell; declared Hlog selection controller and complete native coordinate reader. Same selected sector has L eigenvalue6; its energy instrument is separately required. Conditional on admitted incoming unit, unconditional over route/pump/reader outcomes.',
             falsifier='Failure of source/checkpoint continuity, a full-graph mode residual, coordinate return outside the complete native bound, or a changed-memory control indistinguishable within that bound.',
             measured=dict(fresh_native_routing_calls=293, continued_checkpoints=36,
                           epochs=512, full_position_halfturn_TV_upper='0.000496891',
                           active_dephasing_halfturn_TV_lower='0.035881369',
                           consumed_executed_native_primitives=196316)),
        dict(id='actual_native_signed_pair_loading', rows=[1, 2], packet='pair_loader',
             statement='One continuing native integer source selects actual positive and negative units in disjoint cells. Current-only censuses reconstruct their pair codes; all144 loaded basis states intertwine with the original full native signed-pair action, including contact/vacuum and exterior outputs.',
             domain='Declared confined schedules and local totals+1/-1; identity concatenation on the admitted disjoint-cell chart. Complete quantum sector/AND operations are supplied by the energy consumer, not inferred from classical encoding.',
             falsifier='A current-only pair reconstruction mismatch, failure of charge-conjugate fine-record equality, mismatch of a loaded native action column, or conflation of collision and source vacuum.',
             measured=dict(native_calls=472146, actual_native_pairs=12, selected_controls=0,
                           loaded_native_action_columns=181440, complete_product_words=65536)),
        dict(id='actual_native_common_pair_routing', rows=[1, 2, 7], packet='pair_route',
             statement='Every actual selected signed pair continues without reset to delta0-delta14 through native current-controlled routes, retaining an actual positive0-to1-to0 frame movement. This supplies the canonical common pump input.',
             domain='Twelve admitted signed-pair handoffs, confined native routes with32-event caps per move, explicit timeout outcomes and full current observer. Coherent loading remains a separately implemented operation.',
             falsifier='A discontinuity in the actual source, an observer mismatch, a wrong final pair/frame or a dropped timeout.',
             measured=dict(actual_continuations=12, literal_native_calls=150, hold_control_calls=384,
                           maximum_moves=8, timeout_upper='1-(1-2^-32)^8')),
        dict(id='full_native_logarithmic_source_flow', rows=[3, 4], packet='flow',
             statement='Native processor words implement the complete OLD240 logarithmic transfer flow with all source modes and failure branches. The held-seam source destroys compact persistence, and complete dephasing/reset destroy the source-reference interference contrast.',
             domain='Declared coherent Hlog=-1260logT program with proper-code/load/export/reset/read services. Chosen-ray persistence is distinct from dynamic state selection or physical self-starting oscillation.',
             falsifier='A full-basis or native-primitive certificate failure, excessive complete error, missing changed-source escape, or failure of coherence controls to separate.',
             measured=dict(source_endpoint_census_calls=302400, pulse_applications=8725933,
                           native_matrix_executions=196316, source_dimension=240, padding_labels=16,
                           held_seam_escape_at4=0.15828662859338302, full_clean_artifacts_byte_equal=10)),
        dict(id='native_selected_compact_source_energy', rows=[2], packet='energy',
             statement='Native region/escape/reinjection dynamics select the compact source sectors from actual routed units. The selected states are held under the count-source generator and read by a source-coupled native flag/AND phase instrument. One511keV anchor gives a separated-pair interval containing1022keV and excluding the old fixed-count excess.',
             domain='Declared native coherent controls, source-address preparation memories and disjoint signed-source loading; Hlog preparation and Poisson-count L holding are distinct modes. Complete finite instrument bounds are not finite-shot confidence intervals or an empirical measurement of nature.',
             falsifier='Failure of the full-source contraction, native source/flag coupling or actual signed handoff; a complete pair ratio excluding2; an undiscriminated changed clock; or missing complement and conditioning error.',
             criterion_support='Local scientific C3 passes on the declared construction; the common packet implements the full loading join. Other-group G6 and official award are separate.',
             measured=dict(actual_routed_pairs=12, ideal_single_generator=6, ideal_pair_generator=12,
                           component_complete_error_upper=float(Fraction(energy_routed['complete_join_error_upper'])),
                           pair_to_single_interval=energy_routed['pair_to_single_interval'],
                           pair_keV_at_single511_anchor_interval=energy_routed['pair_keV_at511_anchor_interval'],
                           source_connected_reader_pulse_applications=44970,
                           fresh_native_primitive_matrices=45, literal_native_vacuum_calls=1260)),
        dict(id='common_selected_compact_source', rows=[1, 2, 7], packet='common',
             statement='The actual signed-pair handoffs feed one declared native preparation and a complete broad regional loader before assay choice. Energy and positive-constituent spatial observations consume the same selected source densities. Fine current-history effects are diagonal partial injections, transferring the exact finite-domain stopped-current response to the loaded coherent source.',
             domain='OLD240 disjoint-cell signed-source chart and declared native coherent control architecture. Preparation address memories are distinct from the loaded physical pair. Fresh repetitions implement incompatible assays. Regional charge and source response support the selected-unit route; independent topology, electric identification and a general spinor/exchange law are not inferred.',
             falsifier='An actual handoff or loaded-action mismatch, different preparation densities, a dark-sector filter before the spatial assay, a failed full loader predicate, excessive complete error, nonseparation of the active dephasing control, or a pair interval excluding2.',
             criterion_support='Same-group finite source-derived unit/energy/apolar-field evidence. The full other-group review and exact outgoing submission remain separate gates.',
             measured=dict(actual_routed_pairs=12, loaded_fine_action_columns=181440,
                           density_max_difference=common['field_energy_density_max_difference'],
                           broad_loader_error_upper=float(Fraction(common['complete_broad_loader_error_upper'])),
                           energy_complete_error_upper=float(Fraction(common['energy_complete_error_upper'])),
                           positive_position_halfturn_TV_upper=float(Fraction(common['spatial_halfturn_TV_upper'])),
                           active_dephasing_TV_lower=float(Fraction(common['active_dephasing_TV_lower'])),
                           pair_to_single_interval=common['pair_to_single_interval'],
                           pair_keV_at_single511_anchor_interval=common['pair_keV_at_single511_anchor_interval'],
                           compressed_native_routine_applications=68,
                           fresh_native_calls=0, fresh_fullflow_executions=0)),
    ]
    claims.extend([
        dict(id='record_only_unknown_total_unit_extraction', rows=[1], packet='unknown_total',
             statement='Current-only reconstruction, potential descent and native nearest-opposite-pair cleanup produce the actual canonical signed pair without a supplied unit total. Reversed signs are physically routed; zero input gives a separate no-unit outcome.',
             domain='Two confined twelve-port cells with unknown opposite totals K and -K, |K| <= 12, initial absolute loads <= 64 and initial vacuum datums at 0 and 14. The remainder starts empty. Integer alphabet, neutral capacity, controller and retained count are supplied.',
             falsifier='A reconstructed current state disagrees with the actual source; a strict native repair fails potential descent; the same-reader control reaches the selected class; or a claimed handoff leaves any surplus source or uses a reset.',
             reusable_theorem='On a connected graph, a shortest path for a globally nearest opposite occupied pair has vacant interior. Native unit moves and terminal annihilation clear any neutral binary signed source in at most k(N-1) move requests, with capped failure at most k(N-1)2^-C.',
             measured=dict(native_calls=5576, actual_handoffs=24, physically_reversed_orientations=12,
                           no_unit_outcomes=1, rejected_controls=25,
                           complete_preparation_inconclusive_upper='8501/288230376151711744')),
        dict(id='full_pair_global_spatial_return', rows=[7], packet='aligned_pair',
             statement='One declared injective spatial chart gives a single proper half-turn preserving the full selected signed-pair density and reversing the recorded frame. The complete 65536-bin native position instrument distinguishes this return from active coordinate dephasing with the same chart and preparation.',
             domain='Independent reducing compact preparations from zero selected overlap, the implemented broad loader, supplied target-informed second-body alignment and sixteen native address readers. Finite comparison bounds do not assert exact finite-state return or spinor exchange.',
             falsifier='The full tensor density fails the common physical rotation, any native outcome is omitted, or the same-chart active control overlaps the complete native return bound.',
             reusable_theorem='For independent states (1-w_i)tau_i+sigma_i with positive orthogonal sigma_i, trace distance from the selected product is exactly 1-product_i(1-w_i). Marginal agreement is insufficient without this structure.',
             measured=dict(joint_address_outcomes=65536, terminal_bit_readers=16,
                           full_joint_halfturn_TV_upper=float(Fraction(aligned['full_joint_halfturn_TV_upper'])),
                           active_control_TV_lower=float(Fraction(aligned['active_dephasing_full_joint_TV_lower'])),
                           retained_unaligned_joint_TV='1/3', fresh_native_calls=0)),
        dict(id='source_local_native_count_continuation', rows=[2], packet='source_local_flow',
             statement='Verified native edge blocks determine the entire charged generator. A specialized local-word compiler takes graph edges, time and refinement, without selected spectral values or a compact-sector filter. Native pointer records distinguish held-edge and unprepared-source alternatives; the independently counted pair passes the one-anchor composition test.',
             domain='Admitted positive one-unit calibration chart and distinct signed-pair chart, with supplied coherent processor and auxiliary count clock. Vacuum is a separately typed zero reference; contact comparison is a source-operator norm audit, not a measured fringe. Finite product accuracy is certified on the selected input family.',
             falsifier='A native charged-block identity or full coupled-reflection failure; an excessive complete implementation error; a native reader insensitive to source state/law; inconsistent clock durations; or a nonadditive pair interval.',
             reusable_theorem='For native endpoint projections T_e, H_e=2(I-T_e). The charged pair generator is the collision-Dirichlet compression of L tensor I + I tensor L; on the antisymmetric subspace it is the exterior-square generator. Distinct annihilation records keep their separate probabilities.',
             measured=dict(literal_pinned_source_calls=384324, literal_clock_calls=258,
                           executed_native_primitive_matrices=22, coupled_pair_truth_columns=72273600,
                           represented_pair_pulses_per_setting=15439539200,
                           common_clock_slots_per_time=154395440640, failed_gate_exit_tests=23,
                           held_edge_fringe_separation=local['held_edge_fringe_separation'],
                           unprepared_source_fringe_separation=local['unprepared_source_fringe_separation'],
                           contact_operator_norm_difference=local['contact_tensor_difference'],
                           pair_keV_at_single511_anchor_interval=local['pair_keV_at_single511_anchor_interval'])),
        dict(id='common_record_selected_source_local_generator', rows=[1, 2, 7], packet='common_v2',
             statement='Actual unknown-total native extraction supplies the full signed body for source-local energy and joint spatial assays after identical compact pumps and broad loading. A separately routed actual positive-only calibration uses the same first pump. Complete native pointer, joint outcomes, common count clock and all preparation/failure costs are joined on the declared menu.',
             domain='Twenty-four nonzero members of the bounded neutral family and thirty-six admitted positive-only calibration sources, with separate preparation domains. Supplied integer alphabet, vacuum data, architecture, aligned chart, retained count and one511keV anchor remain explicit. The all-input preparation bound is proved on its declared family.',
             falsifier='Any actual source/interface mismatch, a calibration replaced by a free partial trace, disagreement of the full selected body or native pointers, missing initial/terminal failures, or loss of native/control separation after complete errors.',
             measured=dict(actual_pair_handoffs=24, actual_positive_calibrations=36,
                           full_selected_density_max_difference=common_v2['full_selected_energy_spatial_density_max_difference'],
                           energy_complete_error_upper={k:float(Fraction(v)) for k,v in common_v2['energy_complete_errors_upper'].items()},
                           joint_halfturn_TV_upper=float(Fraction(common_v2['full_joint_halfturn_TV_upper'])),
                           active_dephasing_TV_lower=float(Fraction(common_v2['active_dephasing_full_joint_TV_lower'])),
                           pair_keV_at_single511_anchor_interval=common_v2['pair_keV_at511_anchor'],
                           fresh_native_calls=0, fresh_fullflow_executions=0)),
    ])
    paths = [Path(__file__), PACKAGE/'ORIGINS.json', PACKAGE/'README.md', PACKAGE/'replay.py',
             ROOT/'paper/tex_fragments/NATIVE_SELECTION_AND_SOURCE_RESPONSE.tex']
    for name, origin in origins['packets'].items():
        paths += [PACKAGE/name/'MANIFEST.sha256', PACKAGE/name/origin.get('receipt_path', 'RECEIPT.json')]
    return dict(schema='oph.native_source_repairs.v1', claims=claims,
                evidence_file_counts=file_counts, origins=origins,
                common_source_scope='The common_v2 packet joins actual unknown-total signed-pair extraction, native source-memory selection and broad loading before unfiltered source-local energy or full joint spatial assays. A distinct actual positive-only calibration joins the same first pump and common auxiliary count clock. The runtime-input wording correction is preserved in reviews/common_v2_interpretation.md. Current-history diagonal effects transfer the finite-domain mean response theorem; a classical discrepancy quotient is not promoted to a coherent encoding.',
                physical_predictions_added=0, official_OpenWave_awards_added=0,
                full_other_group_G6=False, submission_ready=False,
                review_group='Codex same-group, with separate implementations and cross-worker reviews',
                retained_failed_proposal='contact_failed',
                paper_fragment='paper/tex_fragments/NATIVE_SELECTION_AND_SOURCE_RESPONSE.tex',
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != text:
        raise SystemExit('Stale source repair claim ledger')
    print('PASS_NATIVE_SOURCE_REPAIR_CLAIMS')


if __name__ == '__main__':
    main()
