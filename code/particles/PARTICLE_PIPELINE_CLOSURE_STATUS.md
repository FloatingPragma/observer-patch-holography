# Particle Pipeline Closure Status

Generated: `2026-06-29T10:22:10Z`

Single closure gate for source-only rows and empirical hadron closure rows.

## Scope

- Scope: `source_only_rows_plus_empirical_hadron_closure_policy`
- Source-only hadrons in local scope: `False`
- Empirical hadron closure surface: `True`
- Chrome workers needed: `False`
- Hadron scope reason: Source-only hadron rows require a working OPH hardware backend such as GLORB/Echosahedron. Empirical hadron closure stays in a separate output class; the e+e- spectral payload has a source registry and schema.

## Receipt Gates

| Receipt label | Closable | Local artifact | Worker policy |
| --- | --- | --- | --- |
| `closed_exact_selected_branch` | `True` | `particles/hierarchy/issue_332_rg_naturality_certificate.json` | not_needed_for_closed_certificate |
| `closed_projection_bridge_with_exact_residual` | `True` | `particles/hierarchy/certificates/R_EW_tick_projection_certificate.json` | not_needed_for_closed_certificate |
| `closed_full_local_global_hierarchy_resonance` | `True` | `particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json` | not_needed_for_closed_certificate |
| `closed_blocker_isolated_endpoint_package` | `True` | `P_derivation/runtime/thomson_endpoint_package_current.json` | not_needed_for_closed_package |
| `closed_blocker_isolated_source_residual_no_go` | `True` | `P_derivation/runtime/thomson_endpoint_contract_current.json` | not_needed_until_source_spectral_measure_payload_exists |
| `closed_canonical_guarded_trunk_adoption` | `True` | `P_derivation/runtime/p_closure_trunk_current.json` | not_needed_for_guarded_codepath_closure |
| `closed_material_sync_no_live_publish` | `True` | `paper/deriving_the_particle_zoo_from_observer_consistency.tex` | not_needed_for_material_sync |
| `closed_declared_convention_contract` | `True` | `P_derivation/runtime/rg_matching_threshold_contract_current.json` | not_needed_for_closed_contract |
| `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| `closed_corpus_limited_charged_end_to_end_no_go` | `True` | `particles/runs/leptons/charged_end_to_end_impossibility_theorem.json` | not_needed_until_new_uncentered_trace_lift_source_exists |
| `closed_corpus_limited_codomain_no_go` | `True` | `particles/runs/calibration/direct_top_bridge_contract.json` | not_needed_until_new_response_kernel_source_exists |
| `closed_provenance_ledger_and_declared_sensitivity_taxonomy` | `True` | `particles/runs/status/blind_prediction_provenance.json` | not_needed_for_closed_provenance_taxonomy |
| `closed_keep_visible_comparison_tension` | `True` | `particles/runs/neutrino/neutrino_lane_closure_contract.json` | not_needed |
| `closed_selected_class_scope_visible` | `True` | `particles/runs/flavor/quark_lane_closure_contract.json` | not_needed |
| `closed_corpus_limited_global_classification_no_go` | `True` | `particles/runs/flavor/quark_class_uniform_public_frame_descent_obstruction.json` | not_needed_until_new_global_public_frame_classifier_source_exists |
| `open_theta_qcd_bar_theta_vanishing_gap` | `False` | `particles/runs/status/particle_derivation_gap_ledger.json` | not_needed_until_a_concrete_strong_cp_packet_exists |

## Companion Claim Boundaries

| Topic | Claim label | Boundary | Gate |
| --- | --- | --- | --- |
| Strong CP | `open_theta_qcd_bar_theta_vanishing_gap` | The selected-class quark support wrapper conditionally carries the running-quark sextet and exact forward Yukawas on the public class f_P, but its sigma datum is still target-derived. The available corpus does not derive theta_QCD, does not emit the physical anomaly-invariant bar(theta), and does not prove that the physical strong-CP phase vanishes. | Keep strong CP explicit as an open branch. Reopen only for a theorem-grade descent from exact quark/Yukawa phase data to the determinant-line phase contribution, together with a theorem fixing the topological-angle contribution and proving the physical strong-CP phase vanishes on the realized branch. |

## Latest Non-Hadron Predictions

| Particle ID | Mass |
| --- | ---: |
| `gluon` | `0.0 GeV` |
| `graviton` | `0.0 GeV` |
| `higgs` | `125.1995304097179 GeV` |
| `photon` | `0.0 GeV` |

## Withheld Non-Prediction Rows

| Particle ID | Claim label | Reason |
| --- | --- | --- |
| `electron` | `exact_target_anchored_current_family_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `muon` | `exact_target_anchored_current_family_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `tau` | `exact_target_anchored_current_family_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `up_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `charm_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `top_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `down_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `strange_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `bottom_quark` | `selected_class_target_anchored_exact_witness` | target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction |
| `electron_neutrino` | `scale_free_weighted_cycle_theorem_with_compare_only_absolute_attachment_candidate` | compare_only_absolute_or_adapter_surface_kept_out_of_public_prediction_table |
| `muon_neutrino` | `scale_free_weighted_cycle_theorem_with_compare_only_absolute_attachment_candidate` | compare_only_absolute_or_adapter_surface_kept_out_of_public_prediction_table |
| `tau_neutrino` | `scale_free_weighted_cycle_theorem_with_compare_only_absolute_attachment_candidate` | compare_only_absolute_or_adapter_surface_kept_out_of_public_prediction_table |

## Finalization Gates

- `nonhadron_prediction_surface_buildable`: `True`
- `source_only_hadrons_suppressed_by_default`: `True`
- `empirical_hadron_closure_policy_documented`: `True`
- `empirical_hadron_spectral_dataset_integrated`: `False`
- `p_trunk_candidate_only`: `True`
- `obstruction_only_worker_result_allowed`: `True`
- `paper_material_sync_complete_without_live_publish`: `True`
- `source_spectral_stage_gate`: `populated source spectral measure payload plus interval certificate`
- `hierarchy_local_global_resonance_closed`: `True`
- `higgs_naturality_defect_closed`: `True`
- `pixel_screen_resonance_summary_closed`: `True`
