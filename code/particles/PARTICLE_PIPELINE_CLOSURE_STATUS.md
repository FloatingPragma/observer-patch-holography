# Particle Pipeline Closure Status

Generated: `2026-05-06T12:38:37Z`

Single simplified closure gate for the non-hadron particle pipeline.

## Scope

- Scope: `nonhadron_particles_plus_candidate_P_root_metadata`
- Hadrons in local scope: `False`
- Chrome workers needed: `False`
- Hadron scope reason: Production hadrons require a working OPH hardware backend such as GLORB/Echosahedron. Issues #153 and #157 are closed as out-of-scope/computationally blocked for the pipeline; local surrogate code and Chrome workers are non-promoting.

## Issue Gates

| Issue | State | Closable | Local next artifact | Chrome policy |
| --- | --- | --- | --- | --- |
| #223 | `closed_blocker_isolated_endpoint_package` | `True` | `P_derivation/runtime/thomson_endpoint_package_current.json` | not_needed_for_closed_package |
| #235 | `closed_blocker_isolated_source_residual_no_go` | `True` | `P_derivation/runtime/thomson_endpoint_contract_current.json` | not_needed_until_source_spectral_measure_payload_exists |
| #224 | `closed_canonical_guarded_trunk_adoption` | `True` | `P_derivation/runtime/p_closure_trunk_current.json` | not_needed_for_guarded_codepath_closure |
| #225 | `closed_material_sync_no_live_publish` | `True` | `paper/deriving_the_particle_zoo_from_observer_consistency.tex` | not_needed_for_material_sync |
| #32 | `closed_declared_convention_contract` | `True` | `P_derivation/runtime/rg_matching_threshold_contract_current.json` | not_needed_for_closed_contract |
| #153 | `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| #157 | `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| #201 | `closed_corpus_limited_charged_end_to_end_no_go` | `True` | `particles/runs/leptons/charged_end_to_end_impossibility_theorem.json` | not_needed_until_new_uncentered_trace_lift_source_exists |
| #207 | `closed_corpus_limited_codomain_no_go` | `True` | `particles/runs/calibration/direct_top_bridge_contract.json` | not_needed_until_new_response_kernel_source_exists |
| #234 | `closed_provenance_ledger_and_declared_sensitivity_taxonomy` | `True` | `particles/runs/status/blind_prediction_provenance.json` | not_needed_for_closed_provenance_taxonomy |
| #117 | `closed_keep_visible_comparison_tension` | `True` | `particles/runs/neutrino/neutrino_lane_closure_contract.json` | not_needed |
| #198 | `closed_selected_class_scope_visible` | `True` | `particles/runs/flavor/quark_lane_closure_contract.json` | not_needed |
| #199 | `closed_corpus_limited_global_classification_no_go` | `True` | `particles/runs/flavor/quark_class_uniform_public_frame_descent_obstruction.json` | not_needed_until_new_global_public_frame_classifier_source_exists |
| #155 | `open_theta_qcd_bar_theta_vanishing_gap` | `False` | `particles/runs/status/particle_derivation_gap_ledger.json` | not_needed_until_a_concrete_strong_cp_packet_exists |

## Companion Status Branches

| Topic | State | Current boundary | Next action |
| --- | --- | --- | --- |
| Strong CP | `open_theta_qcd_bar_theta_vanishing_gap` | The selected-class exact quark theorem emits the PDG 2025 running-quark sextet and exact forward Yukawas on the public class f_P, but the current corpus does not derive theta_QCD, does not emit the physical anomaly-invariant bar(theta), and does not prove that the physical strong-CP phase vanishes. | Keep strong CP explicit as an open branch. Reopen only for a theorem-grade descent from exact quark/Yukawa phase data to the determinant-line phase contribution together with a theorem fixing the topological-angle contribution and proving the physical strong-CP phase vanishes on the realized branch. |

## Latest Non-Hadron Predictions

| Particle ID | Mass |
| --- | ---: |
| `bottom_quark` | `4.182999999999994 GeV` |
| `charm_quark` | `1.2729999999999992 GeV` |
| `down_quark` | `0.004699999999999999 GeV` |
| `electron` | `0.0005109989499999994 GeV` |
| `electron_neutrino` | `0.017454720257976796 eV` |
| `gluon` | `0.0 GeV` |
| `graviton` | `0.0 GeV` |
| `higgs` | `125.1995304097179 GeV` |
| `muon` | `0.10565837550000004 GeV` |
| `muon_neutrino` | `0.019481987935919015 eV` |
| `photon` | `0.0 GeV` |
| `strange_quark` | `0.09349999999999999 GeV` |
| `tau` | `1.7769324651340912 GeV` |
| `tau_neutrino` | `0.05307522145074924 eV` |
| `top_quark` | `172.35235532883115 GeV` |
| `up_quark` | `0.0021600000000000005 GeV` |
| `w_boson` | `80.377 GeV` |
| `z_boson` | `91.18797809193725 GeV` |

## Finalization Gates

- `nonhadron_prediction_surface_buildable`: `True`
- `hadrons_suppressed_by_default`: `True`
- `p_trunk_candidate_only`: `True`
- `obstruction_only_worker_result_allowed`: `True`
- `paper_material_sync_complete_without_live_publish`: `True`
- `source_spectral_stage_gate`: `populated source spectral measure payload plus interval certificate`
