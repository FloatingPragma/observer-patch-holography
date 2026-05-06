# Particle Pipeline Closure Status

Generated: `2026-05-06T04:24:20Z`

Single simplified closure gate for the current non-hadron particle pipeline.

## Scope

- Current scope: `nonhadron_particles_plus_candidate_P_root_metadata`
- Hadrons in current local scope: `False`
- Chrome workers needed now: `False`
- Hadron scope reason: Production hadrons require a working OPH hardware backend such as GLORB/Echosahedron. Issues #153 and #157 are closed as out-of-scope/computationally blocked for the current pipeline; local surrogate code and Chrome workers are non-promoting.

## Issue Gates

| Issue | State | Closable now | Local next artifact | Chrome policy |
| --- | --- | --- | --- | --- |
| #223 | `open_constructive_contract` | `False` | `P_derivation/runtime/thomson_endpoint_contract_current.json` | not_useful_until_source_endpoint_packet_exists |
| #224 | `open_waiting_certified_root` | `False` | `P_derivation/runtime/p_closure_trunk_current.json` | not_useful_until_endpoint_and_interval_gates_close |
| #32 | `open_constructive_contract` | `False` | `P_derivation/runtime/rg_matching_threshold_contract_current.json` | only_after_beta_threshold_packet_is_populated |
| #153 | `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| #157 | `closed_out_of_scope_computationally_blocked` | `True` | `particles/runs/hadron/ward_projected_spectral_measure_contract.json` | do_not_use_for_backend_execution |
| #201 | `closed_current_corpus_charged_end_to_end_no_go` | `True` | `particles/runs/leptons/charged_end_to_end_impossibility_theorem.json` | not_needed_until_new_uncentered_trace_lift_source_exists |
| #207 | `closed_current_corpus_codomain_no_go` | `True` | `particles/runs/calibration/direct_top_bridge_contract.json` | not_needed_until_new_response_kernel_source_exists |
| #234 | `open_provenance_ledger_emitted_sensitivity_open` | `False` | `particles/runs/status/blind_prediction_provenance.json` | only_for_audit_after_convention_sensitivity_sweep_exists |
| #117 | `closed_keep_visible_comparison_tension` | `True` | `particles/runs/neutrino/neutrino_lane_closure_contract.json` | not_needed |
| #198 | `closed_selected_class_scope_visible` | `True` | `particles/runs/flavor/quark_lane_closure_contract.json` | not_needed |
| #199 | `closed_current_corpus_global_classification_no_go` | `True` | `particles/runs/flavor/quark_class_uniform_public_frame_descent_obstruction.json` | not_needed_until_new_global_public_frame_classifier_source_exists |

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
