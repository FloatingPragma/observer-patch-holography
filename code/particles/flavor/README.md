# Flavor Lane

This directory is the `/particles` sandbox for turning the current flavor
continuation into a forward matrix-based branch.

The strict source-only quark-spread audit is closed as a non-identifiability
theorem on the current corpus. The source shape laws leave two independent
positive endpoint spans, giving an exact `(R_{>0})^2` fiber. See
[`derive_quark_sigma_source_nonidentifiability_obstruction.py`](derive_quark_sigma_source_nonidentifiability_obstruction.py).
The target packet also mixes renormalization conventions, and its GeV-valued
matrices are mass textures rather than physical dimensionless Yukawa matrices.
See
[`derive_quark_running_mass_scheme_convention_obstruction.py`](derive_quark_running_mass_scheme_convention_obstruction.py).

The intended chain is:

1. derive a refinement-stable flavor observable from overlap / defect data
2. derive a common sector response object for `u,d,e,nu`
3. derive entrywise suppression and phase tensors
4. build complex forward Yukawa matrices `Y_u` and `Y_d`
5. compute singular values, left diagonalizers, `V_CKM`, and Jarlskog
6. export blind forward artifacts before any compare surface is attached

## How To Read The Active Flavor Files

The live flavor/quark scripts start with a compact derivation summary that
states:

- `Chain role`: how the file fits into the current flavor-to-mass path
- `Mathematics`: the key transport, factorization, or spectral step
- `OPH-derived inputs`: which active `/particles` artifacts it consumes
- `Output`: the artifact it emits and the next residual object if the lane is
  work in progress

For the mass-facing quark path, the active tail is:

- `derive_quark_sector_mean_split.py`
- `derive_quark_sector_descent.py`
- `derive_quark_diagonal_B_odd_source_scalar_evaluator.py`
- `build_forward_yukawas.py`

Current scripts:

- [`derive_family_transport_kernel.py`](derive_family_transport_kernel.py)
- [`derive_generation_bundle_branch_generator.py`](derive_generation_bundle_branch_generator.py)
- [`derive_overlap_edge_line_lift.py`](derive_overlap_edge_line_lift.py)
- [`derive_overlap_edge_transport_cocycle.py`](derive_overlap_edge_transport_cocycle.py)
- [`derive_overlap_flavor_observable.py`](derive_overlap_flavor_observable.py)
- [`derive_sector_transport_pushforward.py`](derive_sector_transport_pushforward.py)
- [`derive_charged_budget_pushforward.py`](derive_charged_budget_pushforward.py)
- [`derive_suppression_phase_tensors.py`](derive_suppression_phase_tensors.py)
- [`derive_charged_dirac_odd_deformation_form.py`](derive_charged_dirac_odd_deformation_form.py)
- [`derive_quark_odd_response_law.py`](derive_quark_odd_response_law.py)
- [`derive_quark_edge_statistics_spread_candidate.py`](derive_quark_edge_statistics_spread_candidate.py)
- [`derive_quark_sigma_source_nonidentifiability_obstruction.py`](derive_quark_sigma_source_nonidentifiability_obstruction.py)
- [`derive_quark_sigma_source_datum_no_target_leak_required.py`](derive_quark_sigma_source_datum_no_target_leak_required.py)
- [`derive_quark_running_mass_scheme_convention_obstruction.py`](derive_quark_running_mass_scheme_convention_obstruction.py)
- [`derive_quark_d12_internal_backread_source_payload.py`](derive_quark_d12_internal_backread_source_payload.py)
- [`derive_quark_d12_internal_backread_yukawa_dictionary.py`](derive_quark_d12_internal_backread_yukawa_dictionary.py)
- [`derive_quark_d12_internal_backread_descent.py`](derive_quark_d12_internal_backread_descent.py)
- [`derive_quark_d12_internal_backread_forward_yukawas.py`](derive_quark_d12_internal_backread_forward_yukawas.py)
- [`derive_quark_sector_descent.py`](derive_quark_sector_descent.py)
- [`build_forward_yukawas.py`](build_forward_yukawas.py)
- [`export_flavor_dictionary_artifact.py`](export_flavor_dictionary_artifact.py)
- [`export_blind_forward_artifact.py`](export_blind_forward_artifact.py)
- [`test_no_ckm_import.py`](test_no_ckm_import.py)
- [`test_flavor_dictionary_disambiguation.py`](test_flavor_dictionary_disambiguation.py)
- [`test_transport_kernel_persistence.py`](test_transport_kernel_persistence.py)
- [`test_sector_pushforward_functoriality.py`](test_sector_pushforward_functoriality.py)
- [`test_sector_residual_factorization.py`](test_sector_residual_factorization.py)
- [`test_charged_budget_partition_invariance.py`](test_charged_budget_partition_invariance.py)
- [`test_shared_budget_refinement_limit.py`](test_shared_budget_refinement_limit.py)
- [`test_shared_budget_reconstruction_identity.py`](test_shared_budget_reconstruction_identity.py)
- [`test_scalarization_label_blindness.py`](test_scalarization_label_blindness.py)
- [`test_sector_local_budget_isolation.py`](test_sector_local_budget_isolation.py)
- [`test_conjugacy_riesz_bound.py`](test_conjugacy_riesz_bound.py)
- [`test_edge_line_lift_boundary.py`](test_edge_line_lift_boundary.py)
- [`test_true_edge_cocycle_identity.py`](test_true_edge_cocycle_identity.py)
- [`test_hermitian_descendant_riesz_margin.py`](test_hermitian_descendant_riesz_margin.py)
- [`test_edge_statistics_nonplaceholder.py`](test_edge_statistics_nonplaceholder.py)
- [`test_cycle_holonomy_from_edge_cocycle.py`](test_cycle_holonomy_from_edge_cocycle.py)
- [`test_observable_certificate_complete.py`](test_observable_certificate_complete.py)
- [`test_quark_placeholder_gate.py`](test_quark_placeholder_gate.py)
- [`test_quark_d12_internal_backread_dictionary.py`](test_quark_d12_internal_backread_dictionary.py)
- [`test_quark_d12_internal_backread_forward_yukawas.py`](test_quark_d12_internal_backread_forward_yukawas.py)
- [`test_quark_descent_requires_projector_action.py`](test_quark_descent_requires_projector_action.py)
- [`test_no_entrywise_quark_amplitudes.py`](test_no_entrywise_quark_amplitudes.py)
- [`test_quark_sector_nonclone.py`](test_quark_sector_nonclone.py)
- [`test_quark_budget_neutrality.py`](test_quark_budget_neutrality.py)
- [`test_quark_noncentrality_witness.py`](test_quark_noncentrality_witness.py)
- [`test_degenerate_splitter_fallback_demotes_proof.py`](test_degenerate_splitter_fallback_demotes_proof.py)
- [`test_quark_odd_response_no_hidden_normalization.py`](test_quark_odd_response_no_hidden_normalization.py)
- [`test_quark_zero_odd_scalar_corollary.py`](test_quark_zero_odd_scalar_corollary.py)
- [`test_quark_edge_statistics_spread_candidate.py`](test_quark_edge_statistics_spread_candidate.py)

These scripts do **not** claim a derived OPH flavor observable. They establish
the artifact boundary and the forward matrix pipeline so the next
math/code work has a concrete home in `/particles`.

Current theorem-shaped local surfaces:

- builder-facing local frontier:
  `source_readback_u_log_per_side_and_source_readback_d_log_per_side`
- broader supported continuation frontier:
  `oph_light_quark_isospin_overlap_defect_selector_law`
  This keeps the recovered-core no-go explicit while exposing the D12
  continuation-level selector shell that would fix the light-sector pure-`B`
  payload pair once its value has a supported emission.
- sharper target-1 theorem package:
  `oph_light_quark_overlap_defect_value_law`
  This is the internalized target-free bridge theorem on the local code
  surface: `Delta_ud_overlap = (1/6) * log(c_d / c_u)`. On the emitted D12
  mass ray it immediately yields the older `quark_d12_t1_value_law` wrapper by
  `t1 = 5 * Delta_ud_overlap`.
- exact current-family light-ratio theorem:
  `oph_quark_current_family_light_ratio_theorem`
  This proves a narrower but valid theorem on the exact same-family witness:
  `ell_ud = log(m_d / m_u)` is fixed exactly by the exact light masses and
  equals `log(g_d / g_u) + E_d_log_exact[0] - E_u_log_exact[0]`. It does not
  promote the public target-free mass theorem `H_mass := log(c_d / c_u)`.
- synthesized target-free bridge theorem package:
  `oph_quark_target_free_bridge_theorem`
  This packages the internalized public bridge theorem constructively:
  `Delta_ud_overlap = (1/6) * log(c_d / c_u)` and equivalently
  `t1 = (5/6) * log(c_d / c_u)`, together with the induced D12/source/transport
  corollaries. On the local code surface this bridge is treated as
  internalized; the remaining full physical quark frontier is the physical-sheet
  lift/readout pair rather than the D12 mass bridge.
- computed current-family target-anchored D12 scalar package:
  `oph_quark_d12_current_family_target_anchored_value_package`
  This computes the missing D12 scalar package explicitly on the exact
  current-family witness: `ell_ud = log(m_d / m_u)`,
  `Delta_ud_overlap = ell_ud / 6`, `t1 = 5 * Delta_ud_overlap`, and the
  induced pure-`B` source and transport package on the theorem-grade sigma
  branch.
- exact D12 transport reduction on any fixed sigma branch:
  `oph_quark_d12_overlap_transport_law`
  This states the transport side exactly: once a sigma branch is fixed, the
  weighted odd transport pair `(tau_u, tau_d)` and `Lambda_ud_B_transport` are
  affine readbacks of the single selector scalar `Delta_ud_overlap`, so the
  remaining scalar burden is the value law itself rather than a free tau-pair.
- constructive continuation-only sidecar emitted:
  `oph_quark_d12_internal_backread_source_payload`
  This materializes the previously named pure-`B` payload pair
  `source_readback_u_log_per_side`, `source_readback_d_log_per_side` together
  with `J_B_source_u`, `J_B_source_d` on the internal D12 backread surface.
- strongest current explicit quark dictionary sidecar:
  `oph_quark_d12_internal_backread_yukawa_dictionary`
  This combines the direct edge-statistics spread bridge with the emitted D12
  source payload and weighted overlap-transport lift to produce one
  continuation-only quark Yukawa dictionary surface.
- strongest current explicit forward-Yukawa sidecar:
  `oph_quark_d12_internal_backread_forward_yukawas`
  This emits actual `Y_u`, `Y_d` and a certified forward matrix surface on
  the D12 internal-backread continuation scope, using the emitted source
  payload and the exact edge-branch weighted transport law.

1. `derive_family_transport_kernel.py` exports a conjugacy-class family kernel
   candidate with refinement intertwiners, conjugacy defects, and three-cluster
   gap certificates.
2. `derive_generation_bundle_branch_generator.py` exports a centered
   compressed branch-generator candidate on the realized three-generation
   charged bundle, together with the simple-spectrum certificate that has become
   the sharp reduced flavor blocker.
3. `derive_overlap_edge_line_lift.py` exports the explicit projective
   polar-Riesz common-refinement eigenline transport as a readout of that
   centered generator candidate. Same-label diagonal transport is tracked
   there; off-diagonal flavor-edge overlaps are downstream induced edge
   objects, not the transport maps themselves.
4. `derive_overlap_edge_transport_cocycle.py` exports the induced overlap-edge
   cocycle candidate, with non-placeholder edge amplitudes, cycle holonomy,
   explicit defect/gap bookkeeping, and the lifted Hermitian-descendant Riesz
   margin that closes the standard-math persistence step on the current
   family.
5. `derive_overlap_flavor_observable.py` exports a persistent-spectral-triple
   candidate with intrinsic labels `f1,f2,f3`, projector certificates, non-floor
   pair suppressions, and cycle holonomy traced to the cocycle artifact.
6. `derive_charged_dirac_odd_deformation_form.py` isolates the remaining odd
   charged-shape / Riesz burden behind the quark response law.
7. `derive_quark_odd_response_law.py` isolates the remaining quark theorem
   burden as an explicit odd-response-law boundary between the persistent flavor
   object and the sector descent.
8. `derive_quark_edge_statistics_spread_candidate.py` exposes the strongest
   current direct bridge from edge-side suppression data to the quark spread
   pair on the realized ordered family, with explicit residuals against the
   active same-family spread pair.
9. `derive_quark_d12_internal_backread_source_payload.py` emits the missing
   pure-`B` quark source payload pair on the D12 internal backread sidecar,
   turning the old named residual object into an explicit artifact.
10. `derive_quark_d12_internal_backread_yukawa_dictionary.py` combines that
   emitted source payload with the edge-statistics spread bridge and the D12
   weighted transport law to produce one continuation-only quark Yukawa
   dictionary sidecar.
11. `derive_quark_d12_internal_backread_descent.py` materializes a
   continuation-only descent artifact with no open pure-`B` payload blocker on
   its own sidecar inputs.
12. `derive_quark_d12_internal_backread_forward_yukawas.py` reuses the live
   forward builder on that descent and emits a certified continuation-only
   forward Yukawa surface with explicit `Y_u` and `Y_d`.
13. `derive_quark_d12_t1_value_law.py` records the reduced exact D12 mass
   boundary more sharply: once `t1` is emitted, the pure-`B` source payload is
   fixed exactly, and the transport side is closed on any fixed sigma
   branch.
14. `derive_quark_sector_descent.py` consumes that boundary and exports a
   projector-resolved odd quark splitter candidate `Xi_q` that separates `u/d`
   in factorized-only mode while keeping theorem status at candidate-only.
15. `build_forward_yukawas.py` blocks silent promotion of dense-amplitude or
   non-projector-resolved quark artifacts explicitly.

The active local chain:

1. normalize a refinement-indexed family transport kernel
2. derive the centered compressed generation-bundle branch-generator candidate
3. derive a projective same-label eigenline transport readout from that generator
4. derive the induced overlap-edge transport cocycle with explicit defect/gap bookkeeping
5. reduce that into family projectors, spectral gaps, pairwise suppressions, and cycle phases
6. derive a common sector response object
7. certify functoriality, cocycle provenance, and residual-factorization guards
8. derive sector suppression/phase tensors plus a shared charged-budget transport artifact
9. derive the charged odd deformation-form boundary
10. derive an explicit odd quark response-law boundary from the full projector algebra
11. expose the strongest current direct edge-statistics candidate for the quark spread pair
12. emit the missing pure-`B` source payload pair on the D12 internal backread sidecar
13. combine the edge-statistics spread bridge with that emitted D12 payload into a continuation-only quark Yukawa dictionary sidecar
14. derive a candidate odd quark sector splitter from the shared boundary on the main builder path
15. build factorized-only forward Yukawas for the quark sectors while preserving any explicit edge-statistics bridge metadata
16. export one blind dictionary artifact that can later feed charged-lepton and neutrino lanes too, including the current constructive quark bridge provenance

## License And Patent Policy

This particle flavor surface is part of the OPH public repository. See the main
[LICENSE](../../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../../PATENTS.md).
