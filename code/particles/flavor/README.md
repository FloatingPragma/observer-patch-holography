# Flavor Lane

This directory is the `/particles` sandbox for turning the current flavor
continuation into a forward matrix-based branch.

The intended chain is:

1. derive a refinement-stable flavor observable from overlap / defect data
2. derive a common sector response object for `u,d,e,nu`
3. derive entrywise suppression and phase tensors
4. build complex forward Yukawa matrices `Y_u` and `Y_d`
5. compute singular values, left diagonalizers, `V_CKM`, and Jarlskog
6. export blind forward artifacts before any compare surface is attached

## How To Read The Active Flavor Files

The live flavor/quark scripts now start with a compact derivation summary that
states:

- `Chain role`: how the file fits into the current flavor-to-mass path
- `Mathematics`: the key transport, factorization, or spectral step
- `OPH-derived inputs`: which active `/particles` artifacts it consumes
- `Output`: the artifact it emits and the next residual object if the lane is
  still open

For the mass-facing quark path, the active tail is:

- `derive_quark_sector_mean_split.py`
- `derive_quark_sector_descent.py`
- `derive_quark_diagonal_B_odd_source_scalar_evaluator.py`
- `build_forward_yukawas.py`

Current scripts:

- [`derive_family_transport_kernel.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_family_transport_kernel.py)
- [`derive_generation_bundle_branch_generator.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_generation_bundle_branch_generator.py)
- [`derive_overlap_edge_line_lift.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_overlap_edge_line_lift.py)
- [`derive_overlap_edge_transport_cocycle.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_overlap_edge_transport_cocycle.py)
- [`derive_overlap_flavor_observable.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_overlap_flavor_observable.py)
- [`derive_sector_transport_pushforward.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_sector_transport_pushforward.py)
- [`derive_charged_budget_pushforward.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_charged_budget_pushforward.py)
- [`derive_suppression_phase_tensors.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_suppression_phase_tensors.py)
- [`derive_charged_dirac_odd_deformation_form.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_charged_dirac_odd_deformation_form.py)
- [`derive_quark_odd_response_law.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_quark_odd_response_law.py)
- [`derive_quark_sector_descent.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/derive_quark_sector_descent.py)
- [`build_forward_yukawas.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/build_forward_yukawas.py)
- [`export_flavor_dictionary_artifact.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/export_flavor_dictionary_artifact.py)
- [`export_blind_forward_artifact.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/export_blind_forward_artifact.py)
- [`test_no_ckm_import.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_no_ckm_import.py)
- [`test_flavor_dictionary_disambiguation.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_flavor_dictionary_disambiguation.py)
- [`test_transport_kernel_persistence.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_transport_kernel_persistence.py)
- [`test_sector_pushforward_functoriality.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_sector_pushforward_functoriality.py)
- [`test_sector_residual_factorization.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_sector_residual_factorization.py)
- [`test_charged_budget_partition_invariance.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_charged_budget_partition_invariance.py)
- [`test_shared_budget_refinement_limit.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_shared_budget_refinement_limit.py)
- [`test_shared_budget_reconstruction_identity.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_shared_budget_reconstruction_identity.py)
- [`test_scalarization_label_blindness.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_scalarization_label_blindness.py)
- [`test_sector_local_budget_isolation.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_sector_local_budget_isolation.py)
- [`test_conjugacy_riesz_bound.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_conjugacy_riesz_bound.py)
- [`test_edge_line_lift_boundary.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_edge_line_lift_boundary.py)
- [`test_true_edge_cocycle_identity.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_true_edge_cocycle_identity.py)
- [`test_hermitian_descendant_riesz_margin.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_hermitian_descendant_riesz_margin.py)
- [`test_edge_statistics_nonplaceholder.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_edge_statistics_nonplaceholder.py)
- [`test_cycle_holonomy_from_edge_cocycle.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_cycle_holonomy_from_edge_cocycle.py)
- [`test_observable_certificate_complete.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_observable_certificate_complete.py)
- [`test_quark_placeholder_gate.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_placeholder_gate.py)
- [`test_quark_descent_requires_projector_action.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_descent_requires_projector_action.py)
- [`test_no_entrywise_quark_amplitudes.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_no_entrywise_quark_amplitudes.py)
- [`test_quark_sector_nonclone.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_sector_nonclone.py)
- [`test_quark_budget_neutrality.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_budget_neutrality.py)
- [`test_quark_noncentrality_witness.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_noncentrality_witness.py)
- [`test_degenerate_splitter_fallback_demotes_proof.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_degenerate_splitter_fallback_demotes_proof.py)
- [`test_quark_odd_response_no_hidden_normalization.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_odd_response_no_hidden_normalization.py)
- [`test_quark_zero_odd_scalar_corollary.py`](/Users/muellerberndt/Projects/oph-meta/particles/flavor/test_quark_zero_odd_scalar_corollary.py)

These scripts do **not** claim the OPH flavor observable is already derived.
They establish the artifact boundary and the forward matrix pipeline so the next
math/code work has a concrete home in `/particles`.

Current theorem-shaped local surfaces:

- builder-facing local frontier:
  `source_readback_u_log_per_side_and_source_readback_d_log_per_side`
- broader honest continuation frontier:
  `oph_light_quark_isospin_overlap_defect_selector_law`
  This keeps the recovered-core no-go explicit while exposing the D12
  continuation-level selector shell that would fix the light-sector pure-`B`
  payload pair once its value is honestly emitted.

1. `derive_family_transport_kernel.py` exports a conjugacy-class family kernel
   candidate with refinement intertwiners, conjugacy defects, and three-cluster
   gap certificates.
2. `derive_generation_bundle_branch_generator.py` now exports a centered
   compressed branch-generator candidate on the realized three-generation
   charged bundle, together with the simple-spectrum certificate that has become
   the sharp reduced flavor blocker.
3. `derive_overlap_edge_line_lift.py` now exports the explicit projective
   polar-Riesz common-refinement eigenline transport as a readout of that
   centered generator candidate. Same-label diagonal transport is tracked
   there; off-diagonal flavor-edge overlaps are downstream induced edge
   objects, not the transport maps themselves.
4. `derive_overlap_edge_transport_cocycle.py` exports the induced overlap-edge
   cocycle candidate, with non-placeholder edge amplitudes, cycle holonomy,
   explicit defect/gap bookkeeping, and the lifted Hermitian-descendant Riesz
   margin that now closes the standard-math persistence step on the current
   family.
5. `derive_overlap_flavor_observable.py` exports a persistent-spectral-triple
   candidate with intrinsic labels `f1,f2,f3`, projector certificates, non-floor
   pair suppressions, and cycle holonomy traced to the cocycle artifact.
6. `derive_charged_dirac_odd_deformation_form.py` isolates the remaining odd
   charged-shape / Riesz burden behind the quark response law.
7. `derive_quark_odd_response_law.py` isolates the remaining quark theorem
   burden as an explicit odd-response-law boundary between the persistent flavor
   object and the sector descent.
8. `derive_quark_sector_descent.py` now consumes that boundary and exports a
   projector-resolved odd quark splitter candidate `Xi_q` that separates `u/d`
   in factorized-only mode while keeping theorem status at candidate-only.
9. `build_forward_yukawas.py` now blocks silent promotion of dense-amplitude or
   non-projector-resolved quark artifacts explicitly.

The active local chain is now:

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
11. derive a candidate odd quark sector splitter from that shared boundary
12. build factorized-only forward Yukawas for the quark sectors
13. export one blind dictionary artifact that can later feed charged-lepton and neutrino lanes too
