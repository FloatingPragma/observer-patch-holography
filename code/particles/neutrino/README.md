# Neutrino Sandbox

This directory holds the `/particles` neutrino sandbox for the local prediction
chain.

## Scientific Status

The intrinsic isotropic ansatz fails its exact atmospheric-splitting cap. The
weighted-cycle branch is a target-informed template
candidate. Its family transport kernel consists of two hand-written matrices,
and its exponent law was promoted after candidate laws had been ranked against
the measured splitting ratio.

NuFIT 6.1 rejects the frozen point on the official correlated
`T23/DCP` profile. The profile value is `20.11955` with the tabulated
atmospheric likelihood and `18.43528` without it; the two-parameter 3σ contour
value is `11.82916`. Exhaustive enumeration of both stored cycle orientations,
all six flavor-row relabelings, and all six mass-column relabelings finds no
normal-ordering-consistent relabeling of the stored candidate that passes both
gates. This does not derive or exhaust source-side physical charged-basis
placements. The stored
weighted-cycle eigensystem is reproducible. The builders export the Majorana
Takagi matrix; this computational correction does not derive its physical
status. Separately, the shared-basis construction defines
`U_nu_shared = U_e_left U_wc` and then recovers `U_wc` by cancellation. That
identity does not derive the physical charged-lepton basis. The failure belongs
to the weighted-cycle continuation candidate, while finite OPH is unaffected
because the source kernel, selector law, and physical charged basis were never
derived.

No neutrino value in this directory has source-only prediction status. The
lane emits no physical PMNS matrix and no absolute neutrino masses.

The purpose is controlled scope. The active forward objects are:

- a local neutrino scale anchor from the D10 core
- a real symmetric family-response tensor
- an explicit Majorana holonomy-lift boundary
- a residual-phase envelope support-check surface
- a selector-point versus selector-law split
- a blind real-first Majorana matrix artifact
- derived ascending singular-value gaps, with physical ordering explicitly unresolved
- a shared-basis PMNS diagnostic whose physical charged-lepton basis remains open

## Header Convention

The active neutrino scripts begin with a short derivation header that says:

- `Chain role`: where the file sits in the local neutrino chain
- `Mathematics`: the main projector, phase, pullback, or spectral step
- `OPH-derived inputs`: which values are coming directly from the local OPH
  chain in `/particles`
- `Output`: which artifact is emitted and what it unlocks downstream

That convention is there to make the neutrino lane readable without mixing the
local standard-math selector closure with the work-in-progress stricter OPH-only
promotion questions.

## Current Working Split

- [`derive_neutrino_scale_anchor.py`](derive_neutrino_scale_anchor.py)
  Derives the active local neutrino scale anchor directly from the current D10
  core via `m_star = v^2 / mu_u`.
- [`derive_family_response_tensor.py`](derive_family_response_tensor.py)
  Builds the real symmetric family-response tensor from the local scale anchor
  and the common sector-response object.
- [`derive_majorana_holonomy_lift.py`](derive_majorana_holonomy_lift.py)
  Holds the congruence-gauge Majorana phase boundary, selector candidates, and
  the current equal-split selector-point closure.
- [`derive_majorana_deformation_bilinear_form.py`](derive_majorana_deformation_bilinear_form.py)
  Exports the OPH-only deformation-bilinear-form boundary. On the current
  isotropic branch it closes the residual metric class up to scale and leaves
  the upstream OPH Hessian origin theorem explicit.
- [`derive_majorana_overlap_defect_hessian.py`](derive_majorana_overlap_defect_hessian.py)
  Exports the exact upstream Hessian-shaped object isolating the remaining
  OPH-only neutrino blocker.
- [`derive_majorana_overlap_defect_action_germ.py`](derive_majorana_overlap_defect_action_germ.py)
  Exports the local quadratic overlap-defect action germ that sits between
  the selector center and the exact scalar-evaluator obstruction.
- [`derive_majorana_phase_pullback_metric.py`](derive_majorana_phase_pullback_metric.py)
  Exports the pullback metric and finite-angle chordal-distortion action behind
  the selector law. In the sandbox this closes the law under a
  standard-math-fixed Hilbert-Schmidt/Frobenius route, while the stricter
  OPH-only ambient-metric derivation remains separate.
- [`build_forward_majorana_matrix.py`](build_forward_majorana_matrix.py)
  Builds the neutrino matrix artifact in explicit `real_seed`,
  `canonical_selector`, or `residual_envelope` mode.
- [`build_majorana_phase_envelope.py`](build_majorana_phase_envelope.py)
  Sweeps the residual Majorana phase plane and exports gap-vs-radius stability
  certificates for the ascending spectrum and collective-mode location; it
  does not select physical normal or inverted ordering.
- [`build_forward_splittings.py`](build_forward_splittings.py)
  Turns the blind matrix artifact into masses, splittings, and explicit
  real-seed versus phase-certified outputs.
- [`build_pmns_from_shared_flavor_basis.py`](build_pmns_from_shared_flavor_basis.py)
  Builds PMNS only if a charged-lepton left basis is supplied.
- [`derive_neutrino_attachment_bridge_scalar_corridor.py`](derive_neutrino_attachment_bridge_scalar_corridor.py)
  Fuses the strongest current compare-only routes for the diagnostic bridge
  scalar sidecar `B_nu`, including the best defect-weighted `mu_e`-assisted
  route, exposes the smaller exact correction scalar `C_nu` above the best
  internal candidate proxy, and records a narrower shortlist-consensus window on the
  rejected candidate surface.
- [`derive_neutrino_bridge_correction_candidate_audit.py`](derive_neutrino_bridge_correction_candidate_audit.py)
  Audits the reduced correction scalar `C_nu` directly and induces a narrower
  measured-reference compare-only `B_nu` window than the older direct bridge
  corridor. This is a diagnostic surface with no prediction status.
- [`derive_neutrino_bridge_correction_invariant_scaffold.py`](derive_neutrino_bridge_correction_invariant_scaffold.py)
  Records the reduced correction geometry `C_nu` below the raw bridge scalar
  `B_nu` and above the internal positive proxy `P_nu` on the scaffold beneath
  the rejected weighted-cycle candidate.
- [`derive_neutrino_bridge_rigidity_theorem.py`](derive_neutrino_bridge_rigidity_theorem.py)
  Retains the compare-only bridge coordinate
  `C_nu = sum_gap^2 * prod_qbar * solar_response_over_mstar^-0.5`. Its source
  audit blocks theorem and prediction promotion.
- [`derive_neutrino_absolute_attachment_theorem.py`](derive_neutrino_absolute_attachment_theorem.py)
  Displays compare-only values for `B_nu`, `lambda_nu`, and an absolute
  weighted-cycle family. These values are not physical mass predictions.
- [`score_neutrino_nufit61.py`](score_neutrino_nufit61.py)
  Scores the frozen scale-free point against hash-pinned official NuFIT 6.1
  profile tables. Overlapping profile projections are reported separately and
  never summed.
- [`audit_neutrino_pmns_conventions.py`](audit_neutrino_pmns_conventions.py)
  Enumerates the phase orientation and every row/column assignment, separating
  genuine convention checks from post-hoc model changes.
- [`build_neutrino_precision_candidate_lock.py`](build_neutrino_precision_candidate_lock.py)
  Freezes the failed point as an immutable anti-rescue record. The lock records
  `prediction_promotion_allowed: false` until a new source-emitted construction
  is frozen before reference data are opened.

## Guards

- [`test_no_pmns_import.py`](test_no_pmns_import.py)
- [`test_no_oscillation_import.py`](test_no_oscillation_import.py)
- [`test_no_dirac_k3_reuse.py`](test_no_dirac_k3_reuse.py)
- [`test_majorana_residual_factorization.py`](test_majorana_residual_factorization.py)
- [`test_shared_flavor_basis_contract.py`](test_shared_flavor_basis_contract.py)
- [`test_no_zero_phase_shortcut.py`](test_no_zero_phase_shortcut.py)
- [`test_majorana_selector_equivariance.py`](test_majorana_selector_equivariance.py)
- [`test_phase_envelope_gap_certificate.py`](test_phase_envelope_gap_certificate.py)
- [`test_isotropic_selector_reason_gate.py`](test_isotropic_selector_reason_gate.py)
- [`test_pullback_action_matches_hs_distortion.py`](test_pullback_action_matches_hs_distortion.py)
- [`test_pullback_metric_finite_difference.py`](test_pullback_metric_finite_difference.py)
- [`test_selector_point_vs_law_status_split.py`](test_selector_point_vs_law_status_split.py)
- [`test_ambient_metric_s3_rigidity.py`](test_ambient_metric_s3_rigidity.py)
- [`test_oph_only_scope_requires_deformation_form.py`](test_oph_only_scope_requires_deformation_form.py)
- [`test_hessian_recovers_deformation_bilinear_form.py`](test_hessian_recovers_deformation_bilinear_form.py)
- [`test_oph_only_hessian_not_frobenius_primitive.py`](test_oph_only_hessian_not_frobenius_primitive.py)
- [`test_overlap_defect_action_germ_class.py`](test_overlap_defect_action_germ_class.py)
- [`test_overlap_defect_action_scale_requires_evaluator.py`](test_overlap_defect_action_scale_requires_evaluator.py)

## Typical Sandbox Flow

```bash
python3 particles/neutrino/derive_neutrino_scale_anchor.py
python3 particles/neutrino/derive_majorana_holonomy_lift.py
python3 particles/neutrino/derive_family_response_tensor.py
python3 particles/neutrino/derive_majorana_overlap_defect_action_germ.py
python3 particles/neutrino/derive_majorana_overlap_defect_hessian.py
python3 particles/neutrino/derive_majorana_deformation_bilinear_form.py
python3 particles/neutrino/derive_majorana_phase_pullback_metric.py
python3 particles/neutrino/build_forward_majorana_matrix.py
python3 particles/neutrino/build_majorana_phase_envelope.py
python3 particles/neutrino/build_forward_splittings.py
python3 particles/neutrino/derive_neutrino_weighted_cycle_repair.py
python3 particles/neutrino/derive_neutrino_bridge_rigidity_theorem.py
python3 particles/neutrino/derive_neutrino_absolute_attachment_theorem.py
python3 particles/neutrino/export_forward_neutrino_closure_bundle.py
python3 particles/neutrino/export_blind_forward_artifact.py
python3 particles/neutrino/test_no_pmns_import.py
python3 particles/neutrino/test_no_oscillation_import.py
python3 particles/neutrino/test_no_dirac_k3_reuse.py
python3 particles/neutrino/test_majorana_residual_factorization.py
python3 particles/neutrino/test_shared_flavor_basis_contract.py
python3 particles/neutrino/test_no_zero_phase_shortcut.py
python3 particles/neutrino/test_majorana_selector_equivariance.py
python3 particles/neutrino/test_phase_envelope_gap_certificate.py
python3 particles/neutrino/test_isotropic_selector_reason_gate.py
python3 particles/neutrino/test_pullback_action_matches_hs_distortion.py
python3 particles/neutrino/test_pullback_metric_finite_difference.py
python3 particles/neutrino/test_selector_point_vs_law_status_split.py
python3 particles/neutrino/test_ambient_metric_s3_rigidity.py
python3 particles/neutrino/test_oph_only_scope_requires_deformation_form.py
python3 particles/neutrino/test_hessian_recovers_deformation_bilinear_form.py
python3 particles/neutrino/test_oph_only_hessian_not_frobenius_primitive.py
python3 particles/neutrino/test_overlap_defect_action_germ_class.py
python3 particles/neutrino/test_overlap_defect_action_scale_requires_evaluator.py
```

The NuFIT scorer reads official compressed tables supplied on the command
line. Their URLs, byte counts, and SHA-256 values are pinned in
[`nufit61_sources.json`](nufit61_sources.json). The tables are not vendored
because NuFIT publishes no explicit redistribution license.

Current closure split:

- the same-label scalar values are numerically complete, but the family kernel
  is a template and the overlap-line lift is candidate-only; their physical
  source-closure gate therefore remains open
- the normalizer and centered scalar evaluator are exact only conditional on
  those source-open inputs
- the real symmetric branch can carry scale, family shape, and real-seed
  splittings
- the phase envelope is the compulsory support check for sharp splittings or
  ordering
- the selector point closes on the current isotropic branch as a principal
  equal split
- the selector law closes locally under the standard-math-fixed
  Hilbert-Schmidt pullback route
- the stricter OPH-only route has a concrete artifact boundary:
  `majorana_overlap_defect_action_germ.json`,
  `majorana_overlap_defect_hessian.json`, and
  `majorana_deformation_bilinear_form.json` separate the local quadratic
  theorem from the exact missing scalar-evaluator object
- the bridge artifact retains the compare-only coordinate
  `C_nu = sum_gap^2 * prod_qbar * solar_response_over_mstar^-0.5`
- the absolute-attachment artifact displays `B_nu`, `lambda_nu`, masses, and
  splittings only as non-promotable comparison coordinates
- the exact positive-segment adapter, bridge corridor, and correction audit are
  retained only as diagnostic sidecars on the rejected candidate surface
- the weighted-cycle calculation emits a measured-data-informed comparison candidate;
  the NuFIT 6.1 correlated profile rejects it
- the shared-basis constructor recovers the weighted-cycle matrix by an
  algebraic cancellation; the physical charged-lepton basis remains open, so
  the result is not a physical PMNS matrix

The outputs stay under `particles/runs/neutrino/`. Prediction promotion
requires a source-emitted family transport kernel, source-derived cycle and
label laws, refinement stability, and a pre-reference hash lock.

## License

This particle neutrino surface is part of the OPH public repository. See the
main [LICENSE](../../../LICENSE).
