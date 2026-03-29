# Neutrino Sandbox

This directory holds the `/particles` neutrino sandbox on the current local
prediction chain.

The current goal is not to overclaim closure. The goal is to keep the active
forward objects explicit:

- a local neutrino scale anchor from the D10 core
- a real symmetric family-response tensor
- an explicit Majorana holonomy-lift boundary
- a residual-phase envelope honesty-check surface
- a selector-point versus selector-law split
- a blind real-first Majorana matrix artifact
- derived splittings and blind ordering
- conditional PMNS only after a shared charged-lepton basis exists

## Header Convention

The active neutrino scripts now begin with a short derivation header that says:

- `Chain role`: where the file sits in the local neutrino chain
- `Mathematics`: the main projector, phase, pullback, or spectral step
- `OPH-derived inputs`: which values are coming directly from the local OPH
  chain in `/particles`
- `Output`: which artifact is emitted and what it unlocks downstream

That convention is there to make the neutrino lane readable without mixing the
local standard-math selector closure with the still-open stricter OPH-only
promotion questions.

## Current Working Split

- [`derive_neutrino_scale_anchor.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_neutrino_scale_anchor.py)
  Derives the active local neutrino scale anchor directly from the current D10
  core via `m_star = v^2 / mu_u`.
- [`derive_family_response_tensor.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_family_response_tensor.py)
  Builds the real symmetric family-response tensor from the local scale anchor
  and the common sector-response object.
- [`derive_majorana_holonomy_lift.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_majorana_holonomy_lift.py)
  Holds the congruence-gauge Majorana phase boundary, selector candidates, and
  the current equal-split selector-point closure.
- [`derive_majorana_deformation_bilinear_form.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_majorana_deformation_bilinear_form.py)
  Exports the OPH-only deformation-bilinear-form boundary. On the current
  isotropic branch it closes the residual metric class up to scale and leaves
  the upstream OPH Hessian origin theorem explicit.
- [`derive_majorana_overlap_defect_hessian.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_majorana_overlap_defect_hessian.py)
  Exports the exact upstream Hessian-shaped object now isolating the remaining
  OPH-only neutrino blocker.
- [`derive_majorana_overlap_defect_action_germ.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_majorana_overlap_defect_action_germ.py)
  Exports the local quadratic overlap-defect action germ that now sits between
  the selector center and the exact scalar-evaluator obstruction.
- [`derive_majorana_phase_pullback_metric.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/derive_majorana_phase_pullback_metric.py)
  Exports the pullback metric and finite-angle chordal-distortion action behind
  the selector law. In the sandbox this now closes the law under a
  standard-math-fixed Hilbert-Schmidt/Frobenius route, while the stricter
  OPH-only ambient-metric derivation remains separate.
- [`build_forward_majorana_matrix.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/build_forward_majorana_matrix.py)
  Builds the neutrino matrix artifact in explicit `real_seed`,
  `canonical_selector`, or `residual_envelope` mode.
- [`build_majorana_phase_envelope.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/build_majorana_phase_envelope.py)
  Sweeps the residual Majorana phase plane and exports gap-vs-radius stability
  certificates for splittings and ordering.
- [`build_forward_splittings.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/build_forward_splittings.py)
  Turns the blind matrix artifact into masses, splittings, and explicit
  real-seed versus phase-certified outputs.
- [`build_pmns_from_shared_flavor_basis.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/build_pmns_from_shared_flavor_basis.py)
  Builds PMNS only if a charged-lepton left basis is supplied.

## Guards

- [`test_no_pmns_import.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_no_pmns_import.py)
- [`test_no_oscillation_import.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_no_oscillation_import.py)
- [`test_no_dirac_k3_reuse.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_no_dirac_k3_reuse.py)
- [`test_majorana_residual_factorization.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_majorana_residual_factorization.py)
- [`test_shared_flavor_basis_contract.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_shared_flavor_basis_contract.py)
- [`test_no_zero_phase_shortcut.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_no_zero_phase_shortcut.py)
- [`test_majorana_selector_equivariance.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_majorana_selector_equivariance.py)
- [`test_phase_envelope_gap_certificate.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_phase_envelope_gap_certificate.py)
- [`test_isotropic_selector_reason_gate.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_isotropic_selector_reason_gate.py)
- [`test_pullback_action_matches_hs_distortion.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_pullback_action_matches_hs_distortion.py)
- [`test_pullback_metric_finite_difference.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_pullback_metric_finite_difference.py)
- [`test_selector_point_vs_law_status_split.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_selector_point_vs_law_status_split.py)
- [`test_ambient_metric_s3_rigidity.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_ambient_metric_s3_rigidity.py)
- [`test_oph_only_scope_requires_deformation_form.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_oph_only_scope_requires_deformation_form.py)
- [`test_hessian_recovers_deformation_bilinear_form.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_hessian_recovers_deformation_bilinear_form.py)
- [`test_oph_only_hessian_not_frobenius_primitive.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_oph_only_hessian_not_frobenius_primitive.py)
- [`test_overlap_defect_action_germ_class.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_overlap_defect_action_germ_class.py)
- [`test_overlap_defect_action_scale_requires_evaluator.py`](/Users/muellerberndt/Projects/oph-meta/particles/neutrino/test_overlap_defect_action_scale_requires_evaluator.py)

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

Current closure split:

- the real symmetric branch can carry scale, family shape, and real-seed
  splittings
- the phase envelope is the compulsory honesty check for sharp splittings or
  ordering
- the selector point now closes on the current isotropic branch as a principal
  equal split
- the selector law now closes locally under the standard-math-fixed
  Hilbert-Schmidt pullback route
- the stricter OPH-only route now has a concrete artifact boundary:
  `majorana_overlap_defect_action_germ.json`,
  `majorana_overlap_defect_hessian.json`, and
  `majorana_deformation_bilinear_form.json` now separate the local quadratic
  theorem from the exact missing scalar-evaluator object
- PMNS stays blocked until the charged-lepton left basis arrives with the same
  family labels and orientation certificate

The outputs stay under `particles/runs/neutrino/` until an actual OPH-derived
family-response and Majorana phase theorem exists.
