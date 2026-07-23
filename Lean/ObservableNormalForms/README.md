# Observation-Determined Normal Forms — Lean artifact

This directory is a standalone, generic Lean 4 formalization accompanying
*Observation-Determined Normal Forms: Stability, Obstructions, and Refinement
in Constraint and Rewrite Systems*. Its theorem modules import only Mathlib
and other modules in this artifact.

The standalone project is pinned by:

- `lean-toolchain`: Lean `v4.29.1`
- `lakefile.lean`: Mathlib `v4.29.1`
- `lake-manifest.json`: pinned, checksummed dependency revisions
- `ObservableNormalForms.lean`: library root

Build this directory directly with:

```sh
lake build
```

The directory can be zipped as-is for a submission artifact. Generated `.lake`
content and build logs are intentionally excluded.

## Contents

- `Exact.lean`: proof-carrying partial normalizer, universal property,
  reachable-fiber theorem, and the corrected empty/singleton/confluence
  alternative.
- `ObserverConfluence.lean`: exact equivalence between boundary identification
  on consistent states and cross-source normal-endpoint uniqueness modulo an
  arbitrary relation, plus weak-normalization existence and a complete
  two-bit positive/separation example.
- `Stability.lean`: the heterogeneous two-output estimate, symmetric form,
  approximate schedule independence, and sensor-enrichment certificate.
- `Refinement.lean`: one-step approximate naturality, exact naturality,
  arbitrary-depth metric telescope, same-level agreement, fine-to-coarse
  comparison from an exact tower receipt, and the anchored/nested metric cores.
- `Repair.lean`: strong-repair projection criterion, no-repair certificate,
  robust margin, and machine-checked counterexamples for two missing
  nonemptiness hypotheses found during the audit.
- `Functional.lean`: synchronous ranked-dependency settling, fixed-point, and
  uniqueness theorems.
- `Stochastic.lean`: finite Markov affine-drift iteration and a one-time
  finite-state tail bound.
- `ConditionalResampling.lean`: finite weighted conditional-resampling
  kernels, fixed-point characterization, idempotence, weighted
  self-adjointness, Pythagorean energy identity, contraction, and the exact
  fiber-support/equal-row/detailed-balance matrix-recognition converse.
- `Examples/Rule90.lean`: standalone width-three kernel, image, readout, and
  reverse-repair obstruction proofs.
- `AxiomAudit.lean`: theorem-level `#print axioms` receipt.
- `PROOF_INDEX.md`: paper-label-to-Lean mapping and explicit formalization gaps.
- `SUBMISSION_MANIFEST.md`: exact archive contents and reproduction commands.
- `BUILD_RECEIPT.md`: pinned local parent/standalone build and axiom-audit result.

All theorem-bearing source files in this directory are intended to be
`sorry`-free. The build and `#print axioms` receipt should be regenerated for
the final submission archive after the manuscript theorem numbering freezes.

## Audited edge conditions encoded here

The current manuscript includes the nonemptiness conditions exposed during
the audit. The artifact machine-checks why they matter:

1. The strong-repair projection equivalence needs a nonempty write space.
   `empty_write_space_counterexample` proves the failure without it.
2. The real-valued robust repair margin needs a nonempty relation `R`.
   Mathlib's ordinary `infDist` to the empty set is zero, as proved by
   `empty_relation_repairMargin_zero`. An alternative would be to formulate
   the result with an extended-real distance.

The current exact product-calculus theorem likewise uses a nonempty finite
index family. That theorem is not yet formalized here.
