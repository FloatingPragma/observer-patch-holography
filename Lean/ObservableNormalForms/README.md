# Observation-Determined Normal Forms: Lean artifact

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
- `MechanismVariants.lean`: domain-neutral variant-indexed trace, protected
  behavior, authority-view, support relation, exact-range encoding, and
  pairwise-versus-family comparison interfaces. `ComparisonPolicy` keeps the
  target, observation, authority view, selected support, cost provenance, and
  family typed together; cost does not enter the behavior order or repair
  width. The exact-range bridge composes actual encoded-collar surjectivity
  with fixed-relation `StrongRepair`. These declarations do not imply
  mechanism reachability, strategy preservation, or full `Repair` width.
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
- `Examples/AmdA2A5Conditional.lean`: an explicitly conditional finite
  mechanism example instantiating Will's shared A2/A5 mechanism-design fixture
  with the OPH calculus. It keeps broad/flash target, raw/class authority,
  outcome/full/hybrid observation, strategy, vacuity, and incomplete-family
  controls as explicit adapters. It supplies separate broad/outcome and
  flash/hybrid policies, an exact-range-to-protected-representative repair
  equivalence, a static fixed-relation missing-support no-go, and one record
  carrying all declared control outcomes simultaneously. `AxiomAudit.lean`
  imports it to print every theorem's dependencies; it remains a separately
  reproducible designated target:
  `lake build ObservableNormalForms.Examples.AmdA2A5Conditional`.
- `AxiomAudit.lean`: theorem-level `#print axioms` receipt for the public
  results and the complete conditional-example theorem set.
- `PROOF_INDEX.md`: paper-label-to-Lean mapping and explicit formalization gaps.
- `SUBMISSION_MANIFEST.md`: exact archive contents and reproduction commands.
- `BUILD_RECEIPT.md`: pinned local parent/standalone build and axiom-audit result.

All theorem-bearing source files in this directory are intended to be
`sorry`-free. The build and `#print axioms` receipt should be regenerated for
the final submission archive after the manuscript theorem numbering freezes.

The conditional module does not claim same-initial executable repair,
all-path or trace completeness, full `Repair` width, a complete support
family, or a family-wide, global, or unique minimum. Its missing-support
result is a static fixed-`R` no-go and has no fabricated exact-range witness.
`FullControlActivation` says only that all declared control outcomes coexist
as fields of one Lean value. Its two comparison-policy fields remain
distinct and do not assert one common target, authority, observation,
support, or family policy.

## Audited edge conditions encoded here

The current manuscript includes the nonemptiness conditions exposed during
the audit. The artifact machine-checks why they matter:

1. The strong-repair projection equivalence needs a nonempty write space.
   `empty_write_space_counterexample` proves the failure without it.
2. The real-valued robust repair margin needs a nonempty relation `R`.
   Mathlib's ordinary `infDist` to the empty set is zero, as proved by
   `empty_relation_repairMargin_zero`. An alternative would be to formulate
   the result with an extended-real distance.

`Refinement.lean` proves the arbitrary-depth metric telescope and the
receipt-composition cores.  The manuscript's family-uniform inverse and
residual moduli are not represented by a separate Lean declaration here.
Finite descriptions of their parametric infinite counterfamilies are instead
recomputed by
`code/consensus/verify_issue_517_proof_obligations.py` in the repository
root: stagewise injectivity without a family inverse modulus, stagewise error
bounds without a family residual modulus, an unbounded family of observation
Lipschitz constants, a nonvanishing tower tail, and a missing Lipschitz
restriction bound, with additional cofinal-limit controls
for completeness, input compatibility, and vanishing solver error.  The
Python artifact checks exact witness constructors and regression instances;
the quantified family argument remains the manuscript proof. This
distinction prevents the executable controls from being reported as a Lean
proof of the infinite-family limit statement.

The issue-517 transactional local-diamond and prepared-acceptance/lock BFT
theorems are outside this Lean package. They have manuscript proofs and finite Python
reference receipts under `code/consensus/`; neither general theorem is claimed
as Lean, TLA+, or protocol-model-checker verified.

The exact product-calculus theorem in this package uses a nonempty finite
index family. The manuscript permits an infinite family only under pairwise
weighted-\(\ell^p\) summability; that generalization is outside this package.
