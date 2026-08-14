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

Python `3.12` or newer is required only for the deterministic finite-model
explorer described below.

The directory can be zipped as-is for a submission artifact. Generated `.lake`
content and build logs are intentionally excluded.

## 60-second review route

1. Read the protected-obstruction result and scope boundary below.
2. Inspect the protected-obstruction table in [`PROOF_INDEX.md`](PROOF_INDEX.md).
3. Review `ObservableNormalForms/ProtectedObstructions.lean`, then
   `ObservableNormalForms/Examples/ProtectedObstructions.lean` for the
   fine/coarse TwoBit separation, the full/missing fixed-relation boundary,
   the proper-target closed-trap and ambiguous-endpoint chains, and the
   common-carrier strict reversal.
4. Inspect `ObservableNormalForms/AxiomAudit.lean`; run `lake build`, the
   designated conditional target, and
   `python tools/verify_protected_obstruction_models.py`.
5. Review the manuscript and claim-registry diff last. Build evidence supports
   reproducibility; it does not replace theorem or assumption review.

## Protected-obstruction result and scope

For a declared finite scheduler, the development separates four nested
questions for each protected behavior: whether its consistent fiber exists,
whether some declared active source reaches it with positive first-hit
probability, whether every active source hits almost surely, and whether the
first-hit endpoint class is unique modulo the declared silent equivalence.
Their complementary cuts form an exact obstruction profile. Explicit
supported-trace adapters recover the four existing `RealizedBehavior` and
`BehaviorCut` coordinates of `MechanismVariants.lean`, and an exact
finite-state morphism transports the profile under the full source, target,
kernel-lumping, and endpoint-quotient hypotheses; observation commutation is
packaged separately in `FrozenMorphism` and is not consumed by the layer
transport.

Finite-state first-hit expansion, closed-class reachability, Bellman
least-fixed-point reasoning, and strong lumpability are standard Markov-chain
ingredients. The contribution here is their source-quantified assembly around
protected observations, the exact OPH behavior-cut adapter, and the verified
native activations and counterexamples. No rate, expected hitting time,
mixing, infinite-state or refinement-tower theorem, scalar ranking, deployed
implementation refinement, or arbitrary-scheduler invariance is claimed.

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
- `ProtectedObstructions.lean`: canonical finite-state first-hit and transport
  development. It builds exact-time mass and fixed-horizon cumulative endpoint
  mass from `FiniteMarkovKernel.pathWeight`, takes bounded monotone limits,
  derives positive-path and no-reachable-closed-trap characterizations, and
  proves the least nonnegative Bellman fixed-point property. Structural
  morphisms carry protected/initial/target reflection and coverage, kernel
  lumping, and quotient exactness; fiber nonemptiness, positive/almost-sure/
  selection exactness, endpoint pushforward, and L0--L3 exactness are derived
  theorems rather than premise fields. The file also contains the T0--T8 public
  adapter, common-carrier family completion, and full-support initial-law
  corollaries.
- `Examples/ProtectedObstructions.lean`: nonidentity transport and finite
  load-bearing fixtures for common-carrier strict reversal, an empty target
  fiber, fine/coarse TwoBit, full/missing fixed relations, and proper-target
  closed-trap and ambiguous-endpoint chains. The proper-target models keep
  equality at the source distinct from positive transition support.
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
- `tools/verify_protected_obstruction_models.py`: deterministic finite-model
  explorer for the protected-obstruction fixtures and their minimum state
  counts.

All theorem-bearing source files in this directory are intended to be
`sorry`-free. The build and `#print axioms` receipt should be regenerated for
the final submission archive after the manuscript theorem numbering freezes.

## Protected-obstruction verification boundary

Lean is canonical for the stochastic protected-obstruction implementation:
the transition kernel, path-weight first-hit law, scheduler correspondence,
behavior-cut adapter, and quantitative endpoint pushforward are all derived
here, and the axiom audit covers every public theorem.

The Bellman fixed-point characterization is a theorem derived from the
canonical path law, not the definition of the law. In particular, the artifact
does not claim uniqueness among all nonnegative Bellman fixed points when a
reachable closed non-target trap exists.

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
