# Proof-to-Paper Index

Mapping between Lean 4 theorems in this project and statements in
*Paradise as Fixed-Point Consensus* (B. Müller, 2026; source in
`paper/paradise_as_fixed_point_consensus.tex`).

## Completion summary

- **Proposition 4.2 (theorem-grade): 0 / 7 obligations formalised → 0%**
  (2 Prop-4.2 statements + 2 OPH-side hypothesis predicates +
  `repair_respects_gauge` congruence + `gaugeEquiv_equivalence` +
  `Φ`-formula match)
- **Definition 4.1 (Public world): 0 / 10 Lean items built → 0%**
  (3 carrier types + 1 repair + 1 potential + 1 gauge relation +
  equivalence proof + congruence proof + 2 quotient/NF construction items)
- Abstract-rewriting skeleton (preliminary): 5 / 5 proofs → 100%
- OPH primitives (declared, sorry-bearing): 0 / 10 discharged → 0%

The headline number is **0% of Proposition 4.2** until the OPH-specific
structures are in place. The 100% skeleton number is preliminary
infrastructure, not progress on Prop 4.2 itself. Counts revised per
math-seat audit (2026-05-19): previous "3+2" undercounted by ≈3×.

## Status legend

- ✅  proven, sorry-free, matches paper statement
- 🟡  proven but generic (does not yet refer to OPH structure)
- 🔸  motivational only (no formal anchor in paper)
- ⬜  declared as `sorry` / not yet formalised
- ❌  stated, contains unintended `sorry` or `admit`

## Definition 4.1 — Public world

> `World = NF(x) / ∼_gauge`, where `NF(x)` is the terminal state reached by
> accepted repair and `∼_gauge` identifies hidden local presentations with
> the same declared observable overlap data.

Paper: `paper/paradise_as_fixed_point_consensus.tex` §4, lines 305–313.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.Records` | `Primitives` | ⬜ | TeX macro line 28; structural content per OPHConsensus. |
| `OPH.Patch` | `Primitives` | ⬜ | TeX macro line 31; structural content per OPHConsensus. |
| `OPH.Obs` | `Primitives` | ⬜ | TeX macro line 30; structural content per OPHConsensus. |
| `OPH.Repair` | `Primitives` | ⬜ | "Built from local recovery moves" (line 297). Not *asynchronous* — paper says local. |
| `OPH.Φ` | `Primitives` | ⬜ | Concrete formula `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))` (line 300). |
| `OPH.gaugeEquiv` (`∼_gauge`) | `Primitives` | ⬜ | Predicate on declared observable overlap data (line 311). |
| `OPH.gaugeEquiv_equivalence` | `Primitives` | ⬜ | `∼_gauge` is an equivalence relation. |
| `OPH.repair_respects_gauge` | `Primitives` | ⬜ | `∼_gauge` is a `Repair`-congruence. Load-bearing for Prop 4.2 sentence 2. |
| `OPH.NF` | (TODO) | ⬜ | Terminal state of accepted repair (local recovery moves). |
| `OPH.World` | (TODO) | ⬜ | Quotient `Records / ∼_gauge` restricted to `NF` representatives. |

## Proposition 4.2 — Fixed-point reading of reality

> `World ∈ Fix(Repair)`, `Repair(World) = World`, and (under confluence +
> completeness) terminal public state is independent of update schedule on
> the physical quotient.

Paper: `paper/paradise_as_fixed_point_consensus.tex` §4, lines 321–328.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.Confluence` | `Primitives` | ⬜ | OPH-specific confluence hypothesis (line 326); per OPHConsensus. |
| `OPH.Completeness` | `Primitives` | ⬜ | OPH-specific completeness hypothesis (line 326); per OPHConsensus. |
| `OPH.world_is_fixedPt` | (TODO) | ⬜ | `Repair(World) = World` on the physical quotient. |
| `OPH.schedule_independence` | (TODO) | ⬜ | Under `Confluence ∧ Completeness`, terminal public state is independent of update schedule on the physical quotient. |

## Abstract-rewriting skeleton (preliminary)

Generic results over `r : X → X → Prop` and `T : X → X`. Verified
sorry-free. These are not Prop 4.2 — they are the abstract layer the OPH
proof will eventually instantiate.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.AbstractRewriting.newman_lemma` | `AbstractRewriting` | 🟡 | Terminating + locally confluent → confluent. Generic. |
| `OPH.AbstractRewriting.unique_normal_form` | `AbstractRewriting` | 🟡 | Confluent → unique normal forms. |
| `OPH.AbstractRewriting.newman_unique_nf` | `AbstractRewriting` | 🟡 | Combined: terminating + locally confluent → unique normal forms. |
| `OPH.AbstractRewriting.deterministic_full` | `AbstractRewriting` | 🟡 | Deterministic op with descent potential → unique fixed-point reached. |
| `OPH.AbstractRewriting.fixedPt_zero_potential` | `AbstractRewriting` | 🔸 | Any fixed point of a repair op has `Φ = 0`. **No paper anchor** — the passage at lines 330–334 is prose, not a labelled corollary; paper makes no formal `Φ(W) = 0` claim. Kept as preliminary ARS decoration only. |

## Gap analysis: skeleton → theorem-grade Prop 4.2

(Revised per math-seat audit 2026-05-19.)

| Skeleton has | Prop 4.2 needs | Paper anchor |
|---|---|---|
| `r : X → X → Prop` (opaque) | `Repair` operating on `Records`, inducing rewriting on `World` representatives | lines 30, 297 |
| `X` (opaque) | `Records` carrier type (TeX macro, structural content per OPHConsensus) | lines 28, 1615 |
| `IsNormalForm` (opaque) | `NF` as terminal state of accepted repair (built from *local* recovery moves — paper does NOT say *asynchronous*) | line 297 |
| `Φ : X → NNReal` (opaque) | Concrete `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))` on patch-overlap data | line 300 |
| (none) | `∼_gauge` is an equivalence relation | line 311 |
| (none) | **`∼_gauge` is a `Repair`-congruence** — descent to the physical quotient depends on this | line 327 ("on the physical quotient") |
| (none) | Quotient construction `World = NF(x) / ∼_gauge` | lines 305–313 |
| `Terminating r` (assumption) | OPH `Completeness` predicate (structurally defined per OPHConsensus) implies termination on `Records` | line 326 |
| `LocallyConfluent r` (assumption) | OPH `Confluence` predicate (structurally defined per OPHConsensus) implies local confluence on patch repairs | line 326 |
| Confluence on `r` | Confluence stated **directly on the physical quotient**, not lifted from below. Prop 4.2 sentence 2 imposes conditions on the quotient. | lines 326–327 |

Closing each row of this table is the work that takes the skeleton to a
theorem-grade Prop 4.2 statement. Several rows depend on the companion
paper *Reality as a Consensus Protocol* (OPHConsensus) — the target is
paper-incomplete, not just Lean-incomplete.

## Update protocol

When adding a Lean statement that targets a paper item:

1. Add a row in the appropriate section above.
2. Cite the paper line range (use
   `paper/paradise_as_fixed_point_consensus.tex` line numbers).
3. Set status (✅ / 🟡 / 🔸 / ⬜ / ❌).
4. Update the completion summary.
5. Tag the PR with the affected proposition/definition for the cross-audit
   pass (Dula via Grok, Ben via ChatGPT-Pro). Bernhard merges.
