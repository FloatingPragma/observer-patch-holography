# Proof-to-Paper Index

Mapping between Lean 4 theorems in this project and statements in
*Paradise as Fixed-Point Consensus* (B. Müller, 2026; source in
`paper/paradise_as_fixed_point_consensus.tex`).

## Completion summary

- **Proposition 4.2 (theorem-grade): 0 / 9 obligation families formalised → 0%**
  (2 Prop-4.2 statements + Lyapunov termination/descent +
  OPH confluence + repair completeness + `repair_respects_gauge`
  quotient descent + `gaugeEquiv_equivalence` + `Φ`-formula match +
  quotient/NF construction)
- **Definition 4.1 (Public world): 0 / 11 Lean items built → 0%**
  (3 carrier types + repair-site/local-step package + 1 global repair +
  1 potential + 1 gauge relation + equivalence proof + congruence proof +
  2 quotient/NF construction items)
- Abstract-rewriting skeleton (preliminary): 5 / 5 proofs → 100%
- OPH primitives (declared, sorry-bearing): 0 / 10 discharged → 0%
- #304 boundary-fiber carrier witness (`Rule90.lean`, PR #385): 5 / 5 theorems, sorry-free → 100%
  (first non-degenerate `Hfib` discharge on a linear information-set carrier +
  `H1`–`H3` local-repair no-go; a carrier-level witness only. It does **not** bear
  on the Prop 4.2 / Def 4.1 counts above, which remain 0%)

The headline number is **0% of Proposition 4.2** until the OPH-specific
structures are in place. The 100% skeleton number is preliminary
infrastructure, not progress on Prop 4.2 itself. Counts revised per
math-seat audit (2026-05-19): previous "3+2" undercounted by ≈3×.

## Status legend

- ✅  proven, sorry-free, matches paper statement
- 🟡  proven but generic (does not yet refer to OPH structure)
- 🔸  motivational only (no formal anchor in paper)
- ⬜  declared as `sorry` or thin placeholder / not yet formalised
- ❌  stated, contains unintended `sorry` or `admit`

## Definition 4.1: Public world

> `World = NF(x) / ∼_gauge`, where `NF(x)` is the terminal state reached by
> accepted repair and `∼_gauge` identifies hidden local presentations with
> the same declared observable overlap data.

Paper: `paper/paradise_as_fixed_point_consensus.tex` §4, lines 305–313.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.Records` | `Primitives` | ⬜ | TeX macro line 28; structural content per OPHConsensus. |
| `OPH.Patch` | `Primitives` | ⬜ | TeX macro line 31; structural content per OPHConsensus. |
| `OPH.Obs` | `Primitives` | ⬜ | TeX macro line 30; structural content per OPHConsensus. |
| `OPH.Site` | `Primitives` | ⬜ | Repair-site carrier for local accepted repair steps; structural content per OPHConsensus. |
| `OPH.Repair`, `OPH.localRepair`, `OPH.acceptedStep` | `Primitives` | ⬜ | "Built from local recovery moves" (line 297), composed under asynchronous schedules in OPHConsensus. |
| `OPH.Φ` | `Primitives` | ⬜ | Concrete formula `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))` (line 300). |
| `OPH.gaugeEquiv` (`∼_gauge`) | `Primitives` | ⬜ | Predicate on declared observable overlap data (line 311). |
| `OPH.gaugeEquiv_equivalence` | `Primitives` | ⬜ | `∼_gauge` is an equivalence relation. |
| `OPH.repair_respects_gauge` | `Primitives` | ⬜ | `∼_gauge` is a `Repair`-congruence. Load-bearing for Prop 4.2 sentence 2. |
| `OPH.NF` | (TODO) | ⬜ | Terminal state of accepted repair (local recovery moves). |
| `OPH.World` | (TODO) | ⬜ | Quotient `Records / ∼_gauge` restricted to `NF` representatives. |

## Proposition 4.2: Fixed-point reading of reality

> `World ∈ Fix(Repair)`, `Repair(World) = World`, and (under confluence +
> completeness) terminal public state is independent of update schedule on
> the physical quotient.

Paper: `paper/paradise_as_fixed_point_consensus.tex` §4, lines 321–328.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.Confluence` | `Primitives` | ⬜ | OPH-specific confluence hypothesis (line 326); per OPHConsensus. |
| `OPH.LyapunovDescent` | `Primitives` | ⬜ | Strict descent of accepted repairs; OPHConsensus uses this plus finite state/value range for termination. |
| `OPH.Termination` | `Primitives` | ⬜ | Termination of the accepted-step relation; not implied by repair completeness. |
| `OPH.Completeness` | `Primitives` | ⬜ | Normal forms are exactly consistent states; per OPHConsensus Assumption `ass:complete`. |
| `OPH.world_is_fixedPt` | (TODO) | ⬜ | `Repair(World) = World` on the physical quotient. |
| `OPH.schedule_independence` | (TODO) | ⬜ | Under `Confluence ∧ Completeness`, terminal public state is independent of update schedule on the physical quotient. |

## Abstract-rewriting skeleton (preliminary)

Generic results over `r : X → X → Prop` and `T : X → X`. Verified
sorry-free. These are the abstract layer the OPH
proof instantiates after the concrete repair layer is supplied.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.AbstractRewriting.newman_lemma` | `AbstractRewriting` | 🟡 | Terminating + locally confluent → confluent. Generic. |
| `OPH.AbstractRewriting.unique_normal_form` | `AbstractRewriting` | 🟡 | Confluent → unique normal forms. |
| `OPH.AbstractRewriting.newman_unique_nf` | `AbstractRewriting` | 🟡 | Combined: terminating + locally confluent → unique normal forms. |
| `OPH.AbstractRewriting.deterministic_full` | `AbstractRewriting` | 🟡 | Deterministic op with descent potential → unique fixed-point reached. |
| `OPH.AbstractRewriting.fixedPt_zero_potential` | `AbstractRewriting` | 🔸 | Any fixed point of a repair op has `Φ = 0`. **No paper anchor**: the passage at lines 330-334 is prose, not a labelled corollary; paper makes no formal `Φ(W) = 0` claim. Kept as preliminary ARS decoration only. |

## #304: Boundary-fiber observer uniqueness: Rule 90 carrier witness

> `Hfib`: within a fixed boundary fiber, consistent states are a gauge-singleton
> the hypothesis binder of `boundary_fiber_observer_unique` (`Primitives.lean`).

Instantiated on a linear (Rule 90) two-patch carrier, the first case where the
`#304` `Hfib` binder is non-degenerate (a proper, failable information set + a
non-trivial gauge), unlike `demoCarrier`'s seed tautology. A carrier-level
witness, **not** a Prop 4.2 / Def 4.1 item. Module:
`ObserverPatchHolography/Rule90.lean` (PR #385). Sorry-free, standard axioms
only (`propext`, `Classical.choice`, `Quot.sound`; no `native_decide`).

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.rule90_Hfib_good` | `Rule90` | ✅ | Reading bottom cells {0,1} discharges `Hfib` on the consistent set, instantiating the `boundary_fiber_observer_unique` binder. |
| `OPH.rule90_Hfib_bad_fails` | `Rule90` | ✅ | Reading {0,2} fails: same carrier, coarser boundary. `Hfib` is about *which* cells `B` reads, not how many. |
| `OPH.rule90_gauge_nontrivial` | `Rule90` | ✅ | Seeds (0,0,0), (1,0,1) observably identical: the gauge contains `ker(Rule90)` (kernel-pair exhibit). |
| `OPH.rule90_no_frustrationFree_repair` | `Rule90` | ✅ | No `H1`–`H3` local repair exists on this carrier. Scope: transactional/multi-patch repair and relaxed `H2` are **not** ruled out. |
| `OPH.rule90t_outer_eq` | `Rule90` | ✅ | Helper: every Rule 90 image has equal outer cells. |

## Gap analysis: skeleton → theorem-grade Prop 4.2

(Revised per math-seat audit 2026-05-19.)

| Skeleton has | Prop 4.2 needs | Paper anchor |
|---|---|---|
| `r : X → X → Prop` (opaque) | `Repair` operating on `Records`, inducing rewriting on `World` representatives | lines 30, 297 |
| `X` (opaque) | `Records` carrier type (TeX macro, structural content per OPHConsensus) | lines 28, 1615 |
| `IsNormalForm` (opaque) | `NF` as terminal state of accepted repair built from local recovery moves, with local steps executed under asynchronous schedules in OPHConsensus | line 297; OPHConsensus §3 |
| `Φ : X → NNReal` (opaque) | Concrete `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))` on patch-overlap data | line 300 |
| (none) | `∼_gauge` is an equivalence relation | line 311 |
| (none) | **`∼_gauge` is a `Repair`-congruence**; descent to the physical quotient depends on this | line 327 ("on the physical quotient") |
| (none) | Quotient construction `World = NF(x) / ∼_gauge` | lines 305–313 |
| `Terminating r` (assumption) | OPH Lyapunov descent plus finite patch-net/value-set control gives termination of `acceptedStep`; repair completeness is a separate normal-form/consistency condition | OPHConsensus Prop. `lyapunov-termination`, Assumption `complete` |
| `LocallyConfluent r` (assumption) | OPH `Confluence` predicate (structurally defined per OPHConsensus) implies local confluence on patch repairs | line 326 |
| Confluence on `r` | Confluence stated **directly on the physical quotient**, not lifted from below. Prop 4.2 sentence 2 imposes conditions on the quotient. | lines 326–327 |

Closing each row of this table is the work that takes the skeleton to a
theorem-grade Prop 4.2 statement. Several rows depend on the companion
paper *Reality as a Consensus Protocol* (OPHConsensus); the target is
paper-incomplete as well as Lean-incomplete.

## Update protocol

When adding a Lean statement that targets a paper item:

1. Add a row in the appropriate section above.
2. Cite the paper line range (use
   `paper/paradise_as_fixed_point_consensus.tex` line numbers).
3. Set status (✅ / 🟡 / 🔸 / ⬜ / ❌).
4. Update the completion summary.
5. Tag the PR with the affected proposition/definition for the cross-audit
   pass (Dula via Grok, Ben via ChatGPT-Pro). Bernhard merges.
