# Proof-to-Paper Index

For the corrected typed Einstein composition theorem (GitHub #578), see the
current statement-by-statement and premise audit in
[`EINSTEIN_BRANCH_INDEX.md`](EINSTEIN_BRANCH_INDEX.md). That development is
sorry-free and formalizes the algebraic/compositional implication while
leaving continuum, asymptotic, physical-identification, and tower-nonemptiness
premises explicit.

The cross-paper confluence, collar, icosahedral, and event-geometry boundaries
are mapped in [`BRIDGE_BOUNDARY_INDEX.md`](BRIDGE_BOUNDARY_INDEX.md). The
finite counterexamples are imported by the umbrella source library, and the
Screen algebra modules are part of the default Lake build.

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
- Part-A coupling algebra (`BridgeEquivalence.lean`, `CapacityFixedPoint.lean`,
  `SeedPi.lean`): 13 / 13 lemmas, sorry-free → 100%
  (algebraic layer of the coupling theorem only; no physical-derivation
  content, the physical identities I1/I2 are outside the formalised set;
  does **not** bear on the Prop 4.2 / Def 4.1 counts)
- Physical A5-to-SM reduct boundary
  (`Screen/PhysicalA5ForcingNoGo.lean`): 6 / 6 declarations, sorry-free.
  The finite theorems prove non-reconstruction of defect, current, and matter
  completions from their stated forgetful reducts, plus the settlement and
  pairing controls. They do not prove nonexistence of a selector from a richer
  operational source packet and do not instantiate any physical producer.
- Icosahedral port-frame Gram identity (`Screen/PortFrameGram.lean`, #568):
  12 / 12 declarations, sorry-free, kernel `decide` only, standard axioms
  (`propext`, `Classical.choice`, `Quot.sound`). Proves on the explicit
  twelve-port icosahedral incidence: 5-regularity, symmetric irreflexive
  adjacency, fixed-point-free antipodal involution, the distance profile
  witnesses (common neighbor for every distance-two pair, no common neighbor
  across an antipodal pair), and the exact Gram identities `G² = 4·G` and
  `tr G = 12` in the integer-scaled `ℤ(√5)` normalization `5·G`. The
  rank-three conclusion (idempotent `G/4` with trace 3) is the one on-paper
  step. It does **not** derive the icosahedral carrier from OPH axioms,
  construct the `A5` action, or bear on the Prop 4.2 / Def 4.1 counts.
- Icosahedral port action (`Screen/A5PortAction.lean`, #568): 9 / 9
  declarations, sorry-free, kernel `decide`, axioms `propext`/`Quot.sound`
  only. Exhibits the sixty proper icosahedral rotations as explicit port
  permutations and proves the listed set is a group of order sixty (closed,
  identity, two-sided inverses in the list, pairwise distinct, each row a
  genuine permutation), that every element preserves the port incidence and
  commutes with the antipode, and that the action is transitive on ports.
  The identification of this order-sixty incidence-automorphism group with
  the abstract `A5` (simplicity or classification) stays on paper; no
  physical current or descent content.
- Universal cap coupling from icosahedral symmetry
  (`Screen/A5CouplingSymmetry.lean`, #576/#568): 9 / 9 declarations,
  sorry-free, standard axioms. Finite facts by kernel `decide` (caps map to
  caps under all sixty rotations, the composition law, the right-translation
  bijection of the group list, transitivity on ports); the quantified
  statements are symbolic over arbitrary rational port data: `A5`-invariant
  readouts have equal cap sums at every port, group-averaged cap sums are
  port-independent for every data assignment, and hence the per-cap ratio of
  any two averaged readouts (entropy over stress flux) is universal with
  zero spread. This converts the Einstein-branch universality clause from a
  measurement into a theorem whose one remaining hypothesis is
  `A5`-equivariance of the source law; checking that hypothesis on an
  implemented law is a finite simulator check. No physical coupling value or
  promotion is implied.
- A5 group/trichotomy receipt corpus (`Screen/A5OPH.lean`, external
  cross-audit contribution, merged 2026-07-23): 18 / 18 audited
  declarations, sorry-free, standard axioms (`propext`, `Classical.choice`,
  `Quot.sound`), no `native_decide`. Six namespaces: `A5Order6` (order-six
  subgroups of `A5`, nonabelian hence `≅ S3`), `A5NoZ6` (no element of
  order six, hence no `ℤ/6`; the odd order-six `S5` element exhibited as
  the sharpness control), `OPHGap` (negative control: a `G`-equivariant
  linear isomorphism does not determine a Lie bracket — guards the
  module-decomposition step against over-reading; no corpus surface makes
  the guarded claim), `OPHSelection` (finite core of the compact-Lie
  trichotomy: every `A5`-action on at most four objects is trivial;
  `11 = 3+8` and `12 = 3+3+3+3` uniquely over the compact-simple dimension
  list `{3,8,10}`), `OPHInner` (the centre is characteristic; the
  inner-action reduction (H2) ⇒ (H*) stated modulo two declared classical
  inputs — exponential surjectivity and `Ad(G) = Inn` on compact connected
  groups; `A5` embeds in no group with a unique involution, hence
  `A5 ⊄ SU(2)`), `OPHTrichotomy` (no compact semisimple algebra in
  dimensions 1, 2, 4, 5, 7; the paper's noncentral-quintet witness
  `[iS, iT] = -2(E₁₂ - E₂₁) ≠ 0`; the screen gluing class
  `Λ₊/(Λ₁ ⊕ Λ₅) ≅ ℤ/6` with proper-rotation invariance and antipodal sign
  reversal — the bare quotient isomorphism is also in `Z6Exact.lean`; this
  module adds the two invariance clauses and the trichotomy context).
  Provenance: authored in the external cross-audit lane (Aristotle build,
  Lean/Mathlib v4.28.0), rebuilt verbatim under this project's toolchain
  (v4.29.1); rebuild record in `DULA/VERIFICATION.md` at the meta-repo
  root. Declared classical inputs stay on paper: compact-simple
  classification, reductive decomposition, exponential surjectivity, the
  dimension-six branch module argument (#604), and the `ℚ(√5)` rationality
  lemma (#605). No physical receipt is formalized. Does **not** bear on the
  Prop 4.2 / Def 4.1 counts.
- #304 boundary-fiber carrier witness (`Rule90.lean`, PR #385): 5 / 5 theorems, sorry-free → 100%
  (first non-degenerate `Hfib` discharge on a linear information-set carrier +
  `H1`–`H3` local-repair no-go; a carrier-level witness only. It does **not** bear
  on the Prop 4.2 / Def 4.1 counts above, which remain 0%)
- #304 application theorem (`BoundaryFiber.lean`): 8 / 8 audited declarations,
  sorry-free → 100% (the declared boundary map `B_OPH` defined concretely as
  the root-packet readback on the verified rooted-tree packet-net domain, its
  injectivity modulo gauge proved class-wide, and the endpoint payoff composed
  through the neutral equivalence; the class restriction is the named
  named premise for nets outside the declared domain. It does **not** bear on
  the Prop 4.2 / Def 4.1 counts)
- Finite event algebras (`EventAlgebra` library, journal-neutral bundle):
  64 / 64 audited declarations, sorry-free → 100% (standard axioms only;
  `chsh_mul_self` needs only `propext` + `Quot.sound`). A self-contained
  neutral-vocabulary development : it does **not** bear on the
  Prop 4.2 / Def 4.1 counts. See "Finite event algebras" below.

The headline number is **0% of Proposition 4.2** until the OPH-specific
structures are in place. The 100% skeleton number is preliminary
infrastructure, not progress on Prop 4.2 itself. Counts revised per
math-seat audit (2026-05-19): previous "3+2" undercounted by ≈3×.

> **Addendum (2026-07-11) : the three declared `Primitives.lean` admissions
> are discharged.** `localRepair` is a constructed single-site
> transactional recovery move (fires iff an incident overlap is broken and
> the site can satisfy all its overlaps at once, repaired state chosen from
> declared overlap data), `Repair` is a choice-canonical asynchronous
> schedule composed to a normal form (well-founded by broken-edge count),
> and `repair_respects_gauge` (Prop 4.2 sentence 2 congruence) is a proven
> theorem, all pinned at `[propext, Classical.choice, Quot.sound]`. The
> file's own `LyapunovDescent` and `Termination` obligations are also
> discharged for the constructed operator (`lyapunovDescent_holds`,
> `termination_holds`), with normal-form/reachability receipts
> (`Repair_normalForm`, `Repair_reachable`). NOT discharged: OPH
> `Confluence` (false in general : `demoCarrier_not_confluent`),
> `Completeness` for the constructed operator (holds exactly on
> frustration-free dynamics via the conditional `H1`–`H3` layer), and the
> `NF`/`World` quotient construction. The tables and counts below predate
> this discharge and need a math-seat re-audit before renumbering; rows for
> `OPH.Repair`/`OPH.localRepair`/`OPH.repair_respects_gauge` should read
> as constructed/proven rather than ⬜-declared.

> **Addendum (2026-07-18) : issue #544 collar chain (layer separation +
> state-side no-gos), 41 / 41 audited declarations, sorry-free.**
> `CollarClause.lean`/`CollarLayer.lean` prove the layer separation (the
> overlap-consistency repair layer factors through the realized constraint
> family; the collar clause is a declared input, not a theorem).
> `CollarStates.lean` (T0) proves `stateSide_axioms_do_not_force`:
> the state-side axioms (density matrices, Gibbs states of the retained
> family, Umegaki relative entropy, family Klein inequality) do not force
> the clause. `CollarStatesT1.lean` (T1) proves `Eflux_does_not_force` and
> `EfluxChannel_deselects_XXC`: the flux conditional expectation exists,
> kills the cross-cut coupling without forcing the clause
> (deselection ≠ exclusion). `CollarModularT2.lean` (T2) proves
> `naive_modular_recast_does_not_exclude` and
> `centralizer_diagonal_strictly_contains_flux`: the naive modular recast
> is vacuous and the corrected recast buys only the diagonal clause.
> `CollarStatesBridge.lean` records the reuse receipts against the
> `EventAlgebra` library (`fluxPartition`,
> `collarReindex_EfluxL_eq_compressed_pinching`,
> `EfluxL_absorbs_pinching`): `EfluxL` is the blockwise normalised-trace
> compression of `pinchingExpectation` : strictly finer than any pinching,
> hence retained rather than replaced. All 41 audited declarations pin at
> `[propext, Classical.choice, Quot.sound]`; none bear on the
> Prop 4.2 / Def 4.1 counts.

## Status legend

- ✅  proven, sorry-free, matches paper statement
- 🟡  generic proof without OPH-specific structure
- 🔸  motivational only (no formal anchor in paper)
- ⬜  declared as `sorry` or thin placeholder / unformalised
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

## Part-A coupling algebra (algebraic layer)

Formalises the ALGEBRAIC layer of the OPH coupling theorem: the bridge
count/tick equivalence, the capacity fixed-point uniqueness schema, and the
CAP-P seed statement. No physical-derivation content is formalised; the
physical identities I1/I2 (which give the two sides of the bridge their
physical readings) are outside the formalised set. The numeric interval
enclosures stay in the Python certificates (`code/capacity_readback/`,
`code/P_derivation/`); no floating-point numerics enter Lean. All 13 lemmas
are sorry-free with standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`). Not a Prop 4.2 / Def 4.1 item.

Modules: `ObserverPatchHolography/BridgeEquivalence.lean`,
`ObserverPatchHolography/CapacityFixedPoint.lean`,
`ObserverPatchHolography/SeedPi.lean`. Status marks refer to the
coupling-algebra statements, not to *Paradise* anchors.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.BridgeEquivalence.six_eq_mRep_div_betaEW` | `BridgeEquivalence` | ✅ | The literal `6` in the bridge exponent is `m_rep / β_EW` with `m_rep = 24`, `β_EW = 4`. |
| `OPH.BridgeEquivalence.mRep_eq_six_mul_betaEW` | `BridgeEquivalence` | ✅ | The literal `24` decomposes as `6 · β_EW`. |
| `OPH.BridgeEquivalence.bridge_equivalence_beta` | `BridgeEquivalence` | ✅ | General-β bridge: for `β > 0`, `P > 0`, `αU > 0`, `N > π`: `N = π·exp((m_rep/β)·π/(P·αU)) ↔ tickProjection αU N = β·P`. |
| `OPH.BridgeEquivalence.bridge_equivalence` | `BridgeEquivalence` | ✅ | Literal form: `N = π·exp(6π/(P·αU)) ↔ (24π)/(αU·log(N/π)) = 4·P`. |
| `OPH.BridgeEquivalence.bridge_equivalence_tick` | `BridgeEquivalence` | ✅ | Tick corollary at `β = 4` with `24 = 6·β` carried as an explicit hypothesis. |
| `OPH.CapacityFixedPoint.averagingMap_isFixedPt_iff` | `CapacityFixedPoint` | ✅ | For `λ ≠ 0`, `x` is a fixed point of `(1−λ)·x + λ·K` iff `x = K`. |
| `OPH.CapacityFixedPoint.averagingMap_unique_fixedPt` | `CapacityFixedPoint` | ✅ | For `0 < λ ≤ 1` the averaging map has exactly one fixed point, `K`. |
| `OPH.CapacityFixedPoint.banach_unique_fixedPt` | `CapacityFixedPoint` | ✅ | `ContractingWith` map on a nonempty complete metric space has exactly one fixed point; wraps Mathlib `ContractingWith.fixedPoint`. |
| `OPH.CapacityFixedPoint.capacity_Icc_unique_fixedPt` | `CapacityFixedPoint` | ✅ | Lipschitz self-map of `Icc a b` with constant `K < 1` has exactly one fixed point in `Icc a b` (interval form used by the capacity certificates). |
| `OPH.SeedPi.capReadback_log_coords` | `SeedPi` | ✅ | In log coordinates `y = log(N/π)` the CAP-P readback `π·(N/π)^s` is the linear map `y ↦ s·y`. |
| `OPH.SeedPi.capReadback_log_coords_averaging` | `SeedPi` | ✅ | The log-coordinate readback equals `averagingMap (1−s) 0`, tying CAP-P to the uniqueness schema. |
| `OPH.SeedPi.capReadback_fixedPt_iff` | `SeedPi` | ✅ | CAP-P seed: for `s ≠ 1`, `N > 0`: `π·(N/π)^s = N ↔ N = π` (map shape from `F_candidate_capP.py`). |
| `OPH.SeedPi.linear_branch_no_positive_fixedPt` | `SeedPi` | ✅ | Additive CAP-P branches `F(N) = s·N`, `s < 1`, have no positive fixed point. |

## #304: Generic theorem and concrete bridge

The standalone proof package at `ObservableNormalForms/` proves the
exact substrate-neutral equivalence

```text
boundary identification on consistent states modulo E
  ↔ cross-source normal-endpoint uniqueness modulo E.
```

`OPH.boundaryIdentifiesModulo_iff_observerEndpointUniqueModuloLR` in
`ObserverPatchHolography/Bridges/ObservableNormalForms.lean`
machine-checks that Jonathan's H1--H3 completeness theorem and boundary
preservation instantiate that equivalence for `acceptedStepLR`.

This completes the generic observer-confluence interface.  The application
obligation it leaves open (a declared physical boundary map `B` and a proof
that any two concrete consistent records with the same `B` value are
`gaugeEquiv`) is discharged on the declared domain by the `BoundaryFiber`
module below.

## #304: Declared boundary map and application theorem (`BoundaryFiber`)

The concrete #304 gate on the verified rooted-tree packet-net domain
(*Reality as a Consensus Protocol*, Definition `def:tree-packet-domain`):
`TreePacketNet` realises the domain as an `OPHCarrier` family, `B_OPH` is the
root-packet readback, and injectivity modulo gauge holds class-wide. Module:
`ObserverPatchHolography/BoundaryFiber.lean`. Sorry-free, standard
axioms only. Companion note: `BOUNDARY_FIBER_APPLICATION.md`. Outside the
declared class the identifiability premise is named and can fail
(`demoCarrier_Hfib_fails`, `rule90_Hfib_bad_fails`). Not a Prop 4.2 /
Def 4.1 item.

| Lean name | Module | Status | Notes |
|---|---|---|---|
| `OPH.TreePacketNet` / `TreePacketNet.carrier` | `BoundaryFiber` | ✅ | The declared domain: finite rooted tree, per-patch state `A × K_i`, interfaces read the packet component only; hidden labels are the declared redundancy data. |
| `OPH.TreePacketNet.BOPH` | `BoundaryFiber` | ✅ | The declared physical boundary/sector map: root-packet readback `x ↦ (x root).1` (the paper's "root-packet map" taxonomy entry). |
| `OPH.TreePacketNet.packet_eq_root_of_edgeConsistent` | `BoundaryFiber` | ✅ | Unique continuation: consistency propagates the root packet through the bulk (strong induction on tree depth). |
| `OPH.TreePacketNet.BOPH_injective_modulo_gauge` | `BoundaryFiber` | ✅ | **The #304 gate**: `u,v ∈ C_OPH` and `B_OPH(u)=B_OPH(v)` imply `u ~gauge v`, for every net in the class. |
| `OPH.TreePacketNet.BOPH_boundaryIdentifiesModulo` | `BoundaryFiber` | ✅ | The gate restated as the neutral `BoundaryIdentifiesModulo` premise. |
| `OPH.TreePacketNet.treeStep_observationPreserving` | `BoundaryFiber` | ✅ | The paper's tree repair preserves `B_OPH` (`ObservationPreserving` hypothesis). |
| `OPH.TreePacketNet.treeStep_completeFor` | `BoundaryFiber` | ✅ | Tree-repair normal forms are exactly the consistent records (`CompleteFor` hypothesis). |
| `OPH.TreePacketNet.BOPH_observerEndpointUnique` | `BoundaryFiber` | ✅ | Endpoint payoff: equal `B_OPH` readback forces gauge-equal normal forms under any maximal schedule; no confluence input. |
| `OPH.fourVertexNet_gauge_nontrivial` | `BoundaryFiber` | ✅ | The paper's exported four-vertex `ℤ₃`/`ℤ₂` instance: distinct same-boundary consistent records that are `gaugeEquiv`, so the gauge quotient is proper. |
| `OPH.fourVHiddenReadback_not_identifying` | `BoundaryFiber` | ✅ | Failure witness on the verified domain: a readback missing the protected packet does not identify. |

## Finite event algebras (`EventAlgebra` library)

A self-contained, sorry-free development of finite-dimensional quantum
event algebras over `Matrix (Fin n) (Fin n) ℂ`, built for journal
submission. **These modules are OPH-vocabulary-free by design**: namespace
`EventAlgebra`, imports Mathlib only, no observer/patch/screen/record
vocabulary anywhere, and no imports from the rest of this repository. Every
lemma's doc comment carries a scope tag : **algebra-only** (pure
`*`-algebra content) or **consumes a tracial state** (content passing
through the trace pairing `(ρ, P) ↦ Tr(ρ P)`), so the split between the
algebraic layer and the state layer is machine-visible. Modules live in
`EventAlgebra/` with umbrella `EventAlgebra.lean`; Mathlib
friction log in `EventAlgebra/MATHLIB_NOTES.md`. All 64 audited
declarations report `[propext, Classical.choice, Quot.sound]` (the ring
identity `chsh_mul_self` even avoids `Classical.choice`).

Definitions: `IsEvent` (Hermitian idempotent), `IsState` (PSD, trace one),
`bornWeight` (`Tr(ρ P)`), `luedersUpdate`, `certainStates`,
`CenterPartition`, `centerExpectation`, `expectation`.

`Basic.lean` : events and Born weights (25 lemmas):

| Lean name | One-line statement |
|---|---|
| `isEvent_zero` | `0` is an event (the impossible event). |
| `isEvent_one` | `1` is an event (the sure event). |
| `IsEvent.compl` | `1 − P` is an event. |
| `IsEvent.orthogonal_symm` | `PQ = 0 → QP = 0` for events (star of the product). |
| `IsEvent.add` | The sum of orthogonal events is an event. |
| `IsEvent.mul_of_commute` | The product of commuting events is an event. |
| `IsEvent.absorb_of_le` | `PQ = P → QP = P` (subevent absorption). |
| `IsEvent.sub_of_le` | `Q − P` is an event when `P` is a subevent of `Q`. |
| `IsEvent.posSemidef` | Every event is positive semidefinite (`P = Pᴴ P`). |
| `IsEvent.one_sub_two_smul_involution` | `1 − 2P` is a selfadjoint involution (dichotomic observable). |
| `bornWeight_add` | Weight is additive in the event argument. |
| `bornWeight_sub` | Weight is subtractive in the event argument. |
| `bornWeight_sum` | Weight commutes with finite sums of events. |
| `bornWeight_smul` | `bornWeight (c•ρ) P = c · bornWeight ρ P`. |
| `trace_sandwich` | `Tr(PρP) = Tr(ρP)` for idempotent `P`. |
| `star_bornWeight` | Reality: the weight is fixed by conjugation (Hermitian `ρ`, `P`). |
| `bornWeight_eq_re` | The weight equals its own real part. |
| `bornWeight_nonneg` | `0 ≤ Tr(ρP)` in the partial order of `ℂ` (PSD `ρ`, event `P`). |
| `bornWeight_re_nonneg` | Real-part form of nonnegativity. |
| `bornWeight_one` | Normalisation: `Tr(ρ·1) = 1` for a state. |
| `bornWeight_add_of_orthogonal` | Additivity packaged with the orthogonality hypothesis. |
| `bornWeight_le_one` | `Tr(ρP) ≤ 1` via the complement event. |
| `bornWeight_re_le_one` | Real-part form of the upper bound. |
| `bornWeight_mono` | `PQ = P → Tr(ρP) ≤ Tr(ρQ)` (monotone under subevents). |
| `bornWeight_ne_zero_iff_re_pos` | Nonzero weight ↔ strictly positive real part. |

`Lueders.lean` : Lüders conditioning (12 lemmas):

| Lean name | One-line statement |
|---|---|
| `luedersUpdate_posSemidef` | The update is PSD (zero matrix in the degenerate case). |
| `trace_luedersUpdate` | The update has trace one when the weight is nonzero. |
| `luedersUpdate_isState` | The update of a state by a nonzero-weight event is a state. |
| `bornWeight_luedersUpdate_self` | Repeatability: the conditioned state gives `P` weight `1`. |
| `luedersUpdate_idem` | Conditioning twice on `P` equals conditioning once. |
| `luedersUpdate_luedersUpdate_of_commute` | Sequential conditioning on commuting events = conditioning on `PQ`. |
| `luedersUpdate_comm` | Order exchange for commuting events. |
| `luedersUpdate_of_commute` | Classical restriction: `ρP = Pρ` ⇒ update is `(Tr ρP)⁻¹ • (ρP)`. |
| `luedersUpdate_mem_certainStates` | One step of conditioning lands in the certainty set of `P`. |
| `mul_eq_self_of_bornWeight_one` | A state certain of `P` is supported on `P` (`σP = σ = Pσ`). |
| `luedersUpdate_eq_self_of_mem_certainStates` | States certain of `P` are fixed points of conditioning. |
| `luedersUpdate_eq_self_iff` | Fixed points of conditioning = states certain of `P` (given nonzero weight). |

`CenterExpectation.lean` : conditional expectation onto a commutative
center (16 lemmas). This is the quantum-native counterpart of the
classical contractive conditional-resampling projector package in
`ObservableNormalForms/ObservableNormalForms/ConditionalResampling.lean`
(fixes-the-fiber-constants / idempotent / selfadjoint-for-the-weighted-inner-product /
Pythagoras / squared-`L²` contraction); the statements were re-proved from
scratch in the matrix event-algebra setting rather than wrapped, so the
neutral bundle stays self-contained with no imports from that package:

| Lean name | One-line statement |
|---|---|
| `CenterPartition.proj_mul_proj` | `Pᵢ Pⱼ = δᵢⱼ Pᵢ` for a partition of unity. |
| `CenterPartition.proj_commute` | Partition members commute (the center is commutative). |
| `centerExpectation_mul_proj` | Right absorption: `𝔼(X) Pᵢ = Pᵢ X Pᵢ`. |
| `proj_mul_centerExpectation` | Left absorption: `Pᵢ 𝔼(X) = Pᵢ X Pᵢ`. |
| `proj_commute_centerExpectation` | `𝔼(X)` lies in the commutant of the partition. |
| `centerExpectation_fixes` | `𝔼` fixes the commutant pointwise. |
| `centerExpectation_idem` | `𝔼 ∘ 𝔼 = 𝔼` (a projector). |
| `conjTranspose_centerExpectation` | `𝔼` commutes with conjugate transposition. |
| `centerExpectation_isState` | `𝔼` maps states to states. |
| `trace_conjTranspose_centerExpectation_mul` | `𝔼` is selfadjoint for the trace inner product. |
| `trace_centerExpectation_pythagoras` | Pythagoras: `‖X‖² = ‖𝔼X‖² + ‖X − 𝔼X‖²` (trace norms). |
| `trace_centerExpectation_contraction` | Squared-`L²` contractivity of `𝔼`. |
| `trace_centerExpectation_mul_central` | Trace compatibility against every central `C`. |
| `centerExpectation_unique` | `𝔼X` is the unique commutant-valued, trace-compatible element. |
| `bornWeight_centerExpectation` | Born weights of central events are `𝔼`-invariant. |
| `centerExpectation_luedersUpdate` | For central `P`: `𝔼 ∘ L_P = L_P ∘ 𝔼` (classical conditioning). |

`StateFromTrace.lean` : the expectation functional (5 lemmas):

| Lean name | One-line statement |
|---|---|
| `bornWeight_eq_expectation` | The Born weight is the expectation functional at an event. |
| `expectation_add` | Additivity in the observable. |
| `expectation_smul` | Homogeneity in the observable. |
| `expectation_one` | Normalisation: `Tr(ρ·1) = 1` for a state. |
| `expectation_nonneg` | Positivity on all PSD observables (spectral-theorem proof). |

`Tsirelson.lean` : the Tsirelson bound, norm form (6 lemmas; abstract
C*-route shipped, matrix corollary included):

| Lean name | One-line statement |
|---|---|
| `chsh_mul_self` | Ring identity `S² = 4·1 − (a₀a₁−a₁a₀)(b₀b₁−b₁b₀)` (no star/order/norm). |
| `norm_eq_one_of_selfAdjoint_involution` | Selfadjoint involutions have norm one in a unital C*-ring. |
| `norm_commutator_le_two` | Commutators of norm-one elements have norm ≤ 2. |
| `tsirelson_bound` | `‖a₀b₀ + a₀b₁ + a₁b₀ − a₁b₁‖ ≤ 2√2` in any nontrivial unital C*-ring. |
| `tsirelson_bound_of_isCHSHTuple` | The same, consuming Mathlib's `IsCHSHTuple` bundle. |
| `matrix_tsirelson_bound` | Instantiation at `Matrix (Fin n) (Fin n) ℂ` with the L2 operator norm. |

No attainment/tightness claim is made anywhere in `Tsirelson.lean`; the
bound is shipped as an inequality only. Mathlib's
`Mathlib/Algebra/Star/CHSH.lean` provides the complementary *order-form*
`tsirelson_inequality`; the norm form here reuses its `IsCHSHTuple`.

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
