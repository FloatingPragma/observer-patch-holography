import Mathlib
import ObserverPatchHolography.AbstractRewriting

/-!
# OPH Primitives — concrete carrier model (partial discharge)

These are the primitives Proposition 4.2 depends on. Where the companion
paper *Reality as a Consensus Protocol* (`OPHConsensus`) pins down concrete
structural content, we now give it: the patch-net carrier, the global state
type `Records`, the declared-overlap observation map, gauge equivalence as
the kernel of that map, and the weighted mismatch potential `Φ`.

The genuinely paper-incomplete asynchronous-schedule / transactional
machinery (`localRepair`, `Repair`, and the congruence
`repair_respects_gauge` that depends on a fully constructed `Repair`)
remains `sorry`-bearing **on purpose**: `lake build` warns on each, and CI
checks that the count stays fixed until they are discharged.

## What is concrete now (from the paper)

* `OPHCarrier` — *Reality* Def 1.1 (finite patch graph `G=(V,E)`; per-patch
  finite state spaces `S_i`; per-edge interface alphabet `I_e` and
  projections `π_{i,e}, π_{j,e}`) + Def 2 (edge weights `w_e > 0` and a
  per-edge distance `d_e` with `d_e(a,b)=0 ↔ a=b`).
* `Records C := (i : C.Patch) → C.State i` — *Reality* Def 1.1 global state
  space `Σ := ∏_{i∈V} S_i`.
* `Obs C` / `obsMap C` — *Paradise* line 311 declared observable overlap
  data: the per-edge exposed projection pair `e ↦ (π_{i,e}(x_i), π_{j,e}(x_j))`.
* `Φ C` — *Reality* Def 2 / *Paradise* line 300:
  `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))`.
* `gaugeEquiv C` — *Paradise* line 311: the kernel `Setoid.ker (obsMap C)`
  (same declared observable overlap data).
* `gaugeEquiv_equivalence` — `∼_gauge` is an equivalence relation (the kernel
  of any map is an equivalence); discharged by the from-first-principles term
  `⟨fun _ => rfl, Eq.symm, Eq.trans⟩` since `gaugeEquiv` unfolds to an `Eq`.
* `consistent_iff_edgeConsistent` — *Reality* Prop 1: `C = Φ⁻¹(0)`, the
  faithfulness witness keeping the `Φ` model from vacuously falsifying
  `Completeness`.
* `Site C` — *Reality* repair-site index (a local move fires at a patch).
* `demoCarrier` / `obsMap_demoCarrier_nonconstant` — an explicit two-patch
  carrier and a proof that its `obsMap` separates two records. This makes the
  non-vacuity of `gaugeEquiv`/`consistent_iff_edgeConsistent` an in-file fact
  (gaugeEquiv is strictly finer than the total relation), not merely an
  argued universal claim. Adds no `sorry`.

## What stays `sorry` (paper-incomplete async machinery)

* `localRepair`, `Repair` — "built from local recovery moves" (line 297),
  composed under asynchronous schedules in `OPHConsensus`; not pinned to a
  constructive operator with a discharged Lyapunov+confluence proof.
* `repair_respects_gauge` — Prop 4.2 sentence 2 congruence; honestly
  unprovable while `Repair` itself is undefined (faking `Repair := id`
  would make `LyapunovDescent` vacuous and the congruence trivial for the
  wrong reason).
-/

namespace OPH

open Relation  -- `ReflTransGen` lives in the `Relation` namespace (cf. AbstractRewriting.lean)

/-- A finite OPH carrier: the patch graph `G=(V,E)` with per-patch state
    spaces, per-edge interface alphabets and projections, edge weights, and
    per-edge distances. Faithful encoding of *Reality* Def 1.1 + Def 2.

    Paper edges are unordered `{i,j}`; here each edge carries a fixed
    representative orientation `(src e, tgt e)`. This is sound: edge
    consistency `π_{i,e}(s_i) = π_{j,e}(s_j)` is symmetric and `Φ` is
    orientation-independent, so no further quotient on edges is needed. -/
structure OPHCarrier where
  /-- Observer patches `V` (vertices of the finite graph `G`). -/
  Patch : Type
  /-- `V` is finite. -/
  [patchFintype : Fintype Patch]
  /-- Patches have decidable equality (needed for, e.g., discrete metrics). -/
  [patchDecEq : DecidableEq Patch]
  /-- Per-patch local state space `S_i`. A genuine `Patch`-indexed family,
      NOT one shared type — faithful to projections out of *different*
      state spaces. -/
  State : Patch → Type
  /-- Interface edges `E` of the finite graph. -/
  Edge : Type
  /-- `E` is finite (so `Φ` is a finite sum). -/
  [edgeFintype : Fintype Edge]
  /-- Chosen source endpoint `i` of edge `e = {i,j}`. -/
  src : Edge → Patch
  /-- Chosen target endpoint `j` of edge `e = {i,j}`. -/
  tgt : Edge → Patch
  /-- Interface alphabet `I_e`. -/
  Iface : Edge → Type
  /-- Interface projection `π_{i,e} : S_i → I_e`. -/
  projSrc : (e : Edge) → State (src e) → Iface e
  /-- Interface projection `π_{j,e} : S_j → I_e`. -/
  projTgt : (e : Edge) → State (tgt e) → Iface e
  /-- Edge weight `w_e`. -/
  weight : Edge → NNReal
  /-- Per-edge distance `d_e` on the interface alphabet. -/
  dist : (e : Edge) → Iface e → Iface e → NNReal
  /-- *Reality* Def 2: weights are strictly positive. -/
  weight_pos : ∀ e : Edge, 0 < weight e
  /-- *Reality* Def 2: `d_e` separates points (`d_e(a,b)=0 ↔ a=b`). -/
  dist_eq_zero : ∀ (e : Edge) (a b : Iface e), dist e a b = 0 ↔ a = b

attribute [instance] OPHCarrier.patchFintype OPHCarrier.patchDecEq OPHCarrier.edgeFintype

variable (C : OPHCarrier)

/-- *Reality* Def 1.1: the global state space `Σ := ∏_{i∈V} S_i` — an
    assignment of a local state to every patch. (`Paradise` macro `\Records`.) -/
def Records : Type := (i : C.Patch) → C.State i

/-- *Paradise* line 311: the type of declared observable overlap data — the
    per-edge exposed projection-pair family. (`Paradise` macro `\Obs`.) -/
def Obs : Type := (e : C.Edge) → C.Iface e × C.Iface e

/-- The declared observable overlap data of a record: on every edge, the
    pair of interface projections it exposes,
    `e ↦ (π_{i,e}(x_i), π_{j,e}(x_j))` (*Paradise* line 311). This is a
    real, generally-non-constant map; `gaugeEquiv` is its kernel. -/
def obsMap (x : Records C) : Obs C :=
  fun e => (C.projSrc e (x (C.src e)), C.projTgt e (x (C.tgt e)))

/-- *Reality* repair-site index: a local accepted repair step fires at a
    patch. A faithful, non-vacuous index type (it does NOT trivialise
    `localRepair`, which remains a genuine `sorry`). -/
def Site : Type := C.Patch

/-- One transactional/local recovery move at a repair site.
    **Paper-incomplete async machinery — honest `sorry`.** -/
noncomputable def localRepair : Site C → Records C → Records C := sorry

/-- The composite confluent repair operator reaching a normal form.
    **Paper-incomplete async machinery — honest `sorry`.** -/
noncomputable def Repair : Records C → Records C := sorry

/-- One accepted asynchronous repair step: some site's local move changes
    the record. This is the relation the generic abstract-rewriting
    skeleton must eventually instantiate. -/
def acceptedStep (x y : Records C) : Prop :=
  ∃ i : Site C, y = localRepair C i x ∧ localRepair C i x ≠ x

/-- *Reality* Def 2 / *Paradise* line 300: the weighted edge-mismatch
    potential `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))`. A finite
    `Finset.sum` over the (finite) edge set, valued in `ℝ≥0`. -/
noncomputable def Φ (x : Records C) : NNReal :=
  ∑ e : C.Edge, C.weight e * C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e)))

/-- A normal form: no accepted repair step applies. -/
def NormalForm (x : Records C) : Prop :=
  ∀ y : Records C, ¬ acceptedStep C x y

/-- Consistency: zero mismatch potential. By `consistent_iff_edgeConsistent`
    this coincides with the paper's `C = Φ⁻¹(0)` (edge-by-edge agreement). -/
def Consistent (x : Records C) : Prop :=
  Φ C x = 0

/-- Edge-consistency (*Reality* Def 1.1): every edge's two projections agree.
    `C := {s : ∀ e, π_{src e}(s) = π_{tgt e}(s)}`. -/
def EdgeConsistent (x : Records C) : Prop :=
  ∀ e : C.Edge, C.projSrc e (x (C.src e)) = C.projTgt e (x (C.tgt e))

/-- *Reality* Prop 1: the model satisfies `C = Φ⁻¹(0)` — `Φ x = 0` holds iff
    `x` is edge-consistent. This is the faithfulness witness for the `Φ`
    model (it is what stops `Φ` from vacuously falsifying `Completeness`);
    it uses both carrier hypotheses `weight_pos` and `dist_eq_zero`. -/
theorem consistent_iff_edgeConsistent (x : Records C) :
    Consistent C x ↔ EdgeConsistent C x := by
  unfold Consistent EdgeConsistent Φ
  -- Use the nonneg-codomain form `sum_eq_zero_iff_of_nonneg`: it needs only
  -- `AddCommMonoid + PartialOrder + AddLeftMono` (all held by `ℝ≥0`) and takes
  -- the pointwise `0 ≤ ·` proof explicitly, so it avoids the `Subsingleton
  -- (AddUnits ·)` instance search that the bare `Finset.sum_eq_zero_iff`
  -- relies on. (`zero_le _` is the canonical `0 ≤ x` on `ℝ≥0`.)
  rw [Finset.sum_eq_zero_iff_of_nonneg (fun i _ => zero_le _)]
  constructor
  · intro h e
    have he := h e (Finset.mem_univ e)
    rcases mul_eq_zero.mp he with hw | hd
    · exact absurd hw (C.weight_pos e).ne'
    · exact (C.dist_eq_zero e _ _).mp hd
  · intro h e _
    have hd : C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e))) = 0 :=
      (C.dist_eq_zero e _ _).mpr (h e)
    rw [hd, mul_zero]

/-- The Lyapunov-descent obligation: every accepted step strictly lowers `Φ`. -/
def LyapunovDescent : Prop :=
  ∀ x y : Records C, acceptedStep C x y → Φ C y < Φ C x

/-- Termination of the accepted-step relation. -/
def Termination : Prop :=
  WellFounded (fun y x : Records C => acceptedStep C x y)

/-- *Paradise* line 311: two records are gauge-equivalent iff they expose the
    same declared observable overlap data. Idiomatically, this is the
    **kernel setoid** `Setoid.ker (obsMap C)`: `gaugeEquiv C x y` unfolds to
    `obsMap C x = obsMap C y`. It is non-vacuous — strictly finer than the
    total relation whenever `obsMap` is non-constant. -/
def gaugeEquiv (x y : Records C) : Prop :=
  (Setoid.ker (obsMap C)).r x y

/-- `∼_gauge` is an equivalence relation. True for the structural reason that
    `gaugeEquiv` is the kernel of `obsMap`: `gaugeEquiv C x y` unfolds (through
    `Setoid.ker` and `Function.onFun`) to the genuine equality
    `obsMap C x = obsMap C y`, whose reflexivity/symmetry/transitivity are
    `rfl`/`Eq.symm`/`Eq.trans`. We discharge it with this from-first-principles
    term rather than `(Setoid.ker (obsMap C)).iseqv` to avoid relying on the
    `.r`-vs-η defeq between `Equivalence (gaugeEquiv C)` and
    `Equivalence ⇑(Setoid.ker (obsMap C))`. -/
theorem gaugeEquiv_equivalence : Equivalence (gaugeEquiv C) :=
  ⟨fun _ => rfl, Eq.symm, Eq.trans⟩

/-- `∼_gauge` is a `Repair`-congruence. Required by Prop 4.2 sentence 2
    (independence on the physical quotient).

    **Honest `sorry`.** This cannot be soundly proved while `Repair` itself
    is a `sorry`: the only `Repair` instances that close it for free are
    degenerate (`Repair := id` / a constant), which would simultaneously
    make `Termination`/`Confluence`/`Completeness`/`LyapunovDescent` vacuous
    or false. The honest content of Prop 4.2 sentence 2 is precisely that
    the real (async) `Repair` factors through `obsMap`; that is discharged
    only once `Repair` is the genuine consensus operator. -/
theorem repair_respects_gauge :
    ∀ x y : Records C, gaugeEquiv C x y → gaugeEquiv C (Repair C x) (Repair C y) :=
  sorry

/-- OPH confluence condition for accepted asynchronous repair steps
    (Prop 4.2 hypothesis; defined per OPHConsensus). -/
def Confluence : Prop :=
  ∀ x y z : Records C, ReflTransGen (acceptedStep C) x y → ReflTransGen (acceptedStep C) x z →
    ∃ w : Records C, ReflTransGen (acceptedStep C) y w ∧ ReflTransGen (acceptedStep C) z w

/-- OPH repair completeness: normal forms are exactly consistent states.
    Termination is a separate Lyapunov/finite-state obligation. -/
def Completeness : Prop :=
  ∀ x : Records C, NormalForm C x ↔ Consistent C x

/-! ## Non-vacuity witness

A concrete two-patch / one-edge carrier exhibiting that `obsMap` is genuinely
non-constant, so `gaugeEquiv` is strictly finer than the total relation and
`consistent_iff_edgeConsistent` is a statement about a model that actually
exists. This is the explicit anti-degeneracy witness (no `sorry`,
`weight_pos`/`dist_eq_zero` discharged for a real instance). -/

/-- A concrete carrier: two patches `Bool`, one edge `()`, interface `Bool`,
    identity projections, unit weight, and the discrete `{0,1}` distance. -/
def demoCarrier : OPHCarrier where
  Patch := Bool
  State := fun _ => Bool
  Edge := Unit
  src := fun _ => false
  tgt := fun _ => true
  Iface := fun _ => Bool
  projSrc := fun _ s => s
  projTgt := fun _ s => s
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases h : a = b
    · rw [if_pos h]; exact ⟨fun _ => h, fun _ => rfl⟩
    · rw [if_neg h]; exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 h⟩

/-- The observation map of `demoCarrier` is non-constant: the all-`false`
    record and the identity record expose different declared overlap data on
    the single edge (they disagree on the target projection). Hence
    `gaugeEquiv demoCarrier` is strictly finer than the total relation. -/
theorem obsMap_demoCarrier_nonconstant :
    obsMap demoCarrier (fun _ => false) ≠ obsMap demoCarrier (fun b => b) := by
  -- Reduce to the single edge `()` and read off the target component:
  -- it is `false` on the all-`false` record and `true` on the identity record.
  -- We extract a *concrete* `Bool` equality (`false = true`) before deciding,
  -- rather than asking for `Decidable` of the function-typed `obsMap` equality.
  intro h
  have hpt : ((false : Bool), (false : Bool)) = ((false : Bool), (true : Bool)) :=
    congrFun h ()
  exact absurd (congrArg Prod.snd hpt) (by decide)

end OPH

/-! ## Global termination & completeness from LOCAL repair laws

This section proves the mathematical content of two of OPH's open *dynamical*
obligations — **Termination** and **Completeness** (cf. the `Termination`/
`Completeness` `def`s above) — **conditionally, for any local repair move
satisfying the local laws `H1`/`H2`/`H3` below**, derived from those explicit,
faithful, single-site properties. It does **not** close the file's own
`Termination`/`Completeness` `def`s (those are stated over the still-placeholder
`sorry`-defined `localRepair`, so cannot be discharged until `localRepair` is
defined); it establishes the theorems for the abstract move `lr` instead. The laws are
satisfiable by a genuine repair (e.g. a two-`Bool`-patch carrier with one
edge, each patch copying its neighbour to snap the edge consistent), so the
result is conditional, not vacuous — and that satisfiability is itself
machine-checked below as `demoCarrier_terminates` (a concrete `(carrier, repair)`
instance discharging `H1`/`H2`/`H3` with a real, non-empty repair step).

It is deliberately **self-contained and axiom-clean**: it does *not* reference
the `sorry`-defined `localRepair`/`Repair`. The repair move and its laws enter
as `section variable`s (`lr`, `H1`, `H2`, `H3`), so each theorem here closes
with `#print axioms` reporting only `[propext, Classical.choice, Quot.sound]`
(no `sorryAx`, no new `axiom`). Because the file's own `Termination`/
`Completeness` are stated over the `sorry`-defined `acceptedStep`, they cannot
be discharged without touching that `sorry`; the honest, axiom-clean statements
are therefore phrased over the hypothesis-bearing move `lr` (`acceptedStepLR`,
`NormalFormLR`) and are mathematically the same theorems for the real operator
once it satisfies `H1`/`H2`/`H3`.

### Hypotheses are LOCAL; conclusions are GLOBAL (no assume-the-conclusion)

The hypotheses are genuine **single-site** statements:
* `H1` (`lr` changes only site `i`): firing at `i` touches patch `i` only.
* `H2` (`lr` fires iff a local edge is broken): the move at `i` changes `x`
  *iff* some edge incident to `i` is currently inconsistent — a purely local
  trigger.
* `H3` (local satisfiability / frustration-freeness): when the move at `i`
  fires it makes *all* of `i`'s own incident edges consistent. This honestly
  restricts to carriers where a single patch *can* satisfy all its overlaps at
  once (frustration-free); it is a local property, **not** the global claim.

The conclusions are **global** dynamical facts about all of `Records C`:
* `termination`: the asynchronous accepted-step relation is `WellFounded`.
* `completeness`: a record is a global normal form *iff* it is globally
  `Consistent` (`Φ = 0`).

None of the forbidden shortcuts is assumed: we never assume `mismatchCount`
decreases, nor `WellFounded`, nor `Termination`, nor `NormalForm ↔ Consistent`.
Those are *proved* from the three local laws (plus the already-discharged
`consistent_iff_edgeConsistent`).

### The Lyapunov / Inter-Basin termination pattern

The proof is the well-founded-measure pattern: every accepted repair strictly
lowers a **structural `ℕ` surrogate** `mismatchCount` (the number of broken
edges), exactly as every SKI reduction strictly lowers `basin_size` in the
Inter-Basin termination theorem. A `ℕ` surrogate is *needed* because the
carrier potential `Φ : ℝ≥0` is **not** `<`-well-founded; `mismatchCount` is the
well-founded shadow of `Φ` that makes asynchronous descent terminate.

### What remains open (honest scoping; no `sorry`)

`Confluence`/`LocallyConfluent` is **not** provided: asynchronous repairs at
different sites need not commute (`lr i (lr j x)` and `lr j (lr i x)` can
differ), so a frustration-free carrier may still reach distinct normal forms
under distinct schedules — schedule independence / unique normal forms is out
of scope for these hypotheses. There is no `sorry`, `admit`, or new `axiom`
anywhere in this section. -/

namespace OPH

open Relation  -- `ReflTransGen` (used by the confluence theorems below)

section LocalRepairDynamics

variable {C : OPHCarrier}

/-- An edge is consistent at `x` when its two interface projections agree.
    A `Prop` (no `DecidableEq (Iface e)` needed). By `dist_eq_zero` this is
    equivalent to the edge's per-edge distance vanishing
    (`edgeConsistentAt_iff_dist`). Definitionally, `EdgeConsistent C x` is
    `∀ e, edgeConsistentAt e x`. -/
def edgeConsistentAt (e : C.Edge) (x : Records C) : Prop :=
  C.projSrc e (x (C.src e)) = C.projTgt e (x (C.tgt e))

/-- Bridge to the decidable surrogate used by `mismatchCount`: an edge is
    consistent iff its per-edge distance is `0` (uses `dist_eq_zero`). -/
theorem edgeConsistentAt_iff_dist (e : C.Edge) (x : Records C) :
    edgeConsistentAt e x ↔
      C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e))) = 0 :=
  (C.dist_eq_zero e _ _).symm

/-- The set of broken edges of `x`: those whose per-edge distance is nonzero.
    This is decidable *without* `DecidableEq (Iface e)`, because `ℝ≥0` has
    `DecidableEq` (from its `LinearOrder`), so `(· ≠ 0)` is a `DecidablePred`. -/
noncomputable def brokenSet (x : Records C) : Finset C.Edge :=
  Finset.univ.filter
    (fun e => C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e))) ≠ 0)

/-- The well-founded `ℕ` surrogate for `Φ`: the number of broken edges.
    (`Φ : ℝ≥0` is not `<`-well-founded; this `ℕ` shadow is.) -/
noncomputable def mismatchCount (x : Records C) : Nat := (brokenSet x).card

theorem mem_brokenSet {x : Records C} {e : C.Edge} :
    e ∈ brokenSet x ↔
      C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e))) ≠ 0 := by
  unfold brokenSet
  rw [Finset.mem_filter]
  exact ⟨fun h => h.2, fun h => ⟨Finset.mem_univ e, h⟩⟩

/-- An edge is broken at `x` exactly when it is *not* consistent there.
    (`mem_brokenSet` composed with the `dist`-bridge `edgeConsistentAt_iff_dist`,
    using `Ne` `=` `¬ (· = ·)` definitionally.) -/
theorem mem_brokenSet_iff_not_consistent {x : Records C} {e : C.Edge} :
    e ∈ brokenSet x ↔ ¬ edgeConsistentAt e x :=
  mem_brokenSet.trans (not_congr (edgeConsistentAt_iff_dist e x)).symm

-- The abstract local-repair move under study (a section variable):
-- `lr i x` applies the recovery move at site `i` to record `x`.
variable (lr : C.Patch → Records C → Records C)

/-- One accepted asynchronous repair step *for the abstract move `lr`*: some
    site's local move changes the record. Self-contained analogue of the file's
    `acceptedStep`, but over the hypothesis-bearing `lr`, so this section never
    touches the `sorry`-defined `localRepair`. -/
def acceptedStepLR (x y : Records C) : Prop :=
  ∃ i : C.Patch, y = lr i x ∧ lr i x ≠ x

/-- A normal form for `lr`: no accepted `lr`-step applies. -/
def NormalFormLR (x : Records C) : Prop :=
  ∀ y : Records C, ¬ acceptedStepLR lr x y

-- H1 (local: changes only site i): firing the move at site i alters the state
-- of patch i only; every other patch keeps its state.
-- H2 (local trigger: fires iff a local edge is broken): the move at i changes x
-- iff some edge incident to i is currently inconsistent.
-- H3 (local satisfiability / frustration-freeness): when the move at i fires it
-- makes all of i's incident edges consistent (restricts to carriers where a
-- single patch can satisfy all its overlaps at once).
variable
  (H1 : ∀ (i : C.Patch) (x : Records C) (j : C.Patch), j ≠ i → (lr i x) j = x j)
  (H2 : ∀ (i : C.Patch) (x : Records C),
      lr i x ≠ x ↔
        ∃ e : C.Edge, (C.src e = i ∨ C.tgt e = i) ∧ ¬ edgeConsistentAt e x)
  (H3 : ∀ (i : C.Patch) (x : Records C),
      lr i x ≠ x →
        ∀ e : C.Edge, (C.src e = i ∨ C.tgt e = i) → edgeConsistentAt e (lr i x))

-- Thread `lr`, `H1`, `H2`, `H3` uniformly through every theorem below, in this
-- fixed order, so cross-references are unambiguous. (Some lemmas don't use all
-- four; the extra hypotheses are harmless and keep call sites uniform.)
include lr H1 H2 H3

/-- A non-incident edge keeps both its endpoint states, hence its broken-ness,
    when site `i` fires (immediate from `H1`). -/
theorem brokenSet_eq_of_not_incident
    {i : C.Patch} {x : Records C} {e : C.Edge}
    (hs : C.src e ≠ i) (ht : C.tgt e ≠ i) :
    (e ∈ brokenSet (lr i x) ↔ e ∈ brokenSet x) := by
  have hsrc : (lr i x) (C.src e) = x (C.src e) := H1 i x (C.src e) hs
  have htgt : (lr i x) (C.tgt e) = x (C.tgt e) := H1 i x (C.tgt e) ht
  rw [mem_brokenSet, mem_brokenSet, hsrc, htgt]

/-- **Key lemma — Lyapunov descent on the `ℕ` surrogate.** Every accepted step
    strictly lowers `mismatchCount`: the broken-edge set strictly shrinks. This
    is the Inter-Basin `basin_size`-strictly-decreases analogue. -/
theorem mismatchCount_lt {x y : Records C}
    (h : acceptedStepLR lr x y) : mismatchCount y < mismatchCount x := by
  obtain ⟨i, rfl, hfire⟩ := h
  -- It suffices to show `brokenSet (lr i x) ⊂ brokenSet x`; then `card_lt_card`.
  -- (1) Subset: an edge broken in `lr i x` cannot be incident to `i` (those are
  -- made consistent by `H3`), and on non-incident edges broken-ness transfers
  -- back to `x` (`brokenSet_eq_of_not_incident`).
  have hsub : brokenSet (lr i x) ⊆ brokenSet x := by
    intro e he
    by_cases hinc : C.src e = i ∨ C.tgt e = i
    · have hcon : edgeConsistentAt e (lr i x) := H3 i x hfire e hinc
      exact absurd hcon (mem_brokenSet_iff_not_consistent.1 he)
    · have hs : C.src e ≠ i := fun h => hinc (Or.inl h)
      have ht : C.tgt e ≠ i := fun h => hinc (Or.inr h)
      exact (brokenSet_eq_of_not_incident lr H1 H2 H3 hs ht).1 he
  -- (2) Strictness: `H2` exhibits an incident broken edge of `x`; it lies in
  -- `brokenSet x` but not in `brokenSet (lr i x)` (incident → consistent there).
  obtain ⟨e₀, hinc₀, hbroken₀⟩ := (H2 i x).1 hfire
  have hmem_x : e₀ ∈ brokenSet x := mem_brokenSet_iff_not_consistent.2 hbroken₀
  have hcon₀ : edgeConsistentAt e₀ (lr i x) := H3 i x hfire e₀ hinc₀
  have hnot_mem : e₀ ∉ brokenSet (lr i x) :=
    fun hm => mem_brokenSet_iff_not_consistent.1 hm hcon₀
  have hssub : brokenSet (lr i x) ⊂ brokenSet x :=
    (Finset.ssubset_iff_of_subset hsub).2 ⟨e₀, hmem_x, hnot_mem⟩
  exact Finset.card_lt_card hssub

/-- **THEOREM — Termination (global).** The accepted asynchronous-repair
    relation is well-founded. Derived purely from the local laws via the `ℕ`
    measure `mismatchCount`, as the inverse image of `<` on `ℕ` and a
    sub-relation thereof. -/
theorem termination :
    WellFounded (fun y x : Records C => acceptedStepLR lr x y) :=
  -- Same idiom as `Finset.lt_wf`: the step relation is a sub-relation of the
  -- inverse image of `<` on `ℕ` under `mismatchCount`, which is well-founded.
  have H : Subrelation (fun y x : Records C => acceptedStepLR lr x y)
      (InvImage (· < ·) mismatchCount) :=
    fun {_ _} hxy => mismatchCount_lt lr H1 H2 H3 hxy
  Subrelation.wf H <| InvImage.wf _ wellFounded_lt

/-- Local characterisation behind completeness: site `i` is quiescent
    (`lr i x = x`) iff every edge incident to `i` is consistent
    (`H2`, contrapositive). -/
theorem lr_fixed_iff_incident_consistent (i : C.Patch) (x : Records C) :
    lr i x = x ↔ ∀ e : C.Edge, (C.src e = i ∨ C.tgt e = i) → edgeConsistentAt e x := by
  constructor
  · intro hfix e hinc
    by_contra hcon
    exact (H2 i x).mpr ⟨e, hinc, hcon⟩ hfix
  · intro hall
    by_contra hfire
    obtain ⟨e, hinc, hcon⟩ := (H2 i x).mp hfire
    exact hcon (hall e hinc)

/-- A record is a normal form iff *no* site fires (`lr i x = x` for all `i`).
    Unfolds `acceptedStepLR`/`NormalFormLR`. -/
theorem normalForm_iff_all_quiescent (x : Records C) :
    NormalFormLR lr x ↔ ∀ i : C.Patch, lr i x = x := by
  constructor
  · intro hnf i
    by_contra hfire
    exact hnf (lr i x) ⟨i, rfl, hfire⟩
  · intro hquiet y hstep
    obtain ⟨i, _, hfire⟩ := hstep
    exact hfire (hquiet i)

/-- **THEOREM — Completeness (global).** A record is a normal form of the
    accepted-step relation iff it is globally `Consistent` (`Φ = 0`). The
    bridge: no site fires ↔ every incident edge of every site is consistent ↔
    every edge is consistent (each edge is incident to its own `src`) ↔
    `EdgeConsistent` ↔ (`consistent_iff_edgeConsistent`) `Consistent`. -/
theorem completeness (x : Records C) :
    NormalFormLR lr x ↔ Consistent C x := by
  rw [normalForm_iff_all_quiescent lr H1 H2 H3 x, consistent_iff_edgeConsistent C x]
  constructor
  · intro hquiet e
    exact (lr_fixed_iff_incident_consistent lr H1 H2 H3 (C.src e) x).1
      (hquiet (C.src e)) e (Or.inl rfl)
  · intro hcons i
    exact (lr_fixed_iff_incident_consistent lr H1 H2 H3 i x).2 (fun e _ => hcons e)

-- ── Boundary-fiber observer-uniqueness (issue #304) ──────────────────────────
-- The boundary / sector map `B` (#304): a coarse invariant the repair PRESERVES.
-- `HB` = repair preserves `B`; `Hfib` = within a fixed boundary fiber, consistent
-- states are a gauge-singleton. These are #304's STATED hypotheses, not proven here.
variable {β : Type} (B : Records C → β)
  (HB : ∀ (i : C.Patch) (x : Records C), B (lr i x) = B x)
  (Hfib : ∀ x y : Records C, B x = B y → Consistent C x → Consistent C y → gaugeEquiv C x y)

include HB in
/-- The boundary map is invariant along an entire accepted-repair reduction. -/
theorem boundary_preserved_reduction {a b : Records C}
    (h : ReflTransGen (acceptedStepLR lr) a b) : B b = B a := by
  induction h with
  | refl => rfl
  | tail _ hstep ih =>
      obtain ⟨i, hc, _⟩ := hstep
      rw [hc, HB]; exact ih

include H1 H2 H3 HB Hfib in
/-- **THEOREM — Boundary-fiber observer-uniqueness (issue #304, observer-facing half).**
    Any two records with the SAME boundary value settle to the same observer-facing
    normal form (`gaugeEquiv`). It needs only `completeness` (normal form ⟹ consistent)
    + boundary preservation (`HB`) + the singleton-consistent-fiber hypothesis (`Hfib`)
    — confluence does NOT enter. Conditional on H1–H3 + HB + Hfib (exactly #304's stated
    hypotheses). A non-vacuous witness needs a boundary-PINNING (directional) repair:
    a single shared edge has two consistent states `(0,0)`/`(1,1)` and the symmetric
    copy-move reaches either, so the singleton fiber fails for it
    (cf. `demoCarrier_not_confluent`). -/
theorem boundary_fiber_observer_unique {x y nfx nfy : Records C}
    (hB : B x = B y)
    (hx : ReflTransGen (acceptedStepLR lr) x nfx) (hxn : NormalFormLR lr nfx)
    (hy : ReflTransGen (acceptedStepLR lr) y nfy) (hyn : NormalFormLR lr nfy) :
    gaugeEquiv C nfx nfy := by
  have hBx : B nfx = B x := boundary_preserved_reduction lr H1 H2 H3 B HB hx
  have hBy : B nfy = B y := boundary_preserved_reduction lr H1 H2 H3 B HB hy
  have hCx : Consistent C nfx := (completeness lr H1 H2 H3 nfx).1 hxn
  have hCy : Consistent C nfy := (completeness lr H1 H2 H3 nfy).1 hyn
  have hBB : B nfx = B nfy := by rw [hBx, hB, ← hBy]
  exact Hfib nfx nfy hBB hCx hCy

-- H4 (GLOBAL commutation): EVERY ordered pair of sites `i, j` commutes on every
-- record (the classical *diamond* condition, stated globally). This is a STRONG
-- hypothesis: it demands even adjacent, edge-sharing sites commute, which the
-- natural copy-repair does NOT satisfy (`demoCarrier_repairs_dont_commute`). So
-- H4 is a SUFFICIENT extra law for global Confluence, not a necessary one, and it
-- has no non-trivial witness in this file — the only carrier, `demoCarrier`,
-- violates it (and is in fact non-confluent, `demoCarrier_not_confluent`). The
-- honest weaker hypothesis would restrict commutation to NON-INCIDENT pairs
-- (sites sharing no edge, expressible via the incidence predicate already used in
-- H2/H3), but that alone does not close the diamond when incident repairs
-- genuinely diverge. H4 is NOT implied by H1–H3.
variable
  (H4 : ∀ (i j : C.Patch) (x : Records C), lr i (lr j x) = lr j (lr i x))

-- H4 is used only inside the proofs below (the statements quantify over the
-- abstract relation `acceptedStepLR lr`), so force its inclusion explicitly.
include H4

-- H1/H2/H3 are unused by the diamond argument (it needs only `H4`); drop them
-- from THIS lemma so its type honestly reads "commuting moves are locally
-- confluent". `confluence_of_commute` below still carries them (for `termination`).
omit H1 H2 H3 in
/-- **Local confluence from single-step commutation** (the diamond condition).
    From two accepted steps at sites `i`, `j`, the common join is
    `lr j (lr i x) = lr i (lr j x)` (by `H4`); each side reaches it in ≤ 1 step
    (zero if that site is already quiescent there). -/
theorem locallyConfluent_of_commute :
    AbstractRewriting.LocallyConfluent (acceptedStepLR lr) := by
  rintro x _ _ ⟨i, rfl, _⟩ ⟨j, rfl, _⟩
  refine ⟨lr j (lr i x), ?_, ?_⟩
  · rcases eq_or_ne (lr j (lr i x)) (lr i x) with h | h
    · rw [h]
    · exact ReflTransGen.single ⟨j, rfl, h⟩
  · rw [← H4 i j x]
    rcases eq_or_ne (lr i (lr j x)) (lr j x) with h | h
    · rw [h]
    · exact ReflTransGen.single ⟨i, rfl, h⟩

/-- **THEOREM — Confluence (Church–Rosser) under commutation.** Termination
    (H1–H3, via `termination`) together with local confluence (H4, via
    `locallyConfluent_of_commute`) yields global confluence, through Newman's
    lemma. With `termination` this further gives UNIQUE normal forms (the repo's
    `AbstractRewriting.newman_unique_nf`) — a schedule-independent "objective
    public reality" — but only the join property `Confluent` is concluded here.
    CAVEAT (read with `demoCarrier_not_confluent`): this is the SUFFICIENT
    direction, conditional on the strong, GLOBAL law `H4`, which has no
    non-trivial witness in this file — the only carrier, `demoCarrier`, provably
    violates `H4` and is in fact non-confluent. So this theorem says precisely
    "IF every pair of repairs commutes, schedules agree"; the witnessed,
    load-bearing fact is the negative one. -/
theorem confluence_of_commute :
    AbstractRewriting.Confluent (acceptedStepLR lr) :=
  AbstractRewriting.newman_lemma (acceptedStepLR lr)
    (termination lr H1 H2 H3) (locallyConfluent_of_commute lr H4)

end LocalRepairDynamics

/-! ## Non-vacuity witness: the local laws are satisfiable by a real repair

`demoCarrier` (two `Bool` patches, one edge) with the neighbour-copy repair
`demoLR` satisfies `H1`/`H2`/`H3` and has a non-empty accepted-step relation, so
`demoCarrier_terminates` is a genuine, non-vacuous instance of `termination` —
not a claim about an unsatisfiable hypothesis set. -/

/-- A genuine local repair on `demoCarrier`: patch `i` copies its neighbour `!i`,
    snapping the single edge consistent. Changes only patch `i`. -/
def demoLR : demoCarrier.Patch → Records demoCarrier → Records demoCarrier :=
  fun i x => Function.update x i (x (!i))

/-- `demoLR` fires (changes the record) exactly when patch `i` disagrees with its
    neighbour. -/
theorem demoLR_eq_self_iff (i : demoCarrier.Patch) (x : Records demoCarrier) :
    demoLR i x = x ↔ x (!i) = x i := by
  constructor
  · intro h
    have hi := congrFun h i
    simp only [demoLR, Function.update_self] at hi
    exact hi
  · intro h
    funext k
    rcases eq_or_ne k i with hk | hk
    · subst hk; simpa only [demoLR, Function.update_self] using h
    · simp only [demoLR, Function.update_of_ne hk]

theorem demoLR_H1 :
    ∀ (i : demoCarrier.Patch) (x : Records demoCarrier) (j : demoCarrier.Patch),
      j ≠ i → (demoLR i x) j = x j := by
  intro i x j hj
  simp only [demoLR, Function.update_of_ne hj]

theorem demoLR_H3 :
    ∀ (i : demoCarrier.Patch) (x : Records demoCarrier),
      demoLR i x ≠ x →
        ∀ e : demoCarrier.Edge,
          (demoCarrier.src e = i ∨ demoCarrier.tgt e = i) →
            edgeConsistentAt e (demoLR i x) := by
  intro i x _ e _
  show (demoLR i x) false = (demoLR i x) true
  cases i <;> rfl

theorem demoLR_H2 :
    ∀ (i : demoCarrier.Patch) (x : Records demoCarrier),
      demoLR i x ≠ x ↔
        ∃ e : demoCarrier.Edge,
          (demoCarrier.src e = i ∨ demoCarrier.tgt e = i) ∧ ¬ edgeConsistentAt e x := by
  intro i x
  rw [ne_eq, demoLR_eq_self_iff]
  have hiff : (x (!i) ≠ x i) ↔ (x false ≠ x true) := by
    cases i
    · exact ne_comm
    · exact Iff.rfl
  constructor
  · intro h
    refine ⟨(), ?_, hiff.mp h⟩
    cases i
    · exact Or.inl rfl
    · exact Or.inr rfl
  · rintro ⟨_, _, hnc⟩
    exact hiff.mpr hnc

/-- The accepted-step relation for `demoLR` is non-empty: the identity record has
    a broken edge (`false ≠ true`), so `demoLR false` fires. -/
theorem demoLR_has_step :
    ∃ x y : Records demoCarrier, acceptedStepLR demoLR x y := by
  refine ⟨(fun b => b), demoLR false (fun b => b), false, rfl, ?_⟩
  rw [ne_eq, demoLR_eq_self_iff]
  show (!false : Bool) ≠ false
  decide

/-- **Non-vacuity payoff.** `termination` instantiated on the real, non-trivial
    witness `(demoCarrier, demoLR)`. -/
theorem demoCarrier_terminates :
    WellFounded (fun y x : Records demoCarrier => acceptedStepLR demoLR x y) :=
  termination demoLR demoLR_H1 demoLR_H2 demoLR_H3

/-- **H4 fails for the natural repair.** On `demoCarrier` the two patches share
    one edge, so their copy-moves can fail to commute. Concretely, from the
    identity record `id = (fun b => b)`, repairing `false` then `true` gives the
    constant `true`, whereas `true` then `false` gives the constant `false` — a
    single record witnessing `lr i (lr j ·) ≠ lr j (lr i ·)`. Hence `demoLR`
    violates `H4`, so `confluence_of_commute` does not apply to it. This is not
    merely a failed sufficient condition — `demoLR` is in fact NON-CONFLUENT
    (`demoCarrier_not_confluent` below): the two firing orders reach two distinct
    normal forms, so on this carrier there is genuinely no unique objective public
    reality. -/
theorem demoCarrier_repairs_dont_commute :
    ∃ x : Records demoCarrier,
      demoLR true (demoLR false x) ≠ demoLR false (demoLR true x) := by
  refine ⟨(fun b => b), fun h => ?_⟩
  have h2 : (true : Bool) = false := congrFun h false
  exact absurd h2 (by decide)

/-- **THE WITNESSED PAYOFF — `demoLR` is genuinely NON-CONFLUENT.** From the
    identity record `id = (fun b => b)`, firing patch `false` reaches the constant
    `true` and firing patch `true` reaches the constant `false`; both are normal
    forms (no patch fires on a constant record) and they differ. So a single
    record has two distinct normal forms — `¬ Confluent (acceptedStepLR demoLR)` —
    the concrete failure of a unique "objective public reality" that `H4` (and
    hence `confluence_of_commute`) rules out by hypothesis. Unlike
    `confluence_of_commute` (conditional on the witness-less global `H4`), THIS
    result holds outright on the concrete `demoCarrier`. Together the three
    theorems give a self-contained picture *on this carrier*: the async copy-repair
    `demoLR` is non-confluent (here); commutation `H4` is a SUFFICIENT condition for
    confluence (`confluence_of_commute`, abstract); and `demoLR` fails it
    (`demoCarrier_repairs_dont_commute`). No claim is made that *every* async repair
    is non-confluent — only this one is exhibited.
    Proof: `AbstractRewriting.unique_normal_form` forces any two normal forms
    reached from one record to coincide; the two we exhibit do not. -/
theorem demoCarrier_not_confluent :
    ¬ AbstractRewriting.Confluent (acceptedStepLR demoLR) := by
  intro hc
  have hfire_f : demoLR false (fun b => b) ≠ (fun b => b) := by
    rw [ne_eq, demoLR_eq_self_iff]; show (!false : Bool) ≠ false; decide
  have hfire_t : demoLR true (fun b => b) ≠ (fun b => b) := by
    rw [ne_eq, demoLR_eq_self_iff]; show (!true : Bool) ≠ true; decide
  -- both single-step results are normal forms: no patch fires on a constant record
  have hnf_y : AbstractRewriting.IsNormalForm (acceptedStepLR demoLR)
      (demoLR false (fun b => b)) := by
    rintro w ⟨i, _, hne⟩; apply hne; rw [demoLR_eq_self_iff]; cases i <;> rfl
  have hnf_z : AbstractRewriting.IsNormalForm (acceptedStepLR demoLR)
      (demoLR true (fun b => b)) := by
    rintro w ⟨i, _, hne⟩; apply hne; rw [demoLR_eq_self_iff]; cases i <;> rfl
  have hy : AbstractRewriting.ReducesToNF (acceptedStepLR demoLR)
      (fun b => b) (demoLR false (fun b => b)) :=
    ⟨ReflTransGen.single ⟨false, rfl, hfire_f⟩, hnf_y⟩
  have hz : AbstractRewriting.ReducesToNF (acceptedStepLR demoLR)
      (fun b => b) (demoLR true (fun b => b)) :=
    ⟨ReflTransGen.single ⟨true, rfl, hfire_t⟩, hnf_z⟩
  -- if confluent, the two normal forms would be equal — but const true ≠ const false
  have heq := AbstractRewriting.unique_normal_form (acceptedStepLR demoLR) hc hy hz
  have h2 : (true : Bool) = false := congrFun heq false
  exact absurd h2 (by decide)

/-! ### The SYMMETRIC half of the #304 dichotomy — `demoCarrier` witnesses `Hfib` failing

`boundary_fiber_observer_unique` shows: IF the repair pins each boundary-fiber to a
single gauge class (`Hfib`), the observer reconstructs a unique public branch — and it
does so WITHOUT confluence. The theorems below exhibit the complementary countermodel:
the symmetric copy-repair `demoLR` makes `Hfib` FALSE, so the same inputs that
`boundary_fiber_observer_unique` consumes hold while its conclusion fails. The witness is
the SAME two normal forms `demoCarrier_not_confluent` already exhibits: `(1,1)` and `(0,0)`.

HONEST SCOPE: this is the FORWARD direction (`Hfib` ⟹ unique; symmetric ⟹ countermodel),
NOT a biconditional — observer-uniqueness is keyed to `Hfib` (a static fiber hypothesis),
which is logically independent of confluence. And `demoBoundary` is the trivial boundary,
legitimate as the COARSEST `B` (the fairest test of whether symmetric repair can pin its
fiber) but carrying no interior/boundary split. -/

/-- The trivial (constant) boundary on `demoCarrier` records — the concrete instance of the
    abstract `B : Records C → β` from `boundary_fiber_observer_unique`. On a single-edge
    carrier the only repair-invariant boundary is the trivial one (the coarsest `B`). -/
def demoBoundary : Records demoCarrier → Unit := fun _ => ()

/-- `demoLR` preserves `demoBoundary` (the `HB` premise of #304), trivially. -/
theorem demoBoundary_HB (i : demoCarrier.Patch) (x : Records demoCarrier) :
    demoBoundary (demoLR i x) = demoBoundary x := rfl

/-- Every constant record is `Consistent` (`Φ = 0`): the single edge's two identity
    projections agree on a constant record. -/
theorem demoCarrier_const_consistent (v : Bool) :
    Consistent demoCarrier (fun _ => v) := by
  rw [consistent_iff_edgeConsistent]; intro e; rfl

/-- The two constant normal forms are NOT gauge-equivalent: their `obsMap`s differ on the
    single edge's source projection. Mirrors `obsMap_demoCarrier_nonconstant`, read on
    `Prod.fst`. (`h` is `gaugeEquiv` = `obsMap _ = obsMap _` definitionally.) -/
theorem demoCarrier_consts_not_gaugeEquiv :
    ¬ gaugeEquiv demoCarrier (fun _ => true) (fun _ => false) := by
  intro h
  have h' : obsMap demoCarrier (fun _ => true) = obsMap demoCarrier (fun _ => false) := h
  have hpt : ((true : Bool), (true : Bool)) = ((false : Bool), (false : Bool)) :=
    congrFun h' ()
  exact absurd (congrArg Prod.fst hpt) (by decide)

/-- **COMPLEMENT THEOREM — `demoCarrier` WITNESSES `Hfib` FAILING (static form).**
    The literal negation, in #304's own vocabulary, of the singleton-consistent-fiber
    hypothesis `Hfib` of `boundary_fiber_observer_unique`, at `B := demoBoundary`:
    two `Consistent` records with equal boundary that are NOT `gaugeEquiv`. -/
theorem demoCarrier_Hfib_fails :
    ¬ (∀ x y : Records demoCarrier,
          demoBoundary x = demoBoundary y →
          Consistent demoCarrier x → Consistent demoCarrier y →
          gaugeEquiv demoCarrier x y) := by
  intro Hfib
  exact demoCarrier_consts_not_gaugeEquiv
    (Hfib (fun _ => true) (fun _ => false) rfl
      (demoCarrier_const_consistent true) (demoCarrier_const_consistent false))

/-- **COMPLEMENT THEOREM (reachability-explicit) — the SYMMETRIC half of the dichotomy.**
    From ONE start (`id`), `demoLR` reaches two normal forms with the SAME boundary that are
    NOT `gaugeEquiv`. The inputs match exactly what `boundary_fiber_observer_unique` consumes
    (reductions + `NormalFormLR` + equal boundary), yet the conclusion `gaugeEquiv` is FALSE —
    because the symmetric repair makes the one premise it does not supply, `Hfib`, fail. -/
theorem demoCarrier_boundary_fiber_not_unique :
    ∃ start nf₁ nf₂ : Records demoCarrier,
      Relation.ReflTransGen (acceptedStepLR demoLR) start nf₁ ∧ NormalFormLR demoLR nf₁ ∧
      Relation.ReflTransGen (acceptedStepLR demoLR) start nf₂ ∧ NormalFormLR demoLR nf₂ ∧
      demoBoundary nf₁ = demoBoundary nf₂ ∧
      Consistent demoCarrier nf₁ ∧ Consistent demoCarrier nf₂ ∧
      ¬ gaugeEquiv demoCarrier nf₁ nf₂ := by
  have hfire_f : demoLR false (fun b => b) ≠ (fun b => b) := by
    rw [ne_eq, demoLR_eq_self_iff]; show (!false : Bool) ≠ false; decide
  have hfire_t : demoLR true (fun b => b) ≠ (fun b => b) := by
    rw [ne_eq, demoLR_eq_self_iff]; show (!true : Bool) ≠ true; decide
  have hnf₁ : NormalFormLR demoLR (demoLR false (fun b => b)) := by
    rw [normalForm_iff_all_quiescent demoLR demoLR_H1 demoLR_H2 demoLR_H3]
    intro i; rw [demoLR_eq_self_iff]; cases i <;> rfl
  have hnf₂ : NormalFormLR demoLR (demoLR true (fun b => b)) := by
    rw [normalForm_iff_all_quiescent demoLR demoLR_H1 demoLR_H2 demoLR_H3]
    intro i; rw [demoLR_eq_self_iff]; cases i <;> rfl
  have heq₁ : demoLR false (fun b => b) = (fun _ => true) := by funext k; cases k <;> rfl
  have heq₂ : demoLR true (fun b => b) = (fun _ => false) := by funext k; cases k <;> rfl
  refine ⟨(fun b => b), demoLR false (fun b => b), demoLR true (fun b => b),
    Relation.ReflTransGen.single ⟨false, rfl, hfire_f⟩, hnf₁,
    Relation.ReflTransGen.single ⟨true, rfl, hfire_t⟩, hnf₂, rfl, ?_, ?_, ?_⟩
  · rw [heq₁]; exact demoCarrier_const_consistent true
  · rw [heq₂]; exact demoCarrier_const_consistent false
  · rw [heq₁, heq₂]; exact demoCarrier_consts_not_gaugeEquiv

/-! ### The two POSITIVE routes to observer-uniqueness (the corrected #304 cut)

`demoCarrier_Hfib_fails` shows the symmetric repair + trivial boundary breaks uniqueness. TWO
DIFFERENT levers restore it, and they are NOT the same theorem:
  ROUTE A — refine the BOUNDARY so the fiber is a gauge-singleton (Hfib holds). A property of B,
            repair-free.
  ROUTE B — use a SELECTING (deterministic) repair, which is CONFLUENT: a unique normal form per
            input, schedule-independent (Church-Rosser; a last-writer-wins / strong-eventual-
            consistency style resolver). A property of the repair dynamics, independent of Hfib. -/

/-- **ROUTE A — a finer boundary makes `Hfib` HOLD.** With `B := obsMap` (the finest sector map)
    the fiber is a gauge-singleton, so the singleton-fiber premise of `boundary_fiber_observer_unique`
    is satisfied — definitionally, since `gaugeEquiv` unfolds to `obsMap`-equality. The positive
    complement of `demoCarrier_Hfib_fails`: `Hfib` is bought by REFINING B (a boundary property),
    not by choosing a repair (which is why "selecting repair proves Hfib" was a category error). -/
theorem demoCarrier_Hfib_holds_finerB :
    ∀ x y : Records demoCarrier, obsMap demoCarrier x = obsMap demoCarrier y →
      Consistent demoCarrier x → Consistent demoCarrier y → gaugeEquiv demoCarrier x y :=
  fun _ _ h _ _ => h

/-- A SELECTING (directional) repair on `demoCarrier`: deterministically snap the edge to
    patch-`false`'s value. Unlike the symmetric `demoLR`, this is a single-valued operator (a
    last-writer-wins style resolver), so its induced rewriting is deterministic. -/
def demoDirT : Records demoCarrier → Records demoCarrier := fun x => fun _ => x false

/-- Descent potential for `demoDirT`: 1 if the edge is broken (the two patches differ),
    else 0. Phrased with `Bool.xor`/`toNat` to stay first-order — no `Decidable` synthesis
    through the dependent (semireducible) `Records demoCarrier` type. -/
def demoDirΦ : Records demoCarrier → ℕ := fun x => (Bool.xor (x true) (x false)).toNat

theorem demoDirΦ_desc (x : Records demoCarrier) :
    demoDirT x ≠ x → demoDirΦ (demoDirT x) < demoDirΦ x := by
  intro hne
  have hb : x true ≠ x false := fun h => hne (by funext k; cases k <;> simp [demoDirT, h])
  have key : Bool.xor (x true) (x false) = true := by
    cases hxt : x true <;> cases hxf : x false <;> simp_all
  simp only [demoDirΦ, demoDirT]
  rw [key]
  cases hxf : x false <;> decide

/-- **ROUTE B — the SELECTING repair is CONFLUENT (positive twin of `demoCarrier_not_confluent`).**
    `demoDirT` is deterministic, so its induced rewriting reaches a UNIQUE normal form per input,
    schedule-independent (Newman's lemma via the `DeterministicRepair` route). This is CONFLUENCE —
    the lever a directional/selecting repair actually buys — and a DIFFERENT theorem from the `Hfib`
    boundary route above. The symmetric `demoLR` provably fails it (`demoCarrier_not_confluent`). -/
theorem demoCarrier_dir_confluent :
    AbstractRewriting.Confluent (AbstractRewriting.stepRel demoDirT) :=
  AbstractRewriting.newman_lemma (AbstractRewriting.stepRel demoDirT)
    (AbstractRewriting.descent_terminating demoDirT demoDirΦ demoDirΦ_desc)
    (AbstractRewriting.deterministic_locally_confluent demoDirT)

/-- **THE CRUX — confluence is NOT observer-uniqueness.** `demoDirT` is confluent
    (`demoCarrier_dir_confluent`: same input → one normal form). But that is a property of
    the REPAIR, not of the boundary. Here are TWO inputs with the SAME (trivial) boundary that
    reach DIFFERENT normal forms `(1,1)` and `(0,0)` under that very confluent repair — so the
    observer, seeing only the boundary, still cannot reconstruct a unique reality. Confluence
    fixes order/races ("same input → same NF"); it does NOT give observer-facing uniqueness
    ("same boundary → same NF"). The lever for the latter is the boundary (`Hfib`), not the
    repair's directionality. -/
theorem demoCarrier_dir_not_observer_unique :
    ∃ x y nfx nfy : Records demoCarrier,
      Relation.ReflTransGen (AbstractRewriting.stepRel demoDirT) x nfx ∧
      AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) nfx ∧
      Relation.ReflTransGen (AbstractRewriting.stepRel demoDirT) y nfy ∧
      AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) nfy ∧
      demoBoundary x = demoBoundary y ∧
      ¬ gaugeEquiv demoCarrier nfx nfy := by
  have hx : demoDirT (fun b => !b) = (fun _ => true) := by funext k; rfl
  have hy : demoDirT (fun b => b) = (fun _ => false) := by funext k; rfl
  have hfpT : demoDirT (fun _ => true) = (fun _ => true) := by funext k; rfl
  have hfpF : demoDirT (fun _ => false) = (fun _ => false) := by funext k; rfl
  have hstepx : AbstractRewriting.stepRel demoDirT (fun b => !b) (fun _ => true) := by
    refine ⟨hx.symm, ?_⟩
    rw [hx]; intro h; simpa using congrFun h true
  have hstepy : AbstractRewriting.stepRel demoDirT (fun b => b) (fun _ => false) := by
    refine ⟨hy.symm, ?_⟩
    rw [hy]; intro h; simpa using congrFun h true
  have hnfT : AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) (fun _ => true) := by
    rintro z ⟨_, hne⟩; exact hne hfpT
  have hnfF : AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) (fun _ => false) := by
    rintro z ⟨_, hne⟩; exact hne hfpF
  exact ⟨fun b => !b, fun b => b, fun _ => true, fun _ => false,
    Relation.ReflTransGen.single hstepx, hnfT,
    Relation.ReflTransGen.single hstepy, hnfF,
    rfl, demoCarrier_consts_not_gaugeEquiv⟩

/-! ### Axiom audit — the reconstruction layer depends only on standard axioms.
The `#print axioms` outputs below confirm that the boundary-fiber reconstruction theorem
and all its concrete witnesses depend ONLY on the standard Lean/Mathlib axioms
(`propext`, `Classical.choice`, `Quot.sound`) and NOT on any of the file's three honest
`sorry`s (`localRepair`, `Repair`, `repair_respects_gauge`) — i.e. the "machine-checked"
claim for observer-reconstruction is sorry-free. -/
#print axioms boundary_fiber_observer_unique
#print axioms boundary_preserved_reduction
#print axioms demoCarrier_Hfib_fails
#print axioms demoCarrier_Hfib_holds_finerB
#print axioms demoCarrier_dir_confluent
#print axioms demoCarrier_dir_not_observer_unique

end OPH
