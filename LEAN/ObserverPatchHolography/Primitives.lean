import Mathlib

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

This section discharges two of OPH's open *dynamical* obligations —
**Termination** and **Completeness** (cf. the `Termination`/`Completeness`
`def`s above) — as genuine theorems, derived from explicit, faithful, **local**
repair properties of an abstract single-site recovery move.

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

end LocalRepairDynamics

end OPH
