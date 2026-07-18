import Mathlib
import ObserverPatchHolography.AbstractRewriting

/-!
# OPH Primitives — concrete carrier model (admission-free)

These are the primitives Proposition 4.2 depends on. Where the companion
paper *Reality as a Consensus Protocol* (`OPHConsensus`) pins down concrete
structural content, this file gives it: the patch-net carrier, the global state
type `Records`, the declared-overlap observation map, gauge equivalence as
the kernel of that map, and the weighted mismatch potential `Φ`.

The asynchronous-schedule / transactional machinery (`localRepair`,
`Repair`, and the congruence `repair_respects_gauge`) — formerly the file's
three declared `sorry`s — is now **constructed and discharged** (see
"The async machinery, now constructed" below). The file, and with it the
observer-reconstruction layer, is admission-free: no `sorry`, no `admit`,
no `native_decide`, no new axiom; CI checks that `lake build` emits zero
`sorry` warnings.

## Concrete content from the paper

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

## The async machinery, now constructed

* `localRepair i` — the single-site transactional recovery move of *Reality*
  line 297 ("built from local recovery moves"): patch `i` fires exactly when
  (a) one of its declared overlaps is broken (`LocalTrigger`) and (b) it can
  transactionally re-satisfy *all* of its declared overlaps at once
  (`LocallySolvable` — frustration-freeness at `(i, x)`), and it then
  replaces its own state with a repaired one chosen (via `Classical.choice`)
  from the *declared overlap data alone*. That last point is the design
  load-bearer: because the trigger, the solvability predicate, and the chosen
  replacement are all functions of `obsMap x` (not of the gauge-hidden
  interior), every lemma needed for Prop 4.2 sentence 2 follows.
* `Repair` — the asynchronous schedule composed to a normal form: repeatedly
  fire a (choice-canonical) firing site until none fires. Well-founded by the
  broken-edge count `mismatchCount` (each accepted move strictly shrinks the
  broken-edge set). `Repair_reachable` and `Repair_normalForm` certify that
  the result is a genuine normal form of the accepted-step relation reached
  by accepted asynchronous steps.
* `repair_respects_gauge` — Prop 4.2 sentence 2, now a theorem: the repair
  dynamics consume only declared overlap data, so `obsMap ∘ Repair` factors
  through `obsMap`.
* **Non-degeneracy receipts** (the failure mode this file always warned
  about — `Repair := id` closing the congruence "for the wrong reason" — is
  ruled out by theorems): `lyapunovDescent_holds` and `termination_holds`
  discharge the file's own `LyapunovDescent` and `Termination` obligations
  for the constructed operator (both would be vacuous-false-or-empty under a
  degenerate repair with a nonempty step relation, and the step relation is
  provably nonempty: `acceptedStep_demoCarrier_nonempty`), and
  `Repair_eq_self_of_consistent` shows repair does nothing on already
  consistent records.

## Honest scope (what "admission-free" does and does not close)

`sorry`-free is **not** "repair theory fully closed". Explicitly:

* `Completeness` for the constructed operator is **not** claimed. The
  operator fires only when a site can satisfy all of its overlaps at once,
  so on *frustrated* carriers (a broken incident edge but no single-site
  fix) a broken record can be a normal form. `Completeness` holds exactly on
  frustration-free dynamics; that is the content of the conditional `H1`–`H3`
  development below, which remains the general statement.
* `Confluence` is **not** claimed — and is false in general:
  `demoCarrier_not_confluent` below exhibits a non-confluent instance. The
  constructed `Repair` is one canonical (choice-selected) schedule; its
  gauge-congruence is what Prop 4.2 sentence 2 requires, not schedule
  independence.
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
    patch. A faithful, non-vacuous index type. -/
def Site : Type := C.Patch

/-! ### The repair kernel

The single-site notions the constructed repair operator is built from: edge
consistency, the broken-edge set and its `ℕ` count (the well-founded shadow
of `Φ` used for termination), the local firing trigger, single-site
transactional solvability, and the gauge-congruence lemmas showing all of
them are functions of the declared overlap data `obsMap` alone. -/

section RepairKernel

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

/-- Edge consistency is a property of the declared overlap data alone: it
    reads exactly the two `obsMap` components of the edge. -/
theorem edgeConsistentAt_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (e : C.Edge) :
    edgeConsistentAt e x ↔ edgeConsistentAt e y := by
  have hx := congrFun h e
  have h1 : C.projSrc e (x (C.src e)) = C.projSrc e (y (C.src e)) :=
    congrArg Prod.fst hx
  have h2 : C.projTgt e (x (C.tgt e)) = C.projTgt e (y (C.tgt e)) :=
    congrArg Prod.snd hx
  unfold edgeConsistentAt
  rw [h1, h2]

/-- Updating gauge-equivalent records at the same patch with the same state
    yields gauge-equivalent records: the declared overlap data of the update
    depends only on the new state and the old overlap data. -/
theorem obsMap_update_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : C.Patch) (s : C.State i) :
    obsMap C (Function.update x i s) = obsMap C (Function.update y i s) := by
  funext e
  have hx := congrFun h e
  have h1 : C.projSrc e (x (C.src e)) = C.projSrc e (y (C.src e)) :=
    congrArg Prod.fst hx
  have h2 : C.projTgt e (x (C.tgt e)) = C.projTgt e (y (C.tgt e)) :=
    congrArg Prod.snd hx
  unfold obsMap
  congr 1
  · rcases eq_or_ne (C.src e) i with hsrc | hsrc
    · subst hsrc
      simp only [Function.update_self]
    · simp only [Function.update_of_ne hsrc]
      exact h1
  · rcases eq_or_ne (C.tgt e) i with htgt | htgt
    · subst htgt
      simp only [Function.update_self]
    · simp only [Function.update_of_ne htgt]
      exact h2

/-- The local firing trigger at site `i`: some edge incident to `i` is
    broken. (Same incidence spelling as the `H2`/`H3` laws below.) -/
def LocalTrigger (i : C.Patch) (x : Records C) : Prop :=
  ∃ e : C.Edge, (C.src e = i ∨ C.tgt e = i) ∧ ¬ edgeConsistentAt e x

/-- `s` is a transactional repair for site `i` at `x`: installing it makes
    every edge incident to `i` consistent at once. -/
def SolvesAt (i : C.Patch) (x : Records C) (s : C.State i) : Prop :=
  ∀ e : C.Edge, (C.src e = i ∨ C.tgt e = i) →
    edgeConsistentAt e (Function.update x i s)

/-- Local satisfiability (frustration-freeness at `(i, x)`): *some*
    replacement state for patch `i` satisfies all of `i`'s declared overlaps
    at once. On frustrated instances this fails and the local move honestly
    does not fire — no single-site move can repair an unsatisfiable
    neighbourhood. -/
def LocallySolvable (i : C.Patch) (x : Records C) : Prop :=
  ∃ s : C.State i, SolvesAt i x s

/-- The firing trigger is a function of the declared overlap data alone. -/
theorem localTrigger_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : C.Patch) :
    LocalTrigger i x ↔ LocalTrigger i y :=
  exists_congr fun e => and_congr_right fun _ =>
    not_congr (edgeConsistentAt_congr h e)

/-- The transactional-repair predicate is a function of the declared overlap
    data alone (pointwise in the candidate state `s`). -/
theorem solvesAt_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : C.Patch) (s : C.State i) :
    SolvesAt i x s ↔ SolvesAt i y s :=
  forall_congr' fun e => imp_congr_right fun _ =>
    edgeConsistentAt_congr (obsMap_update_congr h i s) e

/-- Local satisfiability is a function of the declared overlap data alone. -/
theorem locallySolvable_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : C.Patch) :
    LocallySolvable i x ↔ LocallySolvable i y :=
  exists_congr (solvesAt_congr h i)

/-- `Classical.choose` picks the *same* witness from pointwise-equivalent
    predicates. This is what turns the `_congr` lemmas above into on-the-nose
    equalities of repaired states: gauge-equivalent records feed `choose`
    literally equal predicates, so the chosen repair (and the chosen firing
    site in `Repair`) coincide. -/
private theorem choose_eq_of_pred_iff {α : Sort*} {p q : α → Prop}
    (hpq : ∀ a, p a ↔ q a) (hp : ∃ a, p a) (hq : ∃ a, q a) :
    Classical.choose hp = Classical.choose hq := by
  have hpq' : p = q := funext fun a => propext (hpq a)
  subst hpq'
  rfl

end RepairKernel

open Classical in
/-- One transactional/local recovery move at a repair site (*Reality* line
    297, "built from local recovery moves"). Site `i` fires exactly when an
    incident overlap is broken (`LocalTrigger`) *and* it can transactionally
    satisfy all of its declared overlaps at once (`LocallySolvable`); it then
    installs a repaired state chosen from the declared overlap data. On
    frustrated or already-locally-consistent sites it is the identity. -/
noncomputable def localRepair : Site C → Records C → Records C := fun i x =>
  if h : LocalTrigger i x ∧ LocallySolvable i x then
    Function.update x i (Classical.choose h.2)
  else x

theorem localRepair_of_fire (i : Site C) (x : Records C)
    (h : LocalTrigger i x ∧ LocallySolvable i x) :
    localRepair C i x = Function.update x i (Classical.choose h.2) :=
  dif_pos h

theorem localRepair_of_quiescent (i : Site C) (x : Records C)
    (h : ¬ (LocalTrigger i x ∧ LocallySolvable i x)) :
    localRepair C i x = x :=
  dif_neg h

/-- `H1` for the constructed move: firing at `i` changes patch `i` only. -/
theorem localRepair_apply_of_ne (i : Site C) (x : Records C)
    (j : C.Patch) (hj : j ≠ i) :
    localRepair C i x j = x j := by
  by_cases h : LocalTrigger i x ∧ LocallySolvable i x
  · rw [localRepair_of_fire C i x h, Function.update_of_ne hj]
  · rw [localRepair_of_quiescent C i x h]

/-- `H3` for the constructed move: when the move at `i` fires, it makes all
    of `i`'s incident edges consistent. -/
theorem localRepair_repairs (i : Site C) (x : Records C)
    (hfire : localRepair C i x ≠ x) :
    ∀ e : C.Edge, (C.src e = i ∨ C.tgt e = i) →
      edgeConsistentAt e (localRepair C i x) := by
  by_cases h : LocalTrigger i x ∧ LocallySolvable i x
  · intro e hinc
    rw [localRepair_of_fire C i x h]
    exact Classical.choose_spec h.2 e hinc
  · exact absurd (localRepair_of_quiescent C i x h) hfire

/-- Exact firing characterisation: the move at `i` changes `x` iff an
    incident edge is broken *and* the site is locally solvable. (The `H2`
    trigger law holds in the forward direction unconditionally, and as an
    iff exactly on frustration-free instances — see the honest-scope note in
    the header.) -/
theorem localRepair_ne_iff (i : Site C) (x : Records C) :
    localRepair C i x ≠ x ↔ (LocalTrigger i x ∧ LocallySolvable i x) := by
  constructor
  · intro hne
    by_contra hcond
    exact hne (localRepair_of_quiescent C i x hcond)
  · intro h heq
    obtain ⟨e₀, hinc₀, hbrk₀⟩ := h.1
    apply hbrk₀
    have hcons := Classical.choose_spec h.2 e₀ hinc₀
    rwa [← localRepair_of_fire C i x h, heq] at hcons

/-- Whether site `i` fires is a function of the declared overlap data. -/
theorem localRepair_fire_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : Site C) :
    localRepair C i x ≠ x ↔ localRepair C i y ≠ y := by
  rw [localRepair_ne_iff, localRepair_ne_iff]
  exact and_congr (localTrigger_congr h i) (locallySolvable_congr h i)

/-- The single-site move is a gauge congruence: on gauge-equivalent inputs it
    installs the *same* chosen repair (the choice reads only overlap data),
    so the outputs are gauge-equivalent. This is the single-step engine of
    `repair_respects_gauge`. -/
theorem obsMap_localRepair_congr {x y : Records C}
    (h : obsMap C x = obsMap C y) (i : Site C) :
    obsMap C (localRepair C i x) = obsMap C (localRepair C i y) := by
  by_cases hx : LocalTrigger i x ∧ LocallySolvable i x
  · have hy : LocalTrigger i y ∧ LocallySolvable i y :=
      ⟨(localTrigger_congr h i).mp hx.1, (locallySolvable_congr h i).mp hx.2⟩
    have hs : Classical.choose hx.2 = Classical.choose hy.2 :=
      choose_eq_of_pred_iff (solvesAt_congr h i) hx.2 hy.2
    rw [localRepair_of_fire C i x hx, localRepair_of_fire C i y hy, ← hs]
    exact obsMap_update_congr h i (Classical.choose hx.2)
  · have hy : ¬ (LocalTrigger i y ∧ LocallySolvable i y) := fun hy' =>
      hx ⟨(localTrigger_congr h i).mpr hy'.1, (locallySolvable_congr h i).mpr hy'.2⟩
    rw [localRepair_of_quiescent C i x hx, localRepair_of_quiescent C i y hy]
    exact h

/-- **Lyapunov descent on the `ℕ` surrogate.** Every genuine firing strictly
    shrinks the broken-edge set: incident edges are repaired (`H3`),
    non-incident edges are untouched (`H1`), and the trigger's broken edge
    leaves the set. This is what makes `Repair`'s recursion well-founded. -/
theorem mismatchCount_localRepair_lt (i : Site C) (x : Records C)
    (hfire : localRepair C i x ≠ x) :
    mismatchCount (localRepair C i x) < mismatchCount x := by
  have hsub : brokenSet (localRepair C i x) ⊆ brokenSet x := by
    intro e he
    by_cases hinc : C.src e = i ∨ C.tgt e = i
    · exact absurd (localRepair_repairs C i x hfire e hinc)
        (mem_brokenSet_iff_not_consistent.1 he)
    · have hs : C.src e ≠ i := fun hh => hinc (Or.inl hh)
      have ht : C.tgt e ≠ i := fun hh => hinc (Or.inr hh)
      rw [mem_brokenSet, localRepair_apply_of_ne C i x _ hs,
        localRepair_apply_of_ne C i x _ ht] at he
      exact mem_brokenSet.2 he
  obtain ⟨e₀, hinc₀, hbrk₀⟩ := ((localRepair_ne_iff C i x).1 hfire).1
  have hmem : e₀ ∈ brokenSet x := mem_brokenSet_iff_not_consistent.2 hbrk₀
  have hnot : e₀ ∉ brokenSet (localRepair C i x) := fun hm =>
    mem_brokenSet_iff_not_consistent.1 hm (localRepair_repairs C i x hfire e₀ hinc₀)
  exact Finset.card_lt_card ((Finset.ssubset_iff_of_subset hsub).2 ⟨e₀, hmem, hnot⟩)

open Classical in
/-- The composite repair operator: fire a (choice-canonical) firing site,
    repeat until no site fires. One concrete asynchronous schedule composed
    to a normal form (*Reality* line 297 / `OPHConsensus`), well-founded by
    `mismatchCount` descent. `Repair_reachable`/`Repair_normalForm` below
    certify it reaches a genuine `acceptedStep`-normal form; it is *not*
    claimed to be schedule-independent (`demoCarrier_not_confluent`). -/
noncomputable def Repair (x : Records C) : Records C :=
  if h : ∃ i : Site C, localRepair C i x ≠ x then
    Repair (localRepair C (Classical.choose h) x)
  else x
termination_by mismatchCount x
decreasing_by
  exact mismatchCount_localRepair_lt C (Classical.choose h) x (Classical.choose_spec h)

theorem Repair_of_fire (x : Records C)
    (h : ∃ i : Site C, localRepair C i x ≠ x) :
    Repair C x = Repair C (localRepair C (Classical.choose h) x) := by
  conv_lhs => rw [Repair]
  exact dif_pos h

theorem Repair_of_normal (x : Records C)
    (h : ¬ ∃ i : Site C, localRepair C i x ≠ x) :
    Repair C x = x := by
  conv_lhs => rw [Repair]
  exact dif_neg h

/-- One accepted asynchronous repair step: some site's local move changes
    the record. This is the relation the generic abstract-rewriting
    skeleton instantiates after the concrete repair layer is supplied. -/
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

/-- Engine of `repair_respects_gauge`, by strong induction on the
    broken-edge count: on `obsMap`-equal inputs, either no site fires on
    both (both `Repair`s are the identity), or the *same* site is chosen on
    both (`choose_eq_of_pred_iff` over the fire predicate, which is a
    function of overlap data by `localRepair_fire_congr`), the single steps
    stay `obsMap`-equal (`obsMap_localRepair_congr`), and the count strictly
    drops (`mismatchCount_localRepair_lt`). -/
private theorem obsMap_Repair_congr_aux :
    ∀ (n : ℕ) (x y : Records C), mismatchCount x = n →
      obsMap C x = obsMap C y →
        obsMap C (Repair C x) = obsMap C (Repair C y) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro x y hn h
    by_cases hx : ∃ i : Site C, localRepair C i x ≠ x
    · have hy : ∃ i : Site C, localRepair C i y ≠ y :=
        let ⟨i, hi⟩ := hx
        ⟨i, (localRepair_fire_congr C h i).mp hi⟩
      have hsite : Classical.choose hx = Classical.choose hy :=
        choose_eq_of_pred_iff (fun i => localRepair_fire_congr C h i) hx hy
      rw [Repair_of_fire C x hx, Repair_of_fire C y hy, ← hsite]
      exact ih (mismatchCount (localRepair C (Classical.choose hx) x))
        (hn ▸ mismatchCount_localRepair_lt C (Classical.choose hx) x
          (Classical.choose_spec hx))
        _ _ rfl (obsMap_localRepair_congr C h (Classical.choose hx))
    · have hy : ¬ ∃ i : Site C, localRepair C i y ≠ y := fun hy' =>
        hx (let ⟨i, hi⟩ := hy'; ⟨i, (localRepair_fire_congr C h i).mpr hi⟩)
      rw [Repair_of_normal C x hx, Repair_of_normal C y hy]
      exact h

/-- `∼_gauge` is a `Repair`-congruence — Prop 4.2 sentence 2 (independence
    on the physical quotient), now a theorem.

    The proof is exactly the sentence's content: the constructed async
    `Repair` factors through `obsMap`. Every datum the dynamics consult —
    which edges are broken, whether a site can transactionally repair its
    neighbourhood, which repaired state is installed, which site fires
    next — is a function of the declared overlap data, so gauge-equivalent
    records evolve through pointwise gauge-equivalent trajectories
    (`obsMap_Repair_congr_aux`). This is *not* closed "for the wrong
    reason": the operator genuinely fires (`acceptedStep_demoCarrier_nonempty`),
    strictly descends `Φ` (`lyapunovDescent_holds`), and reaches genuine
    normal forms (`Repair_normalForm`). -/
theorem repair_respects_gauge :
    ∀ x y : Records C, gaugeEquiv C x y → gaugeEquiv C (Repair C x) (Repair C y) :=
  fun x y h => obsMap_Repair_congr_aux C (mismatchCount x) x y rfl h

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

/-! ## Non-degeneracy receipts for the constructed repair

The header of this file always named the failure mode of a fake discharge:
a degenerate `Repair` (`id` / a constant) closes `repair_respects_gauge`
"for the wrong reason" while making the dynamical obligations vacuous. The
theorems below rule that out for the constructed operator:

* `Repair_normalForm` / `Repair_reachable` — `Repair` reaches a genuine
  normal form of `acceptedStep` by a genuine accepted asynchronous run.
* `lyapunovDescent_holds` / `termination_holds` — the file's own
  `LyapunovDescent` and `Termination` obligations, discharged.
* `acceptedStep_demoCarrier_nonempty` — the accepted-step relation is
  provably nonempty (the operator really fires), so none of the above is a
  statement about an empty relation.
* `Repair_eq_self_of_consistent` — repair does nothing on consistent
  records. -/

private theorem Repair_normalForm_aux :
    ∀ (n : ℕ) (x : Records C), mismatchCount x = n → NormalForm C (Repair C x) := by
  intro n
  induction' n using Nat.strong_induction_on with n ih
  intro x hn
  by_cases hx : ∃ i : Site C, localRepair C i x ≠ x
  · rw [Repair_of_fire C x hx]
    exact ih (mismatchCount (localRepair C (Classical.choose hx) x))
      (hn ▸ mismatchCount_localRepair_lt C (Classical.choose hx) x
        (Classical.choose_spec hx))
      _ rfl
  · rw [Repair_of_normal C x hx]
    intro y hstep
    obtain ⟨i, _, hfire⟩ := hstep
    exact hx ⟨i, hfire⟩

/-- **Receipt — `Repair` reaches a genuine normal form.** No accepted repair
    step applies to `Repair C x`. -/
theorem Repair_normalForm (x : Records C) : NormalForm C (Repair C x) :=
  Repair_normalForm_aux C (mismatchCount x) x rfl

private theorem Repair_reachable_aux :
    ∀ (n : ℕ) (x : Records C), mismatchCount x = n →
      ReflTransGen (acceptedStep C) x (Repair C x) := by
  intro n
  induction' n using Nat.strong_induction_on with n ih
  intro x hn
  by_cases hx : ∃ i : Site C, localRepair C i x ≠ x
  · rw [Repair_of_fire C x hx]
    exact ReflTransGen.head
      ⟨Classical.choose hx, rfl, Classical.choose_spec hx⟩
      (ih (mismatchCount (localRepair C (Classical.choose hx) x))
        (hn ▸ mismatchCount_localRepair_lt C (Classical.choose hx) x
          (Classical.choose_spec hx))
        _ rfl)
  · rw [Repair_of_normal C x hx]

/-- **Receipt — `Repair` is an accepted asynchronous run.** `Repair C x` is
    reached from `x` by accepted repair steps, i.e. the operator is the
    composition of local recovery moves under one asynchronous schedule —
    the paper's construction, not an unrelated function that happens to
    satisfy the congruence. -/
theorem Repair_reachable (x : Records C) :
    ReflTransGen (acceptedStep C) x (Repair C x) :=
  Repair_reachable_aux C (mismatchCount x) x rfl

/-- **Receipt — repair fixes consistent records.** On `Φ = 0` records no
    site fires and `Repair` is the identity. -/
theorem Repair_eq_self_of_consistent (x : Records C) (hx : Consistent C x) :
    Repair C x = x := by
  apply Repair_of_normal
  rintro ⟨i, hfire⟩
  obtain ⟨e₀, _, hbrk₀⟩ := ((localRepair_ne_iff C i x).1 hfire).1
  exact hbrk₀ ((consistent_iff_edgeConsistent C x).1 hx e₀)

/-- **Receipt — the file's `LyapunovDescent` obligation, discharged.** Every
    accepted step strictly lowers `Φ`: the fired site's incident edges drop
    to zero mismatch (one of them was strictly positive by
    `weight_pos`/`dist_eq_zero`), all other edges are untouched. A degenerate
    `Repair` was rejected in this file precisely because it would make this
    obligation vacuous; for the constructed operator it is a theorem over a
    provably nonempty step relation (`acceptedStep_demoCarrier_nonempty`). -/
theorem lyapunovDescent_holds : LyapunovDescent C := by
  intro x y hstep
  obtain ⟨i, rfl, hfire⟩ := hstep
  unfold Φ
  apply Finset.sum_lt_sum
  · intro e _
    by_cases hinc : C.src e = i ∨ C.tgt e = i
    · have hcons := localRepair_repairs C i x hfire e hinc
      rw [(edgeConsistentAt_iff_dist e _).1 hcons, mul_zero]
      exact zero_le _
    · have hs : C.src e ≠ i := fun hh => hinc (Or.inl hh)
      have ht : C.tgt e ≠ i := fun hh => hinc (Or.inr hh)
      rw [localRepair_apply_of_ne C i x _ hs, localRepair_apply_of_ne C i x _ ht]
  · obtain ⟨e₀, hinc₀, hbrk₀⟩ := ((localRepair_ne_iff C i x).1 hfire).1
    refine ⟨e₀, Finset.mem_univ e₀, ?_⟩
    have hcons := localRepair_repairs C i x hfire e₀ hinc₀
    rw [(edgeConsistentAt_iff_dist e₀ _).1 hcons, mul_zero]
    have hd : C.dist e₀ (C.projSrc e₀ (x (C.src e₀))) (C.projTgt e₀ (x (C.tgt e₀))) ≠ 0 :=
      fun h0 => hbrk₀ ((edgeConsistentAt_iff_dist e₀ x).2 h0)
    exact mul_pos (C.weight_pos e₀) (pos_iff_ne_zero.mpr hd)

/-- **Receipt — the file's `Termination` obligation, discharged.** The
    accepted asynchronous-repair relation for the constructed `localRepair`
    is well-founded, via the `mismatchCount` measure. -/
theorem termination_holds : Termination C :=
  have H : Subrelation (fun y x : Records C => acceptedStep C x y)
      (InvImage (· < ·) mismatchCount) := fun {y x} hxy => by
    obtain ⟨i, rfl, hfire⟩ := hxy
    exact mismatchCount_localRepair_lt C i x hfire
  Subrelation.wf H (InvImage.wf _ wellFounded_lt)

/-- On `demoCarrier`, site `false` genuinely fires from the identity record:
    the single edge is broken and copying the neighbour's value (`true`)
    transactionally repairs it. Proved through the firing characterisation
    `localRepair_ne_iff` — no need to compute the classical choice. -/
theorem localRepair_demoCarrier_fires :
    localRepair demoCarrier false (fun b => b) ≠ (fun b => b) := by
  rw [localRepair_ne_iff]
  constructor
  · refine ⟨(), Or.inl rfl, ?_⟩
    show ¬ ((false : Bool) = (true : Bool))
    decide
  · refine ⟨true, fun e _ => ?_⟩
    show Function.update (fun b : Bool => b) false true false =
      Function.update (fun b : Bool => b) false true true
    rw [Function.update_self, Function.update_of_ne (by decide : (true : Bool) ≠ false)]

/-- **Receipt — the accepted-step relation is nonempty.** The constructed
    repair genuinely fires, so the dynamical receipts above are not
    statements about an empty relation. -/
theorem acceptedStep_demoCarrier_nonempty :
    ∃ x y : Records demoCarrier, acceptedStep demoCarrier x y :=
  ⟨fun b => b, localRepair demoCarrier false (fun b => b),
    false, rfl, localRepair_demoCarrier_fires⟩

/-! ### Axiom audit — the constructed repair layer is admission-free.
The `#print axioms` outputs below confirm that the constructed
`localRepair`/`Repair`, the discharged `repair_respects_gauge`, and every
non-degeneracy receipt depend only on the standard Lean/Mathlib axioms
(`propext`, `Classical.choice`, `Quot.sound`) — no `sorryAx`, no
`native_decide`, no new axiom. -/
#print axioms localRepair
#print axioms Repair
#print axioms repair_respects_gauge
#print axioms Repair_normalForm
#print axioms Repair_reachable
#print axioms Repair_eq_self_of_consistent
#print axioms lyapunovDescent_holds
#print axioms termination_holds
#print axioms localRepair_demoCarrier_fires
#print axioms acceptedStep_demoCarrier_nonempty

end OPH

/-! ## Global termination & completeness from LOCAL repair laws

This section proves the mathematical content of two of OPH's open *dynamical*
obligations — **Termination** and **Completeness** (cf. the `Termination`/
`Completeness` `def`s above) — **conditionally, for any local repair move
satisfying the local laws `H1`/`H2`/`H3` below**, derived from those explicit,
faithful, single-site properties. It establishes the theorems for the abstract
move `lr`, independently of the constructed `localRepair` above. (The file's
own `Termination` `def` is now discharged outright for the constructed
operator — `termination_holds` — while its `Completeness` `def` is *not*
claimed there: the constructed move fires only on locally solvable sites, so
full `Completeness` remains exactly the frustration-free conditional proved
here.) The laws are
satisfiable by a genuine repair (e.g. a two-`Bool`-patch carrier with one
edge, each patch copying its neighbour to snap the edge consistent), so the
result is conditional, not vacuous — and that satisfiability is itself
machine-checked below as `demoCarrier_terminates` (a concrete `(carrier, repair)`
instance discharging `H1`/`H2`/`H3` with a real, non-empty repair step).

It is deliberately **self-contained and axiom-clean**: it does *not* reference
the constructed `localRepair`/`Repair`. The repair move and its laws enter
as `section variable`s (`lr`, `H1`, `H2`, `H3`), so each theorem here closes
with `#print axioms` reporting only `[propext, Classical.choice, Quot.sound]`
(no `sorryAx`, no new `axiom`). The statements are phrased over the
hypothesis-bearing move `lr` (`acceptedStepLR`, `NormalFormLR`) and are
mathematically the same theorems for any concrete operator that satisfies
`H1`/`H2`/`H3` — including the constructed `localRepair` on frustration-free
carriers.

### Hypotheses are LOCAL; conclusions are GLOBAL (no assume-the-conclusion)

The hypotheses are genuine **single-site** statements:
* `H1` (`lr` changes only site `i`): firing at `i` touches patch `i` only.
* `H2` (`lr` fires iff a local edge is broken): the move at `i` changes `x`
  *iff* some edge incident to `i` is locally inconsistent — a purely local
  trigger.
* `H3` (local satisfiability / frustration-freeness): when the move at `i`
  fires it makes *all* of `i`'s own incident edges consistent. This explicitly
  restricts to carriers where a single patch *can* satisfy all its overlaps at
  once (frustration-free); it is a local property, **not** the global claim.

The conclusions are **global** dynamical facts about all of `Records C`:
* `termination`: the asynchronous accepted-step relation is `WellFounded`.
* `completeness`: a record is a global normal form *iff* it is globally
  `Consistent` (`Φ = 0`).

None of the forbidden shortcuts is assumed: we never assume `mismatchCount`
decreases, nor `WellFounded`, nor `Termination`, nor `NormalForm ↔ Consistent`.
Those are *proved* from the three local laws (plus the discharged
`consistent_iff_edgeConsistent`).

### The Lyapunov / Inter-Basin termination pattern

The proof is the well-founded-measure pattern: every accepted repair strictly
lowers a **structural `ℕ` surrogate** `mismatchCount` (the number of broken
edges), exactly as every SKI reduction strictly lowers `basin_size` in the
Inter-Basin termination theorem. A `ℕ` surrogate is *needed* because the
carrier potential `Φ : ℝ≥0` is **not** `<`-well-founded; `mismatchCount` is the
well-founded shadow of `Φ` that makes asynchronous descent terminate.

### What remains open (explicit scoping; no `sorry`)

`Confluence`/`LocallyConfluent` is **not** provided: asynchronous repairs at
different sites need not commute (`lr i (lr j x)` and `lr j (lr i x)` can
differ), so a frustration-free carrier may reach distinct normal forms
under distinct schedules — schedule independence / unique normal forms is out
of scope for these hypotheses. There is no `sorry`, `admit`, or new `axiom`
anywhere in this section. -/

namespace OPH

open Relation  -- `ReflTransGen` (used by the confluence theorems below)

section LocalRepairDynamics

variable {C : OPHCarrier}

-- `edgeConsistentAt`, `brokenSet`, `mismatchCount` and their membership
-- lemmas now live in the RepairKernel section above (they are shared with
-- the constructed `localRepair`/`Repair`); this section keeps using them.

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
-- iff some edge incident to i is locally inconsistent.
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
    (cf. `demoCarrier_not_confluent`).
    SCOPE (explicit): `HB` and `Hfib` are jointly satisfiable in principle, but NO single carrier in this
    file instantiates BOTH — `demoBoundary` has `HB` (`demoBoundary_HB`) but fails `Hfib`
    (`demoCarrier_Hfib_fails`); `obsMap`/`demoSeedBoundary` give `Hfib` but are not `demoLR`-preserved.
    The two premises are exhibited on separate carriers; a joint witness on a richer multi-edge carrier
    is open modeling task. So this is an explicitly-scoped conditional, not vacuous. -/
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
-- explicit weaker hypothesis would restrict commutation to NON-INCIDENT pairs
-- (sites sharing no edge, expressible via the incidence predicate used in
-- H2/H3), but that alone does not close the diamond when incident repairs
-- genuinely diverge. H4 is NOT implied by H1–H3.
variable
  (H4 : ∀ (i j : C.Patch) (x : Records C), lr i (lr j x) = lr j (lr i x))

-- H4 is used only inside the proofs below (the statements quantify over the
-- abstract relation `acceptedStepLR lr`), so force its inclusion explicitly.
include H4

-- H1/H2/H3 are unused by the diamond argument (it needs only `H4`); drop them
-- from THIS lemma so its type explicitly reads "commuting moves are locally
-- confluent". `confluence_of_commute` below carries them (for `termination`).
omit H1 H2 H3 in
/-- **Local confluence from single-step commutation** (the diamond condition).
    From two accepted steps at sites `i`, `j`, the common join is
    `lr j (lr i x) = lr i (lr j x)` (by `H4`); each side reaches it in ≤ 1 step
    (zero if that site is quiescent there). -/
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
the SAME two normal forms `demoCarrier_not_confluent` exhibits: `(1,1)` and `(0,0)`.

EXPLICIT SCOPE: this is the FORWARD direction (`Hfib` ⟹ unique; symmetric ⟹ countermodel),
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

/-! ### Two DIFFERENT uniqueness notions, two DIFFERENT levers (the corrected #304 cut)

`demoCarrier_Hfib_fails` shows the symmetric repair + trivial boundary breaks OBSERVER-uniqueness
("same boundary → same normal form"). Two levers act, on two DIFFERENT notions — do not conflate them:
  ROUTE A — refine the BOUNDARY so the consistent fiber is a gauge-singleton (`Hfib` holds). This buys
            OBSERVER-uniqueness directly; a property of B, repair-free (`demoCarrier_Hfib_holds_seed`).
  ROUTE B — use a SELECTING (deterministic) repair: it is CONFLUENT — a unique normal form per INPUT,
            schedule-independent (Church-Rosser; last-writer-wins / strong-eventual-consistency). This
            is per-INPUT uniqueness, a property of the repair dynamics — and it is NOT observer-
            uniqueness: under a coarse boundary the SAME confluent repair fails it
            (`demoCarrier_dir_not_observer_unique`). Observer-uniqueness returns only once the boundary
            is also refined (`demoCarrier_dir_observer_unique_under_seed`) — i.e. via Route A. -/

/-- **ROUTE A (DEFINITIONAL endpoint) — at the finest boundary `Hfib` holds verbatim.** With
    `B := obsMap`, `Hfib`'s conclusion `gaugeEquiv x y` IS its hypothesis `obsMap x = obsMap y`
    (`gaugeEquiv` unfolds definitionally to `obsMap`-equality), so the proof is `fun _ _ h _ _ => h`
    and discards BOTH `Consistent` premises. This endpoint therefore CANNOT fail and is NOT independent
    evidence — it is the trivial top of the boundary lattice, recorded only for completeness. The
    load-bearing positive result is `demoCarrier_Hfib_holds_seed` (a strictly coarser boundary that
    genuinely uses consistency). `Hfib` is a property of B, not of the repair (which is why "selecting
    repair proves Hfib" was a category error). -/
theorem demoCarrier_Hfib_holds_finerB :
    ∀ x y : Records demoCarrier, obsMap demoCarrier x = obsMap demoCarrier y →
      Consistent demoCarrier x → Consistent demoCarrier y → gaugeEquiv demoCarrier x y :=
  fun _ _ h _ _ => h

/-- The SEED boundary: read a SINGLE cell (patch `false`). A coarse boundary — strictly between
    the trivial `demoBoundary` (reads nothing) and the full `obsMap` (reads the whole overlap). -/
def demoSeedBoundary : Records demoCarrier → Bool := fun x => x false

/-- **SEED base case (issue #304, Boundary-Fiber Confluence).** `Hfib` holds for the one-cell seed
    boundary: reading patch `false` PLUS the consistency premise pins the whole record, because
    edge-agreement (`consistent_iff_edgeConsistent`) forces the unread cell to follow the read one. So
    `Hfib` does NOT require the full observable — a single seed cell suffices. The proof genuinely
    consumes both `Consistent` hypotheses (unlike `demoCarrier_Hfib_holds_finerB`, which is
    definitional), so this is real singleton-consistent-fiber content. SCOPE (explicit): on this 1-edge
    carrier every `Consistent` record is constant, so the "bulk" is just the one other cell, and on the
    consistent fiber `demoSeedBoundary` has the SAME separating power as the full `obsMap` — genuine
    seed-determines-record, but the multi-edge *propagation* the name evokes is open modeling task, not
    exhibited here. (`demoSeedBoundary` is the boundary-fineness criterion only; it is NOT preserved by
    the symmetric `demoLR` — pairing it with a boundary-preserving repair, the `HB` premise, is a
    separate modeling step. It IS preserved by the selecting `demoDirT`, which is how
    `demoCarrier_dir_observer_unique_under_seed` uses it.) -/
theorem demoCarrier_Hfib_holds_seed :
    ∀ x y : Records demoCarrier, demoSeedBoundary x = demoSeedBoundary y →
      Consistent demoCarrier x → Consistent demoCarrier y → gaugeEquiv demoCarrier x y := by
  intro x y hseed hcx hcy
  simp only [demoSeedBoundary] at hseed
  rw [consistent_iff_edgeConsistent] at hcx hcy
  have hx : x false = x true := hcx ()
  have hy : y false = y true := hcy ()
  have hxy : x = y := by
    funext k
    cases k
    · exact hseed
    · rw [← hx, hseed, hy]
  subst hxy
  rfl

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

/-- **THE CRUX — confluence ALONE is not observer-uniqueness.** `demoDirT` is confluent
    (`demoCarrier_dir_confluent`: same input → one normal form). But that is a property of the REPAIR,
    not of the boundary. Here are TWO inputs with the SAME *trivial* boundary (`demoBoundary : … → Unit`,
    which reads nothing, so equality holds for ALL inputs by `rfl`) that reach DIFFERENT normal forms
    `(1,1)` and `(0,0)` under that very confluent repair — so confluence ("same input → same NF") does
    NOT by itself give observer-facing uniqueness ("same boundary → same NF"). This is a separation,
    exhibited against the COARSEST boundary; it does NOT claim the repair's directionality is irrelevant
    in general. The controlled converse is machine-checked next: hold THIS SAME confluent repair fixed
    and REFINE the boundary to the one-cell `demoSeedBoundary`, and observer-uniqueness is RESTORED
    (`demoCarrier_dir_observer_unique_under_seed`). So the lever that flips reconstructability — with the
    repair held constant — is the boundary's fineness. -/
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

/-- **POSITIVE COMPANION TO THE CRUX — the boundary is the controlling lever, machine-checked.**
    Hold the SAME confluent repair `demoDirT` fixed and REFINE the boundary from the trivial
    `demoBoundary` (under which `demoCarrier_dir_not_observer_unique` shows observer-uniqueness FAILS)
    to the one-cell `demoSeedBoundary`: observer-uniqueness is RESTORED — any two inputs with equal
    seed-boundary reach `gaugeEquiv` normal forms under that very confluent repair. So with the repair
    held constant, varying ONLY the boundary flips reconstructability; the lever is the boundary's
    fineness, not the repair's directionality. Proof: `demoDirT` snaps every record to `(x false, x false)`,
    so it (i) preserves `demoSeedBoundary` along every reduction and (ii) has only `Consistent`
    (edge-constant) normal forms; then equal-seed-boundary consistent records are `gaugeEquiv` by
    `demoCarrier_Hfib_holds_seed`. -/
theorem demoCarrier_dir_observer_unique_under_seed :
    ∀ x y nfx nfy : Records demoCarrier,
      Relation.ReflTransGen (AbstractRewriting.stepRel demoDirT) x nfx →
      AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) nfx →
      Relation.ReflTransGen (AbstractRewriting.stepRel demoDirT) y nfy →
      AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) nfy →
      demoSeedBoundary x = demoSeedBoundary y →
      gaugeEquiv demoCarrier nfx nfy := by
  -- `demoSeedBoundary` is invariant along every `demoDirT` reduction.
  have seed_inv : ∀ z w : Records demoCarrier,
      Relation.ReflTransGen (AbstractRewriting.stepRel demoDirT) z w →
      demoSeedBoundary w = demoSeedBoundary z := by
    intro z w h
    induction h with
    | refl => rfl
    | tail _ hstep ih =>
        obtain ⟨hw, _⟩ := hstep
        rw [hw]
        simpa only [demoSeedBoundary, demoDirT] using ih
  -- A `demoDirT` normal form is a fixed point, hence edge-constant, hence `Consistent`.
  have nf_consistent : ∀ nf : Records demoCarrier,
      AbstractRewriting.IsNormalForm (AbstractRewriting.stepRel demoDirT) nf →
      Consistent demoCarrier nf := by
    intro nf hnf
    have hfix : demoDirT nf = nf := by
      by_contra hne
      exact hnf (demoDirT nf) ⟨rfl, hne⟩
    have hcell : nf false = nf true := by
      have h2 := congrFun hfix true
      simpa only [demoDirT] using h2
    rw [consistent_iff_edgeConsistent]
    intro _; exact hcell
  intro x y nfx nfy hx hxn hy hyn hseed
  have hcx : Consistent demoCarrier nfx := nf_consistent nfx hxn
  have hcy : Consistent demoCarrier nfy := nf_consistent nfy hyn
  have hbb : demoSeedBoundary nfx = demoSeedBoundary nfy := by
    rw [seed_inv x nfx hx, hseed, ← seed_inv y nfy hy]
  exact demoCarrier_Hfib_holds_seed nfx nfy hbb hcx hcy

/-! ### Axiom audit — the reconstruction layer depends only on standard axioms.
The `#print axioms` outputs below confirm that the boundary-fiber reconstruction theorem
and all its concrete witnesses depend ONLY on the standard Lean/Mathlib axioms
(`propext`, `Classical.choice`, `Quot.sound`) — the same footprint as the constructed
repair layer audited above; the file's formerly-declared `sorry`s are discharged, so the
"machine-checked" claim for observer-reconstruction carries no admissions anywhere. -/
#print axioms boundary_fiber_observer_unique
#print axioms boundary_preserved_reduction
#print axioms demoCarrier_Hfib_fails
#print axioms demoCarrier_Hfib_holds_finerB
#print axioms demoCarrier_Hfib_holds_seed
#print axioms demoCarrier_dir_confluent
#print axioms demoCarrier_dir_not_observer_unique
#print axioms demoCarrier_dir_observer_unique_under_seed
#print axioms termination
#print axioms completeness

end OPH
