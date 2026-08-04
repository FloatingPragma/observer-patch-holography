import Mathlib
import ObserverPatchHolography.Primitives
import ObserverPatchHolography.RateNonidentifiability

/-!
# The rate no-go does not reach the OPH repair layer: a located obstruction

`ObserverPatchHolography.RateNonidentifiability` proves an abstract non-identifiability
result: for a deterministic map `T : S → S`, a stutter extension leaves the declared
reading (locked set, accepted-repair reachability, basin) fixed while multiplying the
first-locking count by an arbitrary factor. This file asks the only question that makes
that result repository-facing, **does it instantiate at the repository's own observable
map and accepted-repair relation?**, and answers **no**, with the failure located at a
single hypothesis.

## The attempted bridge

The repository's own definitions, all in `ObserverPatchHolography.Primitives`:

* the carrier `OPHCarrier` and the state type `Records C` (`Primitives.lean:107`, `147`);
* the declared observable map `obsMap C : Records C → Obs C` (`Primitives.lean:157`);
* the accepted asynchronous repair relation `acceptedStep C` (`Primitives.lean:446`);
* the single-site move `localRepair C` it is built from (`Primitives.lean:315`).

These are the definitions the repository's own results use: `acceptedStep` carries
`NormalForm`, `Confluence`, `Completeness`, `LyapunovDescent` and `Termination`
(`Primitives.lean:456`–`499`), and `obsMap` carries `gaugeEquiv` and
`repair_respects_gauge` (`Primitives.lean:506`, `564`). Nothing here is invented for
this file.

## The obstruction, in one line

The abstract construction needs an accepted step that **moves the state while leaving the
declared observable fixed**. The stutter extension has exactly such steps, and
`RateNonidentifiability.stutter_has_projection_fixing_step` says so openly. The
repository's `acceptedStep C` provably has **none**: `acceptedStep_changes_obs`.

The reason is a chain of facts, all of them in tree:

1. `brokenSet` is defined from `C.dist e (C.projSrc e (x (C.src e))) (C.projTgt e (x (C.tgt e)))`
   (`Primitives.lean:194`), which is exactly the pair `obsMap C x e`. So the broken-edge
   count is a function of the declared observable alone (`mismatchCount_congr`).
2. Every accepted step strictly lowers that count, by the in-tree Lyapunov descent
   `mismatchCount_localRepair_lt` (`Primitives.lean:397`).
3. Hence every accepted step strictly changes the declared observable.

So the extension the abstract theorem builds cannot be an OPH repair dynamics: an
observably-idle accepted step is absent from the construction and impossible in the
repository's model.

## What holds instead: the observable *bounds* the rate, and in fact fixes it

The obstruction is not a gap that better plumbing would close. The concrete layer
satisfies the **opposite** of what the abstract no-go concludes:

* `firstLock_le_mismatchCount`: the number of accepted repair steps to first locking is
  at most `mismatchCount x`, a quantity determined by `obsMap C x` alone. The declared
  observable therefore *bounds* the relaxation time. The abstract theorem's family of
  arbitrarily slow systems with one reading has no counterpart here.
* `firstLock_obs_determined`: two records with equal declared observable data have equal
  first-locking counts. The count is a *function of the declared observable*, which is
  precisely the functional the abstract `no_rate_functional` refutes for its own weaker
  reading.

## Why this is not a contradiction

The abstract `DeclaredReading` (`RateNonidentifiability.lean:364`) is a triple of
relations on the state space. The repository's `obsMap` is a *per-state* observable and is
strictly finer: `obsMap` distinguishes records that the abstract reading cannot separate.
The abstract theorem is true of its own reading. What is established here
is that it does **not** transfer to `obsMap`: transferring it is impossible,
because `obsMap` determines the rate under the `stepOnce` schedule.

## Scope

Everything below concerns sense (2) of "speed" in `RateNonidentifiability`'s header:
relaxation time as a count of accepted repair steps. Nothing here is a physical duration,
a clock map, a spectral gap, or a control-parameter law, and no OPH physics claim is
endorsed.

The one-step scheduler `stepOnce` defined below is the body
of the in-tree `Repair` recursion (`Primitives.lean:423`–`426`) stopped after a single
firing. `stepOnce_acceptedStep` proves every move it makes is an `acceptedStep C`, and
`locked_stepOnce_iff_normalForm` proves its fixed points are exactly the in-tree
`NormalForm C`. It is needed only because `FirstLock` counts iterations of a function
while `acceptedStep` is a relation; see the discussion at `stepOnce`.
-/

open Function Relation OPH.AbstractRewriting OPH.RateNonidentifiability

namespace OPH.RateBridgeObstruction

variable {C : OPHCarrier}

/-! ## Step 1: the declared observable determines the broken-edge count -/

/-- The broken-edge set is a function of the declared observable data alone. `brokenSet`
tests `C.dist e` on exactly the two components of `obsMap C x e`, so records with equal
observable data break on the same edges. -/
theorem brokenSet_congr {x y : Records C} (h : obsMap C x = obsMap C y) :
    brokenSet x = brokenSet y := by
  ext e
  rw [mem_brokenSet_iff_not_consistent, mem_brokenSet_iff_not_consistent]
  exact not_congr (edgeConsistentAt_congr h e)

/-- The mismatch count is a function of the declared observable data alone. -/
theorem mismatchCount_congr {x y : Records C} (h : obsMap C x = obsMap C y) :
    mismatchCount x = mismatchCount y :=
  congrArg Finset.card (brokenSet_congr h)

/-! ## Step 2: the obstruction

The abstract stutter construction requires an accepted step that leaves the declared
projection fixed. The repository's accepted step relation admits none. -/

/-- **THE OBSTRUCTION.** Every accepted asynchronous repair step strictly changes the
declared observable data. Combining the in-tree Lyapunov descent
`mismatchCount_localRepair_lt` with `mismatchCount_congr`: a step that left `obsMap`
fixed would leave `mismatchCount` fixed, and no accepted step does. -/
theorem acceptedStep_changes_obs {x y : Records C} (h : acceptedStep C x y) :
    obsMap C x ≠ obsMap C y := by
  obtain ⟨i, rfl, hfire⟩ := h
  intro hobs
  have h1 : mismatchCount x = mismatchCount (localRepair C i x) := mismatchCount_congr hobs
  have h2 : mismatchCount (localRepair C i x) < mismatchCount x :=
    mismatchCount_localRepair_lt C i x hfire
  omega

/-- The same fact stated as the hypothesis that fails. `stutterStep` has reachable steps
that move an unlocked state while fixing the declared projection
(`RateNonidentifiability.stutter_has_projection_fixing_step`); `acceptedStep C` has no such
step, so no stutter extension of the repository's repair dynamics exists. -/
theorem no_observably_idle_acceptedStep (x y : Records C) :
    ¬ (acceptedStep C x y ∧ obsMap C x = obsMap C y) :=
  fun h => acceptedStep_changes_obs h.1 h.2

/-- **The obstruction is not an artefact of the one constructed operator.** It holds for
the hypothesis-bearing relation `acceptedStepLR` that carries the repository's main
dynamics theorems (`Primitives.termination`, `completeness`, `confluence_of_commute`), for
*every* local move satisfying the declared laws H1–H3. The engine is the in-tree
`mismatchCount_lt` (`Primitives.lean:905`). Any OPH repair dynamics obeying the declared
local laws therefore fails the abstract construction's key hypothesis. -/
theorem acceptedStepLR_changes_obs (lr : C.Patch → Records C → Records C)
    (H1 : ∀ (i : C.Patch) (x : Records C) (j : C.Patch), j ≠ i → lr i x j = x j)
    (H2 : ∀ (i : C.Patch) (x : Records C),
      lr i x ≠ x ↔ ∃ e : C.Edge, (C.src e = i ∨ C.tgt e = i) ∧ ¬ edgeConsistentAt e x)
    (H3 : ∀ (i : C.Patch) (x : Records C), lr i x ≠ x →
      ∀ e : C.Edge, (C.src e = i ∨ C.tgt e = i) → edgeConsistentAt e (lr i x))
    {x y : Records C} (h : acceptedStepLR lr x y) :
    obsMap C x ≠ obsMap C y := by
  intro hobs
  have h1 : mismatchCount x = mismatchCount y := mismatchCount_congr hobs
  have h2 : mismatchCount y < mismatchCount x := mismatchCount_lt lr H1 H2 H3 h
  omega

/-! ## Step 3: a one-step scheduler, so that `FirstLock` has something to count

`FirstLock` counts iterations of a *function*; `acceptedStep C` is a *relation*,
nondeterministic in the firing site. The in-tree `Repair C` is not usable either: it runs
the schedule to a normal form in one application, so it locks after at most one step from
anywhere and carries no rate.

`stepOnce` is the in-tree `Repair` recursion body (`Primitives.lean:423`–`426`) stopped
after a single firing: it fires the same choice-canonical site `Repair` would fire, once.
Its faithfulness is discharged rather than assumed: `stepOnce_acceptedStep` (every move is an
`acceptedStep C`) and `locked_stepOnce_iff_normalForm` (its fixed points are exactly the
in-tree `NormalForm C`). -/

open Classical in
/-- One firing of the choice-canonical site: the body of `Repair`'s recursion, run once. -/
noncomputable def stepOnce (C : OPHCarrier) (x : Records C) : Records C :=
  if h : ∃ i : Site C, localRepair C i x ≠ x then
    localRepair C (Classical.choose h) x
  else x

theorem stepOnce_of_fire (C : OPHCarrier) (x : Records C)
    (h : ∃ i : Site C, localRepair C i x ≠ x) :
    stepOnce C x = localRepair C (Classical.choose h) x :=
  dif_pos h

theorem stepOnce_of_normal (C : OPHCarrier) (x : Records C)
    (h : ¬ ∃ i : Site C, localRepair C i x ≠ x) :
    stepOnce C x = x :=
  dif_neg h

/-- `stepOnce` moves exactly when some site fires. -/
theorem stepOnce_ne_iff (C : OPHCarrier) (x : Records C) :
    stepOnce C x ≠ x ↔ ∃ i : Site C, localRepair C i x ≠ x := by
  constructor
  · intro hne
    by_contra hcon
    exact hne (stepOnce_of_normal C x hcon)
  · intro h
    rw [stepOnce_of_fire C x h]
    exact Classical.choose_spec h

/-- **Faithfulness, half one.** Every move `stepOnce` makes is a genuine accepted
asynchronous repair step of the repository's own relation. -/
theorem stepOnce_acceptedStep (C : OPHCarrier) (x : Records C) (hx : stepOnce C x ≠ x) :
    acceptedStep C x (stepOnce C x) := by
  have h : ∃ i : Site C, localRepair C i x ≠ x := (stepOnce_ne_iff C x).1 hx
  exact ⟨Classical.choose h, stepOnce_of_fire C x h, Classical.choose_spec h⟩

/-- **Faithfulness, half two.** The locked states of `stepOnce` are exactly the in-tree
`NormalForm C`, the states from which no accepted repair step applies. -/
theorem locked_stepOnce_iff_normalForm (C : OPHCarrier) (x : Records C) :
    Locked (stepOnce C) x ↔ NormalForm C x := by
  constructor
  · intro hlock y hstep
    obtain ⟨i, _, hfire⟩ := hstep
    exact (stepOnce_ne_iff C x).2 ⟨i, hfire⟩ hlock
  · intro hnf
    by_contra hne
    exact hnf _ (stepOnce_acceptedStep C x hne)

/-- Every `stepOnce` move strictly lowers the mismatch count: the in-tree descent lemma,
transported along `stepOnce`. -/
theorem mismatchCount_stepOnce_lt (C : OPHCarrier) {x : Records C}
    (hx : ¬ Locked (stepOnce C) x) :
    mismatchCount (stepOnce C x) < mismatchCount x := by
  have h : ∃ i : Site C, localRepair C i x ≠ x := (stepOnce_ne_iff C x).1 hx
  rw [stepOnce_of_fire C x h]
  exact mismatchCount_localRepair_lt C (Classical.choose h) x (Classical.choose_spec h)

/-- `stepOnce` terminates, through the in-tree descent lemma of the abstract layer. -/
theorem stepOnce_terminating (C : OPHCarrier) : Terminating (stepRel (stepOnce C)) :=
  descent_terminating (stepOnce C) mismatchCount fun x hne =>
    mismatchCount_stepOnce_lt C (x := x) hne

/-! ## Step 4: the declared observable bounds the relaxation time -/

/-- Each unlocked iterate consumes at least one unit of the mismatch count. -/
theorem mismatchCount_iterate_add_le (C : OPHCarrier) (x : Records C) :
    ∀ k : ℕ, (∀ j < k, ¬ Locked (stepOnce C) ((stepOnce C)^[j] x)) →
      mismatchCount ((stepOnce C)^[k] x) + k ≤ mismatchCount x := by
  intro k
  induction k with
  | zero => intro _; simp
  | succ k ih =>
      intro hpre
      have hk : ¬ Locked (stepOnce C) ((stepOnce C)^[k] x) := hpre k (Nat.lt_succ_self k)
      have hstep : mismatchCount (stepOnce C ((stepOnce C)^[k] x))
          < mismatchCount ((stepOnce C)^[k] x) := mismatchCount_stepOnce_lt C hk
      have hih := ih fun j hj => hpre j (by omega)
      rw [Function.iterate_succ_apply']
      omega

/-- **The declared observable bounds the rate.** The number of accepted repair steps to
first locking is at most `mismatchCount x`, and `mismatchCount` is a function of
`obsMap C x` alone (`mismatchCount_congr`). There is therefore no family of
observably-identical OPH repair systems with unboundedly different relaxation times,
which is exactly what the abstract stutter construction produces for its own reading. -/
theorem firstLock_le_mismatchCount (C : OPHCarrier) {n : ℕ} {x : Records C}
    (h : FirstLock (stepOnce C) n x) : n ≤ mismatchCount x := by
  have := mismatchCount_iterate_add_le C x n h.2
  omega

/-! ## Step 5: the declared observable *determines* the relaxation time -/

/-- `Classical.choose` picks the same witness from pointwise-equivalent predicates. This
restates `Primitives.choose_eq_of_pred_iff`, which is `private` to that file and so cannot
be referenced here; the proof is the same three lines and introduces no additional content. -/
private theorem choose_eq_of_pred_iff {α : Sort*} {p q : α → Prop}
    (hpq : ∀ a, p a ↔ q a) (hp : ∃ a, p a) (hq : ∃ a, q a) :
    Classical.choose hp = Classical.choose hq := by
  have hpq' : p = q := funext fun a => propext (hpq a)
  subst hpq'
  rfl

/-- Locking is a function of the declared observable data: whether any site fires is,
by the in-tree `localRepair_fire_congr`. -/
theorem locked_stepOnce_congr {x y : Records C} (h : obsMap C x = obsMap C y) :
    Locked (stepOnce C) x ↔ Locked (stepOnce C) y := by
  have key : (∃ i : Site C, localRepair C i x ≠ x) ↔ (∃ i : Site C, localRepair C i y ≠ y) :=
    exists_congr fun i => localRepair_fire_congr C h i
  constructor
  · intro hlx
    by_contra hny
    exact ((stepOnce_ne_iff C x).2 (key.2 ((stepOnce_ne_iff C y).1 hny))) hlx
  · intro hly
    by_contra hnx
    exact ((stepOnce_ne_iff C y).2 (key.1 ((stepOnce_ne_iff C x).1 hnx))) hly

/-- `stepOnce` is a gauge congruence: it fires the same site and installs the same repair
on records with equal declared observable data. This is `Primitives.obsMap_Repair_congr`'s
single-step engine, stated for `stepOnce`. -/
theorem obsMap_stepOnce_congr {x y : Records C} (h : obsMap C x = obsMap C y) :
    obsMap C (stepOnce C x) = obsMap C (stepOnce C y) := by
  by_cases hx : ∃ i : Site C, localRepair C i x ≠ x
  · have hy : ∃ i : Site C, localRepair C i y ≠ y := by
      obtain ⟨i, hi⟩ := hx
      exact ⟨i, (localRepair_fire_congr C h i).1 hi⟩
    have hsite : Classical.choose hx = Classical.choose hy :=
      choose_eq_of_pred_iff (fun i => localRepair_fire_congr C h i) hx hy
    rw [stepOnce_of_fire C x hx, stepOnce_of_fire C y hy, ← hsite]
    exact obsMap_localRepair_congr C h (Classical.choose hx)
  · have hy : ¬ ∃ i : Site C, localRepair C i y ≠ y := by
      intro hy'
      obtain ⟨i, hi⟩ := hy'
      exact hx ⟨i, (localRepair_fire_congr C h i).2 hi⟩
    rw [stepOnce_of_normal C x hx, stepOnce_of_normal C y hy]
    exact h

/-- Whole trajectories agree on the declared observable data. -/
theorem obsMap_iterate_congr {x y : Records C} (h : obsMap C x = obsMap C y) :
    ∀ n : ℕ, obsMap C ((stepOnce C)^[n] x) = obsMap C ((stepOnce C)^[n] y) := by
  intro n
  induction n with
  | zero => simpa using h
  | succ n ih =>
      rw [Function.iterate_succ_apply', Function.iterate_succ_apply']
      exact obsMap_stepOnce_congr ih

/-- **The exact negation of the abstract no-go, at the repository's own observable.** The
first-locking count is a function of the declared observable data: records exposing the
same overlap data lock after the same number of accepted repair steps under the `stepOnce`
schedule. A rate functional on
`obsMap` therefore exists for that schedule; the abstract `no_rate_functional` refutes one
for the abstract `DeclaredReading`, and that refutation does not transfer here. -/
theorem firstLock_obs_determined {x y : Records C} (h : obsMap C x = obsMap C y) {n : ℕ}
    (hx : FirstLock (stepOnce C) n x) : FirstLock (stepOnce C) n y := by
  refine ⟨(locked_stepOnce_congr (obsMap_iterate_congr h n)).1 hx.1, ?_⟩
  intro m hm hcon
  exact hx.2 m hm ((locked_stepOnce_congr (obsMap_iterate_congr h m)).2 hcon)

/-- The uniqueness companion: observably-identical records cannot have different
first-locking counts. Contrast `RateNonidentifiability.sameReading_differentFirstLock`,
which produces exactly that for the abstract reading. -/
theorem firstLock_eq_of_obs_eq {x y : Records C} (h : obsMap C x = obsMap C y) {m n : ℕ}
    (hx : FirstLock (stepOnce C) m x) (hy : FirstLock (stepOnce C) n y) : m = n :=
  firstLock_unique (firstLock_obs_determined h hx) hy

end OPH.RateBridgeObstruction

/- Axiom audit. -/

#print axioms OPH.RateBridgeObstruction.mismatchCount_congr
#print axioms OPH.RateBridgeObstruction.acceptedStep_changes_obs
#print axioms OPH.RateBridgeObstruction.no_observably_idle_acceptedStep
#print axioms OPH.RateBridgeObstruction.acceptedStepLR_changes_obs
#print axioms OPH.RateBridgeObstruction.stepOnce_acceptedStep
#print axioms OPH.RateBridgeObstruction.locked_stepOnce_iff_normalForm
#print axioms OPH.RateBridgeObstruction.firstLock_le_mismatchCount
#print axioms OPH.RateBridgeObstruction.obsMap_stepOnce_congr
#print axioms OPH.RateBridgeObstruction.firstLock_obs_determined
#print axioms OPH.RateBridgeObstruction.firstLock_eq_of_obs_eq
