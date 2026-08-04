import Mathlib
import ObserverPatchHolography.AbstractRewriting

/-!
# Rate non-identifiability for finite locking systems (abstract; **not** repository-facing)

This file proves an *abstract* non-identifiability result about deterministic maps on an
arbitrary type. It endorses no physical claim. A no-go is the result.

## SCOPE WARNING — read before quoting anything below

**This result does not apply to the repository's own repair dynamics, and the companion module
`ObserverPatchHolography.RateBridgeObstruction` proves that it cannot.**

An earlier framing of this file described it as machine-checking a *corpus* boundary — that
patch consistency data does not determine how fast a locking transition happens. That framing
was wrong and has been withdrawn. The theorems below are about `T : S → S` for an opaque `S`
and the generic rewriting relation `OPH.AbstractRewriting.stepRel`, which the abstract-rewriting
module itself calls a skeleton that "does not commit to OPH-specific structure". They are *not*
about `OPH.Records`, `OPH.obsMap`, or `OPH.acceptedStep`.

The bridge to those definitions was attempted and **fails at a located hypothesis**. The
construction below needs an accepted step that moves the state while leaving the declared
observable fixed — `stutter_has_projection_fixing_step` supplies exactly such steps. The
repository's `OPH.acceptedStep` has none: every accepted repair step strictly lowers
`OPH.mismatchCount`, which is a function of `OPH.obsMap` alone, so every accepted step strictly
changes the declared observable (`RateBridgeObstruction.acceptedStep_changes_obs`, and for every
local move obeying the declared laws H1–H3, `acceptedStepLR_changes_obs`).

The concrete layer in fact satisfies the **opposite** of the conclusion below: the declared
observable *bounds* the relaxation time (`RateBridgeObstruction.firstLock_le_mismatchCount`) and
indeed *determines* it (`RateBridgeObstruction.firstLock_obs_determined`). Read this file as a
theorem about abstract transition systems and about the weaker `DeclaredReading` defined here,
never as a statement about the OPH repair layer.

## Which sense of "speed"

The word "speed" names four separable quantities. **This file addresses exactly one of them.**

1. *Sharpness in a control parameter* — a critical exponent, a critical-window width. **Static.**
   **Not addressed here.** The systems below carry no control parameter, so no statement here
   can be read as one about a critical exponent or window.
2. *Relaxation time of the operator* — the number of `stepRel T` steps taken before the state
   first locks. **This is the quantity addressed here**, as `FirstLock`. It is a count of
   abstract steps, not of `OPH.acceptedStep` steps.
3. *Finite-size width of the ambiguous band*, and its scaling in a system size `L`.
   **Not addressed here.** No size parameter appears.
4. *Contraction factor of the repair map* — a per-step gap. **Not addressed here**, and for a
   sharper reason than the other two: the construction below does *not* preserve it. See
   `stutter_has_projection_fixing_step` and the conditionality discussion under it.

A result about (2) says nothing about (1). The four are not interchangeable and no theorem
in this file may be quoted across them.

## What is proved

For a finite deterministic operator `T`, locking is `IsFixedPt T`, so locked states are exactly
the normal forms of the in-tree **generic** rewriting relation `OPH.AbstractRewriting.stepRel T`.
`stepRel` is defined over an opaque type and is not the repository's accepted repair relation
`OPH.acceptedStep`; see the scope warning above. Given any positive stutter factor `s + 1`, the
stutter extension `stutterStep T s` on the finite carrier `S × Fin (s + 1)`:

* has exactly the same locked set, fibrewise, in both directions (`locked_stutter_iff`);
* has the same projected `stepRel`-reachability (`proj_reachable`, `lift_reachable`);
* has the same projected basin (`basin_stutter_iff`);
* multiplies the first-locking step count by `s + 1` (`firstLock_stutter`).

Hence no function of the declared reading returns the first-locking count
(`no_rate_functional`), while supplying the stutter factor separately makes it determinate
(`rate_determined_by_clock`). The missing datum is exactly a clock or transition law.

## The projection, stated

The projection used here is `Prod.fst`: forget the counter. It is *stipulated by this file*, not
taken from the repository. The corpus prose that an earlier framing cited in its defence
(`README.md` lines 292--295 and 299--300; `physics-problems/compact_record_transients.md` line
73, "a repair eigenvalue is not a physical duration without a clock map") says only that no
*physical clock* has been supplied. It does not say that an adjoined counter is invisible to the
repository's declared observable, and `RateBridgeObstruction` shows it is not: `OPH.obsMap`
strictly changes on every accepted repair step. Whatever `Prod.fst` is, it is not `OPH.obsMap`.

## Why stuttering is admissible for *this* file's reading

`extra/observable_normal_forms.tex` lines 633--665 lets a kernel sampling rewrite schedules whose
support consists of permitted rewrite steps *or stuttering moves*. That is a statement about a
stochastic kernel in the TeX corpus; it has no Lean counterpart in this repository, and it does
not license stuttering in `OPH.acceptedStep`, which provably admits none. The moves used below
are admissible for the `DeclaredReading` defined in this file and for nothing else in tree.

## Conditionality

The no-go is conditional on no per-step contraction factor being declared observable; see
`stutter_has_projection_fixing_step`. That condition is met by the corpus's own statement at
`extra/observable_normal_forms.tex` lines 1486--1495, that the finite resampling checker
supplies "neither a spectral gap nor a convergence rate".
-/

open Function Relation OPH.AbstractRewriting

namespace OPH.RateNonidentifiability

variable {S : Type*}

/-! ## Locking, basins, and first-locking time -/

/-- A state is locked when the repair operator no longer moves it. Locking is defined as
fixed-point-ness so that it coincides with normal-form-ness of the in-tree *generic* rewriting
relation `stepRel`; see `locked_iff_normalForm`. -/
def Locked (T : S → S) (x : S) : Prop := IsFixedPt T x

theorem locked_iff (T : S → S) (x : S) : Locked T x ↔ T x = x := Iff.rfl

/-- The bridge to the in-tree abstract rewriting layer: locked states are exactly the normal
forms of `OPH.AbstractRewriting.stepRel`. -/
theorem locked_iff_normalForm (T : S → S) (x : S) :
    Locked T x ↔ IsNormalForm (stepRel T) x := by
  constructor
  · intro h y hy
    exact hy.2 h
  · intro h
    exact Classical.not_not.1 fun hne => h (T x) ⟨rfl, hne⟩

/-- The basin: states that eventually lock. -/
def Basin (T : S → S) (x : S) : Prop := ∃ n, Locked T (T^[n] x)

/-- **Sense (2) of "speed".** `FirstLock T n x` says the state `x` first locks after exactly `n`
steps of the operator. It counts `stepRel T` steps, not `OPH.acceptedStep` steps. It is not a
critical exponent, not a window width, and not a per-step contraction factor. -/
def FirstLock (T : S → S) (n : ℕ) (x : S) : Prop :=
  Locked T (T^[n] x) ∧ ∀ m < n, ¬ Locked T (T^[m] x)

theorem firstLock_unique {T : S → S} {m n : ℕ} {x : S}
    (hm : FirstLock T m x) (hn : FirstLock T n x) : m = n := by
  by_contra hne
  rcases Nat.lt_or_ge m n with h | h
  · exact hn.2 m h hm.1
  · exact hm.2 n (by omega) hn.1

theorem basin_of_firstLock {T : S → S} {n : ℕ} {x : S} (h : FirstLock T n x) : Basin T x :=
  ⟨n, h.1⟩

theorem basin_of_iterate {α : Type*} {f : α → α} {x y : α} {m : ℕ}
    (h : f^[m] x = y) (hy : Basin f y) : Basin f x := by
  obtain ⟨n, hn⟩ := hy
  exact ⟨n + m, by rw [Function.iterate_add_apply, h]; exact hn⟩

/-- An iterated run whose every intermediate state actually moves is a run of the accepted
repair relation. -/
theorem reflTransGen_of_iterate {α : Type*} (f : α → α) (p : α) (n : ℕ)
    (h : ∀ i < n, f (f^[i] p) ≠ f^[i] p) :
    ReflTransGen (stepRel f) p (f^[n] p) := by
  induction n with
  | zero => exact ReflTransGen.refl
  | succ n ih =>
      refine (ih fun i hi => h i (by omega)).tail ?_
      exact ⟨by rw [Function.iterate_succ_apply'], h n (Nat.lt_succ_self n)⟩

/-! ## The stutter extension -/

section Stutter

variable [DecidableEq S]

/-- The stutter extension with factor `s + 1`, an arbitrary positive natural. The carrier
`S × Fin (s + 1)` is finite whenever `S` is. A locked state never stutters: the locked branch is
tested first, which is what keeps the fixed-point set exactly preserved. -/
def stutterStep (T : S → S) (s : ℕ) : S × Fin (s + 1) → S × Fin (s + 1) :=
  fun p =>
    if T p.1 = p.1 then p
    else if p.2.val = 0 then (T p.1, Fin.last s)
    else (p.1, ⟨p.2.val - 1, by have := p.2.isLt; omega⟩)

/-- The declared observable projection: forget the counter coordinate. -/
def proj (_T : S → S) (s : ℕ) : S × Fin (s + 1) → S := fun p => p.1

/-! The three computation rules for `stutterStep`, stated once so that no downstream proof has
to reason about the nesting of its branches. -/

/-- A locked state does not move. -/
theorem stutterStep_locked (T : S → S) (s : ℕ) (p : S × Fin (s + 1)) (h : Locked T p.1) :
    stutterStep T s p = p := by
  have h' : T p.1 = p.1 := h
  simp only [stutterStep]
  rw [if_pos h']

/-- With the counter exhausted, one base step fires and the counter reloads. -/
theorem stutterStep_reload (T : S → S) (s : ℕ) {x : S} (hx : ¬ Locked T x)
    (c : Fin (s + 1)) (hc : c.val = 0) :
    stutterStep T s (x, c) = (T x, Fin.last s) := by
  have hx' : ¬ (T x = x) := hx
  simp only [stutterStep]
  rw [if_neg hx', if_pos hc]

/-- Otherwise the extension takes a permitted stuttering move. -/
theorem stutterStep_idle (T : S → S) (s : ℕ) {x : S} (hx : ¬ Locked T x)
    (c : Fin (s + 1)) (hc : c.val ≠ 0) :
    stutterStep T s (x, c) = (x, ⟨c.val - 1, by have := c.isLt; omega⟩) := by
  have hx' : ¬ (T x = x) := hx
  simp only [stutterStep]
  rw [if_neg hx', if_neg hc]

/-- **Failure mode 1, defeated.** The full fixed-point set of the extended system is pinned in
both directions: an extended state is locked exactly when its base state is, so
`Fix (stutterStep T s) = Fix T ×ˢ univ`, the full preimage of `Fix T` under `Prod.fst`. Stated
precisely, because as a bare set claim "the fixed-point set is preserved" would be false — each
locked base state has `s + 1` distinct locked lifts. What is exact is the projected predicate. -/
theorem locked_stutter_iff (T : S → S) (s : ℕ) (p : S × Fin (s + 1)) :
    Locked (stutterStep T s) p ↔ Locked T p.1 := by
  obtain ⟨x, c⟩ := p
  constructor
  · intro h
    have h' : stutterStep T s (x, c) = (x, c) := h
    by_contra hne
    have hne' : ¬ Locked T x := hne
    by_cases hc : c.val = 0
    · rw [stutterStep_reload T s hne' c hc] at h'
      exact hne (congrArg Prod.fst h')
    · rw [stutterStep_idle T s hne' c hc] at h'
      have h2 : c.val - 1 = c.val := congrArg (fun q => (Prod.snd q).val) h'
      omega
  · exact stutterStep_locked T s (x, c)

/-- Every extended step projects either to nothing (a stuttering move) or to exactly one step of
the base repair operator. -/
theorem proj_step (T : S → S) (s : ℕ) (p : S × Fin (s + 1)) :
    (stutterStep T s p).1 = p.1 ∨ (stutterStep T s p).1 = T p.1 := by
  obtain ⟨x, c⟩ := p
  by_cases hl : Locked T x
  · left; rw [stutterStep_locked T s (x, c) hl]
  · by_cases hc : c.val = 0
    · right; rw [stutterStep_reload T s hl c hc]
    · left; rw [stutterStep_idle T s hl c hc]

/-- While the base state is unlocked the extension idles, decrementing the counter. -/
theorem stutter_iterate_hold (T : S → S) (s : ℕ) {x : S} (hx : ¬ Locked T x)
    (c : Fin (s + 1)) :
    ∀ i, i ≤ c.val →
      (stutterStep T s)^[i] (x, c) = (x, ⟨c.val - i, by have := c.isLt; omega⟩) := by
  have hx' : ¬ (T x = x) := hx
  intro i
  induction i with
  | zero => intro _; simp
  | succ i ih =>
      intro hi
      have hne : (⟨c.val - i, by have := c.isLt; omega⟩ : Fin (s + 1)).val ≠ 0 := by
        show c.val - i ≠ 0
        omega
      rw [Function.iterate_succ_apply', ih (by omega), stutterStep_idle T s hx _ hne]
      congr 1

/-- After exactly `c + 1` extended steps the base state has advanced by exactly one repair step
and the counter has been reloaded. -/
theorem stutter_step_next (T : S → S) (s : ℕ) {x : S} (hx : ¬ Locked T x)
    (c : Fin (s + 1)) :
    (stutterStep T s)^[c.val + 1] (x, c) = (T x, Fin.last s) := by
  rw [Function.iterate_succ_apply', stutter_iterate_hold T s hx c c.val le_rfl]
  exact stutterStep_reload T s hx _ (show c.val - c.val = 0 from Nat.sub_self c.val)

/-- One block: `s + 1` extended steps realise one base repair step. -/
theorem stutter_block (T : S → S) (s : ℕ) {x : S} (hx : ¬ Locked T x) :
    (stutterStep T s)^[s + 1] (x, Fin.last s) = (T x, Fin.last s) :=
  stutter_step_next T s hx (Fin.last s)

/-- `k` base repair steps cost exactly `(s + 1) * k` extended steps. -/
theorem stutter_iterate_mul (T : S → S) (s : ℕ) (x : S) :
    ∀ k, (∀ j < k, ¬ Locked T (T^[j] x)) →
      (stutterStep T s)^[(s + 1) * k] (x, Fin.last s) = (T^[k] x, Fin.last s) := by
  intro k
  induction k with
  | zero => intro _; simp
  | succ k ih =>
      intro hpre
      have hk : ¬ Locked T (T^[k] x) := hpre k (Nat.lt_succ_self k)
      have hsplit : (s + 1) * (k + 1) = (s + 1) + (s + 1) * k := by ring
      rw [hsplit, Function.iterate_add_apply, ih fun j hj => hpre j (by omega),
        stutter_block T s hk, Function.iterate_succ_apply']

/-- No extended state locks before the whole block schedule has run. -/
theorem stutter_not_locked_before (T : S → S) (s n : ℕ) (x : S)
    (hpre : ∀ j < n, ¬ Locked T (T^[j] x)) :
    ∀ m < (s + 1) * n,
      ¬ Locked (stutterStep T s) ((stutterStep T s)^[m] (x, Fin.last s)) := by
  intro m hm
  obtain ⟨q, r, hrs, hdm⟩ : ∃ q r, r < s + 1 ∧ m = (s + 1) * q + r :=
    ⟨m / (s + 1), m % (s + 1), Nat.mod_lt _ (Nat.succ_pos s),
      (Nat.div_add_mod m (s + 1)).symm⟩
  have hqn : q < n := by
    by_contra hcon
    have hcon' : n ≤ q := Nat.le_of_not_lt hcon
    have hmul : (s + 1) * n ≤ (s + 1) * q := Nat.mul_le_mul (Nat.le_refl (s + 1)) hcon'
    omega
  have hunlocked : ¬ Locked T (T^[q] x) := hpre _ hqn
  have hidx : m = r + (s + 1) * q := by rw [hdm]; ring
  have hrle : r ≤ (Fin.last s).val := by
    simp only [Fin.val_last]
    omega
  rw [hidx, Function.iterate_add_apply,
    stutter_iterate_mul T s x q fun j hj => hpre j (by omega),
    stutter_iterate_hold T s hunlocked (Fin.last s) r hrle]
  intro hcon
  exact hunlocked ((locked_stutter_iff T s _).1 hcon)

/-- **The rate multiplication.** First locking in `n` base repair steps becomes first locking in
`(s + 1) * n` steps of the extension, for an arbitrary positive stutter factor `s + 1`. -/
theorem firstLock_stutter (T : S → S) (s n : ℕ) (x : S) (h : FirstLock T n x) :
    FirstLock (stutterStep T s) ((s + 1) * n) (x, Fin.last s) := by
  refine ⟨?_, stutter_not_locked_before T s n x h.2⟩
  rw [stutter_iterate_mul T s x n h.2]
  exact (locked_stutter_iff T s _).2 h.1

/-! ## Preservation of the projected reading -/

/-- Every extended run projects to a run of the base `stepRel`. -/
theorem proj_reachable (T : S → S) (s : ℕ) {p q : S × Fin (s + 1)}
    (h : ReflTransGen (stepRel (stutterStep T s)) p q) :
    ReflTransGen (stepRel T) p.1 q.1 := by
  induction h with
  | refl => exact ReflTransGen.refl
  | @tail b c _hpb hbc ih =>
      have hb : ¬ Locked T b.1 := fun hlock =>
        hbc.2 ((locked_stutter_iff T s b).2 hlock)
      have hc : c.1 = (stutterStep T s b).1 := by rw [hbc.1]
      rcases proj_step T s b with h1 | h1
      · rw [hc, h1]
        exact ih
      · exact ih.tail ⟨by rw [hc, h1], hb⟩

/-- Every base run is realised by an extended run. -/
theorem lift_reachable (T : S → S) (s : ℕ) {x y : S}
    (h : ReflTransGen (stepRel T) x y) :
    ReflTransGen (stepRel (stutterStep T s)) (x, Fin.last s) (y, Fin.last s) := by
  induction h with
  | refl => exact ReflTransGen.refl
  | @tail b c _hxb hbc ih =>
      refine ih.trans ?_
      have hb : ¬ Locked T b := hbc.2
      have hmove : ∀ i < s + 1,
          stutterStep T s ((stutterStep T s)^[i] (b, Fin.last s))
            ≠ (stutterStep T s)^[i] (b, Fin.last s) := by
        intro i hi
        have hile : i ≤ (Fin.last s).val := by
          simp only [Fin.val_last]
          omega
        rw [stutter_iterate_hold T s hb (Fin.last s) i hile]
        intro hcon
        exact hb ((locked_stutter_iff T s _).1 hcon)
      have hpath := reflTransGen_of_iterate (stutterStep T s) (b, Fin.last s) (s + 1) hmove
      rw [stutter_block T s hb] at hpath
      rw [hbc.1]
      exact hpath

theorem proj_iterate_exists (T : S → S) (s : ℕ) (p : S × Fin (s + 1)) (n : ℕ) :
    ∃ k, ((stutterStep T s)^[n] p).1 = T^[k] p.1 := by
  induction n with
  | zero => exact ⟨0, rfl⟩
  | succ n ih =>
      obtain ⟨k, hk⟩ := ih
      rw [Function.iterate_succ_apply']
      rcases proj_step T s ((stutterStep T s)^[n] p) with h | h
      · exact ⟨k, by rw [h, hk]⟩
      · exact ⟨k + 1, by rw [h, hk, Function.iterate_succ_apply']⟩

theorem basin_stutter_of_base (T : S → S) (s : ℕ) :
    ∀ (n : ℕ) (p : S × Fin (s + 1)), Locked T (T^[n] p.1) → Basin (stutterStep T s) p := by
  intro n
  induction n with
  | zero => intro p hp; exact ⟨0, (locked_stutter_iff T s p).2 hp⟩
  | succ n ih =>
      intro p hp
      by_cases hlk : Locked T p.1
      · exact ⟨0, (locked_stutter_iff T s p).2 hlk⟩
      · have hstep : (stutterStep T s)^[p.2.val + 1] p = (T p.1, Fin.last s) :=
          stutter_step_next T s hlk p.2
        have hnext : Locked T (T^[n] (T p.1)) := by
          rw [← Function.iterate_succ_apply]
          exact hp
        exact basin_of_iterate hstep (ih (T p.1, Fin.last s) hnext)

/-- The basin is preserved exactly, fibrewise and in both directions. -/
theorem basin_stutter_iff (T : S → S) (s : ℕ) (p : S × Fin (s + 1)) :
    Basin (stutterStep T s) p ↔ Basin T p.1 := by
  constructor
  · rintro ⟨n, hn⟩
    obtain ⟨k, hk⟩ := proj_iterate_exists T s p n
    refine ⟨k, ?_⟩
    rw [← hk]
    exact (locked_stutter_iff T s _).1 hn
  · rintro ⟨n, hn⟩
    exact basin_stutter_of_base T s n p hn

/-! ## The declared reading and the no-go -/

/-- The three pieces of data this file declares readable: the locked set, `stepRel`-reachability,
and the basin. This is *this file's* stipulated reading, not the repository's observable
`OPH.obsMap`, which is strictly finer and — unlike this reading — determines the rate
(`RateBridgeObstruction.firstLock_obs_determined`). Nothing here is a clock and nothing here is
a per-step contraction factor. -/
structure DeclaredReading (S : Type*) where
  LockedSet : S → Prop
  Reach : S → S → Prop
  BasinSet : S → Prop

/-- The reading of the unextended system. -/
def reading (T : S → S) : DeclaredReading S where
  LockedSet := Locked T
  Reach := ReflTransGen (stepRel T)
  BasinSet := Basin T

/-- **Failure mode 2, defeated by naming the projection.** The reading of the extended system as
seen through `proj`, which forgets the counter. Everything an observer without a clock can read
off the stuttered system. -/
def projReading (T : S → S) (s : ℕ) : DeclaredReading S where
  LockedSet := fun x => ∃ c : Fin (s + 1), Locked (stutterStep T s) (x, c)
  Reach := fun x y => ∃ c d : Fin (s + 1),
    ReflTransGen (stepRel (stutterStep T s)) (x, c) (y, d)
  BasinSet := fun x => ∃ c : Fin (s + 1), Basin (stutterStep T s) (x, c)

/-- Every stutter extension reads identically to the unextended system, for every factor. -/
theorem projReading_eq (T : S → S) (s : ℕ) : projReading T s = reading T := by
  have hL : (fun x => ∃ c : Fin (s + 1), Locked (stutterStep T s) (x, c)) = Locked T := by
    funext x
    apply propext
    constructor
    · rintro ⟨c, hc⟩
      exact (locked_stutter_iff T s (x, c)).1 hc
    · intro h
      exact ⟨Fin.last s, (locked_stutter_iff T s (x, Fin.last s)).2 h⟩
  have hR : (fun x y => ∃ c d : Fin (s + 1),
      ReflTransGen (stepRel (stutterStep T s)) (x, c) (y, d)) = ReflTransGen (stepRel T) := by
    funext x y
    apply propext
    constructor
    · rintro ⟨c, d, h⟩
      exact proj_reachable T s h
    · intro h
      exact ⟨Fin.last s, Fin.last s, lift_reachable T s h⟩
  have hB : (fun x => ∃ c : Fin (s + 1), Basin (stutterStep T s) (x, c)) = Basin T := by
    funext x
    apply propext
    constructor
    · rintro ⟨c, hc⟩
      exact (basin_stutter_iff T s (x, c)).1 hc
    · intro h
      exact ⟨Fin.last s, (basin_stutter_iff T s (x, Fin.last s)).2 h⟩
  unfold projReading reading
  simp only [DeclaredReading.mk.injEq]
  exact ⟨hL, hR, hB⟩

/-- Two stutter factors, one reading, two different first-locking counts. -/
theorem sameReading_differentFirstLock (T : S → S) (x : S) (n : ℕ)
    (hlock : FirstLock T n x) (hpos : 0 < n) (s₁ s₂ : ℕ) (hne : s₁ ≠ s₂) :
    projReading T s₁ = projReading T s₂ ∧
      FirstLock (stutterStep T s₁) ((s₁ + 1) * n) (x, Fin.last s₁) ∧
      FirstLock (stutterStep T s₂) ((s₂ + 1) * n) (x, Fin.last s₂) ∧
      (s₁ + 1) * n ≠ (s₂ + 1) * n := by
  refine ⟨by rw [projReading_eq, projReading_eq], firstLock_stutter T s₁ n x hlock,
    firstLock_stutter T s₂ n x hlock, ?_⟩
  intro hcon
  have hsucc : s₁ + 1 = s₂ + 1 := Nat.eq_of_mul_eq_mul_right hpos hcon
  exact hne (by omega)

/-- **The no-go, sense (2), for `DeclaredReading` only.** No function from the reading declared
in this file — locked set, `stepRel`-reachability, basin — returns the first-locking step count.
The reading is constant across the stutter family while the count is not.

**This is not a statement about the OPH repair layer.** For the repository's own observable
`OPH.obsMap` the corresponding functional *does* exist: see
`RateBridgeObstruction.firstLock_obs_determined`, and `firstLock_le_mismatchCount` for the
bound. The gap between the two is that `DeclaredReading` is a triple of relations on the state
space while `obsMap` is a per-state observable that separates strictly more records.

In plain English: locked-set, reachability and basin data alone forbid nothing about how many
abstract steps locking takes. This is a statement about identifiability under one stipulated
reading, not about physics and not about OPH. -/
theorem no_rate_functional (T : S → S) (x : S) (n : ℕ)
    (hlock : FirstLock T n x) (hpos : 0 < n)
    (rate : DeclaredReading S → ℕ)
    (hrate : ∀ s : ℕ, FirstLock (stutterStep T s) (rate (projReading T s)) (x, Fin.last s)) :
    False := by
  have e0 : rate (projReading T 0) = 1 * n :=
    firstLock_unique (hrate 0) (firstLock_stutter T 0 n x hlock)
  have e1 : rate (projReading T 1) = 2 * n :=
    firstLock_unique (hrate 1) (firstLock_stutter T 1 n x hlock)
  rw [projReading_eq] at e0 e1
  omega

/-- **A uniqueness corollary, not a construction.** Given that `m` already *is* a first-locking
count for the extension, it equals `(s + 1) * n`. This is `firstLock_unique` against
`firstLock_stutter` and nothing more: it does not take a clock map as input, does not build a
rate functional, and does not show that the stutter factor plus the declared reading suffices to
recover `n`. An earlier gloss called this "the positive half" proving clock-sufficiency; that
was too strong and is withdrawn. -/
theorem rate_determined_by_clock (T : S → S) (x : S) (n s m : ℕ)
    (hlock : FirstLock T n x)
    (hm : FirstLock (stutterStep T s) m (x, Fin.last s)) :
    m = (s + 1) * n :=
  firstLock_unique hm (firstLock_stutter T s n x hlock)

/-! ## What the construction does not preserve, and hence what the no-go is conditional on -/

/-- **Conditionality of the no-go.** The extension has reachable steps that leave the projection
fixed while the base state is unlocked. A declared per-step contraction factor or spectral gap
would therefore distinguish the stutter family, and the no-go above would not apply to it.

This is why nothing in this file may be read as a statement about sense (4). It is also why the
no-go is conditional on no per-step gap being declared.

**This lemma is the precise point at which the bridge to OPH fails.** It exhibits the
observably-idle step the construction depends on. `RateBridgeObstruction.acceptedStep_changes_obs`
proves the repository's `OPH.acceptedStep` admits no such step, so no stutter extension of the
OPH repair dynamics exists. -/
theorem stutter_has_projection_fixing_step (T : S → S) (s : ℕ) (x : S)
    (hx : ¬ Locked T x) (hs : 0 < s) :
    ∃ p : S × Fin (s + 1),
      ¬ Locked T p.1 ∧ ¬ Locked (stutterStep T s) p ∧ (stutterStep T s p).1 = p.1 := by
  refine ⟨(x, Fin.last s), hx, fun hcon => hx ((locked_stutter_iff T s _).1 hcon), ?_⟩
  have hc : (Fin.last s).val ≠ 0 := by
    show s ≠ 0
    omega
  rw [stutterStep_idle T s hx (Fin.last s) hc]

omit [DecidableEq S] in
/-- By contrast the base system has no such step: an unlocked base state always moves. -/
theorem base_has_no_projection_fixing_step (T : S → S) (x : S) (hx : ¬ Locked T x) :
    T x ≠ x := hx

/-- The stutter extension preserves termination, so it does not manufacture a basin. Stated
through the in-tree descent lemma. -/
theorem stutter_terminating_of_descent (T : S → S) (s : ℕ) (Φ : S → ℕ)
    (hdesc : ∀ x, T x ≠ x → Φ (T x) < Φ x) :
    Terminating (stepRel (stutterStep T s)) := by
  refine descent_terminating (stutterStep T s) (fun p => (s + 1) * Φ p.1 + p.2.val) ?_
  intro p hne
  obtain ⟨x, c⟩ := p
  have hbase : ¬ Locked T x := fun hlock => hne ((locked_stutter_iff T s (x, c)).2 hlock)
  have hdrop : Φ (T x) + 1 ≤ Φ x := hdesc x hbase
  have hmul : (s + 1) * (Φ (T x) + 1) ≤ (s + 1) * Φ x :=
    Nat.mul_le_mul (Nat.le_refl _) hdrop
  rw [Nat.mul_add, Nat.mul_one] at hmul
  by_cases hc : c.val = 0
  · rw [stutterStep_reload T s hbase c hc]
    show (s + 1) * Φ (T x) + s < (s + 1) * Φ x + c.val
    omega
  · rw [stutterStep_idle T s hbase c hc]
    show (s + 1) * Φ x + (c.val - 1) < (s + 1) * Φ x + c.val
    omega

end Stutter

/-! ## An explicit finite witness

Four states, a descent repair operator, one locked state, three unlocked states with distinct
first-locking counts. Stutter factor seven gives a twenty-eight-state extension with an
identical declared reading and a first-locking count of twenty-one instead of three. -/

section Witness

/-- The four-state descent witness: mismatch level three decays to the locked level zero. -/
def wT : Fin 4 → Fin 4 := fun i => if i = 0 then 0 else i - 1

instance decidableLocked {S : Type*} [DecidableEq S] (T : S → S) (x : S) :
    Decidable (Locked T x) :=
  inferInstanceAs (Decidable (T x = x))

instance decidableFirstLock {S : Type*} [DecidableEq S] (T : S → S) (n : ℕ) (x : S) :
    Decidable (FirstLock T n x) :=
  inferInstanceAs (Decidable (Locked T (T^[n] x) ∧ ∀ m < n, ¬ Locked T (T^[m] x)))

/-- The locked set of the witness is exactly `{0}`. -/
theorem wT_locked_iff : ∀ i : Fin 4, Locked wT i ↔ i = 0 := by decide

/-- Three unlocked states with three distinct first-locking counts. -/
theorem wT_firstLock_zero : FirstLock wT 0 0 := by decide

theorem wT_firstLock_one : FirstLock wT 1 1 := by decide

theorem wT_firstLock_two : FirstLock wT 2 2 := by decide

theorem wT_firstLock_three : FirstLock wT 3 3 := by decide

/-- The witness terminates, through the in-tree descent lemma. -/
theorem wT_terminating : Terminating (stepRel wT) :=
  descent_terminating wT (fun i => i.val) (by decide)

/-- Every stutter extension of the witness terminates too. -/
theorem wT_stutter_terminating (s : ℕ) : Terminating (stepRel (stutterStep wT s)) :=
  stutter_terminating_of_descent wT s (fun i => i.val) (by decide)

/-- Stutter factor seven multiplies the first-locking count from three to twenty-one. -/
theorem wT_stutter_firstLock : FirstLock (stutterStep wT 6) 21 (3, Fin.last 6) := by
  have h := firstLock_stutter wT 6 3 3 wT_firstLock_three
  have harith : (6 + 1) * 3 = 21 := by norm_num
  rwa [harith] at h

/-- Every stutter extension of the witness has the same declared reading. -/
theorem wT_readings_agree (s₁ s₂ : ℕ) : projReading wT s₁ = projReading wT s₂ := by
  rw [projReading_eq, projReading_eq]

/-- The witness instance of the no-go. -/
theorem wT_no_rate_functional (rate : DeclaredReading (Fin 4) → ℕ)
    (hrate : ∀ s : ℕ,
      FirstLock (stutterStep wT s) (rate (projReading wT s)) (3, Fin.last s)) :
    False :=
  no_rate_functional wT 3 3 wT_firstLock_three (by norm_num) rate hrate

end Witness

end OPH.RateNonidentifiability

/- Axiom audit. -/

#print axioms OPH.RateNonidentifiability.locked_iff_normalForm
#print axioms OPH.RateNonidentifiability.locked_stutter_iff
#print axioms OPH.RateNonidentifiability.proj_reachable
#print axioms OPH.RateNonidentifiability.lift_reachable
#print axioms OPH.RateNonidentifiability.basin_stutter_iff
#print axioms OPH.RateNonidentifiability.firstLock_stutter
#print axioms OPH.RateNonidentifiability.projReading_eq
#print axioms OPH.RateNonidentifiability.sameReading_differentFirstLock
#print axioms OPH.RateNonidentifiability.no_rate_functional
#print axioms OPH.RateNonidentifiability.rate_determined_by_clock
#print axioms OPH.RateNonidentifiability.stutter_has_projection_fixing_step
#print axioms OPH.RateNonidentifiability.stutter_terminating_of_descent
#print axioms OPH.RateNonidentifiability.wT_firstLock_three
#print axioms OPH.RateNonidentifiability.wT_stutter_firstLock
#print axioms OPH.RateNonidentifiability.wT_no_rate_functional
