import ObserverPatchHolography.YangMillsGap

/-!
# Yang–Mills finite repair-gap — a concrete NON-VACUITY witness for Theorem 7.3

`thm_7_3_finite_gap` (in `ObserverPatchHolography.YangMillsGap`) is a **conditional**:
IF a finite nonempty collar family consists of mutually commuting star projections
whose non-commutative product is `P₀`, and the relaxation rates are strictly positive,
THEN the repair generator dominates `c_* · (I − P₀)` with `c_* > 0` (the finite-stage
representation gap `Δ_rep ≥ c_* > 0`).

A conditional is empty unless its hypothesis bundle is jointly satisfiable. Nothing in
the repo currently rules out that the premises (commuting star projections whose
`noncommProd` equals `P₀`, over a nonempty index set, with positive rates) are jointly
unsatisfiable — which would make the flagship theorem vacuous. This file discharges
that: it exhibits a concrete complete real inner-product space
`W := EuclideanSpace ℝ (Fin 1)` and a concrete two-collar family satisfying EVERY
premise, so `thm_7_3_finite_gap` genuinely **fires** — certifying it is not vacuous.
This mirrors the non-vacuity methodology already used for the reconstruction layer in
`Primitives.lean` (`demoCarrier_terminates`, `demoCarrier_dir_confluent`): pair the
abstract theorem with a machine-checked concrete model.

## HONEST SCOPE

This is a *minimal* witness.

* The two collars use the degenerate collar projection `E_C = 0` (constants space
  `{0}`), so the joint fixed projection is `P₀ = 0` (`wP0_eq_zero`) and the gap operator
  is the full identity `I − P₀ = I` — a genuine **nonzero** positive operator. (Contrast
  the `E_C = 1` choice, which would give the empty gap `I − P₀ = 0`: a vacuous-feeling
  "positive gap on a zero-dimensional complement".)
* The two rates are **distinct** (`1` and `2`), so `c_* = min = 1` and the resulting
  operator bound `1 · (I − P₀) ≤ L_r^rep` is non-reflexive (the generator is `3 · I`).

What it does **not** do, and does not claim: exhibit a *proper* nonzero collar
projection (`0 ≠ E_C ≠ 1`) — a natural strengthening left as future work — nor touch
Assumption 9.2 (the continuum certificate `Δ_YM = Δ_rep`), which remains the open
frontier. Adds no `sorry` and no `axiom`; the `#print axioms` lines below are expected
to report only `[propext, Classical.choice, Quot.sound]`.
-/

namespace ObserverPatchHolography.YangMillsGapWitness

open ObserverPatchHolography.YangMillsGap

/-- A concrete complete real inner-product space: 1-dimensional Euclidean space. -/
abbrev W : Type := EuclideanSpace ℝ (Fin 1)

/-- Two collars, indexed by `Bool`. -/
def wS : Finset Bool := {false, true}

/-- Each collar acts by the (degenerate) zero projection — constants space `{0}`. -/
noncomputable def wEc : Bool → (W →L[ℝ] W) := fun _ => 0

/-- Distinct strictly-positive relaxation rates, so `c_* = min = 1 < 2` and the
    resulting operator bound is non-reflexive. -/
def wRate : Bool → ℝ := fun b => if b then 2 else 1

theorem whne : wS.Nonempty := ⟨false, by unfold wS; decide⟩

theorem whE : ∀ a ∈ wS, IsStarProjection (wEc a) :=
  fun _ _ => IsStarProjection.zero _

theorem whc : (↑wS : Set Bool).Pairwise (Function.onFun Commute wEc) := by
  intro a _ b _ _
  exact Commute.refl (0 : W →L[ℝ] W)

/-- The joint fixed (constants) projection `P₀` of the collar family. -/
noncomputable def wP0 : W →L[ℝ] W := wS.noncommProd wEc whc

/-- With every collar the zero projection, the joint fixed projection is `0`, so the
    gap operator `I − P₀` is the full identity `I` — a genuine nonzero positive operator. -/
theorem wP0_eq_zero : wP0 = 0 := by
  simp [wP0, wS, wEc]

theorem whrate : ∀ a ∈ wS, 0 < wRate a := by
  intro b _
  cases b <;> norm_num [wRate]

/-- **NON-VACUITY WITNESS for Theorem 7.3 (finite representation gap).** The flagship
    finite-repair-gap theorem `thm_7_3_finite_gap` fires on a concrete complete real
    inner-product space: there is a strictly positive `c_*` with
    `c_* · (I − P₀) ≤ L_r^rep`. The premise bundle (nonempty collar family, commuting
    star projections, `noncommProd = P₀`, positive rates) is therefore jointly
    satisfiable — the flagship conditional is not vacuous. By `wP0_eq_zero` the gap
    operator `I − P₀ = I` here, and `c_* = min {1, 2} = 1 > 0`. -/
theorem thm_7_3_finite_gap_nonvacuous :
    ∃ cstar : ℝ, 0 < cstar ∧
      cstar • ((1 : W →L[ℝ] W) - wP0) ≤ repairGenerator wS wEc wRate :=
  thm_7_3_finite_gap wS whne wEc wP0 whE whc rfl wRate whrate

/-- The repair generator of this witness evaluates to `3 · I` (rates `1` and `2`, each
    collar `I − E_C = I − 0 = I`). Together with `wP0_eq_zero` (`I − P₀ = I`) and
    `c_* = min {1, 2} = 1`, this makes the witnessed bound the genuinely **non-reflexive**
    `1 · I ≤ 3 · I` (gap slack `2 · I ⪰ 0`), not a potemkin `I ≤ I`. -/
theorem repairGenerator_eq :
    repairGenerator wS wEc wRate = (3 : ℝ) • (1 : W →L[ℝ] W) := by
  unfold repairGenerator
  simp only [wS, wEc, sub_zero]
  rw [show ({false, true} : Finset Bool) = insert false {true} from rfl,
      Finset.sum_insert (by decide), Finset.sum_singleton, ← add_smul]
  norm_num [wRate]

/-! ## Axiom self-audit (build-log visible) -/

#print axioms thm_7_3_finite_gap_nonvacuous
#print axioms wP0_eq_zero
#print axioms repairGenerator_eq

end ObserverPatchHolography.YangMillsGapWitness
