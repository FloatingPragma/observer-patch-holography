import ObserverPatchHolography.YangMillsGap

/-!
# Yang–Mills finite repair-gap — witnesses and a noncommuting countermodel

`thm_7_3_finite_gap` (in `ObserverPatchHolography.YangMillsGap`) is a
legacy/special **commuting-projection conditional**:
IF a finite nonempty collar family consists of mutually commuting star projections
whose non-commutative product is `P₀`, and the relaxation rates are strictly positive,
THEN the repair generator dominates `c_* · (I − P₀)` with `c_* > 0` (the finite-stage
representation gap `Δ_rep ≥ c_* > 0`).

A conditional is empty unless its hypothesis bundle is jointly satisfiable. Nothing in
the repo currently rules out that the premises (commuting star projections whose
`noncommProd` equals `P₀`, over a nonempty index set, with positive rates) are jointly
unsatisfiable — which would make that special conditional vacuous. This file discharges
that: it exhibits a concrete complete real inner-product space
`W := EuclideanSpace ℝ (Fin 1)` and a concrete two-collar family satisfying EVERY
premise, so `thm_7_3_finite_gap` genuinely **fires** — certifying the
commuting branch is not vacuous.
This mirrors the non-vacuity methodology already used for the reconstruction layer in
`Primitives.lean` (`demoCarrier_terminates`, `demoCarrier_dir_confluent`): pair the
abstract theorem with a machine-checked concrete model.

## HONEST SCOPE

This is a *minimal* witness.

* The two collars use the degenerate collar projection `E_C = 0` (fixed space
  `{0}`), so the joint fixed projection is `P₀ = 0` (`wP0_eq_zero`) and the gap operator
  is the full identity `I − P₀ = I` — a genuine **nonzero** positive operator. (Contrast
  the `E_C = 1` choice, which would give the empty gap `I − P₀ = 0`: a vacuous-feeling
  "positive gap on a zero-dimensional complement".)
* The two rates are **distinct** (`1` and `2`), so `c_* = min = 1` and the resulting
  operator bound `1 · (I − P₀) ≤ L_r^rep` is non-reflexive (the generator is `3 · I`).

The first witness also does **not** instantiate `lemma_7_2` or identify its
zero projection with expectation onto the constants of a hidden fiber.  On
`Fin 1`, the space of constant functions is all of `W`, not `{0}`.  The
witness inhabits only the stripped operator premises that actually occur in
`thm_7_3_finite_gap`; connecting a uniform-fiber relaxation to those collar
projections and rates is a separate, open bridge.  It does **not** touch the
current paper's noncommuting Dobrushin branch or any continuum certificate.
The second section gives a
separate exact three-state rational countermodel: two proper symmetric
idempotent Markov projections fail to commute, and their rate-one repair
generator has explicit eigenvalues `1/2` and `3/2`, rather than the commuting
subset sums `0`, `1`, and `2`.  It is a regression guard against importing the
commuting subset-sum conclusion into the noncommuting setting; it is not a
Dobrushin-gap proof or a continuum model.

Adds no `sorry` and no `axiom`; the `#print axioms` lines below are expected to
report only standard Mathlib axioms where applicable.
-/

namespace ObserverPatchHolography.YangMillsGapWitness

open ObserverPatchHolography.YangMillsGap

/-- A concrete complete real inner-product space: 1-dimensional Euclidean space. -/
abbrev W : Type := EuclideanSpace ℝ (Fin 1)

/-- Two collars, indexed by `Bool`. -/
def wS : Finset Bool := {false, true}

/-- Each collar acts by the (degenerate) zero projection, whose fixed space is `{0}`. -/
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

/-- The joint fixed projection `P₀` of this collar family. -/
noncomputable def wP0 : W →L[ℝ] W := wS.noncommProd wEc whc

/-- With every collar the zero projection, the joint fixed projection is `0`, so the
    gap operator `I − P₀` is the full identity `I` — a genuine nonzero positive operator. -/
theorem wP0_eq_zero : wP0 = 0 := by
  simp [wP0, wS, wEc]

theorem whrate : ∀ a ∈ wS, 0 < wRate a := by
  intro b _
  cases b <;> norm_num [wRate]

/-- **NON-VACUITY WITNESS for the legacy commuting Theorem 7.3.** The
    finite-repair-gap theorem `thm_7_3_finite_gap` fires on a concrete complete real
    inner-product space: there is a strictly positive `c_*` with
    `c_* · (I − P₀) ≤ L_r^rep`. The premise bundle (nonempty collar family, commuting
    star projections, `noncommProd = P₀`, positive rates) is therefore jointly
    satisfiable — the special commuting conditional is not vacuous. By `wP0_eq_zero` the gap
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

/-! ## Exact noncommuting three-state countermodel

`ncE₁` and `ncE₂` are the conditional-expectation matrices for the
partitions `{0} | {1,2}` and `{1} | {0,2}` of a uniform three-state
space.  The receipts below check, over `ℚ`, that each is a proper
symmetric idempotent Markov projection and that they do not commute.

With unit rates, `ncL = (I - ncE₁) + (I - ncE₂)`.  Its displayed
nonconstant eigenvectors have eigenvalues `1/2` and `3/2`.  Since a
commuting two-projection subset-sum argument at unit rates would allow
only `0`, `1`, and `2`, this is an exact finite counterexample to using
that spectral conclusion without commutation.
-/

/-- The conditional expectation onto functions constant on `{1,2}`. -/
def ncE₁ : Matrix (Fin 3) (Fin 3) ℚ :=
  !![1, 0, 0;
     0, 1 / 2, 1 / 2;
     0, 1 / 2, 1 / 2]

/-- The conditional expectation onto functions constant on `{0,2}`. -/
def ncE₂ : Matrix (Fin 3) (Fin 3) ℚ :=
  !![1 / 2, 0, 1 / 2;
     0, 1, 0;
     1 / 2, 0, 1 / 2]

/-- The constant vector on the three-state space. -/
def ncConstant : Fin 3 → ℚ := ![1, 1, 1]

/-- A nonconstant low-mode vector. -/
def ncLow : Fin 3 → ℚ := ![1, -1, 0]

/-- A nonconstant high-mode vector. -/
def ncHigh : Fin 3 → ℚ := ![1, 1, -2]

/-- The rate-one noncommuting repair generator. -/
def ncL : Matrix (Fin 3) (Fin 3) ℚ :=
  !![1 / 2, 0, -1 / 2;
     0, 1 / 2, -1 / 2;
     -1 / 2, -1 / 2, 1]

/-- The displayed matrix is exactly `(I - ncE₁) + (I - ncE₂)` with
unit rates. -/
theorem ncL_eq_rate_one_repairGenerator :
    ncL = ((1 : Matrix (Fin 3) (Fin 3) ℚ) - ncE₁)
      + ((1 : Matrix (Fin 3) (Fin 3) ℚ) - ncE₂) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ncL, ncE₁, ncE₂] <;> norm_num

/-- `ncE₁` is symmetric, idempotent, entrywise nonnegative, unital,
and proper (`ncE₁ ≠ 0, 1`). -/
theorem ncE₁_projection_receipt :
    ncE₁.transpose = ncE₁ ∧ ncE₁ * ncE₁ = ncE₁
      ∧ (∀ i j, 0 ≤ ncE₁ i j) ∧ ncE₁.mulVec ncConstant = ncConstant
      ∧ ncE₁ ≠ 0 ∧ ncE₁ ≠ 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · ext i j
    fin_cases i <;> fin_cases j <;> norm_num [ncE₁]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [ncE₁, Matrix.mul_apply, dotProduct, Fin.sum_univ_succ]
  · intro i j
    fin_cases i <;> fin_cases j <;> norm_num [ncE₁]
  · ext i
    fin_cases i <;>
      norm_num [ncE₁, ncConstant, Matrix.mulVec, dotProduct,
        Fin.sum_univ_succ]
  · intro h
    have h00 := congrArg
      (fun M : Matrix (Fin 3) (Fin 3) ℚ => M 0 0) h
    norm_num [ncE₁] at h00
  · intro h
    have h12 := congrArg
      (fun M : Matrix (Fin 3) (Fin 3) ℚ => M 1 2) h
    change (1 / 2 : ℚ) = 0 at h12
    norm_num at h12

/-- `ncE₂` is symmetric, idempotent, entrywise nonnegative, unital,
and proper (`ncE₂ ≠ 0, 1`). -/
theorem ncE₂_projection_receipt :
    ncE₂.transpose = ncE₂ ∧ ncE₂ * ncE₂ = ncE₂
      ∧ (∀ i j, 0 ≤ ncE₂ i j) ∧ ncE₂.mulVec ncConstant = ncConstant
      ∧ ncE₂ ≠ 0 ∧ ncE₂ ≠ 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · ext i j
    fin_cases i <;> fin_cases j <;> norm_num [ncE₂]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [ncE₂, Matrix.mul_apply, dotProduct, Fin.sum_univ_succ]
  · intro i j
    fin_cases i <;> fin_cases j <;> norm_num [ncE₂]
  · ext i
    fin_cases i <;>
      norm_num [ncE₂, ncConstant, Matrix.mulVec, dotProduct,
        Fin.sum_univ_succ]
  · intro h
    have h00 := congrArg
      (fun M : Matrix (Fin 3) (Fin 3) ℚ => M 0 0) h
    norm_num [ncE₂] at h00
  · intro h
    have h02 := congrArg
      (fun M : Matrix (Fin 3) (Fin 3) ℚ => M 0 2) h
    change (1 / 2 : ℚ) = 0 at h02
    norm_num at h02

/-- The two exact conditional expectations do not commute. -/
theorem ncE₁_ncE₂_noncommuting : ncE₁ * ncE₂ ≠ ncE₂ * ncE₁ := by
  intro h
  have h01 := congrArg
    (fun M : Matrix (Fin 3) (Fin 3) ℚ => M 0 1) h
  norm_num [ncE₁, ncE₂, Matrix.mul_apply, dotProduct,
    Fin.sum_univ_succ] at h01

/-- Constants are zero modes of the rate-one repair generator. -/
theorem ncL_constant_eigen : ncL.mulVec ncConstant = 0 := by
  ext i
  fin_cases i <;>
    norm_num [ncL, ncConstant, Matrix.mulVec, dotProduct,
      Fin.sum_univ_succ]

/-- Exact nonconstant eigenpair with eigenvalue `1/2`. -/
theorem ncL_low_eigen :
    ncL.mulVec ncLow = (1 / 2 : ℚ) • ncLow ∧ ncLow ≠ 0 := by
  constructor
  · ext i
    fin_cases i <;>
      norm_num [ncL, ncLow, Matrix.mulVec, dotProduct,
        Fin.sum_univ_succ]
  · intro h
    have h0 := congrFun h (0 : Fin 3)
    norm_num [ncLow] at h0

/-- Exact nonconstant eigenpair with eigenvalue `3/2`. -/
theorem ncL_high_eigen :
    ncL.mulVec ncHigh = (3 / 2 : ℚ) • ncHigh ∧ ncHigh ≠ 0 := by
  constructor
  · ext i
    fin_cases i <;>
      norm_num [ncL, ncHigh, Matrix.mulVec, dotProduct,
        Fin.sum_univ_succ]
  · intro h
    have h0 := congrFun h (0 : Fin 3)
    norm_num [ncHigh] at h0

/-- The two displayed eigenvalues are outside the rate-one
commuting-projection subset sums `{0, 1, 2}`. -/
theorem nc_eigenvalues_not_subset_sums :
    (1 / 2 : ℚ) ≠ 0 ∧ (1 / 2 : ℚ) ≠ 1 ∧ (1 / 2 : ℚ) ≠ 2
      ∧ (3 / 2 : ℚ) ≠ 0 ∧ (3 / 2 : ℚ) ≠ 1 ∧ (3 / 2 : ℚ) ≠ 2 := by
  norm_num

/-! ## Axiom self-audit (build-log visible) -/

#print axioms thm_7_3_finite_gap_nonvacuous
#print axioms wP0_eq_zero
#print axioms repairGenerator_eq
#print axioms ncE₁_projection_receipt
#print axioms ncE₂_projection_receipt
#print axioms ncE₁_ncE₂_noncommuting
#print axioms ncL_eq_rate_one_repairGenerator
#print axioms ncL_low_eigen
#print axioms ncL_high_eigen

end ObserverPatchHolography.YangMillsGapWitness
