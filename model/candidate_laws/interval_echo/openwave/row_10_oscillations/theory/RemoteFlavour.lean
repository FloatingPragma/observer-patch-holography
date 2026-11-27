import Mathlib.Data.Complex.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! A finite remote receiver with inaccessible records of the launch sector. -/

namespace OPH.OpenWave.NeutralFlavour

open scoped BigOperators

abbrev TipWord := Bool × Bool × Bool

inductive Reflection
  | ab | ac | bc
  deriving DecidableEq

def reflect : Reflection → TipWord → TipWord
  | .ab, (a, b, c) => (b, a, c)
  | .ac, (a, b, c) => (c, b, a)
  | .bc, (a, b, c) => (a, c, b)

theorem fixing_reflection (w : TipWord) : ∃ g, reflect g w = w := by
  obtain ⟨a, b, c⟩ := w
  cases a <;> cases b <;> cases c <;>
    first | exact ⟨.ab, rfl⟩ | exact ⟨.ac, rfl⟩ | exact ⟨.bc, rfl⟩

/-- The charged source interchanges the two conjugate preparations under each
arm transposition. Each individual sharp word therefore has equal weight. -/
theorem conjugate_tip_weights {plus minus : TipWord → ℝ}
    (covariant : ∀ g w, plus w = minus (reflect g w)) :
    ∀ w, plus w = minus w := by
  intro w
  obtain ⟨g, hg⟩ := fixing_reflection w
  simpa only [hg] using covariant g w

section Receiver

variable {S E R : Type*} [Fintype S] [Fintype E]

/-- One matrix element after propagation from launch-sector density matrix ρ.
The residual material E is unread. R describes the accessible guide. -/
noncomputable def reduced (V : E → R → S → ℂ) (ρ : S → S → ℂ) (x y : R) : ℂ :=
  ∑ e, ∑ s, ∑ t, V e x s * ρ s t * star (V e y t)

theorem reduced_diagonal (V : E → R → S → ℂ) (tag : E → S)
    (support : ∀ e x s, s ≠ tag e → V e x s = 0)
    (ρ : S → S → ℂ) (x y : R) :
    reduced V ρ x y =
      ∑ e, V e x (tag e) * ρ (tag e) (tag e) * star (V e y (tag e)) := by
  classical
  unfold reduced
  apply Finset.sum_congr rfl
  intro e _
  rw [Finset.sum_eq_single (tag e)]
  · apply Finset.sum_eq_single (tag e)
    · intro j _ hj
      simp [support e y j hj]
    · simp
  · intro i _ hi
    simp [support e x i hi]
  · simp

/-- Includes arbitrary effects on all guide occupations, even if the receiver
interferes different words or the two arriving inputs have different coherences. -/
theorem remote_probability_eq_of_sector_weights [Fintype R]
    (V : E → R → S → ℂ) (tag : E → S)
    (support : ∀ e x s, s ≠ tag e → V e x s = 0)
    (ρ σ : S → S → ℂ) (same : ∀ s, ρ s s = σ s s) (effect : R → R → ℂ) :
    (∑ x, ∑ y, effect y x * reduced V ρ x y) =
    ∑ x, ∑ y, effect y x * reduced V σ x y := by
  have coeff : ∀ x y, reduced V ρ x y = reduced V σ x y := by
    intro x y
    simp only [reduced_diagonal V tag support, same]
  simp only [coeff]

end Receiver

section Classification

variable {O : Type*} [Fintype O]

/-- The receiver sees probabilities p from label zero, q from each other label.
a is the probability of guessing zero on a particular output. The remaining
two guesses have probabilities adding to 1-a. Priors are all 1/3. -/
noncomputable def success (p q a : O → ℝ) : ℝ :=
  (∑ o, (p o * a o + q o * (1 - a o))) / 3

noncomputable def totalVariation (p q : O → ℝ) : ℝ :=
  (∑ o, |p o - q o|) / 2

theorem optimal_success (p q : O → ℝ)
    (hp : ∑ o, p o = 1) (hq : ∑ o, q o = 1) :
    let opt := (1 + totalVariation p q) / 3
    (∀ a : O → ℝ, (∀ o, 0 ≤ a o ∧ a o ≤ 1) → success p q a ≤ opt) ∧
      ∃ a : O → ℝ, (∀ o, a o = 0 ∨ a o = 1) ∧ success p q a = opt := by
  classical
  have maxsum : (∑ o, max (p o) (q o)) = 1 + totalVariation p q := by
    have per : ∀ o, 2 * max (p o) (q o) = p o + q o + |p o - q o| := by
      intro o
      rcases le_total (p o) (q o) with h | h
      · rw [max_eq_right h, abs_of_nonpos (sub_nonpos.mpr h)]
        ring
      · rw [max_eq_left h, abs_of_nonneg (sub_nonneg.mpr h)]
        ring
    have h := calc
      2 * ∑ o, max (p o) (q o) = ∑ o, 2 * max (p o) (q o) := Finset.mul_sum ..
      _ = ∑ o, (p o + q o + |p o - q o|) := Finset.sum_congr rfl fun o _ => per o
      _ = 1 + 1 + ∑ o, |p o - q o| := by rw [Finset.sum_add_distrib, Finset.sum_add_distrib, hp, hq]
    unfold totalVariation
    linarith
  constructor
  · intro a ha
    have h : ∀ o, p o * a o + q o * (1 - a o) ≤ max (p o) (q o) := by
      intro o
      have h₁ := mul_le_mul_of_nonneg_right (le_max_left (p o) (q o)) (ha o).1
      have h₂ := mul_le_mul_of_nonneg_right (le_max_right (p o) (q o)) (sub_nonneg.mpr (ha o).2)
      linarith
    have := Finset.sum_le_sum (s := Finset.univ) (fun o _ => h o)
    unfold success
    linarith
  · refine ⟨fun o => if q o ≤ p o then 1 else 0, fun o => ?_, ?_⟩
    · by_cases h : q o ≤ p o <;> simp [h]
    · have per : ∀ o, p o * (if q o ≤ p o then 1 else 0) +
          q o * (1 - (if q o ≤ p o then 1 else 0)) = max (p o) (q o) := by
        intro o
        by_cases h : q o ≤ p o
        · simp [h]
        · simp [h, max_eq_right (le_of_lt (lt_of_not_ge h))]
      simp only [success, per, maxsum]

theorem equal_prior_ceiling (p q : O → ℝ) (hp₀ : ∀ o, 0 ≤ p o) (hq₀ : ∀ o, 0 ≤ q o)
    (hp : ∑ o, p o = 1) (hq : ∑ o, q o = 1) :
    (1 + totalVariation p q) / 3 ≤ 2 / 3 := by
  have pointwise : ∀ o, |p o - q o| ≤ p o + q o := by
    intro o
    apply abs_le.mpr
    constructor <;> linarith [hp₀ o, hq₀ o]
  have h := Finset.sum_le_sum (s := Finset.univ) (fun o _ => pointwise o)
  rw [Finset.sum_add_distrib, hp, hq] at h
  unfold totalVariation
  linarith

end Classification

end OPH.OpenWave.NeutralFlavour

#print axioms OPH.OpenWave.NeutralFlavour.fixing_reflection
#print axioms OPH.OpenWave.NeutralFlavour.conjugate_tip_weights
#print axioms OPH.OpenWave.NeutralFlavour.reduced_diagonal
#print axioms OPH.OpenWave.NeutralFlavour.remote_probability_eq_of_sector_weights
#print axioms OPH.OpenWave.NeutralFlavour.optimal_success
#print axioms OPH.OpenWave.NeutralFlavour.equal_prior_ceiling
