import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Analysis.Complex.Basic
import Mathlib.Tactic.Linarith

/-! Marginal interference terms and finite postprocessing of a signed click. -/
namespace OPH.NeutralChargedTags

open scoped BigOperators

variable {D E L R : Type*}

/-- A single common total and different sharp daughter charges force
    disjoint donor support. Charges include every unreported compensator. -/
theorem split_support (a b : D → E → ℂ) (qa qb total : ℤ) (q : E → ℤ)
    (ha : ∀ x e, a x e ≠ 0 → qa + q e = total)
    (hb : ∀ y e, b y e ≠ 0 → qb + q e = total)
    (different : qa ≠ qb) (x y : D) (e : E) :
    a x e * star (b y e) = 0 := by
  by_cases h : a x e = 0
  · simp [h]
  by_cases k : b y e = 0
  · simp [k]
  have p := ha x e h
  have s := hb y e k
  exact (different (by linarith)).elim

/-- Every daughter effect only contracts this marginal cross term. -/
theorem marginal_zero [Fintype E] (a b : D → E → ℂ) (qa qb total : ℤ) (q : E → ℤ)
    (ha : ∀ x e, a x e ≠ 0 → qa + q e = total)
    (hb : ∀ y e, b y e ≠ 0 → qb + q e = total)
    (different : qa ≠ qb) (x y : D) :
    ∑ e, a x e * star (b y e) = 0 := by
  exact Finset.sum_eq_zero fun e _ => split_support a b qa qb total q ha hb different x y e

/-- A detector computed only from the accessible tag cannot increase the
    absolute difference of outcome probabilities, summed over all exits. -/
theorem tag_bound [Fintype L] [Fintype R] (K : L → R → ℝ) (d : L → ℝ) (nonneg : ∀ l r, 0 ≤ K l r)
    (stochastic : ∀ l, ∑ r, K l r = 1) :
    (∑ r, |∑ l, d l * K l r|) ≤ ∑ l, |d l| := by
  calc
    _ ≤ ∑ r, ∑ l, |d l * K l r| := by
      apply Finset.sum_le_sum
      intro r _
      exact Finset.abs_sum_le_sum_abs _ _
    _ = ∑ l, ∑ r, |d l| * K l r := by
      rw [Finset.sum_comm]
      simp only [abs_mul, abs_of_nonneg (nonneg _ _)]
    _ = _ := by simp only [← Finset.mul_sum, stochastic, mul_one]

#print axioms split_support
#print axioms marginal_zero
#print axioms tag_bound

end OPH.NeutralChargedTags
