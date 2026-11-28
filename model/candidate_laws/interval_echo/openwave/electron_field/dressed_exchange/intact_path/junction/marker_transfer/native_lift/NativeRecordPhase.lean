import Mathlib.Analysis.Complex.Basic
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic.LinearCombination

/-!
# Unsigned configuration transport and the recovered phase reference

The sum of sharp-configuration amplitudes is a conserved linear functional
for the unsigned permutations on a native orbit. This is a proof functional.
The phase comparison is a separate controlled-arm measurement.
-/

namespace NativeRecordPhase

open scoped BigOperators

variable {Ω : Type*} [Fintype Ω]

/-- Native unsigned relabeling preserves the configuration sum. -/
theorem perm_preserves_sum (p : Equiv.Perm Ω) (x : Ω → ℂ) :
    ∑ i, x (p i) = ∑ i, x i := by
  exact Fintype.sum_equiv p _ _ (fun _ => rfl)

/-- A native addressed term kills the conserved functional. -/
theorem laplacian_has_zero_sum (p : Equiv.Perm Ω) (x : Ω → ℂ) :
    ∑ i, (x i - x (p i)) = 0 := by
  rw [Finset.sum_sub_distrib, perm_preserves_sum, sub_self]

/-- If both arms preserve the sum, their balanced bright port preserves it. -/
theorem bright_retains_sum (x a b : Ω → ℂ)
    (ha : ∑ i, a i = ∑ i, x i) (hb : ∑ i, b i = ∑ i, x i) :
    ∑ i, ((a i + b i) / 2) = ∑ i, x i := by
  rw [← Finset.sum_div, Finset.sum_add_distrib, ha, hb]
  ring

/-- An exact negative comparison is impossible when the preserved
functional of the actual input is nonzero. -/
theorem no_negative_on_observed_input (x a b : Ω → ℂ)
    (ha : ∑ i, a i = ∑ i, x i) (hb : ∑ i, b i = ∑ i, x i)
    (h : (∑ i, x i) ≠ 0) :
    ¬ (∀ i, a i = -(b i)) := by
  intro hn
  have hz : ∑ i, a i = -(∑ i, b i) := by
    simp only [hn, Finset.sum_neg_distrib]
  rw [ha, hb] at hz
  exact h (by linear_combination hz / 2)

#print axioms perm_preserves_sum
#print axioms laplacian_has_zero_sum
#print axioms bright_retains_sum
#print axioms no_negative_on_observed_input

end NativeRecordPhase
