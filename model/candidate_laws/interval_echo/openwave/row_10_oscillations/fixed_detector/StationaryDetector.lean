import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Linarith

namespace OPH.StationaryNeutralDetector

open scoped BigOperators

variable {n i : Type*} [Fintype n] [DecidableEq n] [Fintype i]

/-- A complete source propagates to the identity under a common unitary.
    The stored effect may be lossy, biased, or higher rank. -/
theorem total_effect (u v e : Matrix n n ℂ) (s : i → Matrix n n ℂ)
    (hu : u * v = 1) (hs : ∑ j, s j = 1) :
    (∑ j, (e * (u * s j * v)).trace) = e.trace := by
  rw [← Matrix.trace_sum, ← Finset.mul_sum, ← Finset.sum_mul, ← Finset.mul_sum,
      hs, mul_one, hu, mul_one]

/-- The normalized mixture of conditional click scores is in their interval. -/
theorem mixture_range (w x : i → ℝ) (lo hi : ℝ)
    (hw : ∀ j, 0 ≤ w j) (hs : ∑ j, w j = 1)
    (bounds : ∀ j, lo ≤ x j ∧ x j ≤ hi) :
    lo ≤ ∑ j, w j * x j ∧ ∑ j, w j * x j ≤ hi := by
  constructor
  · calc
      lo = ∑ j, w j * lo := by rw [← Finset.sum_mul, hs, one_mul]
      _ ≤ _ := Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (bounds j).1 (hw j)
  · calc
      _ ≤ ∑ j, w j * hi := Finset.sum_le_sum fun j _ =>
        mul_le_mul_of_nonneg_left (bounds j).2 (hw j)
      _ = hi := by rw [← Finset.sum_mul, hs, one_mul]

theorem separated (a b : i → ℝ) (wa wb : i → ℝ)
    (l₁ h₁ l₂ h₂ : ℝ) (order : h₁ < l₂)
    (pa : ∀ j, 0 ≤ wa j) (pb : ∀ j, 0 ≤ wb j)
    (sa : ∑ j, wa j = 1) (sb : ∑ j, wb j = 1)
    (ba : ∀ j, l₁ ≤ a j ∧ a j ≤ h₁) (bb : ∀ j, l₂ ≤ b j ∧ b j ≤ h₂) :
    (∑ j, wa j * a j) < ∑ j, wb j * b j :=
  lt_of_le_of_lt (mixture_range wa a l₁ h₁ pa sa ba).2
    (lt_of_lt_of_le order (mixture_range wb b l₂ h₂ pb sb bb).1)

#print axioms total_effect
#print axioms mixture_range
#print axioms separated

end OPH.StationaryNeutralDetector
