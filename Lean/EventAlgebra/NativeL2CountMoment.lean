import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Data.Real.StarOrdered
import Mathlib.Tactic

/-!
The algebraic part of the bounded L=2 source attachment. These are standard identities,
with no additional Born representation claim. Concrete source counts and the actual downstream
Python return are separate source-reading/executable obligations.
-/

namespace EventAlgebra.NativeBinarySource.L2

theorem count_first_ratio (tail total population denominator coefficient : ℚ)
    (hTotal : total ≠ 0) (hPopulation : population ≠ 0)
    (hDenominator : denominator ≠ 0) (hCoefficient : coefficient ≠ 0) :
    (coefficient * tail / (population * denominator^2)) /
      (coefficient * (total / (population * denominator^2))) = tail / total := by
  field_simp

theorem normalized_gram_positive {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (Z : ℝ) (hZ : 0 < Z) : (Z⁻¹ • (C.transpose * C)).PosSemidef := by
  have h : C.conjTranspose = C.transpose := by ext i j; simp
  rw [← h]
  exact (Matrix.posSemidef_conjTranspose_mul_self C).smul (inv_nonneg.mpr hZ.le)

theorem normalized_square_positive {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (hC : C.transpose = C) (Z : ℝ) (hZ : 0 < Z) :
    (Z⁻¹ • (C*C)).PosSemidef := by
  have h : C*C = C.transpose*C := congrArg (fun A => A*C) hC.symm
  rw [h]
  exact normalized_gram_positive C Z hZ

theorem normalized_trace {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (Z : ℝ) (hZ : 0 < Z) (hTrace : (C*C).trace = Z) :
    (Z⁻¹ • (C*C)).trace = 1 := by
  rw [Matrix.trace_smul, hTrace]
  exact inv_mul_cancel₀ (ne_of_gt hZ)

theorem witness_ratio :
    (13/1782000 : ℚ) / ((1/60)*(159/220)) = 13/21465 := by norm_num

theorem diagonal_erasure_mismatch : (0 : ℚ) ≠ 13/21465 := by norm_num

theorem coefficient_omission_mismatch :
    (13/1782000 : ℚ)/(159/220) = 13/1287900 ∧
    (13/1782000 : ℚ)/(159/220) ≠ 13/21465 := by norm_num

#print axioms count_first_ratio
#print axioms normalized_gram_positive
#print axioms normalized_square_positive
#print axioms normalized_trace
#print axioms witness_ratio
#print axioms diagonal_erasure_mismatch
#print axioms coefficient_omission_mismatch

end EventAlgebra.NativeBinarySource.L2
