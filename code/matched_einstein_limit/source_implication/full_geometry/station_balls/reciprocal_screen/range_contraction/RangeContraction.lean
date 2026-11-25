import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

namespace RangeContraction

/-- The observed log-length sum is bounded by the entropy deficit. The three
input equalities are respectively the finite relative-entropy identity, the
affine clock effect at the calibrated reference, and the log clock correction. -/
theorem weighted_log_bound {n m : ℕ} (w z : Fin n → ℝ) (b l : Fin m → ℝ)
    (entropyChange relativeEntropy modularChange : ℝ)
    (hw : ∀ j, 0 ≤ w j) (hD : 0 ≤ relativeEntropy)
    (hS : entropyChange = modularChange - relativeEntropy)
    (hK : modularChange = -(∑ j, w j * (Real.exp (z j) - 1)))
    (hL : (∑ j, w j * z j) = ∑ e, b e * l e) :
    (∑ e, b e * l e) ≤ -entropyChange := by
  have h : (∑ j, w j * z j) ≤ ∑ j, w j * (Real.exp (z j) - 1) := by
    apply Finset.sum_le_sum
    intro j _
    apply mul_le_mul_of_nonneg_left _ (hw j)
    linarith [Real.add_one_le_exp (z j)]
  linarith

/-- With positive total edge weight, an entropy gain strictly contracts
the weighted geometric mean of the reconstructed relative lengths. -/
theorem geometric_mean_contracts (signedLog total gain : ℝ) (hB : 0 < total)
    (hG : 0 < gain) (h : signedLog ≤ -gain) :
    Real.exp (signedLog / total) < 1 := by
  apply Real.exp_lt_one_iff.mpr
  exact div_neg_of_neg_of_pos (by linarith) hB

end RangeContraction
