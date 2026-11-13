import ObserverPatchHolography.EinsteinBranch.Tensor

/-!
# Error transport through the existing finite null-tomography decoder

The tensor decoder is already proved to invert its nine null charges. This
module supplies the finite weighted-error inequality needed to propagate
charge intervals through each decoder row. Physical stress/metric attachment
and a measured error budget remain hypotheses, not consequences of this bound.
-/

namespace OPH.EinsteinBranch

open scoped BigOperators

theorem finite_weighted_error_bound {n : ℕ} (a e : Fin n → ℝ) (ε : ℝ)
    (h : ∀ i, |e i| ≤ ε) :
    |∑ i, a i * e i| ≤ (∑ i, |a i|) * ε := by
  calc
    |∑ i, a i * e i| ≤ ∑ i, |a i * e i| := by
      simpa using Finset.abs_sum_le_sum_abs (fun i => a i * e i) Finset.univ
    _ = ∑ i, |a i| * |e i| := by simp only [abs_mul]
    _ ≤ ∑ i, |a i| * ε := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_left (h i) (abs_nonneg (a i))
    _ = (∑ i, |a i|) * ε := by rw [Finset.sum_mul]

theorem finite_decoder_difference_bound {n : ℕ}
    (a q measured : Fin n → ℝ) (ε : ℝ)
    (h : ∀ i, |measured i - q i| ≤ ε) :
    |(∑ i, a i * measured i) - (∑ i, a i * q i)| ≤
      (∑ i, |a i|) * ε := by
  have heq : (∑ i, a i * measured i) - (∑ i, a i * q i) =
      ∑ i, a i * (measured i - q i) := by
    simp only [mul_sub, Finset.sum_sub_distrib]
  rw [heq]
  exact finite_weighted_error_bound a (fun i => measured i - q i) ε h

#print axioms finite_weighted_error_bound
#print axioms finite_decoder_difference_bound

end OPH.EinsteinBranch
