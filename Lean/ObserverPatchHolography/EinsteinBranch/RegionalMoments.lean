import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall

/-! Finite local stress moments sufficient for the small-ball bulk estimate. -/

namespace OPH.EinsteinBranch.RegionalMoments

open scoped BigOperators

/-- Weights may have signs. The representation error permits unresolved microscopic variations. -/
theorem bulk_error {ι : Type*} [Fintype ι] (w t : ι → ℝ)
    (s m center osc gain variation leakage : ℝ)
    (hrep : |s - ∑ i, w i * t i| ≤ leakage)
    (hm : |(∑ i, w i) - m| ≤ gain)
    (htv : (∑ i, |w i|) ≤ variation)
    (hosc : ∀ i, |t i - center| ≤ osc) (hosc0 : 0 ≤ osc) :
    |s - m * center| ≤ leakage + variation * osc + gain * |center| := by
  have hterms : |∑ i, w i * (t i - center)| ≤ variation * osc := calc
    |∑ i, w i * (t i - center)| ≤ ∑ i, |w i * (t i - center)| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, |w i| * osc := Finset.sum_le_sum fun i _ => by
      rw [abs_mul]; exact mul_le_mul_of_nonneg_left (hosc i) (abs_nonneg _)
    _ = (∑ i, |w i|) * osc := by rw [Finset.sum_mul]
    _ ≤ variation * osc := mul_le_mul_of_nonneg_right htv hosc0
  have hsplit : s - m * center =
      (s - ∑ i, w i * t i) + (∑ i, w i * (t i - center)) + ((∑ i, w i) - m) * center := by
    simp_rw [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul]
    ring
  rw [hsplit]
  calc
    _ ≤ |s - ∑ i, w i * t i| + |∑ i, w i * (t i - center)| + |((∑ i, w i) - m) * center| :=
      (abs_add_le _ _).trans (add_le_add (abs_add_le _ _) (le_refl _))
    _ ≤ leakage + variation * osc + gain * |center| := by
      rw [abs_mul]
      exact add_le_add (add_le_add hrep hterms) (mul_le_mul_of_nonneg_right hm (abs_nonneg _))

/-- A common spatial energy profile gives one forecast, regardless of the spatial weights. -/
theorem shared_profile_error (s₁ s₂ forecast e : ℝ)
    (h₁ : |s₁ - forecast| ≤ e) (h₂ : |s₂ - forecast| ≤ e) :
    |s₁ - s₂| / 2 ≤ e := by
  have h := abs_sub_le s₁ forecast s₂
  rw [abs_sub_comm forecast s₂] at h
  linarith

/-- Opposite feasible perturbations turn a nonnegative conditional information
and a quadratic entropy cost into a signed first-variation error bound. -/
theorem conditional_first_variation (info slope step cost : ℝ) (ht : 0 < step)
    (hp : 0 ≤ info + step * slope + step ^ 2 * cost)
    (hm : 0 ≤ info - step * slope + step ^ 2 * cost) :
    |slope| ≤ info / step + step * cost := by
  have h : |slope| ≤ (info + step ^ 2 * cost) / step := by
    apply (le_div_iff₀ ht).mpr
    rcases le_total 0 slope with h | h
    · rw [abs_of_nonneg h]
      linarith [mul_comm step slope]
    · rw [abs_of_nonpos h, neg_mul]
      linarith [mul_comm step slope]
  convert h using 1
  field_simp

/-- A thermal gain common to two ball radii has at least one third relative
error against a target gain proportional to the radius. The bound is sharp. -/
theorem two_radius_error (x e : ℝ)
    (h₁ : |x - 1| ≤ e) (h₂ : |x / 2 - 1| ≤ e) : (1:ℝ)/3 ≤ e := by
  have hsmall := (abs_le.mp h₁).2
  have hlarge := (abs_le.mp h₂).1
  linarith

assert_no_sorry conditional_first_variation
assert_no_sorry two_radius_error
#print axioms conditional_first_variation
#print axioms two_radius_error
assert_no_sorry bulk_error
assert_no_sorry shared_profile_error

#print axioms bulk_error
#print axioms shared_profile_error

end OPH.EinsteinBranch.RegionalMoments
