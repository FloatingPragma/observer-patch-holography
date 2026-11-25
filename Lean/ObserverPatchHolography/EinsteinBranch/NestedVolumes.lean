import ObserverPatchHolography.EinsteinBranch.RegionalMoments

/-! Matched proper-volume shells and the radius difference reconstructed from their areas. -/

namespace OPH.EinsteinBranch.NestedVolumes

open scoped BigOperators

/-- An effective shell area is the matched proper volume increment divided by
its radial thickness. The independently measured area can have a residual.
The balance residual here is in nats and includes any finite-amplitude remainder. -/
theorem radius_contrast {ι : Type*} [Fintype ι]
    (w a b s e : ι → ℝ) (η r₀ r₁ c : ℝ) (hη : 0 < η)
    (ha : ∀ i, 0 < a i) (hb : ∀ i, 0 < b i)
    (hbalance : ∀ i, s i + η * (b i - a i) = e i)
    (hchart : r₁ - r₀ = c + ∑ i, (w i / b i - w i / a i)) :
    r₁ - r₀ = c + ∑ i, w i * (s i - e i) / (η * a i * b i) := by
  rw [hchart]
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  rw [← hbalance i]
  field_simp [ne_of_gt hη, ne_of_gt (ha i), ne_of_gt (hb i)]
  ring

/-- A surplus on the shells together with a bound on the independent chart
residual gives an experimentally testable lower bound for the change in radius. -/
theorem radius_lower {ι : Type*} [Fintype ι]
    (w a b s e budget : ι → ℝ) (η r₀ r₁ c ε : ℝ) (hη : 0 < η)
    (ha : ∀ i, 0 < a i) (hb : ∀ i, 0 < b i) (hw : ∀ i, 0 ≤ w i)
    (he : ∀ i, e i ≤ budget i) (hc : -ε ≤ c)
    (hbalance : ∀ i, s i + η * (b i - a i) = e i)
    (hchart : r₁ - r₀ = c + ∑ i, (w i / b i - w i / a i)) :
    -ε + ∑ i, w i * (s i - budget i) / (η * a i * b i) ≤ r₁ - r₀ := by
  rw [radius_contrast w a b s e η r₀ r₁ c hη ha hb hbalance hchart]
  apply add_le_add hc
  apply Finset.sum_le_sum
  intro i _
  exact div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left (sub_le_sub_left (he i) _) (hw i))
    (le_of_lt (mul_pos (mul_pos hη (ha i)) (hb i)))

assert_no_sorry radius_contrast
assert_no_sorry radius_lower
#print axioms radius_contrast
#print axioms radius_lower

end OPH.EinsteinBranch.NestedVolumes
