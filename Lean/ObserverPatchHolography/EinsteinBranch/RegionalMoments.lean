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

/-- A signed local stress wholly inside a ball contributes two common moments.
The radial charge is in units with pi absorbed in the moments. -/
noncomputable def twoMomentCharge (r e m : ℝ) : ℝ := r * e - m / r

theorem three_radius_moment_annihilator (r e m : ℝ) (hr : 0 < r) :
    4 * twoMomentCharge r e m - 17 * twoMomentCharge (4*r) e m
      + 4 * twoMomentCharge (16*r) e m = 0 := by
  unfold twoMomentCharge
  have h : r ≠ 0 := ne_of_gt hr
  field_simp
  ring

/-- Flat diamond kernel divided by pi, with d the recorded proper distance. -/
noncomputable def flatDiamondWeight (r d : ℝ) : ℝ :=
  if d ≤ r then r - d^2 / r else 0

private theorem flatDiamondWeight_exterior (r d : ℝ) (hr : 0 < r) (hd : r ≤ d) :
    flatDiamondWeight r d = 0 := by
  rcases eq_or_lt_of_le hd with h | h
  · subst d
    simp [flatDiamondWeight, sq, hr.ne']
  · simp [flatDiamondWeight, not_le.mpr h]

private theorem three_radius_atom (r d : ℝ) (hr : 0 < r) (hd : d ≤ r ∨ 16*r ≤ d) :
    4 * flatDiamondWeight r d - 17 * flatDiamondWeight (4*r) d
      + 4 * flatDiamondWeight (16*r) d = 0 := by
  rcases hd with h | h
  · have h₄ : d ≤ 4*r := by linarith
    have h₁₆ : d ≤ 16*r := by linarith
    simpa [flatDiamondWeight, h, h₄, h₁₆, twoMomentCharge] using
      three_radius_moment_annihilator r 1 (d^2) hr
  · rw [flatDiamondWeight_exterior r d hr (by linarith),
        flatDiamondWeight_exterior (4*r) d (by linarith) (by linarith),
        flatDiamondWeight_exterior (16*r) d (by linarith) h]
    norm_num

/-- Arbitrary signed atoms inside the first ball or outside the last give the
three-radius identity. Atoms on either bounding surface are included. -/
theorem three_radius_finite_stress {ι : Type*} [Fintype ι] (r : ℝ) (d w : ι → ℝ)
    (hr : 0 < r) (support : ∀ i, d i ≤ r ∨ 16*r ≤ d i) :
    4 * (∑ i, w i * flatDiamondWeight r (d i))
      - 17 * (∑ i, w i * flatDiamondWeight (4*r) (d i))
      + 4 * (∑ i, w i * flatDiamondWeight (16*r) (d i)) = 0 := by
  simp only [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  apply Finset.sum_eq_zero
  intro i _
  calc
    _ = w i * (4 * flatDiamondWeight r (d i) - 17 * flatDiamondWeight (4*r) (d i)
          + 4 * flatDiamondWeight (16*r) (d i)) := by ring
    _ = 0 := by rw [three_radius_atom r (d i) hr (support i), mul_zero]

/-- A common entropy slope at least lower * R⁴ cannot be represented by
these charges to error smaller than (9 / 266500) * lower at every radius,
when each error is normalized by that ball's own radius to the fourth. -/
theorem three_radius_normalized_error (r s b₁ b₄ b₁₆ lower error : ℝ) (hr : 0 < r)
    (hs : lower*r^4 ≤ s) (hB : 4*b₁ - 17*b₄ + 4*b₁₆ = 0)
    (h₁ : |s-b₁| ≤ error*r^4) (h₄ : |s-b₄| ≤ error*(4*r)^4)
    (h₁₆ : |s-b₁₆| ≤ error*(16*r)^4) :
    9*lower/266500 ≤ error := by
  have low₁ := (abs_le.mp h₁).1
  have high₄ := (abs_le.mp h₄).2
  have low₁₆ := (abs_le.mp h₁₆).1
  have hscaled : (9*lower)*r^4 ≤ (266500*error)*r^4 := by
    linarith only [hs, hB, low₁, high₄, low₁₆]
  have h := (mul_le_mul_iff_left₀ (pow_pos hr 4)).mp hscaled
  linarith

assert_no_sorry three_radius_moment_annihilator
assert_no_sorry three_radius_finite_stress
assert_no_sorry three_radius_normalized_error
#print axioms three_radius_moment_annihilator
#print axioms three_radius_finite_stress
#print axioms three_radius_normalized_error

assert_no_sorry conditional_first_variation
assert_no_sorry two_radius_error
#print axioms conditional_first_variation
#print axioms two_radius_error
assert_no_sorry bulk_error
assert_no_sorry shared_profile_error

#print axioms bulk_error
#print axioms shared_profile_error

end OPH.EinsteinBranch.RegionalMoments
