import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Data.Real.Sqrt
import Mathlib.Tactic

namespace PolarizedRanges

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The four squared cross distances of two pairs antipodal about one center.
The angle between their axes is arbitrary. -/
theorem four_cross_distances (x y : E) :
    dist x y ^ 2 + dist x (-y) ^ 2 + dist (-x) y ^ 2 + dist (-x) (-y) ^ 2 =
      4 * (‖x‖ ^ 2 + ‖y‖ ^ 2) := by
  have ha : -x - y = -(x + y) := by abel
  have hs : -x - -y = -(x - y) := by abel
  simp only [dist_eq_norm, ha, hs, norm_neg, sub_neg_eq_add]
  linarith [parallelogram_law_with_norm ℝ x y]

/-- In measured exchange units, the center ranges are 1 and 2 and each
cross range is 3. Every individual proper range may move by its own
fractional error bounded by e. The endpoints are genuine tangent
antipodes; no assumption about their angle or RMS distance is needed. -/
theorem antipodal_fit_error (x y : E) {e : ℝ} (he : 0 ≤ e)
    (hx : |‖x‖ - 1| ≤ e) (hy : |‖y‖ - 2| ≤ 2 * e)
    (hpp : |dist x y - 3| ≤ 3 * e) (hpm : |dist x (-y) - 3| ≤ 3 * e)
    (hmp : |dist (-x) y - 3| ≤ 3 * e) (hmm : |dist (-x) (-y) - 3| ≤ 3 * e) :
    (3 - Real.sqrt 5) / (3 + Real.sqrt 5) ≤ e := by
  have hs := Real.sqrt_nonneg 5
  have hs2 : Real.sqrt 5 ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  apply (div_le_iff₀ (by positivity)).mpr
  by_cases h : 1 ≤ e
  · nlinarith [mul_nonneg hs he]
  have hlow : 0 ≤ 3 * (1 - e) := by linarith
  have cross {d : ℝ} (hd : |d - 3| ≤ 3 * e) : (3 * (1 - e)) ^ 2 ≤ d ^ 2 := by
    have hle : 3 * (1 - e) ≤ d := by linarith [(abs_le.mp hd).1]
    exact (sq_le_sq₀ hlow (hlow.trans hle)).mpr hle
  have xr : ‖x‖ ^ 2 ≤ (1 + e) ^ 2 := by
    apply (sq_le_sq₀ (norm_nonneg _) (by positivity)).mpr
    linarith [(abs_le.mp hx).2]
  have yr : ‖y‖ ^ 2 ≤ (2 * (1 + e)) ^ 2 := by
    apply (sq_le_sq₀ (norm_nonneg _) (by positivity)).mpr
    linarith [(abs_le.mp hy).2]
  have bound : (3 * (1 - e)) ^ 2 ≤ 5 * (1 + e) ^ 2 := by
    nlinarith only [xr, yr, cross hpp, cross hpm, cross hmp, cross hmm,
      four_cross_distances x y]
  have rhs : (Real.sqrt 5 * (1 + e)) ^ 2 = 5 * (1 + e) ^ 2 := by rw [mul_pow, hs2]
  rw [← rhs] at bound
  have linear := le_of_sq_le_sq bound (by positivity)
  linarith

/-- Four independent candidate endpoints, translated so their reported
center is zero. Expanding this squared norm tests every one of the ten
distances without assuming that any diameter stays centered. -/
theorem free_endpoints_identity (a b c d : E) :
    12 * (‖a‖^2 + ‖b‖^2) + 6 * (‖c‖^2 + ‖d‖^2) -
      4 * dist a b ^ 2 - dist c d ^ 2 -
      2 * (dist a c ^ 2 + dist a d ^ 2 + dist b c ^ 2 + dist b d ^ 2) =
      ‖(a+b) + (a+b) + (c+d)‖^2 := by
  simp only [dist_eq_norm, norm_add_sq_real, norm_sub_sq_real, inner_add_left,
    inner_add_right, real_inner_self_eq_norm_sq, real_inner_comm b a,
    real_inner_comm c a, real_inner_comm c b, real_inner_comm d a,
    real_inner_comm d b, real_inner_comm d c]
  ring

/-- A position-independent lower bound from all ten observed ranges. -/
theorem free_endpoints_fit_error (a b c d : E) {e : ℝ} (he : 0 ≤ e)
    (ha : |‖a‖ - 1| ≤ e) (hb : |‖b‖ - 1| ≤ e)
    (hc : |‖c‖ - 2| ≤ 2*e) (hd : |‖d‖ - 2| ≤ 2*e)
    (hab : |dist a b - 2| ≤ 2*e) (hcd : |dist c d - 4| ≤ 4*e)
    (hac : |dist a c - 3| ≤ 3*e) (had : |dist a d - 3| ≤ 3*e)
    (hbc : |dist b c - 3| ≤ 3*e) (hbd : |dist b d - 3| ≤ 3*e) :
    (Real.sqrt 13 - 3) / (Real.sqrt 13 + 3) ≤ e := by
  have hs := Real.sqrt_nonneg 13
  have hs2 : Real.sqrt 13 ^ 2 = 13 := Real.sq_sqrt (by norm_num)
  apply (div_le_iff₀ (by positivity)).mpr
  by_cases h : 1 ≤ e
  · nlinarith [mul_nonneg hs (sub_nonneg.mpr h)]
  have hlo {l z : ℝ} (hl : 0 ≤ l) (hz : |z-l| ≤ l*e) :
      (l*(1-e))^2 ≤ z^2 := by
    have nonneg : 0 ≤ l*(1-e) := mul_nonneg hl (by linarith)
    have hle : l*(1-e) ≤ z := by linarith [(abs_le.mp hz).1]
    exact (sq_le_sq₀ nonneg (nonneg.trans hle)).mpr hle
  have hhi {l z : ℝ} (hl : 0 ≤ l) (hzn : 0 ≤ z) (hz : |z-l| ≤ l*e) :
      z^2 ≤ (l*(1+e))^2 := by
    apply (sq_le_sq₀ hzn (by positivity)).mpr
    linarith [(abs_le.mp hz).2]
  have hpa := hhi (show 0 ≤ (1:ℝ) by norm_num) (norm_nonneg a) (by simpa using ha)
  have hpb := hhi (show 0 ≤ (1:ℝ) by norm_num) (norm_nonneg b) (by simpa using hb)
  have hpc := hhi (show 0 ≤ (2:ℝ) by norm_num) (norm_nonneg c) hc
  have hpd := hhi (show 0 ≤ (2:ℝ) by norm_num) (norm_nonneg d) hd
  have hsq : 13*(1-e)^2 ≤ 9*(1+e)^2 := by
    nlinarith only [hpa,hpb,hpc,hpd,hlo (by norm_num) hab,hlo (by norm_num) hcd,
      hlo (by norm_num) hac,hlo (by norm_num) had,hlo (by norm_num) hbc,hlo (by norm_num) hbd,
      free_endpoints_identity a b c d, sq_nonneg ‖(a+b)+(a+b)+(c+d)‖]
  have lhs : (Real.sqrt 13 * (1-e))^2 = 13*(1-e)^2 := by rw [mul_pow, hs2]
  have rhs : (3*(1+e))^2 = 9*(1+e)^2 := by ring
  rw [← lhs, ← rhs] at hsq
  have linear := le_of_sq_le_sq hsq (by positivity)
  linarith

private abbrev Point := ℝ × ℝ × ℝ

private def squaredDistance (u v : Point) : ℝ :=
  (u.1 - v.1)^2 + (u.2.1 - v.2.1)^2 + (u.2.2 - v.2.2)^2

/-- Without exact antipodal placement, a Euclidean survey can move every
central range, diameter and cross range by at most one tenth of its value.
The squared inequalities avoid rounding any of the actual lengths. -/
theorem unconstrained_fit_at_one_tenth :
    let o : Point := (0, 0, 0)
    let a : Point := (9/10, 0, 3/5)
    let b : Point := (-9/10, 0, 3/5)
    let c : Point := (0, 9/5, -6/5)
    let d : Point := (0, -9/5, -6/5)
    ∀ entry ∈ ([(o,a,1), (o,b,1), (o,c,2), (o,d,2), (a,b,2), (c,d,4),
        (a,c,3), (a,d,3), (b,c,3), (b,d,3)] : List (Point × Point × ℝ)),
      (9/10 * entry.2.2)^2 ≤ squaredDistance entry.1 entry.2.1 ∧
      squaredDistance entry.1 entry.2.1 ≤ (11/10 * entry.2.2)^2 := by
  dsimp only
  intro entry h
  simp only [List.mem_cons, List.not_mem_nil, or_false] at h
  rcases h with h|h|h|h|h|h|h|h|h|h <;> subst entry <;> norm_num [squaredDistance]

end PolarizedRanges

#print axioms PolarizedRanges.four_cross_distances
#print axioms PolarizedRanges.antipodal_fit_error
#print axioms PolarizedRanges.free_endpoints_identity
#print axioms PolarizedRanges.free_endpoints_fit_error
#print axioms PolarizedRanges.unconstrained_fit_at_one_tenth
