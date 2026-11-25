import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Util.AssertNoSorry

/-! Between two references the weakening slope of a balanced bit has a
    minimum spread. A bounded geometric spread leaves an endpoint residual. -/

namespace OPH.EinsteinBranch.AreaReadoutSpread

open Set

noncomputable def weakSlope (k s : ℝ) : ℝ :=
  k * ((Real.log (1 + k*s) - Real.log (1 - k*s))/2)

theorem hasDerivAt_weakSlope {k s : ℝ} (h : |k*s| < 1) :
    HasDerivAt (weakSlope k) (k^2 / (1-(k*s)^2)) s := by
  obtain ⟨hl, hr⟩ := abs_lt.mp h
  have hk : HasDerivAt (fun t : ℝ => k*t) k s := by
    simpa using (hasDerivAt_id s).const_mul k
  have hp := (hk.const_add 1).log (by dsimp; linarith : 1 + k*s ≠ 0)
  have hm := (hk.const_sub 1).log (by dsimp; linarith : 1 - k*s ≠ 0)
  have plus : 0 < 1 + k*s := by linarith
  have minus : 0 < 1 - k*s := by linarith
  have hden : 1 - (k*s)^2 ≠ 0 := ne_of_gt (by nlinarith [mul_pos plus minus])
  have hfac : 1 - k^2*s^2 ≠ 0 := by simpa [mul_pow] using hden
  convert (((hp.sub hm).div_const 2).const_mul k) using 1
  simp only [mul_pow]
  field_simp [hfac, (show 1+k*s ≠ 0 by linarith), (show 1-k*s ≠ 0 by linarith)]
  ring

theorem weakSlope_spread {k a b : ℝ} (hab : a ≤ b)
    (h : ∀ s ∈ Icc a b, |k*s| < 1) :
    k^2 * (b-a) ≤ weakSlope k b - weakSlope k a := by
  have hf : ∀ s ∈ Icc a b, HasDerivAt (weakSlope k) (k^2 / (1-(k*s)^2)) s :=
    fun s hs => hasDerivAt_weakSlope (h s hs)
  apply (convex_Icc a b).mul_sub_le_image_sub_of_le_deriv
    (fun s hs => (hf s hs).continuousAt.continuousWithinAt)
    (fun s hs => (hf s (interior_subset hs)).differentiableAt.differentiableWithinAt)
    _ a (left_mem_Icc.mpr hab) b (right_mem_Icc.mpr hab) hab
  intro s hs
  have ht := h s (interior_subset hs)
  rw [(hf s (interior_subset hs)).deriv]
  obtain ⟨hl, hr⟩ := abs_lt.mp ht
  have hden : 0 < 1 - (k*s)^2 := by
    nlinarith [mul_pos (show 0 < 1+k*s by linarith) (show 0 < 1-k*s by linarith)]
  apply (le_div_iff₀ hden).mpr
  nlinarith [mul_nonneg (sq_nonneg k) (sq_nonneg (k*s))]

/-- The area slopes include the common calibration. The same tolerance must
    cover both endpoint residuals. -/
theorem endpoint_residual {k a b areaA areaB budget err : ℝ} (hab : a ≤ b)
    (h : ∀ s ∈ Icc a b, |k*s| < 1) (harea : |areaB-areaA| ≤ budget)
    (ha : |weakSlope k a-areaA| ≤ err) (hb : |weakSlope k b-areaB| ≤ err) :
    (k^2 * (b-a) - budget)/2 ≤ err := by
  have hspread := weakSlope_spread hab h
  linarith [(abs_le.mp harea).2, (abs_le.mp ha).1, (abs_le.mp hb).2]

assert_no_sorry weakSlope_spread
assert_no_sorry endpoint_residual
#print axioms endpoint_residual

end OPH.EinsteinBranch.AreaReadoutSpread
