import Mathlib
import Mathlib.Util.AssertNoSorry

/-! Two readouts of one elapsed-time law on the five-port material clock.
The input is exact agreement of both mean phase deficits with zero. -/
namespace OPH.EinsteinBranch.TwoFrequencyClock

open MeasureTheory

/-- The ratio of the two gaps of the five-cycle is (3 + √5)/2. -/
theorem cycle_ratio_irrational :
    Irrational (((5 + Real.sqrt 5) / 4) / ((5 - Real.sqrt 5) / 4)) := by
  have hs : (Real.sqrt 5) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  have hlo : Real.sqrt 5 < 5 := by nlinarith [Real.sqrt_nonneg 5]
  have hi : Irrational (Real.sqrt 5) := by norm_num
  rw [irrational_iff_ne_rational] at hi ⊢
  intro a b hb heq
  have hb' : (b : ℝ) ≠ 0 := by exact_mod_cast hb
  have hratio : ((5 + Real.sqrt 5) / 4) / ((5 - Real.sqrt 5) / 4)
      = (3 + Real.sqrt 5) / 2 := by
    apply (div_eq_iff (by linarith)).2
    nlinarith
  rw [hratio, div_eq_div_iff (by norm_num) hb'] at heq
  apply hi (2*a-3*b) b hb
  apply (eq_div_iff hb').2
  push_cast
  linarith

/-- Exact cosine agreement fixes each elapsed time in the support of both readings. -/
theorem elapsed_eq_of_phases {low high t x : ℝ}
    (hlow : low ≠ 0) (hirr : Irrational (high / low))
    (ha : Real.cos (low * (x-t)) = 1) (hb : Real.cos (high * (x-t)) = 1) :
    x = t := by
  obtain ⟨a, ha⟩ := (Real.cos_eq_one_iff _).1 ha
  obtain ⟨b, hb⟩ := (Real.cos_eq_one_iff _).1 hb
  by_contra h
  have hd : x-t ≠ 0 := sub_ne_zero.mpr h
  have hzero : a ≠ 0 := by
    intro h0
    subst a
    simp [hlow, hd] at ha
  apply (irrational_iff_ne_rational _).1 hirr b a hzero
  apply (div_eq_div_iff hlow (by exact_mod_cast hzero)).2
  apply mul_right_cancel₀ hd
  calc
    high * (a : ℝ) * (x-t) = (a : ℝ) * (high * (x-t)) := by ring
    _ = (b : ℝ) * ((a : ℝ) * (2 * Real.pi)) := by rw [← hb]; ring
    _ = (b : ℝ) * low * (x-t) := by rw [ha]; ring

private theorem phase_one_ae (μ : Measure ℝ) (g t : ℝ)
    (hi : Integrable (fun x : ℝ => 1 - Real.cos (g*(x-t))) μ)
    (hz : (∫ x, 1 - Real.cos (g*(x-t)) ∂μ) = 0) :
    ∀ᵐ x ∂μ, Real.cos (g*(x-t)) = 1 := by
  have hn : 0 ≤ (fun x : ℝ => 1 - Real.cos (g*(x-t))) :=
    fun x => sub_nonneg.mpr (Real.cos_le_one _)
  filter_upwards [(integral_eq_zero_iff_of_nonneg hn hi).mp hz] with x hx
  change 1 - Real.cos (g*(x-t)) = 0 at hx
  linarith

/-- A common random time with both exact ordinary phases equals the appointment
almost surely. The statement permits unbounded support and supplies no
finite-precision exclusion for arbitrarily late near-recurrences. -/
theorem common_random_time (μ : Measure ℝ) (low high t : ℝ)
    (hl : low ≠ 0) (hr : Irrational (high / low))
    (hil : Integrable (fun x : ℝ => 1 - Real.cos (low*(x-t))) μ)
    (hih : Integrable (fun x : ℝ => 1 - Real.cos (high*(x-t))) μ)
    (hzl : (∫ x, 1 - Real.cos (low*(x-t)) ∂μ) = 0)
    (hzh : (∫ x, 1 - Real.cos (high*(x-t)) ∂μ) = 0) :
    ∀ᵐ x ∂μ, x = t := by
  filter_upwards [phase_one_ae μ low t hil hzl, phase_one_ae μ high t hih hzh] with x hlx hhx
  exact elapsed_eq_of_phases hl hr hlx hhx

/-- Specialization to the two measured modes of the five-port clock. -/
theorem five_port_random_time (μ : Measure ℝ) (t : ℝ)
    (hint : ∀ g ∈ ({(5 - Real.sqrt 5)/4, (5 + Real.sqrt 5)/4} : Finset ℝ),
      Integrable (fun x : ℝ => 1 - Real.cos (g*(x-t))) μ)
    (hzero : ∀ g ∈ ({(5 - Real.sqrt 5)/4, (5 + Real.sqrt 5)/4} : Finset ℝ),
      (∫ x, 1 - Real.cos (g*(x-t)) ∂μ) = 0) : ∀ᵐ x ∂μ, x = t := by
  have hs := Real.sq_sqrt (show (0:ℝ) ≤ 5 by norm_num)
  exact common_random_time μ _ _ t (by nlinarith [Real.sqrt_nonneg 5])
    cycle_ratio_irrational (hint _ (by simp)) (hint _ (by simp))
    (hzero _ (by simp)) (hzero _ (by simp))

assert_no_sorry five_port_random_time
assert_no_sorry cycle_ratio_irrational
assert_no_sorry common_random_time
end OPH.EinsteinBranch.TwoFrequencyClock
