import Mathlib
import Mathlib.Util.AssertNoSorry

/-! Observable inequalities used in the SPARC shape, cosmic budget and cubic
strength calculations. The observational intervals and the three-dimensional
current construction are supplied separately by their respective experiments. -/
namespace OPH.EinsteinBranch.GravityCalibration

/-- A signed combination of upper bounds with nonpositive weights certifies
a lower bound. Exact identification at a feasible point proves optimality. -/
theorem certificate_lower_bound {ι : Type*} (s : Finset ι) (w lhs rhs : ι → ℝ)
    (hw : ∀ i ∈ s, w i ≤ 0) (hrows : ∀ i ∈ s, lhs i ≤ rhs i) :
    ∑ i ∈ s, w i * rhs i ≤ ∑ i ∈ s, w i * lhs i :=
  Finset.sum_le_sum fun i hi => mul_le_mul_of_nonpos_left (hrows i hi) (hw i hi)

/-- If each tracer's outer squared-speed term is at most q times its inner
term, no nonnegative mixture can rise faster. The component terms may be signed. -/
theorem tracer_shape_bound {ι : Type*} (s : Finset ι) (w inner outer : ι → ℝ) (q : ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (hrow : ∀ i ∈ s, outer i ≤ q * inner i) :
    ∑ i ∈ s, w i * outer i ≤ q * ∑ i ∈ s, w i * inner i := by
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun i hi => by
    simpa only [mul_left_comm] using mul_le_mul_of_nonneg_left (hrow i hi) (hw i hi)

/-- The outer excess fraction uses an inward residual at the inner ring and a
positive outer squared speed. Inputs L and U bound the two observed squared speeds. -/
theorem outer_excess_fraction {q inner outer accountedInner accountedOuter L U : ℝ}
    (hq : 0 ≤ q) (hL : 0 < L) (hU : 0 ≤ U) (ho : L ≤ outer) (hi : inner ≤ U)
    (hin : accountedInner ≤ inner) (hshape : accountedOuter ≤ q * accountedInner) :
    1 - q * U / L ≤ (outer - accountedOuter) / outer := by
  have hp : 0 < outer := lt_of_lt_of_le hL ho
  have bound := hshape.trans (mul_le_mul_of_nonneg_left (hin.trans hi) hq)
  have frac : accountedOuter / outer ≤ (q * U) / L :=
    (div_le_div_iff₀ hp hL).2 ((mul_le_mul_of_nonneg_right bound hL.le).trans
      (mul_le_mul_of_nonneg_left ho (mul_nonneg hq hU)))
  rw [sub_div, div_self hp.ne']
  linarith

/-- Finite integration by parts gives the unit flux. Its certified squared
norm yields the lower bound on every grounded center-normalized profile. -/
theorem profile_energy_floor {ι : Type*} (edges : Finset ι) (flow grad : ι → ℝ)
    (hflux : ∑ e ∈ edges, flow e * grad e = 1)
    (hflow : ∑ e ∈ edges, (flow e)^2 ≤ 3/8) :
    8/3 ≤ ∑ e ∈ edges, (grad e)^2 := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq edges flow grad
  rw [hflux, one_pow] at h
  have hg : 0 ≤ ∑ e ∈ edges, (grad e)^2 := Finset.sum_nonneg fun e _ => sq_nonneg _
  linarith [h.trans (mul_le_mul_of_nonneg_right hflow hg)]

/-- With the historical common energy/clock convention, no choice of energy
or length unit changes the ratio between the three physical readings. -/
theorem unit_ratio_invariant (g gap speed u a hbarc G c E : ℝ)
    (hu : u ≠ 0) (ha : a ≠ 0) (hh : hbarc ≠ 0) (hv : speed ≠ 0)
    (hG : G ≠ 0) (hc : c ≠ 0) (hE : E ≠ 0) :
    ((g*a/(u*hbarc))/G) * (gap*u/E)^2 / ((speed*u*a/hbarc)/c)
      = g*gap^2*c/(speed*G*E^2) := by
  field_simp

/-- In the flat, free-ruler branch the exact interval certificate exceeds this
rounded lower dust density. Even the largest allowed baryon and neutrino
terms then require over 1.278 additional dust units per baryon. -/
theorem flat_source_fraction {b ν m f : ℝ}
    (hb : 0 < b) (hbmax : b ≤ 0.02383) (hν : ν ≤ 0.003)
    (hm : 0.05728571 ≤ m) (hstock : m = b*(1+f) + ν) : 1.278 < f := by
  by_contra! hf
  have h₁ := mul_le_mul_of_nonneg_left (show 1+f ≤ 1+1.278 by linarith) hb.le
  have h₂ := mul_le_mul_of_nonneg_right hbmax (show (0:ℝ) ≤ 1+1.278 by norm_num)
  linarith

assert_no_sorry flat_source_fraction
assert_no_sorry certificate_lower_bound
assert_no_sorry tracer_shape_bound
assert_no_sorry outer_excess_fraction
assert_no_sorry profile_energy_floor
assert_no_sorry unit_ratio_invariant
end OPH.EinsteinBranch.GravityCalibration
