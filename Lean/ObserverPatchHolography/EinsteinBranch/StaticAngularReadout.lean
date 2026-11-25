import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall

/-! Directional shell and clock readings on a static region. The analytic
small-ball expansions and physical ruler calibration have separate hypotheses. -/
namespace OPH.EinsteinBranch.StaticAngularReadout

open scoped BigOperators

noncomputable def meanSix (q : Fin 6 → ℝ) : ℝ := (∑ i, q i) / 6

noncomputable def aperture (q : Fin 6 → ℝ) (i : Fin 6) : ℝ := (6 * meanSix q + q i) / 7

theorem mean_aperture (q : Fin 6 → ℝ) : meanSix (aperture q) = meanSix q := by
  unfold aperture meanSix
  simp only [div_eq_mul_inv, ← Finset.sum_mul, Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  ring

theorem decode (q : Fin 6 → ℝ) (i : Fin 6) :
    7 * aperture q i - 6 * meanSix (aperture q) = q i := by
  rw [mean_aperture]
  unfold aperture
  ring

/-- The dimension-three trace is three times the mean of the six normalized
directional Hessian readings. Curvature and matter may be compared separately. -/
theorem static_null (spatial lapse : Fin 6 → ℝ) (lap : ℝ)
    (traceRead : lap = 3 * meanSix lapse) (i : Fin 6) :
    (7 * aperture spatial i - 6 * meanSix (aperture spatial)) + lap -
      (7 * aperture lapse i - 6 * meanSix (aperture lapse)) =
      spatial i + 3 * meanSix lapse - lapse i := by
  rw [decode, decode, traceRead]

assert_no_sorry mean_aperture
assert_no_sorry decode
assert_no_sorry static_null

end OPH.EinsteinBranch.StaticAngularReadout
