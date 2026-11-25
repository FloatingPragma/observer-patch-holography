import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall

/-! Integrated force balance and the mass measured by a static clock gradient
at first order. Coordinate weights, force and geometry have independent meanings. -/
namespace OPH.EinsteinBranch.StaticMassFlux

open scoped BigOperators

/-- Finite summation by parts. The incidence matrix can include oriented
boundary edges, provided the same source and boundary convention is used. -/
theorem coordinate_force {V E : Type*} [Fintype V] [Fintype E]
    (x : V → ℝ) (inc : V → E → ℝ) (f : E → ℝ) :
    (∑ v, x v * ∑ e, inc v e * f e) =
    ∑ e, (∑ v, x v * inc v e) * f e := by
  simp only [Finset.mul_sum, Finset.sum_mul, mul_assoc]
  exact Finset.sum_comm

/-- A force-conserving finite network has zero total stress moment in
each tested coordinate when all supporting forces are included. -/
theorem static_virial {V E : Type*} [Fintype V] [Fintype E]
    (x : V → ℝ) (inc : V → E → ℝ) (f : E → ℝ)
    (balance : ∀ v, ∑ e, inc v e * f e = 0) :
    ∑ e, (∑ v, x v * inc v e) * f e = 0 := by
  rw [← coordinate_force]
  simp [balance]

/-- Integrated trace reversal with a possible measured boundary moment.
The inputs can be differences of stationary, weak-field preparations. -/
theorem mass_flux (flux energy pressure boundary ward field κ : ℝ)
    (virial : pressure = boundary - ward)
    (poisson : flux = κ*(energy + pressure)/2 + field) :
    flux - κ*energy/2 = κ*(boundary - ward)/2 + field := by
  rw [poisson, virial]
  ring

/-- The exterior average lapse shift between spheres, once the pressure
moment and the vacuum-subtracted Gauss law have been supplied. -/
theorem shell_average (r s G energy measured : ℝ)
    (read : measured = G*energy*(1/r - 1/s))
    (hr : 0 < r) (hs : r < s) (hG : 0 < G) (hm : 0 < energy) :
    0 < measured := by
  rw [read]
  exact mul_pos (mul_pos hG hm) (sub_pos.mpr (one_div_lt_one_div_of_lt hr hs))

/-- If one clock response fits both exterior shell averages to accuracy ε,
the integrated mass signal has magnitude at most 2ε, including coverage errors. -/
theorem constant_clock_bound (inner outer clock signal ε : ℝ)
    (hi : |inner - clock| ≤ ε) (ho : |outer - clock| ≤ ε)
    (h : outer - inner = signal) : |signal| ≤ 2*ε := by
  obtain ⟨a,b⟩ := abs_le.mp hi
  obtain ⟨c,d⟩ := abs_le.mp ho
  apply abs_le.mpr
  constructor <;> linarith

assert_no_sorry constant_clock_bound
assert_no_sorry coordinate_force
assert_no_sorry static_virial
assert_no_sorry mass_flux
assert_no_sorry shell_average

end OPH.EinsteinBranch.StaticMassFlux
