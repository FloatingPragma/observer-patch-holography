import Mathlib

/-! Quadratic inequalities for the spectrum of a single vacancy on an
icosahedral federation whose contacts are blank or carry native zero-contrast
marks. Real and imaginary amplitudes obey the same inequalities separately.
The executable checks the native marked-seam columns and the incidence and
private terms on the physical carrier. -/
namespace SourceLaw.SilentSeamVacancy

/-- A free edge pays a lower bound at its leaf by charging the neighboring
vertex. The completed square makes the estimate sharp for one amplitude ratio. -/
theorem leaf_cost {a : ℝ} (ha : 0 < 1 - 2*a) (u v : ℝ) :
    a*u^2 ≤ (u-v)^2/2 + a/(1-2*a)*v^2 := by
  have hid : (1-2*a)*((u-v)^2 + 2*(a/(1-2*a))*v^2 - 2*a*u^2)
      = ((1-2*a)*u-v)^2 := by
    field_simp [ne_of_gt ha]
    ring
  have hn : 0 ≤ (1-2*a)*((u-v)^2 + 2*(a/(1-2*a))*v^2 - 2*a*u^2) := by
    rw [hid]
    positivity
  by_contra bad
  push Not at bad
  have hneg : (u-v)^2 + 2*(a/(1-2*a))*v^2 - 2*a*u^2 < 0 := by linarith
  exact (not_lt_of_ge hn) (mul_neg_of_pos_of_neg ha hneg)

/-- At each internal vertex, the native potential pays the lower bound
and the cost of all its incident leaf edges. Equality is possible at degree two. -/
theorem internal_budget {a b d : ℝ} (hd : 2 ≤ d) (hb : b ≤ 8)
    (hab : 8-2*b=a) : a + b*d ≤ 8*(d-1) := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hb) (sub_nonneg.mpr hd)]

/-- The physical vacancy floor on nontrivial components is the smaller
root of the three-port even-mode characteristic polynomial. -/
theorem sharp_constants :
    let a := (19 - Real.sqrt 297)/4
    0 < a ∧ 0 < 1-2*a ∧ a/(1-2*a) ≤ 8 ∧ 8-2*(a/(1-2*a))=a := by
  have hs := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 297)
  have hp := Real.sqrt_nonneg 297
  have hlo : (17:ℝ) < Real.sqrt 297 := by nlinarith
  have hhi : Real.sqrt 297 < (18:ℝ) := by nlinarith
  dsimp
  have hd : 0 < 1 - 2 * ((19 - Real.sqrt 297)/4) := by linarith
  refine ⟨by linarith, hd, ?_, ?_⟩
  · apply (div_le_iff₀ hd).2
    nlinarith [sq_nonneg (Real.sqrt 297 - 17), sq_nonneg (Real.sqrt 297 - 18)]
  · let a := (19 - Real.sqrt 297)/4
    have hquad : 2*a^2 - 19*a + 8 = 0 := by dsimp [a]; nlinarith
    have hratio : a/(1-2*a) = (8-a)/2 := by
      apply (div_eq_iff (ne_of_gt hd)).2
      change a = ((8-a)/2) * (1-2*a)
      nlinarith [hquad]
    change 8 - 2 * (a/(1-2*a)) = a
    linarith

/-- A finite collection of leaf bounds and core budgets can be added with
equal and opposite payment. The application partitions the quadratic form
by leaves and vertices of degree at least two. -/
theorem add_payments {l c paid energyLeaf energyCore : ℝ}
    (hl : l ≤ energyLeaf + paid) (hc : c + paid ≤ energyCore) :
    l + c ≤ energyLeaf + energyCore := by linarith

#print axioms leaf_cost
#print axioms internal_budget
#print axioms sharp_constants
#print axioms add_payments

end SourceLaw.SilentSeamVacancy
