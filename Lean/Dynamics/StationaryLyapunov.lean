import Dynamics.ProtectedCharge

/-! # A stationary finite law cannot have a strictly negative mean drift

This result reuses the finite record push/pull duality. The source must
separately supply stochasticity and its pointwise drift inequalities.
It selects no thermal energy, bath or physical clock.
-/

namespace OPH.Dynamics

open scoped BigOperators

variable {State : Type*} [Fintype State]

theorem stationary_mean_drift_zero (kernel : State → State → ℝ)
    (p energy : State → ℝ) (hp : kernelPush kernel p = p) :
    (∑ x, p x * ((kernelPull kernel energy) x - energy x)) = 0 := by
  have h := kernelPush_kernelPull_areExpectationDual kernel p energy
  rw [hp] at h
  simp only [chargeExpectation] at h
  simp_rw [mul_sub]
  rw [Finset.sum_sub_distrib, h]
  exact sub_self _

theorem no_full_support_stationary_of_strict_drift
    (kernel : State → State → ℝ) (p energy : State → ℝ)
    (positive : ∀ x, 0 < p x)
    (nonpositive : ∀ x, (kernelPull kernel energy) x - energy x ≤ 0)
    (strict : ∃ x, (kernelPull kernel energy) x - energy x < 0) :
    kernelPush kernel p ≠ p := by
  intro hp
  have hzero := stationary_mean_drift_zero kernel p energy hp
  have hn : ∀ x ∈ (Finset.univ : Finset State),
      p x * ((kernelPull kernel energy) x - energy x) ≤ 0 := by
    intro x _
    exact mul_nonpos_of_nonneg_of_nonpos (le_of_lt (positive x)) (nonpositive x)
  have hall := (Finset.sum_eq_zero_iff_of_nonpos hn).mp hzero
  obtain ⟨x, hx⟩ := strict
  have hz := hall x (Finset.mem_univ x)
  have hlt := mul_neg_of_pos_of_neg (positive x) hx
  linarith

end OPH.Dynamics

#print axioms OPH.Dynamics.stationary_mean_drift_zero
#print axioms OPH.Dynamics.no_full_support_stationary_of_strict_drift
