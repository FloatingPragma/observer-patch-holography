import QFT.LatticeMomentumBalance
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# Quantitative coarse Ward defects

The concrete local CAR identity and the collision-suppressed initial-state
bound are supplied by `model/dilute_momentum_ward`. These finite-sum results
show that spatial coarse-graining sums the certified local error budgets.
They do not turn an approximate Ward column into exact Lorentz stress.
-/

namespace OPH.QFT.DiluteMomentumWard

open OPH.QFT.CoarseEnergyWard

variable {X Y : Type*} [Fintype X] [Fintype Y] [DecidableEq Y]

omit [Fintype Y] in
theorem coarse_error_bound (c : X → Y) (f e : X → ℝ)
    (bound : ∀ x, |f x| ≤ e x) (a : Y) :
    |push c f a| ≤ push c e a := by
  unfold push
  apply (Finset.abs_sum_le_sum_abs _ _).trans
  apply Finset.sum_le_sum
  intro x _
  by_cases h : c x = a
  · simp only [h, ite_true]
    exact bound x
  · simp [h]

theorem coarse_Ward_defect_bound (c : X → Y) (delta force error : X → ℝ)
    (J : X → X → ℝ)
    (balance : ∀ x, delta x + ∑ y, J x y = force x)
    (bound : ∀ x, |force x| ≤ error x) (a : Y) :
    |push c delta a + ∑ b, coarse c J a b| ≤ push c error a := by
  rw [OPH.QFT.LatticeMomentumBalance.coarse_forced_balance c delta force J balance a]
  exact coarse_error_bound c force error bound a

theorem merged_error_bound (f e : X → ℝ) (bound : ∀ x, |f x| ≤ e x) :
    |∑ x, f x| ≤ ∑ x, e x := by
  exact (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun x _ => bound x)

#print axioms coarse_error_bound
#print axioms coarse_Ward_defect_bound
#print axioms merged_error_bound

end OPH.QFT.DiluteMomentumWard
