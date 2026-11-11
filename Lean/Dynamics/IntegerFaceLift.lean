import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.BigOperators.Group.Finset.Sigma
import Mathlib.Data.Int.ModEq
import Mathlib.Tactic.NormNum

/-! # Integer face lifts on a closed incidence carrier

The incidence and update law are declared finite data. The executable
attachment checks them against the actual oriented carrier and accepted
native twists. No hypercharge or particle identity follows from these
integer identities.
-/

namespace OPH.Dynamics.IntegerFaceLift

variable {Face Seam : Type*} [Fintype Face] [Fintype Seam]

def lift (incidence : Face → Seam → ℤ) (link : Seam → ℤ) (f : Face) : ℤ :=
  ∑ s, incidence f s * link s

theorem total_zero (incidence : Face → Seam → ℤ)
    (closed : ∀ s, ∑ f, incidence f s = 0) (link : Seam → ℤ) :
    ∑ f, lift incidence link f = 0 := by
  simp only [lift]
  rw [Finset.sum_comm]
  simp only [← Finset.sum_mul, closed, zero_mul, Finset.sum_const_zero]

omit [Fintype Face] in
theorem lift_add (incidence : Face → Seam → ℤ) (link delta : Seam → ℤ) :
    lift incidence (fun s => link s + delta s) =
      fun f => lift incidence link f + lift incidence delta f := by
  funext f
  simp [lift, mul_add, Finset.sum_add_distrib]

omit [Fintype Face] in
theorem one_twist [DecidableEq Seam] (incidence : Face → Seam → ℤ)
    (link : Seam → ℤ) (seam : Seam) (current : ℤ) (f : Face) :
    lift incidence (fun s => link s + if s = seam then current else 0) f =
      lift incidence link f + incidence f seam * current := by
  simp [lift, mul_add, Finset.sum_add_distrib, mul_ite]

theorem modulo_update (value increment klass : ℤ)
    (hmod : value % 6 = klass) :
    (value + increment) % 6 = (klass + increment % 6) % 6 := by
  rw [Int.add_emod, hmod]

theorem modulo_representative (value winding : ℤ) :
    (value + 6 * winding) % 6 = value % 6 := by
  simp [Int.add_emod]

theorem equal_class_need_not_equal_value :
    (1 : ℤ) % 6 = 7 % 6 ∧ (1 : ℤ) ≠ 7 := by norm_num

end OPH.Dynamics.IntegerFaceLift
