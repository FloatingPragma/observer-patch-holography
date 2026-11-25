import Mathlib.Algebra.Group.Basic
import Mathlib.Tactic.Abel
import Mathlib.Util.AssertNoSorry

/-! Disconnected support cancellation for mixed interaction differences. -/

namespace OPH.EinsteinBranch.MeanForceBoundary

variable {G A : Type*} [AddCommMonoid G] [AddCommGroup A]

/-- Apply one finite difference for each interaction increment, including repeated ones. -/
def mixedDifference : List G → (G → A) → G → A
  | [], f, x => f x
  | h :: rest, f, x => mixedDifference rest f (x + h) - mixedDifference rest f x

theorem mixedDifference_add (steps : List G) (f g : G → A) (x : G) :
    mixedDifference steps (fun y => f y + g y) x =
      mixedDifference steps f x + mixedDifference steps g x := by
  induction steps generalizing x with
  | nil => rfl
  | cons a rest ih =>
    simp only [mixedDifference, ih]
    abel

theorem mixedDifference_translation (steps : List G) (f : G → A) (h : G)
    (hf : ∀ x, f (x + h) = f x) (x : G) :
    mixedDifference steps f (x + h) = mixedDifference steps f x := by
  induction steps generalizing x with
  | nil => exact hf x
  | cons a rest ih =>
    simp only [mixedDifference, add_right_comm x h a, ih]

theorem mixedDifference_inactive (steps : List G) (f : G → A) (h : G)
    (hh : h ∈ steps) (hf : ∀ x, f (x + h) = f x) (x : G) :
    mixedDifference steps f x = 0 := by
  induction steps generalizing x with
  | nil => simp at hh
  | cons a rest ih =>
    rcases List.mem_cons.mp hh with heq | hrest
    · subst a
      exact sub_eq_zero.mpr (mixedDifference_translation rest f h hf x)
    · simp only [mixedDifference, ih hrest, sub_self]

/-- The log of a tensor product is additive after embedding the two factors.
    For an additive split, any mixed difference using an inactive direction of
    each summand vanishes, even when the steps include repetitions. -/
theorem disconnected_cancellation (steps : List G) (f g : G → A) (u v x : G)
    (hu : u ∈ steps) (hv : v ∈ steps)
    (hf : ∀ y, f (y + u) = f y) (hg : ∀ y, g (y + v) = g y) :
    mixedDifference steps (fun y => f y + g y) x = 0 := by
  rw [mixedDifference_add, mixedDifference_inactive steps f u hu hf,
    mixedDifference_inactive steps g v hv hg, add_zero]

assert_no_sorry disconnected_cancellation
#print axioms disconnected_cancellation

end OPH.EinsteinBranch.MeanForceBoundary
