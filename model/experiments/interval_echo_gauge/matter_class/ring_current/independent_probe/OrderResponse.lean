import Mathlib.Tactic.NoncommRing

namespace OPH.NativeWindingProbe

variable {A : Type*} [Ring A]

/-- Alternating signs in two independent encounters, using unnormalized returns. -/
def fourSigns (u x v y : A) : A :=
  (u + x) * (v + y) - (u + x) * (v - y) -
  (u - x) * (v + y) + (u - x) * (v - y)

/-- The sign-independent parts of both encounters disappear exactly. -/
theorem isolate (u x v y : A) : fourSigns u x v y = 4 * (x * y) := by
  unfold fourSigns
  noncomm_ring

/-- Reversing actual visit order extracts the bracket on signed returns. -/
theorem orderedDifference (u x v y : A) :
    fourSigns u x v y - fourSigns v y u x = 4 * (x * y - y * x) := by
  rw [isolate, isolate]
  noncomm_ring

end OPH.NativeWindingProbe
