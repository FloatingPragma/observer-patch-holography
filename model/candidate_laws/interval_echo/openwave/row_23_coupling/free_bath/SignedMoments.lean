import Mathlib.Tactic.NoncommRing

namespace OPH.FreeBathSignedMoments

variable {R : Type*} [Ring R]

/-- No contact meets both signs when their guarded sites are disjoint and nonadjacent. -/
theorem disjoint_mean (h a b : R) :
    (h+a+b) - (h+a-b) - (h-a+b) + (h-a-b) = 0 := by
  noncomm_ring

/-- The second signed spectral moment can survive even when the first vanishes.
No commutation assumption on propagation or the two source terms is used here. -/
theorem disjoint_second (h a b : R) :
    (h+a+b)^2 - (h+a-b)^2 - (h-a+b)^2 + (h-a-b)^2 =
    4 * (a*b+b*a) := by
  noncomm_ring

#print axioms disjoint_mean
#print axioms disjoint_second
end OPH.FreeBathSignedMoments
