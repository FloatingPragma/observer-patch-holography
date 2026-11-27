import Mathlib.Algebra.GroupWithZero.Defs
import Mathlib.Algebra.Group.Basic

namespace IntervalEcho.ChargedFork

/-- A core-only unitary leaves every remote record law invariant if its pulled-back
effect commutes with the pulse. This includes a correlated initial preparation. -/
theorem spectator (M : Type*) [Monoid M] (pulse undo effect : M)
    (inverse : undo * pulse = 1) (independent : effect * pulse = pulse * effect) :
    undo * effect * pulse = effect := by
  rw [mul_assoc, independent, ← mul_assoc, inverse, one_mul]

/-- At a fixed marked vacancy address, every sequence of conserving receiver
interventions has zero transfer into a different terminal address. -/
theorem parked (M : Type*) [MonoidWithZero M] (src dst route : M)
    (disjoint : dst * src = 0) (fixed : route * src = src * route) :
    dst * route * src = 0 := by
  rw [mul_assoc, fixed, ← mul_assoc, disjoint, zero_mul]

end IntervalEcho.ChargedFork
