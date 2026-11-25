import Mathlib

namespace SourceLaw.ObserverTransfer

/-- The energy on a terminal filled seam pulls back to its conserved head
observable on the initial state. F projects onto two filled endpoints. -/
theorem filled_read {R : Type*} [Semigroup R] (h f b u : R)
    (filled : h * f = b * f) (locality : b * f = f * b)
    (flight : b * u = u * b) : h * f * u = f * u * b := by
  rw [filled, locality, mul_assoc, flight, ← mul_assoc]

/-- A parity-even preparation gives zero amplitude in a parity-odd detector
whenever flight commutes with the parity operation. -/
theorem opposite_parity {V : Type*} [AddCommGroup V] [Module ℂ V]
    (c u : V →ₗ[ℂ] V) (read : V →ₗ[ℂ] ℂ) (v : V)
    (source : c v = v) (flight : ∀ x, u (c x) = c (u x))
    (detector : ∀ x, read (c x) = -read x) : read (u v) = 0 := by
  have fixed : c (u v) = u v := by rw [← flight, source]
  have zero := detector (u v)
  rw [fixed] at zero
  linear_combination (1 / 2 : ℂ) * zero

end SourceLaw.ObserverTransfer
