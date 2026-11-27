import Mathlib

/-! A detector orthogonal to a stationary difference of two source words reads a
single bright temporal profile under any diagonal phase preparation. -/

namespace PinnedPair

variable {V : Type*} [AddCommGroup V] [Module ℂ V]

/-- After free evolution, any linear receiver that kills the protected difference
sees only the sum of the two source coefficients. -/
theorem received_amplitude (F : V →ₗ[ℂ] ℂ) (x z : V) (a b : ℂ)
    (hd : F (x - z) = 0) :
    F (a • x + b • z) = (a + b) * F x := by
  have h : F x = F z := sub_eq_zero.mp (by simpa using hd)
  simp only [map_add, map_smul, smul_eq_mul]
  rw [← h, add_mul]

/-- Counts at two waits satisfy the rank-one identity for any two phase kicks.
Both receivers must reject the same protected difference. -/
theorem count_profile (F G : V →ₗ[ℂ] ℂ) (x z : V) (a b c d : ℂ)
    (hf : F (x - z) = 0) (hg : G (x - z) = 0) :
    Complex.normSq (F (a • x + b • z)) * Complex.normSq (G (c • x + d • z)) =
    Complex.normSq (F (c • x + d • z)) * Complex.normSq (G (a • x + b • z)) := by
  rw [received_amplitude F x z a b hf, received_amplitude G x z c d hg,
      received_amplitude F x z c d hf, received_amplitude G x z a b hg]
  simp only [Complex.normSq_mul]
  ring

end PinnedPair
