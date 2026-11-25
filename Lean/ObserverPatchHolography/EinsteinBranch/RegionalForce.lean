import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Util.AssertNoSorry

/-! A finite signed test for reassigning a measured body force to a
time derivative at fixed measured boundary flux. -/
namespace OPH.EinsteinBranch.RegionalForce

variable {ι R : Type*} [Fintype ι] [CommRing R]

/-- Any weight commuting with the generator has zero trace pairing with
a Heisenberg commutator. A scalar such as i can be absorbed into K. -/
theorem stationary_trace (W H K : Matrix ι ι R) (h : W * H = H * W) :
    Matrix.trace (W * (H * K - K * H)) = 0 := by
  rw [mul_sub, Matrix.trace_sub, ← mul_assoc, h,
    Matrix.trace_mul_cycle' W K H, ← mul_assoc, sub_self]

/-- A nonzero measured moment rejects every possible additive dressing
of the finite momentum, including those with nonlocal support. -/
theorem no_commutator [DecidableEq ι] (H F : Matrix ι ι R) (n : ℕ)
    (hf : Matrix.trace (H ^ n * F) ≠ 0) :
    ¬ ∃ K : Matrix ι ι R, F = H * K - K * H := by
  rintro ⟨K, rfl⟩
  exact hf (stationary_trace (H ^ n) H K (by rw [← pow_succ, ← pow_succ']))

assert_no_sorry stationary_trace
assert_no_sorry no_commutator

end OPH.EinsteinBranch.RegionalForce
