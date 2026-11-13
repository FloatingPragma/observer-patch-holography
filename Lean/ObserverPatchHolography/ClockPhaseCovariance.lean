import Mathlib

/-! Algebraic conditional-state invariance for energy-shift instruments.
The result assumes the stated native commutator identities. It supplies no
physical source identification or infinite-volume symmetry-breaking limit. -/
namespace OPH.ClockPhaseCovariance

/-- A single energy-shift branch preserves an invariant density operator.
`L` is the adjoint branch in the operator interpretation; the theorem uses
its explicit opposite-shift identity and therefore applies in any ring. -/
theorem invariant_branch {R : Type*} [Ring R] (H K L rho shift : R)
    (hK : H * K = K * H + K * shift)
    (hL : L * H = H * L + shift * L)
    (hrho : H * rho = rho * H)
    (hshift : shift * rho = rho * shift) :
    H * (K * rho * L) = (K * rho * L) * H := by
  calc
    H * (K * rho * L) = (H * K) * rho * L := by simp only [mul_assoc]
    _ = (K * H + K * shift) * rho * L := by rw [hK]
    _ = K * (H * rho) * L + K * (shift * rho) * L := by noncomm_ring
    _ = K * (rho * H) * L + K * (rho * shift) * L := by rw [hrho, hshift]
    _ = K * rho * (H * L + shift * L) := by noncomm_ring
    _ = K * rho * (L * H) := by rw [hL]
    _ = (K * rho * L) * H := by simp only [mul_assoc]

/-- Invariance survives complete finite coarse graining of invariant branches. -/
theorem invariant_sum {R : Type*} [Ring R] {I : Type*} (S : Finset I)
    (H : R) (branch : I → R)
    (h : ∀ i ∈ S, H * branch i = branch i * H) :
    H * (∑ i ∈ S, branch i) = (∑ i ∈ S, branch i) * H := by
  simp only [Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun i hi => h i hi

#print axioms invariant_branch
#print axioms invariant_sum
end OPH.ClockPhaseCovariance
