import Mathlib.LinearAlgebra.Basis.VectorSpace
import Mathlib.Tactic

/-!
# Regional reconstruction for a canonical scalar split

The local field at two consecutive layers, together with exterior linear
readouts C, determines every local velocity exactly iff ker C ⊆ ker B,
where B is the full exterior-to-region force block. No exterior coupling
is discarded. The theorem applies to the mass-normalized velocity form of
the source scalar action; multiplication by the positive diagonal mass
then recovers canonical momentum.

These are vector-space statements. The Weyl-algebra interpretation and its
Schroedinger representation argument are analytic in the companion paper.
Classical configuration records are not assumed to be quantum observables
or to attest quantum measurement outcomes. Action, regions and access to
the Weyl families are supplied.
-/

set_option autoImplicit false

namespace OPH.QFT.ScalarRegionalTimeSlice

noncomputable section

variable {V E Z : Type*}
  [AddCommGroup V] [Module ℝ V]
  [AddCommGroup E] [Module ℝ E]
  [AddCommGroup Z] [Module ℝ Z]

/-- Position component of the kick-drift-kick canonical split. -/
def nextField (τ : ℝ) (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V)
    (q p : V) (e : E) : V :=
  q + τ • p - (τ ^ 2 / 2) • (L q + B e)

/-- Local initial position, next position and the admitted exterior readouts
separate local velocities. This does not postulate a recovery formula. -/
def DeterminesVelocity (τ : ℝ) (L : V →ₗ[ℝ] V)
    (B : E →ₗ[ℝ] V) (C : E →ₗ[ℝ] Z) : Prop :=
  ∀ q p p' e e', C e = C e' →
    nextField τ L B q p e = nextField τ L B q p' e' → p = p'

/-- Supplying the actual exterior force gives the exact reconstruction. -/
theorem reconstruct_velocity (τ : ℝ) (hτ : τ ≠ 0)
    (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V) (q p : V) (e : E) :
    τ⁻¹ • (nextField τ L B q p e - q) +
      (τ / 2) • (L q + B e) = p := by
  unfold nextField
  have h1 : τ⁻¹ * τ = 1 := inv_mul_cancel₀ hτ
  have h2 : τ⁻¹ * (τ ^ 2 / 2) = τ / 2 := by field_simp
  simp only [smul_sub, smul_add, smul_smul]
  rw [h1, h2, one_smul]
  abel

/-- An invisible exterior perturbation can be cancelled at the next layer
by a nonzero local velocity perturbation. -/
theorem hidden_exterior_perturbation (τ : ℝ)
    (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V) (e : E) :
    nextField τ L B 0 ((τ / 2) • B e) e =
      nextField τ L B 0 0 0 := by
  unfold nextField
  simp only [map_zero, add_zero, zero_add, smul_zero, sub_zero, smul_smul]
  have h : τ * (τ / 2) = τ ^ 2 / 2 := by ring
  rw [h, sub_self]

/-- The necessary and sufficient regional observability criterion. -/
theorem determines_velocity_iff_kernel (τ : ℝ) (hτ : τ ≠ 0)
    (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V) (C : E →ₗ[ℝ] Z) :
    DeterminesVelocity τ L B C ↔ LinearMap.ker C ≤ LinearMap.ker B := by
  constructor
  · intro h e he
    have hCe : C e = C 0 := by simpa using he
    have hp : (τ / 2) • B e = 0 :=
      h 0 ((τ / 2) • B e) 0 e 0 hCe (hidden_exterior_perturbation τ L B e)
    have hn : τ / 2 ≠ 0 := div_ne_zero hτ (by norm_num)
    exact (smul_eq_zero.mp hp).resolve_left hn
  · intro h q p p' e e' hC hnext
    have hCe : e - e' ∈ LinearMap.ker C := by
      simp only [LinearMap.mem_ker, map_sub, hC, sub_self]
    have hBe : B e = B e' := by
      have hb := h hCe
      simpa only [LinearMap.mem_ker, map_sub, sub_eq_zero] using hb
    have hrec := reconstruct_velocity τ hτ L B q p e
    have hrec' := reconstruct_velocity τ hτ L B q p' e'
    rw [hnext, hBe] at hrec
    exact hrec.symm.trans hrec'

/-- The actual force is an admissible and sufficient exterior statistic. -/
theorem force_readout_suffices (τ : ℝ) (hτ : τ ≠ 0)
    (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V) :
    DeterminesVelocity τ L B B :=
  (determines_velocity_iff_kernel τ hτ L B B).2 le_rfl

/-- With no exterior access, exact recovery is possible precisely when the
full exterior coupling vanishes. This is not a no-go for collar readouts. -/
theorem no_collar_iff_zero_coupling (τ : ℝ) (hτ : τ ≠ 0)
    (L : V →ₗ[ℝ] V) (B : E →ₗ[ℝ] V) :
    DeterminesVelocity τ L B (0 : E →ₗ[ℝ] Z) ↔ B = 0 := by
  rw [determines_velocity_iff_kernel τ hτ L B]
  constructor
  · intro h
    ext e
    exact h (by simp)
  · intro h
    simp [h]

#print axioms reconstruct_velocity
#print axioms hidden_exterior_perturbation
#print axioms determines_velocity_iff_kernel
#print axioms force_readout_suffices
#print axioms no_collar_iff_zero_coupling

end

end OPH.QFT.ScalarRegionalTimeSlice
