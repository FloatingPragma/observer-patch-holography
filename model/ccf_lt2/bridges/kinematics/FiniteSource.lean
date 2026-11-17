import QFT.StructuralInheritance
import Dynamics.PrivateInner
import EventAlgebra.PauliMinimalCoupling
import Lean.Elab.Tactic.Omega

/-!
Finite consequences for the CCF-2 matrix extracted by the accompanying audit.
The executable-to-matrix extraction is an explicit external hypothesis here;
this file does not formalize Python, the physical clock, or a continuum limit.
-/

namespace CCF2.Kinematics

open scoped ComplexOrder Matrix

variable {d : Type*} [Fintype d] [DecidableEq d] [Nonempty d]

/-- The complete finite matrix, including all matter, link and geometry states.
No chosen low-energy block or particle dictionary is part of this structure. -/
structure FiniteSource (d : Type*) [Fintype d] [DecidableEq d] where
  hamiltonian : Matrix d d ℂ
  hermitian : hamiltonian.IsHermitian

/-- Finite Gibbs/KMS is available for the same complete Hamiltonian. -/
theorem finite_source_gibbs (S : FiniteSource d) (β : ℝ) :
    0 < OPH.QFT.gibbsPartition S.hamiltonian β ∧
    OPH.QFT.gibbsState S.hamiltonian β 1 = 1 ∧
    (∀ a : Matrix d d ℂ,
      0 ≤ OPH.QFT.gibbsState S.hamiltonian β (aᴴ * a)) ∧
    (∀ a b : Matrix d d ℂ,
      OPH.QFT.gibbsState S.hamiltonian β
        (a * OPH.QFT.heisenbergFlow S.hamiltonian b (Complex.I * β)) =
      OPH.QFT.gibbsState S.hamiltonian β (b * a)) := by
  exact ⟨OPH.QFT.gibbsPartition_pos S.hermitian β,
    OPH.QFT.gibbsState_one S.hermitian β,
    OPH.QFT.gibbsState_star_mul_self_nonneg S.hermitian β,
    OPH.QFT.gibbsState_kms S.hamiltonian β⟩

/-- A stationary density has no changing internal expectations under its own
Hamiltonian. A reference-assisted relative phase is a different experiment. -/
theorem stationary_density_fixed (H ρ : Matrix d d ℂ)
    (h : Commute H ρ) (z : ℂ) : OPH.QFT.heisenbergFlow H ρ z = ρ := by
  have hc : Commute (NormedSpace.exp (z • (Complex.I • H))) ρ :=
    ((h.smul_left Complex.I).smul_left z).exp_left
  unfold OPH.QFT.heisenbergFlow OPH.QFT.matrixInnerFlow
  rw [smul_neg, hc.eq, mul_assoc, OPH.QFT.matrix_exp_mul_exp_neg, mul_one]

/-- A sharp charge in the declared integer-third lattice cannot be half of
an electron. This arithmetic does not establish the physical charge unit. -/
theorem no_half_electron_in_native_charge_lattice (n : ℤ) :
    2 * n ≠ 3 ∧ 2 * n ≠ -3 := by omega

section Unitary
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℂ V]
  [CompleteSpace V] [FiniteDimensional ℂ V]

/-- The finite-matrix operator on its complete Hilbert space has a unitary
real-parameter flow. This theorem supplies no light cone or physical units. -/
theorem finite_source_unitary
    (H : V →L[ℂ] V) (hH : IsSelfAdjoint H) (t : ℝ) :
    NormedSpace.exp (t • ((-Complex.I) • H)) ∈ unitary (V →L[ℂ] V) :=
  OPH.Dynamics.hamiltonianPropagator_mem_unitary H hH t
end Unitary

#print axioms finite_source_gibbs
#print axioms stationary_density_fixed
#print axioms no_half_electron_in_native_charge_lattice
#print axioms finite_source_unitary

end CCF2.Kinematics
