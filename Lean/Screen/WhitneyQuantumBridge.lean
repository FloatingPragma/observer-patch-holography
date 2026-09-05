import Mathlib.Algebra.MvPolynomial.PDeriv
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

set_option autoImplicit false

open scoped BigOperators

namespace OPH.WhitneyQuantumBridge

noncomputable section

/-!
# Canonical quantum oscillator algebra for supplied positive Whitney modes

The domain is the polynomial finite-particle domain in finitely many modes,
not a finite-dimensional truncation. Multiplication and differentiation
realize the bosonic canonical commutation relations. An explicit positive
Planck parameter is an input to canonical quantization, not a derived OPH
constant. The factorial-weight Hilbert completion and self-adjoint closures
are analytic constructions described separately; the identities below are
exact identities on their common invariant polynomial domain.

The normal-mode pullback is conditional on a supplied orthonormal spectral
frame for the actual mass and stiffness pairings. Its frequencies belong to
the ordinary continuous-time semidiscrete action, `omega^2=lambda`; they
are not the principal interpolation frequencies of a discrete prism step.
Gauge/zero-frequency directions, Gauss reduction, physical time, source
selection, and continuum convergence are not selected by this module.
-/

variable {ι : Type*}

abbrev ParticlePolynomial (ι : Type*) := MvPolynomial ι ℂ

def creation (i : ι) (p : ParticlePolynomial ι) : ParticlePolynomial ι :=
  MvPolynomial.X i * p

def annihilation (i : ι) (p : ParticlePolynomial ι) : ParticlePolynomial ι :=
  MvPolynomial.pderiv i p

def numberOperator (i : ι) (p : ParticlePolynomial ι) : ParticlePolynomial ι :=
  creation i (annihilation i p)

/-- The exact noncommutative quantum relation on the invariant domain. -/
theorem annihilation_creation [DecidableEq ι] (i j : ι) (p : ParticlePolynomial ι) :
    annihilation i (creation j p) - creation j (annihilation i p) =
      if i = j then p else 0 := by
  unfold annihilation creation
  rw [MvPolynomial.pderiv_mul]
  by_cases hij : i = j
  · subst j
    simp
  · rw [MvPolynomial.pderiv_X_of_ne (Ne.symm hij)]
    simp [hij]

theorem creation_commute (i j : ι) (p : ParticlePolynomial ι) :
    creation i (creation j p) = creation j (creation i p) := by
  unfold creation
  ring

theorem annihilation_commute (i j : ι) (p : ParticlePolynomial ι) :
    annihilation i (annihilation j p) = annihilation j (annihilation i p) := by
  classical
  by_cases hij : i = j
  · subst j; rfl
  unfold annihilation
  induction p using MvPolynomial.induction_on with
  | C a => simp
  | add p q hp hq => simp [hp, hq]
  | mul_X p k hp =>
    by_cases hi : k = i <;> by_cases hj : k = j <;>
      simp_all [MvPolynomial.pderiv_X] <;> ring

theorem annihilation_vacuum (i : ι) :
    annihilation i (1 : ParticlePolynomial ι) = 0 :=
  MvPolynomial.pderiv_one

/-- Unnormalized occupation monomials have integer number eigenvalues. -/
theorem numberOperator_monomial (i : ι) (n : ι →₀ ℕ) (c : ℂ) :
    numberOperator i (MvPolynomial.monomial n c) =
      (n i : ℂ) • MvPolynomial.monomial n c := by
  have h := MvPolynomial.X_mul_pderiv_monomial (i := i) (m := n) (r := c)
  simpa only [numberOperator, creation, annihilation, Nat.cast_smul_eq_nsmul] using h

/-- Canonical Hamiltonian on the finite-particle domain. -/
def quantumHamiltonian [Fintype ι] (hbar : ℝ) (omega : ι → ℝ)
    (p : ParticlePolynomial ι) : ParticlePolynomial ι :=
  ∑ i, ((hbar * omega i : ℝ) : ℂ) •
    (numberOperator i p + (1 / 2 : ℂ) • p)

def occupationEnergy [Fintype ι] (hbar : ℝ) (omega : ι → ℝ) (n : ι →₀ ℕ) : ℝ :=
  ∑ i, hbar * omega i * ((n i : ℝ) + 1 / 2)

/-- The full occupation spectrum, including the zero-point contribution. -/
theorem quantumHamiltonian_monomial [Fintype ι] (hbar : ℝ) (omega : ι → ℝ)
    (n : ι →₀ ℕ) (c : ℂ) :
    quantumHamiltonian hbar omega (MvPolynomial.monomial n c) =
      (occupationEnergy hbar omega n : ℂ) • MvPolynomial.monomial n c := by
  unfold quantumHamiltonian occupationEnergy
  simp_rw [numberOperator_monomial, ← add_smul, smul_smul]
  rw [← Finset.sum_smul]
  congr 1
  simp

theorem occupationEnergy_nonnegative [Fintype ι] (hbar : ℝ) (hh : 0 < hbar)
    (omega : ι → ℝ) (ho : ∀ i, 0 < omega i) (n : ι →₀ ℕ) :
    0 ≤ occupationEnergy hbar omega n := by
  unfold occupationEnergy
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (mul_pos hh (ho i)).le (by positivity)

/-- Position and momentum use the same real normalization scale. -/
def positionOperator (i : ι) (scale : ℝ) (p : ParticlePolynomial ι) :
    ParticlePolynomial ι :=
  (scale : ℂ) • (creation i p + annihilation i p)

def momentumOperator (i : ι) (omega scale : ℝ) (p : ParticlePolynomial ι) :
    ParticlePolynomial ι :=
  (Complex.I * (omega : ℂ) * (scale : ℂ)) •
    (creation i p - annihilation i p)

/-- Position/momentum CCR before fixing the normalization. -/
theorem position_momentum_commutator (i : ι) (omega scale : ℝ)
    (p : ParticlePolynomial ι) :
    positionOperator i scale (momentumOperator i omega scale p) -
      momentumOperator i omega scale (positionOperator i scale p) =
        (2 * Complex.I * (omega : ℂ) * (scale : ℂ) ^ 2) • p := by
  simp only [positionOperator, momentumOperator, creation, annihilation,
    MvPolynomial.pderiv_mul, MvPolynomial.pderiv_X_self,
    Derivation.map_smul, map_add, map_sub, mul_smul_comm, mul_add, mul_sub,
    one_mul, smul_add, smul_sub, smul_smul]
  module

/-- Exact oscillator quadratic Hamiltonian, before fixing Planck scale. -/
theorem oscillator_quadratic_identity (i : ι) (omega scale : ℝ)
    (p : ParticlePolynomial ι) :
    (1 / 2 : ℂ) •
      (momentumOperator i omega scale (momentumOperator i omega scale p) +
        ((omega ^ 2 : ℝ) : ℂ) • positionOperator i scale (positionOperator i scale p)) =
      ((omega ^ 2 * scale ^ 2 : ℝ) : ℂ) •
        ((2 : ℂ) • numberOperator i p + p) := by
  have hi : (Complex.I * (omega : ℂ) * (scale : ℂ)) *
      (Complex.I * (omega : ℂ) * (scale : ℂ)) =
      -((omega : ℂ) ^ 2 * (scale : ℂ) ^ 2) := by
    calc
      _ = Complex.I ^ 2 * (omega : ℂ) ^ 2 * (scale : ℂ) ^ 2 := by ring
      _ = _ := by rw [Complex.I_sq]; ring
  simp only [positionOperator, momentumOperator, numberOperator, creation, annihilation,
    MvPolynomial.pderiv_mul, MvPolynomial.pderiv_X_self,
    Derivation.map_smul, map_add, map_sub, mul_smul_comm, mul_add, mul_sub,
    one_mul, smul_add, smul_sub, smul_smul, Complex.ofReal_pow, Complex.ofReal_mul, hi]
  module

def positionScale (hbar omega : ℝ) : ℝ := Real.sqrt (hbar / (2 * omega))

theorem positionScale_sq (hbar omega : ℝ) (hh : 0 < hbar) (ho : 0 < omega) :
    positionScale hbar omega ^ 2 = hbar / (2 * omega) := by
  exact Real.sq_sqrt (by positivity)

/-- With canonical normalization the commutator is `i hbar` times identity. -/
theorem canonical_position_momentum (i : ι) (hbar omega : ℝ)
    (hh : 0 < hbar) (ho : 0 < omega) (p : ParticlePolynomial ι) :
    positionOperator i (positionScale hbar omega)
        (momentumOperator i omega (positionScale hbar omega) p) -
      momentumOperator i omega (positionScale hbar omega)
        (positionOperator i (positionScale hbar omega) p) =
      (Complex.I * (hbar : ℂ)) • p := by
  rw [position_momentum_commutator]
  congr 1
  have hs : ((positionScale hbar omega : ℝ) : ℂ) ^ 2 =
      ((hbar / (2 * omega) : ℝ) : ℂ) := by
    exact_mod_cast positionScale_sq hbar omega hh ho
  rw [hs]
  push_cast
  have hone : (omega : ℂ) ≠ 0 := by exact_mod_cast ho.ne'
  field_simp [hone]

/-- The canonical quadratic energy is exactly the occupation Hamiltonian. -/
theorem canonical_oscillator_energy (i : ι) (hbar omega : ℝ)
    (hh : 0 < hbar) (ho : 0 < omega) (p : ParticlePolynomial ι) :
    (1 / 2 : ℂ) •
      (momentumOperator i omega (positionScale hbar omega)
          (momentumOperator i omega (positionScale hbar omega) p) +
        ((omega ^ 2 : ℝ) : ℂ) •
          positionOperator i (positionScale hbar omega)
            (positionOperator i (positionScale hbar omega) p)) =
      ((hbar * omega : ℝ) : ℂ) •
        (numberOperator i p + (1 / 2 : ℂ) • p) := by
  rw [oscillator_quadratic_identity, positionScale_sq hbar omega hh ho]
  have hc : omega ^ 2 * (hbar / (2 * omega)) = hbar * omega / 2 := by
    field_simp
  rw [hc]
  push_cast
  module

/-- One mode's Hamiltonian, acting on the full polynomial domain. -/
def modeHamiltonian (i : ι) (hbar omega : ℝ) (p : ParticlePolynomial ι) :
    ParticlePolynomial ι :=
  ((hbar * omega : ℝ) : ℂ) • (numberOperator i p + (1 / 2 : ℂ) • p)

/-- Exact Heisenberg equation for the position operator. -/
theorem modeHamiltonian_position (i : ι) (hbar omega scale : ℝ)
    (p : ParticlePolynomial ι) :
    modeHamiltonian i hbar omega (positionOperator i scale p) -
      positionOperator i scale (modeHamiltonian i hbar omega p) =
        (-Complex.I * (hbar : ℂ)) • momentumOperator i omega scale p := by
  have hi : (-Complex.I * (hbar : ℂ)) *
      (Complex.I * (omega : ℂ) * (scale : ℂ)) =
      (hbar : ℂ) * (omega : ℂ) * (scale : ℂ) := by
    calc
      _ = -Complex.I ^ 2 * (hbar : ℂ) * (omega : ℂ) * (scale : ℂ) := by ring
      _ = _ := by rw [Complex.I_sq]; ring
  simp only [modeHamiltonian, positionOperator, momentumOperator, numberOperator,
    creation, annihilation, MvPolynomial.pderiv_mul, MvPolynomial.pderiv_X_self,
    Derivation.map_smul, map_add, mul_smul_comm, mul_add,
    one_mul, smul_add, smul_sub, smul_smul, Complex.ofReal_mul, hi]
  module

/-- Exact Heisenberg equation for the momentum operator. -/
theorem modeHamiltonian_momentum (i : ι) (hbar omega scale : ℝ)
    (p : ParticlePolynomial ι) :
    modeHamiltonian i hbar omega (momentumOperator i omega scale p) -
      momentumOperator i omega scale (modeHamiltonian i hbar omega p) =
        (Complex.I * (hbar : ℂ) * (omega : ℂ) ^ 2) • positionOperator i scale p := by
  simp only [modeHamiltonian, positionOperator, momentumOperator, numberOperator,
    creation, annihilation, MvPolynomial.pderiv_mul, MvPolynomial.pderiv_X_self,
    Derivation.map_smul, map_add, map_sub, mul_smul_comm, mul_add, mul_sub,
    one_mul, smul_add, smul_sub, smul_smul, Complex.ofReal_mul]
  module

theorem modeHamiltonian_position_other (i j : ι) (hij : i ≠ j)
    (hbar omega scale : ℝ) (p : ParticlePolynomial ι) :
    modeHamiltonian i hbar omega (positionOperator j scale p) =
      positionOperator j scale (modeHamiltonian i hbar omega p) := by
  have hd : (MvPolynomial.pderiv i) ((MvPolynomial.pderiv j) p) =
      (MvPolynomial.pderiv j) ((MvPolynomial.pderiv i) p) :=
    annihilation_commute i j p
  simp only [modeHamiltonian, positionOperator, numberOperator,
    creation, annihilation, MvPolynomial.pderiv_mul,
    MvPolynomial.pderiv_X_of_ne hij, MvPolynomial.pderiv_X_of_ne (Ne.symm hij),
    Derivation.map_smul, map_add, mul_smul_comm, mul_add,
    zero_mul, zero_add, smul_add, smul_smul, hd]
  rw [mul_left_comm (MvPolynomial.X i) (MvPolynomial.X j)]
  module

theorem modeHamiltonian_momentum_other (i j : ι) (hij : i ≠ j)
    (hbar omega omega' scale : ℝ) (p : ParticlePolynomial ι) :
    modeHamiltonian i hbar omega (momentumOperator j omega' scale p) =
      momentumOperator j omega' scale (modeHamiltonian i hbar omega p) := by
  have hd : (MvPolynomial.pderiv i) ((MvPolynomial.pderiv j) p) =
      (MvPolynomial.pderiv j) ((MvPolynomial.pderiv i) p) :=
    annihilation_commute i j p
  simp only [modeHamiltonian, momentumOperator, numberOperator,
    creation, annihilation, MvPolynomial.pderiv_mul,
    MvPolynomial.pderiv_X_of_ne hij, MvPolynomial.pderiv_X_of_ne (Ne.symm hij),
    Derivation.map_smul, map_add, map_sub, mul_smul_comm, mul_add, mul_sub,
    zero_mul, zero_add, smul_add, smul_sub, smul_smul, hd]
  rw [mul_left_comm (MvPolynomial.X i) (MvPolynomial.X j)]
  module

/-- All other modes commute, so the full Hamiltonian gives the same equation. -/
theorem quantumHamiltonian_position [Fintype ι] (i : ι)
    (hbar : ℝ) (omega : ι → ℝ) (scale : ℝ) (p : ParticlePolynomial ι) :
    quantumHamiltonian hbar omega (positionOperator i scale p) -
      positionOperator i scale (quantumHamiltonian hbar omega p) =
        (-Complex.I * (hbar : ℂ)) • momentumOperator i (omega i) scale p := by
  classical
  have hsum (f : ι → ParticlePolynomial ι) :
      positionOperator i scale (∑ j, f j) = ∑ j, positionOperator i scale (f j) := by
    simp only [positionOperator, creation, annihilation, map_sum, Finset.mul_sum,
      smul_add, Finset.sum_add_distrib, Finset.smul_sum]
  change (∑ j, modeHamiltonian j hbar (omega j) (positionOperator i scale p)) -
    positionOperator i scale (∑ j, modeHamiltonian j hbar (omega j) p) = _
  rw [hsum, ← Finset.sum_sub_distrib]
  rw [Finset.sum_eq_single i]
  · exact modeHamiltonian_position i hbar (omega i) scale p
  · intro j _ hji
    rw [modeHamiltonian_position_other j i hji, sub_self]
  · simp

theorem quantumHamiltonian_momentum [Fintype ι] (i : ι)
    (hbar : ℝ) (omega : ι → ℝ) (scale : ℝ) (p : ParticlePolynomial ι) :
    quantumHamiltonian hbar omega (momentumOperator i (omega i) scale p) -
      momentumOperator i (omega i) scale (quantumHamiltonian hbar omega p) =
        (Complex.I * (hbar : ℂ) * (omega i : ℂ) ^ 2) • positionOperator i scale p := by
  classical
  have hsum (f : ι → ParticlePolynomial ι) :
      momentumOperator i (omega i) scale (∑ j, f j) =
        ∑ j, momentumOperator i (omega i) scale (f j) := by
    simp only [momentumOperator, creation, annihilation, map_sum, Finset.mul_sum,
      smul_sub, Finset.sum_sub_distrib, Finset.smul_sum]
  change (∑ j, modeHamiltonian j hbar (omega j)
    (momentumOperator i (omega i) scale p)) -
    momentumOperator i (omega i) scale (∑ j, modeHamiltonian j hbar (omega j) p) = _
  rw [hsum, ← Finset.sum_sub_distrib]
  rw [Finset.sum_eq_single i]
  · exact modeHamiltonian_momentum i hbar (omega i) scale p
  · intro j _ hji
    rw [modeHamiltonian_momentum_other j i hji, sub_self]
  · simp

section NormalModes

variable [Fintype ι] [DecidableEq ι]
variable {E V : Type*} [AddCommGroup E] [Module ℝ E]
  [AddCommGroup V] [Module ℝ V]

/-- A supplied positive spectral frame must exhaust the constrained sector.
For the Whitney cone, `E` is the 42-edge space, the constraint is `Dᵀ M₁`,
the pairings are `M₁` and `Cᵀ M₂ C`, and there are 30 modes. Existence and
this concrete geometric identification are separate from the algebra here. -/
structure PositiveNormalFrame (mass stiffness : LinearMap.BilinForm ℝ E)
    (constraint : E →ₗ[ℝ] V) where
  vector : ι → E
  omega : ι → ℝ
  positive : ∀ i, 0 < omega i
  mass_orthonormal : ∀ i j, mass (vector i) (vector j) = if i = j then 1 else 0
  eigenmode : ∀ i x, stiffness x (vector i) = omega i ^ 2 * mass x (vector i)
  complete : Submodule.span ℝ (Set.range vector) = LinearMap.ker constraint

def reconstruct (v : ι → E) (q : ι → ℝ) : E := ∑ i, q i • v i

/-- The sum expansion establishing action equality from diagonal pairings. -/
theorem diagonal_pairing (pairing : LinearMap.BilinForm ℝ E) (v : ι → E)
    (weight : ι → ℝ)
    (hdiag : ∀ i j, pairing (v i) (v j) = if i = j then weight i else 0)
    (q r : ι → ℝ) :
    pairing (reconstruct v q) (reconstruct v r) = ∑ i, weight i * q i * r i := by
  simp only [reconstruct, map_sum, LinearMap.sum_apply, map_smul,
    LinearMap.smul_apply, smul_eq_mul]
  apply Finset.sum_congr rfl
  intro i _
  simp [hdiag]
  ring

variable {mass stiffness : LinearMap.BilinForm ℝ E} {constraint : E →ₗ[ℝ] V}

omit [Fintype ι] in
theorem frame_stiffness_diagonal (frame : PositiveNormalFrame (ι := ι)
    mass stiffness constraint) (i j : ι) :
    stiffness (frame.vector i) (frame.vector j) =
      if i = j then frame.omega i ^ 2 else 0 := by
  rw [frame.eigenmode, frame.mass_orthonormal]
  split_ifs with h
  · subst j; simp
  · simp

/-- The reconstructed fields remain in the full Gauss-reduced sector. -/
theorem reconstruct_constrained (frame : PositiveNormalFrame (ι := ι)
    mass stiffness constraint) (q : ι → ℝ) :
    constraint (reconstruct frame.vector q) = 0 := by
  apply LinearMap.mem_ker.mp
  rw [← frame.complete]
  apply Submodule.sum_mem
  intro i _
  exact Submodule.smul_mem _ _ (Submodule.subset_span (Set.mem_range_self i))

/-- Every constrained field is represented; a proper subset of modes is insufficient. -/
theorem reconstruct_surjective (frame : PositiveNormalFrame (ι := ι)
    mass stiffness constraint) (x : E) (hx : constraint x = 0) :
    ∃ q : ι → ℝ, reconstruct frame.vector q = x := by
  have hm : x ∈ Submodule.span ℝ (Set.range frame.vector) := by
    rw [frame.complete]
    exact LinearMap.mem_ker.mpr hx
  exact (Submodule.mem_span_range_iff_exists_fun ℝ).mp hm

/-- Exact pullback of the same semidiscrete Maxwell Lagrangian. -/
theorem same_action_normal_modes (frame : PositiveNormalFrame (ι := ι)
    mass stiffness constraint) (q velocity : ι → ℝ) :
    (mass (reconstruct frame.vector velocity) (reconstruct frame.vector velocity) -
      stiffness (reconstruct frame.vector q) (reconstruct frame.vector q)) / 2 =
      ∑ i, (velocity i ^ 2 - frame.omega i ^ 2 * q i ^ 2) / 2 := by
  rw [diagonal_pairing mass frame.vector (fun _ => 1) frame.mass_orthonormal,
    diagonal_pairing stiffness frame.vector (fun i => frame.omega i ^ 2)
      (frame_stiffness_diagonal frame)]
  simp only [one_mul, ← Finset.sum_sub_distrib, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- Quantizing those same frequencies gives the sum of quadratic energies. -/
theorem same_modes_quantized_energy (frame : PositiveNormalFrame (ι := ι)
    mass stiffness constraint) (hbar : ℝ) (hh : 0 < hbar)
    (p : ParticlePolynomial ι) :
    (∑ i, (1 / 2 : ℂ) •
      (momentumOperator i (frame.omega i) (positionScale hbar (frame.omega i))
          (momentumOperator i (frame.omega i) (positionScale hbar (frame.omega i)) p) +
        ((frame.omega i ^ 2 : ℝ) : ℂ) •
          positionOperator i (positionScale hbar (frame.omega i))
            (positionOperator i (positionScale hbar (frame.omega i)) p))) =
      quantumHamiltonian hbar frame.omega p := by
  unfold quantumHamiltonian
  apply Finset.sum_congr rfl
  intro i _
  exact canonical_oscillator_energy i hbar (frame.omega i) hh (frame.positive i) p

end NormalModes

end

#print axioms annihilation_creation
#print axioms quantumHamiltonian_monomial
#print axioms position_momentum_commutator
#print axioms oscillator_quadratic_identity
#print axioms canonical_position_momentum
#print axioms canonical_oscillator_energy
#print axioms annihilation_commute
#print axioms modeHamiltonian_position
#print axioms modeHamiltonian_momentum
#print axioms same_action_normal_modes
#print axioms quantumHamiltonian_position
#print axioms quantumHamiltonian_momentum
#print axioms reconstruct_surjective
#print axioms same_modes_quantized_energy

end OPH.WhitneyQuantumBridge
