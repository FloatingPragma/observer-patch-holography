import EventAlgebra.PartitionPinchingCP

/-!
# Conditional-expectation generators and their projector GKSL form

For a finite projective partition, this module defines the relaxation
generator

`L(X) = gamma (E(X) - X)`,

where `E` is partition pinching.  It proves that the fixed algebra of `L` is
exactly the partition commutant when `gamma` is nonzero.  For nonnegative real
`gamma`, the same generator is written in the declared finite GKSL form with
the partition projectors as jump operators and common rate `gamma`.

The file formalizes the displayed matrix identity, not a general predicate of
complete positivity or a bundled CPTP semigroup.  In particular, neither the
letters "GKSL" nor explicit jump syntax are used as a substitute for a formal
CP/CPTP interface.
-/

namespace OPH.Dynamics

open Matrix

variable {n k : ℕ}

abbrev CMat (n : ℕ) := Matrix (Fin n) (Fin n) ℂ

/-- Matrix conditional-expectation generator `gamma (E - I)`. -/
noncomputable def conditionalExpectationGenerator
    (gamma : ℝ) (part : EventAlgebra.ProjectivePartition n k) :
    CMat n →ₗ[ℂ] CMat n :=
  (gamma : ℂ) •
    (EventAlgebra.partitionPinchingLinearMap part - LinearMap.id)

@[simp]
theorem conditionalExpectationGenerator_apply
    (gamma : ℝ) (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    conditionalExpectationGenerator gamma part X =
      (gamma : ℂ) • (EventAlgebra.partitionPinching part X - X) := by
  simp [conditionalExpectationGenerator]

/-- One standard GKSL dissipator with jump operator `V`. -/
noncomputable def gkslDissipator (V X : CMat n) : CMat n :=
  V * X * Vᴴ - (2 : ℂ)⁻¹ • (Vᴴ * V * X + X * (Vᴴ * V))

/-- The declared projector-jump GKSL expression with a common real rate. -/
noncomputable def projectorGKSL
    (gamma : ℝ) (part : EventAlgebra.ProjectivePartition n k)
    (X : CMat n) : CMat n :=
  ∑ i, (gamma : ℂ) • gkslDissipator (part.proj i) X

/-- Summing the partition-projector dissipators gives exactly `E(X)-X`.
This is the load-bearing finite matrix identity behind the declared GKSL
form. -/
theorem sum_projector_gkslDissipator
    (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    (∑ i, gkslDissipator (part.proj i) X) =
      EventAlgebra.partitionPinching part X - X := by
  have hterm : ∀ i : Fin k,
      gkslDissipator (part.proj i) X =
        part.proj i * X * part.proj i -
          (2 : ℂ)⁻¹ • (part.proj i * X + X * part.proj i) := by
    intro i
    simp [gkslDissipator, (part.isEvent i).1.eq, (part.isEvent i).2,
      mul_assoc]
  rw [Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_sub_distrib,
    ← Finset.smul_sum, Finset.sum_add_distrib]
  have hleft : (∑ i, part.proj i * X) = X := by
    rw [← Finset.sum_mul, part.complete, one_mul]
  have hright : (∑ i, X * part.proj i) = X := by
    rw [← Finset.mul_sum, part.complete, mul_one]
  rw [hleft, hright]
  change (∑ i, part.proj i * X * part.proj i) - (2 : ℂ)⁻¹ • (X + X) =
    EventAlgebra.partitionPinching part X - X
  rw [show (2 : ℂ)⁻¹ = (1 / 2 : ℂ) by norm_num]
  simp [EventAlgebra.partitionPinching]
  module

/-- Exact declared GKSL identity.  The nonnegativity premise is what makes
`gamma` an admissible common GKSL rate; the matrix equality itself is
algebraic. -/
theorem conditionalExpectationGenerator_eq_projectorGKSL
    (gamma : ℝ) (_hgamma : 0 ≤ gamma)
    (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    conditionalExpectationGenerator gamma part X =
      projectorGKSL gamma part X := by
  rw [conditionalExpectationGenerator_apply, projectorGKSL,
    ← Finset.smul_sum, sum_projector_gkslDissipator]

/-- With nonzero rate, the pointer algebra (the generator kernel) is exactly
the fixed/commutant algebra of partition pinching. -/
theorem conditionalExpectationGenerator_eq_zero_iff_mem_commutant
    {gamma : ℝ} (hgamma : gamma ≠ 0)
    (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    conditionalExpectationGenerator gamma part X = 0 ↔
      X ∈ part.commutant := by
  rw [conditionalExpectationGenerator_apply]
  have hgammaC : (gamma : ℂ) ≠ 0 := by exact_mod_cast hgamma
  constructor
  · intro h
    have hdiff : EventAlgebra.partitionPinching part X - X = 0 :=
      (smul_eq_zero.mp h).resolve_left hgammaC
    exact (EventAlgebra.partitionPinching_eq_self_iff_mem_commutant part X).mp
      (sub_eq_zero.mp hdiff)
  · intro h
    have hfix :=
      (EventAlgebra.partitionPinching_eq_self_iff_mem_commutant part X).mpr h
    simp [hfix]

section MultiCollar

variable {C : Type*} [Fintype C]

/-- Sum of several collar generators acting on the same finite matrix
algebra. -/
noncomputable def multiCollarGenerator
    (gamma : C → ℝ)
    (part : C → EventAlgebra.ProjectivePartition n k) :
    CMat n →ₗ[ℂ] CMat n :=
  ∑ c, conditionalExpectationGenerator (gamma c) (part c)

@[simp]
theorem multiCollarGenerator_apply
    (gamma : C → ℝ)
    (part : C → EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    multiCollarGenerator gamma part X =
      ∑ c, conditionalExpectationGenerator (gamma c) (part c) X := by
  simp [multiCollarGenerator]

/-- Explicit no-cancellation premise for a sum of collar generators.  This
is not automatic for arbitrary linear maps; a future Hilbert-space
dissipativity package may discharge it from positivity and self-adjointness. -/
def MultiCollarNoCancellation
    (gamma : C → ℝ)
    (part : C → EventAlgebra.ProjectivePartition n k) : Prop :=
  ∀ X : CMat n,
    multiCollarGenerator gamma part X = 0 →
      ∀ c, conditionalExpectationGenerator (gamma c) (part c) X = 0

/-- Under explicit nonzero-rate and no-cancellation premises, the stable
algebra of the summed generator is exactly the intersection of the collar
commutants. -/
theorem multiCollarGenerator_eq_zero_iff_stableIntersection
    (gamma : C → ℝ) (hgamma : ∀ c, gamma c ≠ 0)
    (part : C → EventAlgebra.ProjectivePartition n k)
    (hnc : MultiCollarNoCancellation gamma part) (X : CMat n) :
    multiCollarGenerator gamma part X = 0 ↔
      ∀ c, X ∈ (part c).commutant := by
  constructor
  · intro h c
    exact (conditionalExpectationGenerator_eq_zero_iff_mem_commutant
      (hgamma c) (part c) X).mp (hnc X h c)
  · intro h
    rw [multiCollarGenerator_apply]
    apply Finset.sum_eq_zero
    intro c _
    exact (conditionalExpectationGenerator_eq_zero_iff_mem_commutant
      (hgamma c) (part c) X).mpr (h c)

end MultiCollar

end OPH.Dynamics

#print axioms OPH.Dynamics.conditionalExpectationGenerator_eq_projectorGKSL
#print axioms OPH.Dynamics.conditionalExpectationGenerator_eq_zero_iff_mem_commutant
#print axioms OPH.Dynamics.multiCollarGenerator_eq_zero_iff_stableIntersection
