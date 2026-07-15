import EventAlgebra.Lueders

/-!
# Trace-preserving pinching onto the block-diagonal commutant

A finite family of pairwise-orthogonal events summing to the identity (a
`ProjectivePartition`) generates the commutative algebra `D = span{Pᵢ}`.
This module constructs the trace-preserving pinching map

  `pinchingExpectation part X = ∑ i, Pᵢ X Pᵢ`,

whose range is the commutant `N = {X : X Pᵢ = Pᵢ X for every i}`. That
commutant is the block-diagonal algebra `⊕ᵢ M_{rank Pᵢ}(ℂ)`, generally
noncommutative; it is not the commutative algebra `D`. The generating
algebra `D` is precisely the center `Z(N)` of the commutant `N`. The
pinching is the conditional expectation onto `N`, and it satisfies, in the
finite quantum setting, the same package that characterises a classical
conditional expectation:

* **existence**: the pinching lands in the commutant of the partition,
  fixes the commutant pointwise, and is idempotent (a projector);
* **state preservation**: it maps states to states;
* **selfadjointness and contractivity**: it is selfadjoint for the trace
  inner product `⟨X, Y⟩ = Tr(Xᴴ Y)`, satisfies the Pythagorean identity,
  and is a squared-`L²` contraction;
* **uniqueness**: it is the only commutant-valued map compatible with the
  trace against every element of the commutant;
* **compatibility with conditioning**: for an event `P` in the commutant,
  Lüders conditioning commutes with the pinching. When the partition is
  rank one (every `Pᵢ` has rank one) the commutant `N` equals `D` and the
  pinching is the classical conditional expectation, so this recovers
  classical conditioning as the rank-one special case.

## Tagging convention

As in `EventAlgebra.Basic`: each lemma is tagged **algebra-only** or
**trace-dependent** according to whether its content passes through the
trace pairing.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {n k : ℕ}

/-- A **projective partition**: a finite family of pairwise-orthogonal
events summing to the identity. The family generates the commutative
algebra `D = span{Pᵢ}` (see `ProjectivePartition.proj_commute`). -/
structure ProjectivePartition (n k : ℕ) where
  /-- The orthogonal projections constituting the partition. -/
  proj : Fin k → Matrix (Fin n) (Fin n) ℂ
  /-- Every member of the partition is an event. -/
  isEvent : ∀ i, IsEvent (proj i)
  /-- Distinct members are orthogonal. -/
  orthogonal : ∀ i j, i ≠ j → proj i * proj j = 0
  /-- The partition is complete: the members sum to the sure event. -/
  complete : ∑ i, proj i = 1

namespace ProjectivePartition

variable (part : ProjectivePartition n k)

/-- **Algebra-only.** Product formula for partition members:
`Pᵢ Pⱼ = δᵢⱼ Pᵢ`. -/
theorem proj_mul_proj (i j : Fin k) :
    part.proj i * part.proj j = if i = j then part.proj i else 0 := by
  by_cases h : i = j
  · subst h
    simp [(part.isEvent i).2]
  · simp [h, part.orthogonal i j h]

/-- **Algebra-only.** Partition members commute pairwise: the partition
generates the commutative subalgebra `D = span{Pᵢ}`. -/
theorem proj_commute (i j : Fin k) :
    part.proj i * part.proj j = part.proj j * part.proj i := by
  by_cases h : i = j
  · rw [h]
  · rw [part.orthogonal i j h, part.orthogonal j i (Ne.symm h)]

end ProjectivePartition

/-- The **pinching onto the commutant** determined by a partition: the map
`X ↦ ∑ i, Pᵢ X Pᵢ`. -/
noncomputable def pinchingExpectation (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  ∑ i, part.proj i * X * part.proj i

variable (part : ProjectivePartition n k)

/-- **Algebra-only.** Absorption on the right: multiplying the pinching by
a partition member selects the corresponding block. -/
theorem pinchingExpectation_mul_proj (X : Matrix (Fin n) (Fin n) ℂ)
    (i : Fin k) :
    pinchingExpectation part X * part.proj i = part.proj i * X * part.proj i := by
  rw [pinchingExpectation, Finset.sum_mul]
  rw [Finset.sum_eq_single i
    (fun j _ hji => by
      rw [mul_assoc, part.orthogonal j i hji, mul_zero])
    (fun h => absurd (Finset.mem_univ i) h)]
  rw [mul_assoc, (part.isEvent i).2]

/-- **Algebra-only.** Absorption on the left. -/
theorem proj_mul_pinchingExpectation (X : Matrix (Fin n) (Fin n) ℂ)
    (i : Fin k) :
    part.proj i * pinchingExpectation part X = part.proj i * X * part.proj i := by
  rw [pinchingExpectation, Finset.mul_sum]
  rw [Finset.sum_eq_single i
    (fun j _ hji => by
      rw [← mul_assoc, ← mul_assoc, part.orthogonal i j (Ne.symm hji),
        zero_mul, zero_mul])
    (fun h => absurd (Finset.mem_univ i) h)]
  rw [← mul_assoc, ← mul_assoc, (part.isEvent i).2]

/-- **Algebra-only.** The pinching lands in the commutant of the partition:
its output commutes with every partition member. -/
theorem proj_commute_pinchingExpectation (X : Matrix (Fin n) (Fin n) ℂ)
    (i : Fin k) :
    part.proj i * pinchingExpectation part X = pinchingExpectation part X * part.proj i := by
  rw [proj_mul_pinchingExpectation, pinchingExpectation_mul_proj]

/-- **Algebra-only.** The pinching fixes the commutant pointwise: any `X`
commuting with every partition member is left unchanged. -/
theorem pinchingExpectation_fixes {X : Matrix (Fin n) (Fin n) ℂ}
    (h : ∀ i, X * part.proj i = part.proj i * X) :
    pinchingExpectation part X = X := by
  calc pinchingExpectation part X = ∑ i, X * part.proj i :=
        Finset.sum_congr rfl fun i _ => by
          rw [← h i, mul_assoc, (part.isEvent i).2]
    _ = X * ∑ i, part.proj i := (Finset.mul_sum _ _ _).symm
    _ = X := by rw [part.complete, mul_one]

/-- **Algebra-only.** The pinching is idempotent: it is a projector onto
the commutant. -/
theorem pinchingExpectation_idem (X : Matrix (Fin n) (Fin n) ℂ) :
    pinchingExpectation part (pinchingExpectation part X) = pinchingExpectation part X :=
  pinchingExpectation_fixes part fun i =>
    (proj_commute_pinchingExpectation part X i).symm

/-- **Algebra-only.** The pinching commutes with conjugate transposition. -/
theorem conjTranspose_pinchingExpectation (X : Matrix (Fin n) (Fin n) ℂ) :
    (pinchingExpectation part X)ᴴ = pinchingExpectation part Xᴴ := by
  rw [pinchingExpectation, pinchingExpectation, conjTranspose_sum]
  exact Finset.sum_congr rfl fun i _ => by
    simp only [conjTranspose_mul, (part.isEvent i).1.eq, mul_assoc]

/-- **Trace-dependent.** The pinching maps states to states: the pinching
of a density matrix is again a density matrix. -/
theorem pinchingExpectation_isState {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) : IsState (pinchingExpectation part ρ) := by
  constructor
  · exact posSemidef_sum _ fun i _ => by
      have := hρ.1.mul_mul_conjTranspose_same (part.proj i)
      rwa [(part.isEvent i).1.eq] at this
  · rw [pinchingExpectation, trace_sum]
    calc ∑ i, (part.proj i * ρ * part.proj i).trace
        = ∑ i, bornWeight ρ (part.proj i) :=
          Finset.sum_congr rfl fun i _ => trace_sandwich (part.isEvent i).2 ρ
      _ = bornWeight ρ (∑ i, part.proj i) := (bornWeight_sum ρ _ _).symm
      _ = 1 := by rw [part.complete, bornWeight_one hρ]

/-- **Trace-dependent.** Selfadjointness of the pinching for the
trace inner product `⟨X, Y⟩ = Tr(Xᴴ Y)`. -/
theorem trace_conjTranspose_pinchingExpectation_mul
    (X Y : Matrix (Fin n) (Fin n) ℂ) :
    ((pinchingExpectation part X)ᴴ * Y).trace =
      (Xᴴ * pinchingExpectation part Y).trace := by
  have hterm : ∀ i : Fin k,
      (part.proj i * Xᴴ * part.proj i * Y).trace =
        (Xᴴ * (part.proj i * Y * part.proj i)).trace := by
    intro i
    rw [trace_mul_comm (part.proj i * Xᴴ * part.proj i) Y,
      trace_mul_comm Xᴴ (part.proj i * Y * part.proj i),
      show Y * (part.proj i * Xᴴ * part.proj i) =
        (Y * part.proj i * Xᴴ) * part.proj i by simp only [mul_assoc],
      trace_mul_comm (Y * part.proj i * Xᴴ) (part.proj i)]
    simp only [mul_assoc]
  rw [conjTranspose_pinchingExpectation, pinchingExpectation, pinchingExpectation,
    Finset.sum_mul, Finset.mul_sum, trace_sum, trace_sum]
  exact Finset.sum_congr rfl fun i _ => hterm i

/-- **Trace-dependent.** The Pythagorean identity for the
pinching: the squared trace norm splits into the squared norm of the
conditional expectation plus the squared norm of the residual. This is the
quantum counterpart of the classical `L²` energy identity for conditional
expectations. -/
theorem trace_pinchingExpectation_pythagoras (X : Matrix (Fin n) (Fin n) ℂ) :
    (Xᴴ * X).trace =
      ((pinchingExpectation part X)ᴴ * pinchingExpectation part X).trace +
      ((X - pinchingExpectation part X)ᴴ * (X - pinchingExpectation part X)).trace := by
  have hTX : ((pinchingExpectation part X)ᴴ * X).trace =
      (Xᴴ * pinchingExpectation part X).trace :=
    trace_conjTranspose_pinchingExpectation_mul part X X
  have hTT : ((pinchingExpectation part X)ᴴ * pinchingExpectation part X).trace =
      (Xᴴ * pinchingExpectation part X).trace := by
    have := trace_conjTranspose_pinchingExpectation_mul part X
      (pinchingExpectation part X)
    rwa [pinchingExpectation_idem part X] at this
  simp only [conjTranspose_sub, sub_mul, mul_sub, trace_sub]
  rw [hTX, hTT]
  ring

/-- **Trace-dependent.** Squared-`L²` contractivity of the
pinching for the trace norm: `Tr((𝔼X)ᴴ (𝔼X)) ≤ Tr(Xᴴ X)` in the partial
order of `ℂ`. -/
theorem trace_pinchingExpectation_contraction (X : Matrix (Fin n) (Fin n) ℂ) :
    ((pinchingExpectation part X)ᴴ * pinchingExpectation part X).trace ≤
      (Xᴴ * X).trace := by
  rw [trace_pinchingExpectation_pythagoras part X]
  exact le_add_of_nonneg_right
    (posSemidef_conjTranspose_mul_self _).trace_nonneg

/-- **Trace-dependent.** Trace compatibility: against any element
`C` of the commutant of the partition, the pinching is invisible to the
trace pairing. This is the defining property of a conditional
expectation. -/
theorem trace_pinchingExpectation_mul_central (X C : Matrix (Fin n) (Fin n) ℂ)
    (hC : ∀ i, C * part.proj i = part.proj i * C) :
    (pinchingExpectation part X * C).trace = (X * C).trace := by
  have hterm : ∀ i : Fin k,
      ((part.proj i * X * part.proj i) * C).trace =
        (X * (C * part.proj i)).trace := by
    intro i
    rw [trace_mul_comm (part.proj i * X * part.proj i) C,
      show C * (part.proj i * X * part.proj i) =
        (C * part.proj i) * X * part.proj i by simp only [mul_assoc],
      trace_mul_comm ((C * part.proj i) * X) (part.proj i),
      show part.proj i * ((C * part.proj i) * X) =
        ((part.proj i * C) * part.proj i) * X by simp only [mul_assoc],
      ← hC i,
      show ((C * part.proj i) * part.proj i) * X =
        (C * (part.proj i * part.proj i)) * X by simp only [mul_assoc],
      (part.isEvent i).2, trace_mul_comm (C * part.proj i) X]
  have hsum : ∑ i, ((part.proj i * X * part.proj i) * C).trace =
      ∑ i, (X * (C * part.proj i)).trace :=
    Finset.sum_congr rfl fun i _ => hterm i
  rw [pinchingExpectation, Finset.sum_mul, trace_sum, hsum, ← trace_sum,
    ← Finset.mul_sum, ← Finset.mul_sum, part.complete, mul_one]

/-- **Trace-dependent.** Uniqueness of the conditional
expectation: any commutant-valued matrix `Y` that is trace-compatible with
`X` against the whole commutant equals the pinching of `X`. The proof
tests against `C = (Y - 𝔼X)ᴴ` and uses faithfulness of the trace. -/
theorem pinchingExpectation_unique {X Y : Matrix (Fin n) (Fin n) ℂ}
    (hYcomm : ∀ i, Y * part.proj i = part.proj i * Y)
    (hYtr : ∀ C : Matrix (Fin n) (Fin n) ℂ,
      (∀ i, C * part.proj i = part.proj i * C) → (Y * C).trace = (X * C).trace) :
    Y = pinchingExpectation part X := by
  set D := Y - pinchingExpectation part X with hD
  have hDcomm : ∀ i, D * part.proj i = part.proj i * D := by
    intro i
    rw [hD, sub_mul, mul_sub, hYcomm i,
      ← proj_commute_pinchingExpectation part X i]
  have hDHcomm : ∀ i, Dᴴ * part.proj i = part.proj i * Dᴴ := by
    intro i
    have := congrArg conjTranspose (hDcomm i)
    rw [conjTranspose_mul, conjTranspose_mul, (part.isEvent i).1.eq] at this
    exact this.symm
  have htr0 : (D * Dᴴ).trace = 0 := by
    have h1 : (Y * Dᴴ).trace = (X * Dᴴ).trace := hYtr Dᴴ hDHcomm
    have h2 : (pinchingExpectation part X * Dᴴ).trace = (X * Dᴴ).trace :=
      trace_pinchingExpectation_mul_central part X Dᴴ hDHcomm
    rw [hD, sub_mul, trace_sub, h1, h2, sub_self]
  have hD0 : D = 0 := trace_mul_conjTranspose_self_eq_zero_iff.mp htr0
  rw [← sub_eq_zero]
  exact hD0

/-- **Trace-dependent.** For an event `P` in the commutant of the
partition, the Born weight is invariant under the pinching. -/
theorem bornWeight_pinchingExpectation (ρ P : Matrix (Fin n) (Fin n) ℂ)
    (hPc : ∀ i, P * part.proj i = part.proj i * P) :
    bornWeight (pinchingExpectation part ρ) P = bornWeight ρ P :=
  trace_pinchingExpectation_mul_central part ρ P hPc

/-- **Trace-dependent.** **Lüders compatibility**: for an event `P` in the
commutant of the partition, Lüders conditioning commutes with the pinching.
When the partition is rank one (every `Pᵢ` has rank one) the commutant `N`
equals `D` and the pinching is the classical conditional expectation, so
this recovers classical conditioning as the rank-one special case. -/
theorem pinchingExpectation_luedersUpdate (ρ P : Matrix (Fin n) (Fin n) ℂ)
    (hPc : ∀ i, P * part.proj i = part.proj i * P) :
    pinchingExpectation part (luedersUpdate ρ P) =
      luedersUpdate (pinchingExpectation part ρ) P := by
  have hsand : ∀ i : Fin k,
      part.proj i * (P * ρ * P) * part.proj i =
        P * (part.proj i * ρ * part.proj i) * P := by
    intro i
    have h1 : part.proj i * (P * ρ * P) * part.proj i =
        (part.proj i * P) * ρ * (P * part.proj i) := by simp only [mul_assoc]
    have h2 : P * (part.proj i * ρ * part.proj i) * P =
        (P * part.proj i) * ρ * (part.proj i * P) := by simp only [mul_assoc]
    rw [h1, h2, hPc i]
  rw [luedersUpdate, luedersUpdate, bornWeight_pinchingExpectation part ρ P hPc,
    pinchingExpectation, pinchingExpectation]
  simp only [mul_smul_comm, smul_mul_assoc]
  rw [← Finset.smul_sum]
  congr 1
  rw [Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun i _ => hsand i

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms ProjectivePartition.proj_mul_proj
#print axioms ProjectivePartition.proj_commute
#print axioms pinchingExpectation_mul_proj
#print axioms proj_mul_pinchingExpectation
#print axioms proj_commute_pinchingExpectation
#print axioms pinchingExpectation_fixes
#print axioms pinchingExpectation_idem
#print axioms conjTranspose_pinchingExpectation
#print axioms pinchingExpectation_isState
#print axioms trace_conjTranspose_pinchingExpectation_mul
#print axioms trace_pinchingExpectation_pythagoras
#print axioms trace_pinchingExpectation_contraction
#print axioms trace_pinchingExpectation_mul_central
#print axioms pinchingExpectation_unique
#print axioms bornWeight_pinchingExpectation
#print axioms pinchingExpectation_luedersUpdate

end EventAlgebra
