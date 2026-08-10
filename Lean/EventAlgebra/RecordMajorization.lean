import EventAlgebra.PartitionPinching
import Mathlib.Logic.Equiv.Bool

/-!
# B9: the sign-average representation of record pinching

For a projective partition `Pᵢ`, let

`Uₛ = ∑ i, ε(s i) Pᵢ`, where each `ε(s i)` is independently `+1` or `-1`.

Every `Uₛ` is a self-adjoint unitary, and averaging `Uₛ X Uₛᴴ` over all
`2^k` sign assignments is exactly the partition pinching
`∑ i, Pᵢ X Pᵢ`.  This is the algebraic first step requested by B9 before
any spectral-majorization or entropy argument.

The result deliberately does not define matrix entropy or claim
majorization.  Those require the support-aware spectral layer named
separately in B9.  The final two-dimensional control shows that replacing
the independent sign family by the tempting single global sign leaves an
off-diagonal matrix unchanged and therefore is not a pinching argument.
-/

namespace EventAlgebra

open Matrix

variable {n k : ℕ}

/-- An independent sign choice for every block of a `k`-block partition. -/
abbrev RecordSign (k : ℕ) := Fin k → Bool

/-- There are exactly `2^k` independent sign assignments. -/
theorem card_recordSign (k : ℕ) : Fintype.card (RecordSign k) = 2 ^ k := by
  simp

/-- The real uniform weight on the independent sign family. -/
noncomputable def uniformRecordSignWeight (k : ℕ) : ℝ :=
  (2 ^ k : ℝ)⁻¹

theorem uniformRecordSignWeight_pos (k : ℕ) :
    0 < uniformRecordSignWeight k := by
  rw [uniformRecordSignWeight]
  positivity

/-- The uniform sign weights sum to one. -/
theorem sum_uniformRecordSignWeight (k : ℕ) :
    ∑ _s : RecordSign k, uniformRecordSignWeight k = 1 := by
  simp [uniformRecordSignWeight]

/-- The complex sign represented by a Boolean value. -/
def recordSignScalar (b : Bool) : ℂ := if b then -1 else 1

@[simp] theorem recordSignScalar_false : recordSignScalar false = 1 := rfl

@[simp] theorem recordSignScalar_true : recordSignScalar true = -1 := rfl

@[simp] theorem recordSignScalar_not (b : Bool) :
    recordSignScalar (!b) = -recordSignScalar b := by
  cases b <;> simp [recordSignScalar]

@[simp] theorem recordSignScalar_sq (b : Bool) :
    recordSignScalar b * recordSignScalar b = 1 := by
  cases b <;> norm_num [recordSignScalar]

@[simp] theorem star_recordSignScalar (b : Bool) :
    star (recordSignScalar b) = recordSignScalar b := by
  cases b <;> simp [recordSignScalar]

/-- Flip one coordinate of an independent sign assignment. -/
def flipRecordSign (i : Fin k) : RecordSign k ≃ RecordSign k where
  toFun s := Function.update s i (!s i)
  invFun s := Function.update s i (!s i)
  left_inv s := by
    funext j
    by_cases h : j = i
    · subst j
      simp
    · simp [Function.update, h]
  right_inv s := by
    funext j
    by_cases h : j = i
    · subst j
      simp
    · simp [Function.update, h]

@[simp] theorem flipRecordSign_same (s : RecordSign k) (i : Fin k) :
    flipRecordSign i s i = !s i := by
  simp [flipRecordSign]

@[simp] theorem flipRecordSign_other (s : RecordSign k) {i j : Fin k}
    (hij : j ≠ i) : flipRecordSign i s j = s j := by
  simp [flipRecordSign, Function.update, hij]

/-- Orthogonality of the independent sign characters. -/
theorem sum_recordSignScalar_mul (i j : Fin k) :
    (∑ s : RecordSign k, recordSignScalar (s i) * recordSignScalar (s j))
      = if i = j then (2 ^ k : ℂ) else 0 := by
  classical
  by_cases hij : i = j
  · subst j
    simp [recordSignScalar_sq]
  · rw [if_neg hij]
    let f : RecordSign k → ℂ := fun s =>
      recordSignScalar (s i) * recordSignScalar (s j)
    have hflip : ∀ s, f (flipRecordSign i s) = -f s := by
      intro s
      simp only [f, flipRecordSign_same,
        flipRecordSign_other s (show j ≠ i from Ne.symm hij),
        recordSignScalar_not]
      ring
    have hinv : (∑ s, f (flipRecordSign i s)) = ∑ s, f s :=
      (flipRecordSign i).sum_comp f
    have hneg : (∑ s, f (flipRecordSign i s)) = -(∑ s, f s) := by
      simp_rw [hflip, Finset.sum_neg_distrib]
    have hz : (∑ s, f s) = 0 := by
      have heq : (∑ s, f s) = -(∑ s, f s) := hinv.symm.trans hneg
      have htwo : (2 : ℂ) * (∑ s, f s) = 0 := by
        calc
          (2 : ℂ) * (∑ s, f s) = (∑ s, f s) + (∑ s, f s) := by ring
          _ = -(∑ s, f s) + (∑ s, f s) :=
            congrArg (fun z => z + ∑ s, f s) heq
          _ = 0 := neg_add_cancel _
      exact (mul_eq_zero.mp htwo).resolve_left (by norm_num)
    simpa [f] using hz

/-- The signed block reflection associated with one independent sign
assignment. -/
noncomputable def signedRecordUnitary (part : ProjectivePartition n k)
    (s : RecordSign k) : Matrix (Fin n) (Fin n) ℂ :=
  ∑ i, recordSignScalar (s i) • part.proj i

/-- Every signed block reflection is self-adjoint. -/
theorem signedRecordUnitary_star (part : ProjectivePartition n k)
    (s : RecordSign k) :
    (signedRecordUnitary part s)ᴴ = signedRecordUnitary part s := by
  classical
  rw [signedRecordUnitary, Matrix.conjTranspose_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.conjTranspose_smul, (part.isEvent i).1.eq,
    star_recordSignScalar]

/-- Every signed block reflection squares to the identity. -/
theorem signedRecordUnitary_sq (part : ProjectivePartition n k)
    (s : RecordSign k) :
    signedRecordUnitary part s * signedRecordUnitary part s = 1 := by
  classical
  rw [signedRecordUnitary, Finset.sum_mul]
  simp_rw [Finset.mul_sum, smul_mul_smul_comm, part.proj_mul_proj]
  calc
    (∑ i, ∑ j,
        (recordSignScalar (s i) * recordSignScalar (s j)) •
          (if i = j then part.proj i else 0))
        = ∑ i, part.proj i := by
            refine Finset.sum_congr rfl fun i _ => ?_
            rw [Finset.sum_eq_single i]
            · simp [recordSignScalar_sq]
            · intro j _ hji
              simp [Ne.symm hji]
            · intro hi
              exact (hi (Finset.mem_univ i)).elim
    _ = 1 := part.complete

/-- The signed reflections are unitary, in explicit matrix equations. -/
theorem signedRecordUnitary_isUnitary (part : ProjectivePartition n k)
    (s : RecordSign k) :
    (signedRecordUnitary part s)ᴴ * signedRecordUnitary part s = 1 ∧
      signedRecordUnitary part s * (signedRecordUnitary part s)ᴴ = 1 := by
  rw [signedRecordUnitary_star]
  exact ⟨signedRecordUnitary_sq part s, signedRecordUnitary_sq part s⟩

/-- Uniform average over all independent signed conjugations. -/
noncomputable def recordSignAverage (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  ((2 ^ k : ℂ)⁻¹) •
    ∑ s : RecordSign k,
      signedRecordUnitary part s * X * (signedRecordUnitary part s)ᴴ

/-- The displayed average really is a convex combination: it uses the
strictly positive real weight whose sum was proved above. -/
theorem recordSignAverage_eq_uniformCombination
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    recordSignAverage part X =
      ∑ s : RecordSign k, (uniformRecordSignWeight k : ℂ) •
        (signedRecordUnitary part s * X * (signedRecordUnitary part s)ᴴ) := by
  classical
  rw [recordSignAverage, Finset.smul_sum]
  refine Finset.sum_congr rfl fun s _ => ?_
  congr 1
  simp [uniformRecordSignWeight]

/-- The unnormalised signed-conjugation sum.  Character orthogonality kills
all cross-block corners and leaves `2^k` copies of every diagonal block. -/
theorem sum_signedRecordUnitary_conjugations
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    (∑ s : RecordSign k,
        signedRecordUnitary part s * X * (signedRecordUnitary part s)ᴴ)
      = (2 ^ k : ℂ) • partitionPinching part X := by
  classical
  calc
    (∑ s : RecordSign k,
        signedRecordUnitary part s * X * (signedRecordUnitary part s)ᴴ)
        = ∑ s : RecordSign k, ∑ i, ∑ j,
            (recordSignScalar (s i) * recordSignScalar (s j)) •
              (part.proj i * X * part.proj j) := by
            refine Finset.sum_congr rfl fun s _ => ?_
            conv_rhs => rw [Finset.sum_comm]
            rw [signedRecordUnitary_star, signedRecordUnitary,
              Finset.sum_mul]
            simp_rw [smul_mul_assoc]
            rw [Finset.mul_sum]
            refine Finset.sum_congr rfl fun j _ => ?_
            rw [Finset.sum_mul]
            simp_rw [smul_mul_smul_comm]
    _ = ∑ i, ∑ j, ∑ s : RecordSign k,
          (recordSignScalar (s i) * recordSignScalar (s j)) •
            (part.proj i * X * part.proj j) := by
          rw [Finset.sum_comm]
          refine Finset.sum_congr rfl fun i _ => ?_
          rw [Finset.sum_comm]
    _ = ∑ i, (2 ^ k : ℂ) • (part.proj i * X * part.proj i) := by
          refine Finset.sum_congr rfl fun i _ => ?_
          rw [Finset.sum_eq_single i]
          · rw [← Finset.sum_smul, sum_recordSignScalar_mul, if_pos rfl]
          · intro j _ hji
            rw [← Finset.sum_smul, sum_recordSignScalar_mul,
              if_neg (Ne.symm hji), zero_smul]
          · intro hi
            exact (hi (Finset.mem_univ i)).elim
    _ = (2 ^ k : ℂ) • partitionPinching part X := by
          rw [partitionPinching, Finset.smul_sum]

/-- **B9 sign-average theorem.** Arbitrary projective-partition pinching is
the uniform convex average of the `2^k` independently signed unitary
conjugations. -/
theorem recordSignAverage_eq_partitionPinching
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    recordSignAverage part X = partitionPinching part X := by
  rw [recordSignAverage, sum_signedRecordUnitary_conjugations, smul_smul]
  have hpow : (2 ^ k : ℂ) ≠ 0 := by norm_num
  rw [inv_mul_cancel₀ hpow, one_smul]

/-! ## Nondegeneracy control: one global sign is not enough -/

/-- The standard two-coordinate projective partition. -/
noncomputable def binaryCoordinatePartition : ProjectivePartition 2 2 where
  proj i := Matrix.diagonal (Pi.single i 1)
  isEvent i := by
    refine ⟨?_, ?_⟩
    · show (Matrix.diagonal (Pi.single i (1 : ℂ)))ᴴ =
        Matrix.diagonal (Pi.single i 1)
      rw [Matrix.diagonal_conjTranspose]
      congr 1
      funext j
      rcases eq_or_ne j i with h | h <;>
        simp [h, Pi.single_apply]
    · rw [Matrix.diagonal_mul_diagonal]
      congr 1
      funext j
      rcases eq_or_ne j i with h | h <;>
        simp [h, Pi.single_apply]
  orthogonal i j hij := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_zero]
    congr 1
    funext a
    rcases eq_or_ne a i with h | h
    · subst h
      simp [hij]
    · simp [h]
  complete := by
    ext a b
    rw [Matrix.sum_apply]
    rcases eq_or_ne a b with h | h
    · subst h
      simp [Matrix.diagonal_apply_eq, Pi.single_apply,
        Matrix.one_apply_eq]
    · simp [h]

/-- A nonzero off-diagonal matrix erased by coordinate pinching. -/
def binaryOffDiagonal : Matrix (Fin 2) (Fin 2) ℂ := fun i j =>
  if i = 0 ∧ j = 1 then 1 else 0

theorem binaryOffDiagonal_ne_zero : binaryOffDiagonal ≠ 0 := by
  intro h
  have := congr_fun (congr_fun h (0 : Fin 2)) (1 : Fin 2)
  norm_num [binaryOffDiagonal] at this

theorem binaryOffDiagonal_pinching_eq_zero :
    partitionPinching binaryCoordinatePartition binaryOffDiagonal = 0 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    simp [partitionPinching, binaryCoordinatePartition, binaryOffDiagonal,
      Matrix.mul_apply]

/-- Conjugating by only the two global signs `±I` and averaging leaves every
matrix unchanged.  Combined with the preceding off-diagonal witness, this
rules out the degenerate global-sign substitute for independent block signs. -/
theorem globalSignAverage_eq_self (X : Matrix (Fin n) (Fin n) ℂ) :
    ((2 : ℂ)⁻¹) • ((1 : Matrix (Fin n) (Fin n) ℂ) * X * 1 +
      (-1 : Matrix (Fin n) (Fin n) ℂ) * X * (-1)) = X := by
  ext i j
  simp
  ring

theorem globalSignAverage_not_binary_pinching :
    ((2 : ℂ)⁻¹) •
        ((1 : Matrix (Fin 2) (Fin 2) ℂ) * binaryOffDiagonal * 1 +
          (-1 : Matrix (Fin 2) (Fin 2) ℂ) * binaryOffDiagonal * (-1))
      ≠ partitionPinching binaryCoordinatePartition binaryOffDiagonal := by
  rw [globalSignAverage_eq_self, binaryOffDiagonal_pinching_eq_zero]
  exact binaryOffDiagonal_ne_zero

#print axioms EventAlgebra.sum_recordSignScalar_mul
#print axioms EventAlgebra.signedRecordUnitary_isUnitary
#print axioms EventAlgebra.sum_uniformRecordSignWeight
#print axioms EventAlgebra.recordSignAverage_eq_uniformCombination
#print axioms EventAlgebra.recordSignAverage_eq_partitionPinching
#print axioms EventAlgebra.globalSignAverage_not_binary_pinching

end EventAlgebra
