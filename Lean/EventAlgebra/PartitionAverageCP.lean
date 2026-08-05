import EventAlgebra.PartitionAverage

/-!
# Explicit Kraus data for partition averaging

This module constructs a Kraus family for the trace-preserving projection
onto the commutative span of a finite projective partition.  The construction
uses the projected coordinate vectors as a finite Parseval frame.  For a
block projector `P_i` and coordinate matrix units `e_ab`, its Kraus operators
are

`K_iab = (Tr P_i)^(-1/2) P_i e_ab P_i`.

Zero blocks contribute zero.  The module proves, directly in the repository's
finite matrix interface, both load-bearing identities:

* `sum K_iab^H K_iab = 1` (normalization), and
* `partitionAverage X = sum K_iab X K_iab^H` (explicit Kraus form).

Together with the already-proved trace identity for `partitionAverage`, these
are exact Kraus/trace statements.  This file deliberately does not call the
map completely positive or CPTP: no finite complete-positivity predicate or
bundled channel interface is introduced here.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {n k : ℕ}

/-- Coordinate matrix unit `e_ab`. -/
def matrixUnit (a b : Fin n) : Matrix (Fin n) (Fin n) ℂ :=
  Matrix.single a b 1

@[simp]
theorem matrixUnit_apply (a b r s : Fin n) :
    matrixUnit a b r s = if a = r ∧ b = s then 1 else 0 := by
  simp [matrixUnit, Matrix.single]

/-- Summing the two-sided matrix-unit conjugations extracts the trace. -/
theorem sum_matrixUnit_mul_mul_swap (Y : Matrix (Fin n) (Fin n) ℂ) :
    (∑ a : Fin n, ∑ b : Fin n, matrixUnit a b * Y * matrixUnit b a) =
      Y.trace • (1 : Matrix (Fin n) (Fin n) ℂ) := by
  ext r s
  simp only [Matrix.sum_apply]
  have hterm : ∀ a b : Fin n,
      (matrixUnit a b * Y * matrixUnit b a) r s =
        if r = a ∧ s = a then Y b b else 0 := by
    intro a b
    by_cases hra : r = a
    · subst r
      by_cases hsa : s = a
      · subst s
        simp [matrixUnit, Matrix.single, Matrix.mul_apply]
      · simp [matrixUnit, Matrix.single, Matrix.mul_apply, hsa, Ne.symm hsa]
    · simp [matrixUnit, Matrix.single, Matrix.mul_apply, hra, Ne.symm hra]
  rw [Finset.sum_congr rfl fun a _ =>
    Finset.sum_congr rfl fun b _ => hterm a b]
  by_cases hrs : r = s
  · subst s
    simp [Matrix.trace]
  · have hzero : ∀ a : Fin n,
        (∑ b : Fin n, if r = a ∧ s = a then Y b b else 0) = 0 := by
      intro a
      by_cases hra : r = a
      · subst a
        simp [Ne.symm hrs]
      · simp [hra]
    rw [Finset.sum_eq_zero fun a _ => hzero a]
    simp [hrs]

/-- The same trace extractor with the two matrix-unit indices swapped. -/
theorem sum_matrixUnit_swap_mul_mul (Y : Matrix (Fin n) (Fin n) ℂ) :
    (∑ a : Fin n, ∑ b : Fin n, matrixUnit b a * Y * matrixUnit a b) =
      Y.trace • (1 : Matrix (Fin n) (Fin n) ℂ) := by
  rw [Finset.sum_comm]
  exact sum_matrixUnit_mul_mul_swap Y

/-- The real square-root normalization used by one partition block. -/
noncomputable def partitionBlockScale
    (part : ProjectivePartition n k) (i : Fin k) : ℂ :=
  ((Real.sqrt (part.proj i).trace.re)⁻¹ : ℝ)

/-- A positive block trace is represented by its real part. -/
theorem trace_proj_eq_ofReal_re
    (part : ProjectivePartition n k) (i : Fin k) :
    (part.proj i).trace = ((part.proj i).trace.re : ℂ) := by
  have hnonneg : (0 : ℂ) ≤ (part.proj i).trace :=
    (part.isEvent i).posSemidef.trace_nonneg
  have him : (part.proj i).trace.im = 0 :=
    (Complex.nonneg_iff.mp hnonneg).2.symm
  apply Complex.ext
  · simp
  · simp [him]

/-- Squaring the real block scale gives the inverse block trace, including
the zero-block case through totalized inversion. -/
theorem star_partitionBlockScale_mul_self
    (part : ProjectivePartition n k) (i : Fin k) :
    star (partitionBlockScale part i) * partitionBlockScale part i =
      ((part.proj i).trace)⁻¹ := by
  have hnonneg : 0 ≤ (part.proj i).trace.re :=
    (Complex.nonneg_iff.mp (part.isEvent i).posSemidef.trace_nonneg).1
  rw [trace_proj_eq_ofReal_re part i]
  by_cases hzero : (part.proj i).trace.re = 0
  · simp [partitionBlockScale, hzero]
  · have hsqrt : Real.sqrt (part.proj i).trace.re ≠ 0 := by
      exact (Real.sqrt_ne_zero').mpr (lt_of_le_of_ne hnonneg (Ne.symm hzero))
    have hstar :
        star ((((Real.sqrt (part.proj i).trace.re)⁻¹ : ℝ) : ℂ)) =
          (((Real.sqrt (part.proj i).trace.re)⁻¹ : ℝ) : ℂ) := by
      simp
    simp only [partitionBlockScale]
    rw [hstar]
    change
      (((Real.sqrt (part.proj i).trace.re)⁻¹ : ℝ) : ℂ) *
          (((Real.sqrt (part.proj i).trace.re)⁻¹ : ℝ) : ℂ) =
        (((part.proj i).trace.re : ℝ) : ℂ)⁻¹
    norm_cast
    field_simp
    exact (Real.sq_sqrt hnonneg).symm

/-- Parseval-frame Kraus operator for block `i` and coordinate pair `(a,b)`. -/
noncomputable def partitionAverageKraus
    (part : ProjectivePartition n k) (i : Fin k) (a b : Fin n) :
    Matrix (Fin n) (Fin n) ℂ :=
  partitionBlockScale part i •
    (part.proj i * matrixUnit a b * part.proj i)

/-- Conjugate transposition swaps the coordinate indices of a block Kraus
operator. -/
theorem conjTranspose_partitionAverageKraus
    (part : ProjectivePartition n k) (i : Fin k) (a b : Fin n) :
    (partitionAverageKraus part i a b)ᴴ =
      partitionAverageKraus part i b a := by
  simp [partitionAverageKraus, partitionBlockScale, matrixUnit,
    Matrix.conjTranspose_mul, (part.isEvent i).1.eq, mul_assoc]

/-- One block's Parseval-frame Kraus family is complete on that block. -/
theorem partitionAverageKraus_block_complete
    (part : ProjectivePartition n k) (i : Fin k) :
    (∑ a : Fin n, ∑ b : Fin n,
      (partitionAverageKraus part i a b)ᴴ *
        partitionAverageKraus part i a b) = part.proj i := by
  by_cases htrace : (part.proj i).trace = 0
  · have hP : part.proj i = 0 :=
      ((part.isEvent i).posSemidef.trace_eq_zero_iff).mp htrace
    simp [partitionAverageKraus, hP]
  · have hscale := star_partitionBlockScale_mul_self part i
    have hterm : ∀ a b : Fin n,
        (partitionAverageKraus part i a b)ᴴ *
            partitionAverageKraus part i a b =
          ((part.proj i).trace)⁻¹ •
            (part.proj i * matrixUnit b a * part.proj i *
              matrixUnit a b * part.proj i) := by
      intro a b
      rw [conjTranspose_partitionAverageKraus]
      simp only [partitionAverageKraus, smul_mul_assoc, mul_smul_comm,
        smul_smul]
      have hstar : star (partitionBlockScale part i) =
          partitionBlockScale part i := by
        simp [partitionBlockScale]
      have hmul : partitionBlockScale part i * partitionBlockScale part i =
          ((part.proj i).trace)⁻¹ := by
        calc
          partitionBlockScale part i * partitionBlockScale part i =
              star (partitionBlockScale part i) *
                partitionBlockScale part i := by rw [hstar]
          _ = ((part.proj i).trace)⁻¹ := hscale
      rw [hmul]
      congr 1
      calc
        part.proj i * matrixUnit b a * part.proj i *
            (part.proj i * matrixUnit a b * part.proj i) =
          part.proj i * matrixUnit b a *
            (part.proj i * part.proj i) * matrixUnit a b *
              part.proj i := by simp only [mul_assoc]
        _ = part.proj i * matrixUnit b a * part.proj i *
              matrixUnit a b * part.proj i := by
            rw [(part.isEvent i).2]
    rw [Finset.sum_congr rfl fun a _ =>
      Finset.sum_congr rfl fun b _ => hterm a b]
    have hsandwich :
        (∑ a : Fin n, ∑ b : Fin n,
          part.proj i * matrixUnit b a * part.proj i *
            matrixUnit a b * part.proj i) =
          (part.proj i).trace • part.proj i := by
      calc
        (∑ a : Fin n, ∑ b : Fin n,
          part.proj i * matrixUnit b a * part.proj i *
            matrixUnit a b * part.proj i) =
            part.proj i *
              (∑ a : Fin n, ∑ b : Fin n,
                matrixUnit b a * part.proj i * matrixUnit a b) *
              part.proj i := by
                simp only [Finset.mul_sum, Finset.sum_mul, mul_assoc]
        _ = part.proj i *
              ((part.proj i).trace •
                (1 : Matrix (Fin n) (Fin n) ℂ)) * part.proj i := by
              rw [sum_matrixUnit_swap_mul_mul]
        _ = (part.proj i).trace • part.proj i := by
              simp [(part.isEvent i).2]
    calc
      (∑ a : Fin n, ∑ b : Fin n,
        (part.proj i).trace⁻¹ •
          (part.proj i * matrixUnit b a * part.proj i *
            matrixUnit a b * part.proj i)) =
          (part.proj i).trace⁻¹ •
            (∑ a : Fin n, ∑ b : Fin n,
              part.proj i * matrixUnit b a * part.proj i *
                matrixUnit a b * part.proj i) := by
            rw [Finset.smul_sum]
            apply Finset.sum_congr rfl
            intro a _
            rw [Finset.smul_sum]
      _ = part.proj i := by
        rw [hsandwich, smul_smul, inv_mul_cancel₀ htrace, one_smul]

/-- The complete three-index Kraus family is normalized. -/
theorem ProjectivePartition.partitionAverageKraus_complete
    (part : ProjectivePartition n k) :
    (∑ i : Fin k, ∑ a : Fin n, ∑ b : Fin n,
      (partitionAverageKraus part i a b)ᴴ *
        partitionAverageKraus part i a b) = 1 := by
  rw [Finset.sum_congr rfl fun i _ =>
    partitionAverageKraus_block_complete part i]
  exact part.complete

/-- One block of the Kraus sum is exactly the corresponding term in the
partition average. -/
theorem partitionAverageKraus_block_form
    (part : ProjectivePartition n k) (i : Fin k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    (∑ a : Fin n, ∑ b : Fin n,
      partitionAverageKraus part i a b * X *
        (partitionAverageKraus part i a b)ᴴ) =
      (((part.proj i).trace)⁻¹ *
        bornWeight X (part.proj i)) • part.proj i := by
  by_cases htrace : (part.proj i).trace = 0
  · have hP : part.proj i = 0 :=
      ((part.isEvent i).posSemidef.trace_eq_zero_iff).mp htrace
    simp [partitionAverageKraus, hP, bornWeight]
  · have hscale := star_partitionBlockScale_mul_self part i
    have hstar : star (partitionBlockScale part i) =
        partitionBlockScale part i := by
      simp [partitionBlockScale]
    have hmul : partitionBlockScale part i * partitionBlockScale part i =
        ((part.proj i).trace)⁻¹ := by
      calc
        partitionBlockScale part i * partitionBlockScale part i =
            star (partitionBlockScale part i) *
              partitionBlockScale part i := by rw [hstar]
        _ = ((part.proj i).trace)⁻¹ := hscale
    have hterm : ∀ a b : Fin n,
        partitionAverageKraus part i a b * X *
            (partitionAverageKraus part i a b)ᴴ =
          ((part.proj i).trace)⁻¹ •
            (part.proj i * matrixUnit a b * part.proj i * X *
              part.proj i * matrixUnit b a * part.proj i) := by
      intro a b
      rw [conjTranspose_partitionAverageKraus]
      simp only [partitionAverageKraus, smul_mul_assoc, mul_smul_comm,
        smul_smul]
      rw [hmul]
      congr 1
      simp only [mul_assoc]
    rw [Finset.sum_congr rfl fun a _ =>
      Finset.sum_congr rfl fun b _ => hterm a b]
    have hsandwich :
        (∑ a : Fin n, ∑ b : Fin n,
          part.proj i * matrixUnit a b * part.proj i * X *
            part.proj i * matrixUnit b a * part.proj i) =
          bornWeight X (part.proj i) • part.proj i := by
      calc
        (∑ a : Fin n, ∑ b : Fin n,
          part.proj i * matrixUnit a b * part.proj i * X *
            part.proj i * matrixUnit b a * part.proj i) =
            part.proj i *
              (∑ a : Fin n, ∑ b : Fin n,
                matrixUnit a b * (part.proj i * X * part.proj i) *
                  matrixUnit b a) * part.proj i := by
              simp only [Finset.mul_sum, Finset.sum_mul, mul_assoc]
        _ = part.proj i *
              ((part.proj i * X * part.proj i).trace •
                (1 : Matrix (Fin n) (Fin n) ℂ)) * part.proj i := by
              rw [sum_matrixUnit_mul_mul_swap]
        _ = (part.proj i * X * part.proj i).trace • part.proj i := by
              simp [(part.isEvent i).2]
        _ = bornWeight X (part.proj i) • part.proj i := by
              rw [trace_sandwich (part.isEvent i).2 X]
    calc
      (∑ a : Fin n, ∑ b : Fin n,
        (part.proj i).trace⁻¹ •
          (part.proj i * matrixUnit a b * part.proj i * X *
            part.proj i * matrixUnit b a * part.proj i)) =
          (part.proj i).trace⁻¹ •
            (∑ a : Fin n, ∑ b : Fin n,
              part.proj i * matrixUnit a b * part.proj i * X *
                part.proj i * matrixUnit b a * part.proj i) := by
            rw [Finset.smul_sum]
            apply Finset.sum_congr rfl
            intro a _
            rw [Finset.smul_sum]
      _ = (((part.proj i).trace)⁻¹ *
          bornWeight X (part.proj i)) • part.proj i := by
        rw [hsandwich, smul_smul]

/-- Explicit Kraus form of partition averaging. -/
theorem partitionAverage_kraus_form
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part X =
      ∑ i : Fin k, ∑ a : Fin n, ∑ b : Fin n,
        partitionAverageKraus part i a b * X *
          (partitionAverageKraus part i a b)ᴴ := by
  rw [partitionAverage]
  symm
  exact Finset.sum_congr rfl fun i _ =>
    partitionAverageKraus_block_form part i X

end EventAlgebra

#print axioms EventAlgebra.ProjectivePartition.partitionAverageKraus_complete
#print axioms EventAlgebra.partitionAverage_kraus_form
#print axioms EventAlgebra.trace_partitionAverage
