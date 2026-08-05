import EventAlgebra.PartitionAverage

/-!
# Operational invisibility of partition-off-diagonal data

A declared projective partition supplies two finite readouts.  Partition
pinching retains the whole sector-preserving commutant, while partition
averaging retains only the commutative span of the sector projectors.  This
module proves the exact operational quotient for the larger readout: two
matrices have identical trace statistics against every matrix in the
sector-preserving commutant if and only if their pinchings agree. It also proves that a
corner `P_i X P_j` with `i ≠ j` is erased by both readouts.

This is an exact finite superselection statement relative to the supplied
partition and readout.  It does not construct an edge sector, select a
physical partition, or identify a laboratory measurement.
-/

namespace EventAlgebra

open Matrix

variable {n k : ℕ}

/-- Two matrices are operationally equivalent relative to a supplied
partition when every matrix in its sector-preserving commutant gives the same
trace statistic. -/
def PartitionOperationallyEquivalent (part : ProjectivePartition n k)
    (X Y : Matrix (Fin n) (Fin n) ℂ) : Prop :=
  ∀ C : Matrix (Fin n) (Fin n) ℂ, C ∈ part.commutant →
    (X * C).trace = (Y * C).trace

/-- **Trace-dependent.** Equality after pinching is exactly operational
equivalence against the complete sector-preserving commutant.  Thus the
fibres of `partitionPinching` are the operational equivalence classes for
the declared partition. -/
theorem partitionOperationallyEquivalent_iff_pinching_eq
    (part : ProjectivePartition n k)
    (X Y : Matrix (Fin n) (Fin n) ℂ) :
    PartitionOperationallyEquivalent part X Y ↔
      partitionPinching part X = partitionPinching part Y := by
  constructor
  · intro h
    have hYX : partitionPinching part Y = partitionPinching part X :=
      partitionPinching_unique part
        ((ProjectivePartition.mem_commutant_iff part).mp
          (partitionPinching_mem_commutant part Y))
        (fun C hC ↦ by
          calc
            (partitionPinching part Y * C).trace = (Y * C).trace :=
              trace_partitionPinching_mul_commutant part Y C hC
            _ = (X * C).trace :=
              (h C ((ProjectivePartition.mem_commutant_iff part).mpr hC)).symm)
    exact hYX.symm
  · intro h C hC
    have hC' := (ProjectivePartition.mem_commutant_iff part).mp hC
    calc
      (X * C).trace = (partitionPinching part X * C).trace :=
        (trace_partitionPinching_mul_commutant part X C hC').symm
      _ = (partitionPinching part Y * C).trace := by rw [h]
      _ = (Y * C).trace :=
        trace_partitionPinching_mul_commutant part Y C hC'

/-- A matrix is wholly sector-off-diagonal relative to a supplied partition
when its sector-preserving pinching vanishes. -/
def IsPartitionOffDiagonal (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) : Prop :=
  partitionPinching part X = 0

/-- **Trace-dependent.** Every wholly sector-off-diagonal matrix is invisible
to every matrix in the sector-preserving commutant, not only to the central projector
span. -/
theorem trace_mul_eq_zero_of_partitionOffDiagonal
    (part : ProjectivePartition n k) {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : IsPartitionOffDiagonal part X)
    (C : Matrix (Fin n) (Fin n) ℂ) (hC : C ∈ part.commutant) :
    (X * C).trace = 0 := by
  have hC' := (ProjectivePartition.mem_commutant_iff part).mp hC
  rw [← trace_partitionPinching_mul_commutant part X C hC', hX,
    zero_mul, trace_zero]

/-- **Trace-dependent.** Adding wholly sector-off-diagonal data leaves every
sector-preserving trace statistic unchanged. -/
theorem partitionOperationallyEquivalent_add_of_partitionOffDiagonal
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) {D : Matrix (Fin n) (Fin n) ℂ}
    (hD : IsPartitionOffDiagonal part D) :
    PartitionOperationallyEquivalent part (X + D) X := by
  rw [partitionOperationallyEquivalent_iff_pinching_eq]
  change partitionPinchingLinearMap part (X + D) = partitionPinching part X
  rw [(partitionPinchingLinearMap part).map_add, partitionPinchingLinearMap_apply,
    partitionPinchingLinearMap_apply, hD, add_zero]

/-- The `(i,j)` corner of `X` relative to a projective partition. -/
noncomputable def partitionCorner (part : ProjectivePartition n k)
    (i j : Fin k) (X : Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ :=
  part.proj i * X * part.proj j

/-- **Trace-dependent.** A sector-off-diagonal corner has zero trace pairing
with each partition projector. -/
theorem bornWeight_partitionCorner_proj_eq_zero
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X : Matrix (Fin n) (Fin n) ℂ) (l : Fin k) :
    bornWeight (partitionCorner part i j X) (part.proj l) = 0 := by
  unfold bornWeight partitionCorner
  by_cases hjl : j = l
  · subst l
    calc
      ((part.proj i * X * part.proj j) * part.proj j).trace =
          ((part.proj i * X) * part.proj j).trace := by
        rw [mul_assoc, (part.isEvent j).2]
      _ = (part.proj j * (part.proj i * X)).trace :=
        trace_mul_comm (part.proj i * X) (part.proj j)
      _ = ((part.proj j * part.proj i) * X).trace := by
        rw [mul_assoc]
      _ = 0 := by
        rw [part.proj_mul_proj, if_neg (Ne.symm hij), zero_mul, trace_zero]
  · rw [show (part.proj i * X * part.proj j) * part.proj l =
        part.proj i * X * (part.proj j * part.proj l) by
          simp only [mul_assoc]]
    rw [part.proj_mul_proj, if_neg hjl, mul_zero, trace_zero]

/-- **Algebra-only.** Sector pinching erases every cross-sector corner
exactly. -/
theorem partitionPinching_partitionCorner_eq_zero
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionPinching part (partitionCorner part i j X) = 0 := by
  unfold partitionPinching partitionCorner
  apply Finset.sum_eq_zero
  intro l _
  by_cases hli : l = i
  · subst l
    rw [show part.proj i * (part.proj i * X * part.proj j) * part.proj i =
        ((part.proj i * part.proj i) * X) *
          (part.proj j * part.proj i) by simp only [mul_assoc],
      (part.isEvent i).2, part.orthogonal j i (Ne.symm hij), mul_zero]
  · rw [show part.proj l * (part.proj i * X * part.proj j) * part.proj l =
        (((part.proj l * part.proj i) * X) * part.proj j) * part.proj l by
          simp only [mul_assoc],
      part.orthogonal l i hli, zero_mul, zero_mul, zero_mul]

/-- **Trace-dependent.** The declared public partition average erases every
sector-off-diagonal corner exactly. -/
theorem partitionAverage_partitionCorner_eq_zero
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (partitionCorner part i j X) = 0 := by
  unfold partitionAverage
  apply Finset.sum_eq_zero
  intro l _
  rw [bornWeight_partitionCorner_proj_eq_zero part hij X l, mul_zero,
    zero_smul]

/-- **Trace-dependent.** Every public test matrix in the partition span is
operationally insensitive to sector-off-diagonal data. -/
theorem trace_partitionCorner_mul_eq_zero_of_mem_span
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X C : Matrix (Fin n) (Fin n) ℂ) (hC : C ∈ part.span) :
    (partitionCorner part i j X * C).trace = 0 := by
  have hDual := trace_partitionAverage_mul_of_mem part
    (partitionCorner part i j X) hC
  rw [partitionAverage_partitionCorner_eq_zero part hij X, zero_mul,
    trace_zero] at hDual
  exact hDual.symm

/-- **Trace-dependent.** Adding a sector-off-diagonal corner does not change
the declared public average. -/
theorem partitionAverage_add_partitionCorner
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X Y : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (X + partitionCorner part i j Y) =
      partitionAverage part X := by
  rw [partitionAverage_add,
    partitionAverage_partitionCorner_eq_zero part hij Y, add_zero]

/-- **Trace-dependent.** Adding a sector-off-diagonal corner leaves every
public trace statistic unchanged. -/
theorem trace_add_partitionCorner_mul_of_mem_span
    (part : ProjectivePartition n k) {i j : Fin k} (hij : i ≠ j)
    (X Y C : Matrix (Fin n) (Fin n) ℂ) (hC : C ∈ part.span) :
    ((X + partitionCorner part i j Y) * C).trace = (X * C).trace := by
  rw [add_mul, trace_add,
    trace_partitionCorner_mul_eq_zero_of_mem_span part hij Y C hC, add_zero]

#print axioms EventAlgebra.bornWeight_partitionCorner_proj_eq_zero
#print axioms EventAlgebra.partitionOperationallyEquivalent_iff_pinching_eq
#print axioms EventAlgebra.trace_mul_eq_zero_of_partitionOffDiagonal
#print axioms EventAlgebra.partitionOperationallyEquivalent_add_of_partitionOffDiagonal
#print axioms EventAlgebra.partitionPinching_partitionCorner_eq_zero
#print axioms EventAlgebra.partitionAverage_partitionCorner_eq_zero
#print axioms EventAlgebra.trace_partitionCorner_mul_eq_zero_of_mem_span
#print axioms EventAlgebra.partitionAverage_add_partitionCorner
#print axioms EventAlgebra.trace_add_partitionCorner_mul_of_mem_span

end EventAlgebra
