import Mathlib

namespace ObserverPatchHolography.DirectedSeamRepair

/-!
# Directed integer seam balancing

This file isolates the exact local arithmetic of the integer lift of scalar
pair averaging.  A pair with total `s` is sent to the two nearest integers
around `s / 2`.  The two possible directions differ only in which endpoint
receives the upper value.

The results are explicitly conditioned on a fixed total sector.  When the
total is odd, neither directed outcome lies in the exact endpoint-agreement
equalizer.  Averaging the two directions after casting to `ℚ` does recover
the rational pair average.

Nothing in this file selects a seam schedule, proves global convergence, or
identifies this integer kernel with a physical repair law.
-/

/-- The lower member of the nearest-balanced pair with total `a + b`. -/
def lowerLoad (a b : ℤ) : ℤ :=
  (a + b) / 2

/-- The upper member, defined so total preservation is definitional. -/
def upperLoad (a b : ℤ) : ℤ :=
  a + b - lowerLoad a b

/-- The directed outcome which places the lower load at the left endpoint. -/
def balanceForward (a b : ℤ) : ℤ × ℤ :=
  (lowerLoad a b, upperLoad a b)

/-- The reverse directed outcome. -/
def balanceReverse (a b : ℤ) : ℤ × ℤ :=
  (upperLoad a b, lowerLoad a b)

/-- The conserved endpoint total. -/
def pairTotal (q : ℤ × ℤ) : ℤ :=
  q.1 + q.2

/-- A fixed-total sector.  Repair does not identify different sectors. -/
def TotalSector (s : ℤ) : Set (ℤ × ℤ) :=
  {q | pairTotal q = s}

/-- The local quadratic load used to measure strict balancing descent. -/
def quadraticLoad (q : ℤ × ℤ) : ℤ :=
  q.1 ^ 2 + q.2 ^ 2

@[simp]
theorem lower_add_upper (a b : ℤ) :
    lowerLoad a b + upperLoad a b = a + b := by
  simp [upperLoad]

@[simp]
theorem pairTotal_forward (a b : ℤ) :
    pairTotal (balanceForward a b) = a + b := by
  simp [pairTotal, balanceForward]

@[simp]
theorem pairTotal_reverse (a b : ℤ) :
    pairTotal (balanceReverse a b) = a + b := by
  simp [pairTotal, balanceReverse, add_comm]

/-- Both directed outcomes remain in the sector of the input pair. -/
theorem directed_outcomes_preserve_sector (a b : ℤ) :
    balanceForward a b ∈ TotalSector (a + b) ∧
      balanceReverse a b ∈ TotalSector (a + b) := by
  simp [TotalSector]

/-- Sector preservation stated with a separately supplied sector label. -/
theorem directed_outcomes_preserve_named_sector
    {a b s : ℤ} (hsector : (a, b) ∈ TotalSector s) :
    balanceForward a b ∈ TotalSector s ∧
      balanceReverse a b ∈ TotalSector s := by
  have htotal : a + b = s := hsector
  simpa [htotal] using directed_outcomes_preserve_sector a b

/-- The output mismatch is exactly the Euclidean remainder modulo two. -/
theorem upper_sub_lower_eq_emod (a b : ℤ) :
    upperLoad a b - lowerLoad a b = (a + b) % 2 := by
  have h := Int.emod_add_mul_ediv (a + b) 2
  simp only [OfNat.ofNat] at h
  dsimp [upperLoad, lowerLoad]
  omega

/-- The nearest-balanced directed outcomes have residual gap zero or one. -/
theorem output_gap_zero_or_one (a b : ℤ) :
    upperLoad a b - lowerLoad a b = 0 ∨
      upperLoad a b - lowerLoad a b = 1 := by
  rw [upper_sub_lower_eq_emod]
  exact Int.emod_two_eq_zero_or_one (a + b)

/-- Even total gives exact endpoint agreement; odd total leaves gap one. -/
theorem output_agreement_or_unit_residual (a b : ℤ) :
    lowerLoad a b = upperLoad a b ∨
      upperLoad a b = lowerLoad a b + 1 := by
  rcases output_gap_zero_or_one a b with h | h
  · left
    omega
  · right
    omega

/-- Averaging the two directions gives the rational pair average at the
left endpoint. -/
theorem directed_mean_left_eq_rational_pairAverage (a b : ℤ) :
    (((balanceForward a b).1 : ℚ) + (balanceReverse a b).1) / 2 =
      ((a : ℚ) + (b : ℚ)) / 2 := by
  norm_num [balanceForward, balanceReverse]
  have h := lower_add_upper a b
  exact_mod_cast h

/-- Averaging the two directions gives the rational pair average at the
right endpoint as well. -/
theorem directed_mean_right_eq_rational_pairAverage (a b : ℤ) :
    (((balanceForward a b).2 : ℚ) + (balanceReverse a b).2) / 2 =
      ((a : ℚ) + (b : ℚ)) / 2 := by
  norm_num [balanceForward, balanceReverse, add_comm]
  have h := lower_add_upper a b
  exact_mod_cast h

/-- The exact polarization identity for the local quadratic drop.  It needs
only conservation of the pair total and therefore involves no approximation
or asymptotic argument. -/
theorem two_mul_quadratic_drop (a b : ℤ) :
    2 * (quadraticLoad (a, b) -
      quadraticLoad (balanceForward a b)) =
      (a - b) ^ 2 - (upperLoad a b - lowerLoad a b) ^ 2 := by
  simp only [quadraticLoad, balanceForward]
  simp [upperLoad]
  ring

/-- Whenever the incoming squared mismatch exceeds one, nearest balancing
strictly lowers the quadratic load. -/
theorem quadraticLoad_forward_strict (a b : ℤ)
    (hgap : 1 < (a - b) ^ 2) :
    quadraticLoad (balanceForward a b) < quadraticLoad (a, b) := by
  have hdrop := two_mul_quadratic_drop a b
  rcases output_gap_zero_or_one a b with hzero | hone
  · have hsquare :
        (upperLoad a b - lowerLoad a b) ^ 2 = 0 := by
      rw [hzero]
      norm_num
    rw [hsquare] at hdrop
    nlinarith
  · have hsquare :
        (upperLoad a b - lowerLoad a b) ^ 2 = 1 := by
      rw [hone]
      norm_num
    rw [hsquare] at hdrop
    nlinarith

/-- Reversing the direction does not change the quadratic load. -/
theorem quadraticLoad_reverse_eq_forward (a b : ℤ) :
    quadraticLoad (balanceReverse a b) =
      quadraticLoad (balanceForward a b) := by
  simp [quadraticLoad, balanceForward, balanceReverse, add_comm]

/-- The same strict descent holds for the reverse directed outcome. -/
theorem quadraticLoad_reverse_strict (a b : ℤ)
    (hgap : 1 < (a - b) ^ 2) :
    quadraticLoad (balanceReverse a b) < quadraticLoad (a, b) := by
  rw [quadraticLoad_reverse_eq_forward]
  exact quadraticLoad_forward_strict a b hgap

end ObserverPatchHolography.DirectedSeamRepair

/- Axiom audit: exact integer and rational arithmetic with standard axioms. -/

#print axioms ObserverPatchHolography.DirectedSeamRepair.directed_outcomes_preserve_sector
#print axioms ObserverPatchHolography.DirectedSeamRepair.directed_mean_left_eq_rational_pairAverage
#print axioms ObserverPatchHolography.DirectedSeamRepair.quadraticLoad_forward_strict
