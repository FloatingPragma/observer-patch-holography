import Mathlib

namespace OPH.RGRepresentationFrontier

/-!
# Exact representation-index frontier for the RG lane

This file checks the rational arithmetic used by the non-promoting issue #32
frontier.  The finite matter certificate supplies one rank-fifteen chiral
matter copy.  The standard four-dimensional one-loop gauge-beta formula is a
separate imported QFT premise.  Family multiplicity, scalar multiplicity,
continuum kinetic normalization, thresholds, and renormalization schemes are
not selected here.

The sign convention is

  d g_i / d log(mu) = b_i g_i^3 / (16 pi^2).

The opposite asymptotic-freedom convention is `B_i = -b_i`.
-/

/-- Sum of color Dynkin indices over the Weyl states in one matter copy. -/
def su3WeylIndexPerCopy : ℚ :=
  2 * (1 / 2) + 1 / 2 + 1 / 2

/-- Sum of weak Dynkin indices over the Weyl states in one matter copy. -/
def su2WeylIndexPerCopy : ℚ :=
  3 * (1 / 2) + 1 / 2

/-- Hypercharge-squared index over the Weyl states in one matter copy. -/
def u1WeylIndexPerCopy : ℚ :=
  6 * (1 / 6) ^ 2 +
    3 * (-2 / 3) ^ 2 +
    3 * (1 / 3) ^ 2 +
    2 * (-1 / 2) ^ 2 +
    1 ^ 2

/-- One complex color-singlet weak doublet contributes this U(1) index. -/
def u1ScalarIndexPerDoublet : ℚ :=
  2 * (1 / 2) ^ 2

/-- One complex color-singlet weak doublet contributes this SU(2) index. -/
def su2ScalarIndexPerDoublet : ℚ :=
  1 / 2

theorem representation_indices :
    su3WeylIndexPerCopy = 2 ∧
      su2WeylIndexPerCopy = 2 ∧
      u1WeylIndexPerCopy = 10 / 3 ∧
      u1ScalarIndexPerDoublet = 1 / 2 ∧
      su2ScalarIndexPerDoublet = 1 / 2 := by
  norm_num [su3WeylIndexPerCopy, su2WeylIndexPerCopy,
    u1WeylIndexPerCopy, u1ScalarIndexPerDoublet,
    su2ScalarIndexPerDoublet]

/-- Hypercharge coefficient after importing the standard one-loop QFT law. -/
def bY (nG nH : ℚ) : ℚ :=
  (20 / 9) * nG + (1 / 6) * nH

/-- Weak coefficient after importing the standard one-loop QFT law. -/
def b2 (nG nH : ℚ) : ℚ :=
  -22 / 3 + (4 / 3) * nG + (1 / 6) * nH

/-- Color coefficient after importing the standard one-loop QFT law. -/
def b3 (nG : ℚ) : ℚ :=
  -11 + (4 / 3) * nG

/-- The imported one-loop gauge functional written directly in terms of
    fermion and scalar representation indices. -/
def oneLoopGaugeCoefficient (adjointCasimir fermionIndex scalarIndex : ℚ) : ℚ :=
  -(11 / 3) * adjointCasimir +
    (2 / 3) * fermionIndex +
    (1 / 3) * scalarIndex

/-- Conditional arithmetic at the declared `(nG,nH)=(3,1)` completion. -/
theorem declared_three_one_evaluation :
    bY 3 1 = 41 / 6 ∧
      b2 3 1 = -19 / 6 ∧
      b3 3 = -7 := by
  norm_num [bY, b2, b3]

/-- Adding one matter copy shifts all three coefficients by exact amounts. -/
theorem family_copy_shift (nG nH : ℚ) :
    bY (nG + 1) nH - bY nG nH = 20 / 9 ∧
      b2 (nG + 1) nH - b2 nG nH = 4 / 3 ∧
      b3 (nG + 1) - b3 nG = 4 / 3 := by
  constructor
  · dsimp [bY]
    ring
  constructor
  · dsimp [b2]
    ring
  · dsimp [b3]
    ring

/-- Adding one scalar doublet shifts only the U(1) and SU(2) coefficients. -/
theorem scalar_copy_shift (nG nH : ℚ) :
    bY nG (nH + 1) - bY nG nH = 1 / 6 ∧
      b2 nG (nH + 1) - b2 nG nH = 1 / 6 ∧
      b3 nG - b3 nG = 0 := by
  constructor
  · dsimp [bY]
    ring
  constructor
  · dsimp [b2]
    ring
  · ring

/-- A direct summand with zero gauge-representation indices cannot change a
    one-loop gauge coefficient.  This statement does not establish that a
    source-invisible summand has no Yukawa, mass-mixing, or scalar vertex. -/
theorem zero_gauge_index_direct_sum_invariance
    (adjointCasimir fermionIndex scalarIndex : ℚ) :
    oneLoopGaugeCoefficient adjointCasimir (fermionIndex + 0) (scalarIndex + 0) =
      oneLoopGaugeCoefficient adjointCasimir fermionIndex scalarIndex := by
  simp

#print axioms representation_indices
#print axioms declared_three_one_evaluation
#print axioms family_copy_shift
#print axioms scalar_copy_shift
#print axioms zero_gauge_index_direct_sum_invariance

end OPH.RGRepresentationFrontier
