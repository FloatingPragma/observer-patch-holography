import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Algebra.Group.Commute.Defs
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic

/-!
# Finite periodic momentum conservation

The exact CAR attachment in `model/periodic_momentum` supplies commuting
graded translations of the periodic Weyl/contact action. Their finite
Fourier spectral logarithms are scalar polynomials in those translations.
The theorems below transport conservation to these polynomials and record
the integer-coefficient contradiction used after prime cyclotomic reduction.

The CAR source, spectral branch construction, prime minimal-polynomial
argument and local one-particle support reduction are separately proved in
the package. No local stress tensor or Einstein premise is assumed here.
-/

namespace OPH.QFT.PeriodicMomentum

open scoped BigOperators

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem polynomial_conserved (H T : Matrix ι ι ℂ) (h : Commute H T)
    (c : ℕ → ℂ) (n : ℕ) :
    H * (∑ r ∈ Finset.range n, c r • T ^ r) =
      (∑ r ∈ Finset.range n, c r • T ^ r) * H := by
  rw [Finset.mul_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro r _
  rw [Matrix.mul_smul, Matrix.smul_mul, (h.pow_right r).eq]

theorem polynomials_commute (T S : Matrix ι ι ℂ) (h : Commute T S)
    (c d : ℕ → ℂ) (n k : ℕ) :
    (∑ r ∈ Finset.range n, c r • T ^ r) *
      (∑ u ∈ Finset.range k, d u • S ^ u) =
    (∑ u ∈ Finset.range k, d u • S ^ u) *
      (∑ r ∈ Finset.range n, c r • T ^ r) := by
  have hs : Commute S (∑ r ∈ Finset.range n, c r • T ^ r) :=
    polynomial_conserved S T h.symm c n
  exact polynomial_conserved _ S hs.symm d k

/-- Distinct modulo-period spin-trace coefficients cannot become equal by
changing arbitrary integer branches. Prime cyclotomic reduction supplies
the forbidden equality; this elementary conclusion needs only period>2. -/
theorem spin_trace_branches_differ (period a b : ℤ) (hperiod : 2 < period) :
    period * a ≠ 2 + period * b := by
  intro h
  have hm := congrArg (fun z : ℤ => z % period) h
  have two : (2 : ℤ) % period = 2 := Int.emod_eq_of_lt (by omega) hperiod
  simp only [Int.add_emod, Int.mul_emod, Int.emod_self, zero_mul,
    Int.zero_emod, add_zero, two] at hm
  omega

#print axioms polynomial_conserved
#print axioms polynomials_commute
#print axioms spin_trace_branches_differ

end OPH.QFT.PeriodicMomentum
