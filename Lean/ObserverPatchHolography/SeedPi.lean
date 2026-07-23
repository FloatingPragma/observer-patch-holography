import Mathlib
import ObserverPatchHolography.CapacityFixedPoint

/-!
# Coupling-Algebra Part A: CAP-P Seed

Abstract form of the CAP-P seed statement. The multiplicative CAP-P
readback branches in `code/capacity_readback/F_candidate_capP.py` are

```text
F(N) = π · (N/π)^s = π · exp(s · log(N/π)),   s ∈ (0, 1) per branch,
```

and the Python certificate records `N = π` as the unique fixed point
("the unique solution of `(N/π)^s = N/π` with `s ≠ 1`"). This module proves
that exact-form claim:

* `capReadback_fixedPt_iff` — for `s ≠ 1` and `N > 0`,
  `F(N) = N ↔ N = π`.
* `capReadback_log_coords` — in log coordinates `y = log(N/π)` the readback
  is the linear map `y ↦ s·y`, i.e. the `CapacityFixedPoint.averagingMap`
  with `λ = 1 − s` toward target `0` (`capReadback_log_coords_averaging`).
* `linear_branch_no_positive_fixedPt` — the additive CAP-P branches
  `F(N) = s·N` with `s < 1` have no positive fixed point.

**Scope warning.** ALGEBRAIC layer only: the map shape is taken from the
Python certificate and the survival factors `s` (e.g. `e^(−P/24)`,
`1 − P/24`) are abstracted to a real parameter; the certified numeric
interval for `s`, the contraction certificates, and the physical identities
I1/I2 stay outside the formalised set.
-/

namespace OPH.SeedPi

open Function

/-- The CAP-P multiplicative readback map
    `F(N) = π · exp(s · log(N/π)) = π · (N/π)^s`, with the survival factor
    `s` abstracted to a real parameter (matching
    `F_candidate_capP.py`: `F = lambda x: pi * exp(s * log(x / pi))`). -/
noncomputable def capReadback (s N : ℝ) : ℝ :=
  Real.pi * Real.exp (s * Real.log (N / Real.pi))

/-- In log coordinates `y = log(N/π)` the CAP-P readback is the linear map
    `y ↦ s·y`. -/
theorem capReadback_log_coords (s N : ℝ) :
    Real.log (capReadback s N / Real.pi) = s * Real.log (N / Real.pi) := by
  unfold capReadback
  rw [mul_div_cancel_left₀ _ Real.pi_pos.ne', Real.log_exp]

/-- The log-coordinate readback is the `averagingMap` schema of
    `CapacityFixedPoint` with `λ = 1 − s` toward target `0`. -/
theorem capReadback_log_coords_averaging (s N : ℝ) :
    Real.log (capReadback s N / Real.pi) =
      OPH.CapacityFixedPoint.averagingMap (1 - s) 0 (Real.log (N / Real.pi)) := by
  rw [capReadback_log_coords]
  unfold OPH.CapacityFixedPoint.averagingMap
  ring

/-- **CAP-P seed.** For `s ≠ 1` and `N > 0`, the readback fixes `N` iff
    `N = π`: the unique positive fixed point of every multiplicative CAP-P
    branch is the seed `π`. -/
theorem capReadback_fixedPt_iff {s N : ℝ} (hs : s ≠ 1) (hN : 0 < N) :
    capReadback s N = N ↔ N = Real.pi := by
  have hpi : (0 : ℝ) < Real.pi := Real.pi_pos
  have hdivpos : 0 < N / Real.pi := div_pos hN hpi
  constructor
  · intro h
    have hlog : Real.log (N / Real.pi) = s * Real.log (N / Real.pi) := by
      have hcoords := capReadback_log_coords s N
      rw [h] at hcoords
      exact hcoords
    have hzero : Real.log (N / Real.pi) = 0 := by
      have hfac : (s - 1) * Real.log (N / Real.pi) = 0 := by
        linear_combination -hlog
      rcases mul_eq_zero.mp hfac with h' | h'
      · exact absurd (by linarith : s = 1) hs
      · exact h'
    have hone : N / Real.pi = 1 :=
      Real.eq_one_of_pos_of_log_eq_zero hdivpos hzero
    field_simp at hone
    exact hone
  · rintro rfl
    unfold capReadback
    rw [div_self hpi.ne', Real.log_one, mul_zero, Real.exp_zero, mul_one]

/-- The additive CAP-P branches `F(N) = s·N` with `s < 1` have no positive
    fixed point (their unique fixed point `0` lies outside the admissible
    interval). -/
theorem linear_branch_no_positive_fixedPt {s N : ℝ} (hs : s < 1) (hN : 0 < N) :
    s * N ≠ N := by
  have hlt : s * N < N := by nlinarith
  exact hlt.ne

-- Axiom audit: these must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms capReadback_log_coords
#print axioms capReadback_log_coords_averaging
#print axioms capReadback_fixedPt_iff
#print axioms linear_branch_no_positive_fixedPt

end OPH.SeedPi
