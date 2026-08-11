import Variational.DiscreteEulerLagrange

/-!
# Stationary histories outside the Gibbs coverage

Issue B7 (#683).  The committed history projection concentrates on
least-action histories, and the discrete Euler-Lagrange theorem derives
stationarity from minimality.  This module certifies the scoped negative
for the converse coverage question: a history can satisfy the discrete
Euler-Lagrange stationarity condition at its interior junction, with
explicit derivative witnesses, while failing minimality outright, and
the Gibbs weight at every positive multiplier strictly prefers a
non-stationary variation to it.  Stationary-point coverage by the Gibbs
projection therefore fails in general: the projection tracks minimizers
and provably misses this stationary history.

The witness is the three-record zero history under the concave
two-point Lagrangian `L a b = -((b - a)^2) / 2`.  The single-site
variation functional at the interior record is `-x^2`, so the zero
history is a strict local maximum of the variation: stationary with
vanishing adjacent partials, beaten by every nonzero variation, and
carrying strictly smaller Gibbs weight than the unit variation at every
positive multiplier.

**Boundary.**  This is a finite counterexample about the committed
variational interface; no physical action, saddle-point principle,
time, amplitude, or continuum statement enters.  Together with the
committed finite-variation obstruction and the Legendre and multiplier
nonidentifiability receipts, it is the scoped negative exit for the
stationary-saddle item of issue #683, which is closed bounded.
-/

namespace OPH.Variational

noncomputable section

/-- The concave two-point Lagrangian. -/
def saddleL : ℝ → ℝ → ℝ := fun a b => -((b - a) ^ 2) / 2

/-- The three-record zero history. -/
def saddleHistory : Fin 3 → ℝ := fun _ => 0

/-- The zero history is stationary: the two adjacent partial
derivatives at the interior junction exist and sum to zero, in exactly
the form the discrete Euler-Lagrange theorem produces. -/
theorem saddleHistory_stationary :
    ∃ d₂ d₁ : ℝ,
      HasDerivAt (fun x => saddleL (saddleHistory (0 : Fin 2).castSucc) x)
        d₂ (saddleHistory (0 : Fin 2).succ) ∧
      HasDerivAt (fun x => saddleL x (saddleHistory (1 : Fin 2).succ))
        d₁ (saddleHistory (1 : Fin 2).castSucc) ∧
      d₂ + d₁ = 0 := by
  refine ⟨0, 0, ?_, ?_, by norm_num⟩
  · have h : (fun x => saddleL (saddleHistory (0 : Fin 2).castSucc) x) =
        fun x => -(x ^ 2) / 2 := by
      funext x
      simp [saddleL, saddleHistory]
    rw [h]
    have hd : HasDerivAt (fun x : ℝ => -(x ^ 2) / 2) (-(2 * 0) / 2) 0 := by
      have := ((hasDerivAt_pow 2 (0 : ℝ)).neg).div_const 2
      simpa using this
    simpa [saddleHistory] using hd
  · have h : (fun x => saddleL x (saddleHistory (1 : Fin 2).succ)) =
        fun x => -(x ^ 2) / 2 := by
      funext x
      simp [saddleL, saddleHistory]
    rw [h]
    have hd : HasDerivAt (fun x : ℝ => -(x ^ 2) / 2) (-(2 * 0) / 2) 0 := by
      have := ((hasDerivAt_pow 2 (0 : ℝ)).neg).div_const 2
      simpa using this
    simpa [saddleHistory] using hd

/-- The zero history is not a single-site minimizer: the unit variation
strictly lowers the local action. -/
theorem saddleHistory_not_minimal :
    ¬ ∀ x, localAction saddleL saddleHistory ≤
      localAction saddleL
        (Function.update saddleHistory ((0 : Fin 2).succ) x) := by
  intro h
  have h1 := h 1
  have hbase : localAction saddleL saddleHistory = 0 := by
    unfold localAction saddleL saddleHistory
    simp
  have hupd : localAction saddleL
      (Function.update saddleHistory ((0 : Fin 2).succ) 1) = -1 := by
    unfold localAction
    rw [Fin.sum_univ_two]
    have e0 : ((0 : Fin 2).castSucc : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have e2 : ((1 : Fin 2).succ : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have h00 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).castSucc) = 0 := by
      rw [Function.update_of_ne e0]
      rfl
    have h11 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).succ) = 1 := Function.update_self _ _ _
    have h1c : ((1 : Fin 2).castSucc : Fin 3) = (0 : Fin 2).succ := by decide
    have h1m : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).castSucc) = 1 := by
      rw [h1c]
      exact Function.update_self _ _ _
    have h22 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).succ) = 0 := by
      rw [Function.update_of_ne e2]
      rfl
    rw [h00, h11, h1m, h22]
    unfold saddleL
    norm_num
  rw [hbase, hupd] at h1
  norm_num at h1

/-- At every positive multiplier the Gibbs weight strictly prefers the
non-stationary unit variation to the stationary history: the projection
tracks minimizers and misses the stationary point. -/
theorem gibbs_prefers_nonstationary {lam : ℝ} (hlam : 0 < lam) :
    Real.exp (-lam * localAction saddleL saddleHistory) <
      Real.exp (-lam * localAction saddleL
        (Function.update saddleHistory ((0 : Fin 2).succ) 1)) := by
  have hbase : localAction saddleL saddleHistory = 0 := by
    unfold localAction saddleL saddleHistory
    simp
  have hupd : localAction saddleL
      (Function.update saddleHistory ((0 : Fin 2).succ) 1) = -1 := by
    unfold localAction
    rw [Fin.sum_univ_two]
    have e0 : ((0 : Fin 2).castSucc : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have e2 : ((1 : Fin 2).succ : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have h00 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).castSucc) = 0 := by
      rw [Function.update_of_ne e0]
      rfl
    have h11 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).succ) = 1 := Function.update_self _ _ _
    have h1c : ((1 : Fin 2).castSucc : Fin 3) = (0 : Fin 2).succ := by decide
    have h1m : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).castSucc) = 1 := by
      rw [h1c]
      exact Function.update_self _ _ _
    have h22 : Function.update saddleHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).succ) = 0 := by
      rw [Function.update_of_ne e2]
      rfl
    rw [h00, h11, h1m, h22]
    unfold saddleL
    norm_num
  rw [hbase, hupd]
  apply Real.exp_lt_exp.mpr
  nlinarith

end

end OPH.Variational

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.Variational.saddleHistory_stationary
#print axioms OPH.Variational.saddleHistory_not_minimal
#print axioms OPH.Variational.gibbs_prefers_nonstationary
