import Variational.DiscreteEulerLagrange

/-!
# A stationary maximum outside the positive-Gibbs mode set

Issue B7 (#683).  The committed history projection concentrates on
least-action histories, and the discrete Euler-Lagrange theorem derives
stationarity from minimality.  This module certifies a scoped negative for
the converse *mode/minimizer* claim: a history can satisfy the discrete
Euler-Lagrange stationarity condition at its interior junction, with
explicit derivative witnesses, while failing minimality outright, and
the Gibbs weight at every positive multiplier strictly prefers a
non-stationary variation to it.  Thus positive Gibbs weighting does not make
every stationary history a mode or minimizer.

The witness is the three-record zero history under the concave
two-point Lagrangian `L a b = -((b - a)^2) / 2`.  The single-site
variation functional at the interior record is `-x^2`, so the zero
history is a strict local maximum of the variation: stationary with
vanishing adjacent partials, beaten by every nonzero variation, and
carrying strictly smaller Gibbs weight than the unit variation at every
positive multiplier.

**Boundary.**  The witness is a stationary maximum, not a saddle.  The theorem
compares two positive unnormalized Gibbs weights; it does not say that the
stationary history has zero normalized probability or lies outside the
support.  It is also not the committed source action.  Consequently this
counterexample rules out only the universal identification of stationary
histories with positive-Gibbs modes/minimizers.  Genuine saddle recovery,
complex or signed stationary-phase weights, constrained ensembles,
source-selected real enrichments, and refinement limits remain viable.  No
physical action, time, amplitude, or continuum statement enters.
-/

namespace OPH.Variational

noncomputable section

/-- The concave two-point Lagrangian. -/
def concaveControlL : ℝ → ℝ → ℝ := fun a b => -((b - a) ^ 2) / 2

/-- The three-record zero history. -/
def stationaryMaximumHistory : Fin 3 → ℝ := fun _ => 0

/-- The zero history is stationary: the two adjacent partial
derivatives at the interior junction exist and sum to zero, in exactly
the form the discrete Euler-Lagrange theorem produces. -/
theorem stationaryMaximumHistory_stationary :
    ∃ d₂ d₁ : ℝ,
      HasDerivAt (fun x => concaveControlL (stationaryMaximumHistory (0 : Fin 2).castSucc) x)
        d₂ (stationaryMaximumHistory (0 : Fin 2).succ) ∧
      HasDerivAt (fun x => concaveControlL x (stationaryMaximumHistory (1 : Fin 2).succ))
        d₁ (stationaryMaximumHistory (1 : Fin 2).castSucc) ∧
      d₂ + d₁ = 0 := by
  refine ⟨0, 0, ?_, ?_, by norm_num⟩
  · have h : (fun x => concaveControlL (stationaryMaximumHistory (0 : Fin 2).castSucc) x) =
        fun x => -(x ^ 2) / 2 := by
      funext x
      simp [concaveControlL, stationaryMaximumHistory]
    rw [h]
    have hd : HasDerivAt (fun x : ℝ => -(x ^ 2) / 2) (-(2 * 0) / 2) 0 := by
      have := ((hasDerivAt_pow 2 (0 : ℝ)).neg).div_const 2
      simpa using this
    simpa [stationaryMaximumHistory] using hd
  · have h : (fun x => concaveControlL x (stationaryMaximumHistory (1 : Fin 2).succ)) =
        fun x => -(x ^ 2) / 2 := by
      funext x
      simp [concaveControlL, stationaryMaximumHistory]
    rw [h]
    have hd : HasDerivAt (fun x : ℝ => -(x ^ 2) / 2) (-(2 * 0) / 2) 0 := by
      have := ((hasDerivAt_pow 2 (0 : ℝ)).neg).div_const 2
      simpa using this
    simpa [stationaryMaximumHistory] using hd

/-- The zero history is not a single-site minimizer: the unit variation
strictly lowers the local action. -/
theorem stationaryMaximumHistory_not_minimal :
    ¬ ∀ x, localAction concaveControlL stationaryMaximumHistory ≤
      localAction concaveControlL
        (Function.update stationaryMaximumHistory ((0 : Fin 2).succ) x) := by
  intro h
  have h1 := h 1
  have hbase : localAction concaveControlL stationaryMaximumHistory = 0 := by
    unfold localAction concaveControlL stationaryMaximumHistory
    simp
  have hupd : localAction concaveControlL
      (Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1) = -1 := by
    unfold localAction
    rw [Fin.sum_univ_two]
    have e0 : ((0 : Fin 2).castSucc : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have e2 : ((1 : Fin 2).succ : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have h00 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).castSucc) = 0 := by
      rw [Function.update_of_ne e0]
      rfl
    have h11 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).succ) = 1 := Function.update_self _ _ _
    have h1c : ((1 : Fin 2).castSucc : Fin 3) = (0 : Fin 2).succ := by decide
    have h1m : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).castSucc) = 1 := by
      rw [h1c]
      exact Function.update_self _ _ _
    have h22 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).succ) = 0 := by
      rw [Function.update_of_ne e2]
      rfl
    rw [h00, h11, h1m, h22]
    unfold concaveControlL
    norm_num
  rw [hbase, hupd] at h1
  norm_num at h1

/-- At every positive multiplier the unnormalized Gibbs weight of the
non-stationary unit variation is strictly larger than that of the stationary
maximum.  This is an ordering statement, not a zero-mass or support claim. -/
theorem gibbs_prefers_nonstationary {lam : ℝ} (hlam : 0 < lam) :
    Real.exp (-lam * localAction concaveControlL stationaryMaximumHistory) <
      Real.exp (-lam * localAction concaveControlL
        (Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1)) := by
  have hbase : localAction concaveControlL stationaryMaximumHistory = 0 := by
    unfold localAction concaveControlL stationaryMaximumHistory
    simp
  have hupd : localAction concaveControlL
      (Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1) = -1 := by
    unfold localAction
    rw [Fin.sum_univ_two]
    have e0 : ((0 : Fin 2).castSucc : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have e2 : ((1 : Fin 2).succ : Fin 3) ≠ (0 : Fin 2).succ := by decide
    have h00 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).castSucc) = 0 := by
      rw [Function.update_of_ne e0]
      rfl
    have h11 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((0 : Fin 2).succ) = 1 := Function.update_self _ _ _
    have h1c : ((1 : Fin 2).castSucc : Fin 3) = (0 : Fin 2).succ := by decide
    have h1m : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).castSucc) = 1 := by
      rw [h1c]
      exact Function.update_self _ _ _
    have h22 : Function.update stationaryMaximumHistory ((0 : Fin 2).succ) 1
        ((1 : Fin 2).succ) = 0 := by
      rw [Function.update_of_ne e2]
      rfl
    rw [h00, h11, h1m, h22]
    unfold concaveControlL
    norm_num
  rw [hbase, hupd]
  apply Real.exp_lt_exp.mpr
  nlinarith

end

end OPH.Variational

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.Variational.stationaryMaximumHistory_stationary
#print axioms OPH.Variational.stationaryMaximumHistory_not_minimal
#print axioms OPH.Variational.gibbs_prefers_nonstationary
