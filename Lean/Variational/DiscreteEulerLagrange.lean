import Mathlib

namespace OPH.Variational

/-!
# The discrete Euler-Lagrange equation for local path actions

A local path action sums a two-point Lagrangian over consecutive
record pairs. This module proves the discrete Euler-Lagrange equation:
a history that minimizes the local action among all real single-site
interior variations satisfies, at each supplied interior junction, the
vanishing of the sum of the two adjacent partial derivatives.

The interior site is presented as the shared endpoint of two
consecutive segments `k` and `m`, through the index identity
`k.succ = m.castSucc`; such a site is automatically neither the first
nor the last record. The Lagrangian, the path, and the derivative data are
named inputs.  The theorem does not compose with the
finite-state history-projection package, whose path space cannot contain all
real single-site variations; that interface obstruction is recorded in
`Variational.FiniteHistoryBridge`.  No physical action, saddle-point
principle, time, or continuum statement enters.
-/

variable {N : ℕ}

/-- The local path action of a two-point Lagrangian: the sum of
`L` over consecutive record pairs. -/
noncomputable def localAction (L : ℝ → ℝ → ℝ)
    (γ : Fin (N + 1) → ℝ) : ℝ :=
  ∑ n : Fin N, L (γ n.castSucc) (γ n.succ)

/-- Consecutive segments are distinct: the shared-endpoint identity
`k.succ = m.castSucc` forces `k ≠ m`. -/
theorem junction_ne {k m : Fin N} (hkm : k.succ = m.castSucc) :
    k ≠ m := by
  intro hk
  subst hk
  exact absurd hkm.symm (ne_of_lt Fin.castSucc_lt_succ)

/-- Varying one interior record changes exactly the two adjacent
segment terms: the exact difference identity for a single-site
variation at the junction of segments `k` and `m`. -/
theorem localAction_update_diff (L : ℝ → ℝ → ℝ)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N}
    (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction L (Function.update γ k.succ x) - localAction L γ
      = (L (γ k.castSucc) x - L (γ k.castSucc) (γ k.succ))
        + (L x (γ m.succ) - L (γ m.castSucc) (γ m.succ)) := by
  classical
  have hkm' := junction_ne hkm
  unfold localAction
  rw [← Finset.sum_sub_distrib]
  have hvanish : ∀ j ∈ (Finset.univ : Finset (Fin N)),
      j ∉ ({k, m} : Finset (Fin N)) →
      (L ((Function.update γ k.succ x) j.castSucc)
          ((Function.update γ k.succ x) j.succ)
        - L (γ j.castSucc) (γ j.succ)) = 0 := by
    intro j _ hj
    rw [Finset.mem_insert, Finset.mem_singleton, not_or] at hj
    have hcs : j.castSucc ≠ k.succ := by
      rw [hkm]
      exact fun h => hj.2 (Fin.castSucc_injective _ h)
    have hsc : j.succ ≠ k.succ :=
      fun h => hj.1 (Fin.succ_injective _ h)
    rw [Function.update_of_ne hcs, Function.update_of_ne hsc,
      sub_self]
  rw [← Finset.sum_subset (Finset.subset_univ ({k, m} : Finset (Fin N)))
    (fun j hj hjn => hvanish j hj hjn)]
  rw [Finset.sum_pair hkm']
  have h1 : (Function.update γ k.succ x) k.castSucc = γ k.castSucc :=
    Function.update_of_ne (ne_of_lt Fin.castSucc_lt_succ) _ _
  have h2 : (Function.update γ k.succ x) k.succ = x :=
    Function.update_self _ _ _
  have h3 : (Function.update γ k.succ x) m.castSucc = x := by
    rw [← hkm]
    exact Function.update_self _ _ _
  have h4 : (Function.update γ k.succ x) m.succ = γ m.succ :=
    Function.update_of_ne
      (by rw [hkm]; exact (ne_of_lt Fin.castSucc_lt_succ).symm) _ _
  rw [h1, h2, h3, h4]

/-- **The discrete Euler-Lagrange equation.** If the history minimizes
the local action among single-site variations at the interior junction
of segments `k` and `m`, and the two adjacent partial derivatives
exist there, then their sum vanishes:
`∂₂ L(γ(k), γ(k+1)) + ∂₁ L(γ(k+1), γ(k+2)) = 0`. -/
theorem stationary_localAction_discreteEulerLagrange
    (L : ℝ → ℝ → ℝ) (γ : Fin (N + 1) → ℝ) {k m : Fin N}
    (hkm : k.succ = m.castSucc) (d₂ d₁ : ℝ)
    (h₂ : HasDerivAt (fun x => L (γ k.castSucc) x) d₂ (γ k.succ))
    (h₁ : HasDerivAt (fun x => L x (γ m.succ)) d₁ (γ m.castSucc))
    (hmin : ∀ x, localAction L γ
      ≤ localAction L (Function.update γ k.succ x)) :
    d₂ + d₁ = 0 := by
  have hpt : γ m.castSucc = γ k.succ := by rw [hkm]
  have hφmin : IsLocalMin
      (fun x => L (γ k.castSucc) x + L x (γ m.succ)) (γ k.succ) := by
    apply Filter.Eventually.of_forall
    intro x
    have hdiff := localAction_update_diff L γ hkm x
    have := hmin x
    have hself := localAction_update_diff L γ hkm (γ k.succ)
    rw [Function.update_eq_self] at hself
    have hpt' : γ m.castSucc = γ k.succ := hpt
    nlinarith [hdiff, hself]
  have h₁' : HasDerivAt (fun x => L x (γ m.succ)) d₁ (γ k.succ) := by
    rw [← hpt]
    exact h₁
  exact hφmin.hasDerivAt_eq_zero (h₂.add h₁')

end OPH.Variational

#print axioms OPH.Variational.localAction_update_diff
#print axioms OPH.Variational.stationary_localAction_discreteEulerLagrange
