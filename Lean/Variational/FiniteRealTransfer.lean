import Variational.FiniteHistoryBridge

namespace OPH.Variational

open OPH.InformationProjection

/-!
# Finite-to-real variation transfer interface

`Variational.FiniteHistoryBridge` proves that no finite family of
real-valued paths contains every real single-site variation of a fixed
path, so the finite Gibbs package and the real Euler-Lagrange package
cannot be silently identified.  This module supplies the declared
bridge as an interface: an enrichment consisting of a finite path
family, a real embedding, an action-compatibility receipt, and a
declared finite minimizer, together with an undercut receipt stating
that each real single-site variation of the embedded minimizer has
real action at or above the real action of some embedded history.

Under the receipt, the finite minimizer inequality transports to real
local minimality of the embedded minimizer along every single-site
variation at the declared site, and, given the derivative packet, to
the real discrete Euler-Lagrange equation there.  A membership form of
the receipt is impossible for every finite family; the obstruction
theorem below restates `FiniteHistoryBridge` at this interface, which
forces the receipt to be an inequality receipt rather than a closure
receipt.  A Gibbs constructor instantiates the finite side from a
projection packet whose modal history is declared, and a scalar
length-3 witness discharges every premise, with a negative control
embedding on which the receipt and the transferred conclusion both
fail.

The undercut receipt is a named input.  Nothing here derives it from
OPH source data.  The separate committed packet proves existence and a unique
positive exponential parameter at its declared empirical target; issue B7
remains open on prospective source selection of that constraint observable and
level, on the undercut-receipt derivation, and on physical attachment.
-/

variable {N : ℕ}

/-- A finite-to-real enrichment: a finite family of histories, a real
embedding, a finite action that matches the real local action on
embedded paths, and a declared finite minimizer. -/
structure FiniteRealEnrichment (Γ : Type*) [Fintype Γ] (N : ℕ) where
  /-- The two-point Lagrangian of the real local action. -/
  L : ℝ → ℝ → ℝ
  /-- The real embedding of the finite family. -/
  emb : Γ → (Fin (N + 1) → ℝ)
  /-- The finite action, for example a Gibbs effective action. -/
  S : Γ → ℝ
  /-- Action compatibility: the finite action equals the real local
  action of the embedded path. -/
  compat : ∀ g, S g = localAction L (emb g)
  /-- The declared finite minimizer. -/
  g₀ : Γ
  /-- The finite minimizer inequality. -/
  min₀ : ∀ g, S g₀ ≤ S g

/-- **The undercut receipt.**  Every real single-site variation of the
embedded minimizer at the site `k.succ` has real action at or above the
real action of some embedded history.  This is the declared
dense-direction surrogate: the finite family reaches at or below every
real variation without containing it. -/
def UndercutReceipt {Γ : Type*} [Fintype Γ]
    (E : FiniteRealEnrichment Γ N) (k : Fin N) : Prop :=
  ∀ x : ℝ, ∃ g : Γ,
    localAction E.L (E.emb g)
      ≤ localAction E.L (Function.update (E.emb E.g₀) k.succ x)

/-- The finite minimizer inequality transports along the embedding:
the embedded minimizer minimizes the real local action over the
embedded family. -/
theorem embedded_minimum {Γ : Type*} [Fintype Γ]
    (E : FiniteRealEnrichment Γ N) (g : Γ) :
    localAction E.L (E.emb E.g₀) ≤ localAction E.L (E.emb g) := by
  have h := E.min₀ g
  rwa [E.compat, E.compat] at h

/-- **Conditional transfer of minimality.**  Under the undercut
receipt, the finite minimizer inequality transports to real local
minimality of the embedded minimizer along every real single-site
variation at the site `k.succ`.  The universal real quantifier on the
right is exactly the one no finite family can satisfy by membership. -/
theorem transfer_local_min {Γ : Type*} [Fintype Γ]
    (E : FiniteRealEnrichment Γ N) (k : Fin N)
    (hrec : UndercutReceipt E k) :
    ∀ x : ℝ, localAction E.L (E.emb E.g₀)
      ≤ localAction E.L (Function.update (E.emb E.g₀) k.succ x) := by
  intro x
  obtain ⟨g, hg⟩ := hrec x
  exact le_trans (embedded_minimum E g) hg

/-- **Conditional transfer of stationarity.**  Under the undercut
receipt and the derivative packet of `L` at the embedded minimizer,
the finite Gibbs-side minimizer satisfies the real discrete
Euler-Lagrange equation at the declared interior junction:
`∂₂ L + ∂₁ L = 0`.  This composes `transfer_local_min` with
`stationary_localAction_discreteEulerLagrange`. -/
theorem transfer_discreteEulerLagrange {Γ : Type*} [Fintype Γ]
    (E : FiniteRealEnrichment Γ N) {k m : Fin N}
    (hkm : k.succ = m.castSucc) (d₂ d₁ : ℝ)
    (h₂ : HasDerivAt (fun x => E.L (E.emb E.g₀ k.castSucc) x) d₂
      (E.emb E.g₀ k.succ))
    (h₁ : HasDerivAt (fun x => E.L x (E.emb E.g₀ m.succ)) d₁
      (E.emb E.g₀ m.castSucc))
    (hrec : UndercutReceipt E k) :
    d₂ + d₁ = 0 :=
  stationary_localAction_discreteEulerLagrange E.L (E.emb E.g₀) hkm
    d₂ d₁ h₂ h₁ (transfer_local_min E k hrec)

/-- **Membership receipts are impossible.**  For every enrichment and
site, some real single-site variation of the embedded minimizer lies
outside the embedded family, so the undercut receipt cannot be
discharged by closure of the finite family under variations.  This
restates the `FiniteHistoryBridge` obstruction at the transfer
interface and forces the receipt to be an inequality receipt. -/
theorem no_membership_receipt {Γ : Type*} [Fintype Γ]
    (E : FiniteRealEnrichment Γ N) (i : Fin (N + 1)) :
    ¬ ∀ x : ℝ, Function.update (E.emb E.g₀) i x ∈ Set.range E.emb := by
  intro h
  have hsub : Set.range (fun x : ℝ => Function.update (E.emb E.g₀) i x)
      ⊆ Set.range E.emb := by
    rintro _ ⟨x, rfl⟩
    exact h x
  exact (Set.infinite_range_of_injective
      (updateAt_injective (E.emb E.g₀) i)).not_finite
    ((Set.finite_range E.emb).subset hsub)

/-- Build an enrichment from a Gibbs projection packet with a declared
modal history: the finite action is the effective action of the
projected law, and the finite minimizer inequality is
`modal_path_least_action`.  The embedding and its action-compatibility
receipt are named enrichment data. -/
noncomputable def ofGibbsModal {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
    {A : Type*} [Fintype A]
    (L : ℝ → ℝ → ℝ) (τ R : PathSpace Ω N → ℝ)
    (F : A → PathSpace Ω N → ℝ) (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (emb : PathSpace Ω N → (Fin (N + 1) → ℝ))
    (hcompat : ∀ γ, effectiveAction R F lam γ = localAction L (emb γ))
    (γ₀ : PathSpace Ω N) (hmax : ∀ δ, τ δ ≤ τ γ₀) :
    FiniteRealEnrichment (PathSpace Ω N) N where
  L := L
  emb := emb
  S := effectiveAction R F lam
  compat := hcompat
  g₀ := γ₀
  min₀ := modal_path_least_action τ R F lam logZ hτ hlog γ₀ hmax

/-! ## Scalar length-3 witness -/

/-- The free two-point Lagrangian of the witnesses. -/
def witnessL (x y : ℝ) : ℝ := (y - x) ^ 2

/-- The witness family of two scalar chains with 3 records: the affine
chain `0,1,2` and the bent chain `0,0,2`. -/
def witnessEmb : Fin 2 → (Fin 3 → ℝ) := ![![0, 1, 2], ![0, 0, 2]]

/-- Local action of a 3-record chain: the sum of the two segment
terms. -/
theorem localAction_three (L : ℝ → ℝ → ℝ) (p : Fin 3 → ℝ) :
    localAction L p = L (p 0) (p 1) + L (p 1) (p 2) := by
  rw [localAction, Fin.sum_univ_two, Fin.castSucc_zero,
    Fin.succ_zero_eq_one, Fin.castSucc_one, Fin.succ_one_eq_two]

/-- Local action of a 3-record chain after a single-site variation at
the interior record. -/
theorem localAction_update_one (L : ℝ → ℝ → ℝ) (p : Fin 3 → ℝ)
    (x : ℝ) :
    localAction L (Function.update p 1 x) = L (p 0) x + L x (p 2) := by
  rw [localAction_three,
    Function.update_of_ne (by decide : (0 : Fin 3) ≠ 1),
    Function.update_of_ne (by decide : (2 : Fin 3) ≠ 1),
    Function.update_self]

/-- The witness enrichment: the affine chain is the declared finite
minimizer of the free local action over the two-chain family. -/
noncomputable def witnessEnrichment : FiniteRealEnrichment (Fin 2) 2 where
  L := witnessL
  emb := witnessEmb
  S := fun g => localAction witnessL (witnessEmb g)
  compat := fun _ => rfl
  g₀ := 0
  min₀ := by
    intro g
    fin_cases g
    · simp [localAction_three, witnessEmb, witnessL]
    · simp [localAction_three, witnessEmb, witnessL]
      norm_num

/-- The undercut receipt holds on the witness: the affine chain itself
undercuts every real interior variation, since
`x^2 + (2 - x)^2 - 2 = 2 (x - 1)^2 ≥ 0`. -/
theorem witness_receipt : UndercutReceipt witnessEnrichment 0 := by
  intro x
  refine ⟨0, ?_⟩
  show localAction witnessL (witnessEmb 0)
    ≤ localAction witnessL
        (Function.update (witnessEmb 0) (0 : Fin 2).succ x)
  rw [Fin.succ_zero_eq_one, localAction_three, localAction_update_one]
  simp [witnessEmb, witnessL]
  nlinarith [sq_nonneg (x - 1)]

/-- **Transferred witness minimality.**  The embedded witness minimizer
is really locally minimal along every real single-site variation at the
interior record, obtained through `transfer_local_min`. -/
theorem witness_real_local_min (x : ℝ) :
    localAction witnessL (witnessEmb 0)
      ≤ localAction witnessL
          (Function.update (witnessEmb 0) (0 : Fin 2).succ x) :=
  transfer_local_min witnessEnrichment 0 witness_receipt x

/-- **Transferred witness stationarity.**  A nonzero derivative packet
of the free Lagrangian at the embedded witness minimizer discharges
every premise of `transfer_discreteEulerLagrange`, and the discrete
Euler-Lagrange sum vanishes through the transfer theorem. -/
theorem witness_transfer_euler_lagrange :
    ∃ d₂ d₁ : ℝ,
      HasDerivAt
          (fun x => witnessEnrichment.L
            (witnessEnrichment.emb witnessEnrichment.g₀
              (0 : Fin 2).castSucc) x) d₂
          (witnessEnrichment.emb witnessEnrichment.g₀ (0 : Fin 2).succ)
        ∧ HasDerivAt
            (fun x => witnessEnrichment.L x
              (witnessEnrichment.emb witnessEnrichment.g₀
                (1 : Fin 2).succ)) d₁
            (witnessEnrichment.emb witnessEnrichment.g₀
              (1 : Fin 2).castSucc)
        ∧ d₂ ≠ 0
        ∧ d₂ + d₁ = 0 := by
  have h₂ : HasDerivAt
      (fun x => witnessEnrichment.L
        (witnessEnrichment.emb witnessEnrichment.g₀
          (0 : Fin 2).castSucc) x) 2
      (witnessEnrichment.emb witnessEnrichment.g₀ (0 : Fin 2).succ) := by
    have h : HasDerivAt (fun x : ℝ => (x - 0) ^ 2) 2 1 := by
      convert ((hasDerivAt_id (1 : ℝ)).sub_const 0).pow 2 using 1
      norm_num
    show HasDerivAt (fun x => witnessL (witnessEmb 0 0) x) 2
      (witnessEmb 0 1)
    have h0 : witnessEmb 0 0 = (0 : ℝ) := by simp [witnessEmb]
    have h1 : witnessEmb 0 1 = (1 : ℝ) := by simp [witnessEmb]
    rw [h0, h1]
    exact h
  have h₁ : HasDerivAt
      (fun x => witnessEnrichment.L x
        (witnessEnrichment.emb witnessEnrichment.g₀
          (1 : Fin 2).succ)) (-2)
      (witnessEnrichment.emb witnessEnrichment.g₀
        (1 : Fin 2).castSucc) := by
    have h : HasDerivAt (fun x : ℝ => (2 - x) ^ 2) (-2) 1 := by
      convert ((hasDerivAt_id (1 : ℝ)).const_sub 2).pow 2 using 1
      norm_num
    show HasDerivAt (fun x => witnessL x (witnessEmb 0 2)) (-2)
      (witnessEmb 0 1)
    have h2 : witnessEmb 0 2 = (2 : ℝ) := by simp [witnessEmb]
    have h1 : witnessEmb 0 1 = (1 : ℝ) := by simp [witnessEmb]
    rw [h2, h1]
    exact h
  exact ⟨2, -2, h₂, h₁, by norm_num,
    transfer_discreteEulerLagrange witnessEnrichment
      (k := 0) (m := 1) (by decide) 2 (-2) h₂ h₁ witness_receipt⟩

/-! ## Negative control -/

/-- The control family: the bent chain `0,0,2` alone, vacuously its
own finite minimizer. -/
def controlEmb : Fin 1 → (Fin 3 → ℝ) := ![![0, 0, 2]]

/-- The control enrichment: action-compatible and finitely minimal,
yet its embedded minimizer is the wrong real path. -/
noncomputable def controlEnrichment : FiniteRealEnrichment (Fin 1) 2 where
  L := witnessL
  emb := controlEmb
  S := fun g => localAction witnessL (controlEmb g)
  compat := fun _ => rfl
  g₀ := 0
  min₀ := by
    intro g
    rw [Subsingleton.elim g 0]

/-- **Control: the transferred conclusion fails.**  The variation
`x = 1` drops the real action of the embedded control minimizer from
`4` to `2`, so real local minimality along embedded directions is
false for the control embedding. -/
theorem control_local_min_fails :
    ¬ ∀ x : ℝ, localAction witnessL (controlEmb 0)
      ≤ localAction witnessL
          (Function.update (controlEmb 0) (0 : Fin 2).succ x) := by
  intro h
  have h1 := h 1
  rw [Fin.succ_zero_eq_one, localAction_three,
    localAction_update_one] at h1
  simp [controlEmb, witnessL] at h1
  norm_num at h1

/-- **Control: the undercut receipt fails.**  The failure is forced by
`transfer_local_min` through the failed conclusion: were the receipt
available for the control embedding, the transferred minimality would
follow, and it is refuted at `x = 1`.  The receipt is therefore a
substantive obligation of the transfer interface, satisfiable for one
embedding of the same finite data and unsatisfiable for another. -/
theorem control_receipt_fails :
    ¬ UndercutReceipt controlEnrichment 0 := fun hrec =>
  control_local_min_fails (transfer_local_min controlEnrichment 0 hrec)

end OPH.Variational

#print axioms OPH.Variational.embedded_minimum
#print axioms OPH.Variational.transfer_local_min
#print axioms OPH.Variational.transfer_discreteEulerLagrange
#print axioms OPH.Variational.no_membership_receipt
#print axioms OPH.Variational.witness_receipt
#print axioms OPH.Variational.witness_real_local_min
#print axioms OPH.Variational.witness_transfer_euler_lagrange
#print axioms OPH.Variational.control_local_min_fails
#print axioms OPH.Variational.control_receipt_fails
