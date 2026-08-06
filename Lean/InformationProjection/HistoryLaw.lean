import InformationProjection.PathGibbs
import InformationProjection.GlobalObjective

namespace OPH.InformationProjection

open OPH.Thermodynamics

/-!
# The history law as an Axiom-3 information projection

`InformationProjection.PathGibbs` consumes a supplied exponential path
law and proves projection identities and least-action ordering for it.
This module derives the exponential form itself at the representation
level: over any finite history space, for a strictly positive reference
law and an action functional, the exponentially tilted law is the unique
minimizer of relative entropy to the reference plus the
multiplier-weighted mean action, equivalently the unique maximizer of
entropy relative to the reference minus the weighted mean action.  On
the surface of laws sharing the tilted mean action, the tilt is the
unique minimizer of relative entropy to the reference, with the exact
finite Pythagorean identity.  The displayed packet of `PathGibbs`
instantiates this projection: its projected law equals the tilt of its
reference by its multiplier aggregate, its normalization constant is
forced, and its history law is therefore the information projection of
the reference onto the mean-action constraint, with uniqueness.  A
three-history witness realizes the projection with exact expressions in
`exp (-1)`, and a control law with the same mean action and vanishing
mass on one history scores strictly larger relative entropy.

Scope boundary: representation-level derivation.  The reference path
weight, the constraint observables, and the multipliers remain named
inputs; no construction here produces them from OPH source data, and
that production is the open core of issue B7.  The finite-to-real
variation bridge receipts of `Variational.FiniteRealTransfer` are
unchanged by this module.
-/

section GibbsVariational

variable {Γ : Type*} [Fintype Γ] [DecidableEq Γ]

/-- Normalizer of the exponential tilt: the reference-weighted sum of
`exp (-(lam * S g))` over the finite history space. -/
noncomputable def tiltZ (τ S : Γ → ℝ) (lam : ℝ) : ℝ :=
  ∑ g, τ g * Real.exp (-(lam * S g))

/-- The exponential tilt of the reference law `τ` by the action `S` at
multiplier `lam`: `τ g * exp (-(lam * S g))`, normalized. -/
noncomputable def tilt (τ S : Γ → ℝ) (lam : ℝ) (g : Γ) : ℝ :=
  τ g * Real.exp (-(lam * S g)) / tiltZ τ S lam

/-- The mean action of a law: `∑ g, ρ g * S g`. -/
noncomputable def meanAction (ρ S : Γ → ℝ) : ℝ :=
  ∑ g, ρ g * S g

/-- The variational objective of the finite Gibbs principle: relative
entropy to the reference plus the multiplier-weighted mean action. -/
noncomputable def freeEnergy (τ S : Γ → ℝ) (lam : ℝ) (ρ : Γ → ℝ) : ℝ :=
  kl ρ τ + lam * meanAction ρ S

omit [DecidableEq Γ] in
theorem tiltZ_pos [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) : 0 < tiltZ τ S lam :=
  Finset.sum_pos (fun g _ => mul_pos (hτ g) (Real.exp_pos _))
    Finset.univ_nonempty

omit [DecidableEq Γ] in
theorem tilt_pos [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) (g : Γ) : 0 < tilt τ S lam g :=
  div_pos (mul_pos (hτ g) (Real.exp_pos _)) (tiltZ_pos τ S lam hτ)

omit [DecidableEq Γ] in
/-- The tilt is a probability law on the finite history space. -/
theorem tilt_total [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) : ∑ g, tilt τ S lam g = 1 := by
  unfold tilt
  rw [← Finset.sum_div,
    show (∑ g, τ g * Real.exp (-(lam * S g))) = tiltZ τ S lam from rfl]
  exact div_self (ne_of_gt (tiltZ_pos τ S lam hτ))

omit [DecidableEq Γ] in
/-- Exact log form of the tilt: the displayed exponential relation of
the `PathGibbs` packet, with `log (tiltZ)` in the constant slot. -/
theorem tilt_log [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) (g : Γ) :
    Real.log (tilt τ S lam g)
      = Real.log (τ g) - lam * S g - Real.log (tiltZ τ S lam) := by
  unfold tilt
  rw [Real.log_div (ne_of_gt (mul_pos (hτ g) (Real.exp_pos _)))
      (ne_of_gt (tiltZ_pos τ S lam hτ)),
    Real.log_mul (ne_of_gt (hτ g)) (ne_of_gt (Real.exp_pos _)),
    Real.log_exp]
  ring

omit [DecidableEq Γ] in
/-- Pointwise split of one relative-entropy term against the tilt into
the term against the reference plus the linear action contribution. -/
theorem klTerm_tilt_split [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) (r : ℝ) (g : Γ) :
    klTerm r (tilt τ S lam g)
      = klTerm r (τ g) + r * (lam * S g + Real.log (tiltZ τ S lam)) := by
  by_cases hr : r = 0
  · rw [hr, klTerm_zero, klTerm_zero, zero_mul, add_zero]
  rw [klTerm_of_pos hr, klTerm_of_pos hr,
    Real.log_div hr (ne_of_gt (tilt_pos τ S lam hτ g)),
    Real.log_div hr (ne_of_gt (hτ g)), tilt_log τ S lam hτ g]
  ring

omit [DecidableEq Γ] in
/-- **Exact free-energy identity.**  For every normalized law `ρ`,
relative entropy to the tilt equals the variational objective plus the
log normalizer: `D(ρ ‖ w_lam) = D(ρ ‖ τ) + lam * E_ρ[S] + log Z`. -/
theorem kl_tilt_split [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ) (ρ : Γ → ℝ)
    (hτ : ∀ g, 0 < τ g) (hρ1 : ∑ g, ρ g = 1) :
    kl ρ (tilt τ S lam)
      = freeEnergy τ S lam ρ + Real.log (tiltZ τ S lam) := by
  unfold freeEnergy meanAction kl
  rw [Finset.sum_congr rfl fun g _ => klTerm_tilt_split τ S lam hτ (ρ g) g,
    Finset.sum_add_distrib,
    Finset.sum_congr rfl
      (fun g _ => by ring :
        ∀ g ∈ Finset.univ,
          ρ g * (lam * S g + Real.log (tiltZ τ S lam))
            = lam * (ρ g * S g) + Real.log (tiltZ τ S lam) * ρ g),
    Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, hρ1,
    mul_one]
  ring

/-- The tilt scores `-log Z` on the variational objective. -/
theorem freeEnergy_tilt [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (hτ : ∀ g, 0 < τ g) :
    freeEnergy τ S lam (tilt τ S lam) = -Real.log (tiltZ τ S lam) := by
  have h := kl_tilt_split τ S lam (tilt τ S lam) hτ
    (tilt_total τ S lam hτ)
  rw [kl_self] at h
  linarith

/-- **Finite Gibbs variational principle, minimizer half.**  Over
normalized nonnegative laws on the finite history space, the
exponential tilt minimizes relative entropy to the reference plus the
multiplier-weighted mean action. -/
theorem tilt_minimizes_freeEnergy [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g) (hρ1 : ∑ g, ρ g = 1)
    (hτ : ∀ g, 0 < τ g) :
    freeEnergy τ S lam (tilt τ S lam) ≤ freeEnergy τ S lam ρ := by
  have hid := kl_tilt_split τ S lam ρ hτ hρ1
  have hkl : 0 ≤ kl ρ (tilt τ S lam) :=
    kl_nonneg ρ (tilt τ S lam) hρ0
      (fun g => le_of_lt (tilt_pos τ S lam hτ g))
      (fun g hg => absurd hg (ne_of_gt (tilt_pos τ S lam hτ g)))
      (by rw [hρ1, tilt_total τ S lam hτ])
  rw [freeEnergy_tilt τ S lam hτ]
  linarith

/-- **Finite Gibbs variational principle, strict half.**  Every
normalized nonnegative law distinct from the tilt scores strictly
larger on the variational objective; strict positivity of relative
entropy away from equality carries the gap. -/
theorem tilt_strict_freeEnergy [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g) (hρ1 : ∑ g, ρ g = 1)
    (hτ : ∀ g, 0 < τ g) (hne : ρ ≠ tilt τ S lam) :
    freeEnergy τ S lam (tilt τ S lam) < freeEnergy τ S lam ρ := by
  have hid := kl_tilt_split τ S lam ρ hτ hρ1
  have hkl : 0 < kl ρ (tilt τ S lam) :=
    kl_pos_of_ne ρ (tilt τ S lam) hρ0
      (fun g => le_of_lt (tilt_pos τ S lam hτ g))
      (fun g hg => absurd hg (ne_of_gt (tilt_pos τ S lam hτ g)))
      (by rw [hρ1, tilt_total τ S lam hτ]) hne
  rw [freeEnergy_tilt τ S lam hτ]
  linarith

/-- **Uniqueness.**  A normalized nonnegative law attaining the
variational minimum is the tilt. -/
theorem tilt_unique_minimizer [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g) (hρ1 : ∑ g, ρ g = 1)
    (hτ : ∀ g, 0 < τ g)
    (h : freeEnergy τ S lam ρ = freeEnergy τ S lam (tilt τ S lam)) :
    ρ = tilt τ S lam := by
  by_contra hne
  exact absurd h
    (ne_of_gt (tilt_strict_freeEnergy τ S lam ρ hρ0 hρ1 hτ hne))

/-- Entropy reading of the same principle: the tilt maximizes entropy
relative to the reference minus the weighted mean action, since that
score is the negation of the variational objective. -/
theorem tilt_maximizes_entropy_score [Nonempty Γ] (τ S : Γ → ℝ)
    (lam : ℝ) (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g)
    (hρ1 : ∑ g, ρ g = 1) (hτ : ∀ g, 0 < τ g) :
    -kl ρ τ - lam * meanAction ρ S
      ≤ -kl (tilt τ S lam) τ - lam * meanAction (tilt τ S lam) S := by
  have h := tilt_minimizes_freeEnergy τ S lam ρ hρ0 hρ1 hτ
  unfold freeEnergy at h
  linarith

/-- **Constraint form, minimizer half.**  Among normalized nonnegative
laws whose mean action equals the tilted mean action, the tilt
minimizes relative entropy to the reference: the Lagrangian principle
restricted to the constraint surface. -/
theorem tilt_constrained_minimizer [Nonempty Γ] (τ S : Γ → ℝ)
    (lam : ℝ) (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g)
    (hρ1 : ∑ g, ρ g = 1) (hτ : ∀ g, 0 < τ g)
    (hmom : meanAction ρ S = meanAction (tilt τ S lam) S) :
    kl (tilt τ S lam) τ ≤ kl ρ τ := by
  have h := tilt_minimizes_freeEnergy τ S lam ρ hρ0 hρ1 hτ
  unfold freeEnergy at h
  rw [hmom] at h
  linarith

/-- **Constraint form, strict half.**  A law on the constraint surface
distinct from the tilt scores strictly larger relative entropy to the
reference. -/
theorem tilt_constrained_strict [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g) (hρ1 : ∑ g, ρ g = 1)
    (hτ : ∀ g, 0 < τ g)
    (hmom : meanAction ρ S = meanAction (tilt τ S lam) S)
    (hne : ρ ≠ tilt τ S lam) :
    kl (tilt τ S lam) τ < kl ρ τ := by
  have h := tilt_strict_freeEnergy τ S lam ρ hρ0 hρ1 hτ hne
  unfold freeEnergy at h
  rw [hmom] at h
  linarith

/-- **Constraint form, uniqueness.**  A law on the constraint surface
attaining the minimal relative entropy to the reference is the tilt. -/
theorem tilt_constrained_unique [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hρ0 : ∀ g, 0 ≤ ρ g) (hρ1 : ∑ g, ρ g = 1)
    (hτ : ∀ g, 0 < τ g)
    (hmom : meanAction ρ S = meanAction (tilt τ S lam) S)
    (h : kl ρ τ = kl (tilt τ S lam) τ) : ρ = tilt τ S lam := by
  apply tilt_unique_minimizer τ S lam ρ hρ0 hρ1 hτ
  unfold freeEnergy
  rw [hmom, h]

/-- **Pythagorean identity at the tilt.**  On the constraint surface,
relative entropy to the reference splits exactly through the tilt:
`D(ρ ‖ τ) = D(ρ ‖ w_lam) + D(w_lam ‖ τ)`. -/
theorem tilt_pythagorean [Nonempty Γ] (τ S : Γ → ℝ) (lam : ℝ)
    (ρ : Γ → ℝ) (hτ : ∀ g, 0 < τ g) (hρ1 : ∑ g, ρ g = 1)
    (hmom : meanAction ρ S = meanAction (tilt τ S lam) S) :
    kl ρ τ = kl ρ (tilt τ S lam) + kl (tilt τ S lam) τ := by
  have h1 := kl_tilt_split τ S lam ρ hτ hρ1
  have h2 := kl_tilt_split τ S lam (tilt τ S lam) hτ
    (tilt_total τ S lam hτ)
  rw [kl_self] at h2
  unfold freeEnergy at h1 h2
  rw [hmom] at h1
  linarith

end GibbsVariational

/-- A finite type carrying a real family of unit total mass is
inhabited. -/
theorem nonempty_of_sum_eq_one {Γ : Type*} [Fintype Γ] {f : Γ → ℝ}
    (h : ∑ g, f g = 1) : Nonempty Γ := by
  rcases isEmpty_or_nonempty Γ with hemp | hne
  · haveI := hemp
    rw [Finset.univ_eq_empty, Finset.sum_empty] at h
    norm_num at h
  · exact hne

section HistoryPacket

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω] {N : ℕ}
variable {A : Type*} [Fintype A]

omit [Fintype Ω] [DecidableEq Ω] in
/-- Pointwise consequence of the displayed exponential log relation:
the projected path law is the reference times the exponential of minus
the multiplier aggregate, times `exp (-logZ)`. -/
theorem pathGibbs_pointwise
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (γ : PathSpace Ω N) :
    τ γ = R γ * Real.exp (-(∑ a, lam a * F a γ)) * Real.exp (-logZ) := by
  rw [← Real.exp_log (hτ γ), hlog γ,
    show Real.log (R γ) - (∑ a, lam a * F a γ) - logZ
      = Real.log (R γ) + -(∑ a, lam a * F a γ) + -logZ from by ring,
    Real.exp_add, Real.exp_add, Real.exp_log (hR γ)]

omit [DecidableEq Ω] in
/-- Partition identity of the packet: the tilt normalizer of the
reference by the multiplier aggregate cancels `exp (-logZ)` exactly. -/
theorem pathGibbs_partition
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ) :
    tiltZ R (fun γ => ∑ a, lam a * F a γ) 1 * Real.exp (-logZ) = 1 := by
  have hsum : ∑ γ, R γ * Real.exp (-((1 : ℝ) * (∑ a, lam a * F a γ)))
        * Real.exp (-logZ) = ∑ γ, τ γ :=
    Finset.sum_congr rfl fun γ _ => by
      rw [one_mul,
        ← pathGibbs_pointwise τ R F lam logZ hτ hR hlog γ]
  unfold tiltZ
  rw [Finset.sum_mul, hsum, hτ1]

omit [DecidableEq Ω] in
/-- **The history law is the exponential tilt.**  The projected path
law of the `PathGibbs` packet equals the information-projection tilt of
its reference weight by its multiplier aggregate at unit multiplier.
The exponential form is derived from the packet data, and every slot of
the tilt is consumed from the packet. -/
theorem pathGibbs_eq_tilt
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ) :
    τ = tilt R (fun γ => ∑ a, lam a * F a γ) 1 := by
  have hZ := pathGibbs_partition τ R F lam logZ hτ hR hτ1 hlog
  have hZne : tiltZ R (fun γ => ∑ a, lam a * F a γ) 1 ≠ 0 := by
    intro h0
    rw [h0, zero_mul] at hZ
    norm_num at hZ
  funext γ
  unfold tilt
  rw [one_mul, eq_div_iff hZne,
    pathGibbs_pointwise τ R F lam logZ hτ hR hlog γ,
    show R γ * Real.exp (-(∑ a, lam a * F a γ)) * Real.exp (-logZ)
        * tiltZ R (fun γ => ∑ a, lam a * F a γ) 1
      = R γ * Real.exp (-(∑ a, lam a * F a γ))
        * (tiltZ R (fun γ => ∑ a, lam a * F a γ) 1
          * Real.exp (-logZ)) from by ring,
    hZ, mul_one]

omit [DecidableEq Ω] in
/-- **The normalization constant is forced.**  The packet's `logZ` slot
equals the log of the tilt normalizer built from the reference and the
multiplier aggregate, so no freedom remains in the constant. -/
theorem pathGibbs_logZ_forced
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ) :
    logZ = Real.log (tiltZ R (fun γ => ∑ a, lam a * F a γ) 1) := by
  have hZ := pathGibbs_partition τ R F lam logZ hτ hR hτ1 hlog
  have hexp : tiltZ R (fun γ => ∑ a, lam a * F a γ) 1
      = Real.exp logZ := by
    have h := congrArg (fun t => t * Real.exp logZ) hZ
    dsimp only [] at h
    rwa [mul_assoc, ← Real.exp_add, neg_add_cancel, Real.exp_zero,
      mul_one, one_mul] at h
  rw [hexp, Real.log_exp]

/-- **A3 constrained minimizer at the packet.**  Among normalized
nonnegative path laws whose mean multiplier aggregate matches that of
the history law, the history law minimizes relative entropy to the
reference weight. -/
theorem historyLaw_constrained_minimizer
    (P τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hP0 : ∀ γ, 0 ≤ P γ) (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hP1 : ∑ γ, P γ = 1) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (hmom : meanAction P (fun γ => ∑ a, lam a * F a γ)
      = meanAction τ (fun γ => ∑ a, lam a * F a γ)) :
    kl τ R ≤ kl P R := by
  haveI : Nonempty (PathSpace Ω N) := nonempty_of_sum_eq_one hτ1
  have heq := pathGibbs_eq_tilt τ R F lam logZ hτ hR hτ1 hlog
  rw [heq] at hmom ⊢
  exact tilt_constrained_minimizer R (fun γ => ∑ a, lam a * F a γ) 1 P
    hP0 hP1 hR hmom

/-- **A3 constrained uniqueness at the packet.**  A normalized
nonnegative path law on the same mean-aggregate surface attaining the
minimal relative entropy to the reference is the history law itself. -/
theorem historyLaw_constrained_unique
    (P τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hP0 : ∀ γ, 0 ≤ P γ) (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hP1 : ∑ γ, P γ = 1) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (hmom : meanAction P (fun γ => ∑ a, lam a * F a γ)
      = meanAction τ (fun γ => ∑ a, lam a * F a γ))
    (h : kl P R = kl τ R) : P = τ := by
  haveI : Nonempty (PathSpace Ω N) := nonempty_of_sum_eq_one hτ1
  have heq := pathGibbs_eq_tilt τ R F lam logZ hτ hR hτ1 hlog
  rw [heq] at hmom h ⊢
  exact tilt_constrained_unique R (fun γ => ∑ a, lam a * F a γ) 1 P
    hP0 hP1 hR hmom h

/-- **Uniqueness at the committed moment surface.**  Under the full
moment-matching premise of `pathGibbs_minimizer`, a path law attaining
the minimal relative entropy to the reference equals the projected law.
This closes the uniqueness half left open by the committed minimizer
theorem, consuming `pathGibbs_pythagorean` directly. -/
theorem pathGibbs_minimizer_unique
    (P τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hP0 : ∀ γ, 0 ≤ P γ) (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hP1 : ∑ γ, P γ = 1) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (hmom : ∀ a, ∑ γ, P γ * F a γ = ∑ γ, τ γ * F a γ)
    (h : kl P R = kl τ R) : P = τ := by
  have hpyth := pathGibbs_pythagorean P τ R F lam logZ hP0 hτ hR hP1
    hτ1 hlog hmom
  have hzero : kl P τ = 0 := by linarith
  exact eq_of_kl_eq_zero P τ hP0 (fun γ => le_of_lt (hτ γ))
    (fun γ hγ => absurd hγ (ne_of_gt (hτ γ))) (by rw [hP1, hτ1]) hzero

/-- **A3 reading.**  The history law of the `PathGibbs` packet is the
information projection of its reference weight onto the mean-action
constraint surface: it is the exponential tilt of the reference by the
multiplier aggregate, it minimizes relative entropy to the reference
over the surface, and it is the unique law doing so. -/
theorem historyLaw_information_projection
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ) :
    τ = tilt R (fun γ => ∑ a, lam a * F a γ) 1
      ∧ (∀ P : PathSpace Ω N → ℝ, (∀ γ, 0 ≤ P γ) → ∑ γ, P γ = 1 →
          meanAction P (fun γ => ∑ a, lam a * F a γ)
            = meanAction τ (fun γ => ∑ a, lam a * F a γ) →
          kl τ R ≤ kl P R)
      ∧ (∀ P : PathSpace Ω N → ℝ, (∀ γ, 0 ≤ P γ) → ∑ γ, P γ = 1 →
          meanAction P (fun γ => ∑ a, lam a * F a γ)
            = meanAction τ (fun γ => ∑ a, lam a * F a γ) →
          kl P R = kl τ R → P = τ) :=
  ⟨pathGibbs_eq_tilt τ R F lam logZ hτ hR hτ1 hlog,
    fun P hP0 hP1 hmom => historyLaw_constrained_minimizer P τ R F lam
      logZ hP0 hτ hR hP1 hτ1 hlog hmom,
    fun P hP0 hP1 hmom h => historyLaw_constrained_unique P τ R F lam
      logZ hP0 hτ hR hP1 hτ1 hlog hmom h⟩

end HistoryPacket

section Witness

/-- Witness reference law on three histories: uniform masses `1/3`. -/
noncomputable def historyRef : Fin 3 → ℝ := fun _ => 1 / 3

theorem historyRef_apply : ∀ g : Fin 3, historyRef g = 1 / 3 :=
  fun _ => rfl

theorem historyRef_pos : ∀ g, 0 < historyRef g := by
  intro g
  rw [historyRef_apply]
  norm_num

theorem historyRef_sum : ∑ g, historyRef g = 1 := by
  rw [Fin.sum_univ_three, historyRef_apply, historyRef_apply,
    historyRef_apply]
  norm_num

/-- Witness action: history `0` has action `0`; histories `1` and `2`
have action `1`. -/
noncomputable def historyAction : Fin 3 → ℝ :=
  fun g => if g = 0 then 0 else 1

theorem historyAction_apply_0 : historyAction 0 = 0 := by
  unfold historyAction
  rw [if_pos rfl]

theorem historyAction_apply_1 : historyAction 1 = 1 := by
  unfold historyAction
  rw [if_neg (by decide : ¬ (1 : Fin 3) = 0)]

theorem historyAction_apply_2 : historyAction 2 = 1 := by
  unfold historyAction
  rw [if_neg (by decide : ¬ (2 : Fin 3) = 0)]

/-- The witness history law: the tilt of the uniform reference by the
witness action at unit multiplier. -/
noncomputable def historyTilt : Fin 3 → ℝ :=
  tilt historyRef historyAction 1

theorem historyDenom_pos : 0 < 1 + 2 * Real.exp (-1) := by
  have := Real.exp_pos (-1 : ℝ)
  linarith

/-- Exact witness normalizer: `(1 + 2 exp (-1)) / 3`. -/
theorem historyTiltZ_eq :
    tiltZ historyRef historyAction 1 = (1 + 2 * Real.exp (-1)) / 3 := by
  unfold tiltZ
  rw [Fin.sum_univ_three, historyRef_apply, historyRef_apply,
    historyRef_apply, historyAction_apply_0, historyAction_apply_1,
    historyAction_apply_2,
    show -((1 : ℝ) * (0 : ℝ)) = 0 from by norm_num,
    show -((1 : ℝ) * (1 : ℝ)) = -1 from by norm_num, Real.exp_zero]
  ring

theorem historyTilt_apply_0 :
    historyTilt 0 = 1 / (1 + 2 * Real.exp (-1)) := by
  unfold historyTilt tilt
  rw [historyTiltZ_eq, historyRef_apply, historyAction_apply_0,
    show -((1 : ℝ) * (0 : ℝ)) = 0 from by norm_num, Real.exp_zero,
    mul_one]
  have hd : (1 + 2 * Real.exp (-1)) ≠ 0 := ne_of_gt historyDenom_pos
  field_simp

theorem historyTilt_apply_1 :
    historyTilt 1 = Real.exp (-1) / (1 + 2 * Real.exp (-1)) := by
  unfold historyTilt tilt
  rw [historyTiltZ_eq, historyRef_apply, historyAction_apply_1,
    show -((1 : ℝ) * (1 : ℝ)) = -1 from by norm_num]
  have hd : (1 + 2 * Real.exp (-1)) ≠ 0 := ne_of_gt historyDenom_pos
  field_simp

theorem historyTilt_apply_2 :
    historyTilt 2 = Real.exp (-1) / (1 + 2 * Real.exp (-1)) := by
  unfold historyTilt tilt
  rw [historyTiltZ_eq, historyRef_apply, historyAction_apply_2,
    show -((1 : ℝ) * (1 : ℝ)) = -1 from by norm_num]
  have hd : (1 + 2 * Real.exp (-1)) ≠ 0 := ne_of_gt historyDenom_pos
  field_simp

/-- Constraint value of the witness projection: the tilted mean
action. -/
noncomputable def historyControlMass : ℝ :=
  2 * Real.exp (-1) / (1 + 2 * Real.exp (-1))

theorem historyControlMass_pos : 0 < historyControlMass :=
  div_pos (by positivity) historyDenom_pos

theorem historyControlMass_lt_one : historyControlMass < 1 := by
  unfold historyControlMass
  rw [div_lt_one historyDenom_pos]
  linarith

/-- The witness tilt places exactly the constraint value of mean
action. -/
theorem historyTilt_mean :
    meanAction historyTilt historyAction = historyControlMass := by
  unfold meanAction historyControlMass
  rw [Fin.sum_univ_three, historyTilt_apply_0, historyTilt_apply_1,
    historyTilt_apply_2, historyAction_apply_0, historyAction_apply_1,
    historyAction_apply_2]
  have hd : (1 + 2 * Real.exp (-1)) ≠ 0 := ne_of_gt historyDenom_pos
  field_simp
  ring

/-- The control law: the same mean action as the witness tilt, carried
on histories `0` and `1` alone. -/
noncomputable def historyControl : Fin 3 → ℝ :=
  fun g => if g = 0 then 1 - historyControlMass
    else if g = 1 then historyControlMass else 0

theorem historyControl_apply_0 :
    historyControl 0 = 1 - historyControlMass := by
  unfold historyControl
  rw [if_pos rfl]

theorem historyControl_apply_1 :
    historyControl 1 = historyControlMass := by
  unfold historyControl
  rw [if_neg (by decide : ¬ (1 : Fin 3) = 0), if_pos rfl]

theorem historyControl_apply_2 : historyControl 2 = 0 := by
  unfold historyControl
  rw [if_neg (by decide : ¬ (2 : Fin 3) = 0),
    if_neg (by decide : ¬ (2 : Fin 3) = 1)]

theorem historyControl_nonneg : ∀ g, 0 ≤ historyControl g := by
  intro g
  unfold historyControl
  by_cases h0 : g = 0
  · rw [if_pos h0]
    linarith [historyControlMass_lt_one]
  rw [if_neg h0]
  by_cases h1 : g = 1
  · rw [if_pos h1]
    exact le_of_lt historyControlMass_pos
  · rw [if_neg h1]

theorem historyControl_sum : ∑ g, historyControl g = 1 := by
  rw [Fin.sum_univ_three, historyControl_apply_0,
    historyControl_apply_1, historyControl_apply_2]
  ring

/-- The control law sits on the same constraint surface as the witness
tilt. -/
theorem historyControl_mean :
    meanAction historyControl historyAction = historyControlMass := by
  unfold meanAction
  rw [Fin.sum_univ_three, historyControl_apply_0,
    historyControl_apply_1, historyControl_apply_2,
    historyAction_apply_0, historyAction_apply_1,
    historyAction_apply_2]
  ring

theorem historyControl_matches_mean :
    meanAction historyControl historyAction
      = meanAction historyTilt historyAction := by
  rw [historyControl_mean, historyTilt_mean]

/-- The control law differs from the tilt: it vanishes on history `2`,
where the tilt is strictly positive. -/
theorem historyControl_ne_tilt :
    historyControl ≠ tilt historyRef historyAction 1 := by
  intro h
  have h2 := congrFun h 2
  rw [historyControl_apply_2] at h2
  exact absurd h2.symm
    (ne_of_gt (tilt_pos historyRef historyAction 1 historyRef_pos 2))

/-- **Witness minimum.**  On the three-history space, the witness tilt
scores at or below the control law on the variational objective. -/
theorem historyWitness_minimum :
    freeEnergy historyRef historyAction 1 historyTilt
      ≤ freeEnergy historyRef historyAction 1 historyControl :=
  tilt_minimizes_freeEnergy historyRef historyAction 1 historyControl
    historyControl_nonneg historyControl_sum historyRef_pos

/-- **Negative control.**  The control law shares the mean action of
the witness tilt and differs from it, and its relative entropy to the
reference is strictly larger: the constrained projection genuinely
selects the tilt on the witness. -/
theorem historyWitness_negative_control :
    kl historyTilt historyRef < kl historyControl historyRef :=
  tilt_constrained_strict historyRef historyAction 1 historyControl
    historyControl_nonneg historyControl_sum historyRef_pos
    historyControl_matches_mean historyControl_ne_tilt

end Witness

end OPH.InformationProjection

#print axioms OPH.InformationProjection.tilt_total
#print axioms OPH.InformationProjection.kl_tilt_split
#print axioms OPH.InformationProjection.tilt_minimizes_freeEnergy
#print axioms OPH.InformationProjection.tilt_strict_freeEnergy
#print axioms OPH.InformationProjection.tilt_unique_minimizer
#print axioms OPH.InformationProjection.tilt_maximizes_entropy_score
#print axioms OPH.InformationProjection.tilt_constrained_minimizer
#print axioms OPH.InformationProjection.tilt_constrained_strict
#print axioms OPH.InformationProjection.tilt_constrained_unique
#print axioms OPH.InformationProjection.tilt_pythagorean
#print axioms OPH.InformationProjection.pathGibbs_eq_tilt
#print axioms OPH.InformationProjection.pathGibbs_logZ_forced
#print axioms OPH.InformationProjection.historyLaw_constrained_minimizer
#print axioms OPH.InformationProjection.historyLaw_constrained_unique
#print axioms OPH.InformationProjection.pathGibbs_minimizer_unique
#print axioms OPH.InformationProjection.historyLaw_information_projection
#print axioms OPH.InformationProjection.historyWitness_minimum
#print axioms OPH.InformationProjection.historyWitness_negative_control
