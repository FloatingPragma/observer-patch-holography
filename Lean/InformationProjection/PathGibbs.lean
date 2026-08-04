import FiniteConditionalRepair

namespace OPH.InformationProjection

open OPH.Thermodynamics

/-!
# Conditional history projection: path Gibbs law and least action

The third axiom may be instantiated on a finite history type. This module
records the conditional mathematics of that instantiation: when a normalized
positive projected path law has the displayed exponential-tilt form relative
to a positive reference weight and
matches the declared moments, the information-projection identity holds. The
exponent defines an effective path action, and the projected law is antitone
in that action, so a supplied modal history minimizes it. A quantitative
zero-noise bound controls the total mass of every history whose action exceeds
the minimum by a declared gap, and a positive gap yields an actual vanishing
limit as inverse noise tends to infinity.

Scope boundary: the projected weight is a real exponential, and the
statements here concern a statistical selection over finite histories.
No existence theorem constructs the exponential law from OPH source data.
No complex amplitude, interference rule, physical action, clock, or continuum
history space is constructed; the reference path weight, constraint observables,
and moment-matching packet are named inputs, and their source derivation is a
separate obligation.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω] {N : ℕ}

/-- A finite history of length `N`: a record of `N + 1` states. -/
abbrev PathSpace (Ω : Type*) (N : ℕ) := Fin (N + 1) → Ω

section PathProjection

variable {A : Type*} [Fintype A]

/-- **History-level Pythagorean identity.** The state-side
information-projection identity of the thermodynamic package, read on
the finite path space: if the projected path law `τ` is an exponential
tilt of the positive reference weight `R` by the constraint observables, and
`P` matches the constrained moments, then
`D(P ‖ R) = D(P ‖ τ) + D(τ ‖ R)`. -/
theorem pathGibbs_pythagorean
    (P τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hP0 : ∀ γ, 0 ≤ P γ) (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hP1 : ∑ γ, P γ = 1) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (hmom : ∀ a, ∑ γ, P γ * F a γ = ∑ γ, τ γ * F a γ) :
    kl P R = kl P τ + kl τ R :=
  gibbs_pythagorean P τ R F lam logZ hP0 hτ hR hP1 hτ1 hlog hmom

/-- The projected path law minimizes the displayed relative-entropy
functional to the positive reference weight over the constrained moment
surface.  Normalization of `R` is not a premise of this declaration. -/
theorem pathGibbs_minimizer
    (P τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hP0 : ∀ γ, 0 ≤ P γ) (hτ : ∀ γ, 0 < τ γ) (hR : ∀ γ, 0 < R γ)
    (hP1 : ∑ γ, P γ = 1) (hτ1 : ∑ γ, τ γ = 1)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (hmom : ∀ a, ∑ γ, P γ * F a γ = ∑ γ, τ γ * F a γ) :
    kl τ R ≤ kl P R :=
  gibbs_minimizer P τ R F lam logZ hP0 hτ hR hP1 hτ1 hlog hmom

/-- The effective path action selected by the history projection:
minus the log reference weight plus the multiplier-weighted
constraints. -/
noncomputable def effectiveAction (R : PathSpace Ω N → ℝ)
    (F : A → PathSpace Ω N → ℝ) (lam : A → ℝ)
    (γ : PathSpace Ω N) : ℝ :=
  -Real.log (R γ) + ∑ a, lam a * F a γ

omit [Fintype Ω] [DecidableEq Ω] in
/-- The projected path law is the exponential of minus the effective
action, up to the fixed normalization. -/
theorem pathGibbs_eq_exp_neg_effectiveAction
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (γ : PathSpace Ω N) :
    τ γ = Real.exp (-(effectiveAction R F lam γ) - logZ) := by
  rw [← Real.exp_log (hτ γ)]
  congr 1
  rw [hlog γ]
  unfold effectiveAction
  ring

omit [Fintype Ω] [DecidableEq Ω] in
/-- **Exact finite action-order equivalence.** Under the projected path law, probability
order is the reverse of action order: `τ δ ≤ τ γ` exactly when the
effective action of `γ` is at most that of `δ`. -/
theorem prob_le_iff_action_le
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (γ δ : PathSpace Ω N) :
    τ δ ≤ τ γ ↔
      effectiveAction R F lam γ ≤ effectiveAction R F lam δ := by
  rw [pathGibbs_eq_exp_neg_effectiveAction τ R F lam logZ hτ hlog γ,
    pathGibbs_eq_exp_neg_effectiveAction τ R F lam logZ hτ hlog δ,
    Real.exp_le_exp]
  constructor <;> intro h <;> linarith

omit [Fintype Ω] [DecidableEq Ω] in
/-- **The modal history minimizes the effective action.** A history of
maximal projected probability has minimal effective action, as an
exact statement over the finite path space. -/
theorem modal_path_least_action
    (τ R : PathSpace Ω N → ℝ) (F : A → PathSpace Ω N → ℝ)
    (lam : A → ℝ) (logZ : ℝ)
    (hτ : ∀ γ, 0 < τ γ)
    (hlog : ∀ γ, Real.log (τ γ)
      = Real.log (R γ) - (∑ a, lam a * F a γ) - logZ)
    (γ : PathSpace Ω N) (hmax : ∀ δ, τ δ ≤ τ γ) :
    ∀ δ, effectiveAction R F lam γ ≤ effectiveAction R F lam δ :=
  fun δ =>
    (prob_le_iff_action_le τ R F lam logZ hτ hlog γ δ).mp (hmax δ)

end PathProjection

section ZeroNoise

variable {Γ : Type*} [Fintype Γ] [DecidableEq Γ]

/-- The inverse-noise family over a finite history space: the
normalized exponential of minus `β` times the action. -/
noncomputable def noiseFamily (S : Γ → ℝ) (β : ℝ) (γ : Γ) : ℝ :=
  Real.exp (-(β * S γ)) / ∑ δ, Real.exp (-(β * S δ))

omit [DecidableEq Γ] in
/-- **Quantitative baseline-gap bound.** At nonnegative inverse noise `β`,
the total mass of every history whose action exceeds the action of a supplied
baseline `γ₀` by at least `Δ` is bounded by the history count times
`exp (-(β * Δ))`.  No minimality premise is needed for this stronger
pointwise estimate; minimizer semantics enter only in the limit wrappers
below. -/
theorem noiseFamily_concentrates [Nonempty Γ]
    (S : Γ → ℝ) (β Δ : ℝ) (hβ : 0 ≤ β)
    (γ₀ : Γ) :
    ∑ γ ∈ Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ),
        noiseFamily S β γ
      ≤ (Fintype.card Γ : ℝ) * Real.exp (-(β * Δ)) := by
  have hZ : Real.exp (-(β * S γ₀)) ≤ ∑ δ, Real.exp (-(β * S δ)) :=
    Finset.single_le_sum (f := fun δ => Real.exp (-(β * S δ)))
      (fun δ _ => le_of_lt (Real.exp_pos _)) (Finset.mem_univ γ₀)
  have hZpos : 0 < ∑ δ, Real.exp (-(β * S δ)) :=
    lt_of_lt_of_le (Real.exp_pos _) hZ
  have hterm : ∀ γ ∈ Finset.univ.filter
      (fun γ => S γ₀ + Δ ≤ S γ),
      noiseFamily S β γ ≤ Real.exp (-(β * Δ)) := by
    intro γ hγ
    have hgap := (Finset.mem_filter.mp hγ).2
    unfold noiseFamily
    rw [div_le_iff₀ hZpos]
    calc Real.exp (-(β * S γ))
        ≤ Real.exp (-(β * Δ)) * Real.exp (-(β * S γ₀)) := by
          rw [← Real.exp_add]
          apply Real.exp_le_exp.mpr
          nlinarith [mul_le_mul_of_nonneg_left hgap hβ]
      _ ≤ Real.exp (-(β * Δ)) * ∑ δ, Real.exp (-(β * S δ)) := by
          apply mul_le_mul_of_nonneg_left hZ
          exact le_of_lt (Real.exp_pos _)
  calc ∑ γ ∈ Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ),
        noiseFamily S β γ
      ≤ ∑ _γ ∈ Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ),
        Real.exp (-(β * Δ)) := Finset.sum_le_sum hterm
    _ = (Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ)).card
        * Real.exp (-(β * Δ)) := by
        rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ (Fintype.card Γ : ℝ) * Real.exp (-(β * Δ)) := by
        apply mul_le_mul_of_nonneg_right
        · exact_mod_cast Finset.card_filter_le _ _
        · exact le_of_lt (Real.exp_pos _)

omit [DecidableEq Γ] in
/-- **Positive-gap zero-noise limit.** For every fixed `Δ > 0`, the
total probability of histories whose action is at least `Δ` above a
declared minimum tends to zero as the inverse noise tends to infinity.
This is the actual finite-history limit certified by the preceding
exponential tail bound; it does not identify a continuum path measure. -/
theorem noiseFamily_above_gap_mass_tendsto_zero [Nonempty Γ]
    (S : Γ → ℝ) (Δ : ℝ) (hΔ : 0 < Δ)
    (γ₀ : Γ) :
    Filter.Tendsto
      (fun β : ℝ =>
        ∑ γ ∈ Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ),
          noiseFamily S β γ)
      Filter.atTop (nhds 0) := by
  have hnonneg : ∀ β : ℝ,
      0 ≤ ∑ γ ∈ Finset.univ.filter (fun γ => S γ₀ + Δ ≤ S γ),
        noiseFamily S β γ := by
    intro β
    apply Finset.sum_nonneg
    intro γ _
    unfold noiseFamily
    exact div_nonneg (le_of_lt (Real.exp_pos _))
      (Finset.sum_nonneg fun δ _ => le_of_lt (Real.exp_pos _))
  have harg : Filter.Tendsto (fun β : ℝ => -(β * Δ))
      Filter.atTop Filter.atBot := by
    convert Filter.tendsto_id.const_mul_atTop_of_neg
      (neg_lt_zero.mpr hΔ) using 1
    funext β
    dsimp only [id_eq]
    ring
  have hexp : Filter.Tendsto (fun β : ℝ => Real.exp (- (β * Δ)))
      Filter.atTop (nhds 0) :=
    Real.tendsto_exp_atBot.comp harg
  have hbound : Filter.Tendsto
      (fun β : ℝ => (Fintype.card Γ : ℝ) * Real.exp (- (β * Δ)))
      Filter.atTop (nhds 0) := by
    simpa using hexp.const_mul (Fintype.card Γ : ℝ)
  apply squeeze_zero' (Filter.Eventually.of_forall hnonneg) _ hbound
  filter_upwards [Filter.eventually_ge_atTop (0 : ℝ)] with β hβ
  exact noiseFamily_concentrates S β Δ hβ γ₀

omit [DecidableEq Γ] in
/-- Every individual history whose action is strictly above the declared
minimum receives asymptotically zero inverse-noise probability. -/
theorem noiseFamily_strictly_higher_tendsto_zero [Nonempty Γ]
    (S : Γ → ℝ) (γ₀ : Γ)
    (γ : Γ) (hhigh : S γ₀ < S γ) :
    Filter.Tendsto (fun β : ℝ => noiseFamily S β γ)
      Filter.atTop (nhds 0) := by
  let Δ : ℝ := S γ - S γ₀
  have hΔ : 0 < Δ := sub_pos.mpr hhigh
  have htail := noiseFamily_above_gap_mass_tendsto_zero
    S Δ hΔ γ₀
  apply squeeze_zero' _ _ htail
  · apply Filter.Eventually.of_forall
    intro β
    unfold noiseFamily
    exact div_nonneg (le_of_lt (Real.exp_pos _))
      (Finset.sum_nonneg fun δ _ => le_of_lt (Real.exp_pos _))
  · apply Filter.Eventually.of_forall
    intro β
    apply Finset.single_le_sum
    · intro δ _
      unfold noiseFamily
      exact div_nonneg (le_of_lt (Real.exp_pos _))
        (Finset.sum_nonneg fun η _ => le_of_lt (Real.exp_pos _))
    · apply Finset.mem_filter.mpr
      constructor
      · exact Finset.mem_univ γ
      · dsimp only [Δ]
        linarith

omit [DecidableEq Γ] in
/-- **Concentration on the finite minimizer set.** The total
inverse-noise probability of all histories whose action differs from a declared minimum
tends to zero.  Unlike a fixed-gap tail statement, this quantifies over every
nonminimal history; finiteness permits the result to be assembled from the
individual strict-gap limits. -/
theorem noiseFamily_nonminimal_mass_tendsto_zero [Nonempty Γ]
    (S : Γ → ℝ) (γ₀ : Γ) (hmin : ∀ δ, S γ₀ ≤ S δ) :
    Filter.Tendsto
      (fun β : ℝ =>
        ∑ γ ∈ Finset.univ.filter (fun γ => S γ ≠ S γ₀),
          noiseFamily S β γ)
      Filter.atTop (nhds 0) := by
  simpa using tendsto_finset_sum
    (Finset.univ.filter (fun γ => S γ ≠ S γ₀))
    (fun γ hγ => noiseFamily_strictly_higher_tendsto_zero
      S γ₀ γ (lt_of_le_of_ne (hmin γ)
        (Ne.symm (Finset.mem_filter.mp hγ).2)))

omit [DecidableEq Γ] in
/-- The inverse-noise family is a probability law: nonnegative entries
of unit total mass, so the concentration bound controls genuine
probability. -/
theorem noiseFamily_total [Nonempty Γ] (S : Γ → ℝ) (β : ℝ) :
    ∑ γ, noiseFamily S β γ = 1 := by
  unfold noiseFamily
  rw [← Finset.sum_div]
  apply div_self
  have : 0 < ∑ δ, Real.exp (-(β * S δ)) :=
    Finset.sum_pos (fun δ _ => Real.exp_pos _) Finset.univ_nonempty
  exact ne_of_gt this

end ZeroNoise

end OPH.InformationProjection

#print axioms OPH.InformationProjection.pathGibbs_pythagorean
#print axioms OPH.InformationProjection.pathGibbs_minimizer
#print axioms OPH.InformationProjection.pathGibbs_eq_exp_neg_effectiveAction
#print axioms OPH.InformationProjection.prob_le_iff_action_le
#print axioms OPH.InformationProjection.modal_path_least_action
#print axioms OPH.InformationProjection.noiseFamily_concentrates
#print axioms OPH.InformationProjection.noiseFamily_above_gap_mass_tendsto_zero
#print axioms OPH.InformationProjection.noiseFamily_strictly_higher_tendsto_zero
#print axioms OPH.InformationProjection.noiseFamily_nonminimal_mass_tendsto_zero
#print axioms OPH.InformationProjection.noiseFamily_total
