import Mathlib

namespace OPH.Thermodynamics

/-!
# Finite conditional repair and the transition-side second law

This module proves the finite classical core of the OPH thermodynamic
completion program on one regulator stage:

* the pointwise and summed nonnegativity of finite relative entropy;
* the finite log-sum inequality;
* the data-processing inequality for finite stochastic kernels;
* the conditional-resampling theorem: the transition-side information
  projection onto the fibre of a complete repaired visible datum is the
  normalized restriction of the common reference law, and this kernel is
  stochastic, idempotent, reversible, stationary, and fixes every
  fibre-measurable observable;
* the transition-side second law: relative entropy to the reference law
  is nonincreasing under the repair kernel;
* the Gibbs information-projection Pythagorean identity and the unique
  constrained minimizer;
* the finite low-temperature concentration bound and its limit.

Imported premises, not established here:

* the identification of the reference law with a source-derived quotient
  weight shared by the state and transition optimizers (the
  common-reference receipt);
* the identification of a fibre-constant observable with physical energy
  and a calibrated clock (the energy-clock receipt);
* every continuum, refinement-uniform, and physical-realization
  statement.

The deterministic strict-descent normalizer is a different map from the
repair kernel of this file, and no theorem here asserts an entropy
inequality for it; the Python certificate carries the explicit
counterexample.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]

/-- One term of finite relative entropy, with the `0 log 0 = 0`
convention. -/
noncomputable def klTerm (q p : ℝ) : ℝ :=
  if q = 0 then 0 else q * Real.log (q / p)

/-- Finite relative entropy `D(q ‖ p)` over a finite index type. -/
noncomputable def kl (q p : Ω → ℝ) : ℝ :=
  ∑ x, klTerm (q x) (p x)

@[simp] theorem klTerm_zero (p : ℝ) : klTerm 0 p = 0 := by
  simp [klTerm]

theorem klTerm_of_pos {q : ℝ} (hq : q ≠ 0) (p : ℝ) :
    klTerm q p = q * Real.log (q / p) := by
  simp [klTerm, hq]

/-- Pointwise lower bound `klTerm q p ≥ q - p` under absolute
continuity. -/
theorem klTerm_ge {q p : ℝ} (hq : 0 ≤ q) (hp : 0 ≤ p)
    (habs : p = 0 → q = 0) : q - p ≤ klTerm q p := by
  rcases eq_or_lt_of_le hq with hq0 | hqpos
  · simp [← hq0]
    linarith
  rcases eq_or_lt_of_le hp with hp0 | hppos
  · exact absurd (habs hp0.symm) (ne_of_gt hqpos)
  have hlog : Real.log (p / q) ≤ p / q - 1 :=
    Real.log_le_sub_one_of_pos (by positivity)
  have hflip : Real.log (q / p) = -Real.log (p / q) := by
    rw [← Real.log_inv]
    congr 1
    field_simp
  have : 1 - p / q ≤ Real.log (q / p) := by
    rw [hflip]
    linarith
  have hmul := mul_le_mul_of_nonneg_left this (le_of_lt hqpos)
  rw [klTerm_of_pos (ne_of_gt hqpos)]
  have hexp : q * (1 - p / q) = q - p := by
    field_simp
  nlinarith [hmul, hexp]

/-- Relative entropy is nonnegative whenever the two vectors share their
total mass and the reference dominates. -/
theorem kl_nonneg (q p : Ω → ℝ) (hq : ∀ x, 0 ≤ q x) (hp : ∀ x, 0 ≤ p x)
    (habs : ∀ x, p x = 0 → q x = 0)
    (hsum : ∑ x, q x = ∑ x, p x) : 0 ≤ kl q p := by
  have hterm : ∀ x ∈ Finset.univ, q x - p x ≤ klTerm (q x) (p x) :=
    fun x _ => klTerm_ge (hq x) (hp x) (habs x)
  have hsum' := Finset.sum_le_sum hterm
  rw [Finset.sum_sub_distrib] at hsum'
  unfold kl
  linarith [hsum']

/-- Scaling: `klTerm (c q) (c p) = c klTerm q p` for `c ≥ 0`. -/
theorem klTerm_smul {c : ℝ} (hc : 0 ≤ c) (q p : ℝ) :
    klTerm (c * q) (c * p) = c * klTerm q p := by
  rcases eq_or_lt_of_le hc with hc0 | hcpos
  · simp [← hc0]
  by_cases hq : q = 0
  · simp [hq]
  have hcq : c * q ≠ 0 := mul_ne_zero (ne_of_gt hcpos) hq
  rw [klTerm_of_pos hcq, klTerm_of_pos hq]
  have : c * q / (c * p) = q / p := by
    rw [mul_div_mul_left _ _ (ne_of_gt hcpos)]
  rw [this]
  ring

/-- The finite log-sum inequality:
`klTerm (∑ a) (∑ b) ≤ ∑ klTerm (a x) (b x)` for nonnegative entries with
`b`-support dominating `a`-support. -/
theorem logSum (a b : Ω → ℝ) (ha : ∀ x, 0 ≤ a x) (hb : ∀ x, 0 ≤ b x)
    (habs : ∀ x, b x = 0 → a x = 0) :
    klTerm (∑ x, a x) (∑ x, b x) ≤ ∑ x, klTerm (a x) (b x) := by
  set A := ∑ x, a x with hA
  set Bt := ∑ x, b x with hB
  have hApos : 0 ≤ A := Finset.sum_nonneg fun x _ => ha x
  rcases eq_or_lt_of_le hApos with hA0 | hApos'
  · have hzero : ∀ x ∈ Finset.univ, a x = 0 := by
      intro x _
      have := (Finset.sum_eq_zero_iff_of_nonneg
        (fun x _ => ha x)).mp hA0.symm
      exact this x (Finset.mem_univ x)
    have hz : ∑ x, klTerm (a x) (b x) = 0 :=
      Finset.sum_eq_zero fun x hx => by rw [hzero x hx, klTerm_zero]
    rw [← hA0, klTerm_zero, hz]
  have hBpos : 0 < Bt := by
    rcases eq_or_lt_of_le (Finset.sum_nonneg fun x _ => hb x) with h0 | h
    · exfalso
      have hzero := (Finset.sum_eq_zero_iff_of_nonneg
        (fun x _ => hb x)).mp h0.symm
      have : ∀ x ∈ Finset.univ, a x = 0 := fun x hx =>
        habs x (hzero x hx)
      have : A = 0 := Finset.sum_eq_zero this
      exact absurd this (ne_of_gt hApos')
    · exact h
  -- Termwise identity:
  --   klTerm (a x) (b x)
  --     = A * klTerm (a x / A) (b x / Bt) + a x * Real.log (A / Bt).
  have hterm : ∀ x,
      klTerm (a x) (b x)
        = A * klTerm (a x / A) (b x / Bt)
          + a x * Real.log (A / Bt) := by
    intro x
    by_cases hax : a x = 0
    · simp [hax]
    have haxpos : 0 < a x := lt_of_le_of_ne (ha x) (Ne.symm hax)
    have hbx : b x ≠ 0 := fun h => hax (habs x h)
    have hbxpos : 0 < b x := lt_of_le_of_ne (hb x) (Ne.symm hbx)
    have haA : a x / A ≠ 0 :=
      ne_of_gt (div_pos haxpos hApos')
    rw [klTerm_of_pos hax, klTerm_of_pos haA]
    have hratio : a x / A / (b x / Bt) = a x / b x * (Bt / A) := by
      field_simp
    rw [hratio]
    have hlog1 : Real.log (a x / b x * (Bt / A))
        = Real.log (a x / b x) + Real.log (Bt / A) :=
      Real.log_mul (ne_of_gt (div_pos haxpos hbxpos))
        (ne_of_gt (div_pos hBpos hApos'))
    have hlog2 : Real.log (A / Bt) = -Real.log (Bt / A) := by
      rw [← Real.log_inv]
      congr 1
      field_simp
    rw [hlog1, hlog2]
    field_simp
    ring
  have hsum : ∑ x, klTerm (a x) (b x)
      = A * kl (fun x => a x / A) (fun x => b x / Bt)
        + A * Real.log (A / Bt) := by
    rw [Finset.sum_congr rfl fun x _ => hterm x, Finset.sum_add_distrib,
      ← Finset.mul_sum, ← Finset.sum_mul]
    rfl
  have hklnn : 0 ≤ kl (fun x => a x / A) (fun x => b x / Bt) := by
    apply kl_nonneg
    · exact fun x => div_nonneg (ha x) (le_of_lt hApos')
    · exact fun x => div_nonneg (hb x) (le_of_lt hBpos)
    · intro x hx
      have : b x = 0 := by
        by_contra hbx
        have : 0 < b x / Bt :=
          div_pos (lt_of_le_of_ne (hb x) (Ne.symm hbx)) hBpos
        exact absurd hx (ne_of_gt this)
      rw [habs x this]
      simp
    · rw [← Finset.sum_div, ← Finset.sum_div, ← hA, ← hB,
        div_self (ne_of_gt hApos'), div_self (ne_of_gt hBpos)]
  rw [hsum, klTerm_of_pos (ne_of_gt hApos')]
  nlinarith [mul_nonneg (le_of_lt hApos') hklnn]

/-- Pushforward of a finite distribution through a kernel. -/
noncomputable def push (p : Ω → ℝ) (K : Ω → Ω → ℝ) (y : Ω) : ℝ :=
  ∑ x, p x * K x y

/-- **Data processing.** Relative entropy is nonincreasing under every
finite stochastic kernel. -/
theorem kl_push_le (p q : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x) (hq : ∀ x, 0 ≤ q x)
    (habs : ∀ x, q x = 0 → p x = 0)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1) :
    kl (push p K) (push q K) ≤ kl p q := by
  have hstep : ∀ y,
      klTerm (push p K y) (push q K y)
        ≤ ∑ x, klTerm (p x * K x y) (q x * K x y) := by
    intro y
    apply logSum
    · exact fun x => mul_nonneg (hp x) (hK0 x y)
    · exact fun x => mul_nonneg (hq x) (hK0 x y)
    · intro x hx
      rcases mul_eq_zero.mp hx with h | h
      · rw [habs x h]
        ring
      · rw [h]
        ring
  have hsum : ∑ y, ∑ x, klTerm (p x * K x y) (q x * K x y)
      = kl p q := by
    rw [Finset.sum_comm]
    unfold kl
    apply Finset.sum_congr rfl
    intro x _
    have : ∀ y ∈ Finset.univ,
        klTerm (p x * K x y) (q x * K x y)
          = K x y * klTerm (p x) (q x) := by
      intro y _
      rw [mul_comm (p x) (K x y), mul_comm (q x) (K x y)]
      exact klTerm_smul (hK0 x y) (p x) (q x)
    rw [Finset.sum_congr rfl this, ← Finset.sum_mul, hK1 x, one_mul]
  calc kl (push p K) (push q K)
      ≤ ∑ y, ∑ x, klTerm (p x * K x y) (q x * K x y) :=
        Finset.sum_le_sum fun y _ => hstep y
    _ = kl p q := hsum

section ConditionalRepair

variable (π : Ω → ℝ) (b : Ω → B)

/-- Total reference mass of the repaired-visible fibre through `x`. -/
noncomputable def fiberMass (x : Ω) : ℝ :=
  ∑ y ∈ Finset.univ.filter (fun y => b y = b x), π y

/-- The conditional-resampling (heat-bath) repair kernel: resample the
unresolved variables from the reference law inside the fibre of the
complete repaired visible datum. -/
noncomputable def heatBath (x y : Ω) : ℝ :=
  if b y = b x then π y / fiberMass π b x else 0

variable {π b}

theorem fiberMass_pos (hπ : ∀ x, 0 < π x) (x : Ω) :
    0 < fiberMass π b x := by
  unfold fiberMass
  apply Finset.sum_pos'
  · exact fun y _ => le_of_lt (hπ y)
  · exact ⟨x, Finset.mem_filter.mpr ⟨Finset.mem_univ x, rfl⟩, hπ x⟩

theorem fiberMass_congr {x x' : Ω} (h : b x = b x') :
    fiberMass π b x = fiberMass π b x' := by
  unfold fiberMass
  congr 1
  apply Finset.filter_congr
  intro y _
  simp [h]

/-- Row normalization: the repair kernel is stochastic. -/
theorem heatBath_row_sum (hπ : ∀ x, 0 < π x) (x : Ω) :
    ∑ y, heatBath π b x y = 1 := by
  unfold heatBath
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, ← Finset.sum_div]
  exact div_self (ne_of_gt (fiberMass_pos hπ x))

/-- Idempotence: one full fibre resampling reaches conditional
equilibrium. -/
theorem heatBath_idempotent (hπ : ∀ x, 0 < π x) (x y : Ω) :
    ∑ z, heatBath π b x z * heatBath π b z y = heatBath π b x y := by
  unfold heatBath
  by_cases hxy : b y = b x
  · have : ∀ z, (if b z = b x then π z / fiberMass π b x else 0)
        * (if b y = b z then π y / fiberMass π b z else 0)
        = if b z = b x then
            π z / fiberMass π b x * (π y / fiberMass π b x) else 0 := by
      intro z
      by_cases hz : b z = b x
      · have hyz : b y = b z := hz ▸ hxy
        simp only [hz, hyz, if_true]
        rw [fiberMass_congr hz]
      · simp [hz]
    rw [Finset.sum_congr rfl fun z _ => this z, Finset.sum_ite,
      Finset.sum_const_zero, add_zero, if_pos hxy]
    rw [← Finset.sum_mul, ← Finset.sum_div]
    have hMdef : (∑ i with b i = b x, π i) = fiberMass π b x := rfl
    rw [hMdef, div_self (ne_of_gt (fiberMass_pos hπ x)), one_mul]
  · have : ∀ z, (if b z = b x then π z / fiberMass π b x else 0)
        * (if b y = b z then π y / fiberMass π b z else 0) = 0 := by
      intro z
      by_cases hz : b z = b x
      · have : ¬ b y = b z := fun h => hxy (h.trans hz)
        simp [this]
      · simp [hz]
    rw [Finset.sum_congr rfl fun z _ => this z, Finset.sum_const_zero,
      if_neg hxy]

/-- Detailed balance with respect to the reference law. -/
theorem heatBath_detailedBalance (x y : Ω) :
    π x * heatBath π b x y = π y * heatBath π b y x := by
  unfold heatBath
  by_cases h : b y = b x
  · rw [if_pos h, if_pos h.symm, fiberMass_congr h.symm]
    ring
  · rw [if_neg h, if_neg (fun h' => h h'.symm)]
    ring

/-- Stationarity: the reference law is invariant under the repair
kernel. -/
theorem heatBath_stationary (hπ : ∀ x, 0 < π x) (y : Ω) :
    push π (heatBath π b) y = π y := by
  unfold push
  have : ∀ x, π x * heatBath π b x y = π y * heatBath π b y x :=
    fun x => heatBath_detailedBalance x y
  rw [Finset.sum_congr rfl fun x _ => this x, ← Finset.mul_sum,
    heatBath_row_sum hπ y, mul_one]

/-- Every observable measurable through the repaired visible datum is
fixed by the repair kernel. -/
theorem heatBath_fixes_fiberObservable (hπ : ∀ x, 0 < π x)
    (Q : Ω → ℝ) (hQ : ∀ x y, b x = b y → Q x = Q y) (x : Ω) :
    ∑ y, heatBath π b x y * Q y = Q x := by
  unfold heatBath
  have : ∀ y, (if b y = b x then π y / fiberMass π b x else 0) * Q y
      = (if b y = b x then π y / fiberMass π b x else 0) * Q x := by
    intro y
    by_cases h : b y = b x
    · rw [hQ y x h]
    · simp [h]
  rw [Finset.sum_congr rfl fun y _ => this y, ← Finset.sum_mul]
  have := heatBath_row_sum (π := π) (b := b) hπ x
  unfold heatBath at this
  rw [this, one_mul]

/-- **Conditional-resampling theorem.** Among all distributions supported
on the fibre of the repaired visible datum, the heat-bath row minimizes
relative entropy to the reference law: transition-side information
projection selects conditional resampling. -/
theorem heatBath_row_optimal (hπ : ∀ x, 0 < π x) (x : Ω)
    (r : Ω → ℝ) (hr0 : ∀ y, 0 ≤ r y) (hr1 : ∑ y, r y = 1)
    (hsupp : ∀ y, b y ≠ b x → r y = 0) :
    kl (heatBath π b x) π ≤ kl r π := by
  -- Decomposition: for fibre-supported r,
  --   kl r π = kl r (heatBath row) - log (fiberMass).
  have hM := fiberMass_pos (π := π) (b := b) hπ x
  have hdecomp : ∀ s : Ω → ℝ, (∀ y, 0 ≤ s y) →
      (∀ y, b y ≠ b x → s y = 0) →
      kl s π = kl s (heatBath π b x)
        - (∑ y, s y) * Real.log (fiberMass π b x) := by
    intro s hs0 hsupp'
    unfold kl heatBath
    rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro y _
    by_cases hy : b y = b x
    · by_cases hsy : s y = 0
      · simp [hsy]
      rw [klTerm_of_pos hsy, if_pos hy, klTerm_of_pos hsy]
      have hπy := hπ y
      have hspos : 0 < s y := lt_of_le_of_ne (hs0 y) (Ne.symm hsy)
      have : s y / (π y / fiberMass π b x)
          = s y / π y * fiberMass π b x := by
        field_simp
      rw [this, Real.log_mul
        (ne_of_gt (div_pos hspos hπy)) (ne_of_gt hM)]
      ring_nf
    · simp [hsupp' y hy]
  have hrowsupp : ∀ y, b y ≠ b x → heatBath π b x y = 0 := by
    intro y hy
    unfold heatBath
    rw [if_neg hy]
  have hrow0 : ∀ y, 0 ≤ heatBath π b x y := by
    intro y
    unfold heatBath
    by_cases h : b y = b x
    · rw [if_pos h]
      exact le_of_lt (div_pos (hπ y) hM)
    · rw [if_neg h]
  have h1 := hdecomp (heatBath π b x) hrow0 hrowsupp
  have h2 := hdecomp r hr0 hsupp
  rw [heatBath_row_sum hπ x] at h1
  rw [hr1] at h2
  have hself : kl (heatBath π b x) (heatBath π b x) = 0 := by
    unfold kl
    apply Finset.sum_eq_zero
    intro y _
    by_cases h : heatBath π b x y = 0
    · simp [h]
    · rw [klTerm_of_pos h, div_self h, Real.log_one, mul_zero]
  have hklnn : 0 ≤ kl r (heatBath π b x) := by
    apply kl_nonneg
    · exact hr0
    · exact hrow0
    · intro y hy
      by_contra hry
      have hbY : b y = b x := by
        by_contra hb
        exact hry (hsupp y hb)
      unfold heatBath at hy
      rw [if_pos hbY] at hy
      have : 0 < π y / fiberMass π b x := div_pos (hπ y) hM
      exact absurd hy (ne_of_gt this)
    · rw [hr1]
      rw [heatBath_row_sum hπ x]
  rw [h1, h2, hself]
  linarith

/-- **Transition-side second law.** Relative entropy to the reference
law is nonincreasing under the conditional repair kernel. -/
theorem heatBath_secondLaw (hπ : ∀ x, 0 < π x)
    (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    kl (push p (heatBath π b)) π ≤ kl p π := by
  have h := kl_push_le p π (heatBath π b) hp
    (fun x => le_of_lt (hπ x)) (fun x hx => absurd hx (ne_of_gt (hπ x)))
    (fun x y => by
      unfold heatBath
      by_cases h : b y = b x
      · rw [if_pos h]
        exact le_of_lt (div_pos (hπ y) (fiberMass_pos hπ x))
      · rw [if_neg h])
    (heatBath_row_sum hπ)
  have hstat : push π (heatBath π b) = π := by
    funext y
    exact heatBath_stationary hπ y
  rwa [hstat] at h

end ConditionalRepair

section GibbsProjection

/-- **Information-projection Pythagorean identity.** If `log τ` differs
from `log σ` by a linear combination of constraint observables plus a
constant, and `ρ` and `τ` share the constrained moments, then
`D(ρ ‖ σ) = D(ρ ‖ τ) + D(τ ‖ σ)`. -/
theorem gibbs_pythagorean {A : Type*} [Fintype A]
    (ρ τ σ : Ω → ℝ) (Q : A → Ω → ℝ) (lam : A → ℝ) (logZ : ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hτ : ∀ x, 0 < τ x) (hσ : ∀ x, 0 < σ x)
    (hρ1 : ∑ x, ρ x = 1) (hτ1 : ∑ x, τ x = 1)
    (hlog : ∀ x, Real.log (τ x)
      = Real.log (σ x) - (∑ a, lam a * Q a x) - logZ)
    (hmom : ∀ a, ∑ x, ρ x * Q a x = ∑ x, τ x * Q a x) :
    kl ρ σ = kl ρ τ + kl τ σ := by
  -- For any s ≥ 0 with total mass 1 and matching moments,
  --   kl s σ - kl s τ = -(∑ a, lam a * moment a) - logZ.
  have hgap : ∀ s : Ω → ℝ, (∀ x, 0 ≤ s x) → (∑ x, s x = 1) →
      kl s σ - kl s τ
        = -(∑ a, lam a * (∑ x, s x * Q a x)) - logZ := by
    intro s hs0 hs1
    unfold kl
    rw [← Finset.sum_sub_distrib]
    have hterm : ∀ x, klTerm (s x) (σ x) - klTerm (s x) (τ x)
        = s x * (-(∑ a, lam a * Q a x) - logZ) := by
      intro x
      by_cases hsx : s x = 0
      · simp [hsx]
      rw [klTerm_of_pos hsx, klTerm_of_pos hsx]
      have hspos : 0 < s x := lt_of_le_of_ne (hs0 x) (Ne.symm hsx)
      rw [Real.log_div hsx (ne_of_gt (hσ x)),
        Real.log_div hsx (ne_of_gt (hτ x))]
      have := hlog x
      nlinarith [this]
    rw [Finset.sum_congr rfl fun x _ => hterm x]
    have hsplit : ∀ x, s x * (-(∑ a, lam a * Q a x) - logZ)
        = -(∑ a, lam a * (s x * Q a x)) - s x * logZ := by
      intro x
      have h1 : s x * ∑ a, lam a * Q a x
          = ∑ a, lam a * (s x * Q a x) := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun a _ => by ring
      linear_combination -h1
    have hswap : ∑ x, -(∑ a, lam a * (s x * Q a x))
        = -(∑ a, lam a * ∑ x, s x * Q a x) := by
      rw [Finset.sum_neg_distrib]
      congr 1
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun a _ => (Finset.mul_sum _ _ _).symm
    rw [Finset.sum_congr rfl fun x _ => hsplit x,
      Finset.sum_sub_distrib, hswap, ← Finset.sum_mul, hs1, one_mul]
  have h1 := hgap ρ hρ0 hρ1
  have h2 := hgap τ (fun x => le_of_lt (hτ x)) hτ1
  have hmomsum : ∑ a, lam a * (∑ x, ρ x * Q a x)
      = ∑ a, lam a * (∑ x, τ x * Q a x) := by
    apply Finset.sum_congr rfl
    intro a _
    rw [hmom a]
  have hττ : kl τ τ = 0 := by
    unfold kl
    apply Finset.sum_eq_zero
    intro x _
    rw [klTerm_of_pos (ne_of_gt (hτ x)), div_self (ne_of_gt (hτ x)),
      Real.log_one, mul_zero]
  have := h2
  rw [hττ] at this
  rw [hmomsum] at h1
  linarith

/-- The Gibbs state is the unique-in-value constrained minimizer of
relative entropy to the prior. -/
theorem gibbs_minimizer {A : Type*} [Fintype A]
    (ρ τ σ : Ω → ℝ) (Q : A → Ω → ℝ) (lam : A → ℝ) (logZ : ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hτ : ∀ x, 0 < τ x) (hσ : ∀ x, 0 < σ x)
    (hρ1 : ∑ x, ρ x = 1) (hτ1 : ∑ x, τ x = 1)
    (hlog : ∀ x, Real.log (τ x)
      = Real.log (σ x) - (∑ a, lam a * Q a x) - logZ)
    (hmom : ∀ a, ∑ x, ρ x * Q a x = ∑ x, τ x * Q a x) :
    kl τ σ ≤ kl ρ σ := by
  rw [gibbs_pythagorean ρ τ σ Q lam logZ hρ0 hτ hσ hρ1 hτ1 hlog hmom]
  have : 0 ≤ kl ρ τ := by
    apply kl_nonneg
    · exact hρ0
    · exact fun x => le_of_lt (hτ x)
    · exact fun x hx => absurd hx (ne_of_gt (hτ x))
    · rw [hρ1, hτ1]
  linarith

end GibbsProjection

section LowTemperature

variable (E : Ω → ℝ)

/-- Unnormalized Gibbs weight. -/
noncomputable def gibbsWeight (beta : ℝ) (x : Ω) : ℝ :=
  Real.exp (-beta * E x)

/-- Partition function. -/
noncomputable def partitionZ (beta : ℝ) : ℝ :=
  ∑ x, gibbsWeight E beta x

theorem partitionZ_pos (beta : ℝ) [Nonempty Ω] :
    0 < partitionZ E beta :=
  Finset.sum_pos (fun x _ => Real.exp_pos _) Finset.univ_nonempty

/-- Excited-state Gibbs mass above a spectral gap. -/
noncomputable def excitedMass (beta E0 : ℝ) : ℝ :=
  (∑ x ∈ Finset.univ.filter (fun x => E x ≠ E0),
    gibbsWeight E beta x) / partitionZ E beta

/-- **Finite low-temperature concentration.** With ground energy `E0`,
ground degeneracy `g₀`, and gap `Δ > 0`, the excited Gibbs mass obeys
`q_β ≤ (d - g₀)/g₀ · exp (-β Δ)`. -/
theorem excitedMass_le [Nonempty Ω] (beta E0 Δ : ℝ)
    (hβ : 0 ≤ beta) (hΔ : 0 < Δ)
    (hground : ∀ x, E x = E0 ∨ E0 + Δ ≤ E x)
    (hnonempty : (Finset.univ.filter (fun x => E x = E0)).Nonempty) :
    excitedMass E beta E0
      ≤ ((Finset.univ.filter (fun x => E x ≠ E0)).card : ℝ)
        / ((Finset.univ.filter (fun x => E x = E0)).card : ℝ)
        * Real.exp (-beta * Δ) := by
  set G := Finset.univ.filter (fun x => E x = E0) with hG
  set X := Finset.univ.filter (fun x => E x ≠ E0) with hX
  have hGcard : 0 < (G.card : ℝ) := by
    exact_mod_cast Finset.card_pos.mpr hnonempty
  -- Numerator bound: each excited weight is at most exp(-β(E0+Δ)).
  have hnum : ∑ x ∈ X, gibbsWeight E beta x
      ≤ (X.card : ℝ) * Real.exp (-beta * (E0 + Δ)) := by
    have : ∀ x ∈ X, gibbsWeight E beta x
        ≤ Real.exp (-beta * (E0 + Δ)) := by
      intro x hx
      have hEx : E0 + Δ ≤ E x := by
        rcases hground x with h | h
        · exact absurd h (Finset.mem_filter.mp hx).2
        · exact h
      unfold gibbsWeight
      apply Real.exp_le_exp.mpr
      nlinarith
    calc ∑ x ∈ X, gibbsWeight E beta x
        ≤ ∑ _x ∈ X, Real.exp (-beta * (E0 + Δ)) :=
          Finset.sum_le_sum this
      _ = (X.card : ℝ) * Real.exp (-beta * (E0 + Δ)) := by
          rw [Finset.sum_const, nsmul_eq_mul]
  -- Denominator bound: Z is at least g₀ exp(-β E0).
  have hden : (G.card : ℝ) * Real.exp (-beta * E0)
      ≤ partitionZ E beta := by
    unfold partitionZ
    have hsub : ∑ x ∈ G, gibbsWeight E beta x
        = (G.card : ℝ) * Real.exp (-beta * E0) := by
      have : ∀ x ∈ G, gibbsWeight E beta x
          = Real.exp (-beta * E0) := by
        intro x hx
        unfold gibbsWeight
        rw [(Finset.mem_filter.mp hx).2]
      rw [Finset.sum_congr rfl this, Finset.sum_const, nsmul_eq_mul]
    calc (G.card : ℝ) * Real.exp (-beta * E0)
        = ∑ x ∈ G, gibbsWeight E beta x := hsub.symm
      _ ≤ ∑ x, gibbsWeight E beta x := by
          apply Finset.sum_le_sum_of_subset_of_nonneg
            (Finset.filter_subset _ _)
          intro x _ _
          exact le_of_lt (Real.exp_pos _)
  have hZpos := partitionZ_pos E beta
  unfold excitedMass
  rw [div_le_iff₀ hZpos]
  have hchain : ∑ x ∈ X, gibbsWeight E beta x
      ≤ (X.card : ℝ) / (G.card : ℝ) * Real.exp (-beta * Δ)
        * ((G.card : ℝ) * Real.exp (-beta * E0)) := by
    have hfold : (X.card : ℝ) / (G.card : ℝ) * Real.exp (-beta * Δ)
        * ((G.card : ℝ) * Real.exp (-beta * E0))
        = (X.card : ℝ)
          * (Real.exp (-beta * Δ) * Real.exp (-beta * E0)) := by
      field_simp
    have hexp : Real.exp (-beta * Δ) * Real.exp (-beta * E0)
        = Real.exp (-beta * (E0 + Δ)) := by
      rw [← Real.exp_add]
      ring_nf
    rw [hfold, hexp]
    exact hnum
  calc ∑ x ∈ X, gibbsWeight E beta x
      ≤ (X.card : ℝ) / (G.card : ℝ) * Real.exp (-beta * Δ)
        * ((G.card : ℝ) * Real.exp (-beta * E0)) := hchain
    _ ≤ (X.card : ℝ) / (G.card : ℝ) * Real.exp (-beta * Δ)
        * partitionZ E beta := by
        apply mul_le_mul_of_nonneg_left hden
        positivity

/-- Quantitative threshold: beyond `β₀ = log C / Δ` with
`C = (d-g₀)/(g₀ ε)`, the excited mass is below `ε`. Stated through the
bound, so no limit machinery enters. -/
theorem excitedMass_lt_of_beta_large [Nonempty Ω] (beta E0 Δ eps : ℝ)
    (hβ : 0 ≤ beta) (hΔ : 0 < Δ) (heps : 0 < eps)
    (hground : ∀ x, E x = E0 ∨ E0 + Δ ≤ E x)
    (hnonempty : (Finset.univ.filter (fun x => E x = E0)).Nonempty)
    (hthreshold :
      Real.log (((Finset.univ.filter (fun x => E x ≠ E0)).card : ℝ)
        / ((Finset.univ.filter (fun x => E x = E0)).card : ℝ) / eps)
        < beta * Δ) :
    excitedMass E beta E0 < eps := by
  have hbound := excitedMass_le E beta E0 Δ hβ hΔ hground hnonempty
  set C := ((Finset.univ.filter (fun x => E x ≠ E0)).card : ℝ)
    / ((Finset.univ.filter (fun x => E x = E0)).card : ℝ) with hC
  have hCnn : 0 ≤ C := by
    apply div_nonneg
    · exact Nat.cast_nonneg _
    · exact Nat.cast_nonneg _
  have hkey : C * Real.exp (-beta * Δ) < eps := by
    rcases eq_or_lt_of_le hCnn with hC0 | hCpos
    · rw [← hC0, zero_mul]
      exact heps
    · have hlog : Real.log (C / eps) < beta * Δ := hthreshold
      have : C / eps < Real.exp (beta * Δ) := by
        calc C / eps = Real.exp (Real.log (C / eps)) := by
              rw [Real.exp_log (div_pos hCpos heps)]
          _ < Real.exp (beta * Δ) := Real.exp_lt_exp.mpr hlog
      have hexp : 0 < Real.exp (beta * Δ) := Real.exp_pos _
      rw [div_lt_iff₀ heps] at this
      rw [neg_mul, Real.exp_neg]
      rw [mul_inv_lt_iff₀ hexp]
      linarith
  linarith [hbound]

end LowTemperature

section ZerothLaw

variable (E : Ω → ℝ)

/-- Normalized Gibbs state at inverse temperature `beta`. -/
noncomputable def gibbs (beta : ℝ) (x : Ω) : ℝ :=
  gibbsWeight E beta x / partitionZ E beta

/-- **Thermometer theorem (zeroth-law transitivity core).** A system
with two distinct energy levels identifies the inverse temperature: two
Gibbs states that agree as distributions have equal `beta`. Contact of
`A` with `B` and of `B` with `C` therefore forces one common
temperature through any nondegenerate thermometer `B`. -/
theorem gibbs_beta_injective [Nonempty Ω] (beta1 beta2 : ℝ)
    (i j : Ω) (hij : E i ≠ E j)
    (h : gibbs E beta1 = gibbs E beta2) : beta1 = beta2 := by
  have hZ1 := partitionZ_pos E beta1
  have hZ2 := partitionZ_pos E beta2
  have hi := congrFun h i
  have hj := congrFun h j
  unfold gibbs gibbsWeight at hi hj
  have hA := (div_eq_div_iff (ne_of_gt hZ1) (ne_of_gt hZ2)).mp hi
  have hB := (div_eq_div_iff (ne_of_gt hZ1) (ne_of_gt hZ2)).mp hj
  have h3 : Real.exp (-beta1 * E i) * Real.exp (-beta2 * E j)
      * partitionZ E beta2
      = Real.exp (-beta1 * E j) * Real.exp (-beta2 * E i)
      * partitionZ E beta2 := by
    linear_combination Real.exp (-beta2 * E j) * hA
      - Real.exp (-beta2 * E i) * hB
  have key := mul_right_cancel₀ (ne_of_gt hZ2) h3
  rw [← Real.exp_add, ← Real.exp_add] at key
  have hlin := Real.exp_injective key
  have hfac : (beta1 - beta2) * (E j - E i) = 0 := by
    linear_combination hlin
  rcases mul_eq_zero.mp hfac with h' | h'
  · linarith [sub_eq_zero.mp h']
  · exact absurd (sub_eq_zero.mp h').symm hij

end ZerothLaw

section Clausius

/-- Shannon entropy with the `0 log 0 = 0` convention. -/
noncomputable def shannon (p : Ω → ℝ) : ℝ :=
  -∑ x, klTerm (p x) 1

/-- Relative entropy to a faithful reference splits into modular energy
minus entropy: `D(p ‖ τ) = ⟨K⟩_p - S(p)` with `K = -log τ`. -/
theorem kl_eq_energy_sub_shannon (p tau : Ω → ℝ)
    (hτ : ∀ x, 0 < tau x) :
    kl p tau = (∑ x, p x * (-Real.log (tau x))) - shannon p := by
  unfold kl shannon
  rw [← Finset.sum_neg_distrib, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  by_cases hp : p x = 0
  · simp [hp]
  rw [klTerm_of_pos hp, klTerm_of_pos hp]
  rw [Real.log_div hp (ne_of_gt (hτ x))]
  simp [Real.log_one]
  ring

/-- **Modular Clausius inequality.** Under any stochastic kernel that
preserves the faithful reference `τ`, the entropy change dominates the
modular-energy change: `ΔS ≥ Δ⟨K⟩` with `K = -log τ`. On the identified
thermal branch `τ = Z⁻¹ e^{-βH}` this is `ΔS ≥ β Q`. -/
theorem clausius (p tau : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x) (hτ : ∀ x, 0 < tau x)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (hstat : ∀ y, push tau K y = tau y) :
    (∑ x, push p K x * (-Real.log (tau x)))
      - (∑ x, p x * (-Real.log (tau x)))
      ≤ shannon (push p K) - shannon p := by
  have h := kl_push_le p tau K hp (fun x => le_of_lt (hτ x))
    (fun x hx => absurd hx (ne_of_gt (hτ x))) hK0 hK1
  have hstat' : push tau K = tau := funext hstat
  rw [hstat'] at h
  rw [kl_eq_energy_sub_shannon p tau hτ,
    kl_eq_energy_sub_shannon (push p K) tau hτ] at h
  linarith

end Clausius

section Unattainability

variable {π : Ω → ℝ} {b : Ω → B}

/-- Full support is preserved by the repair kernel: the diagonal entry
of every row is positive, so no atom is extinguished. Together with the
rank deficiency of the ground-sector state for `g₀ < d`, finitely many
repair steps cannot reach zero temperature. -/
theorem heatBath_preserves_pos (hπ : ∀ x, 0 < π x)
    (p : Ω → ℝ) (hp : ∀ x, 0 < p x) (y : Ω) :
    0 < push p (heatBath π b) y := by
  unfold push
  have hterm : 0 < p y * heatBath π b y y := by
    apply mul_pos (hp y)
    unfold heatBath
    rw [if_pos rfl]
    exact div_pos (hπ y) (fiberMass_pos hπ y)
  have hnonneg : ∀ x ∈ Finset.univ, 0 ≤ p x * heatBath π b x y := by
    intro x _
    apply mul_nonneg (le_of_lt (hp x))
    unfold heatBath
    by_cases h : b y = b x
    · rw [if_pos h]
      exact le_of_lt (div_pos (hπ y) (fiberMass_pos hπ x))
    · rw [if_neg h]
  exact Finset.sum_pos' hnonneg ⟨y, Finset.mem_univ y, hterm⟩

end Unattainability

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.kl_push_le
#print axioms OPH.Thermodynamics.heatBath_row_optimal
#print axioms OPH.Thermodynamics.heatBath_secondLaw
#print axioms OPH.Thermodynamics.gibbs_minimizer
#print axioms OPH.Thermodynamics.excitedMass_le
#print axioms OPH.Thermodynamics.excitedMass_lt_of_beta_large
#print axioms OPH.Thermodynamics.gibbs_beta_injective
#print axioms OPH.Thermodynamics.clausius
#print axioms OPH.Thermodynamics.heatBath_preserves_pos
