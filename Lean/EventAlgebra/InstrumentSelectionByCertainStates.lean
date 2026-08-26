import EventAlgebra.LuedersPhaseInstrument
import EventAlgebra.FrequencyConcentration
import FiniteConditionalRepair

set_option autoImplicit false

namespace EventAlgebra

/-!
# Selection of the Lüders instrument by certain-state invariance

`EventAlgebra.LuedersPhaseInstrument` proves that a declared effect table
does not determine the instrument (`effect_table_does_not_determine_instrument`):
the swap-twisted instrument carries the same induced effects as the Lüders
instrument and differs from it.  The 2026-08-26 audit records this as
finding F4, "effects do not select a Lüders instrument".  This module
proves that a condition, certain-state invariance, does select the Lüders
map among Kraus-form outcome maps with projective induced effect, in both
directions, and instantiates the selection on the committed
`PhaseInstrument` objects.

## Objects

* `HasKrausFormWithEffect Φ P`: the linear map `Φ` on `Matrix (Fin n) (Fin n) ℂ`
  is written by a finite Kraus family `K : Fin m → CMat n`,
  `Φ X = ∑ c, K c * X * (K c)ᴴ`, whose induced effect `∑ c, (K c)ᴴ * K c`
  is the matrix `P`.  This finite Kraus form is the definition of an outcome
  operation used here.  It implies complete positivity
  (`HasKrausFormWithEffect.isCompletelyPositive`, through the committed
  Kraus criterion) and the induced-effect trace clause
  `(Φ X).trace = bornWeight X P` (`HasKrausFormWithEffect.trace_apply`).
* `FixesCertainStates Φ P`: `Φ σ = σ` for every state `σ` in the committed
  certainty set `certainStates P` (states of Born weight one for `P`).
  This is **certain-state invariance** (no-disturbance on outcome-certain
  states): a state already certain of the outcome is left unchanged by the
  outcome operation.  It is a declared operational hypothesis.
* `NormalizedOutputCertain Φ P`: the component form of the committed
  `PhaseInstrument.Repeatable` (the normalized post-measurement state of
  every state of nonzero Born weight is certain of `P`);
  `PhaseInstrument.repeatable_iff_normalizedOutputCertain` records that the
  committed predicate is definitionally the conjunction of these clauses
  over the committed effect table.

## What is proved

* `kraus_mul_compl_eq_zero`, `kraus_eq_mul_event`: a Kraus family with
  induced effect the event `P` annihilates the complement `1 - P`, so each
  Kraus operator satisfies `K c = K c * P`.
* `rank_one_decomposition_collinear`: if `∑ c, v_c v_cᴴ = ψ ψᴴ` then every
  `v_c` is a scalar multiple of `ψ`.
* `apply_compress_eq_of_fixesCertainStates`: a linear map fixing every state
  certain of the event `P` fixes every compressed matrix `P * X * P`
  (rank-one certain states, then sesquilinear polarization, then the
  matrix-unit expansion of `P * X * P`).
* **Headline** `kraus_component_fixing_certain_states_eq_compress` and
  `instrument_component_fixing_certain_states_is_lueders`: a Kraus-form
  outcome map with projective induced effect `P` that fixes the states
  certain of `P` is the Lüders outcome map `X ↦ P * X * P`.
* `kraus_eq_smul_event_of_fixesCertainStates`: under the same hypotheses
  every Kraus operator is a scalar multiple `λ_c • P` of the event, with
  `∑ c, |λ_c|² = 1` when `P ≠ 0`.
* `instrument_with_projective_effects_fixing_certain_states_is_lueders`: for
  a family of events `P k` and outcome maps `Φ k`, each in Kraus form with
  induced effect `P k` and each fixing `certainStates (P k)`, every `Φ k` is
  the Lüders outcome map of `P k`, and the normalized post-measurement state
  is the committed `luedersUpdate ρ (P k)`.
* `luedersOutcomeMap_hasKrausFormWithEffect`,
  `luedersOutcomeMap_fixesCertainStates`,
  `hasKrausFormWithEffect_and_fixesCertainStates_iff_eq_luedersOutcomeMap`:
  the Lüders map itself satisfies both hypotheses, so the selection is an
  equivalence: `Φ` is a Kraus-form map with induced effect `P` fixing
  `certainStates P` if and only if `Φ = luedersOutcomeMap P`.
* **Comparison with the committed repeatability clause.**
  `normalizedOutputCertain_of_fixesCertainStates`: certain-state invariance
  (with Kraus form) implies the committed-style clause.  The converse
  fails: for a left isometry `U` (`Uᴴ * U = 1`) whose corner preserves the
  range of `P` (`P * (U * P) = U * P`), the twisted map
  `luedersOutcomeMap (U * P)` has Kraus form with induced effect `P`
  (`luedersOutcomeMap_isometry_mul_hasKrausFormWithEffect`) and satisfies
  the committed-style clause
  (`luedersOutcomeMap_isometry_mul_normalizedOutputCertain`), while
  certain-state invariance forces `U * P = λ • P`
  (`isometry_twist_eq_smul_of_fixesCertainStates`).  On `Fin 2` with
  `P = 1` the swap unitary realizes this
  (`swapOutcome_normalizedOutputCertain_one`,
  `committed_repeatability_clause_does_not_select_lueders`): the committed
  repeatability clause does not select the Lüders map; certain-state
  invariance is strictly stronger.
* Scope witnesses on `Fin 2` with `P = 1` (a rank-two projector).  The
  swap-conjugation map `luedersOutcomeMap swapUnitary` has Kraus form with
  induced effect `1`, is not the Lüders map of `1`, and fails certain-state
  invariance (`swapOutcome_hasKrausFormWithEffect_one`,
  `swapOutcome_ne_luedersOutcomeMap_one`,
  `swapOutcome_not_fixesCertainStates_one`); the packaged statement
  `effect_does_not_select_but_certain_states_do` carries this one witness
  together with the selection at `P = 1`.
* **Instantiation on the committed objects.**
  `phaseInstrument_eq_lueders_of_kraus_fixesCertainStates` and
  `phaseInstrument_eq_lueders_iff`: a `PhaseInstrument` whose outcome maps
  have Kraus form with the committed effects `committedEffectPair c i` and
  fix `certainStates (committedEffectPair c i)` is the declared
  `luedersPhaseInstrument`, and conversely.  The committed
  `swapTwistedPhaseInstrument` has Kraus form with the committed effects in
  every context (`swapTwistedPhaseInstrument_hasKrausFormWithEffect`) and
  fails certain-state invariance in the diagonal context
  (`swapTwistedPhaseInstrument_not_fixesCertainStates_diagonal`), so the
  hypothesis bundle is inhabited by exactly the Lüders instrument among the
  two committed instruments
  (`committed_instruments_certain_state_invariance_selects_lueders`).
* **Classical specialization.**
  `classical_subkernel_fixing_certain_laws_is_restriction`: a nonnegative
  finite kernel whose row sums are the indicator of a set `S` and which
  fixes every probability law of mass one on `S` is the restriction kernel
  of `S`, the classical Lüders component.


## Boundary

Finite dimension throughout.  The finite Kraus family is the definition of
an outcome operation; nothing is proved about infinite families or about
maps given only as completely positive.  Certain-state invariance is a
declared operational hypothesis (no-disturbance on outcome-certain states);
it is not derived from any source, run, or consensus rule.  It is strictly
stronger than the committed repeatability predicate
`PhaseInstrument.Repeatable` (normalized output certain of the effect):
that committed clause does not select the Lüders map (corner-isometry
twists satisfy it), so the selection proved here is a selection by
certain-state invariance, not by the committed repeatability clause.
Nothing here source-selects the outcome events `P k`, the preparation, the
readback, or producer identity; the normalization `∑ k, P k = 1` of an
outcome family is not used and not asserted.  The link between
certain-state invariance and consensus repair (a repair leaves
already-agreed configurations fixed) is an inference outside Lean: the
`FiniteConditionalRepair` kernel `heatBath` is proved here not to satisfy
the clause on point laws; the other committed finite kernels and repair
maps (`ObservableNormalForms.conditionalResamplingKernel`,
`OPH.Dynamics.recordKernel` in `Dynamics/PublicMarkov`,
`uniformKernel` in `InformationProjection/ReferenceNormalForm`, and the
quotient repair endomorphisms `publicRepair` / `primitivePublicRepair` in
`Tower/FixedPointEndpoint`) were not examined, and no committed object is
typed here as an instrument component.  Register rows PR-02 (state
surface), PR-03 (valuation additivity) and PR-64 (channel clauses of a
phase-sensitive measurement) are cited as context; none is discharged.
Issues: #730 (primary), #739.

## Tagging convention

As in `EventAlgebra.Basic`: **algebra-only** statements consume no trace
pairing; **trace-dependent** statements pass through `Tr(ρ M)`.
-/

open Matrix
open OPH.QFT
open OPH.Dynamics (CMat IsCompletelyPositive isCompletelyPositive_of_kraus)
open EventAlgebra.FrequencyConcentration (eventW eventW_isEvent)
open scoped ComplexOrder

noncomputable section

variable {n : ℕ}

/-! ## The two hypotheses -/

/-- **Kraus form with induced effect.**  `Φ` is written by a finite Kraus
family `K : Fin m → CMat n` as `Φ X = ∑ c, K c * X * (K c)ᴴ`, and the induced
effect `∑ c, (K c)ᴴ * K c` equals `P`. -/
def HasKrausFormWithEffect (Φ : CMat n →ₗ[ℂ] CMat n) (P : CMat n) : Prop :=
  ∃ (m : ℕ) (K : Fin m → CMat n),
    (∀ X, Φ X = ∑ c, K c * X * (K c)ᴴ) ∧ ∑ c, (K c)ᴴ * K c = P

/-- **Certain-state invariance** (no-disturbance on outcome-certain states).
`Φ` leaves every state certain of `P` unchanged.  A declared operational
hypothesis; strictly stronger than the committed repeatability clause
`PhaseInstrument.Repeatable` (see `NormalizedOutputCertain` below). -/
def FixesCertainStates (Φ : CMat n →ₗ[ℂ] CMat n) (P : CMat n) : Prop :=
  ∀ σ ∈ certainStates P, Φ σ = σ

/-- **Algebra-only.**  A Kraus-form map is completely positive (committed
Kraus criterion). -/
theorem HasKrausFormWithEffect.isCompletelyPositive
    {Φ : CMat n →ₗ[ℂ] CMat n} {P : CMat n}
    (h : HasKrausFormWithEffect Φ P) : IsCompletelyPositive Φ := by
  obtain ⟨m, K, hΦ, -⟩ := h
  exact isCompletelyPositive_of_kraus Φ K hΦ

/-- **Trace-dependent.**  The induced-effect clause: the trace of the
outcome map on any input is the Born weight of the induced effect. -/
theorem HasKrausFormWithEffect.trace_apply
    {Φ : CMat n →ₗ[ℂ] CMat n} {P : CMat n}
    (h : HasKrausFormWithEffect Φ P) (X : CMat n) :
    (Φ X).trace = bornWeight X P := by
  obtain ⟨m, K, hΦ, hK⟩ := h
  rw [hΦ, trace_sum, bornWeight, ← hK, Finset.mul_sum, trace_sum]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [trace_mul_cycle, trace_mul_comm]

/-! ## (1) Kraus operators with projective effect annihilate the complement -/

/-- **Algebra-only.**  If a finite Kraus family has induced effect the event
`P`, every Kraus operator annihilates the complement `1 - P`: compressing
the effect identity by `1 - P` gives a vanishing sum of positive
semidefinite matrices, each of which must vanish. -/
theorem kraus_mul_compl_eq_zero {ι : Type*} [Fintype ι] {K : ι → CMat n}
    {P : CMat n} (hP : IsEvent P) (hK : ∑ c, (K c)ᴴ * K c = P) (c : ι) :
    K c * (1 - P) = 0 := by
  have hQ : IsEvent (1 - P) := hP.compl
  have hsum : ∑ c, (K c * (1 - P))ᴴ * (K c * (1 - P)) = 0 := by
    have hexp : ∑ c, (K c * (1 - P))ᴴ * (K c * (1 - P)) =
        (1 - P)ᴴ * (∑ c, (K c)ᴴ * K c) * (1 - P) := by
      rw [Finset.mul_sum, Finset.sum_mul]
      refine Finset.sum_congr rfl fun c _ => ?_
      rw [conjTranspose_mul]
      simp only [mul_assoc]
    rw [hexp, hK, hQ.1.eq, compl_mul_self_eq_zero hP, zero_mul]
  have hnn : ∀ c ∈ (Finset.univ : Finset ι),
      0 ≤ ((K c * (1 - P))ᴴ * (K c * (1 - P))).trace :=
    fun c _ => (posSemidef_conjTranspose_mul_self _).trace_nonneg
  have htr : ∑ c, ((K c * (1 - P))ᴴ * (K c * (1 - P))).trace = 0 := by
    rw [← trace_sum, hsum, trace_zero]
  have hzero : (K c * (1 - P))ᴴ * (K c * (1 - P)) = 0 :=
    (posSemidef_conjTranspose_mul_self _).trace_eq_zero_iff.mp
      ((Finset.sum_eq_zero_iff_of_nonneg hnn).mp htr c (Finset.mem_univ c))
  exact conjTranspose_mul_self_eq_zero.mp hzero

/-- **Algebra-only.**  Each Kraus operator of a family with projective
induced effect `P` is absorbed by `P` on the right. -/
theorem kraus_eq_mul_event {ι : Type*} [Fintype ι] {K : ι → CMat n}
    {P : CMat n} (hP : IsEvent P) (hK : ∑ c, (K c)ᴴ * K c = P) (c : ι) :
    K c = K c * P := by
  have h := kraus_mul_compl_eq_zero hP hK c
  rw [mul_sub, mul_one, sub_eq_zero] at h
  exact h

/-! ## (2) Rank-one decompositions are collinear -/

/-- **Algebra-only.**  Testing a rank-one decomposition against a vector:
the quadratic forms of both sides at `φ` are sums of squared moduli. -/
theorem sum_star_dotProduct_mul_of_sum_vecMulVec {ι κ : Type*} [Fintype ι]
    [Fintype κ] {v : ι → κ → ℂ} {ψ : κ → ℂ}
    (h : ∑ c, vecMulVec (v c) (star (v c)) = vecMulVec ψ (star ψ))
    (φ : κ → ℂ) :
    ∑ c, star (star φ ⬝ᵥ v c) * (star φ ⬝ᵥ v c) =
      star (star φ ⬝ᵥ ψ) * (star φ ⬝ᵥ ψ) := by
  have hs : ∀ a : κ → ℂ, star a ⬝ᵥ φ = star (star φ ⬝ᵥ a) := fun a => by
    rw [star_dotProduct]
  have := congrArg (fun M : Matrix κ κ ℂ => star φ ⬝ᵥ M *ᵥ φ) h
  simp only [sum_mulVec, dotProduct_sum, vecMulVec_mulVec, dotProduct_smul] at this
  simpa only [hs, mul_comm] using this

/-- **Algebra-only.**  If `∑ c, v_c v_cᴴ = ψ ψᴴ` then every `v_c` is a
scalar multiple of `ψ`: the component of `v_c` orthogonal to `ψ` has
vanishing quadratic form on both sides. -/
theorem rank_one_decomposition_collinear {ι κ : Type*} [Fintype ι]
    [Fintype κ] {v : ι → κ → ℂ} {ψ : κ → ℂ}
    (h : ∑ c, vecMulVec (v c) (star (v c)) = vecMulVec ψ (star ψ)) (c : ι) :
    ∃ lam : ℂ, v c = lam • ψ := by
  -- A test vector orthogonal to `ψ` is orthogonal to every `v c`.
  have hkill : ∀ φ : κ → ℂ, star φ ⬝ᵥ ψ = 0 → star φ ⬝ᵥ v c = 0 := by
    intro φ hφ
    have hsum := sum_star_dotProduct_mul_of_sum_vecMulVec h φ
    rw [hφ, mul_zero] at hsum
    have hnn : ∀ d ∈ (Finset.univ : Finset ι),
        0 ≤ star (star φ ⬝ᵥ v d) * (star φ ⬝ᵥ v d) :=
      fun d _ => star_mul_self_nonneg _
    have hc := (Finset.sum_eq_zero_iff_of_nonneg hnn).mp hsum c (Finset.mem_univ c)
    rcases mul_eq_zero.mp hc with h0 | h0
    · exact star_eq_zero.mp h0
    · exact h0
  by_cases hψ : ψ = 0
  · refine ⟨0, ?_⟩
    have h0 := hkill (v c) (by rw [hψ, dotProduct_zero])
    rw [dotProduct_star_self_eq_zero] at h0
    rw [h0, zero_smul]
  · set lam : ℂ := (star ψ ⬝ᵥ v c) / (star ψ ⬝ᵥ ψ) with hlam
    refine ⟨lam, ?_⟩
    set φ : κ → ℂ := v c - lam • ψ with hφdef
    have hψψ : star ψ ⬝ᵥ ψ ≠ 0 := fun h0 => hψ (dotProduct_star_self_eq_zero.mp h0)
    have hψφ : star ψ ⬝ᵥ φ = 0 := by
      rw [hφdef, dotProduct_sub, dotProduct_smul, smul_eq_mul, hlam,
        div_mul_cancel₀ _ hψψ, sub_self]
    have hφψ : star φ ⬝ᵥ ψ = 0 := by
      rw [star_dotProduct, hψφ, star_zero]
    have hφv : star φ ⬝ᵥ v c = 0 := hkill φ hφψ
    have hφφ : star φ ⬝ᵥ φ = 0 := by
      rw [hφdef, dotProduct_sub, dotProduct_smul, hφv, ← hφdef, hφψ, smul_zero,
        sub_zero]
    have hφ0 : φ = 0 := dotProduct_star_self_eq_zero.mp hφφ
    rw [hφdef] at hφ0
    exact sub_eq_zero.mp hφ0

/-! ## Certain states span the compressed algebra -/

/-- **Trace-dependent.**  The normalized rank-one matrix of a nonzero vector
in the range of the event `P` is a state certain of `P`. -/
theorem smul_vecMulVec_mem_certainStates {P : CMat n} (hP : IsEvent P)
    {ψ : Fin n → ℂ} (hψ : P *ᵥ ψ = ψ) (hne : ψ ≠ 0) :
    (star ψ ⬝ᵥ ψ)⁻¹ • vecMulVec ψ (star ψ) ∈ certainStates P := by
  have hψψ : star ψ ⬝ᵥ ψ ≠ 0 := fun h0 => hne (dotProduct_star_self_eq_zero.mp h0)
  have hpos : 0 < star ψ ⬝ᵥ ψ :=
    lt_of_le_of_ne (dotProduct_star_self_nonneg ψ) (Ne.symm hψψ)
  have hinv : 0 ≤ (star ψ ⬝ᵥ ψ)⁻¹ := (RCLike.inv_pos.mpr hpos).le
  have habsorb : star ψ ᵥ* P = star ψ := by
    have h1 : star ψ ᵥ* Pᴴ = star ψ := by rw [← star_mulVec, hψ]
    rwa [hP.1.eq] at h1
  refine ⟨⟨(posSemidef_vecMulVec_self_star ψ).smul hinv, ?_⟩, ?_⟩
  · rw [trace_smul, trace_vecMulVec, smul_eq_mul, dotProduct_comm ψ (star ψ),
      inv_mul_cancel₀ hψψ]
  · rw [bornWeight, smul_mul_assoc, vecMulVec_mul, habsorb, trace_smul,
      trace_vecMulVec, smul_eq_mul, dotProduct_comm ψ (star ψ), inv_mul_cancel₀ hψψ]

/-- **Trace-dependent.**  A linear map fixing every state certain of `P`
fixes the rank-one matrix `ψ ψᴴ` of every vector in the range of `P`. -/
theorem apply_vecMulVec_self_of_fixesCertainStates {Φ : CMat n →ₗ[ℂ] CMat n}
    {P : CMat n} (hP : IsEvent P) (hfix : FixesCertainStates Φ P)
    {ψ : Fin n → ℂ} (hψ : P *ᵥ ψ = ψ) :
    Φ (vecMulVec ψ (star ψ)) = vecMulVec ψ (star ψ) := by
  by_cases hne : ψ = 0
  · rw [hne, star_zero, vecMulVec_zero, map_zero]
  · have h := hfix _ (smul_vecMulVec_mem_certainStates hP hψ hne)
    rw [map_smul] at h
    have hc : (star ψ ⬝ᵥ ψ)⁻¹ ≠ 0 :=
      inv_ne_zero fun h0 => hne (dotProduct_star_self_eq_zero.mp h0)
    exact smul_right_injective (CMat n) hc h

/-- **Trace-dependent.**  Sesquilinear polarization: a linear map fixing
every state certain of `P` fixes `ψ φᴴ` for all `ψ, φ` in the range of `P`. -/
theorem apply_vecMulVec_of_fixesCertainStates {Φ : CMat n →ₗ[ℂ] CMat n}
    {P : CMat n} (hP : IsEvent P) (hfix : FixesCertainStates Φ P)
    {ψ φ : Fin n → ℂ} (hψ : P *ᵥ ψ = ψ) (hφ : P *ᵥ φ = φ) :
    Φ (vecMulVec ψ (star φ)) = vecMulVec ψ (star φ) := by
  have hdiag : ∀ x : Fin n → ℂ, P *ᵥ x = x →
      Φ (vecMulVec x (star x)) = vecMulVec x (star x) :=
    fun x hx => apply_vecMulVec_self_of_fixesCertainStates hP hfix hx
  -- The two polarization directions `ψ + φ` and `ψ + i φ`.
  have hadd : P *ᵥ (ψ + φ) = ψ + φ := by rw [mulVec_add, hψ, hφ]
  have hI : P *ᵥ (ψ + Complex.I • φ) = ψ + Complex.I • φ := by
    rw [mulVec_add, mulVec_smul, hψ, hφ]
  have e1 := hdiag _ hadd
  have e2 := hdiag _ hI
  have eψ := hdiag ψ hψ
  have eφ := hdiag φ hφ
  set a := Φ (vecMulVec ψ (star φ)) with ha
  set b := vecMulVec ψ (star φ) with hb
  set a' := Φ (vecMulVec φ (star ψ)) with ha'
  set b' := vecMulVec φ (star ψ) with hb'
  have hstarI : star Complex.I = -Complex.I := by
    rw [Complex.star_def, Complex.conj_I]
  -- First direction: `a + a' = b + b'`.
  have h1 : a + a' = b + b' := by
    rw [star_add, vecMulVec_add, add_vecMulVec, add_vecMulVec, map_add, map_add,
      map_add, eψ, eφ] at e1
    rw [← ha, ← hb, ← ha', ← hb'] at e1
    linear_combination (norm := abel) e1
  -- Second direction: `(-i) a + i a' = (-i) b + i b'`.
  have h2 : (-Complex.I) • a + Complex.I • a' = (-Complex.I) • b + Complex.I • b' := by
    rw [star_add, star_smul, hstarI, vecMulVec_add, add_vecMulVec, add_vecMulVec,
      vecMulVec_smul, smul_vecMulVec, smul_vecMulVec, vecMulVec_smul, smul_smul,
      mul_neg, Complex.I_mul_I, neg_neg, one_smul, map_add, map_add, map_add,
      map_smul, map_smul, eψ, eφ] at e2
    rw [← ha, ← hb, ← ha', ← hb'] at e2
    linear_combination (norm := abel) e2
  -- Multiply the second direction by `i`: `a - a' = b - b'`.
  have h3 : a - a' = b - b' := by
    have := congrArg (fun M : CMat n => Complex.I • M) h2
    simp only [smul_add, smul_smul, mul_neg, Complex.I_mul_I, neg_neg, one_smul,
      neg_one_smul] at this
    linear_combination (norm := abel) this
  -- Combine: `2 a = 2 b`.
  have h4 : (2 : ℂ) • a = (2 : ℂ) • b := by
    rw [two_smul, two_smul]
    linear_combination (norm := abel) h1 + h3
  exact smul_right_injective (CMat n) (by norm_num : (2 : ℂ) ≠ 0) h4

/-- **Algebra-only.**  Column `i` of an event lies in its own range. -/
theorem event_mulVec_transpose {P : CMat n} (hP : IsEvent P) (i : Fin n) :
    P *ᵥ Pᵀ i = Pᵀ i := by
  ext a
  have h := congrFun (congrFun hP.2 a) i
  rw [mul_apply] at h
  simpa [mulVec, dotProduct, transpose_apply] using h

/-- **Algebra-only.**  Matrix-unit expansion of the compression `P * X * P`
by rank-one matrices of columns of the Hermitian idempotent `P`. -/
theorem compress_eq_sum_vecMulVec {P : CMat n} (hP : IsEvent P) (X : CMat n) :
    P * X * P = ∑ i, ∑ j, X i j • vecMulVec (Pᵀ i) (star (Pᵀ j)) := by
  ext a b
  have hconj : ∀ j, star (P b j) = P j b := fun j => hP.1.apply j b
  simp only [Matrix.sum_apply, Matrix.smul_apply, vecMulVec_apply, transpose_apply,
    smul_eq_mul, mul_apply, Pi.star_apply, hconj, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring

/-- **Trace-dependent.**  A linear map fixing every state certain of the
event `P` fixes every compressed matrix `P * X * P`. -/
theorem apply_compress_eq_of_fixesCertainStates {Φ : CMat n →ₗ[ℂ] CMat n}
    {P : CMat n} (hP : IsEvent P) (hfix : FixesCertainStates Φ P) (X : CMat n) :
    Φ (P * X * P) = P * X * P := by
  rw [compress_eq_sum_vecMulVec hP X, map_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [map_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [map_smul, apply_vecMulVec_of_fixesCertainStates hP hfix
    (event_mulVec_transpose hP i) (event_mulVec_transpose hP j)]

/-! ## (4) Headline: certain-state invariance selects the Lüders map -/

/-- **Trace-dependent.  Headline (component form).**  A Kraus-form outcome
map with projective induced effect `P` that fixes every state certain of
`P` is the compression `X ↦ P * X * P` on every matrix.  The proof absorbs
each Kraus operator into `P` (`kraus_eq_mul_event`), so `Φ X = Φ (P X P)`,
and then applies the fixed-point extension
`apply_compress_eq_of_fixesCertainStates`. -/
theorem kraus_component_fixing_certain_states_eq_compress {ι : Type*} [Fintype ι]
    {Φ : CMat n →ₗ[ℂ] CMat n} {K : ι → CMat n} {P : CMat n} (hP : IsEvent P)
    (hΦ : ∀ X, Φ X = ∑ c, K c * X * (K c)ᴴ) (hK : ∑ c, (K c)ᴴ * K c = P)
    (hfix : FixesCertainStates Φ P) (X : CMat n) :
    Φ X = P * X * P := by
  have hKP : ∀ c, K c = K c * P := kraus_eq_mul_event hP hK
  have h1 : Φ X = Φ (P * X * P) := by
    rw [hΦ, hΦ]
    refine Finset.sum_congr rfl fun c _ => ?_
    conv_lhs => rw [hKP c]
    rw [conjTranspose_mul, hP.1.eq]
    simp only [mul_assoc]
  rw [h1, apply_compress_eq_of_fixesCertainStates hP hfix]

/-- **Trace-dependent.  Headline.**  A Kraus-form outcome map with
projective induced effect `P` fixing the states certain of `P` is the
Lüders outcome map of `P`. -/
theorem instrument_component_fixing_certain_states_is_lueders
    {Φ : CMat n →ₗ[ℂ] CMat n} {P : CMat n} (hP : IsEvent P)
    (hK : HasKrausFormWithEffect Φ P) (hfix : FixesCertainStates Φ P) :
    Φ = luedersOutcomeMap P := by
  obtain ⟨m, K, hΦ, hKP⟩ := hK
  refine LinearMap.ext fun X => ?_
  rw [luedersOutcomeMap_apply_of_isHermitian hP.1]
  exact kraus_component_fixing_certain_states_eq_compress hP hΦ hKP hfix X

/-! ## (3) Each Kraus operator is a scalar multiple of the event -/

/-- **Algebra-only.**  If `∑ c, K c * X * (K c)ᴴ = P * X * P` for every `X`
and `P` is Hermitian, every `K c` is a scalar multiple of `P`: evaluating on
matrix units flattens the identity to a rank-one decomposition
`∑ c, k_c k_cᴴ = p pᴴ` on `Fin n × Fin n`. -/
theorem kraus_eq_smul_event_of_compress {ι : Type*} [Fintype ι]
    {K : ι → CMat n} {P : CMat n} (hP : P.IsHermitian)
    (h : ∀ X, ∑ c, K c * X * (K c)ᴴ = P * X * P) (c : ι) :
    ∃ lam : ℂ, K c = lam • P := by
  have hterm : ∀ (M : CMat n) (a i b j : Fin n),
      (M * Matrix.single i j (1 : ℂ) * Mᴴ) a b = M a i * star (M b j) := by
    intro M a i b j
    rw [mul_apply]
    simp only [conjTranspose_apply]
    rw [Finset.sum_eq_single j]
    · rw [mul_single_apply_same, mul_one]
    · intro x _ hx
      rw [mul_single_apply_of_ne _ _ _ _ _ hx, zero_mul]
    · intro hj
      exact absurd (Finset.mem_univ j) hj
  have hflat : ∑ d, vecMulVec (fun p : Fin n × Fin n => K d p.1 p.2)
      (star (fun p : Fin n × Fin n => K d p.1 p.2)) =
      vecMulVec (fun p : Fin n × Fin n => P p.1 p.2)
        (star (fun p : Fin n × Fin n => P p.1 p.2)) := by
    ext ⟨a, i⟩ ⟨b, j⟩
    have hij := congrFun (congrFun (h (Matrix.single i j (1 : ℂ))) a) b
    have hR := hterm P a i b j
    rw [hP.eq] at hR
    rw [Matrix.sum_apply, hR] at hij
    simp only [hterm] at hij
    rw [Matrix.sum_apply]
    simp only [vecMulVec_apply, Pi.star_apply]
    exact hij
  obtain ⟨lam, hlam⟩ := rank_one_decomposition_collinear hflat c
  refine ⟨lam, ?_⟩
  ext a b
  have := congrFun hlam (a, b)
  simpa using this

/-- **Trace-dependent.**  Under the headline hypotheses every Kraus
operator is `λ_c • P`, and when `P ≠ 0` the scalars satisfy
`∑ c, |λ_c|² = 1` (read as `∑ c, star λ_c * λ_c = 1`). -/
theorem kraus_eq_smul_event_of_fixesCertainStates {ι : Type*} [Fintype ι]
    {Φ : CMat n →ₗ[ℂ] CMat n} {K : ι → CMat n} {P : CMat n} (hP : IsEvent P)
    (hΦ : ∀ X, Φ X = ∑ c, K c * X * (K c)ᴴ) (hK : ∑ c, (K c)ᴴ * K c = P)
    (hfix : FixesCertainStates Φ P) :
    ∃ lam : ι → ℂ, (∀ c, K c = lam c • P) ∧
      (P ≠ 0 → ∑ c, star (lam c) * lam c = 1) := by
  have hcomp : ∀ X, ∑ c, K c * X * (K c)ᴴ = P * X * P := fun X => by
    rw [← hΦ X]
    exact kraus_component_fixing_certain_states_eq_compress hP hΦ hK hfix X
  choose lam hlam using kraus_eq_smul_event_of_compress hP.1 hcomp
  refine ⟨lam, hlam, fun hP0 => ?_⟩
  have hsum : ∑ c, (K c)ᴴ * K c = (∑ c, star (lam c) * lam c) • P := by
    rw [Finset.sum_smul]
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [hlam c, conjTranspose_smul, smul_mul_smul_comm, hP.1.eq, hP.2]
  rw [hK] at hsum
  have hzero : ((∑ c, star (lam c) * lam c) - 1) • P = 0 := by
    rw [sub_smul, one_smul, ← hsum, sub_self]
  rcases smul_eq_zero.mp hzero with h0 | h0
  · exact sub_eq_zero.mp h0
  · exact absurd h0 hP0

/-! ## The outcome-family form -/

/-- **Trace-dependent.**  For a family of events `P k` with outcome maps
`Φ k`, each in Kraus form with induced effect `P k` and each fixing the
states certain of `P k`, every `Φ k` is the Lüders outcome map of `P k`, and
the normalized post-measurement state on any `ρ` is the committed Lüders
update `luedersUpdate ρ (P k)`.  No normalization `∑ k, P k = 1` of the
family is used. -/
theorem instrument_with_projective_effects_fixing_certain_states_is_lueders
    {κ : Type*} (P : κ → CMat n) (Φ : κ → CMat n →ₗ[ℂ] CMat n)
    (hP : ∀ k, IsEvent (P k)) (hK : ∀ k, HasKrausFormWithEffect (Φ k) (P k))
    (hfix : ∀ k, FixesCertainStates (Φ k) (P k)) :
    (∀ k, Φ k = luedersOutcomeMap (P k)) ∧
      ∀ (k : κ) (ρ : CMat n),
        (bornWeight ρ (P k))⁻¹ • Φ k ρ = luedersUpdate ρ (P k) := by
  have hL : ∀ k, Φ k = luedersOutcomeMap (P k) := fun k =>
    instrument_component_fixing_certain_states_is_lueders (hP k) (hK k) (hfix k)
  refine ⟨hL, fun k ρ => ?_⟩
  rw [hL k, luedersOutcomeMap_normalized (hP k)]

/-! ## (6) The Lüders map satisfies both hypotheses; the selection is an iff -/

/-- **Algebra-only.**  The Lüders outcome map of an event has Kraus form with
the singleton family `{P}` and induced effect `P`. -/
theorem luedersOutcomeMap_hasKrausFormWithEffect {P : CMat n} (hP : IsEvent P) :
    HasKrausFormWithEffect (luedersOutcomeMap P) P :=
  ⟨1, fun _ => P, fun X => by rw [Fin.sum_univ_one, luedersOutcomeMap_apply],
    by rw [Fin.sum_univ_one, hP.1.eq, hP.2]⟩

/-- **Trace-dependent.**  The Lüders outcome map fixes every state certain
of its event (committed absorption `mul_eq_self_of_bornWeight_one`). -/
theorem luedersOutcomeMap_fixesCertainStates {P : CMat n} (hP : IsEvent P) :
    FixesCertainStates (luedersOutcomeMap P) P := by
  intro σ hσ
  obtain ⟨hs, h1⟩ := hσ
  obtain ⟨hr, hl⟩ := mul_eq_self_of_bornWeight_one hs hP h1
  rw [luedersOutcomeMap_apply_of_isHermitian hP.1, hl, hr]

/-- **Trace-dependent.  Selection as an equivalence.**  For an event `P`, a
linear map `Φ` has Kraus form with induced effect `P` and fixes the states
certain of `P` if and only if `Φ` is the Lüders outcome map of `P`.  The
forward direction is the headline; the backward direction is the pair of
lemmas above. -/
theorem hasKrausFormWithEffect_and_fixesCertainStates_iff_eq_luedersOutcomeMap
    {P : CMat n} (hP : IsEvent P) (Φ : CMat n →ₗ[ℂ] CMat n) :
    HasKrausFormWithEffect Φ P ∧ FixesCertainStates Φ P ↔
      Φ = luedersOutcomeMap P := by
  constructor
  · rintro ⟨hK, hfix⟩
    exact instrument_component_fixing_certain_states_is_lueders hP hK hfix
  · rintro rfl
    exact ⟨luedersOutcomeMap_hasKrausFormWithEffect hP,
      luedersOutcomeMap_fixesCertainStates hP⟩

/-! ## The committed repeatability clause is strictly weaker -/

/-- **Normalized-output certainty**, the component form of the committed
`PhaseInstrument.Repeatable`: for every state of nonzero Born weight on
`P`, the normalized output of `Φ` is a state certain of `P`. -/
def NormalizedOutputCertain (Φ : CMat n →ₗ[ℂ] CMat n) (P : CMat n) : Prop :=
  ∀ ρ : CMat n, IsState ρ → bornWeight ρ P ≠ 0 →
    (bornWeight ρ P)⁻¹ • Φ ρ ∈ certainStates P

/-- The committed repeatability predicate is definitionally the conjunction
of the component clauses over the committed effect table. -/
theorem PhaseInstrument.repeatable_iff_normalizedOutputCertain (Φ : PhaseInstrument) :
    Φ.Repeatable ↔ ∀ (c : InstrumentContext) (i : Fin 2),
      NormalizedOutputCertain (Φ.outcomeMap c i) (committedEffectPair c i) :=
  Iff.rfl

/-- **Trace-dependent.**  Certain-state invariance with Kraus form implies
the committed-style clause (through the headline and the committed
`luedersOutcomeMap_normalized_mem_certainStates`). -/
theorem normalizedOutputCertain_of_fixesCertainStates
    {Φ : CMat n →ₗ[ℂ] CMat n} {P : CMat n} (hP : IsEvent P)
    (hK : HasKrausFormWithEffect Φ P) (hfix : FixesCertainStates Φ P) :
    NormalizedOutputCertain Φ P := by
  intro ρ hρ hw
  rw [instrument_component_fixing_certain_states_is_lueders hP hK hfix]
  exact luedersOutcomeMap_normalized_mem_certainStates hρ hP hw

/-- **Algebra-only.**  Twisting the Lüders map of an event by a left
isometry (`Uᴴ * U = 1`) gives a Kraus-form map with the same induced
effect `P` (singleton family `{U * P}`). -/
theorem luedersOutcomeMap_isometry_mul_hasKrausFormWithEffect {U P : CMat n}
    (hU : Uᴴ * U = 1) (hP : IsEvent P) :
    HasKrausFormWithEffect (luedersOutcomeMap (U * P)) P :=
  ⟨1, fun _ => U * P, fun X => by rw [Fin.sum_univ_one, luedersOutcomeMap_apply],
    by rw [Fin.sum_univ_one, conjTranspose_mul, hP.1.eq, mul_assoc, ← mul_assoc Uᴴ, hU,
      one_mul, hP.2]⟩

/-- **Trace-dependent.**  A corner-isometry twist (`Uᴴ * U = 1`, and the
corner `U * P` maps the range of `P` into itself, `P * (U * P) = U * P`)
satisfies the committed-style repeatability clause: every normalized
output is a state certain of `P`. -/
theorem luedersOutcomeMap_isometry_mul_normalizedOutputCertain {U P : CMat n}
    (hU : Uᴴ * U = 1) (hP : IsEvent P) (hUP : P * (U * P) = U * P) :
    NormalizedOutputCertain (luedersOutcomeMap (U * P)) P := by
  intro ρ hρ hw
  have htr : (luedersOutcomeMap (U * P) ρ).trace = bornWeight ρ P :=
    trace_luedersOutcomeMap_mul_of_isometry hU hP ρ
  have hpos : 0 < bornWeight ρ P :=
    lt_of_le_of_ne (bornWeight_nonneg hρ.1 hP) (Ne.symm hw)
  have hinv : 0 ≤ (bornWeight ρ P)⁻¹ := (RCLike.inv_pos.mpr hpos).le
  have hpsd : (luedersOutcomeMap (U * P) ρ).PosSemidef :=
    (luedersOutcomeMap_completelyPositive _).isPositiveMap ρ hρ.1
  -- The output is absorbed by `P` on the right.
  have hconj : P * Uᴴ * P = P * Uᴴ := by
    have h := congrArg conjTranspose hUP
    rw [conjTranspose_mul, conjTranspose_mul, hP.1.eq] at h
    exact h
  have habs : luedersOutcomeMap (U * P) ρ * P = luedersOutcomeMap (U * P) ρ := by
    rw [luedersOutcomeMap_apply, conjTranspose_mul, hP.1.eq]
    calc U * P * ρ * (P * Uᴴ) * P = U * P * ρ * (P * Uᴴ * P) := by
          simp only [mul_assoc]
      _ = U * P * ρ * (P * Uᴴ) := by rw [hconj]
  have hb : bornWeight (luedersOutcomeMap (U * P) ρ) P = bornWeight ρ P := by
    rw [bornWeight, habs, htr]
  refine ⟨⟨hpsd.smul hinv, ?_⟩, ?_⟩
  · rw [trace_smul, htr, smul_eq_mul, inv_mul_cancel₀ hw]
  · rw [bornWeight_smul, hb, inv_mul_cancel₀ hw]

/-- **Trace-dependent.**  If an isometry twist of the Lüders map of `P`
satisfies certain-state invariance, its corner is a scalar multiple of
`P`: the twist is trivial up to a phase. -/
theorem isometry_twist_eq_smul_of_fixesCertainStates {U P : CMat n}
    (hU : Uᴴ * U = 1) (hP : IsEvent P)
    (hfix : FixesCertainStates (luedersOutcomeMap (U * P)) P) :
    ∃ lam : ℂ, U * P = lam • P := by
  obtain ⟨lam, hlam, -⟩ := kraus_eq_smul_event_of_fixesCertainStates
    (Φ := luedersOutcomeMap (U * P)) (K := fun _ : Fin 1 => U * P) hP
    (fun X => by rw [Fin.sum_univ_one, luedersOutcomeMap_apply])
    (by rw [Fin.sum_univ_one, conjTranspose_mul, hP.1.eq, mul_assoc, ← mul_assoc Uᴴ, hU,
      one_mul, hP.2])
    hfix
  exact ⟨lam 0, hlam 0⟩

/-! ## (5) Scope witnesses on `Fin 2` with the rank-two effect `P = 1`

The witnesses reuse the committed `swapUnitary`
(`EventAlgebra.LuedersPhaseInstrument`), the committed witness event
`FrequencyConcentration.eventW = !![1, 0; 0, 0]` with its complement
`1 - eventW`, the committed `sigmaXPlusProj`
(`EventAlgebra.OperationalPhaseInstrument`). -/

/-- **Trace-dependent.**  The committed witness event `eventW` is a state
certain of the sure event `1`. -/
theorem eventW_mem_certainStates_one : eventW ∈ certainStates (1 : CMat 2) := by
  have hs : IsState eventW := ⟨eventW_isEvent.posSemidef, by
    simp [eventW, trace, Fin.sum_univ_two]⟩
  exact ⟨hs, bornWeight_one hs⟩

/-- **Trace-dependent.**  The committed `+X` projector is a state certain
of `1`. -/
theorem sigmaXPlusProj_mem_certainStates_one :
    sigmaXPlusProj ∈ certainStates (1 : CMat 2) := by
  have hs : IsState sigmaXPlusProj := ⟨sigmaXPlusProj_isEvent.posSemidef, by
    simp [sigmaXPlusProj, trace, Fin.sum_univ_two]; norm_num⟩
  exact ⟨hs, bornWeight_one hs⟩

/-- **Algebra-only.**  The swap-conjugation outcome map (Kraus operator
the committed `swapUnitary`) has Kraus form with induced effect `1`. -/
theorem swapOutcome_hasKrausFormWithEffect_one :
    HasKrausFormWithEffect (luedersOutcomeMap swapUnitary) (1 : CMat 2) :=
  ⟨1, fun _ => swapUnitary,
    fun X => by rw [Fin.sum_univ_one, luedersOutcomeMap_apply],
    by rw [Fin.sum_univ_one, swapUnitary_conjTranspose_mul]⟩

/-- **Algebra-only.**  The swap-conjugation map is not the Lüders map of `1`:
it moves `eventW = E₀₀` to `E₁₁`. -/
theorem swapOutcome_ne_luedersOutcomeMap_one :
    luedersOutcomeMap swapUnitary ≠ luedersOutcomeMap (1 : CMat 2) := by
  intro h
  have h00 := congrFun (congrFun (congrArg (fun Φ => Φ eventW) h) 0) 0
  simp [luedersOutcomeMap_apply, swapUnitary, eventW, mul_apply,
    Fin.sum_univ_two, conjTranspose_apply] at h00

/-- **Trace-dependent.**  The swap-conjugation map fails certain-state
invariance for `1`: it does not fix the state `eventW`. -/
theorem swapOutcome_not_fixesCertainStates_one :
    ¬ FixesCertainStates (luedersOutcomeMap swapUnitary) (1 : CMat 2) := by
  intro h
  have h00 := congrFun (congrFun (h _ eventW_mem_certainStates_one) 0) 0
  simp [luedersOutcomeMap_apply, swapUnitary, eventW, mul_apply,
    Fin.sum_univ_two, conjTranspose_apply] at h00

/-- **Trace-dependent.**  The swap-conjugation map satisfies the
committed-style repeatability clause for `1` (the corner-isometry lemma
with `U = swapUnitary`, `P = 1`). -/
theorem swapOutcome_normalizedOutputCertain_one :
    NormalizedOutputCertain (luedersOutcomeMap swapUnitary) (1 : CMat 2) := by
  have h := luedersOutcomeMap_isometry_mul_normalizedOutputCertain
    swapUnitary_conjTranspose_mul isEvent_one (by rw [one_mul])
  rwa [mul_one] at h

/-- **The committed repeatability clause does not select the Lüders map.**
On `Fin 2` with effect `1`, the swap-conjugation map has Kraus form with
induced effect `1`, satisfies the committed-style repeatability clause,
differs from the Lüders map, and fails certain-state invariance. -/
theorem committed_repeatability_clause_does_not_select_lueders :
    ∃ Φ : CMat 2 →ₗ[ℂ] CMat 2, HasKrausFormWithEffect Φ 1 ∧
      NormalizedOutputCertain Φ 1 ∧ Φ ≠ luedersOutcomeMap 1 ∧
      ¬ FixesCertainStates Φ 1 :=
  ⟨luedersOutcomeMap swapUnitary, swapOutcome_hasKrausFormWithEffect_one,
    swapOutcome_normalizedOutputCertain_one, swapOutcome_ne_luedersOutcomeMap_one,
    swapOutcome_not_fixesCertainStates_one⟩

/-- **Algebra-only.**  The committed summed Lüders channel of any event has
Kraus form with the two-element family `P, 1 - P` and induced effect `1`. -/
theorem luedersChannel_hasKrausFormWithEffect_one {P : CMat n} (hP : IsEvent P) :
    HasKrausFormWithEffect (luedersChannel P) 1 := by
  refine ⟨2, ![P, 1 - P], fun X => ?_, ?_⟩
  · rw [Fin.sum_univ_two, luedersChannel_apply]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
  · rw [Fin.sum_univ_two]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
    rw [hP.1.eq, hP.2, hP.compl.1.eq, hP.compl.2, add_sub_cancel]

/-- **Non-vacuity of the scope (swap witness).**  A Kraus-form map with
induced effect `1` on `Fin 2` exists besides the Lüders map and fails
certain-state invariance; every Kraus-form map with induced effect `1`
that satisfies certain-state invariance is the Lüders map.  The selection
theorem is not vacuous and the invariance hypothesis is the selecting
clause within the Kraus-form class. -/
theorem effect_does_not_select_but_certain_states_do :
    (∃ Φ : CMat 2 →ₗ[ℂ] CMat 2, HasKrausFormWithEffect Φ 1 ∧
      Φ ≠ luedersOutcomeMap 1 ∧ ¬ FixesCertainStates Φ 1) ∧
    (∀ Φ : CMat 2 →ₗ[ℂ] CMat 2, HasKrausFormWithEffect Φ 1 →
      FixesCertainStates Φ 1 → Φ = luedersOutcomeMap 1) :=
  ⟨⟨luedersOutcomeMap swapUnitary, swapOutcome_hasKrausFormWithEffect_one,
      swapOutcome_ne_luedersOutcomeMap_one, swapOutcome_not_fixesCertainStates_one⟩,
    fun _ hK hfix =>
      instrument_component_fixing_certain_states_is_lueders isEvent_one hK hfix⟩

/-! ## Instantiation on the committed `PhaseInstrument` objects -/

/-- A phase instrument is determined by its outcome maps (the remaining
fields are propositions). -/
theorem PhaseInstrument.eq_of_outcomeMap_eq {Φ Ψ : PhaseInstrument}
    (h : Φ.outcomeMap = Ψ.outcomeMap) : Φ = Ψ := by
  cases Φ
  cases Ψ
  dsimp only at h
  subst h
  rfl

/-- **Trace-dependent.**  A `PhaseInstrument` whose outcome maps have Kraus
form with the committed effects and satisfy certain-state invariance on
the committed effects is the declared Lüders phase instrument. -/
theorem phaseInstrument_eq_lueders_of_kraus_fixesCertainStates (Φ : PhaseInstrument)
    (hK : ∀ (c : InstrumentContext) (i : Fin 2),
      HasKrausFormWithEffect (Φ.outcomeMap c i) (committedEffectPair c i))
    (hfix : ∀ (c : InstrumentContext) (i : Fin 2),
      FixesCertainStates (Φ.outcomeMap c i) (committedEffectPair c i)) :
    Φ = luedersPhaseInstrument := by
  apply PhaseInstrument.eq_of_outcomeMap_eq
  funext c i
  rw [luedersPhaseInstrument_outcomeMap]
  exact instrument_component_fixing_certain_states_is_lueders
    (committedEffectPair_isEvent c i) (hK c i) (hfix c i)

/-- **Algebra-only.**  The Lüders phase instrument has Kraus form with the
committed effects in every context. -/
theorem luedersPhaseInstrument_hasKrausFormWithEffect (c : InstrumentContext)
    (i : Fin 2) :
    HasKrausFormWithEffect (luedersPhaseInstrument.outcomeMap c i)
      (committedEffectPair c i) := by
  rw [luedersPhaseInstrument_outcomeMap]
  exact luedersOutcomeMap_hasKrausFormWithEffect (committedEffectPair_isEvent c i)

/-- **Trace-dependent.**  The Lüders phase instrument satisfies
certain-state invariance on the committed effects in every context. -/
theorem luedersPhaseInstrument_fixesCertainStates (c : InstrumentContext) (i : Fin 2) :
    FixesCertainStates (luedersPhaseInstrument.outcomeMap c i)
      (committedEffectPair c i) := by
  rw [luedersPhaseInstrument_outcomeMap]
  exact luedersOutcomeMap_fixesCertainStates (committedEffectPair_isEvent c i)

/-- **Trace-dependent.  Selection on the committed objects, as an
equivalence.**  A phase instrument has Kraus-form outcome maps with the
committed effects satisfying certain-state invariance if and only if it is
the declared Lüders phase instrument. -/
theorem phaseInstrument_eq_lueders_iff (Φ : PhaseInstrument) :
    (∀ (c : InstrumentContext) (i : Fin 2),
      HasKrausFormWithEffect (Φ.outcomeMap c i) (committedEffectPair c i) ∧
        FixesCertainStates (Φ.outcomeMap c i) (committedEffectPair c i)) ↔
      Φ = luedersPhaseInstrument := by
  constructor
  · intro h
    exact phaseInstrument_eq_lueders_of_kraus_fixesCertainStates Φ
      (fun c i => (h c i).1) (fun c i => (h c i).2)
  · rintro rfl c i
    exact ⟨luedersPhaseInstrument_hasKrausFormWithEffect c i,
      luedersPhaseInstrument_fixesCertainStates c i⟩

/-- **Algebra-only.**  The committed swap-twisted phase instrument has Kraus
form with the committed effects in every context (singleton family
`{swapUnitary * committedEffectPair c i}`). -/
theorem swapTwistedPhaseInstrument_hasKrausFormWithEffect (c : InstrumentContext)
    (i : Fin 2) :
    HasKrausFormWithEffect (swapTwistedPhaseInstrument.outcomeMap c i)
      (committedEffectPair c i) := by
  rw [swapTwistedPhaseInstrument_outcomeMap]
  exact luedersOutcomeMap_isometry_mul_hasKrausFormWithEffect
    swapUnitary_conjTranspose_mul (committedEffectPair_isEvent c i)

/-- **Trace-dependent.**  The swap-twisted phase instrument fails
certain-state invariance in the diagonal web context, outcome `0`: the
corner `swapUnitary * E₀₀` is not a scalar multiple of `E₀₀`. -/
theorem swapTwistedPhaseInstrument_not_fixesCertainStates_diagonal :
    ¬ FixesCertainStates
      (swapTwistedPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0)
      (committedEffectPair (InstrumentContext.web WebContext.diagonal) 0) := by
  intro h
  rw [swapTwistedPhaseInstrument_outcomeMap] at h
  obtain ⟨lam, hlam⟩ := isometry_twist_eq_smul_of_fixesCertainStates
    swapUnitary_conjTranspose_mul (committedEffectPair_isEvent _ _) h
  have h10 := congrFun (congrFun hlam 1) 0
  rw [committedEffectPair_zero, committedContextEffect_diagonal_eq] at h10
  simp [swapUnitary, mul_apply, Fin.sum_univ_two] at h10

/-- **Selection among the committed instruments.**  Every phase instrument
with Kraus-form outcome maps on the committed effects satisfying
certain-state invariance is the Lüders phase instrument; the committed
swap-twisted instrument has Kraus form with the committed effects and
fails certain-state invariance, so the hypothesis bundle separates the two
committed instruments of `effect_table_does_not_determine_instrument`. -/
theorem committed_instruments_certain_state_invariance_selects_lueders :
    (∀ Φ : PhaseInstrument,
      (∀ (c : InstrumentContext) (i : Fin 2),
        HasKrausFormWithEffect (Φ.outcomeMap c i) (committedEffectPair c i)) →
      (∀ (c : InstrumentContext) (i : Fin 2),
        FixesCertainStates (Φ.outcomeMap c i) (committedEffectPair c i)) →
      Φ = luedersPhaseInstrument) ∧
    (∀ (c : InstrumentContext) (i : Fin 2),
      HasKrausFormWithEffect (swapTwistedPhaseInstrument.outcomeMap c i)
        (committedEffectPair c i)) ∧
    ¬ (∀ (c : InstrumentContext) (i : Fin 2),
      FixesCertainStates (swapTwistedPhaseInstrument.outcomeMap c i)
        (committedEffectPair c i)) :=
  ⟨phaseInstrument_eq_lueders_of_kraus_fixesCertainStates,
    swapTwistedPhaseInstrument_hasKrausFormWithEffect,
    fun h => swapTwistedPhaseInstrument_not_fixesCertainStates_diagonal
      (h (InstrumentContext.web WebContext.diagonal) 0)⟩

/-! ## Classical specialization in the committed finite-kernel vocabulary -/

section ClassicalKernel

open OPH.Thermodynamics (push heatBath fiberMass fiberMass_pos)

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- **Classical selection.**  A nonnegative finite kernel `K` whose row sums
are the indicator of `S` (induced classical effect) and which fixes, under
the committed pushforward `push`, every probability law certain of `S`
(mass one on `S`) is the restriction kernel of `S`: the classical Lüders
component `K x y = if x ∈ S ∧ x = y then 1 else 0`. -/
theorem classical_subkernel_fixing_certain_laws_is_restriction
    (S : Finset Ω) (K : Ω → Ω → ℝ) (hnn : ∀ x y, 0 ≤ K x y)
    (heff : ∀ x, ∑ y, K x y = if x ∈ S then 1 else 0)
    (hfix : ∀ p : Ω → ℝ, (∀ x, 0 ≤ p x) → ∑ x, p x = 1 → ∑ x ∈ S, p x = 1 →
      push p K = p)
    (x y : Ω) : K x y = if x ∈ S ∧ x = y then 1 else 0 := by
  by_cases hx : x ∈ S
  · -- The point law at `x` is certain of `S`; pushing it gives row `x`.
    have hδ : push (Pi.single x (1 : ℝ)) K = Pi.single x 1 := by
      refine hfix _ (fun z => ?_) ?_ ?_
      · rw [Pi.single_apply]
        split_ifs <;> norm_num
      · simp
      · rw [Finset.sum_eq_single x]
        · simp
        · intro z _ hz
          exact Pi.single_eq_of_ne hz 1
        · intro hxS
          exact absurd hx hxS
    have hrow : push (Pi.single x (1 : ℝ)) K y = K x y := by
      unfold push
      rw [Finset.sum_eq_single x]
      · simp
      · intro z _ hz
        rw [Pi.single_eq_of_ne hz, zero_mul]
      · intro hxu
        exact absurd (Finset.mem_univ x) hxu
    rw [← hrow, hδ, Pi.single_apply]
    by_cases hxy : x = y
    · subst hxy
      simp [hx]
    · have hyx : ¬ y = x := fun h => hxy h.symm
      simp [hxy, hyx]
  · -- Row `x` sums to zero with nonnegative entries.
    have hsum : ∑ y, K x y = 0 := by rw [heff x, if_neg hx]
    have h0 := (Finset.sum_eq_zero_iff_of_nonneg (fun y _ => hnn x y)).mp hsum y
      (Finset.mem_univ y)
    rw [h0]
    simp [hx]

/-- **The `FiniteConditionalRepair` kernel is not an instrument component.**
The committed conditional-resampling repair kernel `heatBath π b` (faithful
reference `π`) does not fix the point law at `x` whenever the fibre of `x`
contains a second point `y`: the pushed law has positive mass at `y`.  Only
this kernel is examined; other committed kernels are listed in the module
header as not examined. -/
theorem push_single_heatBath_ne_single {B : Type*} [DecidableEq B]
    {π : Ω → ℝ} {b : Ω → B} (hπ : ∀ x, 0 < π x) {x y : Ω} (hxy : x ≠ y)
    (hb : b y = b x) :
    push (Pi.single x (1 : ℝ)) (heatBath π b) ≠ Pi.single x 1 := by
  intro h
  have hy := congrFun h y
  have hrow : push (Pi.single x (1 : ℝ)) (heatBath π b) y = heatBath π b x y := by
    unfold push
    rw [Finset.sum_eq_single x]
    · simp
    · intro z _ hz
      rw [Pi.single_eq_of_ne hz, zero_mul]
    · intro hxu
      exact absurd (Finset.mem_univ x) hxu
  rw [hrow, Pi.single_eq_of_ne (Ne.symm hxy)] at hy
  unfold heatBath at hy
  rw [if_pos hb] at hy
  exact (div_pos (hπ y) (fiberMass_pos hπ x)).ne' hy

end ClassicalKernel

end

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms HasKrausFormWithEffect.isCompletelyPositive
#print axioms HasKrausFormWithEffect.trace_apply
#print axioms kraus_mul_compl_eq_zero
#print axioms kraus_eq_mul_event
#print axioms rank_one_decomposition_collinear
#print axioms smul_vecMulVec_mem_certainStates
#print axioms apply_vecMulVec_of_fixesCertainStates
#print axioms apply_compress_eq_of_fixesCertainStates
#print axioms kraus_component_fixing_certain_states_eq_compress
#print axioms instrument_component_fixing_certain_states_is_lueders
#print axioms kraus_eq_smul_event_of_compress
#print axioms kraus_eq_smul_event_of_fixesCertainStates
#print axioms instrument_with_projective_effects_fixing_certain_states_is_lueders
#print axioms luedersOutcomeMap_hasKrausFormWithEffect
#print axioms luedersOutcomeMap_fixesCertainStates
#print axioms hasKrausFormWithEffect_and_fixesCertainStates_iff_eq_luedersOutcomeMap
#print axioms PhaseInstrument.repeatable_iff_normalizedOutputCertain
#print axioms normalizedOutputCertain_of_fixesCertainStates
#print axioms luedersOutcomeMap_isometry_mul_hasKrausFormWithEffect
#print axioms luedersOutcomeMap_isometry_mul_normalizedOutputCertain
#print axioms isometry_twist_eq_smul_of_fixesCertainStates
#print axioms swapOutcome_hasKrausFormWithEffect_one
#print axioms swapOutcome_ne_luedersOutcomeMap_one
#print axioms swapOutcome_not_fixesCertainStates_one
#print axioms swapOutcome_normalizedOutputCertain_one
#print axioms committed_repeatability_clause_does_not_select_lueders
#print axioms luedersChannel_hasKrausFormWithEffect_one
#print axioms effect_does_not_select_but_certain_states_do
#print axioms PhaseInstrument.eq_of_outcomeMap_eq
#print axioms phaseInstrument_eq_lueders_of_kraus_fixesCertainStates
#print axioms phaseInstrument_eq_lueders_iff
#print axioms swapTwistedPhaseInstrument_hasKrausFormWithEffect
#print axioms swapTwistedPhaseInstrument_not_fixesCertainStates_diagonal
#print axioms committed_instruments_certain_state_invariance_selects_lueders
#print axioms classical_subkernel_fixing_certain_laws_is_restriction
#print axioms push_single_heatBath_ne_single

end EventAlgebra
