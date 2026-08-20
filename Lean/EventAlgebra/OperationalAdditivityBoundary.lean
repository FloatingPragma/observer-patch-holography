import EventAlgebra.OperationalPhaseAttainment

set_option autoImplicit false

namespace EventAlgebra

/-!
# The additivity boundary of the committed declared-effect count fixture

This module locates register row PR-03 (operational effect additivity)
exactly on the committed static count fixture of
`EventAlgebra.OperationalPhaseAttainment`: the legacy-named inhabitant
`attainedModel` with its eight committed contexts, committed effects, and
generated expected-frequency literals.

**Exact claim.**  Three theorem groups, composed in one receipt
(`operationalAdditivityBoundary_receipt`).

1. **The committed coexistent sums, enumerated exactly.**  Every committed
   outcome effect has trace one, and in dimension two a sum of two trace-one
   effects is an effect precisely when it equals the sure effect.  Hence a
   coexistent effect sum formed within the committed effect set is exactly a
   unit resolution: `IsEffect (E + F)` holds for committed `E`, `F` if and
   only if `E + F = 1` (`committed_coexistent_sum_iff`), and every context
   supplies one such pair (`committed_coexistent_sums_nonvacuous`).

2. **Born form on the static fixture, by extension and by forcing.**  The
   Born functional of the declared diagonal matrix (`bornRunValuation`) is a full
   additive effect valuation in the sense of the committed finite
   Busch-Gleason theorem (`bornRunValuation_isEffectValuation`), and its
   restriction to the committed effect set reproduces every fixture count
   frequency as an exact rational identity (`bornRunValuation_matches`).  On
   the static extension, additivity on all coexistent sums of the full effect
   algebra therefore holds as a theorem of the extension, with no additivity
   premise on the fixture data.  Conversely, full additivity together with
   the fixture counts forces this valuation and no other: every additive
   effect valuation whose committed-effect values are the fixture count
   frequencies agrees with `bornRunValuation` on all effects, and every state
   whose Born weights match the fixture frequencies is the committed
   declared matrix (`additive_valuation_matching_counts_is_born`,
   `born_state_matching_counts_unique`).  The forcing consumes the full
   quantifier of PR-03; nothing here derives that quantifier.

3. **The coverage gap is essential.**  The committed coexistent sums cannot
   carry the PR-03 obligation: any assignment with `v 1 = 1` and per-context
   normalisation is additive on every committed coexistent sum
   (`committed_sums_carry_only_normalisation`), so restricted additivity is
   equivalent to normalisation and forces nothing.  A machine-checked
   countermodel witnesses the gap: `producedCubicValuation` deforms the run
   Born functional by the transverse cubic direction of the committed
   finite-web counterexample, shaped to vanish on the committed effect set.
   It reproduces every fixture count frequency, maps every effect into
   `[0, 1]`, sends the sure effect to one, obeys the complement rule, and is
   additive on every committed coexistent sum, yet it fails additivity at the
   explicit effect pair `halfWitness + halfWitness = witnessSum` with the
   exact rational gap `35/64` against `143/256`
   (`producedCubicValuation_not_additive`), is not an additive effect
   valuation (`producedCubicValuation_not_effectValuation`), and is the Born
   functional of no state (`producedCubicValuation_not_born`).

**Premise consumption, named per the register.**

* PR-02 (committed): the declared diagonal representation `committedRunState` on the
  committed algebra-state surface.
* PR-04 (consumed through the imported attainment module, axiomatize
  disposition): the committed effect set carries the declared phase
  effect; this module adds no effect and no counts.
* PR-03 (registered premise, the subject of this module): consumed only as
  the explicit hypothesis `IsEffectValuation v` of the forcing theorems.
  The countermodel proves that the static fixture does not supply the
  premise: the committed required coexistent sums are the unit resolutions,
  and additivity on them is normalisation.

**Falsifier.** Any failure of the exact identities: a committed effect of
trace other than one, a committed coexistent sum that is an effect without
being the sure effect, a fixture frequency differing from the declared Born
Born weight, a countermodel value leaving `[0, 1]` on some effect, or the
countermodel satisfying full additivity or Born form after all.

**Nonclaims and boundary.**  No instrument, no new effect, no new counts, and no claim
past the committed effect set.  The register row PR-03 is not discharged:
its quantifier ranges over all required coexistent sums of the registered
public frame, the committed static fixture realizes only the unit
resolutions, and the countermodel proves that this restricted family cannot
substitute for the full quantifier.  The scientific status is exact algebraic
conformance of a declared-effect count fixture.  CP outcome maps, a summed
channel, source preparation, and readback provenance remain open separately
from PR-03 additivity.
-/

open Matrix
open OPH.QFT
open scoped ComplexOrder

noncomputable section

/-! ## Trace one on every committed effect -/

/-- Complexification preserves the trace. -/
theorem complexifyRealMatrix_trace (M : Matrix (Fin 2) (Fin 2) ℝ) :
    (complexifyRealMatrix M).trace = ((M.trace : ℝ) : ℂ) := by
  simp [complexifyRealMatrix, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

theorem recordProjector_trace : recordProjector.trace = 1 := by
  norm_num [recordProjector, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

/-- Every conjugated projector has trace one: conjugation by an orthogonal
image preserves the trace of the record projector. -/
theorem conjProjector_trace (g : Fin 6) : (conjProjector g).trace = 1 := by
  unfold conjProjector
  rw [Matrix.trace_mul_comm, ← Matrix.mul_assoc, (gaugeIrrep_unitary g).2,
    one_mul, recordProjector_trace]

theorem sourcePhaseLift_trace : sourcePhaseLift.trace = 1 := by
  rw [sourcePhaseLift_eq_rhoYPlus]
  norm_num [rhoYPlus, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

/-- Every committed outcome-`0` effect has trace one. -/
theorem committedContextEffect_trace (c : InstrumentContext) :
    (committedContextEffect c).trace = 1 := by
  cases c with
  | web wc =>
      show (complexifyRealMatrix (webContextProjector wc)).trace = 1
      cases wc with
      | diagonal =>
          rw [webContextProjector_diagonal, complexifyRealMatrix_trace,
            recordProjector_trace]
          norm_num
      | conjugated g =>
          rw [webContextProjector_conjugated, complexifyRealMatrix_trace,
            conjProjector_trace]
          norm_num
  | phase => exact sourcePhaseLift_trace

/-- Every committed effect of the attained inhabitant, both outcomes in all
eight contexts, has trace one. -/
theorem attainedModel_effect_trace (c : InstrumentContext) (i : Fin 2) :
    (attainedModel.effect c i).trace = 1 := by
  fin_cases i
  · show (committedEffectPair c 0).trace = 1
    simpa [committedEffectPair] using committedContextEffect_trace c
  · show (committedEffectPair c 1).trace = 1
    simp only [committedEffectPair, Matrix.cons_val_one, Matrix.cons_val_fin_one]
    rw [trace_sub, Matrix.trace_one, committedContextEffect_trace c]
    norm_num

/-! ## The committed coexistent sums are exactly the unit resolutions -/

/-- In dimension two, a sum of two trace-one effects is an effect only when
it is the sure effect: the complement of the sum is positive semidefinite
with zero trace, hence zero. -/
theorem effect_sum_of_trace_one {A B : Matrix (Fin 2) (Fin 2) ℂ}
    (hA : A.trace = 1) (hB : B.trace = 1) (hAB : IsEffect (A + B)) :
    A + B = 1 := by
  have htr : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - (A + B)).trace = 0 := by
    rw [trace_sub, trace_add, hA, hB, Matrix.trace_one]
    norm_num
  have hz := hAB.2.trace_eq_zero_iff.mp htr
  have := sub_eq_zero.mp hz
  exact this.symm

/-- **The exact enumeration of the committed coexistent sums.**  A sum of two
committed effects is an effect precisely when it is the sure effect, so the
coexistent sums formed within the committed effect set are exactly the unit
resolutions. -/
theorem committed_coexistent_sum_iff (c c' : InstrumentContext) (i i' : Fin 2) :
    IsEffect (attainedModel.effect c i + attainedModel.effect c' i') ↔
      attainedModel.effect c i + attainedModel.effect c' i' = 1 := by
  constructor
  · exact effect_sum_of_trace_one (attainedModel_effect_trace c i)
      (attainedModel_effect_trace c' i')
  · intro h
    rw [h]
    exact isEffect_one

/-- Nonvacuity: every committed context supplies one coexistent sum, its own
binary resolution. -/
theorem committed_coexistent_sums_nonvacuous (c : InstrumentContext) :
    IsEffect (attainedModel.effect c 0 + attainedModel.effect c 1) := by
  rw [attainedModel.effect_complete c]
  exact isEffect_one

/-! ## The static fixture frequencies and the declared-matrix Born valuation -/

/-- The generated expected frequency of each committed context and outcome:
outcome `0` carries the outcome-`0` frequency, outcome `1` the complementary
frequency, both exact rationals of the committed count literals. -/
def modelFrequency (c : InstrumentContext) : Fin 2 → ℂ :=
  ![((attainedModel.counts c).1 : ℂ) /
      (((attainedModel.counts c).1 + (attainedModel.counts c).2 : ℕ) : ℂ),
    ((attainedModel.counts c).2 : ℂ) /
      (((attainedModel.counts c).1 + (attainedModel.counts c).2 : ℕ) : ℂ)]

/-- The Born weight of every committed effect under the declared diagonal matrix
is the fixture frequency: the stored fit law of the static inhabitant, in
both outcomes. -/
theorem committed_frequency_born (c : InstrumentContext) (i : Fin 2) :
    bornWeight committedRunState (attainedModel.effect c i) =
      modelFrequency c i := by
  fin_cases i
  · have h := attainedModel.born_matches c
    rw [attainedModel_prep_eq] at h
    simpa [modelFrequency, binaryFrequency] using h
  · have h := attainedModel.born_matches_snd c
    rw [attainedModel_prep_eq] at h
    simpa [modelFrequency] using h

/-- The Born functional of the declared diagonal matrix, as a real assignment on
matrices: the candidate valuation extending the static fixture. -/
def bornRunValuation (E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (bornWeight committedRunState E).re

/-- **Additivity as a theorem of the extension.**  The run Born functional
is an additive effect valuation in the sense of the committed finite
Busch-Gleason theorem: on the full effect algebra the PR-03 additivity shape
holds of this extension with no additivity premise. -/
theorem bornRunValuation_isEffectValuation :
    IsEffectValuation bornRunValuation :=
  isEffectValuation_born committedRunState_isState

/-- The run Born valuation restricts on the committed effect set to exactly
the fixture frequencies. -/
theorem bornRunValuation_matches (c : InstrumentContext) (i : Fin 2) :
    ((bornRunValuation (attainedModel.effect c i) : ℝ) : ℂ) =
      modelFrequency c i := by
  rw [bornRunValuation,
    bornWeight_eq_re committedRunState_isState.1.isHermitian
      (attainedModel.effect_isEffect c i).isHermitian,
    committed_frequency_born]

/-- **The static exact-fit values pin the declared matrix.** Every state whose Born
weights on the committed effect set are the fixture frequencies is
the declared matrix: the committed effect set is tomographically
complete on the fixed-trace slice. -/
theorem born_state_matching_counts_unique {ρ : Matrix (Fin 2) (Fin 2) ℂ}
    (hρ : IsState ρ)
    (h : ∀ (c : InstrumentContext) (i : Fin 2),
      bornWeight ρ (attainedModel.effect c i) = modelFrequency c i) :
    ρ = committedRunState := by
  refine ideal_phase_model_completes_tomography attainedModel ?_ ?_
  · rw [hρ.2, committedRunState_isState.2]
  · intro c i
    rw [h c i, committed_frequency_born]

/-- **Forcing under full additivity.**  Every additive effect valuation
whose committed-effect values are the fixture frequencies is the run
Born valuation on all effects.  The hypothesis `IsEffectValuation v` is the
full PR-03 additivity premise; the theorem consumes it and does not derive
it. -/
theorem additive_valuation_matching_counts_is_born
    {v : Matrix (Fin 2) (Fin 2) ℂ → ℝ} (hv : IsEffectValuation v)
    (h : ∀ (c : InstrumentContext) (i : Fin 2),
      ((v (attainedModel.effect c i) : ℝ) : ℂ) = modelFrequency c i) :
    ∀ E, IsEffect E → v E = bornRunValuation E := by
  obtain ⟨ρ, ⟨hρ, hrep⟩, -⟩ := finite_busch_gleason hv
  have hmatch : ∀ (c : InstrumentContext) (i : Fin 2),
      bornWeight ρ (attainedModel.effect c i) = modelFrequency c i := by
    intro c i
    have hE := attainedModel.effect_isEffect c i
    rw [← bornWeight_eq_re hρ.1.isHermitian hE.isHermitian,
      ← hrep _ hE, h c i]
  have hstate := born_state_matching_counts_unique hρ hmatch
  intro E hE
  rw [hrep E hE, hstate, bornRunValuation]

/-! ## The committed coexistent sums carry only normalisation -/

/-- **Restricted additivity is normalisation.**  Any real assignment on
matrices that sends the sure effect to one and is normalised in every
committed context is additive on every coexistent sum formed within the
committed effect set.  The committed sums therefore cannot carry the PR-03
obligation: additivity restricted to them is equivalent to per-context
normalisation. -/
theorem committed_sums_carry_only_normalisation
    (v : Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hone : v 1 = 1)
    (hnorm : ∀ c : InstrumentContext,
      v (attainedModel.effect c 0) + v (attainedModel.effect c 1) = 1)
    {c c' : InstrumentContext} {i i' : Fin 2}
    (hco : IsEffect (attainedModel.effect c i + attainedModel.effect c' i')) :
    v (attainedModel.effect c i + attainedModel.effect c' i') =
      v (attainedModel.effect c i) + v (attainedModel.effect c' i') := by
  have hsum := (committed_coexistent_sum_iff c c' i i').mp hco
  have hA : attainedModel.effect c i = 1 - attainedModel.effect c' i' := by
    rw [← hsum]
    abel
  rw [hsum, hone]
  have hi : i' = 0 ∨ i' = 1 := by omega
  rcases hi with rfl | rfl
  · have hswap : attainedModel.effect c i = attainedModel.effect c' 1 := by
      rw [hA, attainedModel.effect_one c']
    rw [hswap]
    linarith [hnorm c']
  · have hswap : attainedModel.effect c i = attainedModel.effect c' 0 := by
      rw [hA, attainedModel.effect_one c']
      abel
    rw [hswap]
    linarith [hnorm c']

/-! ## The countermodel: a produced-count-compatible cubic valuation

The deformation direction is the transverse cubic of the committed
finite-web counterexample (`EventAlgebra.FiniteWebBornNoGo`), read in the
imaginary off-diagonal coordinate and shaped to vanish at the coordinate
values `0` and `±1/2` realized by the committed effect set. -/

/-- The transverse cubic deformation: an odd polynomial in the imaginary
off-diagonal coordinate that vanishes exactly at `0` and `±1/2`. -/
def cubicDeformation (E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (E 1 0).im * (1 / 4 - (E 1 0).im ^ 2)

/-- The countermodel valuation: the run Born functional deformed by the
shaped transverse cubic. -/
def producedCubicValuation (E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  bornRunValuation E + cubicDeformation E

/-- The run Born valuation in matrix coordinates. -/
theorem bornRunValuation_coords (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornRunValuation E =
      (111 * (E 0 0).re + 68 * (E 1 1).re) / 179 := by
  have h : bornWeight committedRunState E =
      (((111 / 179 : ℝ)) : ℂ) * E 0 0 + (((68 / 179 : ℝ)) : ℂ) * E 1 1 := by
    simp only [bornWeight, committedRunState, Matrix.trace, Matrix.diag,
      Matrix.mul_apply, Fin.sum_univ_two, Matrix.diagonal_apply,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_fin_one]
    push_cast
    norm_num
  rw [bornRunValuation, h]
  simp [Complex.add_re]
  ring

/-- The complement rule holds exactly, on all matrices. -/
theorem producedCubicValuation_compl (E : Matrix (Fin 2) (Fin 2) ℂ) :
    producedCubicValuation (1 - E) = 1 - producedCubicValuation E := by
  have h10 : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - E) 1 0 = - E 1 0 := by
    simp [Matrix.sub_apply, Matrix.one_apply_ne (by decide : (1 : Fin 2) ≠ 0)]
  have hb : bornRunValuation (1 - E) = 1 - bornRunValuation E := by
    unfold bornRunValuation
    rw [bornWeight_sub, bornWeight_one committedRunState_isState,
      Complex.sub_re, Complex.one_re]
  have hd : cubicDeformation (1 - E) = - cubicDeformation E := by
    unfold cubicDeformation
    rw [h10, Complex.neg_im]
    ring
  unfold producedCubicValuation
  rw [hb, hd]
  ring

/-- The countermodel sends the sure effect to one. -/
theorem producedCubicValuation_one :
    producedCubicValuation (1 : Matrix (Fin 2) (Fin 2) ℂ) = 1 := by
  have h10 : ((1 : Matrix (Fin 2) (Fin 2) ℂ)) 1 0 = 0 :=
    Matrix.one_apply_ne (by decide)
  unfold producedCubicValuation cubicDeformation bornRunValuation
  rw [bornWeight_one committedRunState_isState, h10]
  simp

/-- The countermodel sends the null effect to zero. -/
theorem producedCubicValuation_zero :
    producedCubicValuation (0 : Matrix (Fin 2) (Fin 2) ℂ) = 0 := by
  unfold producedCubicValuation cubicDeformation bornRunValuation
  simp [bornWeight]

/-- The shaped cubic vanishes at the three realized coordinate values. -/
theorem cubicDeformation_eq_zero {E : Matrix (Fin 2) (Fin 2) ℂ}
    (h : (E 1 0).im = 0 ∨ (E 1 0).im = 1 / 2 ∨ (E 1 0).im = -(1 / 2)) :
    cubicDeformation E = 0 := by
  unfold cubicDeformation
  rcases h with h | h | h <;> rw [h] <;> norm_num

/-- The imaginary off-diagonal coordinate of every committed effect is `0`,
`1/2`, or `-1/2`. -/
theorem committed_effect_im (c : InstrumentContext) (i : Fin 2) :
    ((attainedModel.effect c i) 1 0).im = 0 ∨
      ((attainedModel.effect c i) 1 0).im = 1 / 2 ∨
      ((attainedModel.effect c i) 1 0).im = -(1 / 2) := by
  have hlift : (sourcePhaseLift 1 0).im = 1 / 2 := by
    rw [sourcePhaseLift_eq_rhoYPlus]
    have h : rhoYPlus 1 0 = Complex.I / 2 := by
      norm_num [rhoYPlus]
    rw [h]
    norm_num [Complex.div_im, Complex.normSq_apply]
  have hzero : ∀ wc : WebContext,
      ((attainedModel.effect (InstrumentContext.web wc) 0) 1 0).im = 0 := by
    intro wc
    rw [attainedModel.web_effect wc]
    exact Complex.ofReal_im _
  have h0 : ∀ c : InstrumentContext,
      ((attainedModel.effect c 0) 1 0).im = 0 ∨
        ((attainedModel.effect c 0) 1 0).im = 1 / 2 := by
    intro c
    cases c with
    | web wc => exact Or.inl (hzero wc)
    | phase =>
        right
        rw [attainedModel_phase_effect]
        exact hlift
  have hi : i = 0 ∨ i = 1 := by omega
  rcases hi with rfl | rfl
  · rcases h0 c with h | h
    · exact Or.inl h
    · exact Or.inr (Or.inl h)
  · have hcompl : ((attainedModel.effect c 1) 1 0).im =
        - ((attainedModel.effect c 0) 1 0).im := by
      rw [attainedModel.effect_one c]
      have h10 : ((1 : Matrix (Fin 2) (Fin 2) ℂ)
          - attainedModel.effect c 0) 1 0 = - (attainedModel.effect c 0) 1 0 := by
        simp [Matrix.sub_apply,
          Matrix.one_apply_ne (by decide : (1 : Fin 2) ≠ 0)]
      rw [h10, Complex.neg_im]
    rcases h0 c with h | h
    · exact Or.inl (by rw [hcompl, h, neg_zero])
    · exact Or.inr (Or.inr (by rw [hcompl, h]))

/-- The deformation vanishes on the whole committed effect set. -/
theorem cubicDeformation_committed (c : InstrumentContext) (i : Fin 2) :
    cubicDeformation (attainedModel.effect c i) = 0 :=
  cubicDeformation_eq_zero (committed_effect_im c i)

/-- The countermodel agrees with the run Born valuation on every committed
effect. -/
theorem producedCubicValuation_committed (c : InstrumentContext) (i : Fin 2) :
    producedCubicValuation (attainedModel.effect c i) =
      bornRunValuation (attainedModel.effect c i) := by
  unfold producedCubicValuation
  rw [cubicDeformation_committed]
  ring

/-- The countermodel reproduces every fixture frequency. -/
theorem producedCubicValuation_matches (c : InstrumentContext) (i : Fin 2) :
    ((producedCubicValuation (attainedModel.effect c i) : ℝ) : ℂ) =
      modelFrequency c i := by
  rw [producedCubicValuation_committed, bornRunValuation_matches]

/-- The countermodel is normalised in every committed context. -/
theorem producedCubicValuation_normalised (c : InstrumentContext) :
    producedCubicValuation (attainedModel.effect c 0) +
      producedCubicValuation (attainedModel.effect c 1) = 1 := by
  rw [attainedModel.effect_one c, producedCubicValuation_compl]
  ring

/-- The countermodel is additive on every coexistent sum formed within the
committed effect set. -/
theorem producedCubicValuation_additive_on_committed_sums
    {c c' : InstrumentContext} {i i' : Fin 2}
    (hco : IsEffect (attainedModel.effect c i + attainedModel.effect c' i')) :
    producedCubicValuation
        (attainedModel.effect c i + attainedModel.effect c' i') =
      producedCubicValuation (attainedModel.effect c i) +
        producedCubicValuation (attainedModel.effect c' i') :=
  committed_sums_carry_only_normalisation producedCubicValuation
    producedCubicValuation_one producedCubicValuation_normalised hco

/-! ## The countermodel takes values in the probability interval -/

/-- Entry bounds of a dimension-two effect: the diagonal real parts lie in
`[0, 1]` and the squared imaginary off-diagonal coordinate is dominated by
both diagonal products, from the determinant nonnegativity of the effect and
of its complement. -/
theorem effect_entry_bounds {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : IsEffect E) :
    0 ≤ (E 0 0).re ∧ (E 0 0).re ≤ 1 ∧ 0 ≤ (E 1 1).re ∧ (E 1 1).re ≤ 1 ∧
      (E 1 0).im ^ 2 ≤ (E 0 0).re * (E 1 1).re ∧
      (E 1 0).im ^ 2 ≤ (1 - (E 0 0).re) * (1 - (E 1 1).re) := by
  have h00 := hE.1.diag_nonneg (i := 0)
  have h11 := hE.1.diag_nonneg (i := 1)
  have hc00 := hE.2.diag_nonneg (i := 0)
  have hc11 := hE.2.diag_nonneg (i := 1)
  have h00' := Complex.le_def.mp h00
  have h11' := Complex.le_def.mp h11
  have hone00 : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - E) 0 0 = 1 - E 0 0 := by
    simp [Matrix.sub_apply]
  have hone11 : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - E) 1 1 = 1 - E 1 1 := by
    simp [Matrix.sub_apply]
  rw [hone00] at hc00
  rw [hone11] at hc11
  have hc00' := Complex.le_def.mp hc00
  have hc11' := Complex.le_def.mp hc11
  have h01 : E 0 1 = star (E 1 0) := by
    conv_lhs => rw [← hE.isHermitian.eq]
    rw [Matrix.conjTranspose_apply]
  have hdet := hE.1.det_nonneg
  rw [Matrix.det_fin_two, h01] at hdet
  have hdet' := Complex.le_def.mp hdet
  have hcdet := hE.2.det_nonneg
  rw [Matrix.det_fin_two] at hcdet
  have hcd01 : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - E) 0 1 = - E 0 1 := by
    simp [Matrix.sub_apply, Matrix.one_apply_ne (by decide : (0 : Fin 2) ≠ 1)]
  have hcd10 : ((1 : Matrix (Fin 2) (Fin 2) ℂ) - E) 1 0 = - E 1 0 := by
    simp [Matrix.sub_apply, Matrix.one_apply_ne (by decide : (1 : Fin 2) ≠ 0)]
  rw [hone00, hone11, hcd01, hcd10, h01] at hcdet
  have hcdet' := Complex.le_def.mp hcdet
  set a := (E 0 0).re
  set b := (E 1 1).re
  set x := (E 1 0).re
  set y := (E 1 0).im
  have hia : (E 0 0).im = 0 := (h00'.2).symm
  have hib : (E 1 1).im = 0 := (h11'.2).symm
  have hstar : (star (E 1 0) * E 1 0).re = x ^ 2 + y ^ 2 := by
    simp [Complex.mul_re]
    ring
  have hdre : (E 0 0 * E 1 1 - star (E 1 0) * E 1 0).re =
      a * b - (x ^ 2 + y ^ 2) := by
    rw [Complex.sub_re, Complex.mul_re, hia, hib, hstar]
    ring
  have hd1 : 0 ≤ a * b - (x ^ 2 + y ^ 2) := by
    have h := hdet'.1
    rw [hdre] at h
    simpa using h
  have hcre : ((1 - E 0 0) * (1 - E 1 1) - - star (E 1 0) * - (E 1 0)).re =
      (1 - a) * (1 - b) - (x ^ 2 + y ^ 2) := by
    have hneg : (- star (E 1 0) * - (E 1 0)) = star (E 1 0) * E 1 0 := by ring
    rw [hneg, Complex.sub_re, Complex.mul_re, hstar, Complex.sub_re,
      Complex.sub_im, Complex.sub_re, Complex.sub_im, Complex.one_re,
      Complex.one_im, hia, hib]
    ring
  have hd2 : 0 ≤ (1 - a) * (1 - b) - (x ^ 2 + y ^ 2) := by
    have h := hcdet'.1
    rw [hcre] at h
    simpa using h
  have ha0 : 0 ≤ a := by simpa using h00'.1
  have hb0 : 0 ≤ b := by simpa using h11'.1
  have ha1 : a ≤ 1 := by
    have h := hc00'.1
    rw [Complex.sub_re, Complex.one_re] at h
    simp only [Complex.zero_re] at h
    linarith
  have hb1 : b ≤ 1 := by
    have h := hc11'.1
    rw [Complex.sub_re, Complex.one_re] at h
    simp only [Complex.zero_re] at h
    linarith
  refine ⟨ha0, ha1, hb0, hb1, ?_, ?_⟩
  · nlinarith [sq_nonneg x]
  · nlinarith [sq_nonneg x]

/-- The countermodel is nonnegative on every effect. -/
theorem producedCubicValuation_nonneg {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : IsEffect E) : 0 ≤ producedCubicValuation E := by
  obtain ⟨ha0, ha1, hb0, hb1, hy1, hy2⟩ := effect_entry_bounds hE
  have hval : producedCubicValuation E =
      (111 * (E 0 0).re + 68 * (E 1 1).re) / 179 +
        (E 1 0).im * (1 / 4 - (E 1 0).im ^ 2) := by
    unfold producedCubicValuation cubicDeformation
    rw [bornRunValuation_coords]
  set a := (E 0 0).re
  set b := (E 1 1).re
  set y := (E 1 0).im
  have h4 : y ^ 2 * y ^ 2 ≤ (a * b) * ((1 - a) * (1 - b)) :=
    mul_le_mul hy1 hy2 (sq_nonneg y) (mul_nonneg ha0 hb0)
  have hy14 : y ^ 2 ≤ 1 / 4 := by
    nlinarith [sq_nonneg (2 * a - 1), sq_nonneg (2 * b - 1), sq_nonneg y,
      sq_nonneg (y ^ 2 - 1 / 4), sq_nonneg (a - b), sq_nonneg (a + b - 1)]
  have h2y : 4 * y ^ 2 ≤ (a + b) ^ 2 := by
    nlinarith [sq_nonneg (a - b)]
  have hab : 0 ≤ a + b := by linarith
  have hylow : -(a + b) ≤ 2 * y := by
    nlinarith [h2y, hab]
  rw [hval]
  rcases le_or_gt 0 y with hy | hy
  · have hq : 0 ≤ 1 / 4 - y ^ 2 := by linarith
    have := mul_nonneg hy hq
    have hfrac : 0 ≤ (111 * a + 68 * b) / 179 := by positivity
    linarith
  · have hq : 0 ≤ 1 / 4 - y ^ 2 := by linarith
    have hcube : y * (1 / 4 - y ^ 2) ≥ y * (1 / 4) := by
      nlinarith [mul_nonneg (neg_nonneg.mpr hy.le) (sq_nonneg y)]
    have hlin : y * (1 / 4) ≥ -(a + b) / 8 := by linarith
    have hdom : (111 * a + 68 * b) / 179 ≥ (a + b) / 8 := by
      nlinarith
    linarith

/-- The countermodel is at most one on every effect. -/
theorem producedCubicValuation_le_one {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : IsEffect E) : producedCubicValuation E ≤ 1 := by
  have h := producedCubicValuation_nonneg hE.compl
  rw [producedCubicValuation_compl] at h
  linarith

/-! ## The explicit additivity failure -/

/-- The half witness: the effect `(1/4) rhoYPlus + (1/8) 1`, with imaginary
off-diagonal coordinate `1/8`, outside the committed effect set. -/
def halfWitness : Matrix (Fin 2) (Fin 2) ℂ :=
  (1 / 4 : ℝ) • rhoYPlus + (1 / 8 : ℝ) • 1

/-- The witness sum: `halfWitness + halfWitness`, with imaginary
off-diagonal coordinate `1/4`, again an effect. -/
def witnessSum : Matrix (Fin 2) (Fin 2) ℂ := halfWitness + halfWitness

/-- The two Pauli-Y projectors resolve the identity. -/
theorem rhoY_add : rhoYPlus + rhoYMinus = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [rhoYPlus, rhoYMinus, Matrix.add_apply, Matrix.one_apply]

theorem halfWitness_isEffect : IsEffect halfWitness := by
  constructor
  · exact (posSemidef_real_smul rhoYPlus_isEvent.posSemidef
      (by norm_num)).add
      (posSemidef_real_smul Matrix.PosSemidef.one (by norm_num))
  · have hminus : rhoYMinus = 1 - rhoYPlus := by
      rw [← rhoY_add]
      abel
    have hkey : (1 : Matrix (Fin 2) (Fin 2) ℂ) - halfWitness =
        (1 / 4 : ℝ) • rhoYMinus + (5 / 8 : ℝ) • 1 := by
      rw [hminus]
      unfold halfWitness
      module
    rw [hkey]
    exact (posSemidef_real_smul rhoYMinus_isEvent.posSemidef
      (by norm_num)).add
      (posSemidef_real_smul Matrix.PosSemidef.one (by norm_num))

theorem witnessSum_isEffect : IsEffect witnessSum := by
  constructor
  · exact (halfWitness_isEffect.1).add (halfWitness_isEffect.1)
  · have hminus : rhoYMinus = 1 - rhoYPlus := by
      rw [← rhoY_add]
      abel
    have hkey : (1 : Matrix (Fin 2) (Fin 2) ℂ) - witnessSum =
        (1 / 2 : ℝ) • rhoYMinus + (1 / 4 : ℝ) • 1 := by
      rw [hminus]
      unfold witnessSum halfWitness
      module
    rw [hkey]
    exact (posSemidef_real_smul rhoYMinus_isEvent.posSemidef
      (by norm_num)).add
      (posSemidef_real_smul Matrix.PosSemidef.one (by norm_num))

/-- Coordinates of the half witness. -/
theorem halfWitness_coords :
    (halfWitness 0 0).re = 1 / 4 ∧ (halfWitness 1 1).re = 1 / 4 ∧
      (halfWitness 1 0).im = 1 / 8 := by
  unfold halfWitness
  have h00 : rhoYPlus 0 0 = 1 / 2 := by norm_num [rhoYPlus]
  have h11 : rhoYPlus 1 1 = 1 / 2 := by norm_num [rhoYPlus]
  have h10 : rhoYPlus 1 0 = Complex.I / 2 := by norm_num [rhoYPlus]
  refine ⟨?_, ?_, ?_⟩ <;>
    simp [Matrix.add_apply, Matrix.smul_apply, h00, h11, h10,
      Complex.real_smul, Complex.mul_re, Complex.mul_im] <;>
    norm_num

/-- The countermodel value of the half witness: `143/512` exactly. -/
theorem producedCubicValuation_halfWitness :
    producedCubicValuation halfWitness = 143 / 512 := by
  obtain ⟨h00, h11, h10⟩ := halfWitness_coords
  unfold producedCubicValuation cubicDeformation
  rw [bornRunValuation_coords, h00, h11, h10]
  norm_num

/-- Coordinates of the witness sum. -/
theorem witnessSum_coords :
    (witnessSum 0 0).re = 1 / 2 ∧ (witnessSum 1 1).re = 1 / 2 ∧
      (witnessSum 1 0).im = 1 / 4 := by
  obtain ⟨h00, h11, h10⟩ := halfWitness_coords
  unfold witnessSum
  refine ⟨?_, ?_, ?_⟩ <;>
    simp [Matrix.add_apply, Complex.add_re, Complex.add_im, h00, h11, h10] <;>
    norm_num

/-- The countermodel value of the witness sum: `35/64` exactly. -/
theorem producedCubicValuation_witnessSum :
    producedCubicValuation witnessSum = 35 / 64 := by
  obtain ⟨h00, h11, h10⟩ := witnessSum_coords
  unfold producedCubicValuation cubicDeformation
  rw [bornRunValuation_coords, h00, h11, h10]
  norm_num

/-- **The explicit additivity failure.**  The countermodel is not additive
at the effect pair `halfWitness, halfWitness`: the sum receives `35/64`
while the summands total `143/256`. -/
theorem producedCubicValuation_not_additive :
    ∃ A B : Matrix (Fin 2) (Fin 2) ℂ, IsEffect A ∧ IsEffect B ∧
      IsEffect (A + B) ∧
      producedCubicValuation (A + B) ≠
        producedCubicValuation A + producedCubicValuation B := by
  refine ⟨halfWitness, halfWitness, halfWitness_isEffect,
    halfWitness_isEffect, witnessSum_isEffect, ?_⟩
  show producedCubicValuation witnessSum ≠ _
  rw [producedCubicValuation_witnessSum, producedCubicValuation_halfWitness]
  norm_num

/-- The countermodel is not an additive effect valuation. -/
theorem producedCubicValuation_not_effectValuation :
    ¬ IsEffectValuation producedCubicValuation := by
  intro hv
  have h := hv.additive halfWitness_isEffect halfWitness_isEffect
    witnessSum_isEffect
  have h' : producedCubicValuation witnessSum =
      producedCubicValuation halfWitness +
        producedCubicValuation halfWitness := h
  rw [producedCubicValuation_witnessSum,
    producedCubicValuation_halfWitness] at h'
  norm_num at h'

/-- The countermodel is the Born functional of no state. -/
theorem producedCubicValuation_not_born :
    ¬ ∃ ρ : Matrix (Fin 2) (Fin 2) ℂ, IsState ρ ∧
      ∀ E, IsEffect E → producedCubicValuation E = (bornWeight ρ E).re := by
  rintro ⟨ρ, hρ, hrep⟩
  have h1 := hrep halfWitness halfWitness_isEffect
  have h2 := hrep witnessSum witnessSum_isEffect
  have hadd : bornWeight ρ witnessSum =
      bornWeight ρ halfWitness + bornWeight ρ halfWitness :=
    bornWeight_add ρ halfWitness halfWitness
  have hre := congrArg Complex.re hadd
  rw [Complex.add_re] at hre
  rw [producedCubicValuation_witnessSum] at h2
  rw [producedCubicValuation_halfWitness] at h1
  rw [← h1, ← h2] at hre
  norm_num at hre

/-- The deformation is live away from the committed set: the countermodel
differs from the run Born valuation at the witness sum, `35/64` against
`1/2`. -/
theorem producedCubicValuation_deviates :
    producedCubicValuation witnessSum ≠ bornRunValuation witnessSum := by
  obtain ⟨h00, h11, -⟩ := witnessSum_coords
  rw [producedCubicValuation_witnessSum, bornRunValuation_coords, h00, h11]
  norm_num

/-! ## The composed receipt -/

/-- **The operational additivity boundary receipt.**  One conjunction:

1. the coexistent sums formed within the committed effect set are exactly
   the unit resolutions, and every context supplies one;
2. the run Born functional is a full additive effect valuation whose
   committed-effect values are the fixture frequencies;
3. full additivity plus the fixture frequencies forces the run Born valuation
   on all effects, and those values pin the matrix to the
   declared diagonal matrix;
4. the countermodel reproduces every fixture frequency, maps every
   effect into `[0, 1]`, sends the sure effect to one, and is additive on
   every committed coexistent sum;
5. the countermodel fails additivity at an explicit effect pair, is not an
   additive effect valuation, and is the Born functional of no state.

Boundary: clause 3 consumes the full PR-03 additivity premise as a
hypothesis; clauses 1, 4, and 5 prove that the committed static fixture
does not supply that premise, because its coexistent sums carry only
normalisation. -/
theorem operationalAdditivityBoundary_receipt :
    (∀ (c c' : InstrumentContext) (i i' : Fin 2),
      IsEffect (attainedModel.effect c i + attainedModel.effect c' i') ↔
        attainedModel.effect c i + attainedModel.effect c' i' = 1) ∧
    (∀ c : InstrumentContext,
      IsEffect (attainedModel.effect c 0 + attainedModel.effect c 1)) ∧
    IsEffectValuation bornRunValuation ∧
    (∀ (c : InstrumentContext) (i : Fin 2),
      ((bornRunValuation (attainedModel.effect c i) : ℝ) : ℂ) =
        modelFrequency c i) ∧
    (∀ (v : Matrix (Fin 2) (Fin 2) ℂ → ℝ), IsEffectValuation v →
      (∀ (c : InstrumentContext) (i : Fin 2),
        ((v (attainedModel.effect c i) : ℝ) : ℂ) = modelFrequency c i) →
      ∀ E, IsEffect E → v E = bornRunValuation E) ∧
    (∀ ρ : Matrix (Fin 2) (Fin 2) ℂ, IsState ρ →
      (∀ (c : InstrumentContext) (i : Fin 2),
        bornWeight ρ (attainedModel.effect c i) = modelFrequency c i) →
      ρ = committedRunState) ∧
    (∀ (c : InstrumentContext) (i : Fin 2),
      ((producedCubicValuation (attainedModel.effect c i) : ℝ) : ℂ) =
        modelFrequency c i) ∧
    (∀ E : Matrix (Fin 2) (Fin 2) ℂ, IsEffect E →
      0 ≤ producedCubicValuation E ∧ producedCubicValuation E ≤ 1) ∧
    producedCubicValuation 1 = 1 ∧
    (∀ (c c' : InstrumentContext) (i i' : Fin 2),
      IsEffect (attainedModel.effect c i + attainedModel.effect c' i') →
      producedCubicValuation
          (attainedModel.effect c i + attainedModel.effect c' i') =
        producedCubicValuation (attainedModel.effect c i) +
          producedCubicValuation (attainedModel.effect c' i')) ∧
    (∃ A B : Matrix (Fin 2) (Fin 2) ℂ, IsEffect A ∧ IsEffect B ∧
      IsEffect (A + B) ∧
      producedCubicValuation (A + B) ≠
        producedCubicValuation A + producedCubicValuation B) ∧
    ¬ IsEffectValuation producedCubicValuation ∧
    ¬ ∃ ρ : Matrix (Fin 2) (Fin 2) ℂ, IsState ρ ∧
      ∀ E, IsEffect E → producedCubicValuation E = (bornWeight ρ E).re := by
  refine ⟨committed_coexistent_sum_iff, committed_coexistent_sums_nonvacuous,
    bornRunValuation_isEffectValuation, bornRunValuation_matches,
    fun v hv h => additive_valuation_matching_counts_is_born hv h,
    fun ρ hρ h => born_state_matching_counts_unique hρ h,
    producedCubicValuation_matches,
    fun E hE => ⟨producedCubicValuation_nonneg hE,
      producedCubicValuation_le_one hE⟩,
    producedCubicValuation_one,
    fun c c' i i' hco =>
      producedCubicValuation_additive_on_committed_sums hco,
    producedCubicValuation_not_additive,
    producedCubicValuation_not_effectValuation,
    producedCubicValuation_not_born⟩

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms complexifyRealMatrix_trace
#print axioms recordProjector_trace
#print axioms conjProjector_trace
#print axioms sourcePhaseLift_trace
#print axioms committedContextEffect_trace
#print axioms attainedModel_effect_trace
#print axioms effect_sum_of_trace_one
#print axioms committed_coexistent_sum_iff
#print axioms committed_coexistent_sums_nonvacuous
#print axioms committed_frequency_born
#print axioms bornRunValuation_isEffectValuation
#print axioms bornRunValuation_matches
#print axioms born_state_matching_counts_unique
#print axioms additive_valuation_matching_counts_is_born
#print axioms committed_sums_carry_only_normalisation
#print axioms bornRunValuation_coords
#print axioms producedCubicValuation_compl
#print axioms producedCubicValuation_one
#print axioms producedCubicValuation_zero
#print axioms cubicDeformation_eq_zero
#print axioms committed_effect_im
#print axioms cubicDeformation_committed
#print axioms producedCubicValuation_committed
#print axioms producedCubicValuation_matches
#print axioms producedCubicValuation_normalised
#print axioms producedCubicValuation_additive_on_committed_sums
#print axioms effect_entry_bounds
#print axioms producedCubicValuation_nonneg
#print axioms producedCubicValuation_le_one
#print axioms rhoY_add
#print axioms halfWitness_isEffect
#print axioms witnessSum_isEffect
#print axioms halfWitness_coords
#print axioms producedCubicValuation_halfWitness
#print axioms witnessSum_coords
#print axioms producedCubicValuation_witnessSum
#print axioms producedCubicValuation_not_additive
#print axioms producedCubicValuation_not_effectValuation
#print axioms producedCubicValuation_not_born
#print axioms producedCubicValuation_deviates
#print axioms operationalAdditivityBoundary_receipt

end EventAlgebra
