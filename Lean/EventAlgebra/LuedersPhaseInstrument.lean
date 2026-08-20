import EventAlgebra.OperationalPhaseAttainment
import Dynamics.ChoiCPTP

set_option autoImplicit false

namespace EventAlgebra

/-!
# A declared Lüders instrument for the committed phase-effect contexts

Register row PR-64 states the channel structure that a phase-sensitive
measurement carries: context- and outcome-indexed completely-positive
trace-nonincreasing state-update maps, a trace-preserving summed channel in
each context, induced effects agreeing with the declared effect table, and
outcome readback agreeing with the recorded public outcome.  The committed
phase packet (`EventAlgebra.OperationalPhaseInstrument`,
`EventAlgebra.OperationalPhaseAttainment`) types a declared effect and a
static Born-table fixture and carries none of those maps.  This module
constructs the Lüders instrument family for the committed eight-context
effect table as exact finite matrix algebra, proves the channel clauses of
PR-64 for it, joins it to the committed run fixture, and proves that the
declared effect table does not determine the instrument.

## What is proved

* General dimension.  For a matrix `P`, `luedersOutcomeMap P` is the linear
  map `X ↦ P * X * Pᴴ`.  For an event `P` (Hermitian idempotent): each
  outcome map is completely positive by a singleton Kraus family
  (`luedersOutcomeMap_completelyPositive`); its trace on any input is the
  Born weight of `P` (`trace_luedersOutcomeMap`), which is the induced-effect
  clause; on positive-semidefinite inputs the trace does not increase
  (`trace_luedersOutcomeMap_le`); the summed channel
  `luedersChannel P = luedersOutcomeMap P + luedersOutcomeMap (1 - P)` is
  trace preserving on every matrix and completely positive, hence CPTP
  (`luedersChannel_isCPTP`); the outcome map on a state of nonzero Born
  weight is the Born weight times the committed Lüders update
  (`luedersOutcomeMap_eq_bornWeight_smul_luedersUpdate`), so the normalized
  post-measurement state is certain of `P`
  (`luedersOutcomeMap_normalized_mem_certainStates`); composition of outcome
  maps is the outcome map of the product (`luedersOutcomeMap_comp`), giving
  idempotence (`luedersOutcomeMap_comp_self`) and outcome exclusivity
  (`luedersOutcomeMap_compl_comp`, `luedersOutcomeMap_comp_compl`).
* The committed contexts.  `PhaseInstrument` is a structure on the eight
  contexts of `InstrumentContext` with outcome-indexed linear maps and the
  PR-64 channel clauses as fields: complete positivity per outcome, trace
  non-increase on positive-semidefinite inputs, trace preservation of the
  summed channel per context, and induced effects equal to the committed
  table `committedEffectPair` (index `0` carries the context effect
  `committedContextEffect c`, index `1` its complement, matching the
  committed indexing).  `PhaseInstrument.Repeatable` is repeatability, the
  normalized post-measurement state being certain of its effect; it is a
  state-update proxy for the readback clause, which is open.
  `luedersPhaseInstrument` inhabits the structure
  with the Lüders maps of the committed effects
  (`committedEffectPair_isEvent` supplies the projector certificates; the
  phase context's effect is `sourcePhaseLift`, whose event certificate is
  `sourcePhaseLift_isEvent`), and it is repeatable
  (`luedersPhaseInstrument_repeatable`); `RepeatablePhaseInstrument`
  bundles the readback field, inhabited by
  `luedersRepeatablePhaseInstrument`.
* Join to the fixture.  On `committedRunState`, the outcome-`0` trace of the
  Lüders instrument in every context equals the fixture frequency of the
  committed static inhabitant (`luedersPhaseInstrument_run_frequency`), with
  the closed forms `111/179` in the diagonal context, `1/2` in the phase
  context, and the full eight-entry table
  (`luedersPhaseInstrument_run_table`).  The closed forms `111/179` and
  `1/2` are also obtained by direct matrix computation through the outcome
  maps (`luedersPhaseInstrument_run_diagonal_direct`,
  `luedersPhaseInstrument_run_phase_direct`), with no use of the fixture's
  Born-fit field; the phase outcome on the run state is half the declared
  projector (`luedersPhaseInstrument_run_phase_output`).  The phase-context
  post-measurement state of the run state is the declared effect itself
  (`luedersUpdate_run_phase`).
* Non-uniqueness.  `swapTwistedPhaseInstrument` conjugates each Lüders
  outcome by the swap unitary.  It satisfies every channel clause
  (`swapTwistedPhaseInstrument : PhaseInstrument`) with the same induced
  effects, it differs from the Lüders instrument on the committed run state
  at an explicit matrix entry (`swapTwisted_ne_lueders_on_run`), and it
  fails the repeatability clause
  (`swapTwistedPhaseInstrument_not_repeatable`).  Hence the declared effect
  table does not determine the instrument
  (`effect_table_does_not_determine_instrument`); the Lüders choice is a
  declared selection.
* Composed receipt.  `luedersPhaseInstrument_receipt` derives the
  conjunction of the clauses above from one typed antecedent bundle
  `PhaseInstrumentDecisionBundle` (the declared context-effect table and
  the committed run state); `phaseInstrument_nonempty` and
  `repeatablePhaseInstrument_nonempty` record nonvacuity.

## What is declared

The effect table `committedContextEffect` (PR-04 for the phase slot, the
committed web projectors for the others), the run state
`committedRunState` (PR-02), and the Lüders selection itself.  The swap
twist is a constructed comparison instrument and carries no selection
status.  Nothing here is source-produced, measured, sampled, validated, or
derived from a run.

## Premise consumption per register row

* PR-04 (consumed, axiomatize disposition): the phase context's effect is
  the declared Pauli +Y projector `sourcePhaseLift`.  The instrument is
  built on top of that declaration and derives nothing about it.
* PR-02 (consumed): the run state is read on the committed algebra-state
  surface.
* PR-64 (channel clauses constructed; readback and selection open): the
  CP, trace-nonincreasing, trace-preserving-sum, and induced-effect clauses
  are constructed and proved for a declared conditional instrument.  The
  clause "outcome readback agrees with the recorded public outcome" is
  open: no public outcome exists in this module, and the repeatability
  predicate is its state-update proxy only.  The selection of the Lüders
  instrument among the instruments compatible with the effect table is a
  declaration, as `effect_table_does_not_determine_instrument` shows.
* PR-65 (open): no source preparation, public outcome, provenance, or
  custody is typed or proved.
* PR-03 (open): the per-context trace-preserving sums give normalization in
  each context only; cross-context valuation additivity is untouched.
* PR-52 (open): physical attachment is untouched.

## What is NOT proved

No measurement, run, sampling, or validation occurs.  The module does not
prove that the Lüders instrument is implemented by any source, that the
committed counts were produced through it, or that any public outcome is
read back through it.  It does not derive the effect table, the run state,
or the selection of Lüders over other compatible instruments.  Uniqueness
of the repeatable instrument compatible with the effect table is not
proved.

## Hypotheses, falsifier, nonclaims

Hypotheses: the committed definitions of `committedContextEffect`,
`committedEffectPair`, `committedRunState`, and the static inhabitant
`attainedModel` with its fixture literals; the CP predicates of
`Dynamics.ChoiCPTP`.

Falsifier: any failure of the exact identities against the committed
objects, for instance a context effect that is not a Hermitian idempotent,
a summed channel whose trace differs from the input trace on some matrix,
an outcome-`0` trace on the run state differing from the fixture
frequency, or equality of the Lüders and swap-twisted outputs on the run
state.

Nonclaims: the instrument is a declared conditional construction; no
emergent or laboratory instrument is claimed, and its existence discharges
no register row.

## Tagging convention

As in `EventAlgebra.Basic`: **algebra-only** statements consume no trace
pairing; **trace-dependent** statements pass through `Tr(ρ M)`.
-/

open Matrix
open OPH.QFT
open OPH.Dynamics (CMat IsCompletelyPositive IsPositiveMap IsTracePreserving
  IsCPTP isCompletelyPositive_of_kraus)
open scoped ComplexOrder

noncomputable section

variable {n : ℕ}

/-! ## Lüders outcome maps in general dimension -/

/-- **Algebra-only.**  The Lüders outcome map of a matrix `P`: the linear map
`X ↦ P * X * Pᴴ`.  For a Hermitian `P` this is the compression `X ↦ P X P`
whose normalization is the committed `luedersUpdate`. -/
def luedersOutcomeMap (P : CMat n) : CMat n →ₗ[ℂ] CMat n where
  toFun X := P * X * Pᴴ
  map_add' X Y := by rw [mul_add, add_mul]
  map_smul' c X := by rw [RingHom.id_apply, mul_smul_comm, smul_mul_assoc]

@[simp]
theorem luedersOutcomeMap_apply (P X : CMat n) :
    luedersOutcomeMap P X = P * X * Pᴴ :=
  rfl

/-- **Algebra-only.**  On a Hermitian matrix the outcome map is the
two-sided compression. -/
theorem luedersOutcomeMap_apply_of_isHermitian {P : CMat n}
    (hP : P.IsHermitian) (X : CMat n) :
    luedersOutcomeMap P X = P * X * P := by
  rw [luedersOutcomeMap_apply, hP.eq]

/-- **Algebra-only.**  The summed Lüders channel of `P`: the sum of the
outcome maps of `P` and of its complement `1 - P`. -/
def luedersChannel (P : CMat n) : CMat n →ₗ[ℂ] CMat n :=
  luedersOutcomeMap P + luedersOutcomeMap (1 - P)

theorem luedersChannel_apply (P X : CMat n) :
    luedersChannel P X = P * X * Pᴴ + (1 - P) * X * (1 - P)ᴴ :=
  rfl

/-- **Algebra-only.**  Every Lüders outcome map is completely positive: the
singleton family `{P}` is a Kraus family for it. -/
theorem luedersOutcomeMap_completelyPositive (P : CMat n) :
    IsCompletelyPositive (luedersOutcomeMap P) :=
  isCompletelyPositive_of_kraus (luedersOutcomeMap P) (fun _ : Fin 1 => P)
    (fun X => by rw [Fin.sum_univ_one, luedersOutcomeMap_apply])

/-- **Algebra-only.**  Every Lüders outcome map is a positive map. -/
theorem luedersOutcomeMap_posSemidef (P : CMat n) {X : CMat n}
    (hX : X.PosSemidef) : (luedersOutcomeMap P X).PosSemidef :=
  hX.mul_mul_conjTranspose_same P

/-- **Trace-dependent.**  **Induced effect.**  For an event `P`, the trace
of the outcome map on any matrix is the Born weight of `P`: the effect
induced by the Lüders outcome map is `P` itself. -/
theorem trace_luedersOutcomeMap {P : CMat n} (hP : IsEvent P) (X : CMat n) :
    (luedersOutcomeMap P X).trace = bornWeight X P := by
  rw [luedersOutcomeMap_apply_of_isHermitian hP.1, trace_sandwich hP.2]

/-- **Trace-dependent.**  The induced effect of the complement outcome is
the complement event. -/
theorem trace_luedersOutcomeMap_compl {P : CMat n} (hP : IsEvent P)
    (X : CMat n) :
    (luedersOutcomeMap (1 - P) X).trace = bornWeight X (1 - P) :=
  trace_luedersOutcomeMap hP.compl X

/-- **Trace-dependent.**  The Born weights of an event and of its complement
sum to the trace of the input, on every matrix. -/
theorem bornWeight_add_compl (X P : CMat n) :
    bornWeight X P + bornWeight X (1 - P) = X.trace := by
  rw [← bornWeight_add, add_sub_cancel, bornWeight, mul_one]

/-- **Trace-dependent.**  **Trace non-increase.**  On a positive-semidefinite
input the outcome map of an event does not increase the trace, in the
partial order of `ℂ`. -/
theorem trace_luedersOutcomeMap_le {P : CMat n} (hP : IsEvent P) {X : CMat n}
    (hX : X.PosSemidef) : (luedersOutcomeMap P X).trace ≤ X.trace := by
  rw [trace_luedersOutcomeMap hP, ← bornWeight_add_compl X P]
  exact le_add_of_nonneg_right (bornWeight_nonneg hX hP.compl)

/-- **Trace-dependent.**  On a state the outcome trace is at most one. -/
theorem trace_luedersOutcomeMap_le_one {P : CMat n} (hP : IsEvent P)
    {ρ : CMat n} (hρ : IsState ρ) : (luedersOutcomeMap P ρ).trace ≤ 1 := by
  rw [← hρ.2]
  exact trace_luedersOutcomeMap_le hP hρ.1

/-- **Trace-dependent.**  **Trace preservation of the summed channel**, on
every matrix. -/
theorem luedersChannel_tracePreserving {P : CMat n} (hP : IsEvent P) :
    IsTracePreserving (luedersChannel P) := by
  intro X
  rw [luedersChannel, LinearMap.add_apply, trace_add, trace_luedersOutcomeMap hP,
    trace_luedersOutcomeMap_compl hP, bornWeight_add_compl]

/-- **Algebra-only.**  The summed channel is completely positive. -/
theorem luedersChannel_completelyPositive (P : CMat n) :
    IsCompletelyPositive (luedersChannel P) :=
  (luedersOutcomeMap_completelyPositive P).add
    (luedersOutcomeMap_completelyPositive (1 - P))

/-- **Trace-dependent.**  The summed Lüders channel of an event is CPTP. -/
theorem luedersChannel_isCPTP {P : CMat n} (hP : IsEvent P) :
    IsCPTP (luedersChannel P) :=
  ⟨luedersChannel_completelyPositive P, luedersChannel_tracePreserving hP⟩

/-- **Trace-dependent.**  **Compatibility with the committed update.**  On a
matrix of nonzero Born weight the outcome map is the Born weight times the
committed Lüders update. -/
theorem luedersOutcomeMap_eq_bornWeight_smul_luedersUpdate {ρ P : CMat n}
    (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    luedersOutcomeMap P ρ = bornWeight ρ P • luedersUpdate ρ P := by
  rw [luedersOutcomeMap_apply_of_isHermitian hP.1, luedersUpdate, smul_smul,
    mul_inv_cancel₀ hw, one_smul]

/-- **Trace-dependent.**  The normalized outcome is the committed Lüders
update. -/
theorem luedersOutcomeMap_normalized {ρ P : CMat n} (hP : IsEvent P) :
    (bornWeight ρ P)⁻¹ • luedersOutcomeMap P ρ = luedersUpdate ρ P := by
  rw [luedersOutcomeMap_apply_of_isHermitian hP.1, luedersUpdate]

/-- **Trace-dependent.**  **Repeatability.**  The normalized post-measurement
state of a state with nonzero Born weight is certain of the event. -/
theorem luedersOutcomeMap_normalized_mem_certainStates {ρ P : CMat n}
    (hρ : IsState ρ) (hP : IsEvent P) (hw : bornWeight ρ P ≠ 0) :
    (bornWeight ρ P)⁻¹ • luedersOutcomeMap P ρ ∈ certainStates P := by
  rw [luedersOutcomeMap_normalized hP]
  exact luedersUpdate_mem_certainStates hρ hP hw

/-- **Algebra-only.**  Composition of outcome maps is the outcome map of the
product. -/
theorem luedersOutcomeMap_comp (P Q : CMat n) :
    luedersOutcomeMap Q ∘ₗ luedersOutcomeMap P = luedersOutcomeMap (Q * P) := by
  apply LinearMap.ext
  intro X
  simp only [LinearMap.comp_apply, luedersOutcomeMap_apply, conjTranspose_mul,
    mul_assoc]

/-- **Algebra-only.**  The outcome map of the zero matrix is the zero map. -/
theorem luedersOutcomeMap_zero : luedersOutcomeMap (0 : CMat n) = 0 := by
  apply LinearMap.ext
  intro X
  simp

/-- **Algebra-only.**  **Idempotence.**  Applying the outcome map of an event
twice equals applying it once. -/
theorem luedersOutcomeMap_comp_self {P : CMat n} (hP : IsEvent P) :
    luedersOutcomeMap P ∘ₗ luedersOutcomeMap P = luedersOutcomeMap P := by
  rw [luedersOutcomeMap_comp, hP.2]

/-- **Algebra-only.**  The complement of an event annihilates it on either
side. -/
theorem compl_mul_self_eq_zero {P : CMat n} (hP : IsEvent P) :
    (1 - P) * P = 0 := by
  rw [sub_mul, one_mul, hP.2, sub_self]

theorem self_mul_compl_eq_zero {P : CMat n} (hP : IsEvent P) :
    P * (1 - P) = 0 := by
  rw [mul_sub, mul_one, hP.2, sub_self]

/-- **Algebra-only.**  **Outcome exclusivity.**  The complement outcome after
the event outcome is the zero map. -/
theorem luedersOutcomeMap_compl_comp {P : CMat n} (hP : IsEvent P) :
    luedersOutcomeMap (1 - P) ∘ₗ luedersOutcomeMap P = 0 := by
  rw [luedersOutcomeMap_comp, compl_mul_self_eq_zero hP, luedersOutcomeMap_zero]

/-- **Algebra-only.**  The event outcome after the complement outcome is the
zero map. -/
theorem luedersOutcomeMap_comp_compl {P : CMat n} (hP : IsEvent P) :
    luedersOutcomeMap P ∘ₗ luedersOutcomeMap (1 - P) = 0 := by
  rw [luedersOutcomeMap_comp, self_mul_compl_eq_zero hP, luedersOutcomeMap_zero]

/-! ## Unitary twists of outcome maps -/

/-- **Trace-dependent.**  Twisting an outcome map by a left isometry
(`Uᴴ * U = 1`) leaves the induced effect unchanged. -/
theorem trace_luedersOutcomeMap_mul_of_isometry {U P : CMat n}
    (hU : Uᴴ * U = 1) (hP : IsEvent P) (X : CMat n) :
    (luedersOutcomeMap (U * P) X).trace = bornWeight X P := by
  rw [luedersOutcomeMap_apply, conjTranspose_mul, hP.1.eq, trace_mul_cycle]
  have hcollapse : P * Uᴴ * (U * P) = P := by
    rw [mul_assoc, ← mul_assoc Uᴴ, hU, one_mul, hP.2]
  rw [hcollapse, bornWeight, trace_mul_comm]

/-! ## The committed effect table as projectors -/

theorem committedEffectPair_zero (c : InstrumentContext) :
    committedEffectPair c 0 = committedContextEffect c :=
  rfl

theorem committedEffectPair_one (c : InstrumentContext) :
    committedEffectPair c 1 = 1 - committedContextEffect c :=
  rfl

/-- **Algebra-only.**  Every entry of the committed effect table is a
projection event: the web entries through `webContextProjector_isEvent`,
the phase entry through `sourcePhaseLift_isEvent`, and the complements
through `IsEvent.compl`. -/
theorem committedEffectPair_isEvent (c : InstrumentContext) (i : Fin 2) :
    IsEvent (committedEffectPair c i) := by
  fin_cases i
  · exact committedContextEffect_isEvent c
  · exact (committedContextEffect_isEvent c).compl

/-- The phase-slot effect of the committed table is the declared lift, and
it is a Hermitian idempotent with the explicit entries of the Pauli +Y
projector. -/
theorem committedContextEffect_phase :
    committedContextEffect InstrumentContext.phase = sourcePhaseLift :=
  rfl

theorem sourcePhaseLift_isHermitian : sourcePhaseLift.IsHermitian :=
  sourcePhaseLift_isEvent.1

theorem sourcePhaseLift_idem : sourcePhaseLift * sourcePhaseLift = sourcePhaseLift :=
  sourcePhaseLift_isEvent.2

theorem sourcePhaseLift_entries :
    sourcePhaseLift = !![(1 / 2 : ℂ), -(Complex.I / 2); Complex.I / 2, 1 / 2] :=
  sourcePhaseLift_eq_rhoYPlus

/-! ## The phase-instrument structure on the committed contexts -/

/-- **A phase instrument on the committed context family.**  Context- and
outcome-indexed linear maps on the committed two-dimensional algebra,
carrying the channel clauses of register row PR-64 as fields: complete
positivity of each outcome map, trace non-increase on positive-semidefinite
inputs, trace preservation of the summed channel in each context, and
induced effects equal to the committed effect table `committedEffectPair`
(index `0` is the context effect, index `1` its complement).  The structure
carries no readback clause; repeatability, a state-update proxy for it, is
the separate predicate `PhaseInstrument.Repeatable`, bundled as a field in
`RepeatablePhaseInstrument`.  The structure types no source, public
outcome, provenance, or custody. -/
structure PhaseInstrument where
  /-- The outcome maps, one linear map per context and outcome. -/
  outcomeMap : InstrumentContext → Fin 2 → (CMat 2 →ₗ[ℂ] CMat 2)
  /-- Complete positivity of every outcome map. -/
  completelyPositive : ∀ (c : InstrumentContext) (i : Fin 2),
    IsCompletelyPositive (outcomeMap c i)
  /-- Trace non-increase of every outcome map on positive-semidefinite
  inputs. -/
  traceNonincreasing : ∀ (c : InstrumentContext) (i : Fin 2) (X : CMat 2),
    X.PosSemidef → (outcomeMap c i X).trace ≤ X.trace
  /-- Trace preservation of the summed channel in every context. -/
  summed_tracePreserving : ∀ c : InstrumentContext,
    IsTracePreserving (outcomeMap c 0 + outcomeMap c 1)
  /-- The induced effect of every outcome is the committed table entry. -/
  inducedEffect : ∀ (c : InstrumentContext) (i : Fin 2) (X : CMat 2),
    (outcomeMap c i X).trace = bornWeight X (committedEffectPair c i)

namespace PhaseInstrument

/-- The summed channel of a context. -/
def summedChannel (Φ : PhaseInstrument) (c : InstrumentContext) :
    CMat 2 →ₗ[ℂ] CMat 2 :=
  Φ.outcomeMap c 0 + Φ.outcomeMap c 1

/-- **Trace-dependent.**  The summed channel of every context is CPTP. -/
theorem summedChannel_isCPTP (Φ : PhaseInstrument) (c : InstrumentContext) :
    IsCPTP (Φ.summedChannel c) :=
  ⟨(Φ.completelyPositive c 0).add (Φ.completelyPositive c 1),
    Φ.summed_tracePreserving c⟩

/-- Every outcome map is a positive map (complete positivity specialized to
the trivial amplification). -/
theorem outcomeMap_posSemidef (Φ : PhaseInstrument) (c : InstrumentContext)
    (i : Fin 2) {X : CMat 2} (hX : X.PosSemidef) :
    (Φ.outcomeMap c i X).PosSemidef :=
  (Φ.completelyPositive c i).isPositiveMap X hX

/-- **Trace-dependent.**  The two outcome traces of a state sum to one. -/
theorem outcome_traces_sum (Φ : PhaseInstrument) (c : InstrumentContext)
    {ρ : CMat 2} (hρ : IsState ρ) :
    (Φ.outcomeMap c 0 ρ).trace + (Φ.outcomeMap c 1 ρ).trace = 1 := by
  have h := Φ.summed_tracePreserving c ρ
  rw [LinearMap.add_apply, trace_add, hρ.2] at h
  exact h

/-- **Repeatability, the state-update proxy for the readback clause.**  For
every context, outcome, and state of nonzero Born weight on the outcome's
effect, the normalized post-measurement state is certain of that effect.
Agreement with a recorded public outcome, the PR-64 readback clause itself,
is open. -/
def Repeatable (Φ : PhaseInstrument) : Prop :=
  ∀ (c : InstrumentContext) (i : Fin 2) (ρ : CMat 2), IsState ρ →
    bornWeight ρ (committedEffectPair c i) ≠ 0 →
      (bornWeight ρ (committedEffectPair c i))⁻¹ • Φ.outcomeMap c i ρ ∈
        certainStates (committedEffectPair c i)

/-- Two phase instruments induce the same effects on every input: the
induced-effect field pins the trace of every outcome map. -/
theorem inducedEffect_eq (Φ Ψ : PhaseInstrument) (c : InstrumentContext)
    (i : Fin 2) (X : CMat 2) :
    (Φ.outcomeMap c i X).trace = (Ψ.outcomeMap c i X).trace := by
  rw [Φ.inducedEffect, Ψ.inducedEffect]

end PhaseInstrument

/-- **A repeatable phase instrument.**  A phase instrument bundled with the
repeatability predicate as a field.  The Lüders instrument inhabits it
(`luedersRepeatablePhaseInstrument`); the swap-twisted comparison
instrument below satisfies the channel clauses and fails this field. -/
structure RepeatablePhaseInstrument extends PhaseInstrument where
  /-- The normalized post-measurement state of every outcome is certain of
  its effect. -/
  repeatable : toPhaseInstrument.Repeatable

/-! ## The Lüders instrument of the committed effect table -/

/-- **The declared Lüders phase instrument.**  Each context and outcome
carries the Lüders outcome map of its committed table entry.  The channel
clauses are the general-dimension theorems specialized through
`committedEffectPair_isEvent`.  This is a declared conditional
construction; no source implements it here. -/
def luedersPhaseInstrument : PhaseInstrument where
  outcomeMap c i := luedersOutcomeMap (committedEffectPair c i)
  completelyPositive c i := luedersOutcomeMap_completelyPositive _
  traceNonincreasing c i X hX :=
    trace_luedersOutcomeMap_le (committedEffectPair_isEvent c i) hX
  summed_tracePreserving c := by
    intro X
    rw [LinearMap.add_apply, trace_add,
      trace_luedersOutcomeMap (committedEffectPair_isEvent c 0),
      trace_luedersOutcomeMap (committedEffectPair_isEvent c 1),
      committedEffectPair_zero, committedEffectPair_one, bornWeight_add_compl]
  inducedEffect c i X := trace_luedersOutcomeMap (committedEffectPair_isEvent c i) X

@[simp]
theorem luedersPhaseInstrument_outcomeMap (c : InstrumentContext) (i : Fin 2) :
    luedersPhaseInstrument.outcomeMap c i =
      luedersOutcomeMap (committedEffectPair c i) :=
  rfl

/-- The summed channel of every context of the Lüders instrument is the
general `luedersChannel` of the context effect. -/
theorem luedersPhaseInstrument_summedChannel (c : InstrumentContext) :
    luedersPhaseInstrument.summedChannel c =
      luedersChannel (committedContextEffect c) :=
  rfl

/-- **Trace-dependent.**  The Lüders instrument is repeatable. -/
theorem luedersPhaseInstrument_repeatable : luedersPhaseInstrument.Repeatable :=
  fun c i _ hρ hw =>
    luedersOutcomeMap_normalized_mem_certainStates hρ
      (committedEffectPair_isEvent c i) hw

/-- The Lüders instrument bundled with its repeatability field. -/
def luedersRepeatablePhaseInstrument : RepeatablePhaseInstrument where
  toPhaseInstrument := luedersPhaseInstrument
  repeatable := luedersPhaseInstrument_repeatable

theorem luedersRepeatablePhaseInstrument_toPhaseInstrument :
    luedersRepeatablePhaseInstrument.toPhaseInstrument = luedersPhaseInstrument :=
  rfl

/-- **Algebra-only.**  Idempotence and exclusivity of the Lüders instrument
in every context. -/
theorem luedersPhaseInstrument_idem_exclusive (c : InstrumentContext) :
    (luedersPhaseInstrument.outcomeMap c 0 ∘ₗ luedersPhaseInstrument.outcomeMap c 0 =
        luedersPhaseInstrument.outcomeMap c 0) ∧
      (luedersPhaseInstrument.outcomeMap c 1 ∘ₗ luedersPhaseInstrument.outcomeMap c 1 =
        luedersPhaseInstrument.outcomeMap c 1) ∧
      luedersPhaseInstrument.outcomeMap c 1 ∘ₗ luedersPhaseInstrument.outcomeMap c 0 = 0 ∧
      luedersPhaseInstrument.outcomeMap c 0 ∘ₗ luedersPhaseInstrument.outcomeMap c 1 = 0 := by
  have hE := committedContextEffect_isEvent c
  simp only [luedersPhaseInstrument_outcomeMap, committedEffectPair_zero,
    committedEffectPair_one]
  exact ⟨luedersOutcomeMap_comp_self hE, luedersOutcomeMap_comp_self hE.compl,
    luedersOutcomeMap_compl_comp hE, luedersOutcomeMap_comp_compl hE⟩

/-- The channel clauses along the committed enumeration of the eight
contexts. -/
theorem luedersPhaseInstrument_enumerated (k : Fin 8) :
    IsCPTP (luedersPhaseInstrument.summedChannel (enumerateContexts k)) :=
  luedersPhaseInstrument.summedChannel_isCPTP _

/-- Nonvacuity of the structure. -/
theorem phaseInstrument_nonempty : Nonempty PhaseInstrument :=
  ⟨luedersPhaseInstrument⟩

/-- Nonvacuity of the bundled repeatable structure. -/
theorem repeatablePhaseInstrument_nonempty : Nonempty RepeatablePhaseInstrument :=
  ⟨luedersRepeatablePhaseInstrument⟩

/-! ## Join to the committed run fixture -/

/-- The outcome-`0` effect of the static inhabitant is the committed table
entry. -/
theorem attainedModel_effect_zero (c : InstrumentContext) :
    attainedModel.effect c 0 = committedEffectPair c 0 :=
  rfl

/-- **Trace-dependent.**  **Fixture join.**  In every context the outcome-`0`
trace of the Lüders instrument on the committed run state is the fixture
frequency of the static inhabitant.  This reuses the committed
`attained_born_consistency`; no literal is recomputed. -/
theorem luedersPhaseInstrument_run_frequency (c : InstrumentContext) :
    (luedersPhaseInstrument.outcomeMap c 0 committedRunState).trace =
      binaryFrequency (attainedModel.counts c) := by
  rw [luedersPhaseInstrument.inducedEffect, ← attainedModel_effect_zero,
    ← attainedModel_prep_eq]
  exact attained_born_consistency c

/-- The phase-context outcome-`0` trace on the run state is exactly `1/2`. -/
theorem luedersPhaseInstrument_run_phase :
    (luedersPhaseInstrument.outcomeMap InstrumentContext.phase 0
      committedRunState).trace = 1 / 2 := by
  rw [luedersPhaseInstrument_run_frequency, attained_phase_frequency]

/-- The diagonal-context outcome-`0` trace on the run state is `111/179`. -/
theorem luedersPhaseInstrument_run_diagonal :
    (luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
      committedRunState).trace = 111 / 179 := by
  rw [luedersPhaseInstrument_run_frequency, attainedModel_diagonal_counts,
    binaryFrequency_diagonal_run]

/-- The rotated fixture frequency in closed form. -/
theorem binaryFrequency_rotated_run : binaryFrequency (315, 401) = 315 / 716 := by
  norm_num [binaryFrequency]

/-- The full eight-entry table of outcome-`0` traces on the run state along
the committed enumeration: `111/179` in the diagonal and record-conjugate
contexts, `315/716` in the rotated contexts, `1/2` in the phase context. -/
theorem luedersPhaseInstrument_run_table (k : Fin 8) :
    (luedersPhaseInstrument.outcomeMap (enumerateContexts k) 0
      committedRunState).trace =
      ![(111 / 179 : ℂ), 111 / 179, 111 / 179, 315 / 716, 315 / 716, 315 / 716,
        315 / 716, 1 / 2] k := by
  fin_cases k
  · exact luedersPhaseInstrument_run_diagonal
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_diagonal_run
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_diagonal_run
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_rotated_run
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_rotated_run
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_rotated_run
  · rw [luedersPhaseInstrument_run_frequency]
    exact binaryFrequency_rotated_run
  · exact luedersPhaseInstrument_run_phase

/-- The committed run state as a literal matrix. -/
theorem committedRunState_eq_literal :
    committedRunState = !![(111 / 179 : ℂ), 0; 0, 68 / 179] := by
  rw [committedRunState, Matrix.diagonal_fin_two]
  rfl

/-- The diagonal-context effect of the committed table as a literal complex
matrix. -/
theorem committedContextEffect_diagonal_eq :
    committedContextEffect (InstrumentContext.web WebContext.diagonal) =
      !![(1 : ℂ), 0; 0, 0] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [committedContextEffect, webContextProjector, complexifyRealMatrix,
      recordProjector]

/-- **Trace-dependent.**  The phase-context Born weight of the run state is
exactly `1/2`, read from the committed fixture theorems. -/
theorem run_phase_bornWeight :
    bornWeight committedRunState sourcePhaseLift = 1 / 2 := by
  have h := attained_born_consistency InstrumentContext.phase
  rw [attained_phase_frequency] at h
  exact h

/-- **Trace-dependent.**  The phase outcome of the run state has a
normalized post-measurement state. -/
theorem run_phase_bornWeight_ne_zero :
    bornWeight committedRunState sourcePhaseLift ≠ 0 := by
  rw [run_phase_bornWeight]
  norm_num

/-- **Trace-dependent.**  The phase-context post-measurement state of the
run state is the declared effect itself: the Lüders update of the
record-diagonal run state by the Pauli +Y projector is that projector. -/
theorem luedersUpdate_run_phase :
    luedersUpdate committedRunState sourcePhaseLift = sourcePhaseLift := by
  rw [luedersUpdate, run_phase_bornWeight, sourcePhaseLift_eq_rhoYPlus,
    committedRunState_eq_literal]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [rhoYPlus, Matrix.mul_apply, Fin.sum_univ_two] <;>
    ring_nf <;> simp [Complex.I_sq] <;> norm_num

/-- The diagonal-context outcome-`0` output on the run state in closed form:
the record weight on the record projector. -/
theorem luedersPhaseInstrument_run_diagonal_output :
    luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
        committedRunState =
      !![(111 / 179 : ℂ), 0; 0, 0] := by
  rw [luedersPhaseInstrument_outcomeMap, committedEffectPair_zero,
    luedersOutcomeMap_apply, committedContextEffect_diagonal_eq,
    committedRunState_eq_literal]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Matrix.mul_apply, Fin.sum_univ_two, Matrix.conjTranspose_apply]

/-- **Trace-dependent.**  The diagonal-context outcome-`0` trace on the run
state, obtained by direct matrix computation through the outcome map with
no use of the fixture's Born-fit field. -/
theorem luedersPhaseInstrument_run_diagonal_direct :
    (luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
      committedRunState).trace = 111 / 179 := by
  rw [luedersPhaseInstrument_run_diagonal_output]
  norm_num [Matrix.trace_fin_two]

/-- The phase-context outcome-`0` output on the run state in closed form:
half the declared projector, by direct matrix computation. -/
theorem luedersPhaseInstrument_run_phase_output :
    luedersPhaseInstrument.outcomeMap InstrumentContext.phase 0 committedRunState =
      (1 / 2 : ℂ) • sourcePhaseLift := by
  rw [luedersPhaseInstrument_outcomeMap, committedEffectPair_zero,
    committedContextEffect_phase,
    luedersOutcomeMap_apply_of_isHermitian sourcePhaseLift_isHermitian,
    sourcePhaseLift_eq_rhoYPlus, committedRunState_eq_literal]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [rhoYPlus, Matrix.mul_apply, Fin.sum_univ_two] <;>
    ring_nf <;> simp [Complex.I_sq] <;> norm_num

/-- **Trace-dependent.**  The phase-context outcome-`0` trace on the run
state is `1/2` by direct matrix computation through the outcome map, with no
use of the fixture's Born-fit field; it agrees with the fixture route
`luedersPhaseInstrument_run_phase`. -/
theorem luedersPhaseInstrument_run_phase_direct :
    (luedersPhaseInstrument.outcomeMap InstrumentContext.phase 0
      committedRunState).trace = 1 / 2 := by
  rw [luedersPhaseInstrument_run_phase_output, trace_smul, sourcePhaseLift_eq_rhoYPlus]
  norm_num [rhoYPlus, Matrix.trace_fin_two, smul_eq_mul]

/-! ## Non-uniqueness: a swap-twisted instrument with the same effects -/

/-- The swap unitary on the committed two-dimensional algebra. -/
def swapUnitary : CMat 2 := !![0, 1; 1, 0]

theorem swapUnitary_conjTranspose_mul : swapUnitaryᴴ * swapUnitary = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [swapUnitary, Matrix.mul_apply, Fin.sum_univ_two]

/-- **The swap-twisted phase instrument.**  Each outcome map is the Lüders
outcome map followed by conjugation with the swap unitary, written as the
outcome map of the Kraus operator `swapUnitary * E`.  It is a constructed
comparison instrument: it satisfies every channel clause with the same
induced effects as the Lüders instrument and carries no selection status. -/
def swapTwistedPhaseInstrument : PhaseInstrument where
  outcomeMap c i := luedersOutcomeMap (swapUnitary * committedEffectPair c i)
  completelyPositive c i := luedersOutcomeMap_completelyPositive _
  traceNonincreasing c i X hX := by
    rw [trace_luedersOutcomeMap_mul_of_isometry swapUnitary_conjTranspose_mul
      (committedEffectPair_isEvent c i), ← bornWeight_add_compl X]
    exact le_add_of_nonneg_right
      (bornWeight_nonneg hX (committedEffectPair_isEvent c i).compl)
  summed_tracePreserving c := by
    intro X
    rw [LinearMap.add_apply, trace_add,
      trace_luedersOutcomeMap_mul_of_isometry swapUnitary_conjTranspose_mul
        (committedEffectPair_isEvent c 0),
      trace_luedersOutcomeMap_mul_of_isometry swapUnitary_conjTranspose_mul
        (committedEffectPair_isEvent c 1),
      committedEffectPair_zero, committedEffectPair_one, bornWeight_add_compl]
  inducedEffect c i X :=
    trace_luedersOutcomeMap_mul_of_isometry swapUnitary_conjTranspose_mul
      (committedEffectPair_isEvent c i) X

@[simp]
theorem swapTwistedPhaseInstrument_outcomeMap (c : InstrumentContext) (i : Fin 2) :
    swapTwistedPhaseInstrument.outcomeMap c i =
      luedersOutcomeMap (swapUnitary * committedEffectPair c i) :=
  rfl

/-- The twisted outcome map is the Lüders outcome map followed by the swap
conjugation. -/
theorem swapTwistedPhaseInstrument_outcomeMap_comp (c : InstrumentContext)
    (i : Fin 2) :
    swapTwistedPhaseInstrument.outcomeMap c i =
      luedersOutcomeMap swapUnitary ∘ₗ luedersPhaseInstrument.outcomeMap c i := by
  rw [swapTwistedPhaseInstrument_outcomeMap, luedersPhaseInstrument_outcomeMap,
    luedersOutcomeMap_comp]

/-- The diagonal-context outcome-`0` output of the twisted instrument on the
run state: the record weight moved onto the companion projector. -/
theorem swapTwistedPhaseInstrument_run_diagonal_output :
    swapTwistedPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
        committedRunState =
      !![0, 0; 0, (111 / 179 : ℂ)] := by
  rw [swapTwistedPhaseInstrument_outcomeMap, committedEffectPair_zero,
    luedersOutcomeMap_apply, committedContextEffect_diagonal_eq,
    committedRunState_eq_literal, swapUnitary]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Matrix.mul_apply, Fin.sum_univ_two, Matrix.conjTranspose_apply]

/-- **Explicit entry inequality.**  The two instruments differ on the run
state in the diagonal context at the record entry: `111/179` against `0`. -/
theorem swapTwisted_ne_lueders_on_run :
    swapTwistedPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
        committedRunState ≠
      luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
        committedRunState := by
  rw [swapTwistedPhaseInstrument_run_diagonal_output,
    luedersPhaseInstrument_run_diagonal_output]
  intro h
  have h00 := congrFun (congrFun h 0) 0
  norm_num at h00

/-- The two instruments are distinct. -/
theorem swapTwistedPhaseInstrument_ne_lueders :
    swapTwistedPhaseInstrument ≠ luedersPhaseInstrument := by
  intro h
  apply swapTwisted_ne_lueders_on_run
  rw [h]

/-- **Trace-dependent.**  The companion-weighted output has vanishing Born
weight on the record projector. -/
theorem bornWeight_companion_output_record :
    bornWeight (!![0, 0; 0, (111 / 179 : ℂ)]) (!![(1 : ℂ), 0; 0, 0]) = 0 := by
  norm_num [bornWeight, Matrix.trace_fin_two, Matrix.mul_apply, Fin.sum_univ_two]

/-- **Trace-dependent.**  The twisted instrument fails the repeatability
clause: in the diagonal context the normalized post-measurement state of
the run state is the companion projector, whose Born weight on the record
projector is `0`. -/
theorem swapTwistedPhaseInstrument_not_repeatable :
    ¬ swapTwistedPhaseInstrument.Repeatable := by
  intro h
  have hw : bornWeight committedRunState
      (committedEffectPair (InstrumentContext.web WebContext.diagonal) 0) = 111 / 179 := by
    rw [← luedersPhaseInstrument.inducedEffect]
    exact luedersPhaseInstrument_run_diagonal
  have hmem := h (InstrumentContext.web WebContext.diagonal) 0 committedRunState
    committedRunState_isState (by rw [hw]; norm_num)
  rw [hw, swapTwistedPhaseInstrument_run_diagonal_output] at hmem
  have h1 := hmem.2
  rw [committedEffectPair_zero, committedContextEffect_diagonal_eq, bornWeight_smul,
    bornWeight_companion_output_record, mul_zero] at h1
  exact zero_ne_one h1

/-- **The declared effect table does not determine the instrument.**  Two
phase instruments with the same induced effects on every input and every
channel clause differ on the committed run state; exactly one of them is
repeatable.  The Lüders instrument is a declared selection among the
instruments compatible with the table. -/
theorem effect_table_does_not_determine_instrument :
    ∃ Φ Ψ : PhaseInstrument,
      (∀ (c : InstrumentContext) (i : Fin 2) (X : CMat 2),
        (Φ.outcomeMap c i X).trace = (Ψ.outcomeMap c i X).trace) ∧
      Φ.Repeatable ∧ ¬ Ψ.Repeatable ∧
      Φ.outcomeMap (InstrumentContext.web WebContext.diagonal) 0 committedRunState ≠
        Ψ.outcomeMap (InstrumentContext.web WebContext.diagonal) 0 committedRunState :=
  ⟨luedersPhaseInstrument, swapTwistedPhaseInstrument,
    PhaseInstrument.inducedEffect_eq _ _, luedersPhaseInstrument_repeatable,
    swapTwistedPhaseInstrument_not_repeatable, swapTwisted_ne_lueders_on_run.symm⟩

/-! ## The antecedent bundle and composed receipt -/

/-- **The antecedent bundle.**  The declared context-effect table (PR-04 in
the phase slot, the committed web projectors elsewhere) and the declared
run state (PR-02).  The bundle carries declared data; it certifies no source
production. -/
structure PhaseInstrumentDecisionBundle where
  /-- The declared per-context effect table. -/
  contextEffect : InstrumentContext → CMat 2
  /-- The table is the committed one. -/
  contextEffect_eq : contextEffect = committedContextEffect
  /-- The declared common preparation. -/
  prep : CMat 2
  /-- The preparation is the committed run state. -/
  prep_eq_run_state : prep = committedRunState

/-- The committed instance of the bundle: nonvacuity of the antecedent. -/
def recordedInstrumentBundle : PhaseInstrumentDecisionBundle where
  contextEffect := committedContextEffect
  contextEffect_eq := rfl
  prep := committedRunState
  prep_eq_run_state := rfl

/-- **Composed receipt.**  From the one antecedent bundle, the conjunction:

1. complete positivity of every Lüders outcome map;
2. trace non-increase on positive-semidefinite inputs;
3. CPTP summed channel in every context;
4. induced effects equal to the bundle's table (outcome `0`) and its
   complement (outcome `1`);
5. repeatability (the normalized post-measurement state is certain of its
   effect);
6. idempotence and exclusivity of the outcome maps in every context;
7. the fixture join on the bundle's preparation: outcome-`0` traces equal
   the static inhabitant's fixture frequencies, with `1/2` in the phase
   context and `111/179` in the diagonal context;
8. the phase post-measurement state of the bundle's preparation is the
   bundle's phase effect;
9. non-uniqueness: a second instrument with the same induced effects that
   is not repeatable and differs on the bundle's preparation;
10. nonvacuity of `PhaseInstrument` and of `RepeatablePhaseInstrument`.

Boundary: the receipt certifies a declared conditional instrument.  It
certifies no source implementation, public outcome, provenance, custody, or
validation, and it discharges no register row. -/
theorem luedersPhaseInstrument_receipt (D : PhaseInstrumentDecisionBundle) :
    (∀ (c : InstrumentContext) (i : Fin 2),
      IsCompletelyPositive (luedersPhaseInstrument.outcomeMap c i)) ∧
    (∀ (c : InstrumentContext) (i : Fin 2) (X : CMat 2), X.PosSemidef →
      (luedersPhaseInstrument.outcomeMap c i X).trace ≤ X.trace) ∧
    (∀ c : InstrumentContext, IsCPTP (luedersPhaseInstrument.summedChannel c)) ∧
    (∀ (c : InstrumentContext) (X : CMat 2),
      (luedersPhaseInstrument.outcomeMap c 0 X).trace =
        bornWeight X (D.contextEffect c) ∧
      (luedersPhaseInstrument.outcomeMap c 1 X).trace =
        bornWeight X (1 - D.contextEffect c)) ∧
    luedersPhaseInstrument.Repeatable ∧
    (∀ c : InstrumentContext,
      (luedersPhaseInstrument.outcomeMap c 0 ∘ₗ luedersPhaseInstrument.outcomeMap c 0 =
          luedersPhaseInstrument.outcomeMap c 0) ∧
        (luedersPhaseInstrument.outcomeMap c 1 ∘ₗ luedersPhaseInstrument.outcomeMap c 1 =
          luedersPhaseInstrument.outcomeMap c 1) ∧
        luedersPhaseInstrument.outcomeMap c 1 ∘ₗ luedersPhaseInstrument.outcomeMap c 0 = 0 ∧
        luedersPhaseInstrument.outcomeMap c 0 ∘ₗ luedersPhaseInstrument.outcomeMap c 1 = 0) ∧
    (∀ c : InstrumentContext,
      (luedersPhaseInstrument.outcomeMap c 0 D.prep).trace =
        binaryFrequency (attainedModel.counts c)) ∧
    (luedersPhaseInstrument.outcomeMap InstrumentContext.phase 0 D.prep).trace = 1 / 2 ∧
    (luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
      D.prep).trace = 111 / 179 ∧
    luedersUpdate D.prep (D.contextEffect InstrumentContext.phase) =
      D.contextEffect InstrumentContext.phase ∧
    (∃ Ψ : PhaseInstrument,
      (∀ (c : InstrumentContext) (i : Fin 2) (X : CMat 2),
        (luedersPhaseInstrument.outcomeMap c i X).trace = (Ψ.outcomeMap c i X).trace) ∧
      ¬ Ψ.Repeatable ∧
      luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
          D.prep ≠
        Ψ.outcomeMap (InstrumentContext.web WebContext.diagonal) 0 D.prep) ∧
    Nonempty PhaseInstrument ∧ Nonempty RepeatablePhaseInstrument := by
  obtain ⟨_, rfl, _, rfl⟩ := D
  refine ⟨luedersPhaseInstrument.completelyPositive,
    luedersPhaseInstrument.traceNonincreasing,
    luedersPhaseInstrument.summedChannel_isCPTP,
    fun c X => ⟨luedersPhaseInstrument.inducedEffect c 0 X,
      luedersPhaseInstrument.inducedEffect c 1 X⟩,
    luedersPhaseInstrument_repeatable,
    luedersPhaseInstrument_idem_exclusive,
    luedersPhaseInstrument_run_frequency,
    luedersPhaseInstrument_run_phase,
    luedersPhaseInstrument_run_diagonal,
    luedersUpdate_run_phase,
    ⟨swapTwistedPhaseInstrument, PhaseInstrument.inducedEffect_eq _ _,
      swapTwistedPhaseInstrument_not_repeatable, swapTwisted_ne_lueders_on_run.symm⟩,
    phaseInstrument_nonempty, repeatablePhaseInstrument_nonempty⟩

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms luedersOutcomeMap_completelyPositive
#print axioms trace_luedersOutcomeMap
#print axioms trace_luedersOutcomeMap_le
#print axioms luedersChannel_isCPTP
#print axioms luedersOutcomeMap_eq_bornWeight_smul_luedersUpdate
#print axioms luedersOutcomeMap_normalized_mem_certainStates
#print axioms luedersOutcomeMap_comp
#print axioms luedersOutcomeMap_comp_self
#print axioms luedersOutcomeMap_compl_comp
#print axioms luedersOutcomeMap_comp_compl
#print axioms trace_luedersOutcomeMap_mul_of_isometry
#print axioms committedEffectPair_isEvent
#print axioms sourcePhaseLift_isHermitian
#print axioms sourcePhaseLift_idem
#print axioms luedersPhaseInstrument
#print axioms luedersPhaseInstrument_repeatable
#print axioms luedersRepeatablePhaseInstrument
#print axioms repeatablePhaseInstrument_nonempty
#print axioms luedersPhaseInstrument_idem_exclusive
#print axioms phaseInstrument_nonempty
#print axioms luedersPhaseInstrument_run_frequency
#print axioms luedersPhaseInstrument_run_phase
#print axioms luedersPhaseInstrument_run_diagonal
#print axioms luedersPhaseInstrument_run_table
#print axioms run_phase_bornWeight
#print axioms luedersUpdate_run_phase
#print axioms luedersPhaseInstrument_run_diagonal_output
#print axioms luedersPhaseInstrument_run_diagonal_direct
#print axioms luedersPhaseInstrument_run_phase_output
#print axioms luedersPhaseInstrument_run_phase_direct
#print axioms swapTwistedPhaseInstrument
#print axioms swapTwistedPhaseInstrument_run_diagonal_output
#print axioms swapTwisted_ne_lueders_on_run
#print axioms swapTwistedPhaseInstrument_ne_lueders
#print axioms swapTwistedPhaseInstrument_not_repeatable
#print axioms effect_table_does_not_determine_instrument
#print axioms recordedInstrumentBundle
#print axioms luedersPhaseInstrument_receipt

end EventAlgebra
