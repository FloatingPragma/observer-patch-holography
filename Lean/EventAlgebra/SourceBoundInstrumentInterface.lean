import EventAlgebra.LuedersPhaseInstrument
import EventAlgebra.OperationalAdditivityBoundary

set_option autoImplicit false

namespace EventAlgebra

/-!
# Typed source-binding interface for the declared Lüders phase instrument

Register row PR-65 states the source-production obligation of the committed
phase packet: one source construction supplies the common preparation used
across the phase and real contexts, selects the phase-sensitive effect or
instrument, produces the public outcomes through that instrument, and binds
preparation, context, operation, readback, and receipt custody to one run.
Register row PR-64 carries the clause that outcome readback agrees with the
recorded public outcome, and register row PR-03 carries cross-context
valuation additivity.  This module types the binding that a source-bound run
must supply against the committed objects of
`EventAlgebra.LuedersPhaseInstrument` and
`EventAlgebra.OperationalAdditivityBoundary`, proves exactly which fields
the committed corpus determines, and proves one exact statement past the
additivity boundary module.

## What is proved

* Interface.  `SourceBoundDeterminedData Φ` types the fields of a source
  binding that the committed corpus determines for a phase instrument `Φ`:
  a common preparation equal to the committed run-state diagonal with the
  committed count literals `111/179` and `68/179`, a public outcome table
  equal to the committed effect table read on that preparation, and a
  readback clause identifying the instrument's induced effects on the
  preparation with the public table.  `SourceBoundInstrumentBinding Φ`
  extends it by the custody data of PR-65: an abstract custody digest and
  one binding proposition per context and outcome, required to hold.
* Attainability of the determined part.  Every phase instrument inhabits
  the determined interface (`determinedDataOfInstrument`,
  `sourceBoundDeterminedData_nonempty`); the explicit partial inhabitant
  for the committed Lüders instrument is
  `determinedDataOfInstrument luedersPhaseInstrument`.  For every
  determined datum of every instrument, the public table is the fixture
  frequency table (`SourceBoundDeterminedData.publicTable_matches`), the
  instrument's outcome traces on the preparation equal those frequencies
  (`SourceBoundDeterminedData.readback_matches`), the phase entry is `1/2`
  (`SourceBoundDeterminedData.publicTable_phase`), and the eight
  enumerated outcome-`0` entries are the literal table with values
  `111/179`, `315/716`, and `1/2`
  (`SourceBoundDeterminedData.publicTable_literal`).
* Externality of the custody data.  The digest is a free parameter over
  the same determined part (`binding_digest_free_parameter`,
  `bindingOfInstrumentDigest_injective`); two bindings agreeing on every
  determined field and differing at the digest exist
  (`committed_corpus_does_not_determine_binding`); and any two binding
  inhabitants have propositionally equal binding assignments
  (`SourceBoundInstrumentBinding.outcomeBinding_eq`), so with
  propositional extensionality the required-to-hold Prop clauses carry no
  separation and the digest is the only datum that distinguishes
  bindings.  The committed corpus therefore determines every field of the
  interface except the custody data, whose value is external input.
* Non-forcing.  The committed swap-twisted comparison instrument inhabits
  the same determined fields with the same preparation and the same
  public table (`determined_interface_does_not_pin_instrument`), and full
  bindings with identical determined data and identical digest exist for
  both instruments (`full_binding_does_not_pin_instrument`); by the
  imported theorems the two instruments are distinct and exactly one of
  them is repeatable, so the interface's determined part does not pin the
  instrument and only the repeatability clause separates the two.
* The additivity rung.  `affineRunValuation t` is the affine
  interpolation `(1 - t) * bornRunValuation + t * producedCubicValuation`
  between the committed run Born valuation and the committed cubic
  countermodel of the additivity boundary module.  For every real `t` it
  reproduces every fixture frequency (`affineRunValuation_matches`),
  sends the sure effect to one (`affineRunValuation_map_one`), obeys the
  complement rule (`affineRunValuation_compl`), and is additive on every
  coexistent sum formed within the committed effect set
  (`affineRunValuation_additive_on_committed_sums`).  The additivity
  residue at the witness pair of the boundary module is the exact linear
  form `-(3/256) * t` (`affineRunValuation_gap`), the family is injective
  in `t` (`affineRunValuation_injective`), and a member is an additive
  effect valuation, equivalently the Born functional of some state,
  precisely when `t = 0` (`affineRunValuation_isEffectValuation_iff`,
  `affineRunValuation_born_iff`).  Between additivity on the committed
  coexistent sums, which holds at every `t`, and the full PR-03
  additivity premise, which holds only at `t = 0`, lies a one-parameter
  continuum of valuations satisfying every committed clause.
* Composed receipt.  `sourceBoundInterface_receipt` derives the
  conjunction of the clauses above.

## What is not proved here

No run, no export, no producer, and no custody artifact exists in the
corpus, and none is constructed.  The interface types what a source-bound
run must supply, and the externality theorems prove that the custody
fields are freely stipulable from inside the corpus: an inhabitant of the
full binding certifies no source production, no provenance, and no
custody.  Register rows PR-65 and PR-03 are open: nothing here supplies a
source construction, a recorded public outcome, a run binding, or an
operational derivation of cross-context additivity, and the full PR-03
premise is consumed only inside the hypothesis of the `t = 0`
characterizations.  A recorded public outcome itself is not typed
anywhere in the corpus; the readback clause identifies induced effects
with a frequency table, which is the strongest identification the
committed objects support.  Instrument selection stays a declared choice:
the determined interface is inhabited by every declared instrument
compatible with the committed table.  The simulator-side export contract
is a separate design document and nothing here consumes it.

## Premise consumption per register row

* PR-02 (consumed): the run state is read on the committed algebra-state
  surface.
* PR-04 (consumed through the imported modules): the phase slot of the
  committed effect table is the declared Pauli +Y projector.
* PR-64 (consumed as the type of the instrument parameter; open): the
  channel clauses enter through `PhaseInstrument`; the readback clause is
  typed against the frequency table only, and no recorded public outcome
  exists.
* PR-65 (subject of the interface; open): the clauses of the row are
  mirrored as fields; the determined fields are attained, and the custody
  fields are proved external to the committed corpus.
* PR-03 (subject of the rung; open): consumed only as the hypothesis
  `IsEffectValuation` in the `t = 0` characterizations; the affine family
  proves that the committed clauses hold on a continuum on which the full
  premise fails everywhere except at the Born member.

## Falsifier

Any failure of the exact identities: a determined public-table entry
differing from the fixture frequency, an enumerated entry differing from
the literal eight-entry table, an additivity residue at the witness pair
differing from `-(3/256) * t` for some `t`, a member of the affine family
with `t` nonzero that is an additive effect valuation or the Born
functional of a state, or two bindings with distinct digests that cannot
be constructed over the same determined part.

## Nonclaims

The binding structures are typed obligations, not certificates: their
inhabitants here are committed-data instances with stipulated custody
fields, and their existence discharges no register row.
-/

open Matrix

noncomputable section

/-! ## The committed effect table of the static inhabitant -/

/-- The effect table of the static fixture inhabitant is the committed
effect pair table, by definitional unfolding of the fixture assembly. -/
theorem attainedModel_effect_eq_pair :
    attainedModel.effect = committedEffectPair :=
  rfl

/-- The phase-context fixture frequency in closed form: the committed
count literals `(179, 179)` give exactly `1/2`. -/
theorem modelFrequency_phase_zero :
    modelFrequency InstrumentContext.phase 0 = 1 / 2 := by
  have h := committed_frequency_born InstrumentContext.phase 0
  rw [← h]
  exact run_phase_bornWeight

/-! ## The determined part of the source-binding interface -/

/-- **The determined part of a source binding.**  The fields of register
row PR-65 that the committed corpus determines, typed against the
committed objects for a phase instrument `Φ`.  The instrument parameter
itself stands for the PR-65 clause "selects the phase-sensitive effect or
instrument"; `determined_interface_does_not_pin_instrument` proves that
the determined fields do not force that selection. -/
structure SourceBoundDeterminedData (Φ : PhaseInstrument) where
  /-- PR-65 clause "one source construction supplies the common
  preparation used across the phase and real contexts": the common
  preparation.  Only the matrix is typed; the source construction is not
  typed in the corpus. -/
  prep : Matrix (Fin 2) (Fin 2) ℂ
  /-- The common preparation is the committed run-state diagonal with the
  committed count literals `111/179` and `68/179` (register row PR-02). -/
  prep_eq_run_state : prep = Matrix.diagonal ![(111 / 179 : ℂ), (68 / 179 : ℂ)]
  /-- PR-65 clause "produces the public outcomes through that
  instrument": the public outcome table, one exact value per context and
  outcome.  Only the frequency surface is typed; a recorded public
  outcome is not typed in the corpus. -/
  publicTable : InstrumentContext → Fin 2 → ℂ
  /-- The public outcome table is the committed effect table read on the
  preparation. -/
  publicTable_eq_born : ∀ (c : InstrumentContext) (i : Fin 2),
    publicTable c i = bornWeight prep (committedEffectPair c i)
  /-- PR-64 clause "outcome readback agrees with the recorded public
  outcome", at the frequency surface: the instrument's induced effects on
  the preparation are identified with the public table. -/
  readback : ∀ (c : InstrumentContext) (i : Fin 2),
    (Φ.outcomeMap c i prep).trace = publicTable c i

/-- The preparation of every determined datum is the committed run
state. -/
theorem SourceBoundDeterminedData.prep_committed {Φ : PhaseInstrument}
    (D : SourceBoundDeterminedData Φ) : D.prep = committedRunState :=
  D.prep_eq_run_state

/-- The public table of every determined datum, of every instrument, is
the fixture frequency table. -/
theorem SourceBoundDeterminedData.publicTable_matches {Φ : PhaseInstrument}
    (D : SourceBoundDeterminedData Φ) (c : InstrumentContext) (i : Fin 2) :
    D.publicTable c i = modelFrequency c i := by
  rw [D.publicTable_eq_born, D.prep_committed]
  have h := committed_frequency_born c i
  rw [attainedModel_effect_eq_pair] at h
  exact h

/-- Every instrument carrying a determined datum has its outcome traces
on the committed preparation pinned to the fixture frequencies, in both
outcomes of every context. -/
theorem SourceBoundDeterminedData.readback_matches {Φ : PhaseInstrument}
    (D : SourceBoundDeterminedData Φ) (c : InstrumentContext) (i : Fin 2) :
    (Φ.outcomeMap c i D.prep).trace = modelFrequency c i := by
  rw [D.readback, D.publicTable_matches]

/-- The phase entry of every determined public table is exactly `1/2`. -/
theorem SourceBoundDeterminedData.publicTable_phase {Φ : PhaseInstrument}
    (D : SourceBoundDeterminedData Φ) :
    D.publicTable InstrumentContext.phase 0 = 1 / 2 := by
  rw [D.publicTable_matches, modelFrequency_phase_zero]

/-- The eight enumerated outcome-`0` entries of every determined public
table are the literal table: `111/179` in the diagonal and
record-conjugate contexts, `315/716` in the rotated contexts, `1/2` in
the phase context. -/
theorem SourceBoundDeterminedData.publicTable_literal {Φ : PhaseInstrument}
    (D : SourceBoundDeterminedData Φ) (k : Fin 8) :
    D.publicTable (enumerateContexts k) 0 =
      ![(111 / 179 : ℂ), 111 / 179, 111 / 179, 315 / 716, 315 / 716,
        315 / 716, 315 / 716, 1 / 2] k := by
  have h := luedersPhaseInstrument_run_table k
  rw [luedersPhaseInstrument.inducedEffect] at h
  rw [D.publicTable_eq_born, D.prep_committed]
  exact h

/-- **The determined part is attained by every phase instrument.**  The
committed run state, the Born values of the committed effect table on it,
and the instrument's induced-effect clause inhabit every determined
field.  At `Φ = luedersPhaseInstrument` this is the explicit partial
inhabitant of the binding interface realized by the committed corpus. -/
def determinedDataOfInstrument (Φ : PhaseInstrument) :
    SourceBoundDeterminedData Φ where
  prep := committedRunState
  prep_eq_run_state := rfl
  publicTable := fun c i => bornWeight committedRunState (committedEffectPair c i)
  publicTable_eq_born := fun _ _ => rfl
  readback := fun c i => Φ.inducedEffect c i committedRunState

/-- Nonvacuity of the determined interface, uniformly in the
instrument. -/
theorem sourceBoundDeterminedData_nonempty (Φ : PhaseInstrument) :
    Nonempty (SourceBoundDeterminedData Φ) :=
  ⟨determinedDataOfInstrument Φ⟩

/-! ## The full binding interface and the externality of its custody data -/

/-- **A source-bound instrument binding.**  The determined data extended
by the custody clauses of register row PR-65: "binds preparation,
context, operation, readback, and receipt custody to one run".  The
custody fields are typed obligations whose content the committed corpus
does not determine (`committed_corpus_does_not_determine_binding`). -/
structure SourceBoundInstrumentBinding (Φ : PhaseInstrument) extends
    SourceBoundDeterminedData Φ where
  /-- PR-65 clause "receipt custody": an abstract digest standing for a
  custody artifact of one run.  No committed object constrains its
  value. -/
  custodyDigest : ℕ
  /-- PR-65 clause "binds preparation, context, operation, readback, and
  receipt custody to one run": one proposition per context and outcome
  asserting that the public outcome is bound to the preparation through
  the instrument.  The content of these propositions is external
  input. -/
  outcomeBinding : InstrumentContext → Fin 2 → Prop
  /-- The binding propositions hold.  Under propositional extensionality
  any two assignments satisfying this field are equal
  (`SourceBoundInstrumentBinding.outcomeBinding_eq`), so the digest
  carries all separation between bindings. -/
  outcomeBinding_holds : ∀ (c : InstrumentContext) (i : Fin 2),
    outcomeBinding c i

/-- The binding with a stipulated digest over the committed determined
data.  The digest and the trivially true binding propositions are
external stipulations; the definition certifies no custody. -/
def bindingOfInstrumentDigest (Φ : PhaseInstrument) (d : ℕ) :
    SourceBoundInstrumentBinding Φ where
  toSourceBoundDeterminedData := determinedDataOfInstrument Φ
  custodyDigest := d
  outcomeBinding := fun _ _ => True
  outcomeBinding_holds := fun _ _ => True.intro

/-- Nonvacuity of the full binding interface, uniformly in the
instrument. -/
theorem sourceBoundInstrumentBinding_nonempty (Φ : PhaseInstrument) :
    Nonempty (SourceBoundInstrumentBinding Φ) :=
  ⟨bindingOfInstrumentDigest Φ 0⟩

/-- **The digest is a free parameter.**  Every natural number is realized
as the custody digest of a binding whose determined part is the committed
one. -/
theorem binding_digest_free_parameter (Φ : PhaseInstrument) (d : ℕ) :
    (bindingOfInstrumentDigest Φ d).toSourceBoundDeterminedData =
        determinedDataOfInstrument Φ ∧
      (bindingOfInstrumentDigest Φ d).custodyDigest = d :=
  ⟨rfl, rfl⟩

/-- Distinct digests give distinct bindings: the stipulated-digest family
is injective. -/
theorem bindingOfInstrumentDigest_injective (Φ : PhaseInstrument) :
    Function.Injective (bindingOfInstrumentDigest Φ) := fun _ _ h =>
  congrArg (fun B : SourceBoundInstrumentBinding Φ => B.custodyDigest) h

/-- **The committed corpus does not determine the source binding.**  Two
bindings of the committed Lüders instrument agree on every determined
field and differ at the custody digest; the custody data is external
input, not derivable from the committed corpus. -/
theorem committed_corpus_does_not_determine_binding :
    ∃ B₁ B₂ : SourceBoundInstrumentBinding luedersPhaseInstrument,
      B₁.toSourceBoundDeterminedData = B₂.toSourceBoundDeterminedData ∧
        B₁.custodyDigest ≠ B₂.custodyDigest ∧ B₁ ≠ B₂ := by
  refine ⟨bindingOfInstrumentDigest luedersPhaseInstrument 0,
    bindingOfInstrumentDigest luedersPhaseInstrument 1, rfl, ?_, ?_⟩
  · exact Nat.zero_ne_one
  · intro h
    exact Nat.zero_ne_one (congrArg
      (fun B : SourceBoundInstrumentBinding luedersPhaseInstrument =>
        B.custodyDigest) h)

/-- **The required-to-hold binding propositions carry no separation.**
Any two bindings of the same instrument have equal binding assignments:
each binding proposition holds in both, so propositional extensionality
identifies them.  The custody digest is therefore the only field that can
distinguish two bindings over the same determined data. -/
theorem SourceBoundInstrumentBinding.outcomeBinding_eq {Φ : PhaseInstrument}
    (B₁ B₂ : SourceBoundInstrumentBinding Φ) :
    B₁.outcomeBinding = B₂.outcomeBinding := by
  funext c i
  exact propext ⟨fun _ => B₂.outcomeBinding_holds c i,
    fun _ => B₁.outcomeBinding_holds c i⟩

/-! ## The determined interface does not pin the instrument -/

/-- **Non-forcing at the determined level.**  The committed Lüders
instrument and the committed swap-twisted comparison instrument inhabit
the determined interface with the same preparation and the same public
table, the two instruments are distinct, and exactly one of them is
repeatable: the determined part of the source-binding interface does not
pin the instrument, and only the repeatability clause separates the
two. -/
theorem determined_interface_does_not_pin_instrument :
    ∃ (D : SourceBoundDeterminedData luedersPhaseInstrument)
      (D' : SourceBoundDeterminedData swapTwistedPhaseInstrument),
      D.prep = D'.prep ∧ D.publicTable = D'.publicTable ∧
        swapTwistedPhaseInstrument ≠ luedersPhaseInstrument ∧
        luedersPhaseInstrument.Repeatable ∧
        ¬ swapTwistedPhaseInstrument.Repeatable :=
  ⟨determinedDataOfInstrument luedersPhaseInstrument,
    determinedDataOfInstrument swapTwistedPhaseInstrument, rfl, rfl,
    swapTwistedPhaseInstrument_ne_lueders, luedersPhaseInstrument_repeatable,
    swapTwistedPhaseInstrument_not_repeatable⟩

/-- **Non-forcing at the full-binding level.**  Full bindings with the
same preparation, the same public table, and the same custody digest
exist for the two distinct instruments; even the complete interface with
stipulated custody data does not pin the instrument. -/
theorem full_binding_does_not_pin_instrument :
    ∃ (B : SourceBoundInstrumentBinding luedersPhaseInstrument)
      (B' : SourceBoundInstrumentBinding swapTwistedPhaseInstrument),
      B.prep = B'.prep ∧ B.publicTable = B'.publicTable ∧
        B.custodyDigest = B'.custodyDigest ∧
        luedersPhaseInstrument.Repeatable ∧
        ¬ swapTwistedPhaseInstrument.Repeatable :=
  ⟨bindingOfInstrumentDigest luedersPhaseInstrument 0,
    bindingOfInstrumentDigest swapTwistedPhaseInstrument 0, rfl, rfl, rfl,
    luedersPhaseInstrument_repeatable,
    swapTwistedPhaseInstrument_not_repeatable⟩

/-! ## The additivity rung: the affine family between Born and the
countermodel -/

/-- Unfolding of the run Born valuation of the additivity boundary
module. -/
theorem bornRunValuation_def (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornRunValuation E = (bornWeight committedRunState E).re :=
  rfl

/-- The complement rule of the run Born valuation, on every matrix. -/
theorem bornRunValuation_compl (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornRunValuation (1 - E) = 1 - bornRunValuation E := by
  unfold bornRunValuation
  rw [bornWeight_sub, bornWeight_one committedRunState_isState,
    Complex.sub_re, Complex.one_re]

/-- The run Born value of the half witness is exactly `1/4`. -/
theorem bornRunValuation_halfWitness : bornRunValuation halfWitness = 1 / 4 := by
  obtain ⟨h00, h11, -⟩ := halfWitness_coords
  rw [bornRunValuation_coords, h00, h11]
  norm_num

/-- The run Born value of the witness sum is exactly `1/2`. -/
theorem bornRunValuation_witnessSum : bornRunValuation witnessSum = 1 / 2 := by
  obtain ⟨h00, h11, -⟩ := witnessSum_coords
  rw [bornRunValuation_coords, h00, h11]
  norm_num

/-- **The affine family.**  The interpolation between the committed run
Born valuation and the committed cubic countermodel of the additivity
boundary module, with real parameter `t`. -/
def affineRunValuation (t : ℝ) (E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (1 - t) * bornRunValuation E + t * producedCubicValuation E

/-- The `t = 0` member is the run Born valuation. -/
theorem affineRunValuation_at_zero : affineRunValuation 0 = bornRunValuation := by
  funext E
  unfold affineRunValuation
  ring

/-- The `t = 1` member is the committed cubic countermodel. -/
theorem affineRunValuation_at_one :
    affineRunValuation 1 = producedCubicValuation := by
  funext E
  unfold affineRunValuation
  ring

/-- Every member reproduces every fixture frequency on the committed
effect set. -/
theorem affineRunValuation_matches (t : ℝ) (c : InstrumentContext) (i : Fin 2) :
    ((affineRunValuation t (attainedModel.effect c i) : ℝ) : ℂ) =
      modelFrequency c i := by
  have h : affineRunValuation t (attainedModel.effect c i) =
      bornRunValuation (attainedModel.effect c i) := by
    unfold affineRunValuation
    rw [producedCubicValuation_committed]
    ring
  rw [h]
  exact bornRunValuation_matches c i

/-- Every member sends the sure effect to one. -/
theorem affineRunValuation_map_one (t : ℝ) :
    affineRunValuation t (1 : Matrix (Fin 2) (Fin 2) ℂ) = 1 := by
  unfold affineRunValuation
  rw [producedCubicValuation_one, bornRunValuation_isEffectValuation.map_one]
  ring

/-- Every member obeys the complement rule, on every matrix. -/
theorem affineRunValuation_compl (t : ℝ) (E : Matrix (Fin 2) (Fin 2) ℂ) :
    affineRunValuation t (1 - E) = 1 - affineRunValuation t E := by
  unfold affineRunValuation
  rw [bornRunValuation_compl, producedCubicValuation_compl]
  ring

/-- Every member is additive on every coexistent sum formed within the
committed effect set. -/
theorem affineRunValuation_additive_on_committed_sums (t : ℝ)
    {c c' : InstrumentContext} {i i' : Fin 2}
    (hco : IsEffect (attainedModel.effect c i + attainedModel.effect c' i')) :
    affineRunValuation t (attainedModel.effect c i + attainedModel.effect c' i') =
      affineRunValuation t (attainedModel.effect c i) +
        affineRunValuation t (attainedModel.effect c' i') := by
  have hb := bornRunValuation_isEffectValuation.additive
    (attainedModel.effect_isEffect c i) (attainedModel.effect_isEffect c' i') hco
  have hp := producedCubicValuation_additive_on_committed_sums hco
  unfold affineRunValuation
  rw [hb, hp]
  ring

/-- The member value at the half witness in closed form. -/
theorem affineRunValuation_halfWitness (t : ℝ) :
    affineRunValuation t halfWitness = 1 / 4 + (15 / 512) * t := by
  unfold affineRunValuation
  rw [bornRunValuation_halfWitness, producedCubicValuation_halfWitness]
  ring

/-- The member value at the witness sum in closed form. -/
theorem affineRunValuation_witnessSum (t : ℝ) :
    affineRunValuation t witnessSum = 1 / 2 + (3 / 64) * t := by
  unfold affineRunValuation
  rw [bornRunValuation_witnessSum, producedCubicValuation_witnessSum]
  ring

/-- **The exact additivity residue.**  At the witness pair of the
additivity boundary module the failure of additivity is the linear form
`-(3/256) * t`: it vanishes exactly at the Born member. -/
theorem affineRunValuation_gap (t : ℝ) :
    affineRunValuation t witnessSum -
        (affineRunValuation t halfWitness + affineRunValuation t halfWitness) =
      -(3 / 256) * t := by
  rw [affineRunValuation_witnessSum, affineRunValuation_halfWitness]
  ring

/-- The family is injective: distinct parameters give distinct
valuations, separated at the witness sum. -/
theorem affineRunValuation_injective : Function.Injective affineRunValuation := by
  intro t₁ t₂ h
  have h' := congrFun h witnessSum
  rw [affineRunValuation_witnessSum, affineRunValuation_witnessSum] at h'
  linarith

/-- **Full additivity singles out the Born member.**  A member of the
affine family is an additive effect valuation, the full PR-03 additivity
premise, precisely when `t = 0`. -/
theorem affineRunValuation_isEffectValuation_iff (t : ℝ) :
    IsEffectValuation (affineRunValuation t) ↔ t = 0 := by
  constructor
  · intro hv
    have h := hv.additive halfWitness_isEffect halfWitness_isEffect
      witnessSum_isEffect
    have h' : affineRunValuation t witnessSum =
        affineRunValuation t halfWitness + affineRunValuation t halfWitness := h
    rw [affineRunValuation_witnessSum, affineRunValuation_halfWitness] at h'
    linarith
  · rintro rfl
    rw [affineRunValuation_at_zero]
    exact bornRunValuation_isEffectValuation

/-- **Born representability singles out the same member.**  A member of
the affine family is the Born functional of a state on all effects
precisely when `t = 0`, and the representing state at `t = 0` is the
committed run state. -/
theorem affineRunValuation_born_iff (t : ℝ) :
    (∃ ρ : Matrix (Fin 2) (Fin 2) ℂ, IsState ρ ∧
        ∀ E, IsEffect E → affineRunValuation t E = (bornWeight ρ E).re) ↔
      t = 0 := by
  constructor
  · rintro ⟨ρ, -, hrep⟩
    have h1 := hrep halfWitness halfWitness_isEffect
    have h2 := hrep witnessSum witnessSum_isEffect
    have hadd : bornWeight ρ witnessSum =
        bornWeight ρ halfWitness + bornWeight ρ halfWitness :=
      bornWeight_add ρ halfWitness halfWitness
    have hre := congrArg Complex.re hadd
    rw [Complex.add_re, ← h1, ← h2] at hre
    rw [affineRunValuation_witnessSum, affineRunValuation_halfWitness] at hre
    linarith
  · rintro rfl
    refine ⟨committedRunState, committedRunState_isState, fun E _ => ?_⟩
    rw [affineRunValuation_at_zero, bornRunValuation_def]

/-! ## The composed receipt -/

/-- **The source-binding interface receipt.**  One conjunction:

1. every phase instrument inhabits the determined interface;
2. every determined datum of every instrument has the fixture frequency
   table as its public table and as its readback traces;
3. the custody digest is a free parameter over the committed determined
   part;
4. two bindings of the committed Lüders instrument agree on every
   determined field and are distinct: the committed corpus does not
   determine the binding;
5. the swap-twisted instrument inhabits the same determined fields with
   the same preparation and public table, and exactly one of the two
   instruments is repeatable;
6. a member of the affine family between the run Born valuation and the
   cubic countermodel is an additive effect valuation precisely at
   `t = 0`, with additivity residue exactly `-(3/256) * t` at the
   witness pair.

Boundary: the receipt certifies typing and delimitation only.  It
certifies no source production, no recorded public outcome, no
provenance, no custody, and no operational additivity, and it discharges
no register row. -/
theorem sourceBoundInterface_receipt :
    (∀ Φ : PhaseInstrument, Nonempty (SourceBoundDeterminedData Φ)) ∧
    (∀ (Φ : PhaseInstrument) (D : SourceBoundDeterminedData Φ)
        (c : InstrumentContext) (i : Fin 2),
      D.publicTable c i = modelFrequency c i ∧
        (Φ.outcomeMap c i D.prep).trace = modelFrequency c i) ∧
    (∀ (Φ : PhaseInstrument) (d : ℕ),
      (bindingOfInstrumentDigest Φ d).toSourceBoundDeterminedData =
          determinedDataOfInstrument Φ ∧
        (bindingOfInstrumentDigest Φ d).custodyDigest = d) ∧
    (∃ B₁ B₂ : SourceBoundInstrumentBinding luedersPhaseInstrument,
      B₁.toSourceBoundDeterminedData = B₂.toSourceBoundDeterminedData ∧
        B₁ ≠ B₂) ∧
    (∃ (D : SourceBoundDeterminedData luedersPhaseInstrument)
        (D' : SourceBoundDeterminedData swapTwistedPhaseInstrument),
      D.prep = D'.prep ∧ D.publicTable = D'.publicTable) ∧
    luedersPhaseInstrument.Repeatable ∧
    ¬ swapTwistedPhaseInstrument.Repeatable ∧
    (∀ t : ℝ, IsEffectValuation (affineRunValuation t) ↔ t = 0) ∧
    (∀ t : ℝ,
      affineRunValuation t witnessSum -
          (affineRunValuation t halfWitness + affineRunValuation t halfWitness) =
        -(3 / 256) * t) := by
  obtain ⟨B₁, B₂, hBeq, -, hBne⟩ := committed_corpus_does_not_determine_binding
  exact ⟨sourceBoundDeterminedData_nonempty,
    fun Φ D c i => ⟨D.publicTable_matches c i, D.readback_matches c i⟩,
    binding_digest_free_parameter,
    ⟨B₁, B₂, hBeq, hBne⟩,
    ⟨determinedDataOfInstrument luedersPhaseInstrument,
      determinedDataOfInstrument swapTwistedPhaseInstrument, rfl, rfl⟩,
    luedersPhaseInstrument_repeatable,
    swapTwistedPhaseInstrument_not_repeatable,
    affineRunValuation_isEffectValuation_iff,
    affineRunValuation_gap⟩

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms attainedModel_effect_eq_pair
#print axioms modelFrequency_phase_zero
#print axioms SourceBoundDeterminedData.prep_committed
#print axioms SourceBoundDeterminedData.publicTable_matches
#print axioms SourceBoundDeterminedData.readback_matches
#print axioms SourceBoundDeterminedData.publicTable_phase
#print axioms SourceBoundDeterminedData.publicTable_literal
#print axioms determinedDataOfInstrument
#print axioms sourceBoundDeterminedData_nonempty
#print axioms bindingOfInstrumentDigest
#print axioms sourceBoundInstrumentBinding_nonempty
#print axioms binding_digest_free_parameter
#print axioms bindingOfInstrumentDigest_injective
#print axioms committed_corpus_does_not_determine_binding
#print axioms SourceBoundInstrumentBinding.outcomeBinding_eq
#print axioms determined_interface_does_not_pin_instrument
#print axioms full_binding_does_not_pin_instrument
#print axioms bornRunValuation_def
#print axioms bornRunValuation_compl
#print axioms bornRunValuation_halfWitness
#print axioms bornRunValuation_witnessSum
#print axioms affineRunValuation
#print axioms affineRunValuation_at_zero
#print axioms affineRunValuation_at_one
#print axioms affineRunValuation_matches
#print axioms affineRunValuation_map_one
#print axioms affineRunValuation_compl
#print axioms affineRunValuation_additive_on_committed_sums
#print axioms affineRunValuation_halfWitness
#print axioms affineRunValuation_witnessSum
#print axioms affineRunValuation_gap
#print axioms affineRunValuation_injective
#print axioms affineRunValuation_isEffectValuation_iff
#print axioms affineRunValuation_born_iff
#print axioms sourceBoundInterface_receipt

end EventAlgebra
