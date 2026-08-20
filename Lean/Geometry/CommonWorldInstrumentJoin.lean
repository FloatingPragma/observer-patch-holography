import Geometry.CommonWorldIslandBridge
import EventAlgebra.LuedersPhaseInstrument
import EventAlgebra.PartitionPinching
import ScaledMaxwellStability
import Geometry.EinsteinBranchInhabitant

set_option autoImplicit false

open scoped BigOperators Matrix

/-!
# CW1 instrument join: four typed joins on the bridged common-world record (issue #740)

WHAT IS PROVED.  The committed bridged witness of
`Geometry/CommonWorldIslandBridge.lean` carries the kinematics island and
the screen island on one record, tied by one metric dictionary.  This
module extends that record by four typed joins, each stated on an object
the record carries, and proves of every inhabitant a composed receipt
containing the ten committed clauses of the bridged receipt and four new
clauses.

The threaded objects, named exactly.  (1) The committed twelve-port,
thirty-seam, twenty-face carrier of the record's unit-step Maxwell bundle
`MaxwellEvolutionBundle` of `Screen/TemporalMaxwellEvolution.lean`, on
which the step-scaled bundle `ScaledMaxwellBundle` of
`Screen/ScaledMaxwellStability.lean` is typed.  (2) The committed
two-dimensional record algebra `Matrix (Fin 2) (Fin 2) ℂ` carrying the
record's partition and state fields, the committed diagonal partition
`diagonalPartition 2` of `QFT/ObserverAccessCut.lean`, the committed
pinching `partitionPinching` of `EventAlgebra/PartitionPinching.lean`, and
the phase instrument `PhaseInstrument` of
`EventAlgebra/LuedersPhaseInstrument.lean` with the committed run state
`committedRunState`.  (3) The record's own mass and frame fields on the
committed Hermitian Lorentz module, through `internalClock` and
`frameWorldline` of `Geometry/InternalClockRestFrequency.lean`.  (4) The
diagonal of the instrument's run state as the reference law of a
`RepairLawData (Fin 2) (Fin 2)`, with a horizon record over it, and the
register constructor `registerOfHorizonRecord` of
`Geometry/EinsteinBranchInhabitant.lean`.

Join 1, certified evolution on the shared carrier.  The extended record
carries one `ScaledMaxwellBundle` whose seam history agrees with the
record's unit-step history at steps zero and one and whose seam current is
the record's seam current; the carrier is shared by type.  For every
inhabitant: the record's unit-step Ampere law is the `h = 1` instance of
the scaled law (`maxwell_is_unit_scaled`, through the committed recovery
theorem); every Courant constant dominates the committed eigenvalues `5`
and `3 + √5` of `CᵀC` (`courantBound_ge_eigenvalue`), so every certified
bundle has step `h < 1` and no certified bundle carries the unit step
(`scaledBundle_step_lt_one`, `scaledBundle_step_ne_one`); the scaled
staggered form is nonnegative at every step, and for zero current the
scaled electric and magnetic energies are bounded at every step by the
committed explicit multiples of the initial form; and the unit-step law on
the same carrier admits the committed unbounded zero-current mode
(`unit_step_instability`).  The record carries the declared unit-step law
and its certified replacement on one carrier, with the contrast as a
clause of the receipt.

Join 2, the instrument island on the shared record algebra.  The
committed witness's partition field is the one-block partition
`oneBlockPartition 2` and its state field is `witnessState`, the pure
state on the record direction; both live on `Matrix (Fin 2) (Fin 2) ℂ`,
the carrier of the instrument.  Two requested identifications are false
and are recorded as obstruction theorems: the record's state is not the
run state (`committedWitness_state_ne_runState`, entry `(0,0)` reads `1`
against `111/179`), and no one-block partition carries the diagonal effect
pair, since every `ProjectivePartition 2 1` has its single projector equal
to `1` (`projectivePartition_one_block_proj`,
`record_partition_proj_ne_record_effect`); the two field values are
recorded definitionally (`committedWitness_partition`,
`committedWitness_state`).  The typed transport that does exist is the
refinement of the record's partition by the committed diagonal partition, a
refinement every partition enjoys because the one-block projector is the
identity, so the record-partition side of the transport carries no content
beyond the shared carrier and the load sits in the pinching identity and
the state clause: the diagonal effect pair of the instrument's diagonal
context is exactly the committed `diagonalPartition 2`
(`committedEffectPair_diagonal`), the instrument's diagonal-context summed
channel is the committed pinching of that partition on every matrix
(`lueders_diagonal_channel_eq_pinching`), the record partition's pinching
is the identity and absorbs the diagonal pinching
(`record_pinching_absorbs_diagonal`), and the record's state and the run
state are fixed points of the instrument's diagonal channel
(`diagonal_pinching_witnessState`, `diagonal_pinching_runState`).  The
extended record carries the pinching identity and the state fixed-point
clause as fields; the swap-twisted instrument of lane A violates the
pinching identity at the run state
(`swapTwisted_diagonal_channel_runState`,
`twistedExtension_breaks_pinching`), so the clause is load-bearing
(`pinching_clause_excludes_swapTwisted`).

Join 3, the declared step-to-clock dictionary.  The extended record
carries one positive real `stepDuration` identifying evolution step `n`
with the proper-time parameter `n * stepDuration` on the record's clock
worldline.  For every inhabitant: the internal clock at step `n` is
`exp (-i · mass · n · stepDuration)` (`stepClock_eq_exp`); the clock
repeats after `N` steps exactly when `mass · N · stepDuration` lies in
`2π ℤ` (`stepClock_periodic_iff`, the forward direction through the
committed exponential criterion `exp_neg_I_mul_real_eq_one_iff` and the
reverse through the committed period classification
`internalClock_periodic_iff`); the worldline events of successive steps are
causally ordered in the record's cone order (`stepEvent_causal_mono`);
distinct steps are distinct events (`stepEvent_injective`); and step zero
is the origin.  Dropping positivity is load-bearing: a zero step collapses
every step to the origin (`zeroStepExtension_collapses`) and a negative
step breaks the causal clause between steps zero and one
(`negativeStepExtension_breaks_causal`).

Join 4, the first-law join at the instrument's preparation.  The extended
record carries a `RepairLawData (Fin 2) (Fin 2)` whose reference law is,
entry by entry, the diagonal of the instrument's preparation, and a
horizon record over it.  The lane-E register `registerOfHorizonRecord` at
that law has the thermodynamic first-law datum as its modular-flow supply
by definition, so the composed first-law clause holds on the simplex
tangent space `MassZero (Fin 2)` at the instrument's preparation
(`register_first_law`).  The explicit inhabitant uses the reference
`(111/179, 68/179)`.

The extended record `InstrumentedCommonWorldArchitecture`, one `Carries*`
predicate per join, the composed receipt `instrumentedCommonWorld_receipt`
from the typed antecedent bundle `InstrumentedWorldAntecedents`, the
explicit inhabitant `instrumentedCommittedWitness` extending
`bridgedCommittedWitness` with `luedersRepeatablePhaseInstrument`, lane
B's `demoScaledBundle` at `h = 1/2`, `stepDuration = 1`, and the run-state
repair law, and nonvacuity are stated at the end, together with
`declared_dictionaries_not_forced`: two inhabitants over the same base
differ in the step duration, two differ in the certified step, two differ
in the horizon split, and the effect table does not determine the
instrument.

DECLARED DATA.  Every base-record field stays a declared hypothesis as in
the committed witness.  The scaled bundle, its step `h`, the instrument,
the decision bundle, the step duration, the repair law, and the horizon
record are declared data constrained by the stated clauses; no committed
theorem selects any of them, and `declared_dictionaries_not_forced` records
this for the step duration, the certified step, the horizon split, and the
instrument.

PREMISE CONSUMPTION PER REGISTER ROW.  PR-15 (clock and energy
calibration) is open: the step-to-clock dictionary is declared, identifies
no physical time, and calibrates no clock.  PR-66 is relocated as in lane
B: the unit-step recurrence is the Euler-Lagrange equation of lane B's
declared discrete action, and the record carries the unit-step law as a
declared field beside its certified replacement.  PR-53 (physical
propagation, frame, and comparison attachment) and PR-54 (source
gauge-field, current, and action attachment) are open.  PR-64 (instrument
channel clauses constructed; readback and selection open), PR-65 (source
preparation and provenance open), and PR-03 (cross-context additivity
open) are carried through lane A unchanged.  PR-04 and PR-02 are consumed
through the decision bundle as in lane A.  PR-07 (repair law) and PR-42
(horizon record) are consumed as declared shapes through join 4.  PR-18
and PR-19 are consumed through the flux packet as in the bridged receipt.
PR-52 (observer-to-spacetime attachment) is open.

FALSIFIER.  The module fails if some Courant constant of a certified bundle
lies below `5`, if some certified bundle carries the unit step, if the
instrument's diagonal channel differs from the diagonal pinching at some
matrix, if the record's state or the run state moves under that channel,
if the swap-twisted instrument satisfies the pinching identity at the run
state, if the step clock fails the exponential or period law, if some pair
of successive step events is causally unordered, if the zero or negative
step fails to break its clause, if the reference law differs from the
run-state diagonal at some entry, if the first-law clause fails on some
mass-zero variation, or if the extended record type is uninhabited over
the committed inhabitant.

WHAT IS NOT PROVED (the joins absent from the record; CW1 stays open).
No common action across carriers exists: the screen action of lane B and
the kinematic flow of the Lorentz module share no variational principle
on one carrier (PR-54 open).  No port, seam, or face of the screen carrier
is mapped to a point, interval, or cone of the Lorentz module (PR-53
open); the step dictionary identifies indices with proper-time parameters
and carries no spatial content, and its link to the Maxwell island is the
shared step index type `ℕ` alone: no clause of the record relates a
Maxwell quantity at step `n` to the clock or the worldline event at step
`n`.  No matter dynamics beyond the committed
structure premises enters any field.  No observer readout exists: the
instrument carries channel clauses and a repeatability proxy, no public
outcome, readback, or provenance (PR-64 readback, PR-65 open).  The
identification of the repair law's reference with the run-state diagonal
is declared, and the register's Einstein clause holds by definition at
lane E's chart as recorded there.  The step duration is a declared
positive real with no unit; no physical time, clock calibration,
spacetime, or readout claim is made.  Nothing here closes issue #740, and
no observation-ledger row is promoted: OL-N1 stays owed.

Axiom audit.  Every proof composes committed receipts with exact finite
mathematics; the module adds no project axiom and uses no native decision
procedure.  The guard lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

namespace OPH.CommonWorldInstrumentJoin

open OPH.C1Lorentz OPH.C2Soldering OPH.CommonWorld OPH.CommonWorldIslandBridge
open OPH.CausalComposition OPH.InverseSquareShellLaw
open OPH.TemporalMaxwellEvolution OPH.ScaledMaxwellStability
open OPH.LocalFaceMaxwellAction OPH.DiscreteCoulombGreen
open OPH.QFT EventAlgebra
open OPH.Dynamics (IsCPTP)
open OPH.Thermodynamics OPH.EinsteinBranch OPH.EinsteinBranchInhabitant

noncomputable section

/-! ## Join 1: the certified scaled law beside the declared unit-step law -/

/-- A Courant constant dominates every eigenvalue of the committed local
normal operator `CᵀC` that carries a nonzero eigenvector: the Courant
hypothesis at the eigenvector reads `lam ‖v‖² ≤ Λ ‖v‖²`. -/
theorem courantBound_ge_eigenvalue {Λ lam : ℝ} (hΛ : CourantBound Λ)
    {v : Fin 30 → ℝ} (hv : v ≠ 0)
    (hev : localMaxwellOperator v = lam • v) : lam ≤ Λ := by
  have hc := hΛ v
  have hE : faceEnergy (faceCurvature v) = lam * realSeamEnergy v := by
    rw [← faceInner_self_eq_energy, faceCurvature_codifferential_adjoint]
    have h1 : faceCodifferential (faceCurvature v) = localMaxwellOperator v :=
      rfl
    rw [h1, hev, seamInner_smul_right, realSeamInner_self_eq_energy]
  rw [hE] at hc
  exact le_of_mul_le_mul_right hc (seamEnergy_pos_of_ne_zero v hv)

/-- Every certified bundle's Courant constant is at least the committed
eigenvalue `5`. -/
theorem scaledBundle_courant_ge_five (S : ScaledMaxwellBundle) : 5 ≤ S.Λ :=
  courantBound_ge_eigenvalue S.courant fiveMode_ne_zero fiveMode_eigen

/-- Every certified bundle's Courant constant is at least the sharp
committed eigenvalue `3 + √5`. -/
theorem scaledBundle_courant_ge_golden (S : ScaledMaxwellBundle) :
    3 + Real.sqrt 5 ≤ S.Λ :=
  courantBound_ge_eigenvalue S.courant goldenMode_ne_zero goldenMode_eigen

/-- Every certified bundle has `h² · 5 < 4`. -/
theorem scaledBundle_step_sq_five_lt (S : ScaledMaxwellBundle) :
    S.h ^ 2 * 5 < 4 := by
  have h5 := scaledBundle_courant_ge_five S
  have hc := S.courant_strict
  have hsq : 0 ≤ S.h ^ 2 := sq_nonneg _
  nlinarith

/-- Every certified bundle has step strictly below one. -/
theorem scaledBundle_step_lt_one (S : ScaledMaxwellBundle) : S.h < 1 := by
  have h := scaledBundle_step_sq_five_lt S
  have hpos := S.h_pos
  by_contra hge
  have hge' : 1 ≤ S.h := le_of_not_gt hge
  nlinarith

/-- No certified bundle carries the unit step. -/
theorem scaledBundle_step_ne_one (S : ScaledMaxwellBundle) : S.h ≠ 1 :=
  (scaledBundle_step_lt_one S).ne

/-- The committed unit-step Ampere law of every Maxwell bundle is the
`h = 1` instance of the scaled law, through lane B's recovery theorem. -/
theorem maxwell_is_unit_scaled (M : MaxwellEvolutionBundle) :
    AmpereEvolutionScaled 1 M.A M.phi M.J :=
  (ampereEvolutionScaled_one_iff M.A M.phi M.J).mpr M.ampere

/-! ## Join 2: the instrument island on the shared record algebra -/

/-- Every one-block projective partition has its single projector equal to
the identity: completeness over `Fin 1` is the single-term sum. -/
theorem projectivePartition_one_block_proj {n : ℕ}
    (P : ProjectivePartition n 1) (i : Fin 1) : P.proj i = 1 := by
  have h := P.complete
  rw [Fin.sum_univ_one] at h
  rw [Subsingleton.elim i 0]
  exact h

/-- The pinching of the record's one-block partition is the identity on
every matrix of the record algebra. -/
theorem record_partition_pinching (W : CommonWorldArchitecture 2 1)
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    partitionPinching W.partition X = X := by
  unfold partitionPinching
  rw [Fin.sum_univ_one, projectivePartition_one_block_proj, one_mul, mul_one]

/-- The committed witness's partition field is the one-block partition,
definitionally. -/
theorem committedWitness_partition :
    committedWitness.partition = oneBlockPartition 2 :=
  rfl

/-- The committed witness's state field is the pure record state,
definitionally. -/
theorem committedWitness_state : committedWitness.state = witnessState :=
  rfl

/-- Obstruction: the committed witness's state is the pure record state,
and it differs from the committed run state at entry `(0,0)`: `1` against
`111/179`. -/
theorem committedWitness_state_ne_runState :
    (committedWitness.state : Matrix (Fin 2) (Fin 2) ℂ) ≠ committedRunState := by
  intro h
  have h00 := congrFun (congrFun h 0) 0
  have hl : (committedWitness.state : Matrix (Fin 2) (Fin 2) ℂ) 0 0 = 1 := by
    show Matrix.diagonal (Pi.single 0 1) 0 0 = 1
    simp
  have hr : committedRunState 0 0 = 111 / 179 := by
    rw [committedRunState_eq_literal]
    simp
  rw [hl, hr] at h00
  norm_num at h00

/-- Obstruction: the single projector of the record's partition is the
identity, and the identity differs from the record effect of the
instrument's diagonal context at entry `(1,1)`. -/
theorem record_partition_proj_ne_record_effect (W : CommonWorldArchitecture 2 1) :
    W.partition.proj 0 ≠
      committedEffectPair (InstrumentContext.web WebContext.diagonal) 0 := by
  rw [projectivePartition_one_block_proj, committedEffectPair_zero,
    committedContextEffect_diagonal_eq]
  intro h
  have h11 := congrFun (congrFun h 1) 1
  simp at h11

/-- The record effect of the instrument's diagonal context is the first
projector of the committed diagonal partition. -/
theorem committedEffectPair_diagonal_zero :
    committedEffectPair (InstrumentContext.web WebContext.diagonal) 0 =
      (diagonalPartition 2).proj 0 := by
  rw [committedEffectPair_zero, committedContextEffect_diagonal_eq]
  show _ = Matrix.diagonal (Pi.single 0 1)
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- The complement effect of the instrument's diagonal context is the
second projector of the committed diagonal partition. -/
theorem committedEffectPair_diagonal_one :
    committedEffectPair (InstrumentContext.web WebContext.diagonal) 1 =
      (diagonalPartition 2).proj 1 := by
  rw [committedEffectPair_one, committedContextEffect_diagonal_eq]
  show _ = Matrix.diagonal (Pi.single 1 1)
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- The diagonal effect pair of the instrument's diagonal context is
exactly the committed diagonal partition. -/
theorem committedEffectPair_diagonal (i : Fin 2) :
    committedEffectPair (InstrumentContext.web WebContext.diagonal) i =
      (diagonalPartition 2).proj i := by
  fin_cases i
  · exact committedEffectPair_diagonal_zero
  · exact committedEffectPair_diagonal_one

/-- The typed transport: the committed diagonal partition refines the
record's one-block partition.  Every partition refines a one-block
partition, whose projector is the identity; the statement records the
shared carrier and carries no further content. -/
theorem diagonal_refines_record (W : CommonWorldArchitecture 2 1) (i : Fin 2) :
    ∃ j : Fin 1,
      W.partition.proj j * (diagonalPartition 2).proj i =
        (diagonalPartition 2).proj i :=
  ⟨0, by rw [projectivePartition_one_block_proj, one_mul]⟩

/-- The record partition's pinching absorbs the diagonal pinching: the
pinching tower of the refinement. -/
theorem record_pinching_absorbs_diagonal (W : CommonWorldArchitecture 2 1)
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    partitionPinching W.partition (partitionPinching (diagonalPartition 2) X) =
      partitionPinching (diagonalPartition 2) X :=
  record_partition_pinching W _

/-- **The instrument-to-pinching identity.**  The Lüders instrument's
diagonal-context summed channel is the committed pinching of the diagonal
partition on every matrix of the shared record algebra. -/
theorem lueders_diagonal_channel_eq_pinching (X : Matrix (Fin 2) (Fin 2) ℂ) :
    luedersPhaseInstrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) X =
      partitionPinching (diagonalPartition 2) X := by
  have hP := committedContextEffect_isEvent (InstrumentContext.web WebContext.diagonal)
  rw [luedersPhaseInstrument_summedChannel, luedersChannel_apply, hP.1.eq,
    hP.compl.1.eq]
  unfold partitionPinching
  rw [Fin.sum_univ_two, ← committedEffectPair_diagonal_zero,
    ← committedEffectPair_diagonal_one, committedEffectPair_zero,
    committedEffectPair_one]

/-- The diagonal pinching fixes every diagonal matrix. -/
theorem diagonal_pinching_fixes_diagonal (d : Fin 2 → ℂ) :
    partitionPinching (diagonalPartition 2) (Matrix.diagonal d) =
      Matrix.diagonal d := by
  apply partitionPinching_fixes
  intro i
  show Matrix.diagonal d * Matrix.diagonal (Pi.single i 1) =
    Matrix.diagonal (Pi.single i 1) * Matrix.diagonal d
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext j
  ring

/-- The committed run state is a fixed point of the diagonal pinching. -/
theorem diagonal_pinching_runState :
    partitionPinching (diagonalPartition 2) committedRunState = committedRunState :=
  diagonal_pinching_fixes_diagonal _

/-- The committed witness state is a fixed point of the diagonal
pinching. -/
theorem diagonal_pinching_witnessState :
    partitionPinching (diagonalPartition 2)
        (witnessState : Matrix (Fin 2) (Fin 2) ℂ) = witnessState :=
  diagonal_pinching_fixes_diagonal _

/-- The record's state is certain of the record effect: the diagonal
outcome-`0` trace of the Lüders instrument on the witness state is one,
against `111/179` on the run state. -/
theorem witnessState_record_certain :
    (luedersPhaseInstrument.outcomeMap (InstrumentContext.web WebContext.diagonal) 0
      (witnessState : Matrix (Fin 2) (Fin 2) ℂ)).trace = 1 := by
  rw [luedersPhaseInstrument.inducedEffect, committedEffectPair_zero,
    committedContextEffect_diagonal_eq]
  show bornWeight (Matrix.diagonal (Pi.single 0 1)) _ = 1
  simp [bornWeight, Matrix.trace_fin_two, Matrix.diagonal_mul]

/-- The swap-twisted summed channel is the Lüders summed channel followed
by swap conjugation. -/
theorem swapTwisted_summedChannel (c : InstrumentContext)
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    swapTwistedPhaseInstrument.summedChannel c X =
      luedersOutcomeMap swapUnitary (luedersPhaseInstrument.summedChannel c X) := by
  unfold PhaseInstrument.summedChannel
  rw [swapTwistedPhaseInstrument_outcomeMap_comp,
    swapTwistedPhaseInstrument_outcomeMap_comp]
  simp only [LinearMap.add_apply, LinearMap.comp_apply, map_add]

/-- The swap-twisted instrument's diagonal channel moves the run state: it
exchanges the two diagonal weights. -/
theorem swapTwisted_diagonal_channel_runState :
    swapTwistedPhaseInstrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) committedRunState =
      !![(68 / 179 : ℂ), 0; 0, 111 / 179] := by
  rw [swapTwisted_summedChannel, lueders_diagonal_channel_eq_pinching,
    diagonal_pinching_runState, luedersOutcomeMap_apply,
    committedRunState_eq_literal, swapUnitary]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Matrix.mul_apply, Fin.sum_univ_two, Matrix.conjTranspose_apply]

/-! ## Join 3: the step-to-clock dictionary, general lemmas -/

/-- The declared proper-time parameter of step `n` at step duration `δ`. -/
def stepTime (δ : ℝ) (n : ℕ) : ℝ := n * δ

theorem stepTime_zero (δ : ℝ) : stepTime δ 0 = 0 := by
  simp [stepTime]

theorem stepTime_one (δ : ℝ) : stepTime δ 1 = δ := by
  simp [stepTime]

theorem stepTime_add (δ : ℝ) (n m : ℕ) :
    stepTime δ (n + m) = stepTime δ n + stepTime δ m := by
  simp [stepTime]
  ring

/-- The frame worldline at parameter zero is the origin. -/
theorem frameWorldline_zero (f : FrameHyperboloid) : frameWorldline f 0 = 0 := by
  show (0 : ℝ) • (f : Herm2) = 0
  exact zero_smul ℝ _

/-- The frame worldline is injective in its parameter: the frame has
positive time component. -/
theorem frameWorldline_injective (f : FrameHyperboloid) :
    Function.Injective (frameWorldline f) := by
  intro τ σ h
  have h1 : (frameWorldline f τ).1 = (frameWorldline f σ).1 := by rw [h]
  have hτ : (frameWorldline f τ).1 = τ * (f : Herm2).1 := rfl
  have hσ : (frameWorldline f σ).1 = σ * (f : Herm2).1 := rfl
  rw [hτ, hσ] at h1
  exact mul_right_cancel₀ f.2.2.ne' h1

/-- A negative step duration breaks the causal clause between steps zero
and one: the displacement has negative time component. -/
theorem negative_step_not_causal (f : FrameHyperboloid) {δ : ℝ} (hδ : δ < 0) :
    ¬ causalLE (frameWorldline f (stepTime δ 0)) (frameWorldline f (stepTime δ 1)) := by
  intro h
  have h' : IsFutureCausal
      (frameWorldline f (stepTime δ 1) - frameWorldline f (stepTime δ 0)) := h
  obtain ⟨_, ht⟩ := h'
  have hdiff : (frameWorldline f (stepTime δ 1) -
      frameWorldline f (stepTime δ 0)).1 = δ * (f : Herm2).1 := by
    show stepTime δ 1 * (f : Herm2).1 - stepTime δ 0 * (f : Herm2).1 =
      δ * (f : Herm2).1
    rw [stepTime_one, stepTime_zero]
    ring
  rw [hdiff] at ht
  have hneg := mul_neg_of_neg_of_pos hδ f.2.2
  linarith

/-! ## Join 4: the run-state repair law and horizon record -/

/-- The reference law read from the diagonal of the committed run state. -/
def runRef (x : Fin 2) : ℝ := if x = 0 then 111 / 179 else 68 / 179

theorem runRef_zero : runRef 0 = 111 / 179 := by simp [runRef]
theorem runRef_one : runRef 1 = 68 / 179 := by simp [runRef]

/-- The repair law (PR-07 shape) whose reference is the run-state diagonal,
with the identity as the repaired visible datum. -/
def runRepairLaw : RepairLawData (Fin 2) (Fin 2) where
  ref := runRef
  ref_pos := fun x => by
    unfold runRef
    split_ifs <;> norm_num
  ref_law := by
    rw [Fin.sum_univ_two, runRef_zero, runRef_one]
    norm_num
  visible := id

/-- The reference law equals the diagonal of the committed run state,
entry by entry. -/
theorem runRepairLaw_ref_eq_runState (x : Fin 2) :
    (runRepairLaw.ref x : ℂ) = committedRunState x x := by
  fin_cases x
  · show (runRef 0 : ℂ) = committedRunState 0 0
    rw [runRef_zero, committedRunState_eq_literal]
    push_cast
    simp
  · show (runRef 1 : ℂ) = committedRunState 1 1
    rw [runRef_one, committedRunState_eq_literal]
    push_cast
    simp

/-- A nonconstant boost charge. -/
def runBoost (x : Fin 2) : ℝ := if x = 0 then 1 else 0

/-- The central charge forced by the split. -/
def runCentral (x : Fin 2) : ℝ := -Real.log (runRef x) - 2 * Real.pi * runBoost x

/-- The horizon record (PR-42 shape) over the run-state repair law. -/
def runHorizonRecord : HorizonRecordData runRepairLaw where
  boostCharge := runBoost
  centralCharge := runCentral
  split := fun x => by
    show -Real.log (runRef x) = 2 * Real.pi * runBoost x + runCentral x
    unfold runCentral
    ring
  central_measurable := fun x y h => by rw [show x = y from h]

/-- A second boost charge over the same law, used to show the split is
declared. -/
def runBoostAlt (x : Fin 2) : ℝ := if x = 0 then 0 else 1

def runCentralAlt (x : Fin 2) : ℝ :=
  -Real.log (runRef x) - 2 * Real.pi * runBoostAlt x

/-- A second horizon record over the same law. -/
def runHorizonRecordAlt : HorizonRecordData runRepairLaw where
  boostCharge := runBoostAlt
  centralCharge := runCentralAlt
  split := fun x => by
    show -Real.log (runRef x) = 2 * Real.pi * runBoostAlt x + runCentralAlt x
    unfold runCentralAlt
    ring
  central_measurable := fun x y h => by rw [show x = y from h]

/-! ## The extended record -/

/-- The instrumented common-world architecture: the committed bridged
record extended by four typed joins.  Join 1 carries a certified scaled
bundle sharing the carrier, the initial seam data, and the seam current of
the record's unit-step bundle.  Join 2 carries a repeatable phase
instrument and a decision bundle with the pinching identity on the shared
record algebra and the fixed-point clause for the record's state.  Join 3
carries a declared positive step duration.  Join 4 carries a repair law
whose reference is the diagonal of the instrument's preparation and a
horizon record over it.  Every base field stays a declared hypothesis
exactly as in the committed witness. -/
structure InstrumentedCommonWorldArchitecture extends
    BridgedCommonWorldArchitecture 2 1 where
  /-- The declared certified scaled bundle on the shared carrier. -/
  scaled : ScaledMaxwellBundle
  /-- Shared initial seam configuration with the unit-step bundle. -/
  scaled_initial_seam : scaled.A 0 = maxwell.A 0 ∧ scaled.A 1 = maxwell.A 1
  /-- Shared seam current with the unit-step bundle. -/
  scaled_current : scaled.J = maxwell.J
  /-- The declared repeatable phase instrument. -/
  instrument : RepeatablePhaseInstrument
  /-- The declared decision bundle (PR-04 and PR-02 as in lane A). -/
  decision : PhaseInstrumentDecisionBundle
  /-- The pinching identity: the instrument's diagonal-context summed
  channel is the committed diagonal pinching on every matrix. -/
  diagonal_channel_pinching : ∀ X : Matrix (Fin 2) (Fin 2) ℂ,
    instrument.toPhaseInstrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) X =
      partitionPinching (diagonalPartition 2) X
  /-- The record's state is a fixed point of the diagonal pinching. -/
  state_diagonal_fixed :
    partitionPinching (diagonalPartition 2) state.1 = state.1
  /-- The declared step duration. -/
  stepDuration : ℝ
  /-- Declared positivity of the step duration. -/
  stepDuration_pos : 0 < stepDuration
  /-- The declared repair law (PR-07 shape). -/
  repairLaw : RepairLawData (Fin 2) (Fin 2)
  /-- The reference law is the diagonal of the instrument's preparation. -/
  repairLaw_ref : ∀ x : Fin 2, (repairLaw.ref x : ℂ) = decision.prep x x
  /-- The declared horizon record (PR-42 shape). -/
  horizon : HorizonRecordData repairLaw

/-! ## Join 1 on the record -/

/-- The certified-evolution clause of the composed receipt.  The final
conjunct, the unit-step instability, is a theorem about the unit-step law
on the carrier type and takes no field of the record; the record's own
unit-step bundle is one instance of that law. -/
def CarriesCertifiedEvolution (W : InstrumentedCommonWorldArchitecture) : Prop :=
  AmpereEvolutionScaled 1 W.maxwell.A W.maxwell.phi W.maxwell.J ∧
  (5 ≤ W.scaled.Λ ∧ 3 + Real.sqrt 5 ≤ W.scaled.Λ) ∧
  W.scaled.h < 1 ∧
  (W.scaled.A 0 = W.maxwell.A 0 ∧ W.scaled.A 1 = W.maxwell.A 1 ∧
    W.scaled.J = W.maxwell.J) ∧
  (∀ n, 0 ≤ fieldEnergyScaled W.scaled.h W.scaled.A W.scaled.phi n) ∧
  (W.maxwell.J = (fun _ ↦ 0) → ∀ n,
    realSeamEnergy (electricFieldScaled W.scaled.h W.scaled.A W.scaled.phi n) ≤
        8 * fieldEnergyScaled W.scaled.h W.scaled.A W.scaled.phi 0 /
          (4 - W.scaled.h ^ 2 * W.scaled.Λ) ∧
      faceEnergy (magneticField W.scaled.A n) ≤
        16 * fieldEnergyScaled W.scaled.h W.scaled.A W.scaled.phi 0 /
          (4 - W.scaled.h ^ 2 * W.scaled.Λ)) ∧
  (∃ A : ℕ → Fin 30 → ℝ, AmpereEvolution A (fun _ ↦ 0) (fun _ ↦ 0) ∧
    ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricField A (fun _ ↦ 0) n))

theorem carriesCertifiedEvolution (W : InstrumentedCommonWorldArchitecture) :
    CarriesCertifiedEvolution W := by
  refine ⟨maxwell_is_unit_scaled W.maxwell,
    ⟨scaledBundle_courant_ge_five W.scaled, scaledBundle_courant_ge_golden W.scaled⟩,
    scaledBundle_step_lt_one W.scaled,
    ⟨W.scaled_initial_seam.1, W.scaled_initial_seam.2, W.scaled_current⟩,
    fun n => fieldEnergyScaled_nonneg W.scaled.h W.scaled.h_pos.ne' W.scaled.Λ
      W.scaled.courant W.scaled.courant_strict.le W.scaled.A W.scaled.phi n,
    fun hJ n => ?_,
    unit_step_instability⟩
  have hAmp0 : AmpereEvolutionScaled W.scaled.h W.scaled.A W.scaled.phi (fun _ ↦ 0) := by
    have h := W.scaled.ampere
    rw [W.scaled_current, hJ] at h
    exact h
  exact stability_certificate W.scaled.h W.scaled.h_pos.ne' W.scaled.Λ
    W.scaled.Λ_nonneg W.scaled.courant W.scaled.courant_strict W.scaled.A
    W.scaled.phi hAmp0 n

/-! ## Join 2 on the record -/

/-- The instrument-join clause of the composed receipt: the pinching
identity, its absorption by the record partition, the refinement, the two
fixed points, the Born readout of the record's state through the
instrument, repeatability, and the CPTP summed channels. -/
def CarriesInstrumentJoin (W : InstrumentedCommonWorldArchitecture) : Prop :=
  (∀ X : Matrix (Fin 2) (Fin 2) ℂ,
    W.instrument.toPhaseInstrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) X =
      partitionPinching (diagonalPartition 2) X) ∧
  (∀ X : Matrix (Fin 2) (Fin 2) ℂ,
    partitionPinching W.partition
        (W.instrument.toPhaseInstrument.summedChannel
          (InstrumentContext.web WebContext.diagonal) X) =
      W.instrument.toPhaseInstrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) X) ∧
  (∀ i : Fin 2, ∃ j : Fin 1,
    W.partition.proj j * (diagonalPartition 2).proj i =
      (diagonalPartition 2).proj i) ∧
  W.instrument.toPhaseInstrument.summedChannel
      (InstrumentContext.web WebContext.diagonal) W.state.1 = W.state.1 ∧
  W.instrument.toPhaseInstrument.summedChannel
      (InstrumentContext.web WebContext.diagonal) W.decision.prep =
    W.decision.prep ∧
  (∀ i : Fin 2,
    (W.instrument.toPhaseInstrument.outcomeMap
        (InstrumentContext.web WebContext.diagonal) i W.state.1).trace =
      bornWeight W.state.1 ((diagonalPartition 2).proj i)) ∧
  W.instrument.toPhaseInstrument.Repeatable ∧
  (∀ c : InstrumentContext, IsCPTP (W.instrument.toPhaseInstrument.summedChannel c))

theorem carriesInstrumentJoin (W : InstrumentedCommonWorldArchitecture) :
    CarriesInstrumentJoin W := by
  refine ⟨W.diagonal_channel_pinching,
    fun X => by rw [W.diagonal_channel_pinching]; exact record_partition_pinching _ _,
    diagonal_refines_record W.toCommonWorldArchitecture,
    by rw [W.diagonal_channel_pinching, W.state_diagonal_fixed],
    by rw [W.diagonal_channel_pinching, W.decision.prep_eq_run_state,
      diagonal_pinching_runState],
    fun i => by
      rw [W.instrument.toPhaseInstrument.inducedEffect, committedEffectPair_diagonal],
    W.instrument.repeatable,
    fun c => W.instrument.toPhaseInstrument.summedChannel_isCPTP c⟩

/-! ## Join 3 on the record -/

/-- The worldline event of step `n`. -/
def stepEvent (W : InstrumentedCommonWorldArchitecture) (n : ℕ) : Herm2 :=
  frameWorldline W.frame (stepTime W.stepDuration n)

/-- The internal clock read at step `n`. -/
def stepClock (W : InstrumentedCommonWorldArchitecture) (n : ℕ) : ℂ :=
  internalClock W.mass W.frame (stepTime W.stepDuration n)

/-- The clock at step `n` is the exact rotation
`exp (-i · mass · n · stepDuration)`. -/
theorem stepClock_eq_exp (W : InstrumentedCommonWorldArchitecture) (n : ℕ) :
    stepClock W n =
      Complex.exp (-Complex.I * ((W.mass * stepTime W.stepDuration n : ℝ) : ℂ)) :=
  internalClock_eq_exp W.mass W.frame _

/-- The clock repeats after `N` steps exactly when
`mass · N · stepDuration` is an integer multiple of `2π`.  The forward
direction passes through the committed exponential criterion
`exp_neg_I_mul_real_eq_one_iff` at step zero; the reverse direction passes
through the committed period classification `internalClock_periodic_iff`. -/
theorem stepClock_periodic_iff (W : InstrumentedCommonWorldArchitecture) (N : ℕ) :
    (∀ n : ℕ, stepClock W (n + N) = stepClock W n) ↔
      ∃ m : ℤ, W.mass * stepTime W.stepDuration N = m * (2 * Real.pi) := by
  constructor
  · intro hper
    have h0 := hper 0
    unfold stepClock at h0
    rw [internalClock_eq_exp, internalClock_eq_exp, Nat.zero_add, stepTime_zero] at h0
    simp only [mul_zero, Complex.ofReal_zero, Complex.exp_zero] at h0
    exact (exp_neg_I_mul_real_eq_one_iff _).mp h0
  · rintro ⟨m, hm⟩ n
    have hper := (internalClock_periodic_iff W.mass W.frame
      (stepTime W.stepDuration N)).mpr ⟨m, hm⟩
    unfold stepClock
    rw [stepTime_add]
    exact hper _

/-- Successive step events are causally ordered in the record's cone
order. -/
theorem stepEvent_causal_mono (W : InstrumentedCommonWorldArchitecture)
    {n m : ℕ} (h : n ≤ m) : causalLE (stepEvent W n) (stepEvent W m) :=
  frameWorldline_causal_mono W.toCommonWorldArchitecture _ _
    (mul_le_mul_of_nonneg_right (Nat.cast_le.mpr h) W.stepDuration_pos.le)

/-- Distinct steps are distinct events. -/
theorem stepEvent_injective (W : InstrumentedCommonWorldArchitecture) :
    Function.Injective (stepEvent W) := by
  intro n m h
  have h1 := frameWorldline_injective W.frame h
  have h2 : (n : ℝ) * W.stepDuration = (m : ℝ) * W.stepDuration := h1
  have h3 := mul_right_cancel₀ W.stepDuration_pos.ne' h2
  exact_mod_cast h3

/-- Step zero is the origin. -/
theorem stepEvent_zero (W : InstrumentedCommonWorldArchitecture) :
    stepEvent W 0 = 0 := by
  unfold stepEvent
  rw [stepTime_zero]
  exact frameWorldline_zero _

/-- The step-dictionary clause of the composed receipt. -/
def CarriesStepDictionary (W : InstrumentedCommonWorldArchitecture) : Prop :=
  (∀ n : ℕ, stepClock W n =
    Complex.exp (-Complex.I * ((W.mass * stepTime W.stepDuration n : ℝ) : ℂ))) ∧
  (∀ N : ℕ, (∀ n : ℕ, stepClock W (n + N) = stepClock W n) ↔
    ∃ m : ℤ, W.mass * stepTime W.stepDuration N = m * (2 * Real.pi)) ∧
  (∀ n m : ℕ, n ≤ m → causalLE (stepEvent W n) (stepEvent W m)) ∧
  Function.Injective (stepEvent W) ∧
  stepEvent W 0 = 0

theorem carriesStepDictionary (W : InstrumentedCommonWorldArchitecture) :
    CarriesStepDictionary W :=
  ⟨stepClock_eq_exp W, stepClock_periodic_iff W,
    fun _ _ h => stepEvent_causal_mono W h, stepEvent_injective W,
    stepEvent_zero W⟩

/-! ## Join 4 on the record -/

/-- The lane-E register at the record's repair law and horizon record. -/
def register (W : InstrumentedCommonWorldArchitecture) :
    EinsteinRegisterData ℕ (MassZero (Fin 2)) TorusPoint :=
  registerOfHorizonRecord W.repairLaw W.horizon

/-- The register's modular-flow supply is the first-law datum at the
record's reference: definitional. -/
theorem register_modularFlowData (W : InstrumentedCommonWorldArchitecture) :
    (register W).modularFlowData =
      thermoFirstLawData (Fin 2) W.repairLaw.ref W.horizon.boostCharge
        W.horizon.centralCharge :=
  rfl

/-- The composed first-law clause at the record's reference, on the simplex
tangent space `MassZero (Fin 2)`. -/
theorem register_first_law (W : InstrumentedCommonWorldArchitecture)
    (δρ : MassZero (Fin 2)) :
    (register W).modularFlowData.entropyDifferential δρ =
      2 * Real.pi * (register W).modularFlowData.bulkPairing δρ +
        (register W).modularFlowData.edgeDifferential δρ :=
  composed_first_law_is_horizon_first_law W.repairLaw W.horizon δρ

/-- The first-law-join clause of the composed receipt. -/
def CarriesFirstLawJoin (W : InstrumentedCommonWorldArchitecture) : Prop :=
  (∀ x : Fin 2, (W.repairLaw.ref x : ℂ) = W.decision.prep x x) ∧
  (∀ x : Fin 2, (W.repairLaw.ref x : ℂ) = committedRunState x x) ∧
  (∑ x : Fin 2, W.repairLaw.ref x) = 1 ∧
  (∀ x : Fin 2, -Real.log (W.repairLaw.ref x) =
    2 * Real.pi * W.horizon.boostCharge x + W.horizon.centralCharge x) ∧
  (register W).modularFlowData =
    thermoFirstLawData (Fin 2) W.repairLaw.ref W.horizon.boostCharge
      W.horizon.centralCharge ∧
  (∀ δρ : MassZero (Fin 2),
    (register W).modularFlowData.entropyDifferential δρ =
      2 * Real.pi * (register W).modularFlowData.bulkPairing δρ +
        (register W).modularFlowData.edgeDifferential δρ)

theorem carriesFirstLawJoin (W : InstrumentedCommonWorldArchitecture) :
    CarriesFirstLawJoin W :=
  ⟨W.repairLaw_ref,
    fun x => by rw [W.repairLaw_ref x, W.decision.prep_eq_run_state],
    W.repairLaw.ref_law, W.horizon.split, register_modularFlowData W,
    register_first_law W⟩

/-! ## The composed receipt on one antecedent bundle -/

/-- The typed antecedent bundle of the instrumented witness: one extended
record and one declared radial flux packet (PR-18 and PR-19, consumed
exactly as in the bridged receipt). -/
structure InstrumentedWorldAntecedents where
  /-- The extended common-world record. -/
  world : InstrumentedCommonWorldArchitecture
  /-- PR-18 and PR-19 as the committed `RadialFluxData` premise packet. -/
  flux : RadialFluxData

/-- **The instrumented common-world receipt (issue #740, four added joins;
CW1 stays open).**  Every antecedent bundle simultaneously carries the ten
committed clauses of the bridged receipt on its one record and the four
join clauses: certified evolution beside the unit-step law on the shared
carrier, the instrument island through the pinching identity on the shared
record algebra, the declared step-to-clock dictionary on the record's own
clock worldline, and the first-law join at the instrument's preparation.
The joins absent from the record are listed in the file header and stay on
the open register rows PR-15, PR-52, PR-53, PR-54, PR-64, and PR-65. -/
theorem instrumentedCommonWorld_receipt (A : InstrumentedWorldAntecedents) :
    CarriesMassShell A.world.toCommonWorldArchitecture ∧
      CarriesInternalClock A.world.toCommonWorldArchitecture ∧
      CarriesFreeEvolution A.world.toCommonWorldArchitecture ∧
      CarriesCausalOrder A.world.toCommonWorldArchitecture ∧
      CarriesKinematicJoins A.world.toCommonWorldArchitecture ∧
      CarriesTemporalMaxwell A.world.toCommonWorldArchitecture ∧
      CarriesElectroweakBreaking A.world.toCommonWorldArchitecture ∧
      CarriesScreenPortMeasure A.world.toCommonWorldArchitecture ∧
      CarriesIslandBridge A.world.toBridgedCommonWorldArchitecture ∧
      CarriesBridgedShellReadout A.world.toBridgedCommonWorldArchitecture A.flux ∧
      CarriesCertifiedEvolution A.world ∧
      CarriesInstrumentJoin A.world ∧
      CarriesStepDictionary A.world ∧
      CarriesFirstLawJoin A.world := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10⟩ :=
    bridgedCommonWorld_receipt ⟨A.world.toBridgedCommonWorldArchitecture, A.flux⟩
  exact ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10,
    carriesCertifiedEvolution A.world, carriesInstrumentJoin A.world,
    carriesStepDictionary A.world, carriesFirstLawJoin A.world⟩

/-! ## The explicit inhabitant extending the bridged inhabitant -/

/-- The explicit inhabitant: the committed bridged witness extended by lane
B's demonstration bundle at `h = 1/2`, lane A's Lüders repeatable
instrument with the recorded decision bundle, unit step duration, and the
run-state repair law with its horizon record.  The parameter values select
nothing. -/
def instrumentedCommittedWitness : InstrumentedCommonWorldArchitecture where
  toBridgedCommonWorldArchitecture := bridgedCommittedWitness
  scaled := demoScaledBundle
  scaled_initial_seam := by
    constructor
    · show demoScaledA (1 / 2) 0 = demoA 0
      rw [demoScaledA_zero]
      exact (show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]).symm
    · show demoScaledA (1 / 2) 1 = demoA 1
      rw [demoScaledA_one]
      exact (show demoA 1 = demoInitial from by simp only [demoA]).symm
  scaled_current := rfl
  instrument := luedersRepeatablePhaseInstrument
  decision := recordedInstrumentBundle
  diagonal_channel_pinching := lueders_diagonal_channel_eq_pinching
  state_diagonal_fixed := diagonal_pinching_witnessState
  stepDuration := 1
  stepDuration_pos := one_pos
  repairLaw := runRepairLaw
  repairLaw_ref := runRepairLaw_ref_eq_runState
  horizon := runHorizonRecord

/-- The extended inhabitant restricts to the bridged inhabitant
definitionally. -/
theorem instrumentedCommittedWitness_extends :
    instrumentedCommittedWitness.toBridgedCommonWorldArchitecture =
      bridgedCommittedWitness :=
  rfl

/-- The extended record is inhabited. -/
theorem instrumentedCommonWorld_inhabited :
    Nonempty InstrumentedCommonWorldArchitecture :=
  ⟨instrumentedCommittedWitness⟩

/-- The explicit antecedent bundle of the composed receipt. -/
def committedInstrumentedAntecedents : InstrumentedWorldAntecedents :=
  ⟨instrumentedCommittedWitness, unitFlux⟩

/-- The composed receipt holds of the one explicit inhabitant and the one
demonstration flux packet. -/
theorem instrumentedCommittedWitness_receipt :
    CarriesMassShell instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesInternalClock instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesFreeEvolution instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesCausalOrder instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesKinematicJoins instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesTemporalMaxwell instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesElectroweakBreaking
        instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesScreenPortMeasure
        instrumentedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesIslandBridge instrumentedCommittedWitness.toBridgedCommonWorldArchitecture ∧
      CarriesBridgedShellReadout
        instrumentedCommittedWitness.toBridgedCommonWorldArchitecture unitFlux ∧
      CarriesCertifiedEvolution instrumentedCommittedWitness ∧
      CarriesInstrumentJoin instrumentedCommittedWitness ∧
      CarriesStepDictionary instrumentedCommittedWitness ∧
      CarriesFirstLawJoin instrumentedCommittedWitness :=
  instrumentedCommonWorld_receipt committedInstrumentedAntecedents

/-- The inhabitant's two screen laws: the certified bundle is bounded at
every step while the unit-step law on the same carrier admits an unbounded
zero-current mode. -/
theorem instrumentedCommittedWitness_contrast :
    (∀ n, realSeamEnergy (electricFieldScaled instrumentedCommittedWitness.scaled.h
        instrumentedCommittedWitness.scaled.A instrumentedCommittedWitness.scaled.phi n) ≤
      8 * fieldEnergyScaled instrumentedCommittedWitness.scaled.h
        instrumentedCommittedWitness.scaled.A instrumentedCommittedWitness.scaled.phi 0 /
        (4 - instrumentedCommittedWitness.scaled.h ^ 2 *
          instrumentedCommittedWitness.scaled.Λ)) ∧
    (∃ A : ℕ → Fin 30 → ℝ, AmpereEvolution A (fun _ ↦ 0) (fun _ ↦ 0) ∧
      ∀ M : ℝ, ∃ n : ℕ, M < realSeamEnergy (electricField A (fun _ ↦ 0) n)) :=
  ⟨fun n => ((carriesCertifiedEvolution instrumentedCommittedWitness).2.2.2.2.2.1
    rfl n).1, unit_step_instability⟩

/-- The inhabitant's record state is certain of the record effect while the
instrument's preparation is not: `1` against `111/179`. -/
theorem instrumentedCommittedWitness_state_contrast :
    (instrumentedCommittedWitness.instrument.toPhaseInstrument.outcomeMap
        (InstrumentContext.web WebContext.diagonal) 0
        instrumentedCommittedWitness.state.1).trace = 1 ∧
    (instrumentedCommittedWitness.instrument.toPhaseInstrument.outcomeMap
        (InstrumentContext.web WebContext.diagonal) 0
        instrumentedCommittedWitness.decision.prep).trace = 111 / 179 :=
  ⟨witnessState_record_certain, luedersPhaseInstrument_run_diagonal⟩

/-! ## Load-bearing receipts: dropping a clause readmits an excluded pair -/

/-- The instrument extension with the pinching clause and repeatability
dropped. -/
structure BareInstrumentedArchitecture extends BridgedCommonWorldArchitecture 2 1 where
  /-- A phase instrument with no pinching clause. -/
  instrument : PhaseInstrument
  /-- A decision bundle. -/
  decision : PhaseInstrumentDecisionBundle

/-- The readmitted pair: the bare extension of the bridged witness by the
swap-twisted instrument. -/
def twistedExtension : BareInstrumentedArchitecture where
  toBridgedCommonWorldArchitecture := bridgedCommittedWitness
  instrument := swapTwistedPhaseInstrument
  decision := recordedInstrumentBundle

/-- The swap-twisted extension violates the dropped pinching clause at the
run state. -/
theorem twistedExtension_breaks_pinching :
    twistedExtension.instrument.summedChannel
        (InstrumentContext.web WebContext.diagonal) twistedExtension.decision.prep ≠
      partitionPinching (diagonalPartition 2) twistedExtension.decision.prep := by
  show swapTwistedPhaseInstrument.summedChannel _ committedRunState ≠
    partitionPinching _ committedRunState
  rw [swapTwisted_diagonal_channel_runState, diagonal_pinching_runState,
    committedRunState_eq_literal]
  intro h
  have h00 := congrFun (congrFun h 0) 0
  norm_num at h00

/-- The retained pinching clause is load-bearing: no extended record
carries the swap-twisted instrument. -/
theorem pinching_clause_excludes_swapTwisted
    (W : InstrumentedCommonWorldArchitecture)
    (h : W.instrument.toPhaseInstrument = swapTwistedPhaseInstrument) : False := by
  have hc := W.diagonal_channel_pinching committedRunState
  rw [h, swapTwisted_diagonal_channel_runState, diagonal_pinching_runState,
    committedRunState_eq_literal] at hc
  have h00 := congrFun (congrFun hc 0) 0
  norm_num at h00

/-- The step dictionary with positivity dropped. -/
structure BareStepDictionary extends BridgedCommonWorldArchitecture 2 1 where
  /-- A declared step duration with no sign clause. -/
  stepDuration : ℝ

/-- The zero-step extension of the bridged witness. -/
def zeroStepExtension : BareStepDictionary where
  toBridgedCommonWorldArchitecture := bridgedCommittedWitness
  stepDuration := 0

/-- A zero step collapses every step to the origin. -/
theorem zeroStepExtension_collapses (n : ℕ) :
    frameWorldline zeroStepExtension.frame
      (stepTime zeroStepExtension.stepDuration n) = 0 := by
  show frameWorldline standardFrame (stepTime 0 n) = 0
  have h : stepTime 0 n = 0 := by simp [stepTime]
  rw [h]
  exact frameWorldline_zero _

/-- The negative-step extension of the bridged witness. -/
def negativeStepExtension : BareStepDictionary where
  toBridgedCommonWorldArchitecture := bridgedCommittedWitness
  stepDuration := -1

/-- A negative step breaks the causal clause between steps zero and one. -/
theorem negativeStepExtension_breaks_causal :
    ¬ causalLE
      (frameWorldline negativeStepExtension.frame
        (stepTime negativeStepExtension.stepDuration 0))
      (frameWorldline negativeStepExtension.frame
        (stepTime negativeStepExtension.stepDuration 1)) :=
  negative_step_not_causal _ (by norm_num : (-1 : ℝ) < 0)

/-! ## The declared dictionaries select nothing -/

/-- A second certified bundle at `h = 1/3`, `Λ = 6`, with the same initial
seam data and sources as lane B's demonstration bundle. -/
def demoScaledBundleThird : ScaledMaxwellBundle where
  h := 1 / 3
  Λ := 6
  h_pos := by norm_num
  Λ_nonneg := by norm_num
  courant := committed_courant
  courant_strict := by norm_num
  A := demoScaledA (1 / 3)
  phi := fun _ ↦ 0
  rho := fun _ ↦ 0
  J := fun _ ↦ 0
  ampere := demoScaled_ampere (1 / 3) (by norm_num)
  gauss_init := by
    rw [electricFieldScaled_temporal_gauge, demoScaledA_one, demoScaledA_zero,
      sub_zero, map_neg, map_smul, demoInitial_boundary, smul_zero, neg_zero]
  continuity := by
    intro n
    simp

/-- The four joins absent from the extended record, as named labels: no
common action across carriers (PR-54), no port, seam, or face to module
carrier map (PR-53), no matter dynamics beyond the committed structure
premises, and no observer readout (PR-64 readback, PR-65).  The record has
no field of any of these kinds. -/
inductive OpenJoin : Type
  /-- One variational principle for the screen action and the kinematic
  flow on one carrier. -/
  | commonAction
  /-- Ports, seams, and faces sent to points, intervals, or cones of the
  declared module. -/
  | carrierMap
  /-- Matter dynamics beyond the committed structure premises. -/
  | matterDynamics
  /-- Public outcome, readback, and provenance of the instrument. -/
  | observerReadout

/-- **The declared dictionaries are not forced.**  Two inhabitants over the
same bridged base differ in the step duration (PR-15 open), two differ in
the certified step (PR-53 open), two differ in the horizon split over the
same law (PR-42 declared), and the effect table does not determine the
instrument (PR-64 selection open, through lane A).  No committed theorem
selects any of the four. -/
theorem declared_dictionaries_not_forced :
    (∃ W₁ W₂ : InstrumentedCommonWorldArchitecture,
      W₁.toBridgedCommonWorldArchitecture = W₂.toBridgedCommonWorldArchitecture ∧
        W₁.stepDuration ≠ W₂.stepDuration) ∧
    (∃ W₁ W₂ : InstrumentedCommonWorldArchitecture,
      W₁.toBridgedCommonWorldArchitecture = W₂.toBridgedCommonWorldArchitecture ∧
        W₁.scaled.h ≠ W₂.scaled.h) ∧
    (∃ W₁ W₂ : InstrumentedCommonWorldArchitecture,
      W₁.toBridgedCommonWorldArchitecture = W₂.toBridgedCommonWorldArchitecture ∧
        W₁.repairLaw = W₂.repairLaw ∧
        W₁.horizon.boostCharge ≠ W₂.horizon.boostCharge) ∧
    (∃ Φ Ψ : PhaseInstrument,
      (∀ (c : InstrumentContext) (i : Fin 2) (X : Matrix (Fin 2) (Fin 2) ℂ),
        (Φ.outcomeMap c i X).trace = (Ψ.outcomeMap c i X).trace) ∧
      Φ.Repeatable ∧ ¬ Ψ.Repeatable ∧
      Φ.outcomeMap (InstrumentContext.web WebContext.diagonal) 0 committedRunState ≠
        Ψ.outcomeMap (InstrumentContext.web WebContext.diagonal) 0 committedRunState) := by
  refine ⟨⟨instrumentedCommittedWitness,
      { instrumentedCommittedWitness with stepDuration := 2, stepDuration_pos := two_pos },
      rfl, by show (1 : ℝ) ≠ 2; norm_num⟩,
    ⟨instrumentedCommittedWitness,
      { instrumentedCommittedWitness with
        scaled := demoScaledBundleThird
        scaled_initial_seam := by
          constructor
          · show demoScaledA (1 / 3) 0 = demoA 0
            rw [demoScaledA_zero]
            exact (show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]).symm
          · show demoScaledA (1 / 3) 1 = demoA 1
            rw [demoScaledA_one]
            exact (show demoA 1 = demoInitial from by simp only [demoA]).symm
        scaled_current := rfl },
      rfl, by show (1 / 2 : ℝ) ≠ 1 / 3; norm_num⟩,
    ⟨instrumentedCommittedWitness,
      { instrumentedCommittedWitness with horizon := runHorizonRecordAlt },
      rfl, rfl, ?_⟩,
    effect_table_does_not_determine_instrument⟩
  intro h
  have h0 : runBoost 0 = runBoostAlt 0 := congrFun h 0
  simp [runBoost, runBoostAlt] at h0

end

end OPH.CommonWorldInstrumentJoin

/- Axiom audit: committed receipts and exact finite mathematics only.
Expected axioms per line: at most `propext`, `Classical.choice`,
`Quot.sound`.  No native decision procedure is used. -/

#print axioms OPH.CommonWorldInstrumentJoin.courantBound_ge_eigenvalue
#print axioms OPH.CommonWorldInstrumentJoin.scaledBundle_courant_ge_five
#print axioms OPH.CommonWorldInstrumentJoin.scaledBundle_courant_ge_golden
#print axioms OPH.CommonWorldInstrumentJoin.scaledBundle_step_lt_one
#print axioms OPH.CommonWorldInstrumentJoin.scaledBundle_step_ne_one
#print axioms OPH.CommonWorldInstrumentJoin.maxwell_is_unit_scaled
#print axioms OPH.CommonWorldInstrumentJoin.projectivePartition_one_block_proj
#print axioms OPH.CommonWorldInstrumentJoin.record_partition_pinching
#print axioms OPH.CommonWorldInstrumentJoin.committedWitness_partition
#print axioms OPH.CommonWorldInstrumentJoin.committedWitness_state
#print axioms OPH.CommonWorldInstrumentJoin.committedWitness_state_ne_runState
#print axioms OPH.CommonWorldInstrumentJoin.record_partition_proj_ne_record_effect
#print axioms OPH.CommonWorldInstrumentJoin.committedEffectPair_diagonal
#print axioms OPH.CommonWorldInstrumentJoin.diagonal_refines_record
#print axioms OPH.CommonWorldInstrumentJoin.record_pinching_absorbs_diagonal
#print axioms OPH.CommonWorldInstrumentJoin.lueders_diagonal_channel_eq_pinching
#print axioms OPH.CommonWorldInstrumentJoin.diagonal_pinching_runState
#print axioms OPH.CommonWorldInstrumentJoin.diagonal_pinching_witnessState
#print axioms OPH.CommonWorldInstrumentJoin.witnessState_record_certain
#print axioms OPH.CommonWorldInstrumentJoin.swapTwisted_diagonal_channel_runState
#print axioms OPH.CommonWorldInstrumentJoin.frameWorldline_injective
#print axioms OPH.CommonWorldInstrumentJoin.negative_step_not_causal
#print axioms OPH.CommonWorldInstrumentJoin.runRepairLaw_ref_eq_runState
#print axioms OPH.CommonWorldInstrumentJoin.carriesCertifiedEvolution
#print axioms OPH.CommonWorldInstrumentJoin.carriesInstrumentJoin
#print axioms OPH.CommonWorldInstrumentJoin.stepClock_eq_exp
#print axioms OPH.CommonWorldInstrumentJoin.stepClock_periodic_iff
#print axioms OPH.CommonWorldInstrumentJoin.stepEvent_causal_mono
#print axioms OPH.CommonWorldInstrumentJoin.stepEvent_injective
#print axioms OPH.CommonWorldInstrumentJoin.carriesStepDictionary
#print axioms OPH.CommonWorldInstrumentJoin.register_modularFlowData
#print axioms OPH.CommonWorldInstrumentJoin.register_first_law
#print axioms OPH.CommonWorldInstrumentJoin.carriesFirstLawJoin
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommonWorld_receipt
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommittedWitness_extends
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommonWorld_inhabited
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommittedWitness_receipt
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommittedWitness_contrast
#print axioms OPH.CommonWorldInstrumentJoin.instrumentedCommittedWitness_state_contrast
#print axioms OPH.CommonWorldInstrumentJoin.twistedExtension_breaks_pinching
#print axioms OPH.CommonWorldInstrumentJoin.pinching_clause_excludes_swapTwisted
#print axioms OPH.CommonWorldInstrumentJoin.zeroStepExtension_collapses
#print axioms OPH.CommonWorldInstrumentJoin.negativeStepExtension_breaks_causal
#print axioms OPH.CommonWorldInstrumentJoin.declared_dictionaries_not_forced
