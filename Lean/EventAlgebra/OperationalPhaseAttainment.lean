import EventAlgebra.PhaseInstrumentDetermination

set_option autoImplicit false

namespace EventAlgebra

/-!
# Declared phase-effect count-table conformance under the PR-04 decision

This module inhabits the committed static fit structure
`IdealPhasePOVMCountModel` of `EventAlgebra.OperationalPhaseInstrument`
with the literals of the committed receipt
`code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json`, and
proves exact algebraic conformance laws for those literals.  Register row
PR-04 carries a recorded decision (dated 2026-08-18, lane issue 730) that
declares the committed exact lift
`I/2 - (2 sqrt 3 / 3) i (QP - PQ)` as a phase-sensitive effect, equal to
the Pauli +Y projector.  A matrix effect is not a quantum operation or an
instrument.  In particular, this file provides no completely-positive
outcome maps, trace-preserving summed channel, postmeasurement update, or
effect/readback compatibility theorem.

**Premise consumption, named per the register.**

* PR-04 (consumed, axiomatize disposition): the phase effect of the static
  fit inhabitant is the declared Pauli +Y projector `sourcePhaseLift`.
  Nothing here derives that effect from source data.
* PR-02 (committed): the record-diagonal state of the committed run
  literals `(111, 68)` is read on the committed algebra-state surface.
* PR-03 remains the separate cross-context valuation-additivity premise.
  Per-context binary POVM completeness does not supply it.
* PR-64 and PR-65 remain open for an operational instrument realization
  and source-produced preparation/readback provenance, respectively.

**Count semantics.**  The count literals are a deterministic semantic
conformance fixture: an external calculator evaluates declared Born weights
and scales the rational values to integer pairs.  They are neither sampled
data nor validation.  Inside Lean, `born_matches` is already a field of
`IdealPhasePOVMCountModel`; the theorem
`attained_counts_cross_multiply_born_fit` only cross-multiplies that assumed
fit equation.  It proves no generation, minimality, custody, producer
identity, or binding to the external receipt.  The independent Python
verifier checks the external fixture and its least-multiple rule separately.

**Exact claim.**  The static inhabitant lives over the committed
algebraic core, carries the declared effect in its phase slot and the
declared diagonal state as its model preparation, reproduces the committed
run
literals in its diagonal context, satisfies the exact Born-consistency
law in every context, satisfies the cross-multiplied fit law in every
context, has every count mass a positive integer multiple of the
committed run mass, and has phase counts `(179, 179)` inside the
committed integer receipt window of
`EventAlgebra.PhaseInstrumentDetermination`.  One composed receipt
(`operationalPhaseAttainment_receipt`) derives all clauses from the
typed phase-effect decision bundle.  The theorem and several definitions
retain historical `Attainment` names for API compatibility only.

**Falsifier.**  Any failure of the exact equalities against the
committed literals: a Born weight differing from the receipt fraction, a
count pair failing the scaled-weight law, a phase count pair outside the
integer window, or a phase effect differing from the committed lift.
The Python verifier of the receipt checks the same equalities
independently from the pinned payload.

**Nonclaims and boundary.**  The inhabitant is not an operational,
emergent, or laboratory instrument, its counts are not a preregistered
statistical run, and the
physical-attachment premises of the register (row PR-52 among them) stay
open.  The synthetic-inhabitant boundary of the base module is
unchanged: the type system cannot distinguish this inhabitant from a
synthetic one, and the declared-effect decision lives outside the type
system.
Register row PR-04 is consumed, never discharged: the row remains on the
register under its axiomatize disposition.  Operational realization and
source provenance remain open under PR-64 and PR-65.
-/

open Matrix
open OPH.QFT
open scoped ComplexOrder

noncomputable section

/-! ## The declared diagonal matrix and generated count literals -/

/-- The declared record-diagonal representation: the diagonal matrix built
from the committed run literals `(111, 68)` at mass `179` on the
committed algebra-state surface (register row PR-02).  Its matrix coincides
with the base module's `syntheticIdealPrep`.  It has the same diagonal shape
as states delimited by `EventAlgebra.SourceReachability`, but no reachability
theorem applies: that module's enumerated carriers do not include `Fin 2`,
and its implication is reachability to diagonality, not the converse. -/
def committedRunState : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![(111 / 179 : ℂ), (68 / 179 : ℂ)]

theorem committedRunState_eq_synthetic :
    committedRunState = syntheticIdealPrep := rfl

theorem committedRunState_isState : IsState committedRunState :=
  syntheticIdealPrep_isState

/-- The committed run state has vanishing upper off-diagonal entry: the
free coordinate of `PhaseInstrumentDetermination.prep_coordinates` is
zero for this preparation. -/
theorem committedRunState_offdiag_zero : committedRunState 0 1 = 0 :=
  Matrix.diagonal_apply_ne _ (by decide)

/-- The phase-context count literals of the committed fixture: the exact
Born weight `1/2` of the declared effect under the declared diagonal
state, scaled to `358 = 2 * 179`, the least positive multiple of the
committed run mass clearing the denominator. -/
def producedPhaseCounts : ℕ × ℕ := (179, 179)

/-- The rotated-context count literals of the committed receipt: labels
`0` and `1` conjugate to the record projector (weight `111/179` at mass
`179`), labels `2` through `5` carry the rotated weight `315/716` at
mass `716 = 4 * 179`. -/
def producedRotatedCounts : Fin 6 → ℕ × ℕ :=
  ![(111, 68), (111, 68), (315, 401), (315, 401), (315, 401), (315, 401)]

theorem producedRotatedCounts_eq_synthetic :
    producedRotatedCounts = syntheticRotatedCounts := rfl

/-- The phase fit equation: the Born weight of the declared effect under
the declared diagonal state is exactly the fixture frequency. -/
theorem produced_born_phase :
    bornWeight committedRunState sourcePhaseLift =
      binaryFrequency producedPhaseCounts := by
  rw [committedRunState_eq_synthetic, sourcePhaseLift_eq_rhoYPlus]
  norm_num [bornWeight, syntheticIdealPrep, rhoYPlus, binaryFrequency,
    producedPhaseCounts, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    Fin.sum_univ_two]

/-! ## The static fit inhabitant (legacy `attained` names) -/

/-- The fit-data package of the recorded decision: the declared diagonal
state, the rotated and phase count literals, and exact Born-fit equations.
The package carries no count provenance or operational semantics. -/
def recordedDecisionFitData : IdealPhaseFitData where
  prep := committedRunState
  prep_isState := committedRunState_isState
  rotatedCounts := producedRotatedCounts
  rotatedCounts_pos := by intro g; fin_cases g <;> decide
  phaseCounts := producedPhaseCounts
  phaseCounts_pos := by decide
  born_diagonal := by
    rw [committedRunState_eq_synthetic]
    exact synthetic_born_diagonal
  born_rotated := by
    intro g
    rw [committedRunState_eq_synthetic, producedRotatedCounts_eq_synthetic]
    exact synthetic_born_rotated g
  born_phase := produced_born_phase

/-- The static fit inhabitant assembled from the recorded decision's data.
The historical name does not assert operational attainment. -/
def attainedModel : IdealPhasePOVMCountModel :=
  modelOfIdealPhaseFitData recordedDecisionFitData

/-! ## Pinned coordinates of the attained inhabitant -/

theorem attainedModel_core :
    attainedModel.core = QuantumSurface.committedAlgebraicPhaseCompletion :=
  rfl

/-- The phase slot carries the declared phase effect: the committed exact
lift, equal to the Pauli +Y projector. -/
theorem attainedModel_phase_effect :
    attainedModel.effect InstrumentContext.phase 0 = sourcePhaseLift := by
  show committedEffectPair InstrumentContext.phase 0 = sourcePhaseLift
  simp [committedEffectPair, committedContextEffect]

theorem attainedModel_prep_eq : attainedModel.prep = committedRunState := rfl

theorem attainedModel_phase_counts :
    attainedModel.counts InstrumentContext.phase = (179, 179) := rfl

theorem attainedModel_rotated_counts (g : Fin 6) :
    attainedModel.counts (InstrumentContext.web (WebContext.conjugated g)) =
      producedRotatedCounts g := rfl

theorem attainedModel_diagonal_counts :
    attainedModel.counts (InstrumentContext.web WebContext.diagonal) =
      (111, 68) :=
  attainedModel.diagonal_counts_run

/-! ## The exact laws of the receipt -/

/-- **Exact Born consistency.**  In every context the fixture frequency of
the static inhabitant equals the Born weight of the context effect
under the committed run state.  This is the structure's fit law
specialized to the receipt literals. -/
theorem attained_born_consistency (c : InstrumentContext) :
    bornWeight attainedModel.prep (attainedModel.effect c 0) =
      binaryFrequency (attainedModel.counts c) :=
  attainedModel.born_matches c

/-- Cross-multiplication of the model's assumed exact Born-fit field.  This
does not establish how the literals were generated, that their masses are
minimal, or that they are bound to an external producer receipt. -/
theorem attained_counts_cross_multiply_born_fit (c : InstrumentContext) :
    ((attainedModel.counts c).1 : ℂ) =
      bornWeight attainedModel.prep (attainedModel.effect c 0) *
        (((attainedModel.counts c).1 + (attainedModel.counts c).2 : ℕ) : ℂ) := by
  rw [attainedModel.born_matches c, binaryFrequency]
  have h : (((attainedModel.counts c).1 + (attainedModel.counts c).2 : ℕ) : ℂ) ≠ 0 := by
    exact_mod_cast (attainedModel.counts_pos c).ne'
  exact (div_mul_cancel₀ _ h).symm

/-- Every count mass of the static inhabitant is a positive integer
multiple of the committed run mass `179`: the diagonal and
record-conjugate contexts carry mass `179`, the rotated contexts carry
`716 = 4 * 179`, and the phase context carries `358 = 2 * 179`. -/
theorem attained_masses (c : InstrumentContext) :
    ∃ k : ℕ, 0 < k ∧
      (attainedModel.counts c).1 + (attainedModel.counts c).2 = k * 179 := by
  cases c with
  | web wc =>
      cases wc with
      | diagonal =>
          exact ⟨1, one_pos, by rw [attainedModel.diagonal_counts_run]; decide⟩
      | conjugated g =>
          fin_cases g
          · exact ⟨1, one_pos, by rw [attainedModel_rotated_counts]; decide⟩
          · exact ⟨1, one_pos, by rw [attainedModel_rotated_counts]; decide⟩
          · exact ⟨4, by norm_num, by rw [attainedModel_rotated_counts]; decide⟩
          · exact ⟨4, by norm_num, by rw [attainedModel_rotated_counts]; decide⟩
          · exact ⟨4, by norm_num, by rw [attainedModel_rotated_counts]; decide⟩
          · exact ⟨4, by norm_num, by rw [attainedModel_rotated_counts]; decide⟩
  | phase =>
      exact ⟨2, by norm_num, by rw [attainedModel_phase_counts]; decide⟩

/-- The phase frequency of the static inhabitant in closed form. -/
theorem attained_phase_frequency :
    binaryFrequency (attainedModel.counts InstrumentContext.phase) = 1 / 2 := by
  rw [attainedModel_phase_counts]
  norm_num [binaryFrequency]

/-- Consistency with the determination module: the closed-form frequency
law `phase_frequency_eq` instantiated at the static inhabitant.  The
committed run state has zero off-diagonal imaginary part, so the law
forces frequency `1/2`, exactly the produced literal. -/
theorem attained_phase_frequency_determined :
    binaryFrequency (attainedModel.counts InstrumentContext.phase) =
      1 / 2 - ((attainedModel.prep 0 1).im : ℂ) :=
  phase_frequency_eq attainedModel rfl

/-- The produced phase counts satisfy the committed integer receipt
window, derived from the general window theorem of the determination
module rather than recomputed. -/
theorem attained_phase_counts_window :
    (32041 : ℤ) * (((attainedModel.counts InstrumentContext.phase).1 : ℤ) -
        ((attainedModel.counts InstrumentContext.phase).2 : ℤ)) ^ 2 ≤
      30192 * (((attainedModel.counts InstrumentContext.phase).1 : ℤ) +
        ((attainedModel.counts InstrumentContext.phase).2 : ℤ)) ^ 2 :=
  phase_counts_receipt_window attainedModel rfl

/-- The same window at the receipt literals, checked by literal
arithmetic: `32041 * 0 ≤ 30192 * 358 ^ 2`. -/
theorem attained_phase_counts_window_literal :
    (32041 : ℤ) * ((179 : ℤ) - 179) ^ 2 ≤ 30192 * ((179 : ℤ) + 179) ^ 2 := by
  norm_num

/-! ## The phase-effect decision bundle and composed conformance theorem -/

/-- **The antecedent bundle of the recorded decision.**  Its fields name
the consumed register rows: the declared effect is the axiomatized
PR-04 phase effect (recorded decision dated 2026-08-18), pinned to
the committed exact lift, and the declared preparation is the
record-diagonal run state on the committed PR-02 algebra-state surface.
The bundle carries the decision's declared data; it certifies no source
production. -/
structure PhaseEffectDecisionBundle where
  /-- The declared phase effect: register row PR-04 under its
  axiomatize disposition. -/
  effect : Matrix (Fin 2) (Fin 2) ℂ
  /-- The recorded decision pins the effect to the committed exact
  lift. -/
  effect_eq_lift : effect = sourcePhaseLift
  /-- The declared common preparation on the committed record surface
  (register row PR-02). -/
  prep : Matrix (Fin 2) (Fin 2) ℂ
  /-- The recorded decision pins the preparation to the committed run
  state. -/
  prep_eq_run_state : prep = committedRunState

/-- The committed instance of the bundle: the recorded decision's
declared data.  Nonvacuity of the antecedent. -/
def recordedDecisionBundle : PhaseEffectDecisionBundle where
  effect := sourcePhaseLift
  effect_eq_lift := rfl
  prep := committedRunState
  prep_eq_run_state := rfl

/-- **The composed static conformance theorem (legacy OL-C5 name).**  From the
one antecedent bundle of the recorded decision, the conjunction:

1. the static inhabitant lives over the committed algebraic core;
2. its phase slot carries the bundle's declared effect;
3. its preparation is the bundle's declared state, with vanishing upper
   off-diagonal entry;
4. its diagonal counts are the committed run literals `(111, 68)`;
5. exact Born consistency holds in every context;
6. the cross-multiplied assumed Born-fit law holds in every context;
7. every count mass is a positive integer multiple of the committed run
   mass `179`;
8. the phase counts are `(179, 179)` and satisfy the committed integer
   receipt window.

Boundary: this legacy-named receipt certifies static algebraic conformance
only.  It does not certify an operational instrument, generation semantics,
receipt custody, or laboratory validation. -/
theorem operationalPhaseAttainment_receipt (D : PhaseEffectDecisionBundle) :
    attainedModel.core = QuantumSurface.committedAlgebraicPhaseCompletion ∧
    attainedModel.effect InstrumentContext.phase 0 = D.effect ∧
    attainedModel.prep = D.prep ∧
    D.prep 0 1 = 0 ∧
    attainedModel.counts (InstrumentContext.web WebContext.diagonal) =
      (111, 68) ∧
    (∀ c : InstrumentContext,
      bornWeight attainedModel.prep (attainedModel.effect c 0) =
        binaryFrequency (attainedModel.counts c)) ∧
    (∀ c : InstrumentContext,
      ((attainedModel.counts c).1 : ℂ) =
        bornWeight attainedModel.prep (attainedModel.effect c 0) *
          (((attainedModel.counts c).1 + (attainedModel.counts c).2 : ℕ) : ℂ)) ∧
    (∀ c : InstrumentContext, ∃ k : ℕ, 0 < k ∧
      (attainedModel.counts c).1 + (attainedModel.counts c).2 = k * 179) ∧
    attainedModel.counts InstrumentContext.phase = (179, 179) ∧
    (32041 : ℤ) * (((attainedModel.counts InstrumentContext.phase).1 : ℤ) -
        ((attainedModel.counts InstrumentContext.phase).2 : ℤ)) ^ 2 ≤
      30192 * (((attainedModel.counts InstrumentContext.phase).1 : ℤ) +
        ((attainedModel.counts InstrumentContext.phase).2 : ℤ)) ^ 2 := by
  refine ⟨rfl, ?_, ?_, ?_, attainedModel.diagonal_counts_run,
    attained_born_consistency, attained_counts_cross_multiply_born_fit,
    attained_masses, attainedModel_phase_counts,
    attained_phase_counts_window⟩
  · rw [D.effect_eq_lift]
    exact attainedModel_phase_effect
  · rw [D.prep_eq_run_state]
    exact attainedModel_prep_eq
  · rw [D.prep_eq_run_state]
    exact committedRunState_offdiag_zero

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms committedRunState_eq_synthetic
#print axioms committedRunState_isState
#print axioms committedRunState_offdiag_zero
#print axioms producedRotatedCounts_eq_synthetic
#print axioms produced_born_phase
#print axioms recordedDecisionFitData
#print axioms attainedModel
#print axioms attainedModel_core
#print axioms attainedModel_phase_effect
#print axioms attainedModel_prep_eq
#print axioms attainedModel_phase_counts
#print axioms attainedModel_rotated_counts
#print axioms attainedModel_diagonal_counts
#print axioms attained_born_consistency
#print axioms attained_counts_cross_multiply_born_fit
#print axioms attained_masses
#print axioms attained_phase_frequency
#print axioms attained_phase_frequency_determined
#print axioms attained_phase_counts_window
#print axioms attained_phase_counts_window_literal
#print axioms recordedDecisionBundle
#print axioms operationalPhaseAttainment_receipt

end EventAlgebra
