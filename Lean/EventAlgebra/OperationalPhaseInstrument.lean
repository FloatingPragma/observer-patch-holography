import EventAlgebra.QuantumAdequacySurface

set_option autoImplicit false

namespace EventAlgebra

/-!
# Typed operational PR-04 instrument semantics (OL-C5, issue #730)

This module types the operational content that register row PR-04 owes and
proves exactly three things about it, all premise-free finite mathematics on
the committed `Fin 2` record algebra:

* **The instrument type** (`OperationalPhaseInstrument`): one structure whose
  fields are the algebraic phase-completion core (defaulting to the committed
  `QuantumSurface.committedAlgebraicPhaseCompletion`), one declared
  preparation shared by every context, the finite context family
  `InstrumentContext` (the committed web contexts of `QFT.SourceContextWeb`
  plus one designated phase context), per-context binary outcome effects
  forming an exact POVM, an explicit phase-sensitivity clause, realized
  outcome counts per context with the diagonal counts pinned to the committed
  run literals `(111, 68)`, and one exact Born-consistency law tying every
  context's realized frequency to the same preparation.
* **Conditional sufficiency** (`operational_instrument_completes_tomography`,
  `preparation_identified_by_realized_data`): every inhabitant of the
  structure yields fixed-trace state identification, by composing the
  instrument fields with the committed `completes` clause of the algebraic
  core; in particular the preparation is the unique state reproducing the
  realized frequencies.  PR-04 discharge therefore reduces to inhabiting the
  structure.
* **Exact necessity** (`diagonal_context_insufficient`,
  `diagonal_receipts_do_not_identify`, `rotated_context_separates_x_pair`,
  `phase_effect_not_diagonal`, `phase_effect_genuinely_complex`,
  `phase_context_necessary`): two distinct certified states
  (`rhoXPlus`, `rhoXMinus`) carry identical Born data on every diagonal
  effect, so the diagonal context alone identifies no state; the committed
  rotated context separates that pair, and the committed Pauli-Y pair (the
  no-real-closure clause `QuantumSurface.real_closure_blind`) defeats even
  the full phase-free real closure, so the phase-sensitivity field of the
  instrument can be met by no diagonal effect and by no complexified real
  effect.  The demanded phase context is necessary, not decorative.

`partialWitness` records exactly which fields are inhabited today from
committed objects: the algebraic core, every context effect matrix, the
POVM certificates, and the diagonal outcome counts — and consumes the
missing data (`MissingPhaseReceipts`: a preparation, rotated and phase
outcome counts, and the Born-consistency receipts) as explicit hypotheses.
`DischargesPR04` is the single reduction Prop, proved equivalent to
`Nonempty MissingPhaseReceipts` (`dischargesPR04_iff_missing_receipts`).

**What is NOT proved.**  PR-04 stays open.  No declaration in this file or
elsewhere in the repository inhabits `OperationalPhaseInstrument`,
`MissingPhaseReceipts`, or `DischargesPR04`: the committed payload supplies
no preparation declaration, no rotated or phase outcome counts, and no
common-preparation validation, and none are constructed here.  An inhabitant
would certify internal coherence of supplied data — one preparation, exact
POVMs, exact Born match — while the provenance of the supplied counts (that
they are producer receipts of a committed simulator run rather than chosen
numbers) is a custody condition the register tracks outside the type system;
Lean cannot enforce it, and the register's discharge discipline prohibits
synthetic counts.  No frame valuation is consumed: register rows PR-02 and
PR-03 do not enter, and no probability reading beyond the committed diagonal
frequencies is asserted.  This file supplies typed operational semantics, a
proved reduction, and a proved necessity — not an operational attainment.
-/

open Matrix
open OPH.QFT
open scoped ComplexOrder

noncomputable section

/-! ## The committed web contexts, repackaged -/

/-- The record-side real projector of each committed web context: the
diagonal context carries `recordProjector`, the conjugated context of gauge
label `g` carries `conjProjector g`.  Committed adapters of
`QFT.SourceContextWeb`, repackaged as one function on `WebContext`. -/
def webContextProjector : WebContext → Matrix (Fin 2) (Fin 2) ℝ
  | WebContext.diagonal => recordProjector
  | WebContext.conjugated g => conjProjector g

theorem webContextProjector_diagonal :
    webContextProjector WebContext.diagonal = recordProjector := rfl

theorem webContextProjector_conjugated (g : Fin 6) :
    webContextProjector (WebContext.conjugated g) = conjProjector g := rfl

/-- Every committed web-context projector is a symmetric real matrix. -/
theorem webContextProjector_symm (c : WebContext) :
    (webContextProjector c)ᵀ = webContextProjector c := by
  cases c with
  | diagonal =>
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [webContextProjector, recordProjector]
  | conjugated g =>
      have h := conjProjector_offdiag_symmetric g
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [webContextProjector, Matrix.transpose_apply, h]

/-- Every committed web-context projector is idempotent. -/
theorem webContextProjector_idem (c : WebContext) :
    webContextProjector c * webContextProjector c = webContextProjector c := by
  cases c with
  | diagonal => exact recordProjector_idem
  | conjugated g => exact conjProjector_idem g

/-- The complexification of a symmetric real matrix is Hermitian. -/
theorem complexifyRealMatrix_isHermitian_of_symmetric
    {E : Matrix (Fin 2) (Fin 2) ℝ} (hE : Eᵀ = E) :
    (complexifyRealMatrix E).IsHermitian := by
  show (complexifyRealMatrix E)ᴴ = complexifyRealMatrix E
  ext i j
  have h : E j i = E i j := by
    have h' := congrFun (congrFun hE i) j
    simpa [Matrix.transpose_apply] using h'
  simp [complexifyRealMatrix, Matrix.conjTranspose_apply,
    Complex.conj_ofReal, h]

/-- The complexified web-context projectors are projection events of the
represented algebra. -/
theorem webContextProjector_isEvent (c : WebContext) :
    IsEvent (complexifyRealMatrix (webContextProjector c)) := by
  refine ⟨complexifyRealMatrix_isHermitian_of_symmetric
    (webContextProjector_symm c), ?_⟩
  rw [← complexifyRealMatrix_mul, webContextProjector_idem]

/-! ## The instrument context family -/

/-- The context family of the typed PR-04 instrument: every committed web
context (the diagonal record/companion context and the six conjugated
contexts of `QFT.SourceContextWeb`) together with one designated
phase-sensitive context. -/
inductive InstrumentContext : Type where
  | web : WebContext → InstrumentContext
  | phase : InstrumentContext
  deriving DecidableEq

/-- The eight contexts of the family, enumerated. -/
def enumerateContexts : Fin 8 → InstrumentContext :=
  ![InstrumentContext.web WebContext.diagonal,
    InstrumentContext.web (WebContext.conjugated 0),
    InstrumentContext.web (WebContext.conjugated 1),
    InstrumentContext.web (WebContext.conjugated 2),
    InstrumentContext.web (WebContext.conjugated 3),
    InstrumentContext.web (WebContext.conjugated 4),
    InstrumentContext.web (WebContext.conjugated 5),
    InstrumentContext.phase]

/-- The context family is finite: eight contexts exhaust it. -/
theorem enumerateContexts_surjective : Function.Surjective enumerateContexts := by
  intro c
  cases c with
  | web wc =>
      cases wc with
      | diagonal => exact ⟨0, rfl⟩
      | conjugated g =>
          fin_cases g
          · exact ⟨1, rfl⟩
          · exact ⟨2, rfl⟩
          · exact ⟨3, rfl⟩
          · exact ⟨4, rfl⟩
          · exact ⟨5, rfl⟩
          · exact ⟨6, rfl⟩
  | phase => exact ⟨7, rfl⟩

/-- The outcome-`0` frequency of a binary count pair, as a complex scalar:
`a / (a + b)`. -/
def binaryFrequency (p : ℕ × ℕ) : ℂ := (p.1 : ℂ) / ((p.1 + p.2 : ℕ) : ℂ)

/-- The committed diagonal run frequency in closed form. -/
theorem binaryFrequency_diagonal_run :
    binaryFrequency (111, 68) = 111 / 179 := by
  norm_num [binaryFrequency]

/-! ## (P1) The typed operational instrument -/

/-- **The typed operational PR-04 instrument.**  One inhabitant of this
structure is exactly what the register row PR-04 declares missing: a
source-attached phase-sensitive instrument with common-preparation outcome
receipts on the committed `Fin 2` record algebra.  The diagonal fields are
inhabitable from committed objects (`partialWitness`); the preparation, the
rotated and phase outcome counts, and every Born-consistency receipt are
exactly what the committed payload does not supply.  No field is asserted to
hold of the committed run; the structure is a demand, not a record. -/
structure OperationalPhaseInstrument where
  /-- (i) The algebraic core: an algebraic phase completion in the sense of
  the PR-04 bundle, defaulting to the committed one.  Its `lift` is the
  declared phase-sensitive projection event; its `separates` and `completes`
  clauses are the committed algebraic receipts this instrument composes
  with. -/
  core : AlgebraicPhaseCompletion := QuantumSurface.committedAlgebraicPhaseCompletion
  /-- (ii) The declared preparation: one state reference shared by every
  context.  The committed payload declares no preparation; this field is
  part of what PR-04 owes. -/
  prep : Matrix (Fin 2) (Fin 2) ℂ
  /-- The declared preparation is a state: positive semidefinite with unit
  trace. -/
  prep_isState : IsState prep
  /-- (iv) Per-context binary outcome effects on the represented
  dimension-two algebra. -/
  effect : InstrumentContext → Fin 2 → Matrix (Fin 2) (Fin 2) ℂ
  /-- POVM nonnegativity: every outcome effect is positive semidefinite. -/
  effect_posSemidef : ∀ (c : InstrumentContext) (i : Fin 2),
    (effect c i).PosSemidef
  /-- POVM normalisation: in every context the two outcome effects sum
  exactly to the sure effect. -/
  effect_complete : ∀ c : InstrumentContext, effect c 0 + effect c 1 = 1
  /-- The web contexts carry the committed complexified projectors: outcome
  `0` of the diagonal context is the record projector, outcome `0` of the
  conjugated context of label `g` is `conjProjector g`.  Outcome `1` is
  forced by `effect_complete` (`effect_one`). -/
  web_effect : ∀ c : WebContext,
    effect (InstrumentContext.web c) 0 =
      complexifyRealMatrix (webContextProjector c)
  /-- The phase context carries the core's declared phase lift as its
  outcome-`0` effect. -/
  phase_effect : effect InstrumentContext.phase 0 = core.lift
  /-- The phase-sensitivity clause, stated so it cannot be met vacuously:
  the phase outcome effect assigns different Born weights to the two
  committed Pauli-Y states.  No diagonal effect and no complexification of
  a real matrix satisfies this clause (`phase_effect_not_diagonal`,
  `phase_effect_genuinely_complex`); with `phase_effect` it follows from
  the core's `separates` clause, and it is kept explicit as the exact
  non-diagonality demand on the phase context. -/
  phase_sensitive :
    bornWeight rhoYPlus (effect InstrumentContext.phase 0) ≠
      bornWeight rhoYMinus (effect InstrumentContext.phase 0)
  /-- (v) Realized outcome counts per context.  Only the diagonal entry has
  a committed value today (`diagonal_counts_run`); the rotated and phase
  entries are the missing receipts named by PR-04. -/
  counts : InstrumentContext → ℕ × ℕ
  /-- Every context carries positive realized outcome mass: no context of
  the family may be left without data. -/
  counts_pos : ∀ c : InstrumentContext, 0 < (counts c).1 + (counts c).2
  /-- The diagonal counts are the committed run literals: the two realized
  cells `(record 18, companion 13) = 111` and `(record 18, companion 3) =
  68` of `QFT.SourceContextWeb` (`realizedOutcomeCounts`). -/
  diagonal_counts_run :
    counts (InstrumentContext.web WebContext.diagonal) = (111, 68)
  /-- (vi) The Born-consistency law and the common-preparation clause in one
  field: for every context of the family, the realized outcome-`0` frequency
  equals the Born weight of that context's outcome-`0` effect under the one
  declared preparation `prep`.  A single `prep` appears for all contexts;
  this is the common-preparation validation PR-04 owes.  Outcome `1` is
  forced (`born_matches_snd`). -/
  born_matches : ∀ c : InstrumentContext,
    bornWeight prep (effect c 0) = binaryFrequency (counts c)

namespace OperationalPhaseInstrument

/-- Outcome `1` of every context is the complement effect, forced by POVM
normalisation. -/
theorem effect_one (I : OperationalPhaseInstrument) (c : InstrumentContext) :
    I.effect c 1 = 1 - I.effect c 0 := by
  have h := I.effect_complete c
  rw [← h]
  abel

/-- Every outcome effect is an effect in the represented-algebra sense:
`0 ≤ E ≤ 1` in the positive-semidefinite order.  The exact POVM fields
compose to the standard effect certificate. -/
theorem effect_isEffect (I : OperationalPhaseInstrument)
    (c : InstrumentContext) : ∀ i : Fin 2, IsEffect (I.effect c i) := by
  rw [Fin.forall_fin_two]
  refine ⟨⟨I.effect_posSemidef c 0, ?_⟩, ⟨I.effect_posSemidef c 1, ?_⟩⟩
  · rw [← I.effect_one c]
    exact I.effect_posSemidef c 1
  · rw [I.effect_one c]
    simpa using I.effect_posSemidef c 0

/-- The diagonal outcome mass of every inhabitant is the committed run
mass `179`. -/
theorem diagonal_mass (I : OperationalPhaseInstrument) :
    (I.counts (InstrumentContext.web WebContext.diagonal)).1 +
      (I.counts (InstrumentContext.web WebContext.diagonal)).2 = 179 := by
  rw [I.diagonal_counts_run]
  norm_num

/-- Every inhabitant's diagonal counts are exactly the committed realized
outcome counts of the run. -/
theorem diagonal_counts_committed (I : OperationalPhaseInstrument) :
    realizedOutcomeCounts WebContext.diagonal =
      some (I.counts (InstrumentContext.web WebContext.diagonal)) := by
  rw [I.diagonal_counts_run]
  rfl

/-- The outcome-`1` Born-consistency law, forced by normalisation: the
realized outcome-`1` frequency equals the Born weight of the complement
effect under the same preparation. -/
theorem born_matches_snd (I : OperationalPhaseInstrument)
    (c : InstrumentContext) :
    bornWeight I.prep (I.effect c 1) =
      ((I.counts c).2 : ℂ) /
        (((I.counts c).1 + (I.counts c).2 : ℕ) : ℂ) := by
  have hsum : (((I.counts c).1 + (I.counts c).2 : ℕ) : ℂ) =
      ((I.counts c).1 : ℂ) + ((I.counts c).2 : ℂ) := by
    push_cast
    ring
  have hs : ((I.counts c).1 : ℂ) + ((I.counts c).2 : ℂ) ≠ 0 := by
    rw [← hsum]
    exact_mod_cast (I.counts_pos c).ne'
  rw [I.effect_one c, bornWeight_sub, bornWeight_one I.prep_isState,
    I.born_matches c, binaryFrequency, hsum]
  field_simp
  ring

end OperationalPhaseInstrument

/-! ## (P2) Conditional sufficiency: an inhabitant completes tomography -/

/-- **Conditional sufficiency.**  Any inhabitant of the typed instrument
yields fixed-trace state identification: two matrices of equal trace that
agree in Born weight on every outcome effect of the instrument are equal.
The proof composes the pinned web and phase effects with the committed
`completes` clause of the algebraic core; the instrument adds no separating
power beyond the core and removes none.  PR-04 discharge reduces to
inhabiting `OperationalPhaseInstrument`; this theorem is the reduction's
sufficiency direction and consumes no committed outcome data beyond the
inhabitant's own fields. -/
theorem operational_instrument_completes_tomography
    (I : OperationalPhaseInstrument)
    {rho sigma : Matrix (Fin 2) (Fin 2) ℂ}
    (htr : Matrix.trace rho = Matrix.trace sigma)
    (h : ∀ (c : InstrumentContext) (i : Fin 2),
      bornWeight rho (I.effect c i) = bornWeight sigma (I.effect c i)) :
    rho = sigma := by
  have hrec := h (InstrumentContext.web WebContext.diagonal) 0
  have hrot := h (InstrumentContext.web (WebContext.conjugated 3)) 0
  have hlift := h InstrumentContext.phase 0
  rw [I.web_effect WebContext.diagonal, webContextProjector_diagonal] at hrec
  rw [I.web_effect (WebContext.conjugated 3),
    webContextProjector_conjugated] at hrot
  rw [I.phase_effect] at hlift
  exact I.core.completes rho sigma htr hrec hrot hlift

/-- **Operational identification of the preparation.**  For any inhabitant,
the declared preparation is the unique state whose Born weights reproduce
the realized outcome frequencies of every context: a state matching the
realized data is the preparation.  This is the exact sense in which an
inhabitant would supply operational state identification from realized
receipts. -/
theorem preparation_identified_by_realized_data
    (I : OperationalPhaseInstrument)
    {sigma : Matrix (Fin 2) (Fin 2) ℂ} (hsigma : IsState sigma)
    (hfreq : ∀ c : InstrumentContext,
      bornWeight sigma (I.effect c 0) = binaryFrequency (I.counts c)) :
    sigma = I.prep := by
  have hAll : ∀ (c : InstrumentContext) (i : Fin 2),
      bornWeight sigma (I.effect c i) = bornWeight I.prep (I.effect c i) := by
    intro c
    rw [Fin.forall_fin_two]
    have h0 : bornWeight sigma (I.effect c 0) =
        bornWeight I.prep (I.effect c 0) :=
      (hfreq c).trans (I.born_matches c).symm
    refine ⟨h0, ?_⟩
    rw [I.effect_one c, bornWeight_sub, bornWeight_sub,
      bornWeight_one hsigma, bornWeight_one I.prep_isState, h0]
  exact operational_instrument_completes_tomography I
    (hsigma.2.trans I.prep_isState.2.symm) hAll

/-! ## (P3) Exact necessity: the diagonal context identifies no state

Two distinct certified states with equal diagonals witness that
diagonal-only outcome families separate nothing: the committed native
context's Born data are blind to every off-diagonal coordinate.  The
committed rotated context separates this pair, and the committed Pauli-Y
pair defeats the entire phase-free real closure, so each escalation step of
the instrument's context family is necessary. -/

/-- Projection onto the `+X` axis: `(1/2)(1 + σ_x)`. -/
def sigmaXPlusProj : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1 / 2, 1 / 2; 1 / 2, 1 / 2]

/-- Projection onto the `-X` axis: `(1/2)(1 - σ_x)`. -/
def sigmaXMinusProj : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1 / 2, -(1 / 2); -(1 / 2), 1 / 2]

theorem sigmaXPlusProj_isEvent : IsEvent sigmaXPlusProj := by
  constructor
  · show sigmaXPlusProjᴴ = sigmaXPlusProj
    ext i j
    fin_cases i <;> fin_cases j <;> simp [sigmaXPlusProj]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [sigmaXPlusProj, Matrix.mul_apply, Fin.sum_univ_two]

theorem sigmaXMinusProj_isEvent : IsEvent sigmaXMinusProj := by
  constructor
  · show sigmaXMinusProjᴴ = sigmaXMinusProj
    ext i j
    fin_cases i <;> fin_cases j <;> simp [sigmaXMinusProj]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [sigmaXMinusProj, Matrix.mul_apply, Fin.sum_univ_two]

/-- Mixed state polarized toward `+X`: `(1/2)(1 + (1/2) σ_x)`. -/
def rhoXPlus : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1 / 2, 1 / 4; 1 / 4, 1 / 2]

/-- Mixed state polarized toward `-X`: `(1/2)(1 - (1/2) σ_x)`. -/
def rhoXMinus : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1 / 2, -(1 / 4); -(1 / 4), 1 / 2]

/-- The `+X`-polarized state is the exact convex mixture
`(3/4) P₊ + (1/4) P₋` of the two `X` projections. -/
theorem rhoXPlus_eq_mix :
    rhoXPlus = (3 / 4 : ℝ) • sigmaXPlusProj + (1 / 4 : ℝ) • sigmaXMinusProj := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    · simp [rhoXPlus, sigmaXPlusProj, sigmaXMinusProj, Matrix.add_apply,
        Complex.real_smul]
      ring

/-- The `-X`-polarized state is the exact convex mixture
`(1/4) P₊ + (3/4) P₋`. -/
theorem rhoXMinus_eq_mix :
    rhoXMinus = (1 / 4 : ℝ) • sigmaXPlusProj + (3 / 4 : ℝ) • sigmaXMinusProj := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    · simp [rhoXMinus, sigmaXPlusProj, sigmaXMinusProj, Matrix.add_apply,
        Complex.real_smul]
      ring

/-- `rhoXPlus` is a certified state: positive semidefinite (as a convex
mixture of projections) with unit trace. -/
theorem rhoXPlus_isState : IsState rhoXPlus := by
  constructor
  · rw [rhoXPlus_eq_mix]
    exact (posSemidef_real_smul sigmaXPlusProj_isEvent.posSemidef
        (by norm_num)).add
      (posSemidef_real_smul sigmaXMinusProj_isEvent.posSemidef (by norm_num))
  · norm_num [rhoXPlus, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

/-- `rhoXMinus` is a certified state. -/
theorem rhoXMinus_isState : IsState rhoXMinus := by
  constructor
  · rw [rhoXMinus_eq_mix]
    exact (posSemidef_real_smul sigmaXPlusProj_isEvent.posSemidef
        (by norm_num)).add
      (posSemidef_real_smul sigmaXMinusProj_isEvent.posSemidef (by norm_num))
  · norm_num [rhoXMinus, Matrix.trace, Matrix.diag, Fin.sum_univ_two]

/-- The two `X`-polarized states differ exactly in their off-diagonal
entry. -/
theorem rhoXPlus_ne_rhoXMinus : rhoXPlus ≠ rhoXMinus := by
  intro h
  have h01 := congrFun (congrFun h 0) 1
  norm_num [rhoXPlus, rhoXMinus] at h01

/-- The two `X`-polarized states assign equal Born weights to every
diagonal matrix: diagonal effects read only the shared diagonal. -/
theorem rhoX_agree_on_diagonal (E : Matrix (Fin 2) (Fin 2) ℂ)
    (hE : E.IsDiag) :
    bornWeight rhoXPlus E = bornWeight rhoXMinus E := by
  have h01 : E 0 1 = 0 := hE (by decide)
  have h10 : E 1 0 = 0 := hE (by decide)
  simp [bornWeight, rhoXPlus, rhoXMinus, Matrix.trace, Matrix.diag,
    Matrix.mul_apply, Fin.sum_univ_two, h01, h10]

/-- The two committed Pauli-Y states also assign equal Born weights to every
diagonal matrix. -/
theorem rhoY_agree_on_diagonal (E : Matrix (Fin 2) (Fin 2) ℂ)
    (hE : E.IsDiag) :
    bornWeight rhoYPlus E = bornWeight rhoYMinus E := by
  have h01 : E 0 1 = 0 := hE (by decide)
  have h10 : E 1 0 = 0 := hE (by decide)
  simp [bornWeight, rhoYPlus, rhoYMinus, Matrix.trace, Matrix.diag,
    Matrix.mul_apply, Fin.sum_univ_two, h01, h10]

/-- **(P3) The diagonal context identifies no state.**  Two distinct
certified fixed-trace states — `(1/2)(1 ± (1/2) σ_x)`, positive semidefinite
by `rhoXPlus_isState` and `rhoXMinus_isState` — carry identical Born data on
every diagonal effect.  No diagonal-only outcome family separates states, so
the phase-sensitive context demanded by `OperationalPhaseInstrument` is
necessary: realized diagonal receipts, however many, close no tomography. -/
theorem diagonal_context_insufficient :
    ∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      IsState rho ∧ IsState sigma ∧ rho ≠ sigma ∧
        ∀ E : Matrix (Fin 2) (Fin 2) ℂ, E.IsDiag →
          bornWeight rho E = bornWeight sigma E :=
  ⟨rhoXPlus, rhoXMinus, rhoXPlus_isState, rhoXMinus_isState,
    rhoXPlus_ne_rhoXMinus, rhoX_agree_on_diagonal⟩

theorem recordProjector_isDiag : recordProjector.IsDiag := by
  intro i j hij
  fin_cases i <;> fin_cases j <;>
    first
      | exact absurd rfl hij
      | simp [recordProjector]

theorem complexifyRealMatrix_isDiag {E : Matrix (Fin 2) (Fin 2) ℝ}
    (hE : E.IsDiag) : (complexifyRealMatrix E).IsDiag := by
  intro i j hij
  simp [complexifyRealMatrix, hE hij]

/-- The insufficiency instantiated on any inhabitant: both outcomes of the
instrument's own diagonal context assign equal Born data to the two distinct
`X`-polarized states.  The diagonal context of every inhabitant is blind to
this pair; the separating power of an inhabitant lives entirely in its
non-diagonal contexts. -/
theorem diagonal_receipts_do_not_identify (I : OperationalPhaseInstrument) :
    ∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      IsState rho ∧ IsState sigma ∧ rho ≠ sigma ∧
        ∀ i : Fin 2,
          bornWeight rho
              (I.effect (InstrumentContext.web WebContext.diagonal) i) =
            bornWeight sigma
              (I.effect (InstrumentContext.web WebContext.diagonal) i) := by
  have hdiag : (I.effect (InstrumentContext.web WebContext.diagonal) 0).IsDiag := by
    rw [I.web_effect WebContext.diagonal, webContextProjector_diagonal]
    exact complexifyRealMatrix_isDiag recordProjector_isDiag
  refine ⟨rhoXPlus, rhoXMinus, rhoXPlus_isState, rhoXMinus_isState,
    rhoXPlus_ne_rhoXMinus, ?_⟩
  rw [Fin.forall_fin_two]
  refine ⟨rhoX_agree_on_diagonal _ hdiag, ?_⟩
  rw [I.effect_one]
  exact rhoX_agree_on_diagonal _ (Matrix.isDiag_one.sub hdiag)

/-- Born weight of the `+X`-polarized state on the committed rotated
context, in closed form. -/
theorem rhoXPlus_rotated_weight :
    bornWeight rhoXPlus (complexifyRealMatrix (conjProjector 3)) =
      1 / 2 - ((sqrt3 : ℝ) : ℂ) / 8 := by
  rw [conjProjector_three]
  simp [bornWeight, rhoXPlus, complexifyRealMatrix, Matrix.trace,
    Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- Born weight of the `-X`-polarized state on the committed rotated
context, in closed form. -/
theorem rhoXMinus_rotated_weight :
    bornWeight rhoXMinus (complexifyRealMatrix (conjProjector 3)) =
      1 / 2 + ((sqrt3 : ℝ) : ℂ) / 8 := by
  rw [conjProjector_three]
  simp [bornWeight, rhoXMinus, complexifyRealMatrix, Matrix.trace,
    Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- The committed rotated context separates the pair the diagonal context
cannot: the two `X`-polarized states receive Born weights differing by
`sqrt 3 / 4`.  The rotated receipt PR-04 owes carries separating content
beyond every diagonal receipt. -/
theorem rotated_context_separates_x_pair :
    bornWeight rhoXPlus (complexifyRealMatrix (conjProjector 3)) ≠
      bornWeight rhoXMinus (complexifyRealMatrix (conjProjector 3)) := by
  rw [rhoXPlus_rotated_weight, rhoXMinus_rotated_weight]
  intro h
  have hs : ((sqrt3 : ℝ) : ℂ) = 0 := by linear_combination (-4 : ℂ) * h
  have hs' : sqrt3 = 0 := by exact_mod_cast hs
  exact sqrt3_ne_zero hs'

/-- The phase-sensitivity clause is unsatisfiable by diagonal effects: the
phase outcome effect of every inhabitant is genuinely non-diagonal. -/
theorem phase_effect_not_diagonal (I : OperationalPhaseInstrument) :
    ¬ (I.effect InstrumentContext.phase 0).IsDiag := fun hdiag =>
  I.phase_sensitive (rhoY_agree_on_diagonal _ hdiag)

/-- The phase-sensitivity clause is unsatisfiable by every complexified real
matrix: a real candidate would be symmetric (by positivity) and the
committed real-symmetric blindness (`rhoY_bornWeight_equal_on_real_symmetric`)
would contradict phase sensitivity.  This connects the instrument demand to
the committed no-real-closure clause: no phase-free real processing
inhabits the phase context. -/
theorem phase_effect_genuinely_complex (I : OperationalPhaseInstrument) :
    ¬ ∃ E : Matrix (Fin 2) (Fin 2) ℝ,
      complexifyRealMatrix E = I.effect InstrumentContext.phase 0 := by
  rintro ⟨E, hE⟩
  have hherm : (complexifyRealMatrix E).IsHermitian := by
    rw [hE]
    exact (I.effect_posSemidef InstrumentContext.phase 0).isHermitian
  have hsym : Eᵀ = E := by
    ext i j
    rw [Matrix.transpose_apply]
    have h := congrFun (congrFun hherm.eq i) j
    rw [Matrix.conjTranspose_apply] at h
    have h' : ((E j i : ℝ) : ℂ) = ((E i j : ℝ) : ℂ) := by
      simpa [complexifyRealMatrix, Complex.star_def,
        Complex.conj_ofReal] using h
    exact_mod_cast h'
  have hblind := rhoY_bornWeight_equal_on_real_symmetric E hsym
  rw [hE] at hblind
  exact I.phase_sensitive hblind

/-- **The necessity ladder in one statement.**  First conjunct: the diagonal
context alone identifies no state (the `X` pair).  Second conjunct: the
committed no-real-closure clause — even the generous phase-free real closure
identifies no state (the committed Pauli-Y pair).  Each escalation of the
instrument's context family, diagonal to rotated to phase, is forced by an
exact countermodel; the phase-sensitive context of (P1) is necessary, not
decorative. -/
theorem phase_context_necessary :
    (∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      IsState rho ∧ IsState sigma ∧ rho ≠ sigma ∧
        ∀ E : Matrix (Fin 2) (Fin 2) ℂ, E.IsDiag →
          bornWeight rho E = bornWeight sigma E) ∧
    (∃ rho sigma : Matrix (Fin 2) (Fin 2) ℂ,
      IsState rho ∧ IsState sigma ∧ rho ≠ sigma ∧
        ∀ E : Matrix (Fin 2) (Fin 2) ℝ,
          RealSourceEffectClosure E →
            bornWeight rho (complexifyRealMatrix E) =
              bornWeight sigma (complexifyRealMatrix E)) :=
  ⟨diagonal_context_insufficient,
    realSourceEffectClosure_not_tomographically_complete⟩

/-! ## (P4) The partial witness: what is inhabited today -/

/-- The committed outcome-`0` effect of each instrument context: web
contexts carry the committed complexified projectors, the phase context
carries the committed algebraic lift.  Every matrix here is a committed
object; no outcome data accompanies any non-diagonal context. -/
def committedContextEffect : InstrumentContext → Matrix (Fin 2) (Fin 2) ℂ
  | InstrumentContext.web c => complexifyRealMatrix (webContextProjector c)
  | InstrumentContext.phase => sourcePhaseLift

/-- Every committed context effect is a projection event. -/
theorem committedContextEffect_isEvent (c : InstrumentContext) :
    IsEvent (committedContextEffect c) := by
  cases c with
  | web wc => exact webContextProjector_isEvent wc
  | phase => exact sourcePhaseLift_isEvent

/-- The committed binary POVM of each context: the committed effect and its
complement. -/
def committedEffectPair (c : InstrumentContext) :
    Fin 2 → Matrix (Fin 2) (Fin 2) ℂ :=
  ![committedContextEffect c, 1 - committedContextEffect c]

/-- **The missing PR-04 data, typed.**  Exactly the fields of
`OperationalPhaseInstrument` that no committed object supplies: one declared
preparation, realized outcome counts for the rotated and phase contexts,
and the Born-consistency receipts tying every realized frequency — the
committed diagonal one included — to that single preparation.  The diagonal
Born receipt is listed here because the committed run realizes diagonal
counts without declaring any preparation they validate.  Nothing in the
repository inhabits this structure; supplying producer receipts for these
fields is the open content of PR-04. -/
structure MissingPhaseReceipts where
  /-- The declared preparation shared by every context. -/
  prep : Matrix (Fin 2) (Fin 2) ℂ
  /-- The declared preparation is a state. -/
  prep_isState : IsState prep
  /-- Realized outcome counts for each conjugated context.  Labels `0` and
  `1` conjugate to the diagonal projector (`conjProjector_one_eq_record`),
  so their receipts may repeat the diagonal data; labels `2`–`5` are the
  two genuinely rotated contexts. -/
  rotatedCounts : Fin 6 → ℕ × ℕ
  /-- Positive realized mass on every rotated context. -/
  rotatedCounts_pos : ∀ g : Fin 6,
    0 < (rotatedCounts g).1 + (rotatedCounts g).2
  /-- Realized outcome counts for the phase context. -/
  phaseCounts : ℕ × ℕ
  /-- Positive realized mass on the phase context. -/
  phaseCounts_pos : 0 < phaseCounts.1 + phaseCounts.2
  /-- The committed diagonal frequency is the Born weight of the record
  projector under the declared preparation: the common-preparation receipt
  for the one context that already has committed counts. -/
  born_diagonal :
    bornWeight prep (complexifyRealMatrix recordProjector) =
      binaryFrequency (111, 68)
  /-- Each rotated frequency is the Born weight of its conjugated projector
  under the same preparation. -/
  born_rotated : ∀ g : Fin 6,
    bornWeight prep (complexifyRealMatrix (conjProjector g)) =
      binaryFrequency (rotatedCounts g)
  /-- The phase frequency is the Born weight of the committed phase lift
  under the same preparation. -/
  born_phase : bornWeight prep sourcePhaseLift = binaryFrequency phaseCounts

/-- **(P4) The partial witness.**  The fields of the typed instrument that
committed objects inhabit today — the committed algebraic core, every
context effect with its POVM certificates, and the committed diagonal
outcome counts — assembled into a function from the missing data to a full
instrument.  No outcome is fabricated: every datum PR-04 owes enters as a
field of `MissingPhaseReceipts`, and this definition constructs no
inhabitant of anything. -/
def partialWitness (m : MissingPhaseReceipts) : OperationalPhaseInstrument where
  core := QuantumSurface.committedAlgebraicPhaseCompletion
  prep := m.prep
  prep_isState := m.prep_isState
  effect := committedEffectPair
  effect_posSemidef := by
    intro c
    rw [Fin.forall_fin_two]
    constructor
    · simpa [committedEffectPair] using
        (committedContextEffect_isEvent c).posSemidef
    · simpa [committedEffectPair] using
        (committedContextEffect_isEvent c).compl.posSemidef
  effect_complete := by
    intro c
    simp only [committedEffectPair, Matrix.cons_val_zero, Matrix.cons_val_one]
    abel
  web_effect := by
    intro c
    simp [committedEffectPair, committedContextEffect]
  phase_effect := by
    simp [committedEffectPair, committedContextEffect,
      QuantumSurface.committedAlgebraicPhaseCompletion]
  phase_sensitive := by
    have h := sourcePhaseLift_distinguishes_Y_states
    simp only [committedEffectPair, Matrix.cons_val_zero,
      committedContextEffect]
    rw [h.1, h.2]
    exact one_ne_zero
  counts := fun c =>
    match c with
    | InstrumentContext.web WebContext.diagonal => (111, 68)
    | InstrumentContext.web (WebContext.conjugated g) => m.rotatedCounts g
    | InstrumentContext.phase => m.phaseCounts
  counts_pos := by
    intro c
    cases c with
    | web wc =>
        cases wc with
        | diagonal =>
            show (0 : ℕ) < 179
            norm_num
        | conjugated g => exact m.rotatedCounts_pos g
    | phase => exact m.phaseCounts_pos
  diagonal_counts_run := rfl
  born_matches := by
    intro c
    cases c with
    | web wc =>
        cases wc with
        | diagonal =>
            simpa [committedEffectPair, committedContextEffect,
              webContextProjector] using m.born_diagonal
        | conjugated g =>
            simpa [committedEffectPair, committedContextEffect,
              webContextProjector] using m.born_rotated g
    | phase =>
        simpa [committedEffectPair, committedContextEffect] using m.born_phase

/-- The partial witness carries the committed diagonal counts: its diagonal
receipt is exactly the run's `realizedOutcomeCounts`. -/
theorem partialWitness_diagonal_committed (m : MissingPhaseReceipts) :
    realizedOutcomeCounts WebContext.diagonal =
      some ((partialWitness m).counts
        (InstrumentContext.web WebContext.diagonal)) :=
  (partialWitness m).diagonal_counts_committed

/-! ## (P5) The reduction Prop -/

/-- **The PR-04 reduction Prop.**  PR-04's operational discharge, typed:
some operational phase instrument exists over the committed algebraic core.
Supplying the missing phase-context receipts is equivalent to inhabiting
this Prop (`dischargesPR04_iff_missing_receipts`).  Nothing in the
repository inhabits it today: the committed payload realizes diagonal
counts only, and this file constructs no inhabitant.  An inhabitant
certifies internal coherence of the supplied data — one preparation, exact
POVMs, exact Born match across all contexts; that the supplied counts are
producer receipts of a committed run is a custody condition the register
tracks outside the type system, and the register's discharge discipline
prohibits meeting this Prop with synthetic counts. -/
def DischargesPR04 : Prop :=
  ∃ I : OperationalPhaseInstrument,
    I.core = QuantumSurface.committedAlgebraicPhaseCompletion

/-- The missing receipts suffice: any inhabitant of `MissingPhaseReceipts`
completes the partial witness into a full instrument over the committed
core. -/
theorem missing_receipts_discharge (m : MissingPhaseReceipts) :
    DischargesPR04 :=
  ⟨partialWitness m, rfl⟩

/-- **The proved reduction.**  Discharging the typed PR-04 Prop is exactly
supplying the missing phase-context receipts: `DischargesPR04` holds if and
only if `MissingPhaseReceipts` is inhabited.  Forward: an instrument over
the committed core yields the preparation, the rotated and phase counts,
and every Born receipt by reading its own fields.  Backward: the partial
witness completes any receipt package into an instrument.  PR-04 stays
open; this equivalence locates its remaining content exactly. -/
theorem dischargesPR04_iff_missing_receipts :
    DischargesPR04 ↔ Nonempty MissingPhaseReceipts := by
  constructor
  · rintro ⟨I, hcore⟩
    refine ⟨⟨I.prep, I.prep_isState,
      fun g => I.counts (InstrumentContext.web (WebContext.conjugated g)),
      fun g => I.counts_pos _,
      I.counts InstrumentContext.phase, I.counts_pos _, ?_, ?_, ?_⟩⟩
    · have h := I.born_matches (InstrumentContext.web WebContext.diagonal)
      rw [I.web_effect WebContext.diagonal, webContextProjector_diagonal,
        I.diagonal_counts_run] at h
      exact h
    · intro g
      have h := I.born_matches (InstrumentContext.web (WebContext.conjugated g))
      rw [I.web_effect (WebContext.conjugated g),
        webContextProjector_conjugated] at h
      exact h
    · have h := I.born_matches InstrumentContext.phase
      rw [I.phase_effect, hcore] at h
      exact h
  · rintro ⟨m⟩
    exact missing_receipts_discharge m

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms webContextProjector_diagonal
#print axioms webContextProjector_conjugated
#print axioms webContextProjector_symm
#print axioms webContextProjector_idem
#print axioms complexifyRealMatrix_isHermitian_of_symmetric
#print axioms webContextProjector_isEvent
#print axioms enumerateContexts_surjective
#print axioms binaryFrequency_diagonal_run
#print axioms OperationalPhaseInstrument.effect_one
#print axioms OperationalPhaseInstrument.effect_isEffect
#print axioms OperationalPhaseInstrument.diagonal_mass
#print axioms OperationalPhaseInstrument.diagonal_counts_committed
#print axioms OperationalPhaseInstrument.born_matches_snd
#print axioms operational_instrument_completes_tomography
#print axioms preparation_identified_by_realized_data
#print axioms sigmaXPlusProj_isEvent
#print axioms sigmaXMinusProj_isEvent
#print axioms rhoXPlus_eq_mix
#print axioms rhoXMinus_eq_mix
#print axioms rhoXPlus_isState
#print axioms rhoXMinus_isState
#print axioms rhoXPlus_ne_rhoXMinus
#print axioms rhoX_agree_on_diagonal
#print axioms rhoY_agree_on_diagonal
#print axioms diagonal_context_insufficient
#print axioms recordProjector_isDiag
#print axioms complexifyRealMatrix_isDiag
#print axioms diagonal_receipts_do_not_identify
#print axioms rhoXPlus_rotated_weight
#print axioms rhoXMinus_rotated_weight
#print axioms rotated_context_separates_x_pair
#print axioms phase_effect_not_diagonal
#print axioms phase_effect_genuinely_complex
#print axioms phase_context_necessary
#print axioms committedContextEffect_isEvent
#print axioms partialWitness
#print axioms partialWitness_diagonal_committed
#print axioms missing_receipts_discharge
#print axioms dischargesPR04_iff_missing_receipts

end EventAlgebra
