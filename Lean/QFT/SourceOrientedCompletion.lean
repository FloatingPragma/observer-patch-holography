import QFT.ConjugationGauge
import RepairCurrentOrientation

/-!
# The source-oriented phase completion and the conditional Born capstone

`QFT/ConjugationGauge.lean` proves that the ambiguity among the two
candidate phase completions of the real source web is exactly a
conjugation orbit: the candidates are conjugation-conjugate, conjugation
preserves states, effects, and every real Born weight, and no
conjugation-fixed effect separates any state from its conjugate.  That
bounds the web's ambiguity from below; the full Pauli-Y coordinate,
magnitude included, is separately invisible to the web, and no theorem
here claims the orbit exhausts the ambiguity.
`Thermodynamics/RepairCurrentOrientation.lean` extracts from the
committed preregistered run an exact orientation bit whose designated
cycle carries a strict forward excess and whose value flips under time
reversal of the counted order.

This module composes the two under one **declared typed convention**:

> when the designated cycle of the committed table carries a strict
> forward excess (kernel-decided for this table), the completion torsor
> is indexed by `Bool` through `phaseCompletion`, and the repair-current
> orientation bit selects the index.

Under that convention the selected completion coincides with the
algebraic phase lift, the time-reversed counted order of the committed
table selects the conjugate completion (a literal receipt of this table,
recorded below), the frame completed by either torsor element is
tomographically complete on states, and the finite Busch–Gleason
representation theorem composes with the oriented frame: every additive
effect valuation is represented by a state that the three oriented
weights pin uniquely.  On a symmetric (detailed-balanced) table the
comparison bit degenerates to the reversal-even constant `false` in both
time directions and the strict-inequality applicability condition fails,
so the convention declares no orientation there; the committed run's
irreversibility is what makes it applicable.

**Boundary.**  The transport of the thermodynamic bit onto the torsor
index is a typed convention, not a theorem; the opposite pairing is
equally admissible, and the coincidence with the phase lift is made by
the declared pairing rather than found.  What is proved: the ambiguity
among the two candidates is one conjugation orbit, the source supplies a
reversal-odd bit satisfying the applicability condition, both torsor
elements complete state tomography (so the convention is outcome-robust
in either branch), and the composition with the representation theorem
is exact.  No y-magnitude instrument, operational additivity derivation,
or preregistered validation of the oriented readout is claimed; those
remain the open B13 gates.
-/

namespace OPH.QFT

open Matrix
open scoped ComplexOrder

noncomputable section

/-! ## The completion torsor and the selection convention -/

/-- The two-element completion torsor: the candidate phase completions
of the real web, indexed by an orientation bit. -/
def phaseCompletion : Bool → Matrix (Fin 2) (Fin 2) ℂ
  | true => rhoYPlus
  | false => rhoYMinus

/-- Conjugation acts on the torsor as the bit flip: the gauge group of
the two-candidate ambiguity is exactly `ℤ₂`. -/
theorem matrixConj_phaseCompletion (b : Bool) :
    matrixConj (phaseCompletion b) = phaseCompletion (!b) := by
  cases b
  · exact matrixConj_rhoYMinus
  · exact matrixConj_rhoYPlus

/-- Both completions are certified states, so the torsor lives inside
the state space and conjugation permutes it there. -/
theorem phaseCompletion_isState (b : Bool) :
    EventAlgebra.IsState (phaseCompletion b) := by
  cases b
  · exact rhoYMinus_isState
  · exact rhoYPlus_isState

/-- The applicability condition of the declared convention: the
designated cycle of the committed table carries a strict forward
excess. -/
def orientationApplicable : Prop :=
  OPH.Thermodynamics.cycleForward OPH.Thermodynamics.repairCounts 5 4 3 <
    OPH.Thermodynamics.cycleForward OPH.Thermodynamics.repairCounts 3 4 5

/-- The committed table satisfies the applicability condition. -/
theorem orientationApplicable_holds : orientationApplicable := by
  unfold orientationApplicable
  decide

/-- **The declared convention**: the source's repair-current orientation
bit selects the completion.  Its use is justified by
`orientationApplicable_holds`; on tables failing the condition no
orientation is declared. -/
def sourceOrientedPhase : Matrix (Fin 2) (Fin 2) ℂ :=
  phaseCompletion OPH.Thermodynamics.repairOrientationBit

/-- The committed run's orientation selects `rhoYPlus`. -/
theorem sourceOrientedPhase_eq_rhoYPlus :
    sourceOrientedPhase = rhoYPlus :=
  (congrArg phaseCompletion
    OPH.Thermodynamics.repairOrientationBit_true).trans rfl

/-- Under the declared pairing, the selected completion coincides with
the algebraic phase lift of the earned pair.  The coincidence is made by
the pairing convention; with the opposite (equally admissible) pairing
the selected completion would be the conjugate lift. -/
theorem sourceOrientedPhase_eq_sourcePhaseLift :
    sourceOrientedPhase = sourcePhaseLift := by
  rw [sourceOrientedPhase_eq_rhoYPlus, sourcePhaseLift_eq_rhoYPlus]

/-! ## Time reversal on the committed table -/

/-- The orientation bit of the time-reversed counted order. -/
def reversedOrientationBit : Bool :=
  decide (OPH.Thermodynamics.cycleForward
      (OPH.Thermodynamics.reversal OPH.Thermodynamics.repairCounts) 5 4 3 <
    OPH.Thermodynamics.cycleForward
      (OPH.Thermodynamics.reversal OPH.Thermodynamics.repairCounts) 3 4 5)

/-- Literal receipt of the committed table: its time-reversed counted
order selects exactly the conjugate completion.  This glues three
decided literals (bit `true`, reversed bit `false`, and the conjugation
of `rhoYPlus`); it is a property of this table, not a general
equivariance, and on a balanced table both bits would coincide. -/
theorem reversal_selects_conjugate :
    phaseCompletion reversedOrientationBit = matrixConj sourceOrientedPhase := by
  have hrev : reversedOrientationBit = false :=
    OPH.Thermodynamics.reversal_flips_orientation
  rw [hrev, sourceOrientedPhase_eq_rhoYPlus, matrixConj_rhoYPlus]
  rfl

/-- On a symmetric (detailed-balanced) table the comparison bit is
`false` for every cycle in either time direction: the bit degenerates to
a reversal-even constant and the strict-inequality applicability
condition of the convention fails identically. -/
theorem symmetric_table_bit_degenerate (C : Fin 8 → Fin 8 → ℕ)
    (hC : ∀ a b : Fin 8, C a b = C b a) (a b c : Fin 8) :
    decide (OPH.Thermodynamics.cycleForward C c b a <
      OPH.Thermodynamics.cycleForward C a b c) = false := by
  have h : OPH.Thermodynamics.cycleForward C c b a =
      OPH.Thermodynamics.cycleForward C a b c := by
    unfold OPH.Thermodynamics.cycleForward
    rw [hC c b, hC b a, hC a c]
    ring
  simp [h]

/-! ## Both torsor elements complete state tomography -/

/-- The tomographic frame completed by a torsor element: the record
projector, one realized rotated web projector, and the chosen
completion. -/
def completionTomography (b : Bool) (rho : Matrix (Fin 2) (Fin 2) ℂ) :
    Fin 3 → ℂ :=
  ![EventAlgebra.bornWeight rho (complexifyRealMatrix recordProjector),
    EventAlgebra.bornWeight rho (complexifyRealMatrix (conjProjector 3)),
    EventAlgebra.bornWeight rho (phaseCompletion b)]

/-- Conjugating the state turns the phase-lifted frame into the starred
minus frame: the two branches transport into each other along the
gauge. -/
theorem sourcePhaseTomography_matrixConj (rho : Matrix (Fin 2) (Fin 2) ℂ) :
    sourcePhaseTomography (matrixConj rho) =
      fun k => star (completionTomography false rho k) := by
  funext k
  fin_cases k <;>
    simp [sourcePhaseTomography, completionTomography, phaseCompletion,
      bornWeight_matrixConj_left, matrixConj_complexify,
      matrixConj_sourcePhaseLift]

/-- **Both branches are tomographically complete.**  The frame completed
by either torsor element identifies every certified state, so the
declared convention is outcome-robust: whichever value the source bit
had taken, the selected frame would separate states. -/
theorem completionTomography_injective_on_states (b : Bool)
    {rho sigma : Matrix (Fin 2) (Fin 2) ℂ}
    (hrho : EventAlgebra.IsState rho) (hsigma : EventAlgebra.IsState sigma)
    (h : completionTomography b rho = completionTomography b sigma) :
    rho = sigma := by
  cases b
  · have hconj : sourcePhaseTomography (matrixConj rho) =
        sourcePhaseTomography (matrixConj sigma) := by
      rw [sourcePhaseTomography_matrixConj, sourcePhaseTomography_matrixConj]
      funext k
      rw [congrFun h k]
    have hcc := sourcePhaseTomography_injective_on_states
      (matrixConj_isState hrho) (matrixConj_isState hsigma) hconj
    have hback := congrArg matrixConj hcc
    rwa [matrixConj_involutive, matrixConj_involutive] at hback
  · have hplus : sourcePhaseTomography rho = sourcePhaseTomography sigma := by
      funext k
      have hk := congrFun h k
      fin_cases k <;>
        simpa [sourcePhaseTomography, completionTomography, phaseCompletion,
          sourcePhaseLift_eq_rhoYPlus] using hk
    exact sourcePhaseTomography_injective_on_states hrho hsigma hplus

/-! ## The oriented frame -/

/-- The oriented tomographic frame: the frame completed by the
source-selected torsor element. -/
def orientedTomography (rho : Matrix (Fin 2) (Fin 2) ℂ) : Fin 3 → ℂ :=
  ![EventAlgebra.bornWeight rho (complexifyRealMatrix recordProjector),
    EventAlgebra.bornWeight rho (complexifyRealMatrix (conjProjector 3)),
    EventAlgebra.bornWeight rho sourceOrientedPhase]

/-- The oriented frame is the bit instance of the two-branch frame. -/
theorem orientedTomography_eq_completionTomography :
    orientedTomography =
      completionTomography OPH.Thermodynamics.repairOrientationBit := rfl

/-- The oriented frame coincides with the phase-lifted frame. -/
theorem orientedTomography_eq_sourcePhaseTomography :
    orientedTomography = sourcePhaseTomography := by
  funext rho
  unfold orientedTomography sourcePhaseTomography
  rw [sourceOrientedPhase_eq_sourcePhaseLift]

/-- **Oriented tomographic completeness.**  Three Born weights against
the oriented frame identify every certified state. -/
theorem orientedTomography_injective_on_states
    {rho sigma : Matrix (Fin 2) (Fin 2) ℂ}
    (hrho : EventAlgebra.IsState rho) (hsigma : EventAlgebra.IsState sigma)
    (h : orientedTomography rho = orientedTomography sigma) :
    rho = sigma := by
  rw [orientedTomography_eq_sourcePhaseTomography] at h
  exact sourcePhaseTomography_injective_on_states hrho hsigma h

/-- The oriented completion separates the conjugation orbit that every
conjugation-fixed frame must conflate: weight `1` on `rhoYPlus`, weight
`0` on its conjugate. -/
theorem sourceOrientedPhase_separates_conj_orbit :
    EventAlgebra.bornWeight rhoYPlus sourceOrientedPhase = 1 ∧
      EventAlgebra.bornWeight (matrixConj rhoYPlus) sourceOrientedPhase = 0 := by
  rw [sourceOrientedPhase_eq_sourcePhaseLift, matrixConj_rhoYPlus]
  exact sourcePhaseLift_distinguishes_Y_states

/-- Without the oriented entry the remaining frame is conjugation-fixed,
so it cannot separate any state from its conjugate: the first two frame
weights agree on every conjugation orbit.  This is the exact sense in
which an orientation datum beyond the conjugation-fixed web is
necessary. -/
theorem unoriented_frame_conj_blind (rho : Matrix (Fin 2) (Fin 2) ℂ) :
    (EventAlgebra.bornWeight (matrixConj rho)
        (complexifyRealMatrix recordProjector)).re =
      (EventAlgebra.bornWeight rho (complexifyRealMatrix recordProjector)).re ∧
    (EventAlgebra.bornWeight (matrixConj rho)
        (complexifyRealMatrix (conjProjector 3))).re =
      (EventAlgebra.bornWeight rho (complexifyRealMatrix (conjProjector 3))).re :=
  ⟨conj_invisible_on_fixed_effects rho _ (matrixConj_complexify _),
    conj_invisible_on_fixed_effects rho _ (matrixConj_complexify _)⟩

/-! ## The conditional Born capstone -/

/-- **Oriented Born capstone.**  Composition of the representation
theorem with oriented tomographic completeness: every additive effect
valuation on the two-by-two event algebra is represented by a state
whose full Born functional it equals, and that representing state is
pinned among all states by its three Born weights against the oriented
frame alone. -/
theorem oriented_born_capstone (v : Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hv : EventAlgebra.IsEffectValuation v) :
    ∃ rho : Matrix (Fin 2) (Fin 2) ℂ,
      (EventAlgebra.IsState rho ∧
        ∀ E, EventAlgebra.IsEffect E →
          v E = (EventAlgebra.bornWeight rho E).re) ∧
      ∀ sigma : Matrix (Fin 2) (Fin 2) ℂ, EventAlgebra.IsState sigma →
        orientedTomography sigma = orientedTomography rho → sigma = rho := by
  obtain ⟨rho, hrho, -⟩ := EventAlgebra.finite_busch_gleason hv
  exact ⟨rho, hrho, fun sigma hsigma h =>
    orientedTomography_injective_on_states hsigma hrho.1 h⟩

/-- Agreement on the three oriented-frame weights forces agreement of
the full Born valuation on every effect: the oriented frame determines
the whole outcome functional of a certified state. -/
theorem oriented_frame_determines_valuation
    {rho sigma : Matrix (Fin 2) (Fin 2) ℂ}
    (hrho : EventAlgebra.IsState rho) (hsigma : EventAlgebra.IsState sigma)
    (h : orientedTomography rho = orientedTomography sigma) :
    ∀ E, EventAlgebra.IsEffect E →
      (EventAlgebra.bornWeight rho E).re =
        (EventAlgebra.bornWeight sigma E).re := by
  intro E _
  rw [orientedTomography_injective_on_states hrho hsigma h]

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.matrixConj_phaseCompletion
#print axioms OPH.QFT.orientationApplicable_holds
#print axioms OPH.QFT.sourceOrientedPhase_eq_sourcePhaseLift
#print axioms OPH.QFT.reversal_selects_conjugate
#print axioms OPH.QFT.symmetric_table_bit_degenerate
#print axioms OPH.QFT.sourcePhaseTomography_matrixConj
#print axioms OPH.QFT.completionTomography_injective_on_states
#print axioms OPH.QFT.orientedTomography_injective_on_states
#print axioms OPH.QFT.sourceOrientedPhase_separates_conj_orbit
#print axioms OPH.QFT.unoriented_frame_conj_blind
#print axioms OPH.QFT.oriented_born_capstone
#print axioms OPH.QFT.oriented_frame_determines_valuation
