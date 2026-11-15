import QFT.SourcePhaseLiftBridge
import EventAlgebra.FiniteBuschGleason

/-!
# The conjugation gauge of the missing Born phase

The B13 no-gos prove that the committed real source web cannot see the
Pauli-Y direction and that its generous real closure cannot generate the
phase lift.  This module upgrades that blindness from a family of
countermodels to a structural statement: entrywise complex conjugation is
a statistics-preserving symmetry of the finite quantum formalism which
fixes every conjugation-fixed effect (in particular the whole
complexified real web) pointwise and exchanges the two candidate phase
completions `rhoYPlus` and `rhoYMinus`.

Consequences proved here, all exact and general (not tied to the two
Pauli-Y states):

* conjugation preserves Hermiticity, positive semidefiniteness, states,
  effects, and events, and equals transposition on Hermitian matrices;
* simultaneous conjugation of state and effect conjugates the Born weight;
  therefore a conjugation-fixed effect gives the same real Born weight on a
  state and its conjugate;
* the two phase completions are conjugation-conjugate, and the phase
  lift's conjugate is the opposite lift.

These statements bound the web's ambiguity from below: the conjugation
orbit is invisible to every conjugation-fixed frame.  They do not assert
that the orbit is the whole ambiguity.  The tomography no-go shows the
real web hides the full Pauli-Y coordinate, magnitude included, and the
y-magnitude readout of an unknown state is exactly what the missing
instrument owes.  What is exact at the orbit level: the ambiguity among
the two constructed candidate completions is one conjugation orbit, so
choosing between them is choosing an orientation of a two-element torsor
acted on by conjugation.  A post-hoc rule applied to the committed repair
counts gives one reversal-odd bit in
`Thermodynamics/RepairCurrentOrientation.lean`, and
`QFT/SourceOrientedCompletion.lean` states the declared, applicability-guarded
transport convention and its consequences.  This does not source-select the
rule or its pairing with the torsor.

**Boundary.**  No instrument, outcome, or physical measurement is
constructed here; the y-magnitude readout of an unknown state remains
absent, and B13 stays open on its instrument, additivity, and validation
gates.
-/

namespace OPH.QFT

open Matrix
open scoped ComplexOrder

noncomputable section

/-- Entrywise complex conjugation of a matrix. -/
def matrixConj (M : Matrix (Fin 2) (Fin 2) ℂ) : Matrix (Fin 2) (Fin 2) ℂ :=
  M.map (starRingEnd ℂ)

@[simp] theorem matrixConj_apply (M : Matrix (Fin 2) (Fin 2) ℂ)
    (i j : Fin 2) : matrixConj M i j = star (M i j) := rfl

theorem matrixConj_involutive (M : Matrix (Fin 2) (Fin 2) ℂ) :
    matrixConj (matrixConj M) = M := by
  ext i j
  simp

theorem matrixConj_add (M N : Matrix (Fin 2) (Fin 2) ℂ) :
    matrixConj (M + N) = matrixConj M + matrixConj N := by
  ext i j
  simp

theorem matrixConj_mul (M N : Matrix (Fin 2) (Fin 2) ℂ) :
    matrixConj (M * N) = matrixConj M * matrixConj N := by
  ext i j
  simp [matrixConj, Matrix.mul_apply]

theorem matrixConj_one :
    matrixConj (1 : Matrix (Fin 2) (Fin 2) ℂ) = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [matrixConj]

theorem matrixConj_sub (M N : Matrix (Fin 2) (Fin 2) ℂ) :
    matrixConj (M - N) = matrixConj M - matrixConj N := by
  ext i j
  simp

/-- Conjugation transposes the trace through `star`. -/
theorem matrixConj_trace (M : Matrix (Fin 2) (Fin 2) ℂ) :
    (matrixConj M).trace = star M.trace := by
  simp [Matrix.trace, matrixConj]

/-- On Hermitian matrices, conjugation is transposition. -/
theorem matrixConj_eq_transpose {M : Matrix (Fin 2) (Fin 2) ℂ}
    (hM : M.IsHermitian) : matrixConj M = Mᵀ := by
  ext i j
  have h := congrFun (congrFun hM j) i
  simp only [Matrix.conjTranspose_apply] at h
  simp [matrixConj, Matrix.transpose_apply, ← h]

theorem matrixConj_isHermitian {M : Matrix (Fin 2) (Fin 2) ℂ}
    (hM : M.IsHermitian) : (matrixConj M).IsHermitian := by
  rw [matrixConj_eq_transpose hM]
  exact hM.transpose

theorem matrixConj_posSemidef {M : Matrix (Fin 2) (Fin 2) ℂ}
    (hM : M.PosSemidef) : (matrixConj M).PosSemidef := by
  rw [matrixConj_eq_transpose hM.isHermitian]
  exact hM.transpose

/-- Conjugation maps states to states. -/
theorem matrixConj_isState {ρ : Matrix (Fin 2) (Fin 2) ℂ}
    (hρ : EventAlgebra.IsState ρ) :
    EventAlgebra.IsState (matrixConj ρ) := by
  refine ⟨matrixConj_posSemidef hρ.1, ?_⟩
  rw [matrixConj_trace, hρ.2]
  simp

/-- Conjugation maps effects to effects. -/
theorem matrixConj_isEffect {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : EventAlgebra.IsEffect E) :
    EventAlgebra.IsEffect (matrixConj E) := by
  refine ⟨matrixConj_posSemidef hE.1, ?_⟩
  have h := matrixConj_posSemidef hE.2
  rw [matrixConj_sub, matrixConj_one] at h
  exact h

/-- Conjugation maps projection events to projection events. -/
theorem matrixConj_isEvent {P : Matrix (Fin 2) (Fin 2) ℂ}
    (hP : EventAlgebra.IsEvent P) :
    EventAlgebra.IsEvent (matrixConj P) := by
  refine ⟨matrixConj_isHermitian hP.1, ?_⟩
  rw [← matrixConj_mul, hP.2]

/-! ## Statistics invariance -/

/-- Conjugating both the state and the effect conjugates the Born
weight. -/
theorem bornWeight_matrixConj (ρ E : Matrix (Fin 2) (Fin 2) ℂ) :
    EventAlgebra.bornWeight (matrixConj ρ) (matrixConj E) =
      star (EventAlgebra.bornWeight ρ E) := by
  unfold EventAlgebra.bornWeight
  rw [← matrixConj_mul, matrixConj_trace]

/-- One-sided form: the conjugated state pairs against any effect as the
original state pairs against the conjugated effect. -/
theorem bornWeight_matrixConj_left (ρ E : Matrix (Fin 2) (Fin 2) ℂ) :
    EventAlgebra.bornWeight (matrixConj ρ) E =
      star (EventAlgebra.bornWeight ρ (matrixConj E)) := by
  conv_lhs => rw [← matrixConj_involutive E]
  exact bornWeight_matrixConj ρ (matrixConj E)

/-- Simultaneously conjugating the state and effect preserves the real part
of their Born weight. -/
theorem bornWeight_re_matrixConj (ρ E : Matrix (Fin 2) (Fin 2) ℂ) :
    (EventAlgebra.bornWeight (matrixConj ρ) (matrixConj E)).re =
      (EventAlgebra.bornWeight ρ E).re := by
  rw [bornWeight_matrixConj]
  exact Complex.conj_re _

/-- **Orbit invisibility.**  If every available effect is
conjugation-fixed, then no outcome statistic separates any state from its
conjugate: the conjugation orbit is statistically invisible.  This is a
short corollary of the trace identity; its content is the bridge from the
web's conjugation-fixedness (`matrixConj_complexify`) to orbit blindness.
It bounds the ambiguity from below and makes no claim that conjugation
exhausts it. -/
theorem conj_invisible_on_fixed_effects (ρ E : Matrix (Fin 2) (Fin 2) ℂ)
    (hE : matrixConj E = E) :
    (EventAlgebra.bornWeight (matrixConj ρ) E).re =
      (EventAlgebra.bornWeight ρ E).re := by
  conv_lhs => rw [← hE]
  exact bornWeight_re_matrixConj ρ E

/-! ## The real web is conjugation-fixed -/

/-- Every complexified real matrix is conjugation-fixed. -/
theorem matrixConj_complexify (R : Matrix (Fin 2) (Fin 2) ℝ) :
    matrixConj (complexifyRealMatrix R) = complexifyRealMatrix R := by
  ext i j
  simp [matrixConj, complexifyRealMatrix]

/-! ## The two completions are conjugation-conjugate -/

theorem matrixConj_rhoYPlus : matrixConj rhoYPlus = rhoYMinus := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [matrixConj, rhoYPlus, rhoYMinus, Complex.ext_iff,
      Complex.conj_ofNat, Complex.div_ofNat_re, Complex.div_ofNat_im]

theorem matrixConj_rhoYMinus : matrixConj rhoYMinus = rhoYPlus := by
  rw [← matrixConj_rhoYPlus, matrixConj_involutive]

/-- The phase lift's conjugate is the opposite completion. -/
theorem matrixConj_sourcePhaseLift :
    matrixConj sourcePhaseLift = rhoYMinus := by
  rw [sourcePhaseLift_eq_rhoYPlus, matrixConj_rhoYPlus]

/-- The specific Pauli-Y blindness of the committed web is one instance
of the gauge theorem: the two Y-states are a conjugation orbit. -/
theorem yStates_are_conj_orbit :
    matrixConj rhoYPlus = rhoYMinus ∧ matrixConj rhoYMinus = rhoYPlus :=
  ⟨matrixConj_rhoYPlus, matrixConj_rhoYMinus⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.matrixConj_isState
#print axioms OPH.QFT.matrixConj_isEffect
#print axioms OPH.QFT.matrixConj_isEvent
#print axioms OPH.QFT.bornWeight_matrixConj_left
#print axioms OPH.QFT.bornWeight_re_matrixConj
#print axioms OPH.QFT.conj_invisible_on_fixed_effects
#print axioms OPH.QFT.matrixConj_complexify
#print axioms OPH.QFT.matrixConj_sourcePhaseLift
