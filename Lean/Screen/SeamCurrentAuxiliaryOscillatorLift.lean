import SeamCurrentDirichletGenerator

namespace OPH.SeamCurrentAuxiliaryOscillatorLift

open OPH.SeamCurrentDirichletGenerator

/-!
# Auxiliary algebraic oscillator lift of the seam-current Dirichlet generator

The canonical completion operator `L = I - P` is nonnegative on its Fourier
characters.  This file gives that mathematical operator an explicitly
auxiliary first-order oscillator lift

`q' = p`, `p' = -L q`.

Rather than postulating formal derivatives, the second application of the
first-order generator defines the algebraic `q''` component.  It satisfies
`q'' + L q = 0` exactly.  Plane waves then have the nonnegative auxiliary
frequency branch `omega = sqrt (completionFourierSymbol k)`, so
`omega^2 = completionFourierSymbol k`.

This construction does not provide a trajectory or identify its mnemonic
prime notation with physical time.  It proves no well-posed flow, inner
product, self-adjointness statement, or conserved energy.  It supplies no
physical-position identification, physical clock, photon or other
field-sector attachment, preferred-frame or boost map, energy unit, SI
conversion, continuum limit, global gluing, or laboratory readout.  Those
identifications remain separate premises.
-/

/-- Product-coordinate point of the response-Gram carrier completion. -/
abbrev CarrierPoint := OPH.SeamCurrentDirichletGenerator.Vec3

/-- Complex scalar functions on the response-Gram carrier completion. -/
abbrev CompletionField := CarrierPoint → ℂ

/-- First and second phase coordinates for the auxiliary oscillator lift. -/
abbrev AuxiliaryPhaseState := CompletionField × CompletionField

/-- First-order generator of the explicitly auxiliary oscillator system
`q' = p`, `p' = -L q`.  The prime notation is mnemonic for this algebraic
lift and does not assert a physical clock or an analytic trajectory. -/
noncomputable def auxiliaryPhaseGenerator
    (state : AuxiliaryPhaseState) : AuxiliaryPhaseState :=
  (state.2, fun x ↦ -completionComplexDirichletGenerator state.1 x)

/-- The algebraic second `q` component obtained by applying the first-order
generator twice. -/
noncomputable def auxiliaryQSecond
    (state : AuxiliaryPhaseState) : CompletionField :=
  (auxiliaryPhaseGenerator (auxiliaryPhaseGenerator state)).1

/-- Exact auxiliary second-order equation `q'' + L q = 0`.  Here `q''` is the
defined algebraic component `auxiliaryQSecond`, not an assumed formal
derivative. -/
theorem auxiliary_q_second_order
    (state : AuxiliaryPhaseState) (x : CarrierPoint) :
    auxiliaryQSecond state x +
        completionComplexDirichletGenerator state.1 x = 0 := by
  simp [auxiliaryQSecond, auxiliaryPhaseGenerator]

/-- Sign reversal of the auxiliary second phase coordinate. -/
def auxiliaryTimeReverse
    (state : AuxiliaryPhaseState) : AuxiliaryPhaseState :=
  (state.1, -state.2)

/-- Auxiliary time reversal is an involution. -/
theorem auxiliaryTimeReverse_involutive
    (state : AuxiliaryPhaseState) :
    auxiliaryTimeReverse (auxiliaryTimeReverse state) = state := by
  simp [auxiliaryTimeReverse]

/-- Conjugating the first-order algebraic generator by reversal of the second
phase coordinate changes its sign.  This is a time-reversal-covariance
identity for the auxiliary algebra.  It does not construct an invertible flow
or assert physical time reversal. -/
theorem auxiliaryPhaseGenerator_time_reversal
    (state : AuxiliaryPhaseState) :
    auxiliaryTimeReverse
        (auxiliaryPhaseGenerator (auxiliaryTimeReverse state)) =
      -auxiliaryPhaseGenerator state := by
  ext x <;> simp [auxiliaryTimeReverse, auxiliaryPhaseGenerator]

/-! ## Exact plane-wave branch -/

/-- Nonnegative auxiliary frequency selected by the canonical completion
Fourier symbol. -/
noncomputable def auxiliaryModeFrequency (k : CarrierPoint) : ℝ :=
  Real.sqrt (completionFourierSymbol k)

theorem auxiliaryModeFrequency_nonnegative (k : CarrierPoint) :
    0 ≤ auxiliaryModeFrequency k := by
  exact Real.sqrt_nonneg _

/-- The auxiliary plane-wave dispersion relation is exactly the canonical
completion symbol. -/
theorem auxiliaryModeFrequency_sq (k : CarrierPoint) :
    auxiliaryModeFrequency k ^ 2 = completionFourierSymbol k := by
  exact Real.sq_sqrt (completionFourierSymbol_nonnegative k)

/-- The negative-imaginary first-order exponent associated with the
nonnegative auxiliary frequency branch. -/
noncomputable def auxiliaryModeExponent (k : CarrierPoint) : ℂ :=
  -(auxiliaryModeFrequency k : ℂ) * Complex.I

/-- Phase state of a completion plane wave on the auxiliary positive-frequency
branch. -/
noncomputable def auxiliaryPlaneWaveState
    (k : CarrierPoint) : AuxiliaryPhaseState :=
  (OPH.PrimitivePortTranslationBridge.planeWave k,
    fun x ↦ auxiliaryModeExponent k *
      OPH.PrimitivePortTranslationBridge.planeWave k x)

/-- The plane-wave phase state is an exact eigenstate of the auxiliary
first-order generator. -/
theorem auxiliaryPhaseGenerator_planeWave (k : CarrierPoint) :
    auxiliaryPhaseGenerator (auxiliaryPlaneWaveState k) =
      auxiliaryModeExponent k • auxiliaryPlaneWaveState k := by
  apply Prod.ext
  · funext x
    rfl
  · funext x
    change
      -completionComplexDirichletGenerator
          (OPH.PrimitivePortTranslationBridge.planeWave k) x =
        auxiliaryModeExponent k *
          (auxiliaryModeExponent k *
            OPH.PrimitivePortTranslationBridge.planeWave k x)
    rw [completionComplexDirichletGenerator_planeWave]
    rw [← auxiliaryModeFrequency_sq]
    unfold auxiliaryModeExponent
    push_cast
    have hI : Complex.I ^ 2 = -(1 : ℂ) := by
      rw [pow_two, Complex.I_mul_I]
    ring_nf
    rw [hI]
    ring

/-- On the plane-wave branch the `q` component obeys
`q'' + omega^2 q = 0`, with the same `omega^2` as the canonical completion
Fourier eigenvalue. -/
theorem auxiliaryPlaneWave_q_equation (k x : CarrierPoint) :
    auxiliaryQSecond (auxiliaryPlaneWaveState k) x +
        (auxiliaryModeFrequency k ^ 2 : ℂ) *
          (auxiliaryPlaneWaveState k).1 x = 0 := by
  have hfrequency :
      (auxiliaryModeFrequency k : ℂ) ^ 2 =
        (completionFourierSymbol k : ℂ) := by
    exact_mod_cast auxiliaryModeFrequency_sq k
  rw [hfrequency]
  have hsecond :=
    auxiliary_q_second_order (auxiliaryPlaneWaveState k) x
  have hmode :
      completionComplexDirichletGenerator
          (auxiliaryPlaneWaveState k).1 x =
        (completionFourierSymbol k : ℂ) *
          (auxiliaryPlaneWaveState k).1 x := by
    exact completionComplexDirichletGenerator_planeWave k x
  rw [hmode] at hsecond
  exact hsecond

/-! ## Axiom audit

The lift uses the previously defined equal-source-counting completion
generator and ordinary algebra.  A2 and A3 select that generator only under
the feasibility, naturality, objective, and unique-minimizer premises carried
by `A2A3DirectedSeamProjection`.  The plane-wave results inherit the disclosed
native-decision dependency used for the exact finite seam count.  The file
introduces no new project axiom.  Physical interpretation requires the open
position, trajectory, clock, field-sector, frame, boost, energy-scale, SI,
continuum, gluing, and readout attachments.
-/

#print axioms auxiliary_q_second_order
#print axioms auxiliaryTimeReverse_involutive
#print axioms auxiliaryPhaseGenerator_time_reversal
#print axioms auxiliaryModeFrequency_nonnegative
#print axioms auxiliaryModeFrequency_sq
#print axioms auxiliaryPhaseGenerator_planeWave
#print axioms auxiliaryPlaneWave_q_equation

end OPH.SeamCurrentAuxiliaryOscillatorLift
