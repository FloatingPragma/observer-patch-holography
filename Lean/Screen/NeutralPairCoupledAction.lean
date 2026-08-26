import SeamChargeContinuity
import Geometry.TransportedChargeForceLaw

set_option autoImplicit false

open scoped BigOperators

namespace OPH.NeutralPairCoupledAction

open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.ChargeFixedInteraction
open OPH.CommonWorldJointAction
open OPH.PortChargeMinimalCoupling
open OPH.WorldlineHopTransport
open OPH.TransportedChargeForceLaw
open OPH.SeamChargeContinuity
open OPH.C1Lorentz (lorentzQ)

/-!
# One neutral-pair field--worldline action on the committed carrier

This module closes a narrowly scoped composition gap.  Two declared
seam-step worldlines carry opposite charges `q` and `-q` on the same
twelve-port, thirty-seam carrier.  `neutralPairCoupledAction` is one real
finite-window action: the source-free local Maxwell window action, the two
minimal-coupling terms, and the two committed clock terms.

The field variations of this one action are exactly the existing sourced
interior finite-window Ampere equations and window Gauss equations at the sum
of the two hopping currents and loads.  A closed two-step variation of either
worldline, with the other worldline and the fields fixed, has exactly the
existing clock-plus-force-law action difference.  The action is gauge
invariant under the exact endpoint condition stated below.  Thus the field
and path variation statements no longer live on actions with independently
declared background sources.

The pair is neutral at every step and its sources obey the committed
continuity equation.  The explicit crossings of seams `0` and `29` recover
the nonzero neutral-pair source and field history of `SeamChargeContinuity`,
so the field-stationarity theorem has an inhabited sourced instance.

All of this is conditional finite mathematics.  The charge `q`, step `h`,
clock units, kinetic and clock terms, minimal coupling, carrier map, and
variation class remain declared.  No theorem selects a joint stationary
worldline, identifies the step with laboratory time, or establishes a
continuum Lorentz-force or Maxwell limit.
-/

noncomputable section

/-! ## Opposite-charge sources -/

/-- Total port load of two declared carrier worldlines with charges `q` and
`-q`. -/
def neutralPairLoad (q : ℝ) (wPos wNeg : SeamStepWorldline) : ℕ → Fin 12 → ℝ :=
  hoppingLoad q (hoppingPath wPos) + hoppingLoad (-q) (hoppingPath wNeg)

/-- Total seam current of the same opposite-charge pair. -/
def neutralPairCurrent (q h : ℝ) (wPos wNeg : SeamStepWorldline) : ℕ → Fin 30 → ℝ :=
  hoppingCurrent q h (hoppingPath wPos) + hoppingCurrent (-q) h (hoppingPath wNeg)

/-- The opposite-charge load is neutral at every step, independently of the
two paths. -/
theorem neutralPairLoad_total (q : ℝ) (wPos wNeg : SeamStepWorldline) (n : ℕ) :
    (∑ p : Fin 12, neutralPairLoad q wPos wNeg n p) = 0 := by
  unfold neutralPairLoad
  simp only [Pi.add_apply, Finset.sum_add_distrib, hoppingLoad_total]
  ring

/-- The opposite-charge sources satisfy the committed scaled continuity
equation identically. -/
theorem neutralPair_continuity (q h : ℝ) (hh : h ≠ 0)
    (wPos wNeg : SeamStepWorldline) (n : ℕ) :
    neutralPairLoad q wPos wNeg (n + 1) - neutralPairLoad q wPos wNeg n +
      h • realBoundary (neutralPairCurrent q h wPos wNeg n) = 0 := by
  have hPos := hopping_continuity q h hh (hoppingPath wPos) n
  have hNeg := hopping_continuity (-q) h hh (hoppingPath wNeg) n
  unfold neutralPairLoad neutralPairCurrent
  simp only [Pi.add_apply, map_add, smul_add]
  calc
    (hoppingLoad q (hoppingPath wPos) (n + 1) +
          hoppingLoad (-q) (hoppingPath wNeg) (n + 1)) -
        (hoppingLoad q (hoppingPath wPos) n +
          hoppingLoad (-q) (hoppingPath wNeg) n) +
        (h • realBoundary (hoppingCurrent q h (hoppingPath wPos) n) +
          h • realBoundary (hoppingCurrent (-q) h (hoppingPath wNeg) n)) =
      (hoppingLoad q (hoppingPath wPos) (n + 1) -
          hoppingLoad q (hoppingPath wPos) n +
          h • realBoundary (hoppingCurrent q h (hoppingPath wPos) n)) +
        (hoppingLoad (-q) (hoppingPath wNeg) (n + 1) -
          hoppingLoad (-q) (hoppingPath wNeg) n +
          h • realBoundary (hoppingCurrent (-q) h (hoppingPath wNeg) n)) := by
            abel
    _ = 0 := by rw [hPos, hNeg, add_zero]

/-! ## The single action -/

/-- One finite-window action for a pair of opposite hopping charges and one
carrier field.  The field part is source-free before the two explicit
minimal-coupling terms are added. -/
def neutralPairCoupledAction (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : ℝ :=
  windowAction h N A φ 0 0 +
    interactionA q h N A φ (hoppingPath wPos) +
    interactionA (-q) h N A φ (hoppingPath wNeg) +
    clockAction (N + 1) (generatedPath tauPos wPos) +
    clockAction (N + 1) (generatedPath tauNeg wNeg)

/-- Source pairing is additive in a pair of port loads and seam currents. -/
theorem sourcePairing_add (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ₁ ρ₂ : ℕ → Fin 12 → ℝ)
    (J₁ J₂ : ℕ → Fin 30 → ℝ) :
    sourcePairing h N A φ (ρ₁ + ρ₂) (J₁ + J₂) =
      sourcePairing h N A φ ρ₁ J₁ + sourcePairing h N A φ ρ₂ J₂ := by
  unfold sourcePairing
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun n _ ↦ ?_
  simp only [Pi.add_apply, seamInner_add_left, portInner_add_left]
  ring

/-- The explicit pair action is exactly the committed window action at the
total pair sources, plus the two clock terms.  There is no background source
and no endpoint remainder. -/
theorem neutralPairCoupledAction_eq_window (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) :
    neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      windowAction h N A φ (neutralPairLoad q wPos wNeg)
          (neutralPairCurrent q h wPos wNeg) +
        clockAction (N + 1) (generatedPath tauPos wPos) +
        clockAction (N + 1) (generatedPath tauNeg wNeg) := by
  unfold neutralPairCoupledAction neutralPairLoad neutralPairCurrent interactionA
  have hBase :
      windowAction h N A φ (hoppingLoad q (hoppingPath wPos))
          (hoppingCurrent q h (hoppingPath wPos)) =
        windowAction h N A φ 0 0 +
          sourcePairing h N A φ (hoppingLoad q (hoppingPath wPos))
            (hoppingCurrent q h (hoppingPath wPos)) := by
    simpa only [zero_add] using
      (windowAction_augment_sourcePairing h N A φ 0
        (hoppingLoad q (hoppingPath wPos)) 0
        (hoppingCurrent q h (hoppingPath wPos)))
  rw [windowAction_augment_sourcePairing h N A φ
      (hoppingLoad q (hoppingPath wPos)) (hoppingLoad (-q) (hoppingPath wNeg))
      (hoppingCurrent q h (hoppingPath wPos))
      (hoppingCurrent (-q) h (hoppingPath wNeg)), hBase]

/-- The joined action is invariant under the committed time-dependent gauge
transformation for gauge functions vanishing at the two window endpoints.
This is endpoint gauge invariance: the endpoint condition is the exact
finite-window boundary condition, not a claim of unrestricted gauge
invariance or a claim about asymptotic spacetime. -/
theorem neutralPairCoupledAction_gauge_invariant (q h tauPos tauNeg : ℝ)
    (hh : h ≠ 0) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ χ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) (h0 : χ 0 = 0) (hN : χ (N + 1) = 0) :
    neutralPairCoupledAction q h tauPos tauNeg N (gaugeTransformA A χ)
        (gaugeTransformPhiScaled h φ χ) wPos wNeg =
      neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg := by
  have hWindow :
      windowAction h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) 0 0 =
        windowAction h N A φ 0 0 := by
    unfold windowAction
    refine Finset.sum_congr rfl fun n _ ↦ ?_
    unfold stepLagrangian
    rw [electricFieldScaled_gauge_invariant]
    have hLocal :=
      (localSourcedAction_gauge_invariant_iff (0 : Fin 30 → ℝ)).mpr
        (map_zero realBoundary)
        (A (n + 1)) (χ (n + 1))
    simp only [Pi.zero_apply, portInner_zero_left, mul_zero, add_zero]
    unfold gaugeTransformA
    rw [hLocal]
  have hPos := interactionA_gauge_invariant q h hh N A φ χ (hoppingPath wPos) h0 hN
  have hNeg := interactionA_gauge_invariant (-q) h hh N A φ χ (hoppingPath wNeg) h0 hN
  unfold neutralPairCoupledAction
  rw [hWindow, hPos, hNeg]

/-! ## Field variations -/

/-- **Field Euler--Lagrange equations of the neutral-pair action.**
Stationarity under endpoint-fixed seam-potential variations and arbitrary
port-potential variations is equivalent to the existing sourced scaled
Ampere equations at the interior indices `m < N` and Gauss equations at the
window indices `n < N + 1`, with the *same two worldlines* as sources. -/
theorem neutralPair_field_equations (q h tauPos tauNeg : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) :
    ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        neutralPairCoupledAction q h tauPos tauNeg N (A + a) φ wPos wNeg =
          neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
            quadraticRemainder h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        neutralPairCoupledAction q h tauPos tauNeg N A (φ + f) wPos wNeg =
          neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
            quadraticRemainder h N 0 f)) ↔
    ((∀ m, m < N →
        electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
          h • (localMaxwellOperator (A (m + 1)) -
            neutralPairCurrent q h wPos wNeg m)) ∧
      (∀ n, n < N + 1 →
        realBoundary (electricFieldScaled h A φ n) =
          neutralPairLoad q wPos wNeg n)) := by
  let rho := neutralPairLoad q wPos wNeg
  let current := neutralPairCurrent q h wPos wNeg
  have hActionA :
      (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        neutralPairCoupledAction q h tauPos tauNeg N (A + a) φ wPos wNeg =
          neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
            quadraticRemainder h N a 0) ↔
      (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        windowAction h N (A + a) φ rho current =
          windowAction h N A φ rho current + quadraticRemainder h N a 0) := by
    constructor
    · intro hs a ha0 haN
      have hs' := hs a ha0 haN
      rw [neutralPairCoupledAction_eq_window, neutralPairCoupledAction_eq_window] at hs'
      linarith
    · intro hs a ha0 haN
      have hs' := hs a ha0 haN
      rw [neutralPairCoupledAction_eq_window, neutralPairCoupledAction_eq_window]
      linarith
  have hActionPhi :
      (∀ f : ℕ → Fin 12 → ℝ,
        neutralPairCoupledAction q h tauPos tauNeg N A (φ + f) wPos wNeg =
          neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
            quadraticRemainder h N 0 f) ↔
      (∀ f : ℕ → Fin 12 → ℝ,
        windowAction h N A (φ + f) rho current =
          windowAction h N A φ rho current + quadraticRemainder h N 0 f) := by
    constructor
    · intro hs f
      have hs' := hs f
      rw [neutralPairCoupledAction_eq_window, neutralPairCoupledAction_eq_window] at hs'
      linarith
    · intro hs f
      have hs' := hs f
      rw [neutralPairCoupledAction_eq_window, neutralPairCoupledAction_eq_window]
      linarith
  constructor
  · rintro ⟨hStatA, hStatPhi⟩
    have hResidual :=
      (action_stationary_A_iff_ampere h hh N A φ rho current).mp
        (hActionA.mp hStatA)
    have hGauss :=
      (action_stationary_phi_iff_gauss h hh N A φ rho current).mp
        (hActionPhi.mp hStatPhi)
    constructor
    · intro m hm
      have hr := hResidual m hm
      unfold ampereResidual at hr
      exact sub_eq_zero.mp hr
    · exact hGauss
  · rintro ⟨hAmpere, hGauss⟩
    have hResidual : ∀ m, m < N → ampereResidual h A φ current m = 0 := by
      intro m hm
      unfold ampereResidual
      rw [sub_eq_zero]
      exact hAmpere m hm
    exact ⟨hActionA.mpr
        ((action_stationary_A_iff_ampere h hh N A φ rho current).mpr hResidual),
      hActionPhi.mpr
        ((action_stationary_phi_iff_gauss h hh N A φ rho current).mpr hGauss)⟩

/-! ## Path variations of the same action -/

/-- Varying the positive-charge path while fixing the fields and the
negative-charge path gives exactly the existing discrete force-law
difference. -/
theorem positivePath_closed_variation (q h tauPos tauNeg : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation wPos j) (hj : j + 1 < N + 1) :
    neutralPairCoupledAction q h tauPos tauNeg N A φ v.worldline wNeg -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      clockDifference wPos v + interactionDifference q h A φ wPos v := by
  unfold neutralPairCoupledAction
  have hInteraction := interactionA_closed_variation q h hh N A φ wPos v hj
  have hClock := clockAction_closed_variation tauPos wPos v (N + 1) hj
  linear_combination hInteraction + hClock

/-- Varying the negative-charge path while fixing the fields and the
positive-charge path gives the same force-law formula at charge `-q`. -/
theorem negativePath_closed_variation (q h tauPos tauNeg : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation wNeg j) (hj : j + 1 < N + 1) :
    neutralPairCoupledAction q h tauPos tauNeg N A φ wPos v.worldline -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      clockDifference wNeg v + interactionDifference (-q) h A φ wNeg v := by
  unfold neutralPairCoupledAction
  have hInteraction := interactionA_closed_variation (-q) h hh N A φ wNeg v hj
  have hClock := clockAction_closed_variation tauNeg wNeg v (N + 1) hj
  linear_combination hInteraction + hClock

/-- **Electric force contact for the positive charge.**  On the same pair
action, delaying a forward crossing of the positive worldline by one rest
changes the action by `q h E` on the crossed seam. -/
theorem positivePath_exchange_response (q h tauPos tauNeg : ℝ) (hh : h ≠ 0)
    (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1)
    (e : Fin 30) (hf : wPos.steps j = .forward e)
    (hr : wPos.steps (j + 1) = .rest) :
    neutralPairCoupledAction q h tauPos tauNeg N A φ
          (exchangeVariation wPos j e hf hr).worldline wNeg -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      q * (h * electricFieldScaled h A φ (j + 1) e) := by
  have hSingle :=
    (exchange_action_difference q h tauPos hh N A φ 0 0 wPos j hj e hf hr).1
  unfold neutralPairCoupledAction
  unfold transportedAction monopoleCoupledAction at hSingle
  linear_combination hSingle

/-- **Electric force contact for the negative charge.**  The same forward
exchange on the negative worldline changes the identical pair action by
`-q h E`. -/
theorem negativePath_exchange_response (q h tauPos tauNeg : ℝ) (hh : h ≠ 0)
    (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1)
    (e : Fin 30) (hf : wNeg.steps j = .forward e)
    (hr : wNeg.steps (j + 1) = .rest) :
    neutralPairCoupledAction q h tauPos tauNeg N A φ wPos
          (exchangeVariation wNeg j e hf hr).worldline -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      -(q * (h * electricFieldScaled h A φ (j + 1) e)) := by
  have hSingle :=
    (exchange_action_difference (-q) h tauNeg hh N A φ 0 0 wNeg j hj e hf hr).1
  unfold neutralPairCoupledAction
  unfold transportedAction monopoleCoupledAction at hSingle
  linear_combination hSingle

/-! ## One typed receipt -/

/-- Field stationarity predicate for the one neutral-pair action. -/
def NeutralPairFieldStationary (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : Prop :=
  (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
    neutralPairCoupledAction q h tauPos tauNeg N (A + a) φ wPos wNeg =
      neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
        quadraticRemainder h N a 0) ∧
  (∀ f : ℕ → Fin 12 → ℝ,
    neutralPairCoupledAction q h tauPos tauNeg N A (φ + f) wPos wNeg =
      neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg +
        quadraticRemainder h N 0 f)

/-- The sourced scaled Ampere/Gauss packet at the pair's own sources. -/
def NeutralPairFieldEquations (q h : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : Prop :=
  (∀ m, m < N →
    electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
      h • (localMaxwellOperator (A (m + 1)) - neutralPairCurrent q h wPos wNeg m)) ∧
  (∀ n, n < N + 1 →
    realBoundary (electricFieldScaled h A φ n) = neutralPairLoad q wPos wNeg n)

/-- Exact closed two-step law for the positive-charge path of the same
action. -/
def PositivePathVariationLaw (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : Prop :=
  ∀ (j : ℕ) (v : ClosedTwoStepVariation wPos j), j + 1 < N + 1 →
    neutralPairCoupledAction q h tauPos tauNeg N A φ v.worldline wNeg -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      clockDifference wPos v + interactionDifference q h A φ wPos v

/-- Exact closed two-step law for the negative-charge path of the same
action. -/
def NegativePathVariationLaw (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : Prop :=
  ∀ (j : ℕ) (v : ClosedTwoStepVariation wNeg j), j + 1 < N + 1 →
    neutralPairCoupledAction q h tauPos tauNeg N A φ wPos v.worldline -
        neutralPairCoupledAction q h tauPos tauNeg N A φ wPos wNeg =
      clockDifference wNeg v + interactionDifference (-q) h A φ wNeg v

/-- One receipt tying neutrality, continuity, both field Euler--Lagrange
equations, and both path force-law differences to the identical action. -/
structure NeutralPairActionReceipt (q h tauPos tauNeg : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) : Prop where
  neutral : ∀ n, (∑ p : Fin 12, neutralPairLoad q wPos wNeg n p) = 0
  continuity : ∀ n, neutralPairLoad q wPos wNeg (n + 1) -
    neutralPairLoad q wPos wNeg n +
      h • realBoundary (neutralPairCurrent q h wPos wNeg n) = 0
  fieldEL : NeutralPairFieldStationary q h tauPos tauNeg N A φ wPos wNeg ↔
    NeutralPairFieldEquations q h N A φ wPos wNeg
  positivePath : PositivePathVariationLaw q h tauPos tauNeg N A φ wPos wNeg
  negativePath : NegativePathVariationLaw q h tauPos tauNeg N A φ wPos wNeg

/-- The general neutral pair carries the complete joined receipt whenever
`h ≠ 0`. -/
theorem neutralPairCoupledAction_receipt (q h tauPos tauNeg : ℝ) (hh : h ≠ 0)
    (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wPos wNeg : SeamStepWorldline) :
    NeutralPairActionReceipt q h tauPos tauNeg N A φ wPos wNeg where
  neutral := neutralPairLoad_total q wPos wNeg
  continuity := neutralPair_continuity q h hh wPos wNeg
  fieldEL := by
    unfold NeutralPairFieldStationary NeutralPairFieldEquations
    exact neutralPair_field_equations q h tauPos tauNeg hh N A φ wPos wNeg
  positivePath := by
    intro j v hj
    exact positivePath_closed_variation q h tauPos tauNeg hh N A φ wPos wNeg v hj
  negativePath := by
    intro j v hj
    exact negativePath_closed_variation q h tauPos tauNeg hh N A φ wPos wNeg v hj

/-! ## The committed sourced inhabitant -/

/-- Hopping paths are determined by their port sequence; the `hop` field is
proof data. -/
theorem hoppingPath_eq_of_port_eq (γ δ : HoppingPath) (h : γ.port = δ.port) :
    γ = δ := by
  cases γ with
  | mk p hp =>
    cases δ with
    | mk r hr =>
      dsimp at h
      subst r
      rfl

/-- The transported path of the one-step crossing worldline is exactly the
committed crossing hopping path. -/
theorem hoppingPath_crossingWorldline (e : Fin 30) :
    hoppingPath (crossingWorldline e) = crossingPath e := by
  apply hoppingPath_eq_of_port_eq
  funext n
  exact crossingWorldline_port e n

/-- The general pair load specializes to the existing explicit pair load. -/
theorem committedNeutralPairLoad_eq (n : ℕ) :
    neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) n = pairLoad n := by
  unfold neutralPairLoad pairLoad familyLoad
  rw [Fin.sum_univ_two, hoppingPath_crossingWorldline, hoppingPath_crossingWorldline]
  rfl

/-- The general pair current specializes to the existing explicit pair
current. -/
theorem committedNeutralPairCurrent_eq (h : ℝ) (n : ℕ) :
    neutralPairCurrent 1 h (crossingWorldline 0) (crossingWorldline 29) n =
      pairCurrent h n := by
  unfold neutralPairCurrent pairCurrent familyCurrent
  rw [Fin.sum_univ_two, hoppingPath_crossingWorldline, hoppingPath_crossingWorldline]
  rfl

/-- **Non-vacuous joined field stationarity.**  For every finite window,
the nonzero neutral-pair history already constructed at `h = 1/2` is a
stationary point of this same neutral-pair action in both field slots.
The clock units are arbitrary because they are constant under field
variation. -/
theorem committedNeutralPair_field_stationary (tauPos tauNeg : ℝ) (N : ℕ) :
    (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N
          (neutralPairBundle.A + a) neutralPairBundle.phi
          (crossingWorldline 0) (crossingWorldline 29) =
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N
          neutralPairBundle.A neutralPairBundle.phi
          (crossingWorldline 0) (crossingWorldline 29) +
        quadraticRemainder (1 / 2) N a 0) ∧
    (∀ f : ℕ → Fin 12 → ℝ,
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N
          neutralPairBundle.A (neutralPairBundle.phi + f)
          (crossingWorldline 0) (crossingWorldline 29) =
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N
          neutralPairBundle.A neutralPairBundle.phi
          (crossingWorldline 0) (crossingWorldline 29) +
        quadraticRemainder (1 / 2) N 0 f) := by
  rw [neutralPair_field_equations 1 (1 / 2) tauPos tauNeg (by norm_num) N]
  constructor
  · intro m _
    rw [committedNeutralPairCurrent_eq]
    exact neutralPairBundle.ampere m
  · intro n _
    rw [committedNeutralPairLoad_eq]
    exact neutralPairBundle_gauss n

/-- The joined field-stationary inhabitant is genuinely sourced at the first
seam and carries the explicit Coulomb-started electric field. -/
theorem committedNeutralPair_nonvacuous :
    neutralPairCurrent 1 (1 / 2) (crossingWorldline 0) (crossingWorldline 29) 0 0 = -2 ∧
      neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) 0 0 = 1 ∧
      electricFieldScaled neutralPairBundle.h neutralPairBundle.A
        neutralPairBundle.phi 0 0 = -1 / 6 := by
  rw [committedNeutralPairCurrent_eq, committedNeutralPairLoad_eq]
  exact neutralPairBundle_nonvacuous

/-- At the explicit declared clock unit `3`, both initial carrier crossings
of the committed pair are timelike in the existing Lorentz-module metric.
The unit `3` is independent of the field step `h = 1/2`; identifying the two
would instead make a crossing spacelike.  This is a carrier statement, and
`3` has no laboratory calibration here. -/
theorem committedNeutralPair_initial_crossings_timelike :
    0 < lorentzQ (generatedPath 3 (crossingWorldline 0) (0 + 1) -
      generatedPath 3 (crossingWorldline 0) 0) ∧
    0 < lorentzQ (generatedPath 3 (crossingWorldline 29) (0 + 1) -
      generatedPath 3 (crossingWorldline 29) 0) := by
  constructor
  · exact seam_step_timelike_of_two_lt 3 (by norm_num) (crossingWorldline 0) 0
      (by simp [crossingWorldline])
  · exact seam_step_timelike_of_two_lt 3 (by norm_num) (crossingWorldline 29) 0
      (by simp [crossingWorldline])

/-- At the field step h = 1/2 itself, both initial crossings are spacelike.
This makes the separation between the independently declared carrier clock
unit τ = 3 and the field step machine-checked in the same module. -/
theorem committedNeutralPair_initial_crossings_spacelike_at_field_step :
    lorentzQ (generatedPath (1 / 2) (crossingWorldline 0) (0 + 1) -
      generatedPath (1 / 2) (crossingWorldline 0) 0) < 0 ∧
    lorentzQ (generatedPath (1 / 2) (crossingWorldline 29) (0 + 1) -
      generatedPath (1 / 2) (crossingWorldline 29) 0) < 0 := by
  constructor
  · rw [lorentzQ_generated_step]
    norm_num [crossingWorldline, stepNormSq, seamNormSq]
  · rw [lorentzQ_generated_step]
    norm_num [crossingWorldline, stepNormSq, seamNormSq]

/-! ## The committed paths are not jointly stationary -/

/-- The committed Coulomb-started history has electric field `5/6` on seam
`0` at step `1`.  The initial potential is a pure coboundary, so its local
Maxwell operator vanishes; the value follows from the sourced Ampere step,
the committed initial field `-1/6`, and current `-2`. -/
theorem committedNeutralPair_electric_one_seam_zero :
    electricFieldScaled neutralPairBundle.h neutralPairBundle.A
      neutralPairBundle.phi 1 0 = 5 / 6 := by
  have hL : localMaxwellOperator (neutralPairBundle.A 1) = 0 := by
    change localMaxwellOperator (-((1 / 2 : ℝ) • realCoulombField (pairLoad 0))) = 0
    unfold realCoulombField localMaxwellOperator
    simp only [LinearMap.coe_comp, Function.comp_apply, map_neg, map_smul,
      faceCurvature_coboundary, smul_zero, neg_zero, map_zero]
  have hE0 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi 0 0 = -1 / 6 :=
    neutralPairBundle_nonvacuous.2.2
  have hAmp := congrFun (neutralPairBundle.ampere 0) 0
  change electricFieldScaled (1 / 2) neutralPairBundle.A neutralPairBundle.phi 1 0 -
      electricFieldScaled (1 / 2) neutralPairBundle.A neutralPairBundle.phi 0 0 =
      (1 / 2) * (localMaxwellOperator (neutralPairBundle.A 1) 0 -
        neutralPairBundle.J 0 0) at hAmp
  rw [hL, neutralPairBundle_nonvacuous.1, hE0] at hAmp
  change electricFieldScaled (1 / 2) neutralPairBundle.A
    neutralPairBundle.phi 1 0 = 5 / 6
  norm_num at hAmp ⊢
  linarith

/-- The canonical closed exchange that delays the positive crossing of seam
`0` from step `0` to step `1`, keeping its two-step endpoint fixed. -/
def committedPositiveExchange : ClosedTwoStepVariation (crossingWorldline 0) 0 :=
  exchangeVariation (crossingWorldline 0) 0 0
    (by simp [crossingWorldline]) (by simp [crossingWorldline])

/-- **Exact path nonstationarity of the advertised sourced packet.**  For a
window containing both steps (`1 ≤ N`), keep the committed fields and the
negative path fixed and delay the positive seam-`0` crossing by one rest.  The
same neutral-pair action changes by exactly `5/12`, and hence the advertised
field-stationary packet fails this allowed positive-path exchange. -/
theorem committedNeutralPair_positive_exchange_nonstationary (N : ℕ) (hN : 1 ≤ N) :
    (neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
          neutralPairBundle.phi committedPositiveExchange.worldline
          (crossingWorldline 29) -
        neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12) ∧
      neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
          neutralPairBundle.phi committedPositiveExchange.worldline
          (crossingWorldline 29) ≠
        neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) := by
  have hj : 0 + 1 < N + 1 := by omega
  have hResponse := positivePath_exchange_response 1 (1 / 2) 3 3 (by norm_num) N
    neutralPairBundle.A neutralPairBundle.phi (crossingWorldline 0)
    (crossingWorldline 29) 0 hj 0 (by simp [crossingWorldline])
      (by simp [crossingWorldline])
  have hE1 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi (0 + 1) 0 = 5 / 6 := by
    simpa only [zero_add] using committedNeutralPair_electric_one_seam_zero
  rw [hE1] at hResponse
  norm_num at hResponse
  have hDifference :
      neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
            neutralPairBundle.phi committedPositiveExchange.worldline
            (crossingWorldline 29) -
          neutralPairCoupledAction 1 (1 / 2) 3 3 N neutralPairBundle.A
            neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12 := by
    simpa only [committedPositiveExchange] using hResponse
  refine ⟨hDifference, ?_⟩
  intro hEqual
  rw [hEqual, sub_self] at hDifference
  norm_num at hDifference

/-- Within this fixed quadratic clock-action form, the same path obstruction
is independent of both declared clock units.  For every pair of `tauPos` and
`tauNeg` values, the positive seam-0 delay changes the action by `5/12`.
Varying only these two units therefore cannot make this particular
field-stationary packet jointly stationary.  This does not quantify over
different clock functionals or additional path-dependent terms. -/
theorem committedNeutralPair_positive_exchange_nonstationary_any_clock_unit
    (tauPos tauNeg : ℝ) (N : ℕ) (hN : 1 ≤ N) :
    (neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi committedPositiveExchange.worldline
          (crossingWorldline 29) -
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12) ∧
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi committedPositiveExchange.worldline
          (crossingWorldline 29) ≠
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) := by
  have hj : 0 + 1 < N + 1 := by omega
  have hResponse := positivePath_exchange_response 1 (1 / 2) tauPos tauNeg
    (by norm_num) N neutralPairBundle.A neutralPairBundle.phi
    (crossingWorldline 0) (crossingWorldline 29) 0 hj 0
      (by simp [crossingWorldline]) (by simp [crossingWorldline])
  have hE1 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi (0 + 1) 0 = 5 / 6 := by
    simpa only [zero_add] using committedNeutralPair_electric_one_seam_zero
  rw [hE1] at hResponse
  norm_num at hResponse
  have hDifference :
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi committedPositiveExchange.worldline
            (crossingWorldline 29) -
          neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12 := by
    simpa only [committedPositiveExchange] using hResponse
  refine ⟨hDifference, ?_⟩
  intro hEqual
  rw [hEqual, sub_self] at hDifference
  norm_num at hDifference

/-- At step `0`, the negative carrier contributes current `2` on its crossed
seam `29`, the charge-conjugate counterpart of the committed current `-2` on
seam `0`. -/
theorem committedNeutralPair_current_zero_seam_29 : neutralPairBundle.J 0 29 = 2 := by
  change pairCurrent (1 / 2) 0 29 = 2
  unfold pairCurrent
  rw [familyCurrent_apply, Fin.sum_univ_two]
  have h0 : pairPath 0 = crossingPath 0 := rfl
  have h1 : pairPath 1 = crossingPath 29 := rfl
  have hq0 : pairCharge 0 = 1 := rfl
  have hq1 : pairCharge 1 = -1 := rfl
  rw [h0, h1, hq0, hq1,
    hoppingCurrent_forward 1 (1 / 2) (crossingPath 0) 0 0
      (crossingPath_zero 0) (crossingPath_one 0),
    hoppingCurrent_forward (-1) (1 / 2) (crossingPath 29) 0 29
      (crossingPath_zero 29) (crossingPath_one 29)]
  have h29 : (29 : Fin 30) ≠ 0 := by decide
  norm_num [h29]

/-- The committed initial Coulomb field is `1/6` on the negative carrier's
seam `29`. -/
theorem committedNeutralPair_electric_zero_seam_29 :
    electricFieldScaled neutralPairBundle.h neutralPairBundle.A
      neutralPairBundle.phi 0 29 = 1 / 6 := by
  show electricFieldScaled (1 / 2)
      (coulombSourcedHistory (1 / 2) (pairLoad 0) (pairCurrent (1 / 2)))
      (fun _ ↦ 0) 0 29 = 1 / 6
  rw [coulombSourcedHistory_electric_zero (1 / 2) (by norm_num),
    pairCoulomb_eq_table, seam_29_endpoints.1, seam_29_endpoints.2]
  simp [pairPotentialZ]
  norm_num

/-- The committed sourced Ampere step gives electric field `-5/6` on seam
`29` at step `1`. -/
theorem committedNeutralPair_electric_one_seam_29 :
    electricFieldScaled neutralPairBundle.h neutralPairBundle.A
      neutralPairBundle.phi 1 29 = -5 / 6 := by
  have hL : localMaxwellOperator (neutralPairBundle.A 1) = 0 := by
    change localMaxwellOperator (-((1 / 2 : ℝ) • realCoulombField (pairLoad 0))) = 0
    unfold realCoulombField localMaxwellOperator
    simp only [LinearMap.coe_comp, Function.comp_apply, map_neg, map_smul,
      faceCurvature_coboundary, smul_zero, neg_zero, map_zero]
  have hAmp := congrFun (neutralPairBundle.ampere 0) 29
  have hE0 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi 0 29 = 1 / 6 :=
    committedNeutralPair_electric_zero_seam_29
  change electricFieldScaled (1 / 2) neutralPairBundle.A neutralPairBundle.phi 1 29 -
      electricFieldScaled (1 / 2) neutralPairBundle.A neutralPairBundle.phi 0 29 =
      (1 / 2) * (localMaxwellOperator (neutralPairBundle.A 1) 29 -
        neutralPairBundle.J 0 29) at hAmp
  rw [hL, committedNeutralPair_current_zero_seam_29, hE0] at hAmp
  change electricFieldScaled (1 / 2) neutralPairBundle.A
    neutralPairBundle.phi 1 29 = -5 / 6
  norm_num at hAmp ⊢
  linarith

/-- The canonical closed exchange that delays the negative crossing of seam
`29` from step `0` to step `1`, keeping its two-step endpoint fixed. -/
def committedNegativeExchange : ClosedTwoStepVariation (crossingWorldline 29) 0 :=
  exchangeVariation (crossingWorldline 29) 0 29
    (by simp [crossingWorldline]) (by simp [crossingWorldline])

/-- **Exact charge-conjugate path nonstationarity.**  For `1 ≤ N`, keep the
committed fields and positive path fixed and delay the negative seam-`29`
crossing by one rest.  Its charge and field signs both reverse relative to
the positive carrier, so for every pair of clock parameters in the fixed
clock-action shape the same action again changes by exactly `5/12` and is not
stationary under this allowed discrete exchange. -/
theorem committedNeutralPair_negative_exchange_nonstationary_any_clock_unit
    (tauPos tauNeg : ℝ) (N : ℕ) (hN : 1 ≤ N) :
    (neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0)
          committedNegativeExchange.worldline -
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12) ∧
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0)
          committedNegativeExchange.worldline ≠
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) := by
  have hj : 0 + 1 < N + 1 := by omega
  have hResponse := negativePath_exchange_response 1 (1 / 2) tauPos tauNeg
    (by norm_num) N neutralPairBundle.A neutralPairBundle.phi
    (crossingWorldline 0) (crossingWorldline 29) 0 hj 29
      (by simp [crossingWorldline]) (by simp [crossingWorldline])
  have hE1 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi (0 + 1) 29 = -5 / 6 := by
    simpa only [zero_add] using committedNeutralPair_electric_one_seam_29
  rw [hE1] at hResponse
  norm_num at hResponse
  have hDifference :
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi (crossingWorldline 0)
            committedNegativeExchange.worldline -
          neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 12 := by
    simpa only [committedNegativeExchange] using hResponse
  refine ⟨hDifference, ?_⟩
  intro hEqual
  rw [hEqual, sub_self] at hDifference
  norm_num at hDifference

/-- **The specific simultaneous neutral-pair exchange.**  For `1 ≤ N`,
delay both committed crossings by one rest while holding the committed fields
fixed.  For arbitrary `tauPos` and `tauNeg` in the fixed clock-action shape,
the two exact `5/12` responses add and the same action changes by exactly
`5/6`.  This concerns only these two canonical closed exchanges; it is not a
no-go for other path pairs, fields, or action shapes. -/
theorem committedNeutralPair_simultaneous_exchange_difference_any_clock_unit
    (tauPos tauNeg : ℝ) (N : ℕ) (hN : 1 ≤ N) :
    neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi committedPositiveExchange.worldline
          committedNegativeExchange.worldline -
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
          neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) = 5 / 6 := by
  have hj : 0 + 1 < N + 1 := by omega
  have hPositiveResponse := positivePath_exchange_response 1 (1 / 2) tauPos tauNeg
    (by norm_num) N neutralPairBundle.A neutralPairBundle.phi
    (crossingWorldline 0) committedNegativeExchange.worldline 0 hj 0
      (by simp [crossingWorldline]) (by simp [crossingWorldline])
  have hE1 : electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi (0 + 1) 0 = 5 / 6 := by
    simpa only [zero_add] using committedNeutralPair_electric_one_seam_zero
  rw [hE1] at hPositiveResponse
  norm_num at hPositiveResponse
  have hPositive :
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi committedPositiveExchange.worldline
            committedNegativeExchange.worldline -
          neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg N neutralPairBundle.A
            neutralPairBundle.phi (crossingWorldline 0)
            committedNegativeExchange.worldline = 5 / 12 := by
    simpa only [committedPositiveExchange] using hPositiveResponse
  have hNegative :=
    (committedNeutralPair_negative_exchange_nonstationary_any_clock_unit
      tauPos tauNeg N hN).1
  linear_combination hPositive + hNegative

/-- **Explicit timelike sourced field-stationary packet.**  At `h = 1/2`
and declared carrier clock unit `3`, the same action has the nonzero
Coulomb-started neutral-pair field stationary in both field slots, and both
initial source crossings are timelike in the carrier Lorentz module.  For
every window with `1 ≤ N`, the exchange theorems above prove that this packet
fails the displayed discrete stationarity condition in both its positive- and
negative-path slots. -/
theorem committedNeutralPair_timelike_sourced_field_stationary (N : ℕ) :
    NeutralPairFieldStationary 1 (1 / 2) 3 3 N neutralPairBundle.A
        neutralPairBundle.phi (crossingWorldline 0) (crossingWorldline 29) ∧
      (0 < lorentzQ (generatedPath 3 (crossingWorldline 0) (0 + 1) -
          generatedPath 3 (crossingWorldline 0) 0) ∧
        0 < lorentzQ (generatedPath 3 (crossingWorldline 29) (0 + 1) -
          generatedPath 3 (crossingWorldline 29) 0)) ∧
      (neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
            (crossingWorldline 29) 0 0 = -2 ∧
        neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) 0 0 = 1 ∧
        electricFieldScaled neutralPairBundle.h neutralPairBundle.A
          neutralPairBundle.phi 0 0 = -1 / 6) := by
  refine ⟨?_, committedNeutralPair_initial_crossings_timelike,
    committedNeutralPair_nonvacuous⟩
  unfold NeutralPairFieldStationary
  exact committedNeutralPair_field_stationary 3 3 N

end

end OPH.NeutralPairCoupledAction

#print axioms OPH.NeutralPairCoupledAction.neutralPairLoad_total
#print axioms OPH.NeutralPairCoupledAction.neutralPair_continuity
#print axioms OPH.NeutralPairCoupledAction.neutralPairCoupledAction_eq_window
#print axioms OPH.NeutralPairCoupledAction.neutralPairCoupledAction_gauge_invariant
#print axioms OPH.NeutralPairCoupledAction.neutralPair_field_equations
#print axioms OPH.NeutralPairCoupledAction.positivePath_closed_variation
#print axioms OPH.NeutralPairCoupledAction.negativePath_closed_variation
#print axioms OPH.NeutralPairCoupledAction.positivePath_exchange_response
#print axioms OPH.NeutralPairCoupledAction.negativePath_exchange_response
#print axioms OPH.NeutralPairCoupledAction.neutralPairCoupledAction_receipt
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_field_stationary
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_nonvacuous
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_initial_crossings_spacelike_at_field_step
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_electric_one_seam_zero
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_positive_exchange_nonstationary
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_positive_exchange_nonstationary_any_clock_unit
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_current_zero_seam_29
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_electric_zero_seam_29
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_electric_one_seam_29
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_negative_exchange_nonstationary_any_clock_unit
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_simultaneous_exchange_difference_any_clock_unit
#print axioms OPH.NeutralPairCoupledAction.committedNeutralPair_timelike_sourced_field_stationary
