import ScaledMaxwellStability

set_option autoImplicit false

namespace OPH.CertifiedScaledStepInstrument

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability

/-!
# Certified scaled-step instrument at the committed carrier threshold

WHAT IS PROVED.  Exact real arithmetic over the committed sharp carrier
constant `Λ = 3 + √5` of `ScaledMaxwellStability`, on the committed
twelve-port, thirty-seam, twenty-oriented-face complex fixed by the
imported modules.

(1) Step certificate.  The declared step `h = 4/5` satisfies the strict
stability condition with a rational margin:
`(4/5)² (3 + √5) ≤ 84/25` and `84/25 < 4` (`stepCertificate_margin`),
hence `(4/5)² (3 + √5) < 4` (`stepCertificate_strict`), via the
square-root bound `√5 < 9/4` from `5 < 81/16`
(`sqrt_five_lt_nine_fourths`).  No numeric axiom is used.

(2) Certified-step receipts.  The stability and energy theorems of
`ScaledMaxwellStability` specialized at the literal step `4/5`:
nonnegativity of the scaled staggered form for every history and every
step count (`certifiedStep_form_nonneg`); exact conservation of the form
for zero current (`certifiedStep_energy_conserved`); the exact per-step
balance with the literal coefficient `2/5`
(`certifiedStep_energy_balance`); the uniform electric and magnetic
energy bounds for every step count (`certifiedStep_energy_bounds`); and
the rational-coefficient corollaries `‖E n‖² ≤ (25/2) 𝓔₀` and
`‖B n‖² ≤ 25 𝓔₀` (`certifiedStep_energy_bounds_rational`), obtained from
the margin `4 - (4/5)² (3 + √5) ≥ 16/25`.

(3) Instrument packaging.  The structure `CertifiedStepInstrument`
carries a positive step, its strict stability certificate against
`Λ = 3 + √5`, and the resulting form-positivity, energy-bound, and
conservation receipts as fields over the committed carrier.  The general
constructor `instrumentOfStrictStep` produces an instrument from any
positive step strictly below the threshold.  `fourFifthsInstrument` is
the explicit inhabitant at the committed step `4/5`
(`fourFifthsInstrument_step`); `halfStepInstrument` is a second
inhabitant at `1/2`, so the committed choice is a declared selection
among qualifying steps, not forced (`step_selection_not_forced`).  The
committed step also inhabits the parent antecedent bundle:
`fourFifthsBundle` is a `ScaledMaxwellBundle` at `h = 4/5`,
`Λ = 3 + √5`, whose history is the committed nonstatic demonstration
history (`fourFifthsBundle_nonstatic`, `fourFifthsBundle_step_matches`),
and the rational energy bound applies to that history
(`fourFifths_receipt_on_demo`).

(4) Design-only handoff interface.  The structure `InstrumentRunHandoff`
lists what a preregistered run of the instrument lane must bind: the
stability-certified instrument, the carrier Courant constant together
with its Courant proof, a declared readout map from field histories to
per-step rational summaries, the declared control histories (port load
and seam current), and a nonempty rational kill band.
`InstrumentRunHandoff.killTriggered` states the kill condition as a
proposition about a history.  `demoHandoff` is an explicit inhabitant at
the committed step with kill band `[-1, 1]` (`demoHandoff_step`,
`demoHandoff_band_consistent`).

WHAT IS NOT PROVED.  The step `4/5` is a declared selection strictly
below the proved threshold `h² (3 + √5) = 4`; no optimality or
uniqueness of the selection is claimed, and `halfStepInstrument`
witnesses that other steps qualify.  No run is authorized or performed
by this module: the handoff interface carries no seed, no execution, and
no frozen prediction, and its readout map, control histories, and kill
band are declared data with no semantics attached.  No physical clock,
units, source coupling, or laboratory readout is claimed.  The
calibration and physical-attachment premises PR-15, PR-53, and PR-54 of
the parent module are open and are not consumed here.

Axiom audit.  Every proof composes the receipts of
`ScaledMaxwellStability` with exact real arithmetic; the module adds no
project axiom and uses no native decision procedure.  The audit lines at
the end of the file show at most `propext`, `Classical.choice`, and
`Quot.sound`.
-/

noncomputable section

/-! ## (1) The step certificate at `h = 4/5` -/

/-- Clean square-root bound feeding the certificate: `√5 < 9/4`,
from `5 < 81/16`. -/
theorem sqrt_five_lt_nine_fourths : Real.sqrt 5 < 9 / 4 := by
  have h : (5 : ℝ) < (9 / 4) ^ 2 := by norm_num
  calc Real.sqrt 5 < Real.sqrt ((9 / 4) ^ 2) :=
        Real.sqrt_lt_sqrt (by norm_num) h
    _ = 9 / 4 := Real.sqrt_sq (by norm_num)

/-- Rational margin of the committed step below the sharp threshold:
`(4/5)² (3 + √5) ≤ 84/25` and `84/25 < 4`. -/
theorem stepCertificate_margin :
    (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5) ≤ 84 / 25 ∧ (84 / 25 : ℝ) < 4 := by
  constructor
  · have h := sqrt_five_lt_nine_fourths
    nlinarith [Real.sqrt_nonneg 5]
  · norm_num

/-- **Strict stability certificate for the committed step.**
`(4/5)² (3 + √5) < 4`. -/
theorem stepCertificate_strict :
    (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5) < 4 :=
  lt_of_le_of_lt stepCertificate_margin.1 stepCertificate_margin.2

/-! ## (2) The certified-step receipts -/

/-- **Certified-step positivity.**  At the step `4/5` the scaled
staggered form is nonnegative for every history and every step count. -/
theorem certifiedStep_form_nonneg (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    0 ≤ fieldEnergyScaled (4 / 5) A φ n :=
  fieldEnergyScaled_nonneg (4 / 5) (by norm_num) (3 + Real.sqrt 5)
    committed_courant_sharp stepCertificate_strict.le A φ n

/-- **Certified-step conservation.**  At the step `4/5` the scaled
staggered form is exactly conserved along every zero-current solution. -/
theorem certifiedStep_energy_conserved (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled (4 / 5) A φ (fun _ ↦ 0)) (n : ℕ) :
    fieldEnergyScaled (4 / 5) A φ n = fieldEnergyScaled (4 / 5) A φ 0 :=
  energy_conserved_scaled (4 / 5) (by norm_num) A φ hAmp n

/-- **Certified-step balance.**  At the step `4/5` the exact per-step
balance of the scaled staggered form carries the literal coefficient
`2/5`. -/
theorem certifiedStep_energy_balance (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled (4 / 5) A φ J) (n : ℕ) :
    fieldEnergyScaled (4 / 5) A φ (n + 1) =
      fieldEnergyScaled (4 / 5) A φ n -
        (2 / 5) * realSeamInner
          (electricFieldScaled (4 / 5) A φ n +
            electricFieldScaled (4 / 5) A φ (n + 1)) (J n) := by
  have h := energy_balance_scaled (4 / 5) (by norm_num) A φ J hAmp n
  rw [show (4 / 5 : ℝ) / 2 = 2 / 5 by norm_num] at h
  exact h

/-- **Certified-step energy bounds.**  At the step `4/5`, for every
zero-current solution, the electric seam energy and the magnetic face
energy are uniformly bounded at every step count by explicit multiples of
the initial staggered form. -/
theorem certifiedStep_energy_bounds (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled (4 / 5) A φ (fun _ ↦ 0)) (n : ℕ) :
    realSeamEnergy (electricFieldScaled (4 / 5) A φ n) ≤
        8 * fieldEnergyScaled (4 / 5) A φ 0 /
          (4 - (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5)) ∧
      faceEnergy (magneticField A n) ≤
        16 * fieldEnergyScaled (4 / 5) A φ 0 /
          (4 - (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5)) :=
  stability_certificate (4 / 5) (by norm_num) (3 + Real.sqrt 5)
    (by positivity) committed_courant_sharp stepCertificate_strict A φ hAmp n

/-- **Rational form of the certified-step energy bounds.**  The margin
`4 - (4/5)² (3 + √5) ≥ 16/25` turns the bounds into rational multiples of
the initial staggered form: `‖E n‖² ≤ (25/2) 𝓔₀` and `‖B n‖² ≤ 25 𝓔₀`. -/
theorem certifiedStep_energy_bounds_rational (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled (4 / 5) A φ (fun _ ↦ 0)) (n : ℕ) :
    realSeamEnergy (electricFieldScaled (4 / 5) A φ n) ≤
        (25 / 2) * fieldEnergyScaled (4 / 5) A φ 0 ∧
      faceEnergy (magneticField A n) ≤
        25 * fieldEnergyScaled (4 / 5) A φ 0 := by
  obtain ⟨hE, hB⟩ := certifiedStep_energy_bounds A φ hAmp n
  have hE0 : 0 ≤ fieldEnergyScaled (4 / 5) A φ 0 :=
    certifiedStep_form_nonneg A φ 0
  have hm := stepCertificate_margin.1
  have hpos : (0 : ℝ) < 4 - (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5) := by
    linarith
  have hprod : 0 ≤ fieldEnergyScaled (4 / 5) A φ 0 *
      (4 - (4 / 5 : ℝ) ^ 2 * (3 + Real.sqrt 5) - 16 / 25) :=
    mul_nonneg hE0 (by linarith)
  constructor
  · refine hE.trans ?_
    rw [div_le_iff₀ hpos]
    nlinarith
  · refine hB.trans ?_
    rw [div_le_iff₀ hpos]
    nlinarith

/-! ## (3) The certified-step instrument -/

/-- A stability-certified scaled-step instrument on the committed
carrier: a positive step, its strict stability certificate against the
sharp carrier constant `Λ = 3 + √5`, and the resulting positivity,
energy-bound, and conservation receipts.  The step value is a declared
selection; the certificate fields are proved propositions over the
committed carrier. -/
structure CertifiedStepInstrument where
  /-- The declared step. -/
  step : ℝ
  /-- The step is positive. -/
  step_pos : 0 < step
  /-- Strict stability certificate against the sharp carrier constant. -/
  stability_strict : step ^ 2 * (3 + Real.sqrt 5) < 4
  /-- Positivity of the scaled staggered form at the declared step. -/
  form_nonneg : ∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ),
    0 ≤ fieldEnergyScaled step A φ n
  /-- Uniform electric and magnetic energy bounds for every step count
  along every zero-current solution at the declared step. -/
  energy_bounded : ∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
    AmpereEvolutionScaled step A φ (fun _ ↦ 0) → ∀ n : ℕ,
      realSeamEnergy (electricFieldScaled step A φ n) ≤
          8 * fieldEnergyScaled step A φ 0 /
            (4 - step ^ 2 * (3 + Real.sqrt 5)) ∧
        faceEnergy (magneticField A n) ≤
          16 * fieldEnergyScaled step A φ 0 /
            (4 - step ^ 2 * (3 + Real.sqrt 5))
  /-- Exact conservation of the scaled staggered form for zero current
  at the declared step. -/
  energy_conserved : ∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
    AmpereEvolutionScaled step A φ (fun _ ↦ 0) → ∀ n : ℕ,
      fieldEnergyScaled step A φ n = fieldEnergyScaled step A φ 0

/-- **General constructor.**  Any positive step strictly below the sharp
threshold yields a certified instrument; every receipt field is
discharged by the corresponding theorem of `ScaledMaxwellStability` at
the sharp carrier constant. -/
def instrumentOfStrictStep (h : ℝ) (h_pos : 0 < h)
    (hstab : h ^ 2 * (3 + Real.sqrt 5) < 4) : CertifiedStepInstrument where
  step := h
  step_pos := h_pos
  stability_strict := hstab
  form_nonneg := fun A φ n ↦
    fieldEnergyScaled_nonneg h h_pos.ne' (3 + Real.sqrt 5)
      committed_courant_sharp hstab.le A φ n
  energy_bounded := fun A φ hAmp n ↦
    stability_certificate h h_pos.ne' (3 + Real.sqrt 5) (by positivity)
      committed_courant_sharp hstab A φ hAmp n
  energy_conserved := fun A φ hAmp n ↦
    energy_conserved_scaled h h_pos.ne' A φ hAmp n

/-- The committed instrument at the declared step `4/5`. -/
def fourFifthsInstrument : CertifiedStepInstrument :=
  instrumentOfStrictStep (4 / 5) (by norm_num) stepCertificate_strict

/-- The committed instrument carries the literal step `4/5`. -/
theorem fourFifthsInstrument_step : fourFifthsInstrument.step = 4 / 5 := rfl

/-- A second inhabitant at the step `1/2`, witnessing that the committed
choice `4/5` is a selection among qualifying steps. -/
def halfStepInstrument : CertifiedStepInstrument :=
  instrumentOfStrictStep (1 / 2) (by norm_num) (by
    have h := sqrt_five_lt_nine_fourths
    nlinarith [Real.sqrt_nonneg 5])

/-- **The step selection is declared, not forced.**  Two certified
instruments with distinct steps inhabit the structure. -/
theorem step_selection_not_forced :
    fourFifthsInstrument.step ≠ halfStepInstrument.step := by
  show (4 / 5 : ℝ) ≠ 1 / 2
  norm_num

/-- The committed step inhabits the parent antecedent bundle: the
`ScaledMaxwellBundle` at `h = 4/5`, `Λ = 3 + √5`, driven by the committed
demonstration history in the gauge `φ = 0` with zero sources. -/
def fourFifthsBundle : ScaledMaxwellBundle where
  h := 4 / 5
  Λ := 3 + Real.sqrt 5
  h_pos := by norm_num
  Λ_nonneg := by positivity
  courant := committed_courant_sharp
  courant_strict := stepCertificate_strict
  A := demoScaledA (4 / 5)
  phi := fun _ ↦ 0
  rho := fun _ ↦ 0
  J := fun _ ↦ 0
  ampere := demoScaled_ampere (4 / 5) (by norm_num)
  gauss_init := by
    rw [electricFieldScaled_temporal_gauge, demoScaledA_one, demoScaledA_zero,
      sub_zero, map_neg, map_smul, demoInitial_boundary, smul_zero, neg_zero]
  continuity := by
    intro n
    simp

/-- The bundle inhabitant and the instrument inhabitant carry the same
committed step. -/
theorem fourFifthsBundle_step_matches :
    fourFifthsBundle.h = fourFifthsInstrument.step := rfl

/-- The bundle history at the committed step is nonstatic: the first two
seam configurations differ. -/
theorem fourFifthsBundle_nonstatic :
    fourFifthsBundle.A 1 ≠ fourFifthsBundle.A 0 := by
  show demoScaledA (4 / 5) 1 ≠ demoScaledA (4 / 5) 0
  rw [demoScaledA_one, demoScaledA_zero]
  have h := demo_nonstatic
  rw [show demoA 1 = demoInitial from by simp only [demoA],
    show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]] at h
  exact h

/-- The rational certified-step bound applied to the nonstatic
demonstration history: the receipts are jointly satisfied on a concrete
evolving solution. -/
theorem fourFifths_receipt_on_demo (n : ℕ) :
    realSeamEnergy
        (electricFieldScaled (4 / 5) (demoScaledA (4 / 5)) (fun _ ↦ 0) n) ≤
      (25 / 2) *
        fieldEnergyScaled (4 / 5) (demoScaledA (4 / 5)) (fun _ ↦ 0) 0 :=
  (certifiedStep_energy_bounds_rational (demoScaledA (4 / 5)) (fun _ ↦ 0)
    (demoScaled_ampere (4 / 5) (by norm_num)) n).1

/-! ## (4) Design-only handoff interface for the instrument lane -/

/-- What a preregistered run of the instrument lane must bind, as typed
data and proved propositions: the stability-certified instrument (which
fixes the committed step), the carrier Courant constant with its Courant
proof (the carrier itself is the committed twelve-port, thirty-seam,
twenty-oriented-face complex fixed by the imported modules), a declared
readout map from field histories to per-step rational summaries, the
declared control histories, and a nonempty rational kill band.  The
structure carries no run, no seed, and no frozen prediction. -/
structure InstrumentRunHandoff where
  /-- The stability-certified instrument whose step the run must use. -/
  instrument : CertifiedStepInstrument
  /-- The carrier Courant constant the certificate is stated against. -/
  carrierConstant : ℝ
  /-- The carrier satisfies the Courant hypothesis at that constant. -/
  carrier_courant : CourantBound carrierConstant
  /-- Declared readout map: field histories to per-step rational
  summaries.  Design-only; no readout semantics are attached. -/
  readout : (ℕ → Fin 30 → ℝ) → (ℕ → Fin 12 → ℝ) → ℕ → ℚ
  /-- Declared control history for the port load. -/
  controlLoad : ℕ → Fin 12 → ℝ
  /-- Declared control history for the seam current. -/
  controlCurrent : ℕ → Fin 30 → ℝ
  /-- Lower edge of the declared kill band. -/
  killLower : ℚ
  /-- Upper edge of the declared kill band. -/
  killUpper : ℚ
  /-- The declared kill band is nonempty. -/
  kill_band_nonempty : killLower ≤ killUpper

/-- The kill condition a preregistered run commits to: the declared
readout of a history leaves the declared band at some step. -/
def InstrumentRunHandoff.killTriggered (H : InstrumentRunHandoff)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) : Prop :=
  ∃ n : ℕ, H.readout A φ n < H.killLower ∨ H.killUpper < H.readout A φ n

/-- Explicit inhabitant of the handoff interface at the committed step:
the `4/5` instrument, the sharp carrier constant `3 + √5`, the zero
readout and zero controls as declared placeholders, and the kill band
`[-1, 1]`. -/
def demoHandoff : InstrumentRunHandoff where
  instrument := fourFifthsInstrument
  carrierConstant := 3 + Real.sqrt 5
  carrier_courant := committed_courant_sharp
  readout := fun _ _ _ ↦ 0
  controlLoad := fun _ ↦ 0
  controlCurrent := fun _ ↦ 0
  killLower := -1
  killUpper := 1
  kill_band_nonempty := by norm_num

/-- The handoff inhabitant binds the committed step `4/5`. -/
theorem demoHandoff_step : demoHandoff.instrument.step = 4 / 5 := rfl

/-- The placeholder readout of the handoff inhabitant never leaves its
kill band, so the interface fields are jointly satisfiable. -/
theorem demoHandoff_band_consistent (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) : ¬ demoHandoff.killTriggered A φ := by
  rintro ⟨n, h | h⟩
  · rw [show demoHandoff.readout A φ n = 0 from rfl,
      show demoHandoff.killLower = -1 from rfl] at h
    norm_num at h
  · rw [show demoHandoff.readout A φ n = 0 from rfl,
      show demoHandoff.killUpper = 1 from rfl] at h
    norm_num at h

end

end OPH.CertifiedScaledStepInstrument

/- Axiom audit: composition of the `ScaledMaxwellStability` receipts with
exact real arithmetic; only the standard axioms appear. -/
#print axioms OPH.CertifiedScaledStepInstrument.sqrt_five_lt_nine_fourths
#print axioms OPH.CertifiedScaledStepInstrument.stepCertificate_margin
#print axioms OPH.CertifiedScaledStepInstrument.stepCertificate_strict
#print axioms OPH.CertifiedScaledStepInstrument.certifiedStep_form_nonneg
#print axioms OPH.CertifiedScaledStepInstrument.certifiedStep_energy_conserved
#print axioms OPH.CertifiedScaledStepInstrument.certifiedStep_energy_balance
#print axioms OPH.CertifiedScaledStepInstrument.certifiedStep_energy_bounds
#print axioms OPH.CertifiedScaledStepInstrument.certifiedStep_energy_bounds_rational
#print axioms OPH.CertifiedScaledStepInstrument.fourFifthsInstrument_step
#print axioms OPH.CertifiedScaledStepInstrument.step_selection_not_forced
#print axioms OPH.CertifiedScaledStepInstrument.fourFifthsBundle_step_matches
#print axioms OPH.CertifiedScaledStepInstrument.fourFifthsBundle_nonstatic
#print axioms OPH.CertifiedScaledStepInstrument.fourFifths_receipt_on_demo
#print axioms OPH.CertifiedScaledStepInstrument.demoHandoff_step
#print axioms OPH.CertifiedScaledStepInstrument.demoHandoff_band_consistent
