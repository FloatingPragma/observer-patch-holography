import Time.ClockReadout
import Geometry.InternalClockRestFrequency
import Geometry.CommonWorldInstrumentJoin
import CertifiedScaledStepInstrument

set_option autoImplicit false

/-!
# PR-15 physical calibration import: SI anchors and the tick dictionary

WHAT IS PROVED.  This module types the PR-15 physical calibration import
one level further than the premise register row: it fixes the exact SI
anchor dictionary, packages the single empirical number of the import as
a declared tick duration, proves the resulting conversion dictionary
affine-exact against the committed affine clock algebra of
`Time/ClockReadout.lean`, and states dimensionful composition receipts
against three committed corpus surfaces.

(1) SI anchor layer.  The 2019 SI fixes the cesium hyperfine frequency,
the speed of light, and the Planck constant by definition.  The structure
`SIAnchors` carries the three values as rational literals together with
their defining equalities.  `siAnchors_unique` proves the structure is a
subsingleton: the anchor dictionary has exactly one inhabitant, so no
measurement uncertainty and no free parameter enters the dictionary
itself.  The empirical content of PR-15 is not a number in this layer; it
is the identification of the corpus clock with the cesium transition,
carried as the propositional premise field of `CalibrationImport`.

(2) Tick dictionary.  `ClockCalibration` carries one declared positive
tick duration `tau` in SI seconds, the single empirical number of the
import.  Derived exact conversions: internal duration `t` corresponds to
`t * tau` laboratory seconds (`labSeconds`); an internal angular rate `m`
per unit internal parameter corresponds to the laboratory frequency
`m / (2 * pi * tau)` hertz (`labFrequency`); an internal energy gap
`E` corresponds to `E * hbarOverTau` joules with the exact-symbolic
factor `(h / (2 * pi)) / tau`.  The dictionary is affine-exact:
`toAffineMap` embeds every calibration into the committed
`PositiveAffineClockMap` algebra, composition of calibrations is a
calibration and commutes with the committed affine composition
(`toAffineMap_comp`), conversion round-trips are identities, and the
second-times-hertz cancellations hold as exact real algebra
(`tau_mul_labFrequency`, `labFrequency_mul_labPeriod`).

(3) Composition receipts, each conditional on a declared tick.
(a) Step dictionary: under a declared tick equal to the record's step
duration, `stepTime` of `Geometry/CommonWorldInstrumentJoin.lean` is the
laboratory time `n * tau` seconds (`stepTime_eq_labSeconds`), the
record's clock period criterion `mass * N * delta` in `2 * pi * ZZ`
becomes a criterion on the laboratory period `N * tau` seconds
(`stepClock_periodic_iff_labSeconds`), and the step clock reads the
laboratory phase of the frequency attached to the record's mass
(`stepClock_lab_frequency`).  (b) Certified scaled step: the committed
step `4/5` of `Screen/CertifiedScaledStepInstrument.lean` is a
dimensionless Courant fraction; under any declared tick the per-step
laboratory duration is `tau` seconds independent of the step value
(`instrumentStepLabDuration_step_independent`,
`certified_step_not_a_duration`), which blocks a wrong seconds-reading
of `4/5`.  (c) Internal-clock rest rate: an internal rate `m` attaches
to the laboratory frequency `m / (2 * pi * tau)` hertz exactly; the
internal clock of `Geometry/InternalClockRestFrequency.lean` read at
laboratory time `t` seconds is the unit-circle rotation at that
frequency (`internalClock_lab_reading`), and its committed internal
period corresponds to the reciprocal laboratory frequency
(`internalClock_lab_period_receipt`).

(4) Non-forcing.  Calibration is an import, never derived: for a nonzero
internal rate, two distinct declared ticks give distinct laboratory
frequencies (`labFrequency_not_forced`), witnessed by the exact
two-instance receipt `tick_declaration_not_forced`.

(5) Worked example.  The tick declared equal to one cesium period,
`tau = 1 / 9192631770` seconds, gives the internal unit angular rate the
laboratory frequency `9192631770 / (2 * pi)` hertz exact-symbolically
(`cesiumTick_unit_rate_frequency`).  This is a worked example of the
dictionary, not a physical identification of the OPH step.

WHAT IS NOT PROVED HERE.  PR-15 keeps its import disposition.  The
identification of the corpus step with any laboratory standard is the
open empirical content of the import; every receipt above is conditional
on a declared tick, and no theorem selects one.  No source-derived clock
exists in the corpus.  No measurement uncertainty is modeled: the SI
anchors are exact by definition and the tick is a declared number, not
an estimated one.  The worked cesium tick is a demonstration of the
dictionary and identifies no physical realization of the internal step.
Nothing in this module calibrates any prediction surface, and nothing
here is a registered, frozen, or scored prediction.
-/

namespace OPH.PhysicalCalibrationImport

noncomputable section

open OPH.C1Lorentz
open OPH.C2Soldering
open OPH.CommonWorldInstrumentJoin
open OPH.CertifiedScaledStepInstrument
open OPH.ScaledMaxwellStability
open OPH.D1Time

/-! ## (1) The SI anchor layer, exact by definition -/

/-- The three 2019 SI anchor values used by the calibration import, as
exact rationals with their defining equalities.  The 2019 SI defines the
cesium-133 ground-state hyperfine transition frequency as exactly
`9192631770` hertz, the speed of light as exactly `299792458` meters per
second, and the Planck constant as exactly `6.62607015e-34` joule
seconds, typed here as `662607015 / 10 ^ 42`.  The values are
definitional, so the structure carries no measurement uncertainty; its
empirical content is zero.  The empirical content of PR-15, the
identification of the corpus clock with the cesium transition, is
carried as the premise field of `CalibrationImport`, not as a number in
this structure. -/
structure SIAnchors where
  /-- The cesium-133 hyperfine transition frequency in hertz. -/
  cesiumHyperfine : ℚ
  /-- The speed of light in meters per second. -/
  lightSpeed : ℚ
  /-- The Planck constant in joule seconds. -/
  planckConstant : ℚ
  /-- Defining equality: the cesium anchor is the exact SI literal. -/
  cesium_exact : cesiumHyperfine = 9192631770
  /-- Defining equality: the light-speed anchor is the exact SI literal. -/
  light_exact : lightSpeed = 299792458
  /-- Defining equality: the Planck anchor is the exact SI literal. -/
  planck_exact : planckConstant = 662607015 / 10 ^ 42

/-- The canonical anchor dictionary: the unique inhabitant of
`SIAnchors`. -/
def siAnchors : SIAnchors where
  cesiumHyperfine := 9192631770
  lightSpeed := 299792458
  planckConstant := 662607015 / 10 ^ 42
  cesium_exact := rfl
  light_exact := rfl
  planck_exact := rfl

/-- **Anchor exactness.**  The anchor dictionary is a subsingleton: any
two inhabitants are equal, because every field is forced to its rational
literal by the defining equalities.  No measurement uncertainty and no
free parameter enters the dictionary itself. -/
theorem siAnchors_unique (a b : SIAnchors) : a = b := by
  obtain ⟨ac, al, ap, hac, hal, hap⟩ := a
  obtain ⟨bc, bl, bp, hbc, hbl, hbp⟩ := b
  subst hac; subst hal; subst hap; subst hbc; subst hbl; subst hbp
  rfl

/-- The anchors carry no freedom: uniqueness of the inhabitant together
with the three exact literal values of the canonical dictionary. -/
theorem anchors_carry_no_freedom :
    (∀ a b : SIAnchors, a = b) ∧
      siAnchors.cesiumHyperfine = 9192631770 ∧
      siAnchors.lightSpeed = 299792458 ∧
      siAnchors.planckConstant = 662607015 / 10 ^ 42 :=
  ⟨siAnchors_unique, rfl, rfl, rfl⟩

/-- The real cast of the Planck anchor is positive. -/
theorem planckConstant_real_pos (a : SIAnchors) :
    (0 : ℝ) < (a.planckConstant : ℝ) := by
  rw [a.planck_exact]
  push_cast
  positivity

/-- The real cast of the cesium anchor is positive. -/
theorem cesiumHyperfine_real_pos (a : SIAnchors) :
    (0 : ℝ) < (a.cesiumHyperfine : ℝ) := by
  rw [a.cesium_exact]
  push_cast
  positivity

/-! ## (2) The calibration dictionary -/

/-- A clock calibration: one declared tick duration `tau > 0` in SI
seconds.  This single number is the empirical import of PR-15: the
identification of one internal step with a laboratory duration.  The
value is declared, never derived; `labFrequency_not_forced` records that
no theorem forces it. -/
structure ClockCalibration where
  /-- The declared tick duration in SI seconds. -/
  tau : ℝ
  /-- Declared positivity of the tick. -/
  tau_pos : 0 < tau

namespace ClockCalibration

variable (cal : ClockCalibration)

/-- Internal duration `t` (in steps) corresponds to `t * tau` laboratory
seconds. -/
def labSeconds (t : ℝ) : ℝ := t * cal.tau

@[simp] theorem labSeconds_def (t : ℝ) : cal.labSeconds t = t * cal.tau := rfl

/-- Laboratory duration `t` seconds corresponds to `t / tau` internal
steps. -/
def internalOfSeconds (t : ℝ) : ℝ := t / cal.tau

@[simp] theorem internalOfSeconds_def (t : ℝ) :
    cal.internalOfSeconds t = t / cal.tau := rfl

/-- An internal angular rate `m` (radians per internal step) corresponds
to the laboratory frequency `m / (2 * pi * tau)` hertz. -/
def labFrequency (m : ℝ) : ℝ := m / (2 * Real.pi * cal.tau)

@[simp] theorem labFrequency_def (m : ℝ) :
    cal.labFrequency m = m / (2 * Real.pi * cal.tau) := rfl

/-- The laboratory period of an internal angular rate `m`: the internal
period `2 * pi / m` converted to seconds. -/
def labPeriod (m : ℝ) : ℝ := cal.labSeconds (2 * Real.pi / m)

/-- Conversion round-trip on durations: seconds of steps of seconds. -/
theorem labSeconds_internalOfSeconds (t : ℝ) :
    cal.labSeconds (cal.internalOfSeconds t) = t := by
  simp only [labSeconds_def, internalOfSeconds_def]
  exact div_mul_cancel₀ t cal.tau_pos.ne'

/-- Conversion round-trip on durations: steps of seconds of steps. -/
theorem internalOfSeconds_labSeconds (t : ℝ) :
    cal.internalOfSeconds (cal.labSeconds t) = t := by
  simp only [labSeconds_def, internalOfSeconds_def]
  exact mul_div_cancel_right₀ t cal.tau_pos.ne'

/-- **Second-times-hertz cancellation.**  One tick duration times the
laboratory frequency of an internal rate `m` is the dimensionless phase
per step over `2 * pi`, as exact real algebra: the seconds unit of `tau`
cancels the hertz unit of the frequency. -/
theorem tau_mul_labFrequency (m : ℝ) :
    cal.tau * cal.labFrequency m = m / (2 * Real.pi) := by
  simp only [labFrequency_def]
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- **Frequency-times-period cancellation.**  For a nonzero internal
rate the laboratory frequency times the laboratory period is exactly
one: hertz times seconds cancels as exact real algebra. -/
theorem labFrequency_mul_labPeriod {m : ℝ} (hm : m ≠ 0) :
    cal.labFrequency m * cal.labPeriod m = 1 := by
  simp only [labFrequency_def, labPeriod, labSeconds_def]
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- The embedding of a calibration into the committed affine clock
algebra of `Time/ClockReadout.lean`: the positive affine map with scale
`tau` and offset zero. -/
def toAffineMap : PositiveAffineClockMap where
  scale := cal.tau
  offset := 0
  scalePositive := cal.tau_pos

/-- The affine embedding applies as the seconds conversion. -/
theorem toAffineMap_map (t : ℝ) :
    cal.toAffineMap.map t = cal.labSeconds t := by
  simp only [toAffineMap, PositiveAffineClockMap.map, labSeconds_def]
  ring

/-- Composition of two calibrations: the tick durations multiply.  The
composite converts through the first dictionary and then rescales
through the second. -/
def comp (second first : ClockCalibration) : ClockCalibration where
  tau := second.tau * first.tau
  tau_pos := mul_pos second.tau_pos first.tau_pos

@[simp] theorem comp_tau (second first : ClockCalibration) :
    (comp second first).tau = second.tau * first.tau := rfl

/-- **Affine exactness of composition.**  The embedding into the
committed affine clock algebra sends composition of calibrations to the
committed affine composition: the dictionary composes with the existing
calibration algebra instead of duplicating it. -/
theorem toAffineMap_comp (second first : ClockCalibration) :
    (comp second first).toAffineMap =
      PositiveAffineClockMap.comp second.toAffineMap first.toAffineMap := by
  simp only [toAffineMap, comp_tau, PositiveAffineClockMap.comp]
  simp

/-- **Affine round-trip.**  Through the committed affine inverse, every
conversion round-trips to the identity. -/
theorem toAffineMap_roundtrip (t : ℝ) :
    cal.toAffineMap.inverse.map (cal.toAffineMap.map t) = t :=
  PositiveAffineClockMap.inverse_map_map cal.toAffineMap t

end ClockCalibration

/-- The exact-symbolic energy conversion factor `(h / (2 * pi)) / tau`
in joules per internal energy unit: the reduced Planck anchor divided by
the declared tick. -/
def hbarOverTau (a : SIAnchors) (cal : ClockCalibration) : ℝ :=
  ((a.planckConstant : ℝ) / (2 * Real.pi)) / cal.tau

/-- The factor as a single quotient. -/
theorem hbarOverTau_eq (a : SIAnchors) (cal : ClockCalibration) :
    hbarOverTau a cal = (a.planckConstant : ℝ) / (2 * Real.pi * cal.tau) :=
  div_div _ _ _

/-- The energy conversion factor is positive. -/
theorem hbarOverTau_pos (a : SIAnchors) (cal : ClockCalibration) :
    0 < hbarOverTau a cal :=
  div_pos (div_pos (planckConstant_real_pos a)
    (mul_pos two_pos Real.pi_pos)) cal.tau_pos

/-- An internal energy gap `E` corresponds to `E * hbarOverTau`
joules. -/
def energyJoules (a : SIAnchors) (cal : ClockCalibration) (E : ℝ) : ℝ :=
  E * hbarOverTau a cal

/-- A laboratory energy `EJ` joules corresponds to `EJ / hbarOverTau`
internal energy units. -/
def internalOfJoules (a : SIAnchors) (cal : ClockCalibration) (EJ : ℝ) : ℝ :=
  EJ / hbarOverTau a cal

/-- Conversion round-trip on energies: internal of joules of
internal. -/
theorem internalOfJoules_energyJoules (a : SIAnchors) (cal : ClockCalibration)
    (E : ℝ) : internalOfJoules a cal (energyJoules a cal E) = E :=
  mul_div_cancel_right₀ E (hbarOverTau_pos a cal).ne'

/-- Conversion round-trip on energies: joules of internal of joules. -/
theorem energyJoules_internalOfJoules (a : SIAnchors) (cal : ClockCalibration)
    (EJ : ℝ) : energyJoules a cal (internalOfJoules a cal EJ) = EJ :=
  div_mul_cancel₀ EJ (hbarOverTau_pos a cal).ne'

/-! ## The import package: anchors, tick, and the premise field -/

/-- The full PR-15 import package: the anchor dictionary, one declared
tick, and the empirical identification premise.  The premise is carried
as a propositional field supplied by the declarer, not as a number: a
physical declaration would supply the laboratory statement identifying
the corpus step with a realized standard, and no theorem in the corpus
supplies or discharges that statement.  PR-15 keeps its import
disposition. -/
structure CalibrationImport where
  /-- The anchor dictionary; unique by `siAnchors_unique`. -/
  anchors : SIAnchors
  /-- The declared tick. -/
  clock : ClockCalibration
  /-- The declarer's empirical identification statement, carried as a
  proposition, not a number. -/
  identification : Prop
  /-- The declaration that the premise is asserted. -/
  identification_declared : identification

/-- The anchor component of every import package is the canonical
dictionary: the import's freedom sits entirely in the tick and the
premise field. -/
theorem calibrationImport_anchors_unique (imp : CalibrationImport) :
    imp.anchors = siAnchors :=
  siAnchors_unique imp.anchors siAnchors

/-! ## (3a) Receipt: the common-world step dictionary -/

/-- Under a tick declared equal to the step duration, the committed step
time `stepTime tau n` of `Geometry/CommonWorldInstrumentJoin.lean` is
the laboratory time `n * tau` seconds. -/
theorem stepTime_eq_labSeconds (cal : ClockCalibration) (n : ℕ) :
    stepTime cal.tau n = cal.labSeconds (n : ℝ) := rfl

/-- **Laboratory reading of the clock period criterion.**  Conditional
on the declared tick `W.stepDuration = cal.tau`: the record's step clock
repeats after `N` steps exactly when the record's mass times the
laboratory period `N * tau` seconds is an integer multiple of
`2 * pi`.  The committed criterion `mass * stepTime delta N` in
`2 * pi * ZZ` becomes a criterion on a laboratory period; the tick
identification itself stays a declared premise. -/
theorem stepClock_periodic_iff_labSeconds
    (W : InstrumentedCommonWorldArchitecture) (cal : ClockCalibration)
    (tick : W.stepDuration = cal.tau) (N : ℕ) :
    (∀ n : ℕ, stepClock W (n + N) = stepClock W n) ↔
      ∃ m : ℤ, W.mass * cal.labSeconds (N : ℝ) = m * (2 * Real.pi) := by
  rw [stepClock_periodic_iff, tick, stepTime_eq_labSeconds]

/-- Under the dictionary `delta = tau` a laboratory angular rate
`omega` (radians per second) advances by `omega * tau` radians per step;
the dictionary attaches that internal per-step rate to the laboratory
frequency `omega / (2 * pi)` hertz, with the tick cancelling exactly. -/
theorem labFrequency_of_angular_rate (cal : ClockCalibration) (omega : ℝ) :
    cal.labFrequency (omega * cal.tau) = omega / (2 * Real.pi) := by
  simp only [ClockCalibration.labFrequency_def]
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- **Laboratory phase of the step clock.**  Conditional on the declared
tick `W.stepDuration = cal.tau`: the record's worldline parameter is in
laboratory seconds, so the record's mass field is a laboratory angular
rate, its per-step angular advance is `mass * tau` radians, and at step
`n` the clock reads the unit-circle rotation of laboratory phase
`2 * pi * nu * t` with `nu = labFrequency (mass * tau) = mass / (2 * pi)`
hertz and `t = n * tau` the laboratory time in seconds. -/
theorem stepClock_lab_frequency (W : InstrumentedCommonWorldArchitecture)
    (cal : ClockCalibration) (tick : W.stepDuration = cal.tau) (n : ℕ) :
    stepClock W n =
      Complex.exp (-Complex.I *
        ((2 * Real.pi * cal.labFrequency (W.mass * cal.tau) *
          cal.labSeconds (n : ℝ) : ℝ) : ℂ)) := by
  rw [stepClock_eq_exp]
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  have h : W.mass * stepTime W.stepDuration n =
      2 * Real.pi * cal.labFrequency (W.mass * cal.tau) *
        cal.labSeconds (n : ℝ) := by
    rw [tick, stepTime_eq_labSeconds, labFrequency_of_angular_rate]
    simp only [ClockCalibration.labSeconds_def]
    field_simp
  rw [h]

/-! ## (3b) Receipt: the certified scaled step is dimensionless -/

/-- The laboratory duration of one evolution step of a certified-step
instrument under a declared tick.  The instrument argument is
deliberately unused: the dictionary assigns the duration from the tick
alone. -/
def instrumentStepLabDuration (cal : ClockCalibration)
    (_I : CertifiedStepInstrument) : ℝ :=
  cal.labSeconds 1

/-- The laboratory duration of one evolution step of a scaled Maxwell
bundle under a declared tick.  The bundle argument is deliberately
unused. -/
def bundleStepLabDuration (cal : ClockCalibration)
    (_S : ScaledMaxwellBundle) : ℝ :=
  cal.labSeconds 1

/-- The per-step laboratory duration is the tick. -/
theorem instrumentStepLabDuration_eq_tau (cal : ClockCalibration)
    (I : CertifiedStepInstrument) :
    instrumentStepLabDuration cal I = cal.tau := one_mul cal.tau

/-- The per-step laboratory duration of a bundle is the tick. -/
theorem bundleStepLabDuration_eq_tau (cal : ClockCalibration)
    (S : ScaledMaxwellBundle) :
    bundleStepLabDuration cal S = cal.tau := one_mul cal.tau

/-- **Step-value independence.**  Under one declared tick, every two
certified-step instruments have the same per-step laboratory duration:
the certified step value is a dimensionless Courant fraction and carries
no seconds. -/
theorem instrumentStepLabDuration_step_independent (cal : ClockCalibration)
    (I I' : CertifiedStepInstrument) :
    instrumentStepLabDuration cal I = instrumentStepLabDuration cal I' := rfl

/-- **The certified step is not a duration.**  The committed step `4/5`
is calibration-invariant, while the per-step laboratory duration changes
with the declared tick: for two ticks with distinct durations the
dimensionless step agrees and the laboratory durations differ.  This
blocks the wrong reading of `4/5` as a number of seconds. -/
theorem certified_step_not_a_duration (cal cal' : ClockCalibration)
    (hne : cal.tau ≠ cal'.tau) :
    fourFifthsInstrument.step = fourFifthsInstrument.step ∧
      instrumentStepLabDuration cal fourFifthsInstrument ≠
        instrumentStepLabDuration cal' fourFifthsInstrument := by
  refine ⟨rfl, ?_⟩
  rw [instrumentStepLabDuration_eq_tau, instrumentStepLabDuration_eq_tau]
  exact hne

/-- The committed bundle at the step `4/5` has the same per-step
laboratory duration as every other bundle under the same tick, equal to
the tick itself. -/
theorem fourFifthsBundle_lab_step_duration (cal : ClockCalibration) :
    bundleStepLabDuration cal fourFifthsBundle = cal.tau :=
  bundleStepLabDuration_eq_tau cal fourFifthsBundle

/-! ## (3c) Receipt: the internal-clock rest rate -/

/-- **Laboratory reading of the internal clock.**  Conditional on the
declared tick: the committed internal clock of
`Geometry/InternalClockRestFrequency.lean` with internal rate `m`, read
at laboratory time `t` seconds (internal parameter `t / tau`), is the
unit-circle rotation of laboratory phase `2 * pi * nu * t` with
`nu = m / (2 * pi * tau)` hertz.  An internal rate attaches to a
laboratory frequency exactly. -/
theorem internalClock_lab_reading (cal : ClockCalibration) (m : ℝ)
    (frame : FrameHyperboloid) (t : ℝ) :
    internalClock m frame (cal.internalOfSeconds t) =
      Complex.exp (-Complex.I *
        ((2 * Real.pi * cal.labFrequency m * t : ℝ) : ℂ)) := by
  rw [internalClock_eq_exp]
  have hτ : cal.tau ≠ 0 := cal.tau_pos.ne'
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  have h : m * cal.internalOfSeconds t =
      2 * Real.pi * cal.labFrequency m * t := by
    simp only [ClockCalibration.internalOfSeconds_def,
      ClockCalibration.labFrequency_def]
    field_simp
  rw [h]

/-- **Laboratory period of the internal clock.**  Conditional on the
declared tick: for a nonzero internal rate the committed internal period
`2 * pi / m` holds, and its laboratory duration in seconds is the exact
reciprocal of the attached laboratory frequency. -/
theorem internalClock_lab_period_receipt (cal : ClockCalibration)
    {m : ℝ} (hm : m ≠ 0) (frame : FrameHyperboloid) :
    Function.Periodic (internalClock m frame) (2 * Real.pi / m) ∧
      cal.labFrequency m * cal.labPeriod m = 1 :=
  ⟨internalClock_periodic hm frame, cal.labFrequency_mul_labPeriod hm⟩

/-! ## (4) Non-forcing: the tick is an import, never derived -/

/-- **Calibration non-forcing, general form.**  For a nonzero internal
rate, two calibrations with distinct declared ticks give distinct
laboratory frequencies: no internal statement fixes the laboratory
frequency of an internal rate without a declared tick. -/
theorem labFrequency_not_forced (cal cal' : ClockCalibration)
    {m : ℝ} (hm : m ≠ 0) (hne : cal.tau ≠ cal'.tau) :
    cal.labFrequency m ≠ cal'.labFrequency m := by
  simp only [ClockCalibration.labFrequency_def]
  have h2π : (2 * Real.pi : ℝ) ≠ 0 :=
    mul_ne_zero two_ne_zero Real.pi_ne_zero
  have hden : 2 * Real.pi * cal.tau ≠ 0 :=
    mul_ne_zero h2π cal.tau_pos.ne'
  have hden' : 2 * Real.pi * cal'.tau ≠ 0 :=
    mul_ne_zero h2π cal'.tau_pos.ne'
  intro h
  have h2 : m * (2 * Real.pi * cal'.tau) = m * (2 * Real.pi * cal.tau) :=
    (div_eq_div_iff hden hden').mp h
  have h3 := mul_left_cancel₀ hm h2
  have h4 := mul_left_cancel₀ h2π h3
  exact hne h4.symm

/-- The declared unit tick: one internal step per SI second.  A worked
declaration, not a physical identification. -/
def unitTick : ClockCalibration where
  tau := 1
  tau_pos := one_pos

@[simp] theorem unitTick_tau : unitTick.tau = 1 := rfl

/-- The declared double tick: one internal step per two SI seconds.  A
worked declaration, not a physical identification. -/
def doubleTick : ClockCalibration where
  tau := 2
  tau_pos := two_pos

@[simp] theorem doubleTick_tau : doubleTick.tau = 2 := rfl

/-- **Two-instance non-forcing receipt.**  The unit tick and the double
tick are distinct declarations, and they give the internal unit rate
distinct laboratory frequencies.  Calibration is an import: the tick is
declared, never derived. -/
theorem tick_declaration_not_forced :
    unitTick.tau ≠ doubleTick.tau ∧
      unitTick.labFrequency 1 ≠ doubleTick.labFrequency 1 := by
  have hne : unitTick.tau ≠ doubleTick.tau := by norm_num
  exact ⟨hne, labFrequency_not_forced unitTick doubleTick one_ne_zero hne⟩

/-! ## (5) Worked exact example: the cesium-period tick -/

/-- The tick declared equal to one cesium period,
`tau = 1 / 9192631770` seconds.  This is a worked example of the
dictionary, not a physical identification of the OPH step: no theorem
identifies the corpus step with the cesium transition, and PR-15 keeps
its import disposition. -/
def cesiumTick : ClockCalibration where
  tau := 1 / 9192631770
  tau_pos := by norm_num

@[simp] theorem cesiumTick_tau : cesiumTick.tau = 1 / 9192631770 := rfl

/-- **Worked frequency receipt.**  Under the cesium-period tick, the
internal unit angular rate has the exact-symbolic laboratory frequency
`9192631770 / (2 * pi)` hertz.  The `2 * pi` divisor records that the
internal rate is angular: an internal rate of `2 * pi` radians per step
would attach to the anchor frequency itself. -/
theorem cesiumTick_unit_rate_frequency :
    cesiumTick.labFrequency 1 = 9192631770 / (2 * Real.pi) := by
  simp only [ClockCalibration.labFrequency_def, cesiumTick_tau,
    mul_one_div, one_div_div]

/-- The worked frequency written through the anchor dictionary: the
laboratory frequency of the internal unit rate under the cesium-period
tick is the cesium anchor over `2 * pi`. -/
theorem cesiumTick_unit_rate_frequency_anchor (a : SIAnchors) :
    cesiumTick.labFrequency 1 = (a.cesiumHyperfine : ℝ) / (2 * Real.pi) := by
  rw [cesiumTick_unit_rate_frequency, a.cesium_exact]
  norm_num

/-- The cesium-period tick and the cesium anchor multiply to one
exactly: one anchor period per tick. -/
theorem cesiumTick_matches_anchor (a : SIAnchors) :
    (a.cesiumHyperfine : ℝ) * cesiumTick.tau = 1 := by
  rw [a.cesium_exact, cesiumTick_tau]
  push_cast
  norm_num

/-- **Worked energy receipt.**  Under the cesium-period tick, the
internal unit energy gap corresponds to the exact-symbolic laboratory
energy `(662607015 / 10 ^ 42) * 9192631770 / (2 * pi)` joules.  The
value is kept symbolic; no numerical evaluation enters. -/
theorem cesiumTick_energy_quantum (a : SIAnchors) :
    hbarOverTau a cesiumTick =
      (662607015 / 10 ^ 42 : ℝ) * 9192631770 / (2 * Real.pi) := by
  rw [hbarOverTau_eq, a.planck_exact, cesiumTick_tau]
  have hπ : Real.pi ≠ 0 := Real.pi_ne_zero
  push_cast
  field_simp

/-- The worked import package: canonical anchors, the cesium-period
tick, and an internal restatement standing in the premise field.  The
proposition here is an internal equality about the declared tick, marked
as a stand-in: a physical declaration would carry the laboratory
identification of the corpus step, which no theorem supplies. -/
def workedImport : CalibrationImport where
  anchors := siAnchors
  clock := cesiumTick
  identification := cesiumTick.tau = 1 / 9192631770
  identification_declared := rfl

/-- A second import package at the unit tick, witnessing at the package
level that the tick declaration is free. -/
def workedImportUnit : CalibrationImport where
  anchors := siAnchors
  clock := unitTick
  identification := unitTick.tau = 1
  identification_declared := rfl

/-- The two worked packages share the forced anchors and differ only in
the declared tick: the import's freedom is exactly the tick and the
premise field. -/
theorem workedImports_differ_only_in_tick :
    workedImport.anchors = workedImportUnit.anchors ∧
      workedImport.clock.tau ≠ workedImportUnit.clock.tau := by
  refine ⟨rfl, ?_⟩
  show (1 / 9192631770 : ℝ) ≠ 1
  norm_num

/-! ## Axiom audit -/

#print axioms siAnchors_unique
#print axioms anchors_carry_no_freedom
#print axioms planckConstant_real_pos
#print axioms cesiumHyperfine_real_pos
#print axioms ClockCalibration.labSeconds_def
#print axioms ClockCalibration.internalOfSeconds_def
#print axioms ClockCalibration.labFrequency_def
#print axioms ClockCalibration.labSeconds_internalOfSeconds
#print axioms ClockCalibration.internalOfSeconds_labSeconds
#print axioms ClockCalibration.tau_mul_labFrequency
#print axioms ClockCalibration.labFrequency_mul_labPeriod
#print axioms ClockCalibration.toAffineMap_map
#print axioms ClockCalibration.comp_tau
#print axioms ClockCalibration.toAffineMap_comp
#print axioms ClockCalibration.toAffineMap_roundtrip
#print axioms hbarOverTau_eq
#print axioms hbarOverTau_pos
#print axioms internalOfJoules_energyJoules
#print axioms energyJoules_internalOfJoules
#print axioms calibrationImport_anchors_unique
#print axioms stepTime_eq_labSeconds
#print axioms stepClock_periodic_iff_labSeconds
#print axioms labFrequency_of_angular_rate
#print axioms stepClock_lab_frequency
#print axioms instrumentStepLabDuration_eq_tau
#print axioms bundleStepLabDuration_eq_tau
#print axioms instrumentStepLabDuration_step_independent
#print axioms certified_step_not_a_duration
#print axioms fourFifthsBundle_lab_step_duration
#print axioms internalClock_lab_reading
#print axioms internalClock_lab_period_receipt
#print axioms labFrequency_not_forced
#print axioms unitTick_tau
#print axioms doubleTick_tau
#print axioms tick_declaration_not_forced
#print axioms cesiumTick_tau
#print axioms cesiumTick_unit_rate_frequency
#print axioms cesiumTick_unit_rate_frequency_anchor
#print axioms cesiumTick_matches_anchor
#print axioms cesiumTick_energy_quantum
#print axioms workedImports_differ_only_in_tick

end

end OPH.PhysicalCalibrationImport
