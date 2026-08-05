import Time.ObserverHistory

/-!
# D1: affine clock comparison algebra

This module isolates the exact algebra available before a clock has a physical
unit.  Positive affine maps preserve ordering, compose, and are invertible.
Two ordered comparison events determine a unique positive-affine interpolation
of their four supplied readout values.  A third event whose source and target
readings differ from the anchors gives a nondegenerate consistency check
without another fit parameter.

The module does not assert that every pair of physical clocks is affinely
related.  Event correspondence and physical calibration remain explicit D1
receipts.
-/

namespace OPH.D1Time

noncomputable section

open OPH.TimeOrderLedger

universe u

/-- Any strictly increasing numerical regrading preserves the bare clock
interface.  Therefore monotonicity alone leaves more freedom than affine
scale and origin. -/
def regradeClock {Record : Type u} {records : ObserverRecordOrder Record}
    (clock : ClockReadout records) (regrade : ℝ → ℝ)
    (regradeStrict : StrictMono regrade) : ClockReadout records where
  read := fun record => regrade (clock.read record)
  strictlyMonotone := fun hBefore =>
    regradeStrict (clock.strictlyMonotone hBefore)

@[simp] theorem regradeClock_read {Record : Type u}
    {records : ObserverRecordOrder Record} (clock : ClockReadout records)
    (regrade : ℝ → ℝ) (regradeStrict : StrictMono regrade) (record : Record) :
    (regradeClock clock regrade regradeStrict).read record =
      regrade (clock.read record) := rfl

/-- Difference of two readings of one supplied clock. -/
def clockIncrement {Record : Type u} {records : ObserverRecordOrder Record}
    (clock : ClockReadout records) (first second : Record) : ℝ :=
  clock.read second - clock.read first

/-- A declared record precedence has positive clock increment. -/
theorem clockIncrement_pos {Record : Type u}
    {records : ObserverRecordOrder Record} (clock : ClockReadout records)
    {first second : Record} (hBefore : records.precedes first second) :
    0 < clockIncrement clock first second :=
  sub_pos.mpr (clock.strictlyMonotone hBefore)

/-- Affine regraduation multiplies intervals by its positive scale and removes
the origin offset. -/
theorem affineGauge_clockIncrement {Record : Type u}
    {records : ObserverRecordOrder Record} (clock : ClockReadout records)
    (scale offset : ℝ) (scalePositive : 0 < scale) (first second : Record) :
    clockIncrement (clock.affineGauge scale offset scalePositive) first second =
      scale * clockIncrement clock first second := by
  simp [clockIncrement, ClockReadout.affineGauge_read]
  ring

/-- A positive affine map between numerical clock readouts.  Application is
kept explicit through `.map`; there is deliberately no function coercion. -/
structure PositiveAffineClockMap where
  scale : ℝ
  offset : ℝ
  scalePositive : 0 < scale

namespace PositiveAffineClockMap

/-- Apply a declared clock comparison map. -/
def map (calibration : PositiveAffineClockMap) (reading : ℝ) : ℝ :=
  calibration.scale * reading + calibration.offset

/-- Identity clock comparison. -/
def identity : PositiveAffineClockMap where
  scale := 1
  offset := 0
  scalePositive := by norm_num

/-- Compose two declared clock comparisons. -/
def comp (second first : PositiveAffineClockMap) : PositiveAffineClockMap where
  scale := second.scale * first.scale
  offset := second.scale * first.offset + second.offset
  scalePositive := mul_pos second.scalePositive first.scalePositive

/-- Invert a positive affine clock comparison. -/
def inverse (calibration : PositiveAffineClockMap) : PositiveAffineClockMap where
  scale := calibration.scale⁻¹
  offset := -(calibration.scale⁻¹ * calibration.offset)
  scalePositive := inv_pos.mpr calibration.scalePositive

@[simp] theorem identity_map (reading : ℝ) : identity.map reading = reading := by
  simp [identity, map]

@[simp] theorem comp_map (second first : PositiveAffineClockMap)
    (reading : ℝ) :
    (comp second first).map reading = second.map (first.map reading) := by
  simp [comp, map]
  ring

@[simp] theorem inverse_map_map (calibration : PositiveAffineClockMap)
    (reading : ℝ) : calibration.inverse.map (calibration.map reading) = reading := by
  have hScale : calibration.scale ≠ 0 := ne_of_gt calibration.scalePositive
  simp [inverse, map]
  field_simp [hScale]
  ring

@[simp] theorem map_inverse_map (calibration : PositiveAffineClockMap)
    (reading : ℝ) : calibration.map (calibration.inverse.map reading) = reading := by
  have hScale : calibration.scale ≠ 0 := ne_of_gt calibration.scalePositive
  simp [inverse, map]
  field_simp [hScale]
  ring

/-- The unique positive affine candidate through two strictly ordered pairs of
readings. -/
def throughTwoPoints (sourceFirst sourceSecond targetFirst targetSecond : ℝ)
    (sourceOrdered : sourceFirst < sourceSecond)
    (targetOrdered : targetFirst < targetSecond) : PositiveAffineClockMap where
  scale := (targetSecond - targetFirst) / (sourceSecond - sourceFirst)
  offset := targetFirst -
    ((targetSecond - targetFirst) / (sourceSecond - sourceFirst)) * sourceFirst
  scalePositive := div_pos (sub_pos.mpr targetOrdered)
    (sub_pos.mpr sourceOrdered)

theorem throughTwoPoints_map_first
    (sourceFirst sourceSecond targetFirst targetSecond : ℝ)
    (sourceOrdered : sourceFirst < sourceSecond)
    (targetOrdered : targetFirst < targetSecond) :
    (throughTwoPoints sourceFirst sourceSecond targetFirst targetSecond
      sourceOrdered targetOrdered).map sourceFirst = targetFirst := by
  simp [throughTwoPoints, map]

theorem throughTwoPoints_map_second
    (sourceFirst sourceSecond targetFirst targetSecond : ℝ)
    (sourceOrdered : sourceFirst < sourceSecond)
    (targetOrdered : targetFirst < targetSecond) :
    (throughTwoPoints sourceFirst sourceSecond targetFirst targetSecond
      sourceOrdered targetOrdered).map sourceSecond = targetSecond := by
  have hSource : sourceSecond - sourceFirst ≠ 0 :=
    ne_of_gt (sub_pos.mpr sourceOrdered)
  simp [throughTwoPoints, map]
  field_simp [hSource]
  ring

/-- Two distinct source readings determine an affine clock map uniquely. -/
theorem eq_of_map_eq_at_two_points (first second : PositiveAffineClockMap)
    {sourceFirst sourceSecond targetFirst targetSecond : ℝ}
    (sourceNe : sourceFirst ≠ sourceSecond)
    (firstAtFirst : first.map sourceFirst = targetFirst)
    (firstAtSecond : first.map sourceSecond = targetSecond)
    (secondAtFirst : second.map sourceFirst = targetFirst)
    (secondAtSecond : second.map sourceSecond = targetSecond) :
    first = second := by
  have hScale : first.scale = second.scale := by
    have hProduct :
        (first.scale - second.scale) * (sourceSecond - sourceFirst) = 0 := by
      dsimp [map] at firstAtFirst firstAtSecond secondAtFirst secondAtSecond
      nlinarith
    rcases mul_eq_zero.mp hProduct with hEqual | hSource
    · linarith
    · exfalso
      apply sourceNe
      linarith
  have hOffset : first.offset = second.offset := by
    dsimp [map] at firstAtFirst secondAtFirst
    rw [hScale] at firstAtFirst
    linarith
  cases first
  cases second
  simp_all

/-- The two-point construction is the unique positive affine comparison that
matches both ordered event readings. -/
theorem throughTwoPoints_unique
    (candidate : PositiveAffineClockMap)
    (sourceFirst sourceSecond targetFirst targetSecond : ℝ)
    (sourceOrdered : sourceFirst < sourceSecond)
    (targetOrdered : targetFirst < targetSecond)
    (candidateFirst : candidate.map sourceFirst = targetFirst)
    (candidateSecond : candidate.map sourceSecond = targetSecond) :
    candidate = throughTwoPoints sourceFirst sourceSecond targetFirst targetSecond
      sourceOrdered targetOrdered := by
  apply eq_of_map_eq_at_two_points candidate _ (ne_of_lt sourceOrdered)
    candidateFirst candidateSecond
  · exact throughTwoPoints_map_first _ _ _ _ sourceOrdered targetOrdered
  · exact throughTwoPoints_map_second _ _ _ _ sourceOrdered targetOrdered

end PositiveAffineClockMap

/-! ## Non-affine negative control -/

/-- Cubing the canonical three-tick clock is another valid monotone clock. -/
def cubicThreeTickClock : ClockReadout (finRankedHistory 3).order :=
  regradeClock (finRankedHistory 3).clock (fun reading : ℝ => reading ^ 3)
    (by exact (show Odd 3 by decide).strictMono_pow)

@[simp] theorem cubicThreeTickClock_read (record : Fin 3) :
    cubicThreeTickClock.read record = (record.1 : ℝ) ^ 3 := rfl

/-- The cubic three-tick clock is not any positive affine regraduation of the
canonical clock.  Thus an additive/metric calibration receipt is required to
reduce the bare monotone gauge to the affine subgroup. -/
theorem cubicThreeTickClock_not_affine :
    ¬ ∃ calibration : PositiveAffineClockMap,
      ∀ record : Fin 3,
        calibration.map ((finRankedHistory 3).clock.read record) =
          cubicThreeTickClock.read record := by
  rintro ⟨calibration, hCalibration⟩
  have hZero := hCalibration (0 : Fin 3)
  have hOne := hCalibration (1 : Fin 3)
  have hTwo := hCalibration (2 : Fin 3)
  norm_num [PositiveAffineClockMap.map] at hZero hOne hTwo
  nlinarith

#print axioms affineGauge_clockIncrement
#print axioms PositiveAffineClockMap.throughTwoPoints_unique
#print axioms cubicThreeTickClock_not_affine

end

end OPH.D1Time
