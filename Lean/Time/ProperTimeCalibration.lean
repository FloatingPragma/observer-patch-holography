import Time.WorldlineRealization
import Time.ClockReadout

/-!
# D1: proper-time calibration on a supplied affine history

An `AffineTimelikeHistory` already contains the substantive calibration
premise: its clock parameter multiplies a future-unit timelike direction in a
supplied event chart.  Under exactly that premise, the positive clock
increment is the positive square root selected by the Lorentz interval, is
chart independent wherever the endpoints overlap, and is additive along a
record chain.

This is a conditional calibration theorem.  It does not construct the atlas,
history, unit direction, affine coordinate law, refinement limit, SI unit, or
physical clock instrument.
-/

namespace OPH.D1Time

noncomputable section

open OPH.TimeOrderLedger
open OPH.C1Lorentz
open OPH.C2Soldering

universe u v w

/-- A pair of records whose order is certified by one declared history. -/
structure OrderedRecordPair {Record : Type u}
    (records : ObserverRecordOrder Record) where
  first : Record
  second : Record
  before : records.precedes first second

namespace AffineTimelikeHistory

variable {Record : Type u} {Event : Type v} {Chart : Type w}
    (history : AffineTimelikeHistory Record Event Chart)

/-- On a unit-timelike affine history, the Lorentz square of a displacement is
the square of the supplied clock increment. -/
theorem lorentzQ_displacement_eq_clockIncrement_sq (first second : Record) :
    lorentzQ (history.atlas.displacement history.chart
      (history.eventAt first) (history.eventAt second)) =
      clockIncrement history.clock first second ^ 2 := by
  rw [history.displacement_eq first second, lorentzQ_smul,
    history.velocity.2.1]
  simp [clockIncrement]

/-- The calibrated proper-time interval selected by a declared precedence. -/
def properTimeBetween (pair : OrderedRecordPair history.records) : ProperTime where
  value := clockIncrement history.clock pair.first pair.second
  nonnegative := le_of_lt (clockIncrement_pos history.clock pair.before)

@[simp] theorem properTimeBetween_value
    (pair : OrderedRecordPair history.records) :
    (history.properTimeBetween pair).value =
      clockIncrement history.clock pair.first pair.second := rfl

/-- The calibrated interval is the positive Lorentz length of the selected
displacement, expressed without introducing an analytic square-root API. -/
theorem properTimeBetween_sq_eq_interval
    (pair : OrderedRecordPair history.records) :
    (history.properTimeBetween pair).value ^ 2 =
      lorentzQ (history.atlas.displacement history.chart
        (history.eventAt pair.first) (history.eventAt pair.second)) := by
  rw [history.lorentzQ_displacement_eq_clockIncrement_sq]
  rfl

/-- The same calibrated squared interval is obtained in every overlapping
chart where both endpoint events are visible. -/
theorem properTimeBetween_sq_eq_interval_in_chart
    (pair : OrderedRecordPair history.records) (other : Chart)
    (hFirst : history.atlas.visible other (history.eventAt pair.first))
    (hSecond : history.atlas.visible other (history.eventAt pair.second)) :
    (history.properTimeBetween pair).value ^ 2 =
      lorentzQ (history.atlas.displacement other
        (history.eventAt pair.first) (history.eventAt pair.second)) := by
  rw [history.atlas.interval_overlap
    (history.visible pair.first) (history.visible pair.second) hFirst hSecond]
  exact history.properTimeBetween_sq_eq_interval pair

/-- Proper-time increments add exactly along three ordered records on the
same supplied affine history. -/
theorem properTimeBetween_add {first second third : Record}
    (firstBeforeSecond : history.records.precedes first second)
    (secondBeforeThird : history.records.precedes second third) :
    (history.properTimeBetween
      ⟨first, third, history.records.trans firstBeforeSecond secondBeforeThird⟩).value =
      (history.properTimeBetween
        ⟨first, second, firstBeforeSecond⟩).value +
      (history.properTimeBetween
        ⟨second, third, secondBeforeThird⟩).value := by
  simp [properTimeBetween, clockIncrement]

/-- Expose interval calibration as an explicit named realization map.  The
map exists only after an `AffineTimelikeHistory` carrying the unit-speed
coordinate law has been supplied. -/
def properTimeCalibration :
    NamedRealizationMap (OrderedRecordPair history.records) ProperTime where
  label := "ordered-record-pair-to-affine-proper-time"
  realize := history.properTimeBetween

end AffineTimelikeHistory

#print axioms AffineTimelikeHistory.properTimeBetween_sq_eq_interval
#print axioms AffineTimelikeHistory.properTimeBetween_sq_eq_interval_in_chart
#print axioms AffineTimelikeHistory.properTimeBetween_add

end

end OPH.D1Time
