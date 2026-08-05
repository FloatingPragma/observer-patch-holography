import Time.ProperTimeCalibration

/-!
# D1: public comparison of clocks at shared events

Two histories can compare readings only after records on both histories are
identified with the same public events.  One shared event fixes neither rate
nor origin.  Two ordered shared events determine the unique positive-affine
interpolation of the four supplied readings.  A third shared event supplies an
exact cross-multiplication equation.  Distinction from both anchor events and
readings makes that equation a nondegenerate no-new-fit-parameter check.
Calling it held out additionally requires a predesignation and custody
protocol, which this module does not construct.

This packet does not assert that arbitrary physical clocks have an affine
relation.  It makes that hypothesis testable and keeps the shared-event
correspondence as supplied evidence.
-/

namespace OPH.D1Time

noncomputable section

open OPH.TimeOrderLedger
open OPH.C1Lorentz
open OPH.C2Soldering

universe uLeft uRight v wLeft wRight

namespace PositiveAffineClockMap

/-- A positive affine clock map of any chosen rate through one comparison
event. -/
def throughOnePoint (source target rate : ℝ) (ratePositive : 0 < rate) :
    PositiveAffineClockMap where
  scale := rate
  offset := target - rate * source
  scalePositive := ratePositive

@[simp] theorem throughOnePoint_map (source target rate : ℝ)
    (ratePositive : 0 < rate) :
    (throughOnePoint source target rate ratePositive).map source = target := by
  simp [throughOnePoint, map]

/-- One shared event cannot determine both affine clock parameters. -/
theorem onePoint_not_unique (source target : ℝ) :
    ∃ first second : PositiveAffineClockMap,
      first ≠ second ∧ first.map source = target ∧ second.map source = target := by
  refine ⟨throughOnePoint source target 1 (by norm_num),
    throughOnePoint source target 2 (by norm_num), ?_, by simp, by simp⟩
  intro hEqual
  have hScale := congrArg PositiveAffineClockMap.scale hEqual
  norm_num [throughOnePoint] at hScale

end PositiveAffineClockMap

/-- Records on two supplied affine histories that represent two common public
events, together with their local precedence certificates. -/
structure TwoSharedEventComparison
    {LeftRecord : Type uLeft} {RightRecord : Type uRight}
    {Event : Type v} {LeftChart : Type wLeft} {RightChart : Type wRight}
    (left : AffineTimelikeHistory LeftRecord Event LeftChart)
    (right : AffineTimelikeHistory RightRecord Event RightChart) where
  leftFirst : LeftRecord
  leftSecond : LeftRecord
  rightFirst : RightRecord
  rightSecond : RightRecord
  firstShared : left.eventAt leftFirst = right.eventAt rightFirst
  secondShared : left.eventAt leftSecond = right.eventAt rightSecond
  leftBefore : left.records.precedes leftFirst leftSecond
  rightBefore : right.records.precedes rightFirst rightSecond

namespace TwoSharedEventComparison

variable {LeftRecord : Type uLeft} {RightRecord : Type uRight}
    {Event : Type v} {LeftChart : Type wLeft} {RightChart : Type wRight}
    {left : AffineTimelikeHistory LeftRecord Event LeftChart}
    {right : AffineTimelikeHistory RightRecord Event RightChart}
    (comparison : TwoSharedEventComparison left right)

/-- The unique positive affine comparison fitted to the two shared events. -/
def calibration : PositiveAffineClockMap :=
  PositiveAffineClockMap.throughTwoPoints
    (left.clock.read comparison.leftFirst)
    (left.clock.read comparison.leftSecond)
    (right.clock.read comparison.rightFirst)
    (right.clock.read comparison.rightSecond)
    (left.clock.strictlyMonotone comparison.leftBefore)
    (right.clock.strictlyMonotone comparison.rightBefore)

@[simp] theorem calibration_map_first :
    comparison.calibration.map (left.clock.read comparison.leftFirst) =
      right.clock.read comparison.rightFirst :=
  PositiveAffineClockMap.throughTwoPoints_map_first _ _ _ _
    (left.clock.strictlyMonotone comparison.leftBefore)
    (right.clock.strictlyMonotone comparison.rightBefore)

@[simp] theorem calibration_map_second :
    comparison.calibration.map (left.clock.read comparison.leftSecond) =
      right.clock.read comparison.rightSecond :=
  PositiveAffineClockMap.throughTwoPoints_map_second _ _ _ _
    (left.clock.strictlyMonotone comparison.leftBefore)
    (right.clock.strictlyMonotone comparison.rightBefore)

/-- No other positive affine comparison matches both shared-event readings. -/
theorem calibration_unique (candidate : PositiveAffineClockMap)
    (candidateFirst :
      candidate.map (left.clock.read comparison.leftFirst) =
        right.clock.read comparison.rightFirst)
    (candidateSecond :
      candidate.map (left.clock.read comparison.leftSecond) =
        right.clock.read comparison.rightSecond) :
    candidate = comparison.calibration :=
  PositiveAffineClockMap.throughTwoPoints_unique candidate _ _ _ _
    (left.clock.strictlyMonotone comparison.leftBefore)
    (right.clock.strictlyMonotone comparison.rightBefore)
    candidateFirst candidateSecond

end TwoSharedEventComparison

/-- A third pair of records certified to represent the same public event. -/
structure ThirdSharedEvent
    {LeftRecord : Type uLeft} {RightRecord : Type uRight}
    {Event : Type v} {LeftChart : Type wLeft} {RightChart : Type wRight}
    {left : AffineTimelikeHistory LeftRecord Event LeftChart}
    {right : AffineTimelikeHistory RightRecord Event RightChart}
    (comparison : TwoSharedEventComparison left right) where
  leftRecord : LeftRecord
  rightRecord : RightRecord
  shared : left.eventAt leftRecord = right.eventAt rightRecord

namespace ThirdSharedEvent

variable {LeftRecord : Type uLeft} {RightRecord : Type uRight}
    {Event : Type v} {LeftChart : Type wLeft} {RightChart : Type wRight}
    {left : AffineTimelikeHistory LeftRecord Event LeftChart}
    {right : AffineTimelikeHistory RightRecord Event RightChart}
    {comparison : TwoSharedEventComparison left right}
    (third : ThirdSharedEvent comparison)

/-- The third shared event agrees with the two-event affine calibration. -/
def IsAffineConsistent : Prop :=
  comparison.calibration.map (left.clock.read third.leftRecord) =
    right.clock.read third.rightRecord

/-- The third comparison event and both of its clock readings are distinct
from the two anchors.  Event distinction alone is insufficient because the
supplied event maps and clocks need not be injective on incomparable records.
The cross-multiplication theorem does not need these premises.  They make the
check algebraically nondegenerate; they do not certify predesignation, custody,
or non-use in constructing either history. -/
def IsNondegenerateThirdPoint : Prop :=
  left.eventAt third.leftRecord ≠ left.eventAt comparison.leftFirst ∧
    left.eventAt third.leftRecord ≠ left.eventAt comparison.leftSecond ∧
    left.clock.read third.leftRecord ≠
      left.clock.read comparison.leftFirst ∧
    left.clock.read third.leftRecord ≠
      left.clock.read comparison.leftSecond ∧
    right.clock.read third.rightRecord ≠
      right.clock.read comparison.rightFirst ∧
    right.clock.read third.rightRecord ≠
      right.clock.read comparison.rightSecond

/-- Exact third-event cross-multiplication equation.  With
`IsNondegenerateThirdPoint`, this compares two interval ratios without a new
fitted parameter.  A prospective held-out interpretation needs an additional
predesignation and custody protocol. -/
theorem affineConsistent_iff_crossMultiplication :
    third.IsAffineConsistent ↔
      (right.clock.read third.rightRecord -
          right.clock.read comparison.rightFirst) *
        (left.clock.read comparison.leftSecond -
          left.clock.read comparison.leftFirst) =
      (right.clock.read comparison.rightSecond -
          right.clock.read comparison.rightFirst) *
        (left.clock.read third.leftRecord -
          left.clock.read comparison.leftFirst) := by
  have hLeft :
      left.clock.read comparison.leftSecond -
        left.clock.read comparison.leftFirst ≠ 0 :=
    ne_of_gt (sub_pos.mpr
      (left.clock.strictlyMonotone comparison.leftBefore))
  dsimp [IsAffineConsistent, TwoSharedEventComparison.calibration,
    PositiveAffineClockMap.throughTwoPoints,
    PositiveAffineClockMap.map]
  field_simp [hLeft]
  constructor <;> intro h <;> nlinarith

end ThirdSharedEvent

/-! ## Executable three-event controls -/

/-- A one-chart algebraic atlas whose coordinate of each public event is its
supplied clock reading along the standard future-unit direction. -/
def threeRecordAtlasOfClock
    (clock : ClockReadout (finRankedHistory 3).order) :
    EventGermAtlas (Fin 3) Unit where
  visible := fun _ _ => True
  coordinate := fun _ record => clock.read record • (standardFrame : Herm2)
  overlap := identityOverlapCocycle Unit
  coordinate_overlap := by
    intro i j event _ _
    simp [LorentzOverlapCocycle.act, identityOverlapCocycle,
      OrientedLorentzEquiv.refl]

/-- Every supplied three-tick clock yields a finite affine control history in
its own declared atlas.  The atlas is constructed test data, not a source or
physical clock realization. -/
def threeRecordHistoryOfClock
    (clock : ClockReadout (finRankedHistory 3).order) :
    AffineTimelikeHistory (Fin 3) (Fin 3) Unit where
  records := (finRankedHistory 3).order
  clock := clock
  atlas := threeRecordAtlasOfClock clock
  chart := ()
  eventAt := id
  origin := 0
  velocity := standardFrame
  visible := by simp [threeRecordAtlasOfClock]
  coordinate_eq := by simp [threeRecordAtlasOfClock]

@[simp] theorem threeRecordHistoryOfClock_read
    (clock : ClockReadout (finRankedHistory 3).order) (record : Fin 3) :
    (threeRecordHistoryOfClock clock).clock.read record = clock.read record :=
  rfl

/-- A nontrivial positive affine regrading of the canonical three-tick
clock, used by the passing algebraic control. -/
def affineThreeTickClock : ClockReadout (finRankedHistory 3).order :=
  (finRankedHistory 3).clock.affineGauge 2 5 (by norm_num)

def affineControlComparison :
    TwoSharedEventComparison
      (threeRecordHistoryOfClock (finRankedHistory 3).clock)
      (threeRecordHistoryOfClock affineThreeTickClock) where
  leftFirst := 0
  leftSecond := 1
  rightFirst := 0
  rightSecond := 1
  firstShared := rfl
  secondShared := rfl
  leftBefore := by
    change (0 : Fin 3) < 1
    decide
  rightBefore := by
    change (0 : Fin 3) < 1
    decide

/-- The final record is a public event distinct from both fitted anchors. -/
def affineControlThird : ThirdSharedEvent affineControlComparison where
  leftRecord := 2
  rightRecord := 2
  shared := rfl

theorem affineControlThird_isNondegenerate :
    affineControlThird.IsNondegenerateThirdPoint := by
  norm_num [ThirdSharedEvent.IsNondegenerateThirdPoint, affineControlThird,
    affineControlComparison, threeRecordHistoryOfClock, affineThreeTickClock,
    RankedHistory.clock_read]
  all_goals decide

theorem affineControlThird_isConsistent :
    affineControlThird.IsAffineConsistent := by
  rw [affineControlThird.affineConsistent_iff_crossMultiplication]
  norm_num [affineControlThird, affineControlComparison,
    affineThreeTickClock, RankedHistory.clock_read]

/-- The same anchors can be fitted when the cubic clock supplies the second
readout, while the distinct third event rejects the affine relation. -/
def nonAffineControlComparison :
    TwoSharedEventComparison
      (threeRecordHistoryOfClock (finRankedHistory 3).clock)
      (threeRecordHistoryOfClock cubicThreeTickClock) where
  leftFirst := 0
  leftSecond := 1
  rightFirst := 0
  rightSecond := 1
  firstShared := rfl
  secondShared := rfl
  leftBefore := by
    change (0 : Fin 3) < 1
    decide
  rightBefore := by
    change (0 : Fin 3) < 1
    decide

def nonAffineControlThird : ThirdSharedEvent nonAffineControlComparison where
  leftRecord := 2
  rightRecord := 2
  shared := rfl

theorem nonAffineControlThird_isNondegenerate :
    nonAffineControlThird.IsNondegenerateThirdPoint := by
  norm_num [ThirdSharedEvent.IsNondegenerateThirdPoint, nonAffineControlThird,
    nonAffineControlComparison, threeRecordHistoryOfClock,
    cubicThreeTickClock, regradeClock,
    RankedHistory.clock_read]
  all_goals decide

theorem nonAffineControlThird_isRejected :
    ¬ nonAffineControlThird.IsAffineConsistent := by
  rw [nonAffineControlThird.affineConsistent_iff_crossMultiplication]
  norm_num [nonAffineControlThird, nonAffineControlComparison,
    cubicThreeTickClock, regradeClock, RankedHistory.clock_read]

#print axioms PositiveAffineClockMap.onePoint_not_unique
#print axioms TwoSharedEventComparison.calibration_unique
#print axioms ThirdSharedEvent.affineConsistent_iff_crossMultiplication
#print axioms affineControlThird_isConsistent
#print axioms nonAffineControlThird_isRejected

end

end OPH.D1Time
