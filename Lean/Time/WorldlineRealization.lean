import Geometry.EventGermDisplacement
import Geometry.ObserverFrameHyperboloid
import Time.ObserverHistory

/-!
# D1: bounded affine-timelike history realization

This file proves a conditional finite realization theorem on the C2 event-germ
interface.  If a supplied record clock parametrizes event coordinates along a
supplied future-unit timelike direction in one supplied chart, then every
precedence is represented by a future-timelike displacement.  The resulting
map inhabits A1's explicit `WorldlineRealization` interface.

The atlas, events, visibility, monotone clock, affine coordinate law, and unit
timelike direction are premises.  The theorem neither constructs them from
observer repair data nor promotes the algebraic `Herm2` carrier to physical
spacetime.
-/

namespace OPH.D1Time

noncomputable section

open OPH.TimeOrderLedger
open OPH.C1Lorentz
open OPH.C2Soldering

universe u v w

/-! ## A nonvacuous algebraic control model -/

/-- Identity affine-overlap data, used only for the explicit bounded controls
below. -/
def identityOverlapCocycle (Chart : Type*) : LorentzOverlapCocycle Chart where
  lorentz := fun _ _ => OrientedLorentzEquiv.refl
  translation := fun _ _ => 0
  lorentz_self := by simp [OrientedLorentzEquiv.refl]
  lorentz_cocycle := by simp [OrientedLorentzEquiv.refl]
  translation_self := by simp
  translation_cocycle := by simp [OrientedLorentzEquiv.refl]

/-- The one-chart atlas on the intrinsic `Herm2` carrier.  This is an
algebraic test model, not a source-produced event population. -/
def identityEventAtlas : EventGermAtlas Herm2 Unit where
  visible := fun _ _ => True
  coordinate := fun _ event => event
  overlap := identityOverlapCocycle Unit
  coordinate_overlap := by
    intro i j event _ _
    simp [identityOverlapCocycle, LorentzOverlapCocycle.act,
      OrientedLorentzEquiv.refl]

/-- Algebraic future timelikeness in the declared C1 Lorentz module. -/
def IsFutureTimelike (displacement : Herm2) : Prop :=
  0 < lorentzQ displacement ∧ 0 < displacement.1

/-- The Lorentz square of a scalar multiple. -/
theorem lorentzQ_smul (scale : ℝ) (vector : Herm2) :
    lorentzQ (scale • vector) = scale ^ 2 * lorentzQ vector := by
  calc
    lorentzQ (scale • vector) =
        lorentzB (scale • vector) (scale • vector) :=
      (lorentzB_self _).symm
    _ = scale * lorentzB vector (scale • vector) :=
      lorentzB_smul_left scale vector (scale • vector)
    _ = scale * (scale * lorentzB vector vector) := by
      rw [lorentzB_smul_right]
    _ = scale ^ 2 * lorentzQ vector := by
      rw [lorentzB_self]
      ring

/-- A positive multiple of a future-unit timelike vector is future timelike. -/
theorem positive_smul_futureTimelike (scale : ℝ) (velocity : FrameHyperboloid)
    (scalePositive : 0 < scale) :
    IsFutureTimelike (scale • (velocity : Herm2)) := by
  constructor
  · rw [lorentzQ_smul, velocity.2.1]
    nlinarith
  · change 0 < scale * velocity.1.1
    exact mul_pos scalePositive velocity.2.2

/-- Future-timelike reachability in one explicitly selected C2 chart. -/
def ChartFutureTimelike {Event Chart : Type*}
    (atlas : EventGermAtlas Event Chart) (chart : Chart)
    (first second : Event) : Prop :=
  atlas.visible chart first ∧ atlas.visible chart second ∧
    IsFutureTimelike (atlas.displacement chart first second)

/-- A supplied finite or infinite record history whose selected chart
coordinates follow one affine future-unit timelike line, parametrized by a
supplied monotone clock. -/
structure AffineTimelikeHistory (Record : Type u) (Event : Type v)
    (Chart : Type w) where
  records : ObserverRecordOrder Record
  clock : ClockReadout records
  atlas : EventGermAtlas Event Chart
  chart : Chart
  eventAt : Record → Event
  origin : Herm2
  velocity : FrameHyperboloid
  visible : ∀ record, atlas.visible chart (eventAt record)
  coordinate_eq : ∀ record,
    atlas.coordinate chart (eventAt record) =
      origin + clock.read record • (velocity : Herm2)

/-- Every canonical finite record chain has a literal future-unit affine
realization in the one-chart algebraic control model. -/
def finAffineTimelikeHistory (n : ℕ) :
    AffineTimelikeHistory (Fin n) Herm2 Unit where
  records := (finRankedHistory n).order
  clock := (finRankedHistory n).clock
  atlas := identityEventAtlas
  chart := ()
  eventAt record :=
    (finRankedHistory n).clock.read record • (standardFrame : Herm2)
  origin := 0
  velocity := standardFrame
  visible := by simp [identityEventAtlas]
  coordinate_eq := by simp [identityEventAtlas]

namespace AffineTimelikeHistory

variable {Record : Type u} {Event : Type v} {Chart : Type w}
    (history : AffineTimelikeHistory Record Event Chart)

/-- Exact coordinate displacement along the affine history. -/
theorem displacement_eq (first second : Record) :
    history.atlas.displacement history.chart
        (history.eventAt first) (history.eventAt second) =
      (history.clock.read second - history.clock.read first) •
        (history.velocity : Herm2) := by
  rw [EventGermAtlas.displacement, history.coordinate_eq,
    history.coordinate_eq]
  rw [sub_smul]
  abel

/-- Every declared record precedence has positive clock increment. -/
theorem clock_increment_pos {first second : Record}
    (hBefore : history.records.precedes first second) :
    0 < history.clock.read second - history.clock.read first :=
  sub_pos.mpr (history.clock.strictlyMonotone hBefore)

/-- Every declared precedence is future timelike in the selected chart. -/
theorem displacement_futureTimelike {first second : Record}
    (hBefore : history.records.precedes first second) :
    IsFutureTimelike
      (history.atlas.displacement history.chart
        (history.eventAt first) (history.eventAt second)) := by
  rw [history.displacement_eq first second]
  exact positive_smul_futureTimelike _ history.velocity
    (history.clock_increment_pos hBefore)

/-- The affine chain supplies A1's explicit order-preserving event map, with
causality restricted to the selected chart's future-timelike relation. -/
def worldline : WorldlineRealization history.records
    (ChartFutureTimelike history.atlas history.chart) where
  eventAt := history.eventAt
  orderPreserving := by
    intro first second hBefore
    exact ⟨history.visible first, history.visible second,
      history.displacement_futureTimelike hBefore⟩

/-- Timelikeness of the same displacement is preserved in any second chart
where both endpoint events are visible.  This is derived from the C2 overlap
law and the transported future-unit direction. -/
theorem displacement_futureTimelike_in_chart {first second : Record}
    (hBefore : history.records.precedes first second) (other : Chart)
    (hFirst : history.atlas.visible other (history.eventAt first))
    (hSecond : history.atlas.visible other (history.eventAt second)) :
    IsFutureTimelike
      (history.atlas.displacement other
        (history.eventAt first) (history.eventAt second)) := by
  let scale := history.clock.read second - history.clock.read first
  have hScale : 0 < scale := history.clock_increment_pos hBefore
  have hOverlap := history.atlas.displacement_overlap
    (history.visible first) (history.visible second) hFirst hSecond
  rw [hOverlap, history.displacement_eq first second,
    OrientedLorentzEquiv.map_smul]
  exact positive_smul_futureTimelike scale
    ((history.atlas.overlap.lorentz history.chart other).mapFrame
      history.velocity) hScale

end AffineTimelikeHistory

/-- The three-record control realizes two successive future-timelike
displacements. -/
theorem threeRecord_worldline_control :
    ChartFutureTimelike identityEventAtlas ()
        ((finAffineTimelikeHistory 3).eventAt 0)
        ((finAffineTimelikeHistory 3).eventAt 1) ∧
      ChartFutureTimelike identityEventAtlas ()
        ((finAffineTimelikeHistory 3).eventAt 1)
        ((finAffineTimelikeHistory 3).eventAt 2) := by
  constructor
  · exact ⟨by simp [identityEventAtlas], by simp [identityEventAtlas],
      (finAffineTimelikeHistory 3).displacement_futureTimelike
        (show (0 : Fin 3) < 1 by decide)⟩
  · exact ⟨by simp [identityEventAtlas], by simp [identityEventAtlas],
      (finAffineTimelikeHistory 3).displacement_futureTimelike
        (show (1 : Fin 3) < 2 by decide)⟩

/-! ## Nonvacuous finite two-chart control -/

/-- Three events with identical coordinates in two declared charts.  This is
a bounded witness for the conditional interface, not a source-produced event
atlas. -/
def threeRecordTwoChartAtlas : EventGermAtlas (Fin 3) (Fin 2) where
  visible := fun _ _ => True
  coordinate := fun _ record => ((record.1 : ℝ), 0)
  overlap := identityOverlapCocycle (Fin 2)
  coordinate_overlap := by
    intro i j event _ _
    simp [LorentzOverlapCocycle.act, identityOverlapCocycle,
      OrientedLorentzEquiv.refl]

/-- A literal three-record affine unit-timelike history visible in both
charts.  All geometric data are explicit finite control data. -/
def threeRecordAffineHistory :
    AffineTimelikeHistory (Fin 3) (Fin 3) (Fin 2) where
  records := (finRankedHistory 3).order
  clock := (finRankedHistory 3).clock
  atlas := threeRecordTwoChartAtlas
  chart := 0
  eventAt := id
  origin := 0
  velocity := standardFrame
  visible := by simp [threeRecordTwoChartAtlas]
  coordinate_eq := by
    intro record
    ext <;> simp [threeRecordTwoChartAtlas, RankedHistory.clock,
      standardFrame]

/-- The first two records of the explicit control have a future-timelike
displacement even when read in the second chart. -/
theorem threeRecord_twoChart_futureTimelike :
    IsFutureTimelike
      (threeRecordAffineHistory.atlas.displacement 1
        (threeRecordAffineHistory.eventAt 0)
        (threeRecordAffineHistory.eventAt 1)) := by
  apply threeRecordAffineHistory.displacement_futureTimelike_in_chart
    (first := 0) (second := 1)
  · change (0 : Fin 3) < 1
    decide
  · simp [threeRecordAffineHistory, threeRecordTwoChartAtlas]
  · simp [threeRecordAffineHistory, threeRecordTwoChartAtlas]

#print axioms positive_smul_futureTimelike
#print axioms AffineTimelikeHistory.displacement_eq
#print axioms AffineTimelikeHistory.worldline
#print axioms AffineTimelikeHistory.displacement_futureTimelike_in_chart
#print axioms threeRecord_worldline_control
#print axioms threeRecord_twoChart_futureTimelike

end

end OPH.D1Time
