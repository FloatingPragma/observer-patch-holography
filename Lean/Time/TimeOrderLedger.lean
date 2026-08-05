import Mathlib

/-!
# Time-and-order type ledger

This module enforces the ontology boundary required by completion-plan issue
`#675`.  Universe closure, repair scheduling, observer record order, modular
parameter, worldline realization, clock readout, proper time, and an optional
global time function are separate types.  None has a `Coe` instance to any
other ledger type.  A downstream construction that relates two layers must
provide an explicit `NamedRealizationMap`.

The structures are interfaces, not existence theorems.  In particular this
file constructs no worldline, clock calibration, proper-time limit, global
time function, or physical identification of modular flow.  The
`WorldlineRealization`, `ClockReadout`, and `GlobalTimeFunction` fields state
what a downstream construction must prove when it supplies an inhabitant.
-/

namespace OPH.TimeOrderLedger

universe u v

/-! ## Pairwise-distinct ledger entries -/

/-- Audit tags for the eight concepts kept separate by this module.  These
tags are documentation and machine-checkable inventory; the corresponding
Lean structures below carry the actual type boundary. -/
inductive LedgerKind where
  | universeClosure
  | repairOrder
  | observerRecordOrder
  | modularParameter
  | worldlineRealization
  | clockReadout
  | properTime
  | globalTimeFunction
  deriving DecidableEq, Repr

/-- The canonical A1 inventory, in the order used by the completion plan. -/
def canonicalLedgerKinds : List LedgerKind :=
  [.universeClosure, .repairOrder, .observerRecordOrder,
    .modularParameter, .worldlineRealization, .clockReadout,
    .properTime, .globalTimeFunction]

/-- The ledger inventory contains eight different tags. -/
theorem canonicalLedgerKinds_pairwise :
    canonicalLedgerKinds.Pairwise (fun left right => left ≠ right) := by
  decide

/-! ## Timeless closure and order data -/

/-- A supplied fixed point of a closure operation.  An inhabitant is a
closure witness, not a duration, execution counter, or clock.  This module
does not construct such a witness for the OPH public-world endpoint. -/
structure UniverseClosure (State : Type u) (close : State -> State) where
  state : State
  isClosed : close state = state

/-- A finite execution schedule for primitive repair moves.  Its position
index records only the order in which moves are attempted; no physical unit
or duration is attached. -/
structure RepairOrder (Move : Type u) where
  length : Nat
  moveAt : Fin length -> Move

namespace RepairOrder

/-- The empty repair schedule. -/
def empty (Move : Type u) : RepairOrder Move where
  length := 0
  moveAt := Fin.elim0

end RepairOrder

/-- A strict order on observer-visible records.  It is deliberately not
installed as an `LT` instance: code using the order must name the supplied
record-order object. -/
structure ObserverRecordOrder (Record : Type u) where
  precedes : Record -> Record -> Prop
  irrefl : forall record, ¬ precedes record record
  trans : forall {first second third},
    precedes first second -> precedes second third -> precedes first third

namespace ObserverRecordOrder

/-- The discrete observer record order, useful for boundary tests. -/
def discrete (Record : Type u) : ObserverRecordOrder Record where
  precedes := fun _ _ => False
  irrefl := by simp
  trans := by simp

end ObserverRecordOrder

/-! ## Local parameters and geometric/clock realizations -/

/-- The real parameter of a local modular automorphism family.  The wrapper
adds no clock, worldline, causal, or unit interpretation. -/
structure ModularParameter where
  value : Real

/-- An explicit order-preserving realization of an observer record history
as causal events.  Supplying an inhabitant is the downstream D1 obligation;
this A1 module supplies only the type. -/
structure WorldlineRealization {Record : Type u} {Event : Type v}
    (records : ObserverRecordOrder Record)
    (causal : Event -> Event -> Prop) where
  eventAt : Record -> Event
  orderPreserving : forall {first second}, records.precedes first second ->
    causal (eventAt first) (eventAt second)

/-- A real-valued record observable strictly increasing on the supplied
observer record order.  It is an ordering coordinate, not yet proper time. -/
structure ClockReadout {Record : Type u}
    (records : ObserverRecordOrder Record) where
  read : Record -> Real
  strictlyMonotone : forall {first second}, records.precedes first second ->
    read first < read second

/-- A nonnegative proper-time interval.  Producing this value from a clock
readout requires a separate calibration map and a geometric receipt. -/
structure ProperTime where
  value : Real
  nonnegative : 0 <= value

/-- A global time function on an explicitly supplied causal branch.  The
strict-causality field makes this stronger than an arbitrary event label.
No inhabitant is assumed; availability is represented by
`OptionalGlobalTimeFunction`. -/
structure GlobalTimeFunction (Event : Type u)
    (causal : Event -> Event -> Prop) where
  read : Event -> Real
  strictlyCausal : forall {first second}, causal first second ->
    read first < read second

/-- The optional stable-causality branch.  `none` is the default logical
possibility; A1 does not choose either constructor. -/
def OptionalGlobalTimeFunction (Event : Type u)
    (causal : Event -> Event -> Prop) : Type u :=
  Option (GlobalTimeFunction Event causal)

/-! ## Named maps are the only bridge API -/

/-- An explicit, labelled realization or calibration map.  There is no
`CoeFun` instance: applying the map also remains an explicit `.realize`
operation. -/
structure NamedRealizationMap (Source : Type u) (Target : Type v) where
  label : String
  realize : Source -> Target

namespace NamedRealizationMap

/-- Explicit composition of two named layer bridges. -/
def comp {Source : Type u} {Middle : Type v} {Target : Type*}
    (second : NamedRealizationMap Middle Target)
    (first : NamedRealizationMap Source Middle) :
    NamedRealizationMap Source Target where
  label := first.label ++ " ; " ++ second.label
  realize := second.realize ∘ first.realize

@[simp] theorem comp_apply {Source : Type u} {Middle : Type v} {Target : Type*}
    (second : NamedRealizationMap Middle Target)
    (first : NamedRealizationMap Source Middle) (source : Source) :
    (comp second first).realize source =
      second.realize (first.realize source) := rfl

end NamedRealizationMap

/-- Expose a supplied worldline realization as an explicit named map from
records to events. -/
def WorldlineRealization.asNamedMap {Record : Type u} {Event : Type v}
    {records : ObserverRecordOrder Record} {causal : Event -> Event -> Prop}
    (worldline : WorldlineRealization records causal) :
    NamedRealizationMap Record Event where
  label := "observer-record-to-event"
  realize := worldline.eventAt

/-- Package a supplied modular-to-worldline interpretation.  The function is
an argument, so no physical interpretation is inferred from the modular
parameter alone. -/
def namedModularRealization {Event : Type u}
    (realize : ModularParameter -> Event) :
    NamedRealizationMap ModularParameter Event where
  label := "modular-parameter-to-realized-event"
  realize := realize

/-- Package a supplied clock calibration.  This is the only A1 API from an
ordering-coordinate readout to a proper-time interval. -/
def namedProperTimeCalibration {Record : Type u}
    {records : ObserverRecordOrder Record}
    (calibrate : ClockReadout records -> ProperTime) :
    NamedRealizationMap (ClockReadout records) ProperTime where
  label := "clock-readout-to-proper-time"
  realize := calibrate

/-! ## Exact affine-gauge boundary for an uncalibrated clock -/

namespace ClockReadout

/-- Positive affine regraduation preserves strict record monotonicity. -/
def affineGauge {Record : Type u} {records : ObserverRecordOrder Record}
    (clock : ClockReadout records) (scale offset : Real)
    (scalePositive : 0 < scale) : ClockReadout records where
  read := fun record => scale * clock.read record + offset
  strictlyMonotone := by
    intro first second hBefore
    simpa [add_comm] using
      add_lt_add_right
        (mul_lt_mul_of_pos_left (clock.strictlyMonotone hBefore) scalePositive)
        offset

@[simp] theorem affineGauge_read {Record : Type u}
    {records : ObserverRecordOrder Record}
    (clock : ClockReadout records) (scale offset : Real)
    (scalePositive : 0 < scale) (record : Record) :
    (affineGauge clock scale offset scalePositive).read record =
      scale * clock.read record + offset := rfl

/-- If at least one record exists, the same record order admits a different
clock origin.  Therefore record order plus monotonicity cannot calibrate a
unique clock. -/
theorem offsetGauge_ne {Record : Type u}
    {records : ObserverRecordOrder Record}
    (clock : ClockReadout records) (record : Record) (offset : Real)
    (offsetNonzero : offset ≠ 0) :
    affineGauge clock 1 offset (by norm_num) ≠ clock := by
  intro hEqual
  have hAtRecord := congrArg (fun candidate => candidate.read record) hEqual
  have hOffsetZero : offset = 0 := by
    dsimp [affineGauge] at hAtRecord
    linarith
  exact offsetNonzero hOffsetZero

end ClockReadout

/- ## Compile-time non-coercion receipts

Lean's theorem language is open-world: downstream code can always declare a
new `Coe` instance.  The correct closed-source audit is therefore to make the
dangerous implicit conversions fail during elaboration in this module.  The
tests below cover the ontology shortcuts most likely to be introduced by
accident.  All supported bridges above require an explicit named map.
-/

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe (UniverseClosure Unit (fun state => state)) (RepairOrder Unit)

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe (RepairOrder Unit) (ObserverRecordOrder Unit)

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe (ObserverRecordOrder Unit) ModularParameter

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe
  (WorldlineRealization (ObserverRecordOrder.discrete Unit)
    (fun _ _ : Unit => False))
  (ClockReadout (ObserverRecordOrder.discrete Unit))

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe ModularParameter ProperTime

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe (ClockReadout (ObserverRecordOrder.discrete Unit)) ProperTime

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth Coe ProperTime
  (GlobalTimeFunction Unit (fun _ _ => False))

/-! ### Exhaustive transitive-coercion matrix

The samples above guard the most likely direct Coe mistakes.  This exhaustive
receipt checks all 56 ordered distinct pairs through Lean's transitive
coercion class.  It therefore also catches an implicit conversion assembled
through an intermediate type. -/

private abbrev ClosureTest := UniverseClosure Unit (fun state => state)
private abbrev RepairTest := RepairOrder Unit
private abbrev RecordTest := ObserverRecordOrder Unit
private abbrev ModularTest := ModularParameter
private abbrev WorldlineTest :=
  WorldlineRealization (ObserverRecordOrder.discrete Unit)
    (fun _ _ : Unit => False)
private abbrev ClockTest := ClockReadout (ObserverRecordOrder.discrete Unit)
private abbrev ProperTest := ProperTime
private abbrev GlobalTest := GlobalTimeFunction Unit (fun _ _ => False)

-- Source: universe closure.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClosureTest GlobalTest

-- Source: repair order.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RepairTest GlobalTest

-- Source: observer record order.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT RecordTest GlobalTest

-- Source: modular parameter.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ModularTest GlobalTest

-- Source: worldline realization.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT WorldlineTest GlobalTest

-- Source: clock readout.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest ProperTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ClockTest GlobalTest

-- Source: proper time.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT ProperTest GlobalTest

-- Source: global time function.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest ClosureTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest RepairTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest RecordTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest ModularTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest WorldlineTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest ClockTest

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeHTCT GlobalTest ProperTest

-- A named realization map also remains an explicit record, not a function.

/-- error: failed to synthesize -/
#guard_msgs (substring := true) in
#synth CoeFun (NamedRealizationMap Unit Unit) (fun _ => Unit → Unit)


/-! ## Axiom audit -/

#print axioms canonicalLedgerKinds_pairwise
#print axioms ClockReadout.offsetGauge_ne

end OPH.TimeOrderLedger
