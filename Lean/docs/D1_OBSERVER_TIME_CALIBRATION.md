# D1 bounded observer histories and clock calibration

The five modules under `Time/` implement the bounded part of V2 issue #691.
They separate what follows from record order, what follows only after a metric
calibration, and what still requires a source-produced physical clock.

## Exact finite results

- `ObserverHistory.lean` turns an explicitly ranked strict history into a
  monotone real readout. The three-record chain is nontrivial; a constant clock
  on a discrete order is the vacuity control.
- `ClockReadout.lean` proves that every strictly increasing regrading preserves
  the bare clock interface. The cubic three-tick clock (0,1,2\mapsto0,1,8)
  is not any positive affine regraduation. Record order therefore supplies
  only order data and an order-compatible scalar readout. Positive affine maps form the explicit calibrated
  comparison algebra, and two ordered pairs determine one interpolation of
  their four supplied readout values uniquely.
- `WorldlineRealization.lean` proves that a supplied affine unit-timelike
  history maps every precedence to a future-timelike displacement, including
  in overlapping charts. A literal three-record, two-chart control inhabits
  the interface.
- `ProperTimeCalibration.lean` proves, on that supplied unit-speed branch, that
  the positive clock increment along the same supplied affine history is
  additive and its square is the chart-invariant Lorentz quadratic interval.
  The result is dimensionless until a physical scale is attached.
- `ClockComparison.lean` proves that one shared event leaves rate and origin
  underdetermined. Two ordered shared events select the unique positive-affine
  interpolation of their four readings; affine consistency at any third
  shared event is equivalent to an exact cross-multiplication equation.
  Distinction of the third event and both readings from both anchors makes
  this a nondegenerate check with no new fit parameter. A held-out test also
  requires a separate predesignation and custody protocol. Two explicit
  algebraic controls use atlases synthesized separately from their clocks.
  The readout `(5,7,9)` passes after fitting `(0,1)` to its first two values,
  while the cubic readout `(0,1,8)` fits the same two anchors and fails at the
  distinct third point. These comparison controls are separate from the
  three-record, two-chart geometric witness.

## Premise and realization boundary

`AffineTimelikeHistory` supplies its event atlas, visibility, event map,
monotone clock, future-unit direction, and affine coordinate law. The theorems
do not derive those fields from A1--A3. Shared-event equality is likewise a
supplied public receipt. No module constructs a source history, refinement
transport, physical instrument, SI unit, preferred foliation, global time
function, modular-time identity, common-atlas physical clock pair, or
predesignation/custody protocol. V2 issue #703 owns source-realized clocks,
physical calibration, and coherence for networks of three or more clocks on
the common tower.

## Verification

From `Lean/`:

```sh
lake build OPHTime
lake build ObserverPatchHolography
```

The printed axiom audits contain only Mathlib's standard classical and quotient
axioms and no project axiom or `sorryAx`.
