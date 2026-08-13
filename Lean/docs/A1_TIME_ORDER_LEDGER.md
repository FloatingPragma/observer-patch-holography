# A1 time-and-order type ledger

`Time/TimeOrderLedger.lean` implements a type boundary, not a chronology
theorem.

## Exact inventory

| Layer | Lean type | What an inhabitant means | What A1 does not infer |
|---|---|---|---|
| Universe endpoint | `UniverseClosure State close` | a supplied fixed point of `close` | a duration, a schedule, or existence of the OPH public-world endpoint |
| Repair execution | `RepairOrder Move` | a finite ordered repair schedule | elapsed time or a physical rate |
| Observer history | `ObserverRecordOrder Record` | a supplied strict order on records | a clock or causal-event realization |
| Modular flow | `ModularParameter` | a real local automorphism parameter | a global chronology or proper time |
| Event realization | `WorldlineRealization records causal` | an explicit order-preserving record-to-event map | existence, timelikeness, continuum convergence, or metric calibration |
| Clock coordinate | `ClockReadout records` | a real record observable strictly monotone on the supplied order | a fixed origin, scale, units, or proper-time identity |
| Proper interval | `ProperTime` | a nonnegative real interval | a metric derivation or clock calibration |
| Optional causal branch | `GlobalTimeFunction Event causal` | a supplied real function strictly increasing on the supplied causal relation | existence on the base theory or an absolute master clock |

The optional branch is represented by
`OptionalGlobalTimeFunction Event causal`. No constructor is selected.

## Coercion policy

No ledger type has a `Coe`, `CoeT`, or `CoeFun` instance to another ledger
type. The module checks all 56 ordered distinct pairs through `CoeHTCT`, so a
transitive coercion assembled through an intermediate ledger type also fails.
Seven direct `Coe` guards make the most dangerous shortcuts visible:

- universe closure to repair order;
- repair order to record order;
- record order to modular parameter;
- worldline realization to clock readout;
- modular parameter to proper time;
- clock readout to proper time;
- proper time to global time function.

Lean is open-world, so no theorem can prevent a downstream file from adding a
new instance. The enforceable in-module contract is the absence of such
instances plus compile-time checks that the conversions fail. Supported
cross-layer translations use `NamedRealizationMap`; it deliberately has no
`CoeFun` instance, so application remains the explicit `.realize` operation.

The module provides only three named packaging operations:

- `WorldlineRealization.asNamedMap` for a worldline already constructed;
- `namedModularRealization` for a modular interpretation supplied as data;
- `namedProperTimeCalibration` for a calibration supplied as data.

None constructs its input map.

## Exact clock-gauge result

`ClockReadout.affineGauge` proves that every positive affine regraduation
preserves strict record monotonicity. If a record exists, theorem
`ClockReadout.offsetGauge_ne` proves that every nonzero offset gives a
different readout on the same record order. Thus record order and monotonicity
do not select a clock origin. This theorem does not say that a physical clock
exists or that all physical calibration freedom is affine; it pins the
specific pre-calibration ambiguity used by the plan.

## Verification

From `Lean/`:

```sh
lake env lean Time/TimeOrderLedger.lean
lake build ObserverPatchHolography
```

The module is imported by `ObserverPatchHolography.lean`, so the ordinary
umbrella build checks it. Its printed axiom audit contains no project axiom or
`sorryAx`. The affine theorem has only Mathlib's standard classical/quotient
axioms.

## Downstream boundary

A1 supplies names and type separation only. The bounded D1 packet in
`D1_OBSERVER_TIME_CALIBRATION.md` adds conditional affine-timelike realization,
dimensionless proper-time algebra, and shared-event comparison, while proving
that bare record order retains arbitrary monotone gauge. Source production,
refinement transport, a physical clock and SI unit remain open downstream.
The stable-causality/global-time branch remains optional.
