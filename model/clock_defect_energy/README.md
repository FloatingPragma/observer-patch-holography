# Continuing clocks, localized excitations and repair energy

These source-bound diagnostics address clock, defect and energy joins on the
native C/E carrier. They build on the recorded-history transport of
[PR #1149](https://github.com/FloatingPragma/observer-patch-holography/pull/1149),
pinned at `3d282b3d59763fa1a5171ad9a4fbc09fc0abf35a`. The [first contract](campaign.json) fixes the diagnostic
domain. A separate [theorem-join contract](joins_campaign.json) admits the
counted extension of the compiled clock inventory. Both retain the native
repair law and accept no integrated physical gate.

| Owner | Retained result | Next missing join |
| --- | --- | --- |
| [#1136](https://github.com/FloatingPragma/observer-patch-holography/issues/1136), [clock services](CLOCK_SERVICES.md) | All 66 transported-core words execute and certify at parameter 4. Reference metadata survives FIFO export, with policy states, interpreter storage and source pulses counted. | Realize the supplied preparation, memory, coherent-control and delivery services, including physical error budgets. |
| [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138), [joint clock and field](JOINT_CLOCK.md) | One native incident count advances the transported clock and protected field; the complete joint law distinguishes a control with identical separate marginals. The [source-pulse Ramsey service](RAMSEY_SERVICE.md) executes the controlled density-sector readout and retains complete outputs and continuation. The [clock-count inference](CLOCK_COUNT_INFERENCE.md) combines the source detector and field bit to reduce uncertainty about their shared incident count on C and E, with a matched-marginal control. | Establish the regional physical clock/field dictionary and whole-wave service boundary. |
| [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139), [exchange history](DEFECT_HISTORY.md) | The native configuration complex has exact integral first homology `Z/2`; complete C/E records realize its exchange class. The [sharp-sector packet](../exchange_sector/README.md) tests live relational charges, exterior distance and active recovery. | Obtain an autonomously protected localized sector and coherent exchange response. |
| [#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140), [mode energy](MODE_ENERGY.md) | Native current records give a positive Dirichlet form. The product mode has a separate rest/kinetic comparison. The [reference-balance test](REFERENCE_BALANCE.md) finds fixed-clock drift, a stationary source/core measure and zero gradient cost for its relational charge. | Attach one rest action, moving response and conserved stress to a localized excitation with a calibrated clock. |

In the [completion plan](../../docs/COMPLETION_PROGRAM.md), #1135 owns the
shared domain and primitives, #1136 the continuing services and #1137 the
theorem interfaces.
The [theorem application map](THEOREM_JOINS.md) states which premises are
checked on the model and which are conditional. The
[continuation](CONTINUATION.md), [localization](LOCALIZED.md) and
[activity-energy](ENERGY.md) diagnostics retain their results: the
four-word inventory is incomplete, the adjacent pair has mean first exit 15,
and its excess bare movement is negative. The positive results of the
theorem-join contract use their own explicit service and observation contracts.

## Clockwork correction and reuse

Affine equality of a clock generator and a field generator is **not required**
for them to share operational time. Existing clock comparisons concern event
readings, reference/algebra maps and calibrated parameters. The finite Gibbs
clock identities, native-return `exp(-arccos(K))` reference and conditional
clocked-field constructions retain their positive results.

The new three-charge comparison is a narrower reference-substitution test.
Its unequal spectrum and noncommutation with the wave transfer cannot be
reported as a general clock-field obstruction. The
[clock audit](CLOCK_WAVE.md) maps the relevant existing constructions and
finite theorem interfaces before stating this limited negative result.

## Common domain and evidence

The source is the binary total-two shell of the twelve-port carrier, dimension
66, with explicit restriction from the dimension-78 complete instrument where
needed. C and E have the same population transition here and different raw
records and coherence behavior. Both controls are retained.

The transported faithful reference extension is not a posterior corner or
the stationary source/FIFO reference. Its core label is classical metadata;
transporting it does not implement port routing or establish a physical
source charge. Attempt count, actual repair count, modular parameter, wave
parameter and any calibrated physical clock are named separately.

The retained receipts bind their implementations, independent checks and their
named contracts. No empirical comparison is added. Preparation and coherent
services retain their parent domains and failure budgets. The expanded
clock inventory counts every added word and its source execution.
The reference-balance packet distinguishes fixed-reference drift from a
stationary co-moving source/core measure with counted classical metadata.
The finite first-exit calculation is an offline path diagnostic, not an added
measurement of the source.

## Replay

From the repository root with the normal Python requirements installed:

```sh
python -m pytest -q -W error model/clock_defect_energy
python -m model.clock_defect_energy.continuation --check
python -m model.clock_defect_energy.clock_wave --check
python -m model.clock_defect_energy.localized --check
python -m model.clock_defect_energy.energy --check
python -m model.clock_defect_energy.clock_services --check
python -m model.clock_defect_energy.joint_clock --check
python -m model.clock_defect_energy.defect_history --check
python -m model.clock_defect_energy.mode_energy --check
python -m model.clock_defect_energy.reference_balance --check
python -m model.clock_defect_energy.ramsey_service --check
python -m model.clock_defect_energy.clock_count_inference --check
```

`--write` regenerates the named packet's receipt after a reviewed source
change. The all-core source-word replay takes several minutes. The frozen
parent receipts and campaign are never rewritten.
