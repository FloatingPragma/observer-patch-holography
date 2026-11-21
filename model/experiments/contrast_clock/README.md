# A finite clock for native events

A fixed Hamiltonian executes contrast-memory events and two observer
readings on one twelve-port cell. Its moving clock register carries the
source through the program, increments local counters and retains both
reader blocks. No controller switches the gates during execution.
The construction includes all 78 two-load states and the reader's
disagreement and unavailable outcomes.

The [proof](PROOF.md) gives a finite observation interval, complete
instruments with arbitrary quantum references, and the exact return that
undoes the internal records. The [protocol](protocol.json) fixes the
construction and controls before its confirmatory replay. The
[receipt](results.json) retains every clock position, every nominal prefix
archive law, and the complete terminal outcomes for every source basis
input, with competing source and apparatus controls.

The terminal no-pair probability is 0.00319 for the contrast law and 0.42647
for exchange alone. Reversing the contrast sign also preserves the pair
in this finite test (0.00279 no-pair probability), while changing the joint
reader distribution. Persistence alone does not select the sign of the law.

This is an experimental apparatus extension of
[A-contrast-memory-spherical-v1](../../candidate_laws/contrast_memory_spherical_v1/SPEC.md). The clock couplings, initial pair,
counter registers and gate-controlled rail interactions are supplied
architecture. Their derivation from the contrast-memory law, indefinite
record storage and a physical time normalization are separate problems.

The clock register is a sequencer of the computation. It steps one declared
program of events through one time-independent law, and it is hardware of
the calculation. No observer reads it, and it is no clock that the observers
share. The clock of a port observer is its own counter of the events of its
seams, as in section 8 of the specification, and the counters of different
ports hold different values along the program.

## Replay

From the repository root with `model/requirements.txt` installed:

```sh
python -B -O -m model.experiments.contrast_clock
python -B -m pytest -q -W error model/experiments/contrast_clock
```

`--verify PATH` checks a supplied receipt against fresh execution and its
scientific invariants. `--record PATH` writes to an unused path. Source
pins cover the consumed runtime and fixed protocol. A changed law or
protocol requires distinct evidence; the command does not refresh pins.

The source gate is the actual contrast-memory exponential, including its
local diagonal phase. The symmetric gate program is a declared event
schedule; continuous-time binding bounds do not transfer to it by name.
Its two readings occur in the same run and retain measurement backaction.
The reported clock action is dimensionless and the bounded clock
Hamiltonian has its own energy convention.
