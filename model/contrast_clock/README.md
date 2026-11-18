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

This is an experimental apparatus extension of
[polynomial-sv](../polynomial_sv/SPEC.md). The clock couplings, initial pair,
counter registers and gate-controlled rail interactions are supplied
architecture. Their derivation from the contrast-memory law, indefinite
record storage and a physical time normalization are separate problems.

## Replay

From the repository root with `model/requirements.txt` installed:

```sh
python -B -O -m model.contrast_clock
python -B -m pytest -q -W error model/contrast_clock
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
Hamiltonian has its own energy convention. No eigenvalue is identified
with a particle mass, and no OpenWave physical row is awarded.
