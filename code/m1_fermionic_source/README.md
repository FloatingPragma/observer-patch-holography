# Local fermionic source: preparation, clocks, reads and interaction

The construction uses observer-like self-reading patches with private local
state, boundary ports, retained readback, feedback and repair operations.
The public receipt binds their preparation, clocks and read instruments.

This package closes the vacuum/one-particle restriction of
[`m1_source_realization`](../m1_source_realization/README.md). It implements
the full parity-preserving fermionic process in that same coherent-code
source, including actual many-particle signs, deterministic code preparation,
measurements, a density interaction and charged physical time. The prior
Bravyi--Kitaev superfast encoding is explicitly credited; its source compiler
and complete resource/noise account are the result here.

Read the [objective, deliverables and exit](CONTRACT.md) and the
[all-size proof](../../docs/research/FERMIONIC_SOURCE_CLOCKS.md).

| Deliverable | Result and check |
| --- | --- |
| Full even CAR algebra | Both parity blocks, signed loop constraints, complete isometries and all-sector gates checked against independent occupation signs and exterior minors. |
| Deterministic preparation | Every small-graph measurement branch retained; local loop basis and causal prefix decoder checked for all syndromes by an independent GF(2) inverse. Preparation wall time scales with diameter. |
| Actual native operations | M6 one-code pulses and its native pair-code CZ implement signed Pauli rotations, both QND outcomes and a quartic density phase. Unknown entangled inputs survive the explicit code-growth instrument. |
| Massive clock | Reflecting two-bank compiler, fermionic flights, exact blank restoration, degree bound, one processor per owner and a conflict-free positive-duration schedule. |
| Reads and interaction | Tree-compiled complex CAR number reads; bounded Fock click effect; non-Slater two-particle output with stabilizer entanglement witness 2 sqrt(2) against separable bound sqrt(2). |
| Valid noise budget | Physical edge dephasing damages even the encoded vacuum; the bound is a full-interface bound; an occupation-only estimate is invalid. A finite analytic clock has positive sufficient noise rate and retained error outcomes. |

## Reproduce

From the repository root with `PYTHONPATH=code`:

```text
python -m m1_fermionic_source.build
python -m m1_fermionic_source.verify
python -m pytest -q code/m1_fermionic_source
```

The build replays its candidate before writing a receipt. The verifier does
not import this package's producer modules. It reconstructs small CAR
operators from occupation signs; tensor-Pauli matrices from I,X,Z; Fock lifts
from determinants; spatial walks from a direct reflecting permutation; loop
corrections from binary Gaussian elimination; and timing/noise budgets with
70-digit arithmetic with one-sided checks on the finite clock bounds.
Weighted detector reads include the native helper reset, conditional pulse
and all four retained histories, with their quantum outputs checked rather
than only the aggregate click probability. It also replays the parent's
source certificate.

The compact receipt contains small matrices, native gate programs and hashes
of exactly regenerated sparse graph/decoder data. No external raw archive is
needed. The verifier rejects missing parity blocks, omitted negative outcomes,
bad loop signs, ordinary swaps, changed complex reads, decoder corruption,
false vacuum-noise claims, zero time/rate budgets, slightly wrong-sided
numerical bounds, incomplete weighted instruments, altered claim/source pins,
extra or duplicate JSON fields and nonfinite numbers. Hostile physics tests
bypass custody to ensure that hashes are not the only defence. CLI rejection
also runs under `python -O` with producer imports disabled.

## Scope and practical limits

The general theorems are analytic. Executions cover five small encoding
graphs, both parity sectors, spatial sides 1--3 and local preparation sides
1--4. Complete spatial matrices and blank-bank columns are replayed, not only
a fitted dispersion curve or a selected successful trajectory. The enormous
q=2^41 clock is a **finite analytic witness**, not an executed many-body
simulation. Its conservative wall-time factor is about 2.05 x 10^7; its
global accounting-error bound demands a sufficient noise rate about
1.13 x 10^-102 in the declared units. These are consistency bounds, not
practical hardware specifications or measured capabilities.

The all-size schedule also yields a controlled refinement family at fixed
physical volume and observation time. With the declared drive/event scaling,
choosing a positive noise rate lambda_q=lambda_0 q^-8 makes the sufficient
state error O(q^-5) and additional accounting error O(q^-1). This is a
construction with a specified noise capability, not its physical selection
or an optimal noise threshold.

CCG and the declared cofinal drive/event capability retain their existing
status. There is no new transport axiom and no freely supplied fermionic
gate set, pure A3 vacuum, instantaneous distributed report or hidden
postselection. Central classical records are assumed reliable in this stated
quantum dephasing model and their operations are counted. As in the parent,
the fixed classical program is distributed in finite charged prehistory
before quantum preparation; the displayed exposure bound starts with that
preparation. Quantum hardware is O(q^3), while this conservative retained
ledger uses O(q^4 log q) central slots at fixed physical volume and horizon.
The construction
does not select fermionic statistics, empirical masses/couplings, a physical
energy standard or a unique microscopic source. The free clock theorem is
not extended to nonzero density coupling. Odd CAR fields are not identified
with commuting local qubit observables.

The three claim rows are `OPH-SOURCE-LOCAL-FERMIONIC-REALIZATION`,
`OPH-SOURCE-ENCODED-FERMION-CLOCK` and `OPH-SOURCE-FERMIONIC-VACUUM-NOISE`.
Their mathematical-model status is retained in the registry and necessity
inventory. Physical-selection premises are not marked discharged.
