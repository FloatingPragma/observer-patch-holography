# The native phase interface from existing source pulses

The balanced mixer used by the [native interface](../native_interface/README.md)
can be implemented by **34 pulses of the existing source compiler**, using
three local seam qubits and leaving six spectator qubits untouched. The
balanced gate is implemented through the compiler's coherent proper-code
services. The balanced
angle is a chosen experimental operation; this does not select it from
A1--A3 or produce coherence using committed native repairs alone.

One additional source `Z` pulse gives the imaginary phase quadrature. The
resulting native coin/current read distinguishes two conjugate phase mixtures
with ideal plus probabilities **one and zero**; the real-only apparatus
gives one half for each. These are finite conditional model predictions.
There are no laboratory observations in this packet.

## Why the same source suffices

The existing [nine-bit native encoding](../observer/source_discrimination/circuit.py)
writes a total-two carrier configuration as

```text
address = seam_state + 8 * spectator_label.
```

For the first seam, all ten native orientation pairs have ordered seam
addresses `(2,1)` and identical spectator labels. They are the configurations
`|uw>` and `|vw>`. Thus the requested mixer is the same local eight-dimensional
operation for every spectator:

```text
B|2> = (|2>+|1>)/sqrt(2),
B|1> = (|2>-|1>)/sqrt(2),
B|j> = |j>                       for j=0,3,4,5,6,7.
```

The full target is `I64 tensor B`, including all 512 source addresses. It
preserves the 78-state native code and its 66-state binary sector exactly.
It also mixes the corresponding pairs in the padding sector; its extension
there is explicitly **not** asserted to be identity. No spectator routing,
extra carrier, state-dependent basis change or arbitrary numerical gate is
inserted into execution.

Let `C=CNOT(0,1)`, `R=exp(-i*pi*Y0/8)` and
`D=exp(-i*pi*|3><3|)`. The exact local identity is

```text
B = C R D R* C.
```

`C` sends `(2,1)` to `(2,3)`. On that two-dimensional corner `D=Z0` and
`R Z0 R*=(X0+Z0)/sqrt(2)`; elsewhere the two rotations cancel. The diagonal
reflection is its complete Walsh expansion

```text
|3><3| = (1/8) sum_mask (-1)^popcount(3 & mask) Z_mask.
```

The scalar term is retained. Cyclic Gray walks share the commuting parity
computations and restore each chart. This needs eight phase pulses and six
source CNOTs. Each phase-exact CNOT uses three existing primitive pulses.
Including the two rotations and two outer CNOTs gives **34 pulses**, comprising
26 single-code and eight pair-code pulses, with **68 load/export service
events**. The unchanged generic orthogonal compiler emits 1,046 pulses for
the same operation and is used as an independent comparison; that count is
asserted in `test_construction.py` and is not retained in a receipt.

Every emitted pulse calls the existing full six-level response evaluator.
Its retained twelve-port control coefficients are checked in `Q(sqrt(5))`
against the source response and realization receipts. `H` controls follow
from the available `X` and `Z` controls. No target matrix is used by
the pulse executor. The target above is only the independent oracle.

## Original binary arrays, references and padding

`local_matrix()` returns the raw eight-by-eight source execution.
`full_matrix()` lifts it to all 512 addresses using the fixed spectator
wiring. No polar correction, normalization, amplitude cutoff or code
projection is performed. All off-target binary entries remain present.
The existing compiler contracts each processor pulse to its analytically
invariant two- or four-dimensional proper code; this packet's certificate concerns
the returned encoded array, not a numerical simulation of all processor
load/export environments.

`certify()` converts every original binary component to an exact rational
number, compares it with the complete independent target in
`Q(sqrt(2))+i Q(sqrt(2))`, and encloses the Frobenius norm using directed
rational square roots. With `epsilon >= ||A-B||op`,

```text
delta = epsilon*(2+epsilon)/2
g     = (1+epsilon)^2
```

bound half the diamond distance of the single-map outputs and the actual
binary channel norm respectively. These follow from
`AXA* - BXB* = (A-B)XA* + BX(A*-B*)` and Schatten norm
submultiplicativity. They apply to arbitrary operators and every untouched
finite reference. Tensoring the map with spectator identity does not increase
its operator norm, so the same bound covers all 512 output coordinates.

For two uses separated by exact CPTP processing, the complete instrument
error is at most `delta*(1+g)`, including binary norm amplification. Adaptive
use on successful branches and identity on failed branches obeys the same
bound. All three retained operations (balanced mixer, one-pulse phase, and
35-pulse imaginary-quadrature analyzer) have portable reported single- and
two-use bounds at most **`1e-12`**. These are numerical implementation bounds,
not measured physical gate fidelities.

The API deliberately provides no compressed “actual 66D source mixer.” A
66D calculation may use `ideal_native_matrix()` as an explicitly ideal
surrogate and add the full-source replacement error. To do so its other
instruments need a stated CPTP extension on the code complement, all padding
and failure records must remain represented in the comparison, and any
subsequent floating density arithmetic needs its own error budget. The
source certificate does not silently grant those conditions. An untouched
reference, all 78 carrier columns, all padding columns, wrong phases,
wrong orientation, off-support defects and hostile array subclasses are
covered by the checks.

## The missing imaginary quadrature costs one pulse

Write `Z_orientation=P_u-P_v`, `X_orientation=sum_w(|uw><vw|+|vw><uw|)` and

```text
Y_orientation = sum_w (-i|uw><vw|+i|vw><uw|).
```

The native coin/current pointer has effects
`D_±=(I±Z_orientation)/2`. The mixer gives
`B*D_±B=(I±X_orientation)/2`. On each active pair, the existing source
`Z` pulse on seam bit zero, with angle `pi/4`, restricts to
`V=exp(i*pi*Z_orientation/4)`. Its other diagonal phases do not affect
the inactive pointer, whose effects are `I/2`. Applying `V` before `B`
therefore gives exactly

```text
(B V)* D_± (B V) = (I±Y_orientation)/2.
```

The complete target unitaries and these effects are checked entry by entry
in the exact number field, including the inactive sector. The two uniform
mixtures of `(|uw>±i|vw>)/sqrt(2)` have plus probabilities one and zero.
Preparation uses the same services: after a successful native orientation
herald, apply `B` and then the inverse source phase pulse. Both herald signs
remain available with their recorded interpretation.

For states known to occupy the active orientation corner, the three reads
determine the reduced orientation qubit after tracing out its spectator.
They do not determine spectator coherences, arbitrary 66-dimensional states
or an unknown active-sector weight. There is no claim of full-carrier
tomography.

## Results used, and what is supplied

The [coherent source-code construction](../../docs/research/COHERENT_SOURCE_CLOCKS.md)
supplies the two ingredients: exact realization of single-code
rotations and phase-exact pair entanglers from the declared twelve-port
response, and identity code transport conditional on preservation of two
complementary read tables. One-basis agreement permits complete
dephasing. This packet adds the concrete three-qubit word and its native
encoding map; it does not derive coherent gluing from classical agreement.

[`SourcePhaseLiftBridge.lean`](../../Lean/QFT/SourcePhaseLiftBridge.lean)
proves `rhoY_bornWeight_equal_on_real_symmetric` and the corresponding
real-effect closure obstruction. Its positive phase-lift and qubit
tomography results require access to the complex effect. Here the explicit
source pulse supplies that missing operation on the reduced orientation
qubit under the existing services. Algebra membership by itself would not
have supplied it. The [native partial-injection argument](../native_interface/README.md#the-exact-reachability-boundary)
and `dephase_kraus` in
[`SourceRepairInstrument.lean`](../../Lean/Geometry/SourceRepairInstrument.lean)
retain their separate committed-repair scope. This is a direct finite-matrix
attachment, with no claim of a Lean proof.

The carrier encoding, chosen balanced measurement, source control response,
coherent two-/four-dimensional proper-code services, apparatus preparation,
storage and readout are premises. Their physical errors and timings are
unbounded here. The result reduces an additional gate-availability premise
within that common construction; it does not select a physical action,
repair candidate, universal scheduler or continuum field theory.

From the repository root:

```sh
python -m pytest -q -W error model/compiled_native_interface
python -m model.compiled_native_interface --check model/compiled_native_interface/results.json
```

The receipt is freshly recomputed; strict comparison preserves exact source
hashes, scalar types, outcomes and outward guarantees. `--write PATH` refuses
to overwrite an existing file. The package manifest binds code, tests,
documentation and retained results.
