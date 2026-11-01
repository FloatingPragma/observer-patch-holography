# Reusable finite-precision records across captured carriers

This packet implements repeated classical record reads between carriers 0 and
3 of the captured W12 support. Two prepared records each encode one of four
levels, `{-3/2,-1/2,1/2,3/2}`, around baseline 2. The request sequence is
`[0,1,0]`: a record is reread after the same transport bus carries a different
record. All 16 initial payload combinations have complete fixed-point tapes.

The exact protocol uses scalar pair means. The finite implementation rounds
their outputs to a `2^-20` grid with ties to even. Its rounding error is proved
and charged separately. Preparation, accepted schedule, isolation, addresses,
version labels, comparator and additional error limits are supplied.

## Captured support embedding

A global port is `12*carrier + local_port`. The register map is injective:

| Abstract cell | Positive global port | Negative global port | Location |
| --- | ---: | ---: | --- |
| Record 0 | 0 | 9 | Carrier 0 |
| Record 1 | 7 | 4 | Carrier 0 |
| Bus 0 | 1 | 5 | Carrier 0 |
| Bus 1 | 14 | 40 | Carriers 1 and 3 |
| Bus 2 | 23 | 38 | Carriers 1 and 3 |
| Bus 3, receiver | 45 | 39 | Carrier 3 |

The positive path is `1 → 14 → 23 → 45`; the negative path is
`5 → 40 → 38 → 39`. They have three edges each and no shared port. Every mean
is checked against [the captured support](../source_routing/support_w12_l3.json),
including the cross-carrier seams `1–14`, `23–45` and `5–40`. The source and
receiver carriers are glued neighbours via `5–40`. The internal paired cells
are mathematical coordinates spanning different carriers; the
protocol uses no readout or reset between those separated rails. Source and
receiver operations use rails on their respective carriers. This is one
embedded route on the declared capture, not a general support compiler or a
source-selection theorem.
The 12 counted registers are the active subsystem of a 15,360-port support.
Other ports are spectators under the supplied schedule; their preparation,
storage and the physical realization of isolation are outside this accounting.

## Native export, transport and cleanup

For balanced rails `[b+a,b-a]`, two matching means copy a source cell into a
bus cell: both amplitudes become their average. Three such paired transfers
move information along the bus. Every intermediate load is retained. A read
uses eight means: two for export and six for forward transport.

The only reset edge used is the actual source-side bus edge `1–5`. Averaging
that balanced pair clears its amplitude. A forward sweep then transforms
the four bus amplitudes into

```
[a1/2, a1/4+a2/2, a1/8+a2/4+a3/2, a1/8+a2/4+a3/2].
```

This seven-mean cleanup sweep reduces the maximum ideal bus amplitude to at
most `7/8` of its value, while leaving both source records untouched. Eighty sweeps
bound every bus residual by

```
epsilon = (3/2)*(7/8)^80 < 0.000034415.
```

There is no assigned zero or free ancillary refresh. Each read is followed
by all 80 sweeps, including the final cleanup and any equal-input means.
Those identity evaluations count as operations; strict repair progress or
fair scheduling is not inferred.

## Payload and error bounds

For a source amplitude `a`, its kth export leaves amplitude
`a/2^k + e`, with `|e| <= epsilon*(1-2^-k)` in the ideal protocol. At the
receiver the reference amplitude is `a/2^(k+3)`, with ideal error at most
`epsilon`. The independent verifier reconstructs exact linear maps from both
arbitrary source amplitudes to every read; their coefficient bounds cover
all amplitudes in `[-3/2,3/2]`, not only the 16 retained histories.

On the integer grid, nearest/even rounding satisfies
`|2*roundedHalf(z)-z| <= 1`, hence per-mean scalar error at most `1/(2Q)` for
`Q=2^20`. Pair means do not increase maximum-coordinate errors. At read K,
an implementation satisfying the declared one-step error bound therefore has
receiver half-contrast error at most

```
epsilon + E0 + K*(delta + 1/(2Q)) + rho.
```

`E0=2^-18` bounds initial error, `delta=2^-28` is an additional per-operation
bound on every rail, and `rho=2^-16` bounds terminal error per sampled rail.
These are assumed allowances, not measured hardware limits or executed
additional disturbances. The retained integer histories execute rounding.
One rounded mean can change total scalar load by a grid unit. It is an
approximation to the exact mean, not another exact native conservation law.
The general perturbation theorem applies to any implementation satisfying
the one-step bounds; no physical mechanism combining these errors is derived.
`program_readout_error` composes native compilation, per-step perturbations,
two local readout errors and the ideal transfer estimate in one Lean theorem.
The cleanup contraction applies to balanced ideal amplitudes. A common offset
on both rails is preserved by means and must remain in the implementation
error budget; cleanup does not restore an arbitrary raw baseline error.

Four-level decoding is unique when that total is strictly below half the
level spacing `1/2^(k+4)`. All three reads have a positive certified margin;
the minimum exceeds `0.0150217` in raw amplitude units. The receiver uses its
two local values and the declared read count to compare against scaled
thresholds. This comparison supplies no physical amplifier or unit-gain write.
At two fractional bits, a test-time control loses the first record and
produces an ambiguous receiver value. Exact midpoint controls also reject
unique decoding.

## Retained histories and costs

The 16 production histories retain **27,264 means**. Every event stores its
two register indices, both input integers, both consumed writer identities
and its output. Negative writer indices identify prepared registers; a
nonnegative index identifies a preceding mean. Receiver checkpoints contain
only local samples, their writers and the declared record/use identifiers.

Each complete history uses the following resources:

| Quantity | Bound or exact count |
| --- | ---: |
| Scalar registers, including both records and all bus rails | 12 |
| Preparation writes | 12 |
| Means, including all cleanup | 1704 |
| Cross-carrier means | 729 |
| Scalar reads by means, plus receiver samples | 3408 + 6 |
| Scalar writes, including preparation | 3420 |
| Scalar state storage on the declared grid | 276 bits |
| Mean-adder / receiver arithmetic width | 24 / 29 bits |
| Expanded mean schedule, two four-bit register indices per mean | 13,632 bits |
| Global-port map | 168 bits |
| Program counter / two use counters | 11 / 4 bits |
| Read-position, record and use metadata | 42 bits |
| Receiver addresses / request record list | 8 / 3 bits |

The bit counts describe specific scalar and control encodings. They do not
include a physical controller implementation, clock, comparator circuitry,
communications framing, codebook/calibration representation, software runtime
or external audit storage. The retained tapes and exact rational verifier
workspace are external evidence, not memory created by the source dynamics.
No full physical-machine capacity bound is claimed.
The storage and arithmetic widths apply to the retained integer model. The
receiver width bounds signed integer subtraction and threshold scaling, not
the Python producer's rational objects or the abstract real perturbation model.
The verifier derives traffic and encoding counts from the validated tapes,
ports and read metadata; the receipt retains observed arithmetic maxima too.

The ideal reread has a nonzero coefficient from the other record. Approximate
cleanup does not erase that raw dependence. Some rounded states can become
exactly blank; this does not make the ideal cleanup exact. Writer ancestry is
retained through every cleanup operation regardless of its output. Logical
intervention comparisons concern the decoded four-level payload under the
margin, not absence of physical causal influence. The protocol creates no
new versions after preparation and performs no arbitrary record writes.

## Proof, verification and scope

[SourceReusableBus.lean](../../Lean/Geometry/SourceReusableBus.lean) contains
25 theorems: compilation into native means, the finite cleanup formula and
counts, contraction and composition, retention and receiver errors, readout
separation, and integer rounding bounds. The downstream audit checks every
theorem transitively against `propext`, `Classical.choice` and `Quot.sound`,
and rejects `sorryAx` and compiler-trusted Boolean evaluation.

`verify.py` imports no producer. It constructs the schedule by event index,
checks captured edges, independently rounds exact rational means, reconstructs
the full ideal linear maps and scrub matrix, and checks all custody and costs.
Source hashes bind the capture, proofs, contract, toolchain and replay code.
The command-line verifier rejects noncanonical control and receipt bytes.

```
python code/source_reusable_bus/verify.py
python -m pytest -q code/source_reusable_bus
cd Lean
lake build Geometry.SourceReusableBusAxiomAudit
```

Intentional regeneration uses `build.py`, followed by
`verify.py --write-receipt`. Both generated artifacts must match their bytes
on Windows and Linux. Dedicated CI executes the controls on Linux;
Lean CI explicitly includes the downstream audit for dependency-only changes.
Adversarial controls maximize signed error at each of the three reads using
exact backward coefficient maps, perturb idle registers, and verify the
common-mode boundary. These controls exercise the declared abstract bounds;
they are not production histories or measurements of physical noise.

The [contract](CONTRACT.md) delimits this result. Native
existence of this finite route does not select the code, scheduler, source
population or feedback law. Scalar amplitudes attenuate and errors accumulate;
unlimited reuse, immutable scalar archives, versioned commits and general
metric-neighbour routing are not derived. The q=13/q=21 execution and the
physical clock/count and common-world obligations are not discharged.
