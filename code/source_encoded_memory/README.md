# Native encoded copy, reset and reuse toward M1

A pair of scalar loads can retain and copy an unknown **logical sign** using
only the existing pair-mean law, then reset used storage without injecting
fresh zeros. The price is decreasing amplitude, hence decreasing error
tolerance. This is a concrete restricted-code construction. It does not
implement M1's arbitrary scalar archive service.

The [contract](CONTRACT.md) fixes the objective, deliverables and bounded exit.
The [pull-request record](https://github.com/FloatingPragma/observer-patch-holography/pull/915) holds the verification history.
`SourceEncodedMemory.lean` has 25 theorems, transitively checked against only
`propext`, `Classical.choice` and `Quot.sound`. It imports the existing scalar
law `pairAverage` from `ScalarSeamRepair.lean` directly.

## Native protocol and what it preserves

Prepare a cell as `[b+a,b-a]`, with `b >= |a|`; blank is `[b,b]`.
The logical decoder compares the two rails. It distinguishes positive and
negative `a`; `a=0` is blank and carries neither bit. The decoder is a
mathematical observation map, not a derived physical comparator.

For any unknown amplitude `a`, two native means give

```
[b+a, b-a ; b, b] -> [b+a/2, b-a ; b+a/2, b]
                  -> [b+a/2, b-a/2 ; b+a/2, b-a/2].
```

The first mean couples the positive rails, the second the negative rails.
Every output and changed source rail is included. The target preparation is
independent of `a`; neither the schedule nor any auxiliary input knows its
sign. Both cells decode the original sign, although neither retains the
original analog amplitude. No factor-two write or protected export occurs.

One mean between either cell's own rails then gives `[b,b]` exactly:

- Clear the target after its copied state exists: the original sign remains
  available for another read into the same blank target.
- Clear the source: the record moves to the target and the old cell is blank.
  It can be visited again later without fresh preparation.

The single-copy and clear words are proved on the full encoded state, allowing
other cells to hold unrelated amplitudes. During the first mean of a copy the
two touched cells are not balanced individually; the second mean must complete
before the specified cleanup. This is a fixed serial protocol, not an arbitrary
interleaving theorem. Protocol writes after initialization are exclusively
pair means; all their actual consumed writers and values are retained.

After `n` copy-and-cleanup cycles or `n` moving hops, the remaining amplitude
is exactly `a/2^n`. The words have `3n` means. Repeated read cycles use four
scalar registers; a declared route through `k` cells uses `2k`. Revisits are
allowed when consecutive cells differ and the required native edges exist.
Lean's `GoodWalk` checks only that consecutive cells differ: it is not an
adjacency or physical-support certificate. The route theorem starts with one
occupied cell and all other cells blank; support admissibility is separately
checked for each retained word. All cells except the current one are blank
after a moving hop. For repeated reads, the copied target state is
retained in the audit tape **before** its cleanup; no passive observation
archive is claimed to have been created.

The exact raw-state formulas use the same baseline `b` in source and target.
Equal target rails at another baseline still represent a logical blank, but
their common-mode offset survives copy/cleanup. They do not provide the
specified reset to `[b,b]`. The noise bound accounts for such offsets rather
than treating logical blankness as a precision refresh.

These logical cycles do not restore the complete analog state. About the
baseline, its quadratic sum starts at `2a^2` and ends at `2a^2/4^n`; the
difference is the sum of the native losses `(x_u-x_v)^2/2`. The verifier checks
that identity after every operation. Q is not identified with physical energy.

## Finite uncertainty and the sharp boundary

The formal error model permits initial error at most `E0` on every rail and
an arbitrary signed disturbance at most `delta` on every rail after each
mean. This includes idle-register drift. Native averaging is nonexpansive
in the maximum-coordinate norm, so after `m` means the state error is at
most `E0 + m*delta`. No cleanup is assumed to refresh this error budget.

With terminal per-rail observation error at most `rho`, a sufficient condition
for a depth-`n` moving record to decode correctly is

```
E0 + 3*n*delta + rho < |a|/2^n.
```

`robust_route_positive` formally composes the native walk, preparation error,
all disturbances and terminal readout. The negative-sign readout follows
from the same state bound and `robust_negative_read`. The inequality is strict:
at terminal error radius `|a|/2^n`, opposite bits both admit the same observed
pair `[b,b]`. No decoder using only that pair distinguishes them. This sharp
terminal ambiguity does not assert that the conservative accumulated bound
is attained by every hardware noise model.

For the declared example `|a|=3/2`, `E0=delta=2^-20`, `rho=2^-16`, the sufficient
test passes at depth 14 and fails at depth 15. A failed sufficient test
means **uncertified**, not an executed bit error. The margin rows are analytic
controls, separately labelled; they are not noisy execution tapes. Tests also
execute bounded noisy reuse, with errors on untouched rails.

Recovering an analog amplitude by multiplying the final deviation by `2^n`
multiplies its error by `2^n` exactly. This mathematical inverse is not a
native amplifier. Constant register count therefore supplies neither uniform
precision nor indefinitely reliable storage. The retained scalar-value bit
tally grows from 14 initially to a peak of 530 for the four-rail, 64-hop
control. That tally counts canonical numerator/denominator representations;
it excludes signs/serialization framing, addresses, versions, controller,
comparator, arithmetic workspaces and the audit archive. No complete machine
memory or physical bit-capacity bound is asserted.

## Locality and retained executions

Eight exact rational executions retain 626 means and 209 complete copied
states. Initial loads have `b=2` and `a=+3/2` or `-3/2`.

| Experiment | Copies | Native means | Scalar registers |
| --- | ---: | ---: | ---: |
| Empty route | 0 | 0 | 4 |
| One copy, both outputs retained | 1 | 2 | 4 |
| Repeated read and target cleanup, each sign | 16 each | 48 each | 4 each |
| Repeated moving reuse, each sign | 64 each | 192 each | 4 each |
| Declared line route, each sign | 24 each | 72 each | 50 each |

The captured controls map their two cells to W12 level-three ports `(0,1)`
and `(5,9)` in **one carrier**. The verifier checks both rungs and both copy
edges in the pinned capture. Sixty-four back-and-forth hops reuse these four
ports; they are not a spatial distance of 64 or an inter-carrier route. Other
captured ports are outside this subprogram and are untouched, not silently
initialized or counted as a full-population execution.

The 24-hop controls use an explicitly declared 25-cell ladder: two parallel
rail paths and one reset rung per cell. They are not a captured W12 embedding.
This support requirement is substantial: the capture has no matching square
with two local cells in two carriers and both copy edges directly between
those carriers. That finite fact does not rule out longer encodings, routed
gadgets or other layouts. No full support compiler is supplied.

Each native mean costs two scalar reads and two writes. Initialization costs
one write per allocated rail. The audit snapshots and mathematical comparisons
are external verification work; their cost is not included in this native
word count. No hardware readout or clock calibration follows from it.

## Reproduction and audit

From the repository root:

```sh
python code/source_encoded_memory/verify.py
python -m pytest -q code/source_encoded_memory
cd Lean
lake build Geometry.SourceEncodedMemoryAxiomAudit
```

The producer runs a direct state machine. The independent verifier instead
composes rational averaging matrices, applies them to the initial state, and
checks the whole experiment against a separately specified schedule. It checks
every read, writer, native output, intermediate copied state, final rail,
allocation, preparation, bit tally, total load and quadratic loss. No producer
import is used. Strict JSON types and canonical fractions reject type aliases
and ambiguous serializations. Eleven content pins bind the proof, toolchain,
support and execution/checking code.

Intentional regeneration, after reviewing changes to those inputs:

```sh
python code/source_encoded_memory/build.py
python code/source_encoded_memory/verify.py --write-receipt
python -m pytest -q code/source_encoded_memory
```

The package tests replay the controls on Linux and check that
the Lean CI target list names the axiom-audit module; the guard reads
executable shell words, so a commented-out target is rejected. Additional
test-time controls cover three-cell revisits, unrelated records, premature
cleanup, common-mode blank offsets and strict readout margins; they do not add
executions to the retained 626-mean receipt.

## Boundary toward M1

This derives native logical copy and logical reset for a supplied balanced
code and schedule. It gives a possible finite-resolution classical precursor;
it does not replace M1's scalar-valued immutable archives, unit-gain exports,
integer accumulators or versioned commits. Arbitrary protected-memory
invariance is not claimed: averaging the stored cell's own rails erases its
bit, and copying into an occupied opposite cell can erase both bits.

The undischarged obligations include selecting/preparing the code and support,
implementing a physical comparator and address/version channel, compiling the
protocol across actual carriers, deriving regeneration or budgeting finite
precision, supporting arbitrary payloads and intervened histories, and joining
population, causal order, physical clock and count/volume semantics on one
family. Native cleanup consumes the old writers even when its exact encoded
output is constant; their causal ancestry must not be removed. Existing
feedback-order and count-volume results do not transfer by changing a label.

No scientific premise, observation row or common-world join is closed by
this bounded packet.
