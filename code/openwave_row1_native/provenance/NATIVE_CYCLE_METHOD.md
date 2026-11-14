# A native cycle meter for the spatial charge class

The primary source is the unchanged one-load twelve-port native process.
Its initial occupation G is retained, and its current tape starts at zero.
After each actual seam event, the source and tape obey

\[
n'=n-jB_e,\qquad F'=F+j\mathbf e_e,\qquad BF=G-n.
\]

The spatial adapter identifies G−n with the boundary degree of its dual
cell circle map. The meter below reads that degree through another actual
native carrier. It uses committed coin/current records, finite memory and
a current-adaptive seam schedule. It needs no coherent code-transfer,
analog phase pulse or ideal projective detector.

## Finite memory on the primary source

Each port counter has six states. Initialize all counters to zero with the
actual zero-tape source preparation, then update only the two addressed
counters by

\[
c_a'=c_a+j\pmod6,\qquad c_b'=c_b-j\pmod6.
\]

For every regional cut A, the actual relative charge
`sum_A(G−n)` is one of −1, 0, 1 because both G and n contain exactly one
unit. Its class is `sum_A c mod6`. The map `0→0, 1→1, 5→−1` therefore
recovers that regional charge exactly. No unbounded history lift is needed
for this specified observable and preparation domain.

The source state retains G. Across the twelve preparations there are 144
pairs `(G,n)` and 133 distinct charge images. In particular, charge zero
does not identify the background port and cannot replace it as the full
native source state. The executable checks all 8,640 native coin events
on those source pairs and all 589,824 regional cuts. A charge-seven control
aliases class one and lies outside this domain. The finite counter theorem
does not extend silently to arbitrary integer loads.

## Reset from an unknown native pointer position

The pointer is another one-load twelve-port carrier. Its initial port is
arbitrary and is not available to the controller. Fix R retries per requested
edge. Scan all thirty carrier seams, taking R native attempts at each seam
until a nonzero current appears. The first nonzero current identifies the
endpoint reached. A zero current does not disclose the hidden load.

If all discovery currents are zero, the pointer never moved. Every port
has degree five, so exactly 5R of the discovery attempts touched the load.
Each such attempt holds it with probability one half. Discovery timeout
therefore has probability exactly `2^(-5R)`, for every initial position.
After discovery, route the known position to port zero along a shortest
carrier path. Its length is at most three. Each requested move succeeds
within R attempts with probability `s=1−2^(-R)`. All failed attempts are
retained no-move records; exhausting the retry budget is a timeout.

For initial port x, order its five incident seams as they occur in the
fixed scan. Let y_k be the other endpoint of the k-th such seam, starting
at k=0, and let d(y_k,0) be the routing length. The exact reset probability is

\[
p_x=\sum_{k=0}^{4}2^{-Rk}s^{1+d(y_k,0)}.
\]

The receipt checks this formula against a separate calculation using both
literal native outcomes of every attempted seam. A uniform timeout bound
is `2^(-5R)+3*2^(-R)`. The reset costs at most 33R attempts.

Each conclusive fine reset word has one initial source column and ends at
port zero. Its native Kraus map is a scalar multiple of `|0><x|`.
Consequently the complete success/timeout construction also makes sense
with arbitrary untouched references. Conditioning on successful reset may
change a correlated reference marginal. The complete instrument retains
the timeout branch and makes no independence claim from that conditioning.

## A winding read through literal native moves

The fixed carrier cycle is

`0 → 1 → 8 → 6 → 2 → 10 → 0`.

All six edges are actual carrier seams. Counter class one commands one
positive traversal, class five one negative traversal, and class zero no
traversal. The controller requests the next seam and waits for its actual
moving current, retrying up to R times. It advances its cursor only after
that current. The retained signed step count divided by six is the pointer
winding. Every successful command is a closed pointer path with the exact
decoded spatial degree. Timeout leaves the pointer at its actual partial
position and supplies no guessed winding.

For the bounded primary source, at most six moves are requested. The
complete success probability is `p_x*s^6` for either nonzero class and
`p_x` for zero. Thus all native no-moves, discovery failures, route failures
and measurement failures are included. The attempt cap is 39R. At R=16 it
is 624 attempts. The receipt gives the exact rational worst timeout and
complete successful and timed-out replay tapes. Every fine tape has
probability `2^(-number of executed attempts)`; the executable controller
defines the entire finite outcome tree.

Unused counter labels 2, 3 and 4 produce explicit admission-refusal records.
They are neither interpreted as a primary-source charge nor silently dropped.
The longer unreduced contour mode is checked against all 132 existing
zero-tape source histories and all twelve contour vertices and pointer
initial positions, covering 19,008 cases. Its bound is separate from the
six-move compressed primary meter.

The pointer is an apparatus for the spatial class supplied by the actual
current records. Its cycle is not an independently evolved field and is
not used as a substitute proof of the spatial cohomology theorem. This
construction establishes a native operational readout of that theorem's
observable. It supplies no electromagnetic strength, particle mass or
automatic platform status.

Run from the workspace root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge/native_cycle_meter.py
```

`NATIVE_CYCLE_METER.json` contains source pins, exact probabilities, full
finite-controller semantics, reset and meter tapes, and the bounded-domain
counter checks.
