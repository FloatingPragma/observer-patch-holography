# Causal dataflow of the actual native controller

The record-only adaptive selector can communicate between disjoint native
seams. This package authenticates that dependency in actual
[`native_program`](../native_program/interpreter.py) executions. It retains
every raw branch and distinguishes endpoint, controller and storage edges.
The full dataflow of this single-FIFO interpreter is a chain. This is an
implementation fact, not an identification of physical spacetime order.
Its read/write footprint is a sound operational dependency
overapproximation, not a test that every edge carries physical influence.

## Source and domain

The input is an admitted local H66 or H78 program, a classical total-two
basis preparation, and a complete execution within the existing interpreter
limits. The repair law, fair or biased coin weights, C/E or mixture choice,
record-only selector, FIFO policy and raw alphabets are those of the actual
program. No alternative scheduler or geometric placement is introduced.

`execute` calls the existing interpreter and keeps every returned raw
history, including histories with exactly zero probability. On a positive
branch, `replay_trace` calls the interpreter's actual `_outcomes` at every
step and requires the next record to be a unique positive-weight outcome.
The native endpoint function independently checks its load and current
payload. Normalizing this one branch is used only to authenticate its next
conditional transition; all original branch weights remain in the complete
output. The final probability law is never postselected or renormalized.

Zero branches retain their labels and zero weights, but receive no claimed
executed value log. A trace of one arbitrary coherent superposition cannot
be assigned classical endpoint values this way, so this adapter does not
claim a coherent-state provenance theorem. Federated carriers, source pulse
services, intrinsic timeouts and host resource stops are also outside this
first complete-run adapter. The interpreter retains those outcomes; the
adapter rejects rather than silently reinterprets such an execution.

## Authenticated values, full footprints and finite order

The log carries twelve native load registers, the FIFO, and the export
archive. Every positive event retains its before/after values, last-writer
versions, read/write sets and dependency categories:

* A native repair reads and writes its two endpoint registers, even when
  its current is zero. A hold does not access native loads.
* `last_current` and `observer` selectors read the latest FIFO record's
  current. The selected seam is checked against the actual `_outcomes`
  branch. The controller is not allowed to inspect the hidden density.
* Every completed repair or hold calls `_commit`, reads and updates the
  shared FIFO, and reads the archive even when it is still empty. Overflow
  appends the oldest record and stamps a new archive writer.

The selector read is retained even where its two alternatives coincide or
the realized outcome is a hold: this is the interpreter's executed read
support, not a claim of minimal counterfactual influence. Immutable program
parameters and externally specified random choices are source inputs, not
invented earlier events. This package does not localize the auxiliary
random-bit service of [`repair_scheduler`](../repair_scheduler/README.md).

An earlier event is an authenticated parent precisely when it wrote a
register in the current read support, its identifier is the current
last-writer version, and its written value is the value now read. These are
the clauses of `ReadAfterWriteAt` in
[`SemanticEventProvenance.lean`](../../Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean).
All writers initially have no event identifier. Event `j` receives the fresh
identifier `j` only after its reads. Thus every parent identifier is strictly
smaller than its child's. The transitive closure is irreflexive; adjoining
equality gives a partial order, and longest-parent-path height is well
defined. This supplies the concrete execution hypotheses of that finite
order construction. The new Python adapter is independently checked; it
is not presented as a new Lean formalization.

Every event after the first reads the FIFO written by its predecessor.
Consequently the **full** logged order is the total emitted order, with
height `j` at event `j`. Endpoint edges alone can still form a sparse DAG;
that projection omits controller and storage dependencies and is a
different relation.

Read/write footprints give a conservative implementation-dependency
relation. Serialization of an append-only record alone is not evidence of
physical signalling between the recorded events. A physical causal order
could use an observer-visible quotient, but removing such edges requires a
semantics and record-equivalence argument. The adaptive countercontrol below
is stronger than a storage-order observation: it changes an actual remote
outcome probability.

The existing
[`HistoryCausalInvariance.lean`](../../Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean)
proves invariance under swaps whose write sets are disjoint and neither
write set meets the other's certified support. `independent` implements
that predicate on the full footprints. Two disjoint endpoint updates
commute on native loads, as independently checked over all 78 native
states and both coins. Their complete single-FIFO commits fail the full
independence predicate because both write the FIFO. The swap theorem
therefore cannot be used to erase that storage dependency.
This failure of one sufficient swap criterion does not prove physical
influence. A semantics that identifies record-equivalent executions could
remove mere logging order, but its equivalence and observer-visible
causal relation would need a separate proof.

## An exact adaptive countercontrol

Run two C repairs with unbiased native coins and these actual commands:

```
first:  fixed seam (0,1)
second: (2,3) if the latest current is positive, otherwise (2,10)
```

Both possible second seams are disjoint from the first. Compare the binary
preparations with occupied ports `{0,2}` and `{1,2}`. Their remote loads at
ports 2, 3 and 10 agree. In the first preparation, the first current is
positive with probability `1/2`; in that case the second repair moves the
unit from 2 to 3 with probability `1/2`. In the second preparation, the
first current is never positive. Hence the final receiver probability is

```
Pr(port 3 nonempty | occupied {0,2}) = 1/4,
Pr(port 3 nonempty | occupied {1,2}) = 0.
```

Both executable distributions have all 625 raw two-step histories, four
with positive probability, and total probability exactly one. Every positive
second event has no endpoint parent but does have the authenticated first
event as its controller and storage parent. Replacing the adaptive second
selector by the fixed seam `(2,3)` gives probability `1/2` for both inputs,
isolating the role of the controller read.
Unlike the storage chain alone, this change in an actual receiver
probability is a direct witness of influence through the adaptive selector.

This does not violate native endpoint locality. It shows that the complete
controller is an additional communication path. Applying the fixed-schedule
endpoint cone of [`common_source/CAUSALITY.md`](../common_source/CAUSALITY.md)
to this adaptive program would omit that path.

## The spacetime premises still needed

The finite authenticated order is useful without a spacetime interpretation.
The existing exact source Gram has rank three, and the finite Lorentz-carrier
construction adjoins a time axis. It still needs a physical spatial readback,
a clock, an edge-speed bound, and the converse from cone membership to
actual ancestry. The small native diamonds in
[`physical_space_time`](../physical_space_time/causal_placement.py) supply
explicit examples, while arbitrary tapes fail that converse.

The stronger golden-population order-limit theorems use their specified
three-dimensional population, all-neighbour read rule and layer duration.
They are not yet a theorem that this native controller realizes those
premises. A count-volume relation, local refinements, actual observer
communication costs, and identification with the same field propagator's
clock and metric are also required. In particular, the S2/null-ray
equivalence alone does not establish those conditions. No dimension,
Lorentz cone, metric, bulk field, or physical signal speed is assigned by
the new dataflow log.

There is a direct obstruction to identifying the full single-FIFO order
with an open Lorentzian region: every pair of its events is comparable,
whereas any open Lorentzian region of spatial dimension at least one
contains spacelike pairs. An order-reflecting image of this chain cannot
contain such a pair. Nor can such images approximate both members of every
fixed pair in an open flat region at common levels with the same closed
cone relation: approximants to a strictly spacelike pair remain spacelike
when sufficiently close. A physical local-event order would therefore need an explicitly
justified change of interface or treatment of the centralized storage
events. It excludes only that order-preserving and order-reflecting
identification of the exact full chain. Bookkeeping serialization alone
does not exclude a physical Lorentzian limit, an observer-visible causal
quotient, the separately declared golden population, a distributed
observer, or another source construction.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_causal_dataflow
```
