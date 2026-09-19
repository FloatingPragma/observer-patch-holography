# Full-family record reads on W12 support

This package executes every declared q=13 and q=21 metric read through
actual glued support seams. It retains the **M1 read/feedback law as a supplied
premise**. It proves a conditional finite compiler result; it does not derive
that law from A1--A3 or exclude every other admissible routing. This resolves
#777 through its explicit declared-M1 exit: the full routing deliverables are
complete and the derivation obligation is transferred to #779/#740.

The complete tapes and their manifests are in
[`evidence/source_net_causal_poset/routed_read_law/`](../../evidence/source_net_causal_poset/routed_read_law/).
Small controls live in `controls/`. The exact contract, including the transfer
of M1's derivation obligation to #779/#740 through PR-52, is
[`specification.json`](specification.json).

| Quantity, per baseline or source-intervention run | q=13 | q=21 |
| --- | ---: | ---: |
| Source sites | 2,197 | 9,261 |
| W12 host carriers | 5,120 (L4) | 20,480 (L5) |
| Declared rounds | 4 | 5 |
| Logical reads, including self | 1,176,764 | 14,559,225 |
| Executed multicast seam hops | 2,972,512 | 46,288,315 |
| Retained events, including preparation/accumulation/commit | 19,113,548 | 292,722,053 |
| Scalar register reads | 23,178,688 | 353,229,265 |
| Scalar register writes | 22,086,060 | 339,010,368 |
| Protected registers, including every retained relay | 3,001,799 | 46,419,927 |
| Mutable registers, including accumulators | 63,637 | 255,021 |
| Maximum shortest glued-carrier distance of a requested read | 44 | 91 |
| Sum of requested read distances, including all rounds | 12,599,152 | 294,030,780 |

## Construction and proof

`capture_support.py` joins precisely the triangular cells sharing a vertex
and uses the public simulator's geometric port-assignment function. Those
public primitives reproduce every face and glued port pair of the old L3
fixture before L4/L5 are captured. Each larger mesh has twelve degree-five
vertices, all other vertices degree six, and `6*C-30` glued seams. The sixty
unused port slots at pentagonal cells are retained. The unpublished
`oph_exact/support_wiring.py` is not a dependency. Geometric floating-point
assignment is captured as declared finite combinatorics; exact replay does
not treat those coordinates as a physical geometry or an interval proof.

Golden-coordinate Morton order assigns the q^3 source sites injectively to
the first q^3 host cells. This placement is supplied. Every host initially
has the twelve-port signed antipodal split of its prepared golden record;
six protected address scalars are also initialized. Repair can change live
loads while the prepared address remains protected, as in the bounded #778
construction. Unoccupied hosts are relays. Source-payload values are separate
immutable registers, initially `site+1`.

For each source version, an ascending-neighbour BFS tree is pruned to the
union of paths to its complete metric-recipient set. Each used edge performs
export, receiver reset, pair mean, factor-two capture, source reset, receiver
reset. All three feedback operations read only their owner carrier's state.
A carrier may access its own different ports through its protected local
memory; that local access is part of M1 and is charged as an export/capture,
not silently treated as an extra intra-carrier mean. Intermediate protected
versions branch to all needed children. A receiver accumulates each required
version once. Only after every read is complete are new immutable layer
versions committed. No route or read interface depends on payload values.

The exact real hop and sharp additive error bound are inherited from
`Lean/Geometry/SourceFeedbackTransport.lean`. The new
`Lean/Geometry/SourceReadRouting.lean` additionally evaluates the concrete
six-operation word in integer half-units, proving source preservation,
capture, cleanup and its 6-event/7-read/7-write cost. This word discharges the
transport hypothesis in `six_event_compiledHistory_eq`, which proves
full-history value equality for every initial integer state and integer
local-commit intervention. It also proves equality of induced logical reachability
from the local-edge and lifted-read certificates. The native verifier checks
the order hypotheses against every consumed writer. The C++ verifier itself
is reviewed executable code, not a Lean-verified interpreter. This covers the complete
order without enumerating billions of logical pairs. A separate Python
traversal checks all 6,561 logical pairs in each of four q=3 controls.

The full q=13/q=21 `source` runs change the centre's initial payload by +1;
every completed value is checked against an independent exact recurrence.
The q=3 controls additionally change a later commit and the initial scratch
loads. Immutable registers are never overwritten in any replay. The earlier
4,928-event `source_feedback_transport` receipt additionally executes old
version rereads after newer versions and intersecting scratch reuse; its
31-test suite remains part of validation. This is an inherited transport
control, not a newly claimed q=21 old-version experiment.

The closed-register sum invariant gives a limited obstruction: no finite
composition of total-load-preserving operations can implement a reset that
changes that total. The new Lean module checks this finite-word statement.
It explains why this compiler's reset is additional feedback. It does not
exclude tomography, other decoders, changed ancillary state, other laws, or
every operational realization of the metric order.

## Resource and custody model

For C hosts, n sites, K rounds, R logical reads and H multicast hops:

- Preparation uses `13*C+8*n` events: twelve scratch ports and one zero per
  host, one accumulator per site, six protected address scalars per site,
  and one initial payload version per site.
- Transport uses `6*H` events, `7*H` scalar reads and `7*H` scalar writes.
- Layer start/commit uses `2*K*n` events; accumulation uses R events and
  `2*R` scalar reads. Each start reads its own zero, and each commit its own
  accumulator. All writes and read-from edges are counted.
- Protected storage is `C+6*n+(K+1)*n+H` scalar registers; mutable storage is
  `12*C+n`. No relay archive is discarded from that accounting.
- A requested depth-d read traverses `6*d` transport events plus its local
  accumulation; shared multicast prefixes are charged once in H. Total
  serial duration in the declared unit-event schedule equals the total
  event count. It is not the distinguished layer count K or a physical clock.

Each scalar in the retained executions has an exact signed 64-bit half-unit
representation. Independent replay checks mean exactness, arithmetic bounds
and every result; all supported runs fit. The arbitrary-intervention theorem
uses unbounded integers: it does not promise that arbitrary interventions or
inputs fit the native binary. Each explicit
event row occupies 64 bytes, including opcode, owner, input/output register
identifiers, both consumed writer identifiers and the result. Two-output
mean writes share the recorded result. Absent reads have a reserved sentinel;
event IDs are row offsets. Register versions are their checked successive
writer IDs, and immutable cells have only one writer.

Input bytes include the complete read menu, port wiring, host mapping and
golden addresses. `implementation_resources` records actual native vector
allocation requests (including simultaneous old/new register allocations),
auxiliary input/BFS buffers, the 16-byte wide arithmetic temporary, the
64-byte event row, and measured native peak RSS including runtime/allocator
overhead. These are host-software measurements, not a physical patch-memory
certificate. The producer's register structure is 24 bytes; the verifier's
extra origin certificate makes its structure 32 bytes. Codec work is in
65,536-event blocks, with a 4 MiB raw block and at most three differential
dependency levels for the archived controls. Stored manifests and compressed
part byte sizes account for the archive metadata too; compression is not a
reduction in executed events or scalar accesses. The raw-block size is not
the codec's total memory use; the RSS observation covers the native producer,
not the whole Python preparation/codec pipeline.

Routing computation is also charged separately. The current implementation
recomputes each source BFS in each round: `K*n*C` vertex dequeues,
`2*K*n*(6*C-30)` directed-edge examinations, `K*n*(C-1)` discoveries and
tree-mark examinations, and `2*K*n*C` parent/mark initializations. Pruning
traverses H parent links and makes R recipient-path queries. These control
counts are recorded in `routing_control_work`; at q=21 there are
11,377,138,500 directed-edge examinations. No optimal compiler-cost claim is
made. Metric-menu preparation is a separate declared setup calculation;
the independent input verifier checks all `n*n` metric decisions. These
control operations and host storage/compression work are not assigned a
physical clock by the unit-event repair schedule.

Later layers store literal field-by-field uint64 differences from the first
layer's explicit tape; intervention phases store differences from their
corresponding baseline phase. Addition modulo 2^64 reconstructs every field.
Byte-plane compression and fixed 32 MiB file splitting are lossless storage
operations. Neither the decoder nor the verifier regenerates omitted events
from a routing recipe. Whole-stream decoded SHA-256 values bind the same
ordered rows emitted during execution. No external custody/signature or
laboratory attestation is inferred.

Semantic read-from order, serial control order, resource/write hazards and
archive custody remain distinct. Constant resets do not consume the previous
scratch payload, but they remain physical resource operations in the tape.
Routing events are not assigned four-volume by #782. Physical clock units,
header/address transmission, source selection, and quantum copying or
instruments are outside this finite classical implementation.

The two production levels are separate finite regulators. No cross-regulator
state/provenance refinement or uniform physical clock/capacity limit is
claimed. Glued-carrier distance counts inter-carrier means with the supplied
local archive-to-port access; it is not the old port-graph distance when
intra-carrier motion is restricted to means alone.

## Reproduce

Python uses the repository's pinned dependencies. Full native replay requires
a Linux C++17 compiler; on Windows, these runs were executed through WSL
Ubuntu. The small Python controls run on both supported platforms.

```sh
python -m pytest -q code/source_read_routing/test_read_routing.py
g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/verify.cpp -o /tmp/verify-routing
python code/source_read_routing/native_negative_controls.py --binary /tmp/verify-routing
python code/source_read_routing/verify.py evidence/source_net_causal_poset/routed_read_law/q13_baseline.json --binary /tmp/verify-routing
python code/source_read_routing/verify.py evidence/source_net_causal_poset/routed_read_law/q13_source.json --binary /tmp/verify-routing
python code/source_read_routing/verify.py evidence/source_net_causal_poset/routed_read_law/q21_baseline.json --binary /tmp/verify-routing
python code/source_read_routing/verify.py evidence/source_net_causal_poset/routed_read_law/q21_source.json --binary /tmp/verify-routing
cd Lean
lake env lean Geometry/SourceReadRouting.lean
```

To execute anew, compile `produce.cpp` and perform execution, accounting and
packing in that order. Package each baseline before its intervention:

```sh
g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/produce.cpp -o /tmp/produce-routing
python code/source_read_routing/run.py --q 13 --variant baseline --binary /tmp/produce-routing --output temp/routing
python code/source_read_routing/account.py temp/routing/q13_baseline.json --binary /tmp/produce-routing --work temp/accounting
python code/source_read_routing/pack.py temp/routing/q13_baseline.json temp/routing-packed --work temp/packing
python code/source_read_routing/verify.py temp/routing-packed/q13_baseline.json --binary /tmp/verify-routing
```

Repeat for the other q/variant combinations. `account.py` uses GNU time and
requires a byte-identical reexecution before adding host memory accounting.
The large intermediate tapes and first-layer differential workspace require
several gigabytes of temporary disk. Final files are split below GitHub's
per-file size limit. For a small end-to-end reproduction of all four controls:

```sh
python code/source_read_routing/check_reproduction.py --producer /tmp/produce-routing --verifier /tmp/verify-routing
```

The `Source Read Routing` CI workflow executes the small controls and both
inherited transport suites on Linux and Windows. Linux also builds the
producer and reproduces all four small controls through execution, accounting,
packing, both independent verifiers and comparison with the retained hashes.
Its replay jobs run the full native replay for both production levels, including semantic negative
controls. The existing mandatory runner is an input to the invariant-mining
campaign's byte-frozen source projection; its bytes are preserved.

## Objective, deliverables, exit and audit corrections

The selected exit is the declared-law alternative explicitly permitted by
[#777](https://github.com/FloatingPragma/observer-patch-holography/issues/777).
The declaration remains M1, and its source-derivation obligation is transferred
to #779/#740 through PR-52. Closing this routing issue does not discharge that
premise or either receiving issue. No A1--A3 derivation or universal routing
impossibility is asserted.

| Acceptance requirement | Delivered evidence |
| --- | --- |
| Establish a sufficient routing rule on the wired federation and its cost | Six-event local feedback hop; fixed pruned-BFS multicast compiler on captured W12 L4/L5; exact hop proof in `SourceFeedbackTransport.lean`; full-layer induction in `SourceReadRouting.lean`; cost formulas and executed census above. The rule is explicitly supplied as M1. |
| Theorem with the finite part formalized | `six_event_hop_exact` and `six_event_hop_cost` evaluate the concrete local word; `six_event_compiledHistory_eq` instantiates full-history intervention preservation. `exact_induced_order` proves both order inclusions. Native replay checks every consumed writer and the local-edge/read-path hypotheses; the independent small oracle exhausts all logical pairs. |
| Full routed receipts at q=13 and q=21 | `q13_baseline.json`, `q13_source.json`, `q21_baseline.json`, `q21_source.json` and every referenced explicit event segment in `evidence/source_net_causal_poset/routed_read_law/`. The independent verifier rebuilds every metric read decision and compares every completed logical value. |
| Retain every intermediate event and account for resources | Every preparation, export, reset, mean, archive, accumulation and commit remains in the decoded tape. Receipts and replay count reads, writes, registers, distances and controller work; host memory and archive sizes are recorded separately. Serial unit-event duration includes all primitives; no physical clock is inferred. |
| Record M1 and transfer its derivation obligation explicitly | `specification.json` lists the supplied read, memory, feedback, placement and control law; `claims/assumption_dictionary.md` records `M1_full_family_read_feedback_and_routing`; PR-52 assigns derivation/selection of that entire law to the open #779/#740 lanes. |

The earlier audit kept #777 open while the q=13/q=21 rows were only a census
and the full-family refinement and resource accounting were missing. Those
deliverables are now actual executions and checked proofs. The closeout uses
the completed construction together with the issue's M1 transfer clause.
The closed-sum reset lemma remains a narrowly scoped supporting obstruction;
it is not the non-derivability/no-go exit.

Relevant prior corrections reviewed:

- [#518](https://github.com/FloatingPragma/observer-patch-holography/issues/518):
  independently computed witnesses and non-identifiability must not be
  replaced by back-solved identities or promoted into physical producers.
- [#507](https://github.com/FloatingPragma/observer-patch-holography/issues/507):
  tests must actually execute in CI, with optional external tools identified.
- [#778](https://github.com/FloatingPragma/observer-patch-holography/issues/778):
  prepared/protected addresses do not prove population production; a
  class-specific invariant is not an obstruction to every allowed law.
- [#900](https://github.com/FloatingPragma/observer-patch-holography/issues/900):
  live derivation obligations must not disappear onto closed issue lanes.
- [3d5188ad](https://github.com/FloatingPragma/observer-patch-holography/commit/3d5188ada2940e076fcaff9a963b6925d26a8468):
  receipt paths use canonical forward slashes, with cross-platform controls.

Source/control interventions use identical interfaces; stale writers,
unsupported seams, wrong laws, missing resets, old-version substitutions,
resource undercounts, corrupt compression, cyclic dependencies and duplicate
JSON keys have rejecting controls. Every negative mutation must actually
change its input; a no-op mutation is rejected by the test itself.

## Local validation (2026-09-19)

The audit regression run passed 138 tests: 96 routing/inherited
transport tests and 42 premise-register/cross-surface tests. The
native checker passed its positive control and rejected all ten semantic
mutations. All four final packed production histories were fully replayed:
623,671,202 explicit events in total. Each run's complete logical values
matched the independent recurrence. All 64 files in the parent archive
inventory passed their byte-size and SHA-256 checks. The public simulator
checkout reproduced the L3 wiring and both new supports byte for byte.

The new Lean module was kernel checked with the repository's Lean 4.29.1,
`autoImplicit=false` and `relaxedAutoImplicit=false`. The order theorem uses
no axioms; the printed dependencies of the other checked claims are only
`propext` and, for value-history equality, `Quot.sound`. This is a standalone
check of the new module, not a claim of a fresh whole-repository Lean build.

All 143 standard mandatory steps passed across resumed runs before the final
exit-contract wording update. The targeted tests above cover that update;
the final full collection imported 5,207 tests. Execution required restoration of stale
Windows CRLF checkout files to their verified exact HEAD bytes and a local
`python3` alias to the installed Python. The frozen mandatory runner was
preserved; the added routing controls execute in their dedicated CI job.
The broader `--full` heavy steps were not run locally. Hosted CI on the
pre-audit commit `5556050` passed both full production replay jobs, the Linux
and Windows controls, the whole Lean build and its two certificate jobs,
all mandatory shards and collection gates, claim validation and paper preview.
See [AUDIT.md](AUDIT.md) for the detailed
contract review and the distinction between that CI evidence and the audit
follow-up. No issue merge or closure is claimed.
