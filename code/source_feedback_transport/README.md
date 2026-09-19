# Reusable local record transport

Bounded software patches exchange classical immutable versions over the
captured W12 seam support. Each patch has mutable scratch ports, protected
records, local readback and feedback. Every operation has public read/write
provenance. The destination can reread an old version after a newer version
has used the same channel.

The supplied local grammar adds three operations to canonical pair averaging:
export a local protected record into a local scratch port; reset that port
from its local protected zero; and archive twice a local port readout. These
are explicit classical feedback operations, already of the type used by the
software observer instruments. They are not selected by the pair-mean law,
and they do not implement quantum copying or a quantum instrument.

One hop exports the payload, resets the receiver, averages the adjacent pair,
captures twice the receiver value, and resets both scratch ports. It uses six
events, seven scalar reads and seven scalar writes. Exact arithmetic restores
both scratch ports to zero and preserves the archived version. A fixed route
uses one hop per seam. Only immutable local records are exported; no decoder
reads a remote baseline, global state or generator reference.

```sh
python3 code/source_feedback_transport/build_transport.py
python3 code/source_feedback_transport/verify_transport.py
python3 -m pytest -q code/source_feedback_transport/test_transport.py
cd Lean
lake env lean Geometry/SourceFeedbackTransport.lean
```

The 16 complete episodes use path depths 1, 4, 12 and 24, with repeated and
intersecting routes. They contain 740 hops and 4,928 retained events. Source,
branch and initial-scratch interventions preserve the fixed implementation.
A branch intervention changes its descendants while leaving later reads of
an older version unchanged. The verifier imports no producer. It checks the
complete support, exact primitive laws, locality, immutable versions, every
consumed writer and the full induced order on logical commits. A coherently
relabelled different route is rejected even when all its local operations
are valid.

The projection concerns semantic data read-from order. Constant resets consume
only a local zero record, so overwritten payloads do not become data parents.
Every reset and every intermediate mean remains in the event bundle and cost
count. Physical resource ordering, write-after-write constraints and hash-chain
order are separate; no physical causal equivalence or clock follows.
Version identity is fixed by the compiled schedule and authenticated writer
chain. The scalar wire does not separately implement a physical address or
version-header channel.

For per-hop export, zero-reset, receiver-mean, readback and decoder errors
bounded by E, Z, M, R and D, respectively, a depth-d route has error at most
E0 + d*(E+Z+2*M+2*R+D), where E0 is its initial archive error. This sharp
additive bound is conditional on the stated noise limits; no hardware noise
limits are measured. Export error includes retained-archive drift/read/write
error. Cleanup errors do not feed the next payload because fresh export and
reset precede every mean. The exact receipt itself executes no noisy gates.

At depth 24, one episode has 108 hops, 706 total events, 141 protected scalar
registers and 25 mutable scratch registers. All archived relay records are
retained in that count. Written-scalar bit counts exclude identifiers,
versions, hashes, event metadata, decoder temporaries and machine overhead;
they do not certify physical memory capacity. Prepared address records remain
separate from the live port loads changed by feedback.

The q=13 and q=21 rows are a checked census from a pinned declared-family
receipt and a unicast compiler cost calculation. They are not new full-family
executions. Full routing at those levels, source selection of copy/reset and
read rules, physical costs, and quantum preparation/readout remain separate.

The subsequent [full-family compiler](../source_read_routing/README.md) now
retains complete q=13/q=21 routed executions and their interventions under
an explicitly retained M1 law. This earlier receipt and its census remain
unchanged; its small episodes are not relabelled as those new executions.
