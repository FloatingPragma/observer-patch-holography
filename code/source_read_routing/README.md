# Full-family record reads on W12 support

This package routes every declared q=13 and q=21 metric read across captured
glued support seams. A six-operation local transport word preserves source
records and implements the declared recurrence and its semantic order,
conditional on supplied memory, feedback, placement and control rules.
Deriving those rules from the source architecture is a separate obligation.

The [contract](CONTRACT.md) gives the construction, proof boundary and costs.
The [production receipts](../../evidence/source_net_causal_poset/routed_read_law/)
retain input, event and phase commitments. The four small q=3 controls retain
complete compressed tapes in `controls/`.

## Complete replay without a large archive

Production streams are regenerated from pinned inputs and source, then fed
directly into an independent native verifier. Every original decoded phase
and whole-stream SHA-256 must match. The verifier reconstructs all metric
decisions, checks every operation and writer, and compares every completed
logical value against an independent exact recurrence. No production binary
files, cloud credentials or multi-gigabyte raw spool are required.

Python uses the repository's pinned dependencies. Native execution requires
a Linux C++17 compiler; Windows users can run these commands inside WSL.

```sh
python -m pip install -r requirements.txt
python -m pytest -q code/source_read_routing code/source_feedback_transport/test_transport.py code/source_routing/test_routing.py
python code/source_read_routing/support_comparison.py
g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/produce.cpp -o /tmp/produce-routing
g++ -std=c++17 -O3 -Wall -Wextra code/source_read_routing/verify.cpp -o /tmp/verify-routing
python code/source_read_routing/native_negative_controls.py --binary /tmp/verify-routing
python code/source_read_routing/check_reproduction.py --producer /tmp/produce-routing --verifier /tmp/verify-routing
for q in 13 21; do
  for variant in baseline source; do
    python code/source_read_routing/regenerate.py evidence/source_net_causal_poset/routed_read_law/q${q}_${variant}.json --producer /tmp/produce-routing --verifier /tmp/verify-routing
  done
done
cd Lean
lake env lean Geometry/SourceReadRouting.lean
```

Each q13 history has 19,113,548 events; each q21 history has 292,722,053.
The four complete histories total 623,671,202 events. Native register storage
scales with the full history: the producer uses 24-byte cells and the verifier
32-byte cells, in addition to their routing buffers and runtime overhead.
The receipts' peak resident memory values describe the recorded producer
runs, not the combined regeneration pipeline or physical patch memory.

The `Source Read Routing` workflow runs Python controls, fresh q3
production and packing on Linux, and four separate complete
production regeneration/replay jobs. A missing phase, changed event, extra
suffix, failed producer, resource discrepancy or checksum mismatch fails.

## Optional materialization

To retain compressed production tapes outside Git, use execution, accounting
and packing in that order. Package each baseline before its intervention:

```sh
python code/source_read_routing/run.py --q 13 --variant baseline --binary /tmp/produce-routing --output temp/routing
python code/source_read_routing/account.py temp/routing/q13_baseline.json --binary /tmp/produce-routing --work temp/accounting
python code/source_read_routing/pack.py temp/routing/q13_baseline.json temp/routing-packed --work temp/packing
python code/source_read_routing/verify.py temp/routing-packed/q13_baseline.json --binary /tmp/verify-routing
```

Repeat for the other production combinations. Accounting uses GNU time and
requires a byte-identical reexecution before recording host memory use.
Packing needs several gigabytes of temporary disk. The retained phase
manifests describe the original compressed parts; their decoded commitments
are the reference used by regeneration.

## Diagnostic boundary

The captured L3/L4/L5 faces and glued ports agree with the
[support-wiring diagnostic](../../evidence/support_wiring_776/README.md).
Full q13 baseline replay also checks its complete read menu and 10,985
logical values. The q13 lag-three interval is interior; the lag-four interval
is clipped. Routing preserves that distinction and does not turn auxiliary
operations into spacetime-volume events. The finite theorem concerns semantic
read-from order, with physical control ancestry and clock identification
outside its conclusion.
