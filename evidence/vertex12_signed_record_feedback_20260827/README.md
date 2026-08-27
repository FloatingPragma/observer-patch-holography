# Literal signed-record feedback receipt

This archive contains the exact finite receipt for a literal
record-conditioned feedback transaction on a twelve-port carrier system. It
was produced by
[`oph-physics-sim` at `ce17921`](https://github.com/muellerberndt/oph-physics-sim/tree/ce17921eb7504106fef1ba445e1349b5367aa676).

Eight carriers each commit a full twelve-coordinate bounded integer record.
For every carrier-port pair, the diagnostic adds a literal (+1) probe to the
working coordinate, reads the corresponding committed integer, and applies

\[
z_{i,p}\leftarrow z_{i,p}+(b_{i,p}-z_{i,p}).
\]

All 96 transactions restore the working coordinate exactly and leave the
record unchanged. Removing the record read and feedback write leaves the
(+1) displacement in place. Increasing the read record coordinate by one
while holding the other inputs fixed changes the feedback increment from
(-1) to zero, so the later write depends on the literal record value.

The receipt covers 720 (A_5) covariance squares, comprising 60 group
elements and 12 source ports on the reference record. It also checks 96
idempotence cases and all 528 disjoint port-pair commutations across the eight
carriers. The producer-free verifier at the pinned simulator commit rebuilds
the records and events from the parent source data and reports `PASS`; its
source and output are retained here.

The standalone archive verifier checks every byte, the receipt and parent
hashes, each literal probe/read/write transition, the ablation and record
counterfactual, the full covariance row census, the idempotence and
commutation rows, and the independent verifier output:

```bash
python3 verify_archive.py
```

To rerun the separate producer-free implementation in a clean simulator
checkout:

```bash
git checkout ce17921eb7504106fef1ba445e1349b5367aa676
python3 -m oph_fpe.dynamics.verify_vertex12_signed_record_feedback_independent \
  data/repair_closure/vertex12_signed_record_feedback_receipt.json
```

This receipt establishes a bounded internal software observer-like
self-reading component. It has not been integrated with the 81,920-row
protected consensus run. Its normal-form statement covers the standalone
literal record-reset transactions. The parent endpoint-repair confluence,
axiomatic selection of the signed source law, spatial translation, laboratory
realization, physical-sector identification, and physical prediction remain
outside the result.
