# Finite source operator join

This package reproduces and independently checks a finite operator join of
the existing observer 86/88 and 86/247 carriers. It consumes the pinned Lean
source paths and transition tables, retains all 32 rows and 31 actual adjacent
transitions, and checks the counted marginal states and chronological
transports. Sparse matrix-unit controls test the noncommutative hinge,
spectator locality and full-algebra generation without dense 2366-by-2366
matrices. The universal statements and actual source-generator generation
are proved in `Lean/QFT/TripleCarrierOperatorJoin.lean`.

From the repository root:

```bash
python3 code/source_operator_join/join.py --output /tmp/operator_join_packet.json
python3 code/source_operator_join/verify.py /tmp/operator_join_packet.json
python3 -m pytest -q code/source_operator_join/test_join.py
cd Lean
lake env lean QFT/TripleCarrierOperatorJoin.lean
```

`operator_join_packet.json` is a deterministic mathematical evidence packet,
not an observed measurement export. The independent verifier does not import
the producer; immutable parent hashes are pinned inside the verifier. It
rejects changed parents even if a packet refreshes its own hashes. Parent
changes require scientific review and a new version, not automatic pin updates.
The Python verdict verifies this source-bound packet and its finite controls;
it does not run the Lean kernel or prove the universal operator theorems.
Its result explicitly marks both of those limits. Run the separate Lean
command for the general intersection, generation and readout proofs.

The operator embedding adjoins an identity spectator. Its reverse normalized
partial trace is a linear retraction, not a multiplicative homomorphism. The
empirical diagonal state uses unnormalized partial traces for its marginals
and is distinct from the supplied uniform tower states. The shared hinge
is the full 13-by-13 matrix algebra; the two overlapping pair algebras do not
commute. Only the exclusive observer-88 and observer-247 factors commute.

The observer-like structure is the retained finite patches, local label
states, transition operators, field-projector readback, checkpoint records
and common-carrier feedback transport. Tensor-slot assembly and path-read
transpositions are supplied postprocessing, not source-selected physical
operations. This package supplies no physical region or clock, quantum
preparation or outcome record, continuum refinement, universal algebraic
pushout, or identification with the separate 64-mode scalar instrument.
