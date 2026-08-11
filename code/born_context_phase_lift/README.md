# B13 algebraic phase-lift audit

`BORN_CONTEXT_WEB_PAYLOAD.v1.json` is a byte-identical, hash-pinned copy of
the simulator payload used by the exact verifier. The default verification
path is hermetic. An optional custody test rehashes the larger sibling
simulator inputs when that checkout is available.

The payload is retained byte-for-byte, including its historical prose. Its
phrase that the web “earns noncommuting effect contexts” is superseded by the
audited boundary: the run realizes gauge labels and diagonal counts; the
two-dimensional representation and projector orbit are declared algebraic
adapters. The phase lift
`I/2 - (2*sqrt(3)/3)*i*(Q*P-P*Q)` is an exact operator-algebra target, not a
source-produced instrument, outcome, or Born-law validation receipt.

Run the hermetic gates from the repository root:

```bash
python3 code/born_context_phase_lift/verify_source_phase_lift.py
python3 -m pytest -q code/born_context_phase_lift/test_source_phase_lift.py
```

Use `--verify-source-hashes` only when the sibling `oph-physics-sim` checkout
is present and the full cross-repository custody replay is desired.
