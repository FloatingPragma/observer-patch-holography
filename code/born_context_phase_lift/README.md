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

`SOURCE_PHASE_SELECTION_PACKET.v1.json` is the committed derived packet for
the declared representation, record projector, commutator order, and
complexification adapter. Its replay yields exactly the unordered pair of
Pauli-Y projectors. The stable label order currently puts `pair:00-02` first
and that pair is `+Y`, but this sign is not source-forced: reversing the event
list does not select `-Y`, whereas reversing the commutator order for a fixed
pair exchanges the two transpose effects.

`SOURCE_PHASE_SELECTION_SEMANTICS_RECEIPT.v1.json` records the exact replay
and mutation checks. These two files establish artifact reproducibility only;
they do not authenticate a producer, run, public outcome, or physical phase
orientation.

Run the hermetic gates from the repository root:

```bash
python3 code/born_context_phase_lift/verify_source_phase_lift.py
python3 -m pytest -q code/born_context_phase_lift/test_source_phase_lift.py
python3 tools/test_source_phase_selection.py \
  --payload code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json \
  --module-dir code/born_context_phase_lift
```

Use `--verify-source-hashes` only when the sibling `oph-physics-sim` checkout
is present and the full cross-repository custody replay is desired.
