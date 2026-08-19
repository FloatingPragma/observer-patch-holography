# PR-04 phase-operation receipt producer

Register row PR-04 carries a recorded decision dated 2026-08-18: the exact
phase lift `I/2 - (2*sqrt(3)/3)*i*(Q*P - P*Q)` of
`code/born_context_phase_lift` is a declared architecture operation under the
row's axiomatize disposition. This directory holds the deterministic producer
for that decision's outcome table, the committed receipt
`PHASE_OPERATION_RECEIPT.v1.json`, and an independent verifier.

The counts follow an explicitly exhaustive deterministic semantics, one of
the two validation routes the committed boundary text names as acceptable:
each context count pair is the exact Born weight of the context effect under
the committed record-diagonal run state `diag(111/179, 68/179)`, scaled to
the least positive integer multiple of the committed run mass 179 that makes
both counts integers. All arithmetic is exact over `Q(sqrt(3), i)`; no
sampling, randomness, or floating point enters. The Lean inhabitant of these
literals is `Lean/EventAlgebra/OperationalPhaseAttainment.lean`.

Run from the repository root:

```bash
python3 code/phase_operation_producer/produce_phase_operation_receipt.py --check
python3 code/phase_operation_producer/verify_phase_operation_receipt.py
python3 -m pytest -q code/phase_operation_producer/test_phase_operation_receipt.py
```

Boundary: the operation is declared, never source-produced, and the counts
are produced, never measured. No laboratory or emergent-instrument claim
follows; the physical-attachment premises of the register stay open.
