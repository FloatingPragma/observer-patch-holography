# PR-04 phase-effect static conformance calculator

Register row PR-04 carries a recorded decision dated 2026-08-18: the exact
phase lift `I/2 - (2*sqrt(3)/3)*i*(Q*P - P*Q)` of
`code/born_context_phase_lift` is a declared phase-sensitive effect under the
row's axiomatize disposition. This directory retains legacy filenames while
holding the deterministic calculator for that effect table, the committed receipt
`PHASE_OPERATION_RECEIPT.v1.json`, and an independent verifier.

The integer pairs are generated expected-frequency numerators: each pair is
the exact Born weight of a declared context effect under the declared matrix
`diag(111/179, 68/179)`, scaled to
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

Boundary: an effect is not a state transition or an instrument. These values
are calculated, never measured, and provide semantic conformance rather than
validation. CP outcome maps and a trace-preserving summed channel remain
PR-64; common source preparation, public outcomes, provenance, and custody
remain PR-65; cross-context additivity remains PR-03.
