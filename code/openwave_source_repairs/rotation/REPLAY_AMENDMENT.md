# Canonical replay correction

The first clean-copy replay found an artifact-checking defect. The original
entry point compared a parsed JSON receipt (lists) against the in-memory
result (which includes tuples). Those representations are not equal as Python
objects even when their canonical serialized bytes are identical. The failure
is retained in `ATTEMPT_02_FAILURE.json`. No scientific assertion is waived.

The original `headline.py` is preserved byte for byte as
`source_experiment_v1.py`; its digest still matches the `headline.py` entry in
the original `PREREG.json`. The new `headline.py` verifies that complete
original source map with this one explicit filename mapping, then executes
the same `execute()` function and compares canonical JSON bytes and complete
decompressed outcome bytes. The new wrapper is frozen separately in
`REPLAY_PREREG.json` before its first execution. It also sets Python's integer
serialization cap explicitly for these fixed bounded inputs.

No source formula, inputs, record, gate, calculation or prior receipt changed.
This amendment changes reproduction/custody mechanics only. The original
registration and failed attempts remain available; the amended replay is not
described as a new prospective scientific discovery.
