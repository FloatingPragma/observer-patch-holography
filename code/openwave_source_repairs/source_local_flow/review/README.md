# Bounded independent source-local energy check

From this review directory, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B check_math.py ../../work/2026-11-16_codex_source_local_count_flow
python3 -B check_clock.py ../../work/2026-11-16_codex_source_local_count_flow
```

Python 3 and NumPy are sufficient. The checker imports no producer module and
does not fetch anything. It verifies existing results by default; use a copy
for an explicitly new `--record` run. The exact producer source/observer
receipts and relevant freezes are bound in `PREREG.json` and `RESULT.json`.
`SCIENTIFIC_REVIEW.md` and `ROUNDING_CHECK.md` state the mathematics, physical
scope and what the bounded replay does not reexecute. All review is same-group.

`CLOCK_RESULT.json` independently binds the common count-clock construction.
`REPLAY_AUDIT.json` records temporary clean-copy verification, unchanged saved
bytes, and rejection of deliberately false named scientific gates. Its
`audit_replay.py` first-record command is intentionally separate from the
read-only headline commands above.
