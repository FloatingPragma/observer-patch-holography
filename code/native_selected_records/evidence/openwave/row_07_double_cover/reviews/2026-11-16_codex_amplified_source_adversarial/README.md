# Amplified source: independent same-group review

Read [REVIEW.md](REVIEW.md) for the scientific verdict, hypotheses, repairs and
scope. No external G6 or official award is granted. [ALGEBRA_REVIEW.md](ALGEBRA_REVIEW.md)
gives the source-word derivation; the preliminary review is retained as history.

The full 256-case audit passed, including all retained states, flags and costs.
It took 264.555s. The bounded complete-workspace and signed-clock interface check
took 23.264s. These counts and runtimes are different from the producer's small
compact headline or its 172.041s full run.

From this directory, using Python 3 and NumPy, the following commands are
read-only by default; they compare saved review evidence and never overwrite it:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B phase_tests.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B workspace_bridge.py ../../work/2026-11-16_codex_stationary_readout_repair_v2/CONFIRMATION
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B audit_records.py ../../work/2026-11-16_codex_stationary_readout_repair_v2/CONFIRMATION
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B mutation_check.py ../../work/2026-11-16_codex_stationary_readout_repair_v2/CONFIRMATION
```

`audit_catalog.py` independently reconstructs complete circuits from the two
header artifacts; `microtrace.py` and `pair_source_check.py` are separately
registered component checks. Their CLI help gives explicit inputs. New first
receipts require an explicit exclusive `--record-to` path. Existing paths fail.
All guards are unconditional under `-O`.

`RECORD_AUDIT.json`, `CATALOG_AUDIT.json`, `WORKSPACE_BRIDGE.json`,
`MUTATIONS.json`, `PHASE_TESTS.json`, `MICROTRACE.json` and
`PAIR_SOURCE_REVIEW.json` retain the full review results. The independent
reader imports only reviewer helpers. Full source custody is bound by
`BINDINGS.json`; `MANIFEST.sha256` seals this review when available.
