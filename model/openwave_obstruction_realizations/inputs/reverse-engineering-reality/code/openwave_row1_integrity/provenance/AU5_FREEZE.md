# FREEZE: AU5 validators

Frozen 2026-11-14. Directory `openwave-push/audits_row1/AU5_validators/`. The two modules use the standard
library only, read nothing from disk and import nothing of the other team.

## Pinned object (checked equal before and after every run, and again when this file was written)

| File under `openwave-tuning/next/ROW1_NATIVE/` | SHA-256 |
| --- | --- |
| `routes/multiload/multiload_readout.py` | `442d697747ee4adbbcafea5c0a5676e848cbb94173c652042893964145d94eef` |
| `routes/multiload/MULTILOAD_READOUT.json` | `b70871dbefcb731cb7b956a7e8a8b89c722fbfcc49781f571c92a167b73f8e1f` |
| `routes/auditor_bridge/native_cycle_meter.py` (pointer policy) | `4db8d707879c416ca5002b782a70e07aa71d4d9a7fe7111b76a64828c1a05329` |

## Files of this directory

| File | SHA-256 | Bytes |
| --- | --- | --- |
| `PROOF.md` | `b52a6e90393b197fa652d14900ebb45ecb75d543a38dba075a77d63bfcbe162d` | 10040 |
| `RESULTS.json` | `3fdf033339cd93f4d1f4bb5769a5e51127358d1c5c24eb56469d45476d585eb6` | 13354 |
| `RESULTS_ZERO_ECHO.json` | `9e1b25160e43941050e98f0d82af0a7a2f55ac99d0d340db6a9f8ed835aed615` | 614 |
| `REVIEW_FOR_CODEX.md` | `4a16190a359c048270bf095f94f2c06a8b58e23e5b70528c1758129ef8ad4709` | 7126 |
| `closing_sweep.py` | `0f1ff01daf8aa2afa1c351eb93cdc6243bb1cbb6cced82b9be950baff47155ce` | 11150 |
| `logs/run.log` | `fa27ce675bf4211bf8031525cfe61e6fc77f2e9f5b5aaebf6f38b5d5f5052d02` | 8836 |
| `logs/run_zero_echo.log` | `64ecc90b54ec768abbb147f447f90efc88f37e1017ad3bb7daf37c41a6ecbd32` | 622 |
| `run.py` | `d9f2d6db709694eedf44787f7773aeddde25afd6a3b75a774834d3d8ef57ed9f` | 36135 |
| `run_zero_echo.py` | `639046e97e2e785de309877376d341f52f1bed6a7cab05a98f9c494f8ef85167` | 4159 |
| `strict_tape_validator.py` | `6c2f1501a4596ba1a09e1e203ea9c0f2411f33286e8b77cc8467aaeacd0c576a` | 21737 |

## Callable signatures

```python
# strict_tape_validator.py
validate(result, inputs=None) -> (ok: bool, reasons: list[str])
    # result: dictionary returned by multiload_readout.execute, with result['inputs'] or `inputs` =
    #         {'G', 'schedule', 'source_coins', 'region', 'pointer_initial', 'retries'}
    # each reason starts with the rule that failed: R0 ... R7
with_inputs(result, G, schedule, source_coins, region, pointer_initial, retries=16) -> dict
validate_pointer_tape(meter, pointer_initial, commanded, retries) -> (reasons: list[str], replay: dict | None)
native(loads, seam, coin) -> (loads, current)
route_to_zero(start) -> tuple        reset_success(port, retries) -> Fraction        attempt_cap(retries) -> int

# closing_sweep.py
closing_sweep(G, retained_tape, source_state, order, coins, *, channel=None, allow_partial=False) -> dict
    # keys: records, comparison, source_before, source_after, reconstruction_before, reconstruction_after,
    #       reconstruction_from, mismatch_positions, first_mismatch, non_native_records, swept_seams,
    #       swept_graph_connected_and_spanning, sweep_moving_events, source_changed_by_sweep, verdict,
    #       latch, audit_only
    # verdict: 'consistent' | 'reconstruction_error_detected' | 'sweep_record_not_native'
    #          | 'inconclusive_sweep_not_spanning' (only with allow_partial=True)
    # channel (tests only): function(record, reconstruction) -> received record; None = faithful reception
sweep_step(source, reconstruction, seam, coin, position=0, channel=None)
    -> (record, source_after, reconstruction_after, mismatch, native_ok, reconstructed_contrast)
reconstruct(G, retained_tape) -> tuple           tape_per_seam(retained_tape) -> list
native_event(loads, seam, coin) -> (loads, contrast, current)
order_family() -> dict[str, list[int]]           seeded_order(seed) -> list[int]
connected_and_spanning(seams) -> bool
```

## Commands

```sh
cd /Users/muellerberndt/Projects/oph-meta/openwave-push/audits_row1/AU5_validators
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY=/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/.venv/bin/python
$PY run.py > logs/run.log 2>&1; echo "exit=$?" >> logs/run.log                        # 53.8 s CPU, RESULTS.json
$PY run_zero_echo.py > logs/run_zero_echo.log 2>&1; echo "exit=$?" >> logs/run_zero_echo.log   # RESULTS_ZERO_ECHO.json
```

Both scripts stop at once if the pinned files differ from the hashes above. They read the pinned modules at
their paths and write only in this directory. Interpreter: Python 3.12.6. No configuration of total load four
and no reserved or acceptance configuration is executed.

## Use by another party

```python
import strict_tape_validator as sv, closing_sweep as cs
ok, reasons = sv.validate(sv.with_inputs(result, G, schedule, source_coins, region, pointer_initial, retries))
report = cs.closing_sweep(G, retained_tape, source_state, cs.order_family()["identity"], coins)
```
