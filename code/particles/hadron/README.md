# Hadron Lane

This directory contains the canonical OPH hadron derivation path inside
`reverse-engineering-reality/code/particles`.

## Current Meaning

The hadron lane is not blocked on a missing symbolic theorem alone. Its honest
frontier is production execution:

- the seeded `N_f = 2+1`, QED-off family is fixed
- the runtime receipt exists
- the payload and evaluation schemas exist
- the missing object is the real production backend correlator dump

So the first promotable outputs on the current branch are isospin-symmetric
`pi_iso` and `N_iso`, not fully physical `p`, `n`, or `pi^0`.

## Core Files

- runtime receipt:
  [derive_runtime_schedule_receipt_n_therm_and_n_sep.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py)
- cfg/source payload:
  [derive_stable_channel_cfg_source_measure_payload.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/derive_stable_channel_cfg_source_measure_payload.py)
- sequence evaluation:
  [derive_stable_channel_sequence_evaluation.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/derive_stable_channel_sequence_evaluation.py)
- groundstate readout:
  [derive_stable_channel_groundstate_readout.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/derive_stable_channel_groundstate_readout.py)
- lane audit:
  [derive_current_hadron_lane_audit.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/derive_current_hadron_lane_audit.py)
- production validator:
  [validate_production_hadron_closure.py](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/hadron/validate_production_hadron_closure.py)

## Active Artifacts

- [runtime_schedule_receipt_N_therm_and_N_sep.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json)
- [stable_channel_cfg_source_measure_payload.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/hadron/stable_channel_cfg_source_measure_payload.json)
- [stable_channel_sequence_evaluation.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/hadron/stable_channel_sequence_evaluation.json)
- [current_hadron_lane_audit.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/hadron/current_hadron_lane_audit.json)
- [hadron_production_closure_validation_report.json](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/runs/hadron/hadron_production_closure_validation_report.json)
- [HADRON_SYSTEMATICS_STATUS.md](/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/code/particles/HADRON_SYSTEMATICS_STATUS.md)

## Canonical Execution Path

From `reverse-engineering-reality/code/particles`:

```bash
python3 hadron/derive_runtime_schedule_receipt_n_therm_and_n_sep.py \
  --n-therm 2048 \
  --n-sep 512 \
  --schedule-provenance "external_runtime_commit" \
  --output runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json

python3 hadron/derive_stable_channel_cfg_source_measure_payload.py \
  --runtime-receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --output runs/hadron/stable_channel_cfg_source_measure_payload.json

python3 hadron/generate_backend_export_bundle_skeleton.py \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --out-dir /path/to/backend_bundle

python3 hadron/normalize_backend_export_bundle.py \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --backend-bundle /path/to/backend_bundle \
  --output runs/hadron/backend_correlator_dump.production.json \
  --manifest-output runs/hadron/oph_hadron_production_backend_manifest.json

python3 hadron/derive_stable_channel_cfg_source_measure_payload.py \
  --runtime-receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --backend-dump runs/hadron/backend_correlator_dump.production.json \
  --output runs/hadron/stable_channel_cfg_source_measure_payload.json

python3 hadron/derive_stable_channel_sequence_evaluation.py \
  --runtime-receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --cfg-source-payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --backend-dump runs/hadron/backend_correlator_dump.production.json \
  --output runs/hadron/stable_channel_sequence_evaluation.json

python3 hadron/validate_production_hadron_closure.py \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --dump runs/hadron/backend_correlator_dump.production.json \
  --evaluation runs/hadron/stable_channel_sequence_evaluation.json \
  --output runs/hadron/hadron_production_closure_validation_report.json

python3 hadron/derive_current_hadron_lane_audit.py
python3 hadron/render_hadron_systematics_status.py
```

## Focused Checks

```bash
python3 -m pytest \
  hadron/test_runtime_schedule_receipt_n_therm_and_n_sep.py \
  hadron/test_stable_channel_sequence_evaluation.py \
  hadron/test_current_hadron_lane_audit.py
```
