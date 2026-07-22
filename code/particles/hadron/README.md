# Hadron Lane

This directory contains the canonical OPH hadron derivation path inside
`reverse-engineering-reality/code/particles`.

## Current Meaning

The hadron lane is not blocked on a missing symbolic theorem alone. Its supported
frontier is production execution:

- the seeded `N_f = 2+1`, QED-off family is fixed
- the runtime receipt exists
- the payload and evaluation schemas exist
- the missing object is one production backend export bundle with publication-complete manifest provenance and real correlator arrays

So the first promotable outputs on the current branch are isospin-symmetric
`pi_iso` and `N_iso`, not fully physical `p`, `n`, or `pi^0`.

The public output policy for QCD-limited rows is documented in
[`../../../docs/HADRON.md`](../../../docs/HADRON.md). That policy separates source-only
OPH values from OPH plus empirical hadron closure values. The empirical
fine-structure closure uses measured `e+e- -> hadrons` data as the strongest
available substitute for a production OPH hadronic spectral payload.

## Core Files

- runtime receipt:
  [derive_runtime_schedule_receipt_n_therm_and_n_sep.py](derive_runtime_schedule_receipt_n_therm_and_n_sep.py)
- cfg/source payload:
  [derive_stable_channel_cfg_source_measure_payload.py](derive_stable_channel_cfg_source_measure_payload.py)
- sequence evaluation:
  [derive_stable_channel_sequence_evaluation.py](derive_stable_channel_sequence_evaluation.py)
- groundstate readout:
  [derive_stable_channel_groundstate_readout.py](derive_stable_channel_groundstate_readout.py)
- lane audit:
  [derive_current_hadron_lane_audit.py](derive_current_hadron_lane_audit.py)
- production validator:
  [validate_production_hadron_closure.py](validate_production_hadron_closure.py)
- one-shot production writeback:
  [run_production_backend_writeback.py](run_production_backend_writeback.py)
- local diagnostic backend:
  [run_local_diagnostic_backend.py](run_local_diagnostic_backend.py)
- empirical `e+e- -> hadrons` source registry:
  [empirical_ee_hadrons_sources.yaml](empirical_ee_hadrons_sources.yaml)

## Active Artifacts

- [runtime_schedule_receipt_N_therm_and_N_sep.json](../runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json)
- [stable_channel_cfg_source_measure_payload.json](../runs/hadron/stable_channel_cfg_source_measure_payload.json)
- [stable_channel_sequence_evaluation.json](../runs/hadron/stable_channel_sequence_evaluation.json)
- [current_hadron_lane_audit.json](../runs/hadron/current_hadron_lane_audit.json)
- [hadron_production_closure_validation_report.json](../runs/hadron/hadron_production_closure_validation_report.json)
- [hadron_production_readiness_report.json](../runs/hadron/hadron_production_readiness_report.json)
- [HADRON_SYSTEMATICS_STATUS.md](../HADRON_SYSTEMATICS_STATUS.md)

## Canonical Execution Path

From `reverse-engineering-reality/code/particles`:

One-shot path once the backend export exists:

```bash
python3 hadron/run_production_backend_writeback.py \
  --sequence-population runs/hadron/stable_channel_sequence_population.json \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --backend-bundle /path/to/backend_bundle_or_inline_backend_input.json \
  --n-therm 2048 \
  --n-sep 512 \
  --schedule-provenance "external_runtime_commit"
```

Local diagnostic backend path for exercising the full writeback/evaluation
stack without claiming production QCD:

```bash
python3 hadron/run_local_diagnostic_backend.py \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --output /tmp/oph-local-diagnostic-backend.json

python3 hadron/run_production_backend_writeback.py \
  --sequence-population runs/hadron/stable_channel_sequence_population.json \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --backend-bundle /tmp/oph-local-diagnostic-backend.json \
  --dump-output /tmp/backend_correlator_dump.production.json \
  --manifest-output /tmp/oph_hadron_production_backend_manifest.json \
  --evaluation-output /tmp/stable_channel_sequence_evaluation.json \
  --closure-output /tmp/hadron_production_closure_validation_report.json \
  --readiness-output /tmp/hadron_production_readiness_report.json \
  --n-therm 2048 \
  --n-sep 512 \
  --schedule-provenance "local_diagnostic_backend"
```

The local diagnostic backend is target-anchored and carries
`execution_class=diagnostic_surrogate`, so the readiness gate keeps
`publication_bundle_ready=false` even when the numeric pipeline is populated.
Replace that backend bundle with real RHMC/HMC or OPH-hardware correlators,
with `execution_class=production`, before any source-only hadron row can be
promoted.

Backend-side readiness report from the current local state:

```bash
python3 hadron/derive_hadron_production_readiness_report.py
```

Expanded manual path:

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

# If h5py is unavailable on the handoff machine, emit only the production
# manifest and dataset index. The backend owner must create correlators.h5 and
# fill every listed dataset path with real production arrays before ingestion,
# or fill the listed array_file text files and assemble an inline backend export.
python3 hadron/generate_backend_export_bundle_skeleton.py \
  --receipt runs/hadron/runtime_schedule_receipt_N_therm_and_N_sep.json \
  --payload runs/hadron/stable_channel_cfg_source_measure_payload.json \
  --out-dir runs/hadron/production_backend_bundle_request \
  --manifest-only

python3 hadron/assemble_backend_export_from_array_files.py \
  --manifest runs/hadron/production_backend_bundle_request/backend_run_manifest.json \
  --dataset-index runs/hadron/production_backend_bundle_request/dataset_index.json \
  --array-dir runs/hadron/production_backend_bundle_request \
  --output runs/hadron/production_backend_bundle_request/backend_export_inlined.json

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

## License

This particle hadron surface is part of the OPH public repository. See the main
[LICENSE](../../../LICENSE).
