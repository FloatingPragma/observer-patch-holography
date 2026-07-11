# Issue 509 blinded IBM Quantum benchmark

This is the public evidence bundle for the production-random
`issue509-20260711-v2` run completed on 11 July 2026. The frozen blind verdict
is `passes_frozen_reduced_repair_kernel_gate`.

The run used a bounded four-qubit observer-like patch: three qubits carry local
sector state, a measured decision boundary creates a finite record, conditional
feedback applies or withholds repair, final readback closes the patch, and
append-only job/count receipts expose the complete execution. No shots were
dropped or postselected.

## Result

- IBM Runtime jobs: 14, all `DONE`
- Circuits: 3,104 (`3,072` dynamic plus `32` diagnostic calibration circuits)
- Shots: 606,208
- Provider-reported quantum time: 196 seconds
- Development slot: `ibm_fez`, layout `[10, 18, 94, 124]`
- Held-out slot: `ibm_kingston`, layout `[21, 47, 50, 107]`
- Calibration, leakage, completeness, and shot-conservation gates: all pass
- Primary S3 endpoint: pass on both backends and pooled
- Global label/layout test: the committed reference is the unique winner among
  256 components

The table reports conditional log-likelihood ratios for record-gated repair
over each frozen null. The strict thresholds are `log(10)` per backend and
`log(100)` pooled. Exponentiated ratios overflow for the large entries; the
label-mixture ratio remains finite at approximately `256`. Every frozen
contamination-sensitivity lower bound also passes.

| Null | `ibm_fez` | `ibm_kingston` | Pooled |
|---|---:|---:|---:|
| lazy heat | 136,738.537 | 141,443.849 | 278,182.386 |
| delayed record | 77,145.222 | 78,954.396 | 156,099.617 |
| shuffled record | 93,252.706 | 96,396.875 | 189,649.582 |
| inverted record | 275,461.330 | 284,142.447 | 559,603.778 |
| state preparation only | 99,891.886 | 103,598.342 | 203,490.228 |
| calibrated noise | 171,961.767 | 180,409.442 | 352,371.209 |
| 256-component label/layout mixture | 5.545 | 5.545 | 5.545 |

The 99% global simultaneous prediction envelope passed with zero violations
across 55,296 primary cells. Pooled Cayley leakage was `3.362%` on `ibm_fez`
and `1.189%` on `ibm_kingston`, below the frozen 10% limits. The worst
individual circuit was `19.792%`, below the separate 25% catastrophic limit.
The largest observed four-bit basis error was `3.516%` on `ibm_fez` and
`2.930%` on `ibm_kingston`, below the 15% diagnostic gate.

The secondary Metropolis--Hastings family behaved as designed: S3, A4, and the
seeded nonabelian decoy favored the frozen `kappa = 1` spectrum over `kappa = 0`
and `kappa = 2`; the Z3 and Z5 abelian negative controls returned exact zero
log-likelihood separation for those comparisons. These secondary results do
not select the primary verdict.

## Blind chronology and commitments

The public catalog, analysis code, thresholds, nulls, two anonymous backend
slots, and 3,104 normalized-OpenQASM logical-circuit digests were committed
before submission. The development role completed before any held-out job was
submitted. Raw results were then sealed into one blinded packet and analyzed
without opening the reveal.

The blind report was frozen at `2026-07-11T06:17:04Z` with:

- catalog precommitment: `ae335e03ca9a70dc6b611f3e95151af8dc4e5b63de1836625d2bce60dc5f0e9c`
- manifest: `78868d2ea896e1fad2b032114fd2b103ad7564ca2ceff072bac7444574bc7ff3`
- analysis lock: `d5658535a03a8f89159554013256f7e0e024f175a4c605d0393872d7a27577fb`
- data packet: `16e1272ff6c76eac6f0c7bf8d3c141245744820cf5409bfa6bc75796e6f63691`
- blind report: `d6c7373ec7a21ac1194144e976864122fbddcc10e9abf47df2855a0e7539b49e`
- blind-report file: `9d60d6939a41f59184b2dee859a5fa5b223eff6d6309a223a1f16ab25e8317ab`

Only after that freeze did the verifier open the reveal. Its secret commitment,
analysis document, backend/layout map, opaque IDs, candidate recipes, and every
logical-circuit digest verified. The unique blind winner
`lbl_0625bffe5d250e6ca1d50df25e68a991` was exactly the committed reference
component.

An earlier `v1` seal is documented in `v1_pre_submission_abort.json`. It was
aborted before any provider submission after a fresh-process check proved that
QPY bytes were not a canonical cross-process logical identity. The production
`v2` seal instead commits normalized OpenQASM 3; QPY hashes are retained only
as per-execution transport receipts.

## Files

- `public_manifest.json`: pre-submission public catalog and thresholds
- `analysis_lock.json.gz`: anonymous frozen candidate tables and analysis lock
- `reveal.json.gz`: post-analysis commit/reveal payload and circuit recipes
- `runtime_*_submission_events.ndjson.gz`: hash-chained submission receipts
- `runtime_harvest_events.ndjson.gz`: provider metadata and raw joined counts
- `blinded_hardware_data_packet.json.gz`: complete fail-closed analysis packet
- `blind_hardware_report.json.gz`: full frozen blind report
- `blind_freeze_receipt.json`: pre-unblinding file-hash receipt
- `job_receipts.json`: concise IBM job and usage index
- `result_summary.json`: machine-readable result and claim boundary
- `CHECKSUMS.json`: compressed and original-file SHA-256 values

The JSON and NDJSON evidence is compressed with deterministic `gzip -9 -n`.
The submission journals bind every opaque ID to its normalized-OpenQASM logical
digest and execution-only compiled-QPY digest. The reveal contains the recipes
needed to rebuild the logical circuits; the harvest journal contains all raw
joined counts and IBM Runtime job metadata.

## Verification

From the repository root, first check the compressed artifacts:

```bash
gzip -t code/ibm_quantum_cloud/runs/issue_509_20260711_v2/*.gz
```

Extract the committed documents and journals into a scratch directory:

```bash
mkdir -p /tmp/oph-issue509-v2
gzip -dc code/ibm_quantum_cloud/runs/issue_509_20260711_v2/analysis_lock.json.gz > /tmp/oph-issue509-v2/analysis_lock.json
gzip -dc code/ibm_quantum_cloud/runs/issue_509_20260711_v2/reveal.json.gz > /tmp/oph-issue509-v2/reveal.json
gzip -dc code/ibm_quantum_cloud/runs/issue_509_20260711_v2/runtime_development_submission_events.ndjson.gz > /tmp/oph-issue509-v2/development.ndjson
gzip -dc code/ibm_quantum_cloud/runs/issue_509_20260711_v2/runtime_heldout_submission_events.ndjson.gz > /tmp/oph-issue509-v2/heldout.ndjson
gzip -dc code/ibm_quantum_cloud/runs/issue_509_20260711_v2/runtime_harvest_events.ndjson.gz > /tmp/oph-issue509-v2/harvest.ndjson
```

With the pinned optional environment from the parent README, verify the reveal
and rebuild all logical circuits:

```bash
PYTHONPATH=code/ibm_quantum_cloud/programs \
  /tmp/oph-ibm-quantum/bin/python \
  code/ibm_quantum_cloud/programs/blind_preregister.py verify \
  --public code/ibm_quantum_cloud/runs/issue_509_20260711_v2/public_manifest.json \
  --reveal /tmp/oph-issue509-v2/reveal.json \
  --analysis-document /tmp/oph-issue509-v2/analysis_lock.json \
  --rebuild-circuits
```

Reseal the raw journals into an independent packet and rerun the frozen blind
analysis:

```bash
PYTHONPATH=code/ibm_quantum_cloud/programs \
  /tmp/oph-ibm-quantum/bin/python \
  code/ibm_quantum_cloud/programs/runtime_analysis_packet.py \
  --manifest code/ibm_quantum_cloud/runs/issue_509_20260711_v2/public_manifest.json \
  --analysis-lock /tmp/oph-issue509-v2/analysis_lock.json \
  --reveal /tmp/oph-issue509-v2/reveal.json \
  --development-submission-journal /tmp/oph-issue509-v2/development.ndjson \
  --heldout-submission-journal /tmp/oph-issue509-v2/heldout.ndjson \
  --harvest-journal /tmp/oph-issue509-v2/harvest.ndjson \
  --created-utc 2026-07-11T06:14:46.204442+00:00 \
  --json-out /tmp/oph-issue509-v2/replayed_packet.json

PYTHONPATH=code/ibm_quantum_cloud/programs \
  /tmp/oph-ibm-quantum/bin/python \
  code/ibm_quantum_cloud/programs/cayley_blind_likelihood_analysis.py \
  --lock /tmp/oph-issue509-v2/analysis_lock.json \
  --data /tmp/oph-issue509-v2/replayed_packet.json \
  --out /tmp/oph-issue509-v2/replayed_report.json
```

The replayed internal packet and report hashes must equal the values above.

## Claim boundary

This is hardware evidence for an engineered finite self-reading repair
instrument against the frozen reduced controller nulls. It demonstrates an
observer-like implementation with local state, ports or boundaries, readback,
records, feedback, and a public evidence bundle on two IBM QPUs.

It is not evidence against unrestricted quantum mechanics and does not
distinguish OPH from standard quantum theory. Standard quantum mechanics
predicts the programmed dynamic circuit. The result does not establish OPH
claims about relativity, the Standard Model, cosmology, or ontology. The
reported evidence values are conditional likelihood ratios, not Bayes factors,
because calibration uncertainty was stress-tested by frozen sensitivity bounds
rather than fully marginalized.
