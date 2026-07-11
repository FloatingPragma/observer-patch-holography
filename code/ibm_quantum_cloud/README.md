# IBM Quantum Cloud Bundle

This directory is the public code-and-data bundle for the IBM Quantum Cloud OPH experiments.

It contains two things:

- `programs/`: runnable preflight, sealing, execution, harvesting, and analysis scripts
- `qc_data/`: representative raw JSON outputs from the completed hardware runs
- `protocols/`: frozen or pre-hardware experimental protocol specifications
- `tests/`: mandatory finite-kernel tests and optional Qiskit/Aer circuit tests

## Scope

These experiments are best read as IBM hardware consistency benchmarks for
OPH-inspired reduced-sector states. Existing `qc_data/` receipts and the issue
`#509` blinded generative benchmark are separate evidence surfaces. The latter
was completed on 11 July 2026, and its sealed execution and result receipts are
published under `runs/issue_509_20260711_v2/`.

They show that:

- the Stage 1 recoverability benchmark behaves in the OPH-favored direction on two real backends;
- the `Z_3` reduced-sector sanity check passes cleanly;
- `Z_5` repeatedly lands near the OPH golden-ratio target on `ibm_marrakesh`;
- `S_3` exposes a real layout-dependent hardware bias, and the corrected layout restores the nonabelian target ratio near `2`.

They do **not** amount to a standalone confirmation of the full OPH framework.
In particular, standard quantum mechanics predicts every programmed dynamic
measurement-and-feedback circuit in the generative benchmark.

## Layout

- `programs/ibm_runtime_common.py`
  Shared IBM Runtime helper utilities.
- `programs/check_ibm_setup.py`
  Auth and backend-discovery verifier.
- `programs/discrete_heatkernel_test.py`
  `Z_3`, `Z_5`, and `S_3` reduced-sector runner.
- `programs/s3_diagnostic_bundle.py`
  `S_3` layout and readout diagnostic bundle.
- `programs/stage1_markov_fingerprint.py`
  3-qubit Markov / recoverability benchmark.
- `programs/generative_repair_kernel.py`
  Pure finite-matrix implementation of the record-gated Cayley repair kernel,
  the matched open-loop heat null, and the dimension-exponent diagnostic.
- `programs/record_gated_cayley_circuits.py`
  Local dynamic-circuit preflight for the self-reading repair instrument. It
  contains no IBM submission path.
- `programs/generative_repair_circuits.py`
  Secondary Metropolis--Hastings `kappa = 0, 1, 2` diagnostic circuits,
  including mandatory abelian non-identifiability controls.
- `programs/blind_preregister.py`
  Opaque circuit catalog, private commit/reveal payload, and digest verifier.
- `programs/blind_preregister_orchestrate.py`
  Two-phase production seal: catalog precommitment followed by an analysis-lock
  binding and atomically committed public/private artifacts.
- `programs/cayley_blind_likelihood_analysis.py`
  Frozen joint-process conditional likelihood ratios, contamination
  sensitivity, calibration/leakage gates, and blinded verdict logic.
- `programs/record_gated_cayley_runtime.py`
  Explicit role-scoped IBM plan, submit, and harvest path with append-only
  hash-chained journals.
- `programs/runtime_analysis_packet.py`
  Fail-closed conversion of the complete two-role runtime journals into the
  blinded analysis packet.
- `protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md`
  Claim boundary, exact predictions, blinding contract, nulls, and failure rule
  for the issue `#509` benchmark.
- `runs/issue_509_20260711_v2/`
  Public manifest, committed reveal, raw provider count journals, sealed blind
  packet/report, job receipts, checksums, and the completed hardware result.
- `qc_data/README.md`
  Index of the public result artifacts.

## Blinded generative repair benchmark

The issue `#509` benchmark freezes five Cayley arms: contemporaneous
`record_gated` repair, matched `open_loop_heat`, and the `delayed_record`,
`shuffled_record`, and `inverted_record` interventions. The secondary
Metropolis--Hastings family compares `kappa = 1` with `kappa = 0` and
`kappa = 2` across identifiable nonabelian spectra and abelian negative
controls.

Every balanced Cayley variant and every Metropolis--Hastings edge circuit uses
`192` shots per backend slot. Separate four-bit basis diagnostic circuits use
`512` shots per prepared code per slot. The two backend/layout roles execute
sequentially: development must be complete before the held-out role can be
submitted, and both roles remain bound to the same sealed catalog and analysis
lock.

The primary evidence measure is a conditional multinomial likelihood ratio,
not a Bayes factor. The primary noise treatment is a frozen contamination
model with mandatory sensitivity channels; four-bit basis calibration is a
diagnostic-only control and is not treated as a factorized causal correction
for feedback errors. A Bonferroni-corrected joint basis-code control against a
`15%` maximum error fraction, a `10%` pooled backend/family leakage ceiling,
and a `25%` catastrophic individual-circuit leakage ceiling are validity gates.
Failure of any makes the run invalid, and no invalid code or shot may be
postselected away.

Before submission, the production-random catalog is committed first and the
hardened anonymous analysis lock is bound second. The public manifest closes
the catalog, normalized-OpenQASM 3 logical-circuit, role/layout,
candidate-table, analysis, private-reveal, and manifest hash chain. Compiled
QPY digests remain transport receipts rather than logical identities; the
mode-`0600` reveal remains outside the checkout until the two-chip blind
analysis is frozen.

### Completed two-QPU result

The production-random `issue509-20260711-v2` run submitted `3,104` circuits
and `606,208` shots in 14 completed IBM Runtime jobs. The development slot was
`ibm_fez` on physical layout `[10, 18, 94, 124]`; only after its completion was
the held-out `ibm_kingston` slot submitted on `[21, 47, 50, 107]`.

The frozen blind report passed every calibration, leakage, completeness, and
shot-conservation gate. On both backends, record-gated repair exceeded the
strict per-backend likelihood-ratio threshold against every frozen controller,
state-preparation, calibrated-noise, and label/layout null. The pooled tests
also exceeded their threshold. The preregistered reference label/layout was
the unique best component among 256 possibilities. After the blind report was
file-hashed, the commitment and all `3,104` normalized-OpenQASM logical-circuit
hashes verified against the reveal.

See [`runs/issue_509_20260711_v2/README.md`](runs/issue_509_20260711_v2/README.md)
for the result table, raw receipts, hashes, and replay commands. This pass is
evidence that the engineered bounded patch performs self-reading repair: it has
local state, a measured boundary record, conditional feedback, final readback,
append-only records, and a public evidence bundle. It is not evidence against
unrestricted quantum mechanics; standard quantum mechanics predicts this
programmed dynamic instrument.

### Local preflight

Create the pinned optional environment:

```bash
python3 -m venv /tmp/oph-ibm-quantum
/tmp/oph-ibm-quantum/bin/python -m pip install -r code/ibm_quantum_cloud/requirements-ibm.txt
```

Run the pure finite-kernel and optional circuit tests:

```bash
python3 -m pytest -q code/ibm_quantum_cloud/tests/test_generative_repair_kernel.py

/tmp/oph-ibm-quantum/bin/python -m pytest -q \
  code/ibm_quantum_cloud/tests/test_generative_repair_kernel.py \
  code/ibm_quantum_cloud/tests/test_record_gated_cayley_circuits.py
```

Generate local validation receipts without credentials or hardware submission:

```bash
python3 code/ibm_quantum_cloud/programs/generative_repair_kernel.py \
  --json-out /tmp/oph-generative-kernel-validation.json

/tmp/oph-ibm-quantum/bin/python \
  code/ibm_quantum_cloud/programs/record_gated_cayley_circuits.py \
  --model s3 \
  --protocol record_gated \
  --shots 192 \
  --json-out /tmp/oph-s3-record-gated-preflight.json
```

These local receipts validate an engineered record-conditioned process against
frozen open-loop and stationary-law nulls. The hardware workflow additionally
requires the two-phase seal, role-scoped sequential execution, complete
hash-chained journals, blind conditional-likelihood report, and verified
reveal required by the issue `#509` acceptance protocol.

A passing hardware result supports only the realization of a blinded finite
self-reading repair instrument against its frozen controller, label, noise,
and decoy nulls. Standard quantum mechanics predicts the programmed instrument,
so the result cannot distinguish OPH from unrestricted quantum mechanics and
does not establish an OPH vacuum, relativity, Standard Model structure, or
cosmology.

## Re-running

The credential file is not part of this repo. The scripts accept either a raw
IBM Cloud API key or a local note file with a line of the form:

```text
IBM cloud API key: <YOUR_API_KEY>
```

Pass it explicitly with `--credentials-file` if it is not in the current working directory.

Example commands:

```bash
python3 code/ibm_quantum_cloud/programs/check_ibm_setup.py --credentials-file /path/to/IBM_cloud.txt

python3 code/ibm_quantum_cloud/programs/discrete_heatkernel_test.py \
  --group z3 \
  --mode hardware \
  --backend ibm_marrakesh \
  --shots 4096 \
  --outdir code/ibm_quantum_cloud/qc_data/z3/rerun_example \
  --credentials-file /path/to/IBM_cloud.txt

python3 code/ibm_quantum_cloud/programs/s3_diagnostic_bundle.py \
  --mode hardware \
  --backend ibm_marrakesh \
  --t-value 0.60 \
  --shots 8192 \
  --use-best-pair-layout \
  --outdir code/ibm_quantum_cloud/qc_data/s3/rerun_example \
  --credentials-file /path/to/IBM_cloud.txt
```

## Notes

- The Stage 1 JSON exports keep the original job IDs and measured outputs, but account-specific instance CRNs are redacted in this public copy.
- The raw JSON is intentionally preserved in its original schema so readers can inspect the measured values directly.

## License And Patent Policy

This code and data bundle is part of the OPH public repository. See the main
[LICENSE](../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../PATENTS.md).
