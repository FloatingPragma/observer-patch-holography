# IBM Quantum Cloud Bundle

This directory is the public code-and-data bundle for the IBM Quantum Cloud OPH experiments.

It contains two things:

- `programs/`: the final runnable scripts used for the public hardware benchmark bundle
- `qc_data/`: representative raw JSON outputs from the completed hardware runs
- `protocols/`: frozen or pre-hardware experimental protocol specifications
- `tests/`: mandatory finite-kernel tests and optional Qiskit/Aer circuit tests

## Scope

These experiments are best read as IBM hardware consistency benchmarks for OPH-inspired reduced-sector states.

They show that:

- the Stage 1 recoverability benchmark behaves in the OPH-favored direction on two real backends;
- the `Z_3` reduced-sector sanity check passes cleanly;
- `Z_5` repeatedly lands near the OPH golden-ratio target on `ibm_marrakesh`;
- `S_3` exposes a real layout-dependent hardware bias, and the corrected layout restores the nonabelian target ratio near `2`.

They do **not** yet amount to a complete standalone confirmation of the full OPH framework.

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
- `protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md`
  Claim boundary, exact predictions, blinding contract, nulls, and failure rule
  for the issue `#509` benchmark.
- `qc_data/README.md`
  Index of the public result artifacts.

## Generative repair preflight

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
  --shots 64 \
  --json-out /tmp/oph-s3-record-gated-preflight.json
```

These receipts validate an engineered record-conditioned process against
frozen open-loop and stationary-law nulls. Standard quantum mechanics predicts
the programmed dynamic circuit, so a passing result is not evidence against
ordinary quantum mechanics.

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
