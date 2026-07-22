# Frozen Legacy IBM Quantum Cloud Derived Summaries

This directory contains representative derived JSON result summaries from the
original IBM Quantum Cloud OPH runs. These files preserve provider job IDs,
run settings, measured summary statistics, and selected counts where present.
They are not a complete archive of raw provider counts, and they are frozen as
engineering consistency records rather than OPH-versus-QM evidence.

For the later blinded benchmark, including complete raw joined counts for every
circuit, append-only job journals, hashes, and exact replay instructions, see
[`../runs/issue_509_20260711_v2/`](../runs/issue_509_20260711_v2/).

## Stage 1 Markov / Recoverability

- `stage1/ibm_marrakesh_summary.json`
  `ibm_marrakesh`, job `d6t4da6sh9gc73di7720`
- `stage1/ibm_fez_seed17_summary.json`
  `ibm_fez`, job `d6t4ejngtkcc73cm8l6g`

These files contain the derived tomography analysis: CMI, Petz-recovery
metrics, selected Pauli expectations, exact reference values, and the
qualitative fingerprint checks. They do not contain the complete 27-basis
tomography count dictionaries or reconstructed density matrices. Petz recovery
was calculated offline from the reconstruction; it was not executed as a QPU
recovery circuit.

## Z3

- `z3/ibm_marrakesh_summary.json`
  `ibm_marrakesh`, job `d6t4ic790okc73et0n4g`

This is the reduced-sector sanity check. The main outputs are `t_from_q1`,
`t_from_q2`, leakage, and their agreement across the three prepared `t`
points. The JSON preserves derived probabilities rather than complete provider
count dictionaries.

## Z5

- `z5/ibm_marrakesh_summary.json`
  `ibm_marrakesh`, job `d6t4iiv90okc73et0nbg`
- `z5/ibm_marrakesh_seed17_summary.json`
  `ibm_marrakesh`, job `d6t4rkjbjfas73fp65m0`
- `z5/ibm_fez_seed17_summary.json`
  `ibm_fez`, job `d6t4ipvgtkcc73cm8qg0`
- `z5/ibm_marrakesh_t09_highshot_summary.json`
  `ibm_marrakesh`, job `d6t50e790okc73et186g`

These files expose the measured `delta2_over_delta1` values, `phi^2` target,
leakage, bootstrap intervals, and the extracted `t` consistency checks. They
are derived summaries rather than complete provider count archives.

## S3

- `s3/ibm_marrakesh_layout_diagnostic_summary.json`
  `ibm_marrakesh`, jobs `d6t5c47gtkcc73cm9pd0` and `d6t5cdvgtkcc73cm9pr0`
- `s3/ibm_marrakesh_reversed_seed17_summary.json`
  `ibm_marrakesh`, job `d6t5csmsh9gc73di8e8g`

These files show the layout-sensitive `S_3` behavior directly. The base layout
sits below `2`, while the reversed layout reaches the target neighborhood near
`2`. The layout diagnostic retains target-circuit counts and fitted assignment
matrices, but not every calibration count dictionary; the seed-17 file is a
derived run summary.

## Account Status Snapshot

- `account_status/usage_snapshot_after_s3_cleanup.json`

This is the IBM usage snapshot captured after the main `S_3` cleanup run set.

## License

This public data index is part of the OPH public repository. See the main
[LICENSE](../../../LICENSE).
