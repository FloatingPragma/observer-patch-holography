# Finite-field calculation inputs

These matrices, candidate eigenframes and numerical contractions concern finite charged, gauge and photon constructions. Active exact producers are in model/field_record_addenda and model/native_charged_cluster_minmax. In particular, the latter certifies a rank-two min–max interval but marks the perturbative many-body transfer allowance unproved. Numerical candidate frames and exploratory energies in this collection have no independent certificate.

`input_manifest.json` gives the relative path, source identifier, byte count
and SHA-256 of every input. JSON documents keep their producer-defined
parameter and result fields. JSONL stores one such record per line. NumPy
containers keep named arrays with their original shapes and dtypes; their
schema can be inspected without enabling object deserialization. The file
and directory names distinguish ensemble, geometry, seed and preparation
where the original producer supplies them. Missing metadata is not inferred.

Source text has a `.raw` suffix and is inert input, avoiding accidental test
collection or execution. It is available for independent inspection and
reuse; path constants and external imports in a comparison producer require
an explicit adaptation before execution. Active verified producers are
linked above and run from the canonical checkout. Source comments and
numerical outputs are not an acceptance certificate.

Run `python -m model.finite_repair_comparison_inputs.verify` from the
repository root to check all comparison-input manifests, lengths and hashes
and inspect JSON/NumPy container syntax. This validates input integrity,
without asserting that every numerical interpretation is correct.

Machine-specific path prefixes and compute-host labels are replaced by
symbolic labels where present. The manifest retains both the original
source digest and the public digest for each adapted file; scientific
numbers and arrays are unchanged. Some comparison JSON documents retain
explicit nonstandard nonfinite constants. The verifier counts those cases
and never turns an undefined numerical value into a finite prediction.

An input containing nonfinite constants is marked `raw_nonfinite_values`
in the manifest, with its exact count. Such an input is ineligible for an
automatic prediction or a passing numerical receipt. Integrity verification
preserves those bytes and reports this classification; it does not admit
nonfinite measurements as valid scientific numbers.
