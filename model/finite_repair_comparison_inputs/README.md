# Finite repair and wave comparison inputs

These finite geometry, repair, coherent-walk and detector experiments retain their numerical arrays and source definitions. The SIM/cf1 data distinguish declared coherent propagation, gauge ensembles, generation-register controls and observer comparisons. Source preparation, address count, clock and frame rules belong to each producer; equal numerical labels do not identify their observables. Model/repair_sector_bridge supplies the exact addressed-current observer and mixing theorem, and model/observer supplies the declared native implementation.

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

The controls directory retains finite native/pump, gauge-apparatus,
state-geometry and protection/reference countercontrols. The seeds retain
small source checks. Their exact construction and preparation are part of
the source text; numerical agreement alone supplies no source-law transfer.

An input containing nonfinite constants is marked `raw_nonfinite_values`
in the manifest, with its exact count. Such an input is ineligible for an
automatic prediction or a passing numerical receipt. Integrity verification
preserves those bytes and reports this classification; it does not admit
nonfinite measurements as valid scientific numbers.
