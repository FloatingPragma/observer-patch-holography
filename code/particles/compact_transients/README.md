# Compact-Transient Receipts

This folder mirrors the compact-record transient simulator contract inside the
paper-stack code tree. It is a receipt scaffold, not a particle-spectrum
prediction and not an OPH confirmation of any transient catalog.

Run:

```bash
python3 particles/compact_transients/build_compact_transient_receipts.py \
  --output particles/runs/compact_transients/receipt_scaffold
```

The generated bundle freezes the CR0-CR4 claim ladder, history schema,
quotient/source/kernel/packet/detector/censoring receipts, FRB control family,
black-hole genealogy and no-generation-leakage guard, refinement/accuracy
contract, and promotion audit.

The default claim is `CR1_QUOTIENT_DIAGNOSTIC`. The generated schema files do
not certify their own readiness gates.

## NON_CLAIMS

The checker verifies presence, byte-level hash integrity, schema conformance,
declared provenance, freshness against recorded inputs, and non-self-authorship.
It does not verify that any value is physically meaningful, correctly derived,
or scientifically true. Requiring an external signature or a reference to an
independently retrievable dataset is an obvious hardening direction that would
make forgery cost more than doing the science, but this checker does not do
either.

CR2 or CR3 evidence may be supplied with:

```bash
python3 particles/compact_transients/build_compact_transient_receipts.py \
  --output particles/runs/compact_transients/receipt_scaffold \
  --evidence-manifest /path/to/external/manifest.json
```

The manifest must use schema
`compact-transient-evidence-manifest-v1` and map each offered receipt to an
external artifact path and SHA-256 content hash. Each artifact declares its
receipt, kind, external producer, current input paths and hashes, assumptions,
units, parameters, outputs, and receipt-specific checked results. The builder
reads and validates that content; filenames and presence alone do not count.
Malformed, stale, wrong-kind, self-authored, and output-directory artifacts
fail closed. See `--help` and the tests in this directory for executable
examples.

The promotion audit records `verdict_checking_tier` as
`SCHEMA_AND_PROVENANCE_ONLY` and `scientific_validation_performed` as `false`,
so downstream readers do not need to inspect this source to distinguish the
implemented checks from scientific validation.

Promotion to `CR3_FROZEN_PHYSICAL_PREDICTION` requires the same evidence
mechanism for frozen controls, refinement stability, and frozen hashes.
Promotion to `CR4_SOURCE_ONLY_OPH_PREDICTION` remains blocked in this builder;
the compact source action, emission microphysics, physical clock, old-host FRB
source theorem, and black-hole genealogy prior are not derived here.
