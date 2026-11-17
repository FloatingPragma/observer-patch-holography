# Verification amendment ledger

1. Parent `/root` reports a static coverage gap in the sealed admission reader:
   exact column/input multiplicities are not required. No producer/source edit
   is authorized or made.
2. `GAP_PREREG.json` freezes the old reader/evidence/manifest and new mutation
   generator before running it. SHA-256:
   `f6b4a8ceadde26ced2dca66c6fe6e16f8cfcaf0ef18e6d0238d720b85f3b2007`.
3. `GAP_PROBE.json` records four malformed catalogs/domains accepted by the old
   reader and one uncompensated dropped row rejected. All original bytes stay
   fixed. The five-check probe takes 36.856 seconds including disposable JSON
   serialization and subprocess work. This confirms a verifier gap, not a
   failed native source theorem.
4. The additive reader derives case identity from actual census/coin records,
   requires exact catalog/domain multiplicities and preserves the original
   source-law/prefix checks. `AMENDMENT_PREREG.json` freezes the repair before
   its first execution. SHA-256:
   `4150b655c8391eccdf24f4f7e9a50554502e4cf6b5d03896f354dd3d04e1e267`.
5. `VALIDATION.json` retains all amended-check dispositions. `HEADLINE.json`
   records the standalone read-only execution and runtime; `CUSTODY.json`
   binds the unchanged ancestor/source and the amendment difference.

All times and file hashes in the JSON registrations are authoritative.
Same-group Codex work only; no external review or official mark is inferred.
