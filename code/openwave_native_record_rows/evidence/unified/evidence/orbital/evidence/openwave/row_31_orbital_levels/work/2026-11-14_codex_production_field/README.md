# Native global production field: corrected and fully retained

The amended run executes347533 native calls and produces outer occupation
fields, native band readouts, four-channel curl and stopping clocks on the
same histories. The population operator derives from native first-hit
feedback; it is not installed in place of the source.

The original producer checks and byte replay passed, but independent review
found a real band-decoder error: in-place normalization of a NumPy column
view modified the stored projector. `failed_initial/` preserves that entire
attempt, its producer-pass result and the failed review finding. The
amendment uses an explicit copy and independently verifies projector
symmetry, idempotence, trace, resolution and polynomial identity before
reading records. Source, seeds, caps and scientific gates are unchanged.

CORRECTION_COMPARISON.json proves every raw event history, occupation,
clock and constant/curl readout stayed identical, while all eight band
readouts were corrected. The receipt's only change is its raw archive
content hash; population consequences are unchanged. Current PREREG.json
is the pre-retry amendment; the original is retained under failed_initial.
Read both attempts. A local corrected pass is not an official cell award.
