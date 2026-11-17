# Row 1 interface amendment v2

The [Fable-B review](../../../../openwave-push/reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B/REVIEW.md)
found genuine interface defects in frozen v1. This additive candidate keeps
the old files and review intact.

- `phase_reader.py` reads an independently supplied field. The record
  sets resolution only. Closure, gluing, integrality and unit-circle checks
  are restored. A reversed field can now fail the comparison.
- `closing_reader.py` obtains the regional charge from thirty closing
  contrast/current/coin records alone, without any preparation history,
  initial filling or terminal snapshot. It joins that reading to the
  unchanged native pointer and retains every unsuccessful outcome.
- `PROOF.md` defines the phase regularity assumption and proves the
  connected-graph reconstruction and complete conditional instrument.
- `RESULT.json`, `EXACT_LAWS.json` and `EXAMPLES.json` record executed
  controls. `FINAL_CODE_RECHECK.json` checks the final finite-alphabet guard.

Run `python3 verify.py` in this directory. The initial full verification
covers 144 declared-field reads, 72 independent reversed-field reads,
72 extra-turn comparison rejections and 72 increment-guard rejections.
It deliberately retains 72 indistinguishable high-frequency aliases:
the reader does not authenticate arbitrary fields from finite samples.

Closing-only tests cover all 133 source relative states in 266 independently
executed sweeps, 798 region reads and 2394 complete source/pointer effects.
All four pointer outcomes occur. Exact dynamic programming retains all
2^30 closing coin histories for every source start; 36 exact pointer
aggregate laws sum to one. These are development controls, not reserved
acceptance runs or an independent review of the new code.

The full instrument is conditional on the state at sweep start. It does
not claim a preparation probability for an omitted prefix. It reads the
last recorded time and cannot detect unrecorded later movement. A physical
electric-charge identification, common field coupling and platform award
remain unresolved. No public filing or official row is changed.
