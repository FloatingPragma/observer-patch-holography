# Mutable current-index scan

Read-only scan on 2026-11-15 during final PR-723 integration. No tracker, row
status or root index was changed by this reviewer. The following findings were
sent to the integrating worker. Line numbers refer to the observed snapshot;
the worker may have corrected them after this record.

| Current file | Observed issue | Requested correction |
| --- | --- | --- |
| `README.md:111,112,117` | Rows 1/2/7 still show M10.2 as open for review | Regenerate after the merged-review status update |
| `notes/ALL_ROWS_STATUS.json:70,724,3369` | Generated cells still show M10.2 open for review; copied next actions retain old adjudication-pending text | Regenerate from current STATUS inputs |
| `notes/onboarding/README.md:9` | Calls PR #723 separate and awaiting review | Record merged arithmetic review with physical cells retained |
| `docs/CANDIDATES.md:36` | Current candidate table calls admission candidate A the pinned candidate | Distinguish admission ancestry from selected A-static-signed-source-v1 |
| `notes/ROWS_01_10_RESEARCH.md:6` | Unqualified statement that no upstream source pin or platform task has changed | Scope statement to the new private follow-up candidates |
| `row_01_charge_quantization/STATUS.json:464` | Top-level platform_pin_moved is false while platform_candidate names adopted static source | Separate adopted platform pin from unpinned new research candidate |
| `row_02_electron_rest_energy/STATUS.json:394` | Same ambiguous false platform_pin_moved, without top-level selected candidate | Add explicit adopted platform source and unpinned research distinction |
| `row_01_charge_quantization/STATUS.json:482` and `row_07_double_cover/STATUS.json:466` | checked_at remains 02:53Z in submission records now describing the 15:17Z merge | Update review observation time or label original snapshot time |

Provenance clarification: the three current `submission.platform_commit` fields
and `notes/UPSTREAM_TRACKER.md:26` retain `72060630...`. If that denotes an earlier
submitted platform head, label it as such and add the final PR head/merge
separately; do not silently rewrite historical evidence. The tracker already
records merge commit `d3029c5e...` and adopted author source `c330c6e7...`.

Exempt from findings: sealed prior packets, explicitly superseded/historical
sections and `prior_active_cases`; claims that the newly implemented private
candidates are unpinned; unresolved physical dictionaries; G6/G9 limitations.
The main upstream tracker, M10.2 filing note and current row README introductions
already correctly distinguish merged arithmetic acceptance from unchanged 🚧
physical statuses.
