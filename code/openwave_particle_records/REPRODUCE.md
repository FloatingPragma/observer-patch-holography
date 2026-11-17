# Evidence verification and scientific replay

Run the package verifier from the RER root. It checks exact archived bytes,
membership in each original manifest and critical receipt distinctions.
`--workspace ..` additionally checks all original manifest leaves. Neither
command executes native repairs or grants physical acceptance.

Scientific replays require the complete sibling OpenWave workspace or the
specified standalone tree. Original artifacts are immutable. Use a disposable
copy for any script which writes receipts or timings.

| Evidence | Scientific entry point relative to workspace | Execution rule |
| --- | --- | --- |
| Common static source | `openwave/shared/static_candidate_join/publication/` | Copy the entire standalone tree and follow its README from the copied root |
| Compact internal clock | `openwave/shared/compact_internal_clock_snapshot_adapter/replay.py` | Use the reader's exact theorem snapshots |
| C/CAR compact join | `openwave/shared/compact_field_spin_join/existing_audit/` | Follow the sealed packet and root-review commands |
| Current and loop clocks | `openwave/shared/native_current_clock/`, `openwave/shared/native_field_cycle/` | Use the documented headline/snapshot readers and retained amendments |
| First return | `openwave/shared/native_field_cycle/first_return/` | Preserve both the amended neutral-null test and original failed control |
| Neutral flavour | `openwave/shared/neutral_flavour_join/` | Headline is read-only; run the receipt-writing apparatus audit in a complete disposable copy |
| Conditional boundary field | `openwave/shared/neutral_source_completion/headline_linked.py` | Complete linked input verification is part of the reader |

The source proofs, receipts and implementation extracts in `evidence/` preserve
their original bytes and terminology. They are archival evidence, not mutable
package instructions. The unexecuted logarithmic-flow, symmetric-source and
contact-energy proposals are listed separately in `ASSESSMENT.json`; no
working snapshot is treated as an experimental seal.
