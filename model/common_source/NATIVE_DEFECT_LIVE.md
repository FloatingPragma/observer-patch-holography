# Native defect receipt custody

The canonical `native_defect_decision_results.json` and the ledger projection
`native_defect_decision_live_results.json` use the unchanged
`native_defect_decision.report()` calculation. Their scientific fields equal
the byte-authenticated receipt under
`../completion/historical/r2063_live_inputs/model/common_source/`.
That archived receipt retains its documentary pin; it is not a claim of
current source custody. The frozen clock campaign explicitly permits the
canonical control's one documentation digest and its completion-campaign
reference to renew. Its scientific parameters and executable sources are fixed.

Run `python3 tools/renew_native_control_canonical.py --renew` for the canonical
control, admitted campaign reference and quantum-defect child. The command
executes the unchanged producers and rejects every non-custody difference.
Run `python3 tools/renew_native_defect_live.py --renew` for the ledger projection.
Both commands accept `--check` for producer replay without writing. Consumers
of changed canonical receipt hashes require their own producer-backed renewal.
