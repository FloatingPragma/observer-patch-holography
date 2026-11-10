# Native defect receipt custody

The ledger consumes `native_defect_decision_live_results.json`, produced by
the unchanged `native_defect_decision.report()` calculation. Every scientific
field equals the immutable `native_defect_decision_results.json`; only source
hashes can differ. The retained receipt, completion campaign and quantum-defect
receipt keep their recorded provenance. Their historical documentation pins are
not represented as hashes of the working tree.

After editing a source document, regenerate the live receipt from the repository
root with `python3 tools/renew_native_defect_live.py --renew`. Run
`python3 tools/renew_native_defect_live.py --check` for a fresh producer replay
and exact payload comparison. The tool rejects changed scientific fields,
source coverage or retained evidence.
