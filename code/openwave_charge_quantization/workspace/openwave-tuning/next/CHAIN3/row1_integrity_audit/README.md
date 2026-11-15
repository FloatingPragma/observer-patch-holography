# Independent row1 integrity audit

Run with Python 3 and its standard library:

```sh
python3 openwave-tuning/next/CHAIN3/row1_integrity_audit/audit.py
python3 openwave-tuning/next/CHAIN3/row1_integrity_audit/branch_tree_counts.py
python3 openwave-tuning/next/CHAIN3/row1_integrity_audit/seal.py
```

The command verifies frozen inputs and independently checks the native
kernel, complete probability laws, actual pointer records, source latches,
and adversarial controls. It writes `RESULT.json`, `SUMMARY.json`, and
`PINNED_INPUTS.json` only in this directory. `PROOF.md` states the method and
the conditional scope. The second command independently verifies all sixty
reported quotient-tree leaf counts using a path-count recurrence and writes
`EXACT_TREE_COUNTS.json`; it does not re-execute every producer leaf tape.
It also retains every independently compared native execution, including
all timeout records, in `NATIVE_POINTER_RECORDS.json`.
The final command verifies cross-artifact consistency, executes the portable
producer and facade in an isolated standard-library process, and seals the
audit manifest.

Use the additive `observer_only.admit_record_only(claim, G, region, retries,
producer=integrity)` for record-only admission. The original producer's
required pointer final-load metadata narrows timeout candidate columns.
The facade ignores that supplied snapshot, supplies synthetic schema
placeholders, and retains every initial column compatible with the tape.
Frozen producer and AU5 code remain unchanged; no official row is awarded.
