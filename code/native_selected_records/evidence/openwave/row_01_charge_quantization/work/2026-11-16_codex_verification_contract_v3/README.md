# Verification contract V3

This packet repairs execution interfaces and tests rejection behavior. It adds
no physical result or award. The sealed b6df source and outgoing V2 files are
unchanged. Read [REVIEW.md](REVIEW.md) for findings and exact scope;
[FINAL_VERDICT.json](FINAL_VERDICT.json) binds the retained process evidence.

The five complete, inspectable calculations are in `descendants/scripts/`.
Their scientific calculation ASTs are unchanged after normalizing the explicit
assert-to-unconditional-check replacement and final output interface. They
require Python 3; energy/common additionally require NumPy. Charge/spin/pair
use only the standard library. No network or installation is performed.

From this directory, default verification reads and compares the complete saved
JSON result without changing it:

```sh
python3 -B descendants/scripts/m10_4_charge.py
python3 -B descendants/scripts/m10_4_energy.py
python3 -B descendants/scripts/m10_4_spin.py
python3 -B descendants/scripts/m10_4_common.py
python3 -B descendants/scripts/m10_4_pair.py
```

Each script also accepts `--verify PATH`. A missing or differing saved result
fails. `--record-to NEW_FILE` is an explicit first-record operation and refuses
overwrite. Scientific checks execute unconditionally even under `python -O`;
the final record operation is reached only after those checks pass. A new
record is reproduction of known results, never retrospective discovery.

Repaired direct producers, including their calculation dependencies and original
saved evidence, are available separately:

```sh
python3 -B descendants/producers/reservoir/experiment.py
python3 -B descendants/producers/pair/pair_source.py
```

Their default mode verifies the saved complete evidence. `--verify DIR` changes
the read-only reference directory. `--record-to NEW_DIR` explicitly creates a
new record directory and refuses overwrite. The reservoir producer retains an
explicitly requested failed scientific record and exits nonzero; it never
turns a failed gate into a successful process. The pair comparison excludes
only runtime from the full source result and separately checks saved receipt
metadata, source digest and column count.

The new historical-source harness is a full archive replay, **not** the small
scientific inspection payload. Its source argument is the `source/` directory
of the preserved b6df export:

```sh
python3 -B descendants/verify_source.py --source /PATH/TO/b6df/source
```

It runs the 28 declared historical jobs, compares every returned scientific
field with the existing clean-replay record, checks current verdicts using
explicit job-specific paths, and verifies the complete source tree has unchanged
bytes, sizes, mtimes and inodes. Only prebound runtime paths are ignored. It
forces assertion-enabled Python in every child and descendant subprocess,
including when the caller uses optimization. Each job has a 60-second limit;
the complete suite took 92.56 seconds wall time and is not represented as a
single under-60-second headline. `--job INDEX` selects an explicitly scoped
subset. `--record-to NEW_FILE` may create a new replay receipt outside the saved
source; it never creates or replaces scientific expectations.

The schemas and expected values in `descendants/SOURCE_JOB_CONTRACTS.json` are
copied from the existing clean b6df replay. They are not fresh predictions.
New V3 scientific work requires its own registered code, data and contracts;
do not replace these historical expectations with the descendant's output.

To verify this audit's retained results without rerunning the science or
writing the verdict:

```sh
python3 -B adjudicate.py
```

The audit-generation programs are retained for full inspection and repeatable
test construction. They create new receipt directories exclusively and therefore
must be run in a fresh disposable copy. They are not the public scientific
verification commands. No command should be run in a way that overwrites this
sealed packet. Full replay and targeted mutation receipts, including original
failures, are retained under the `validation_*`, `original_processes` and
`wrapper_processes` directories. Original computation inputs for counterfactual
gate tests are retained in the original source; gate changes do not alter their
generated scientific records. Tests credit execution/integrity only.

Final integration needs a new source export and a renewed affected checklist.
Historical source producers and V2 interfaces remain historical; they have not
been silently repaired. Same-group validation is not external G6 or release
authorization.
