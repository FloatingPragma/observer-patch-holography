# DELIVERABLE: what every agent produces, and what OpenWave acceptance is

One page. Binding for every agent of every party. Details: [ACCEPTANCE.md](ACCEPTANCE.md) (the platform's
terms with quotations), [FRAMEWORK.md](FRAMEWORK.md) (gates and collaboration).

## 1. Acceptance, in one paragraph

A row is **accepted** when the OpenWave maintainers set its cell in `MODELS.md` to ✅. They do that, one cell at
a time, for a registered task of column M10 whose script they can run, when the row's observable is **computed
from the records of observers running OPH's pinned native process**, every clause of the row's simplest test
holds in the full sense, and nothing that the criterion asks to follow from the model was put in by hand
(no installed field, no fitted or measured value, no supplied strength). If something is still supplied and
named, the cell is ⚠️. If the records do not show the criterion, the cell is ❌. If the result stands on an
installed or declared source, the cell stays 🚧 whatever the result. No agent of ours can accept a row; our job
is to produce the packet that leaves the maintainers nothing to object to.

## 2. What you produce: one row packet

A directory `openwave/row_NN_<name>/work/<date>_<party>_<short-name>/`. It maps one to one onto what OpenWave's
pull requests carry (task document, run, method note, small script, receipt, cell).

| File | Must contain | OpenWave counterpart |
| --- | --- | --- |
| `DESIGN.md` | The clause table: every clause of the simplest test verbatim; for each the observers, the records they hold, the observable as a function of records, the description to fit, the control, everything supplied. The `ESTABLISHED.md` entries and open premises it builds on | task document, part 1 |
| `TASK.md` + `PREREG.txt` | The candidate (law, memory, state, decoder, from `CANDIDATES.md`); gates as numbers; inputs; development and held-out seeds; controls; every value of nature used (or "none"). sha256 and time posted to the ledger **before** the run | `research/tasks/m10_<n>_task_details.md` (gates and inputs before the run) |
| code + `RECEIPT.json` | The run on the filed inputs, calling the pinned event function; every gate with its number and pass/fail; every outcome kept, failures and timeouts included | the run record at a pinned commit; `research/data/m10_<n>_receipt.json` |
| `CONTROLS.md` | (a) a changed law that fails on the row's own object; (b) a planted effect that the instrument detects; (c) the list of statements that are identities of the reader and carry nothing | the self-checks rule: a printed PASS must be able to fail |
| `headline.py` | Regenerates the headline number in under 60 seconds, exit code 0 or 1 | "a small script that regenerates the headline number" |
| `METHOD_NOTE.md` | Equations first; every term mapped to file and line; each result beside its gate; what was not computed; one **Scope** line | `dev_docs/METHOD_NOTE.md` |
| `CLAUSES.md` | Each clause: full / weaker / absent / negative, with the number | the cell's caveats |
| `CELL.md` | Cell text of at most 65 words, the expected mark (✅ ⚠️ ❌), every supplied piece named | the proposed `MODELS.md` cell |
| `MANIFEST.sha256` | sha256 of every file, written last | custody |

Then: post `SEAL row NN <manifest sha256> <path>` to the ledger, update the row's `STATUS.json`, add one line to
`ESTABLISHED.md`, and run `python3 openwave/tools/gate.py check`.

## 3. The packet fails if any of these is true

1. A field law, Hamiltonian, coherent step, coupling or clock was written down and evolved in place of the
   native process (it may only be the description that the records are fitted to).
2. The observable is read from the simulator's internal state and not from observers' records.
3. A gate holds for every law (an identity of the alphabet, the graph or the reader), or can only pass.
4. The changed-law control does not touch the object the row is about.
5. A negative was read from an instrument that misses a planted effect.
6. The terms were not posted before the run.
7. A measured value of nature entered where the criterion asks for a derivation; an anchor was not reported
   with what it then predicts.
8. A clause is borrowed from a run on another candidate.
9. The text claims more than the records show, or uses "accepted", "validated" or "complete" for our own work.
10. An `ESTABLISHED.md` entry was redone without a stated reason.

## 4. What is not a deliverable

A summary, a status report, a plan, a recount of rows, a run on an installed field presented as a row result,
a review by your own party presented as independent. These may be useful; they move no row.

## 5. Shared studies

Work on an open premise (P1 to P6 of `ESTABLISHED.md`) is produced the same way in `openwave/shared/<id>/`,
with the theorem to be proved or the execution to be run stated in `DESIGN.md`. A theorem is delivered as a
Lean file that builds, or a written proof with a reviewer's check, plus the statement of which premise it
discharges.
