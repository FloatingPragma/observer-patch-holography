# Progress tracking

OPH tracks progress in two places with two different questions.

| | Model and law evaluation | Global completion ledger |
| --- | --- | --- |
| Question | How far does one candidate law satisfy the requirements on the executable model? | Which observed properties of the Universe does the current material demonstrate, and which theorem premises must the model confirm? |
| Unit | One criterion per candidate law | One observable or one premise |
| Lives in | [`model/evaluation/`](../../model/evaluation/README.md): [scorecard](../../model/evaluation/SCORECARD.md), criteria files, each candidate's `EVALUATION.json` | [`OBSERVABLES.json`](OBSERVABLES.json) in this folder |
| Updated by | `python -m model evaluate --write` | Editing the ledger; `python tools/check_observables.py` validates it |
| Shown at | The scorecard in the repository | [status.floatingpragma.ai](https://status.floatingpragma.ai/) |

The scorecard grades laws. The ledger grades the program: it combines the finite theorems, the papers
and the executable results of the main candidate (Interval Echo) into one status per observable. A
scorecard result enters the ledger as evidence through the `model` links of the rows it supports.

## The ledger

[`OBSERVABLES.json`](OBSERVABLES.json) has three keys: `schema`, `statuses` (the binding definition of
each status) and `rows`. Every row has the same fields:

| Field | Content |
| --- | --- |
| `id` | Stable identifier. Observable ids are short slugs; premise ids are `PR-nn`. Never rename or reuse an id; the status page links to rows by id. |
| `kind` | `observable` or `premise` |
| `name`, `category` | Display name and group |
| `status` | `Demonstrated`, `Conditional`, `Partial` or `Open` |
| `what` | One sentence: what the observable is, or what the premise assumes |
| `why` | One sentence: why it matters, or which results rely on the premise |
| `premises` | Premise ids the row waits on |
| `lean`, `model`, `papers` | Repository paths: Lean files; `model/` and `code/` evidence; papers and other documents |

Premise rows also carry `statement` (the exact assumption, which `model/evaluation` reads) and
`premise_kind`: `declared` (an architecture choice of the construction), `derivable` (a physical premise
the model must produce), `calibration` (a value taken from data) or `mathematics` (an established result).

The observables are the target list of the program. The premises are the assumptions of OPH's finite
theorems beyond the axioms A1 to A3, the declared basis and established mathematics. A premise that
states the same thing as an observable is not listed separately.

## Statuses

- **Demonstrated**: the current material shows the whole observable. Either a theorem or derivation chain
  (Lean or paper proof) whose premises need no confirmation from the executable model, or a committed
  executable result of the model whose observer records exhibit the observable. A premise is Demonstrated
  when the executable model satisfies it, or when it is established mathematics.
- **Conditional**: proved, with at least one premise that needs confirmation from the executable model
  that the current material does not supply. The row lists those premises.
- **Partial**: part of the observable is shown: some clauses, a restricted domain or preparation, a
  candidate other than the main one, or a retrospective comparison.
- **Open**: no result in the current material.

A status rests on a committed file that states the result: a theorem, a test, a receipt. A README sentence
or a file name is not enough. Agreement with values used for calibration is at most Partial. Particle mass
values are never stated in the ledger.

Judging a premise:

- The `what` of a premise states the assumption and stays fixed when evidence arrives. The status reports
  how much of the `statement` is shown. `what` and `why` never describe an experiment or its outcome; the
  evidence goes in the links.
- A `derivable` premise is Demonstrated when the law of the main candidate produces the `statement` as
  written. A result that holds on a supplied apparatus, routing, schedule, cadence or preparation is
  Partial, and so is a finite instance of a general statement.
- A `declared` premise is Demonstrated when the construction declares the choice and the executable
  implements it.
- When a premise loses or gains Demonstrated, re-judge the rows that list it. A Demonstrated row lists no
  premise that is not Demonstrated.

## Completion score

Each observable scores Demonstrated 1, Partial 1/4, Open 0. A Conditional observable scores 1/2 plus half
the mean score of the premises it waits on, with premises scored on the same scale. The completion score
is the mean over all observables. Demonstrating a premise raises every observable that waits on it. The
status page and `tools/check_observables.py` compute the same number.

## Updating the ledger

When a result lands (a Lean theorem, a model experiment, a paper section):

1. Find the rows it affects: search the ledger for the observable, and for premise ids the result
   discharges.
2. Set the status by the rules above. Add the new evidence paths to `lean`, `model` or `papers`.
3. When a premise becomes Demonstrated, re-judge every Conditional row that waits on it.
4. When a new theorem uses an assumption that has no row, add a premise row with the next free `PR-nn`.
5. Keep `what` and `why` to one sentence each. Do not add status narration or long paragraphs.
6. Run `python tools/check_observables.py` and commit. The status page reads the file from `main`.

A row is removed when its subject leaves the research program, for example a test that only distinguished
retired laws, and a duplicate is merged into the row it repeats. Removed ids are not reused.

A daily cloud routine reads the commits on `main` since the last ledger change, re-judges the rows they
affect by these rules, runs the check and pushes the update.
