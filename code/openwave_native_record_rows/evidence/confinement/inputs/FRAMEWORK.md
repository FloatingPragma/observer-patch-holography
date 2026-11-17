# FRAMEWORK: how a row gets accepted

Binding for every agent that works on the OpenWave rows of OPH. Start with [DELIVERABLE.md](DELIVERABLE.md) (one page: what
you produce, what acceptance is). Beside it:
[AGENTS.md](../../openwave-tuning/AGENTS.md) (the owner's premise and the words) and
[MANUAL.md](MANUAL.md) (the reasons) and [ACCEPTANCE.md](ACCEPTANCE.md) (what counts as acceptance and which evidence is produced). Set down on 2026-11-14 on the
owner's order: "Our plan will be to escalate every single row to accepted, in the sense that the result and
evidence would be fully accepted by OpenWave. We need to build a rigorous framework around this that allows
multiple agents to collaborate, and has quality gates." The owner corrects any line by editing this file.

## 1. The goal and what "accepted" means

**Mission.** Show, row by row, that the effective descriptions of physics (fields, particles, forces,
spacetime) emerge as fits to the records of observers running the native process of OPH, and fit them well
enough that OpenWave accepts each of the 31 rows.

The goal is a validated cell (✅) from the OpenWave maintainers on each of the 31 rows. A cell is their act.
What this framework controls is whether a row is **ready**: its evidence meets the terms the maintainers wrote
(`openwave-push/results_nov14/platform/m10_theory_canonical.md` section 8; `METHOD_NOTE.md`; `CONTRIBUTING.md`),
so that filing it is a formality of process.

| Word | Meaning here |
| --- | --- |
| ready | all gates G0 to G8 of section 4 are passed, each confirmed by a party of another group |
| partial-ready | the same with gate G8 passed at clause level only, every supplied piece named in the cell text; a milestone, never the end state |
| accepted | the maintainers have awarded the cell; recorded under gate G9 by the owner |

No file of this directory uses "accepted", "validated", "closed" or "complete" for anything else.

## 2. The order of work

**native process → observers' records → observable → fitted description → comparison with nature.**

Only work of class N (records of observers on the native process, apparatus supplied) or N+ (the same with an
added register, memory, preparation, routing, schedule or decoder: a new candidate) advances a row. Work of
class F (an installed field law or coherent step), S (a separately installed source) or C (a measured value
put in) is a target description: it says what the records would have to show. It is kept, linked from the
row, and passes no gate.

## 3. Candidates and the pin

OpenWave scores one executable model at a pinned commit. A changed law, memory, state or decoder is a new
candidate. Every candidate in use is registered in [CANDIDATES.md](CANDIDATES.md) with its four words (law,
memory, state, decoder), what an observer has to hold for it, its class and its pin state.

- A row is ready **on a candidate**, and only after the nodes it depends on in [DAG.md](DAG.md) are done. Rows that are ready on different candidates are not added up; the
  README reports the count per candidate.
- A new candidate is registered before its first run. Which candidate is pinned is the owner's decision.
- Identical `native.py` bytes do not make two protocols one candidate: the scheduler, the preparation, the
  memory, the read and reset operations and the decoder belong to the candidate.

## 4. The gates

A gate is `open`, `claimed` (the producer states it holds), `pass` (a party of another group has confirmed it;
G9: the owner) or `fail`. The state of every gate of a row is in its `STATUS.json`;
[../tools/gate.py](../tools/gate.py) refuses a pass that the producer's own group confirmed, a pass without
evidence, a pass before an earlier gate, and a changed manifest.

| Gate | Passed when | File in the row's work directory | Confirmed by |
| --- | --- | --- | --- |
| G0 design | Every clause of the row's simplest test is listed. For each: the observers, the records they hold, the observable as a function of the records, the description to be fitted, the control. The candidate is registered, everything supplied is listed, and the sentence of the theory that says the records should show the clause is quoted with file and line | `DESIGN.md` | design review by another group |
| G1 terms before the run | Gates with numbers, inputs, seeds (development apart from held-out) and controls are written and hashed, and the hash is posted to the ledger before the run | `TASK.md`, `PREREG.txt`, ledger entry | reviewer checks the order from the ledger |
| G2 execution | The run uses the filed inputs. The pinned event function is byte-identical to the pin and executes the events (a compressed sampler states the share of events executed by the pinned function and proves the compression). Every outcome is kept, failures and timeouts included | `RECEIPT.json`, logs | reviewer replays a subset |
| G3 controls | (a) a changed-law control fails on the row's own object for the stated reason; (b) a planted effect is detected by the instrument; (c) no gate is an identity of the alphabet, the graph or the decoder, and none can only pass; (d) every statistical statement names its interval, its selection and its power | `RECEIPT.json`, `CONTROLS.md` | reviewer runs own changed laws |
| G4 reproduction | A script regenerates the headline number in under 60 seconds; a replay is byte-identical; the manifest covers every file | `headline.py`, `MANIFEST.sha256` | reviewer runs both |
| G5 method note and cell | Method note in the platform's form (equations first, map from every term to code, scope line, what was not computed). Cell text of at most 65 words with the expected mark and every supplied piece | `METHOD_NOTE.md`, `CELL.md` | reviewer checks word by word |
| G6 other-party review | A full adversarial review, bound to the manifest, by a party of another group, with no blocking finding open. A review of a component is recorded and does not pass this gate | `reviews/<id>/REVIEW.md`, `REVIEW.json`, `MANIFEST.sha256` | the reviewing party |
| G7 one candidate | The packet runs on the candidate that the column is to be pinned to, or the pin-move statement (what changed in law, memory, state, decoder) is written and the owner has decided. No clause is borrowed from another candidate | `PIN.md` | integrator and owner |
| G8 full criterion | Every clause of the simplest test is carried in the full sense on that candidate; nothing that the criterion asks to follow from the model is supplied or fitted (rows 20, 26, 28: the strength). For partial-ready: a clause table with each clause full, weaker or absent | `CLAUSES.md` | reviewer of another group |
| G9 filed | Pin move, registered task with gates, in-platform run, method note and script are upstream and the maintainers have decided | owner's record in `STATUS.json` | the owner only |

Stages, derived by the checker: `open`, `designed`, `registered`, `executed`, `sealed` (G1 to G5 claimed),
`revise` (a review found blocking findings), `reviewed`, `ready`, `filed`. A sealed negative (the records do not
show the row's law) is evidence, recorded in the row with its number; the row's stage stays `open` and the
next design starts from it.

## 5. Parties, roles and independence

| Party | Group |
| --- | --- |
| Codex CRT (record team), Codex MFR (measurement-first root) | Codex |
| Fable-B and its workers | Fable |
| the external auditor and its lanes | Auditor |
| the owner | Owner |

- **Producer**: one party per row at a time; it writes the design, the task, the run and the packet.
- **Reviewer**: a party of another group. A check by a worker or lane of the producer's own group is useful
  and is labelled same-group; it confirms no gate.
- **Integrator**: keeps `README.md`, `CANDIDATES.md` and the checker in order and runs `gate.py check` after
  every change. The owner names the integrator; until then the party that wrote this file does it.
- **Owner**: decides pins, priorities and reassignments, files upstream, answers the maintainers in his own
  words. Nothing goes upstream without his click.

## 6. How agents work together

0. **Read the memory first.** Before any design, read [ESTABLISHED.md](ESTABLISHED.md) and the row's
   `THEOREMS.md`. Gate G0 is passed only if `DESIGN.md` cites the entries and open premises it builds on and
   states why it does not repeat an entry. Every sealed result, review and decision appends one line to
   ESTABLISHED.md.
1. **Claim.** Set `producer` and `claim` (party, machine-clock time) in the row's `STATUS.json` and post one
   line to the ledger. One producer per row. A claim with no ledger entry for 24 hours is free again.
2. **Ledger.** `openwave-tuning/CODEX_FABLE_SYNC.md` stays the message log, append-only. Entries that move a
   gate start with one of: `CLAIM row N`, `PREREG row N <sha256> <path>`, `SEAL row N <manifest sha256> <path>`,
   `REVIEW row N <bound manifest> <review manifest> <verdict>`, `DECISION`, `STOP`. Times come from the machine
   clock (`date`), never from an estimate.
3. **Review assignment.** A `SEAL` entry is answered by a party of another group with a review; reviews are
   reciprocal across groups.
4. **Findings travel.** A result that bears on another row is posted with the row numbers in its first line
   and added to that row's `README.md` under "From other rows".
5. **Shared studies** (an experiment that serves several rows) live in `shared/<id>/` and pass the same
   gates; each row that uses one names it in its `STATUS.json`.
6. **Disagreement.** A producer who disputes a blocking finding writes the dispute into the row's README under
   "Open disagreement" with both arguments. The owner decides.
7. **Stop.** A `STOP` entry of the owner halts every party. Each writes `INTERRUPTED.md` into unfinished work
   with the exact resume point. Nothing resumes without the owner.

## 7. Layout

```
openwave/
  README.md            the front door: status of the 31 rows, what the records show, decisions for the owner
  docs/
    MANUAL.md          the owner's manual: observations first
    ACCEPTANCE.md      what counts as OpenWave acceptance, and the evidence a row has to produce
    FRAMEWORK.md       this file: gates, roles, collaboration
    DAG.md             what each row needs before it can be accepted
    CANDIDATES.md      the candidate registry
    MATH.md            the mathematics, with the file behind each statement
    INVENTORY.md       what the two older directories hold
  tools/gate.py        the checker: check, render, show
  shared/              studies that serve several rows (the foundations of DAG.md)
  row_NN_<name>/
    README.md          the row: test, clauses, what the records show, target descriptions, what is owed, evidence
    STATUS.json        gates, packets, reviews, blocking findings (read by the checker)
    work/<date>_<party>_<name>/      new packets, one directory each, files as in section 4
    reviews/<date>_<party>_<name>/   reviews of this row's packets
```


The evidence of 2026-11-13 and 2026-11-14 stays where it is, under `openwave-push/` and `openwave-tuning/`:
its manifests pin those paths. Rows link to it. New work goes into the row directories. Data above 100 MB
stays out of the repository and is named by path and hash.

## 8. Rules learned from the reviews of 2026-11-14

Each was found in a sealed packet that day. A reviewer checks every one; any of them blocks G3 or G8.

1. **A reader's identity is no result.** A statement that holds for arbitrary laws (integer degree of an
   integer tape, the period of a counter, half-odd labels of a lifted rotation, conservation modulo 6) is
   marked as an identity and carries no clause.
2. **A gate must be able to fail.** No literal `False` or `True`, no comparison of a call with itself, no
   state found by solving the equation that the gate then tests.
3. **A control must touch the row's object.** A changed law that differs from the native one only away from
   the object tested (for example only on doubly loaded ports when the object is one unit) is no control.
4. **A detector is validated before a null is read.** A gate that returns false on a planted effect makes
   its negative uninformative.
5. **The pinned function runs.** A sampler that draws from a compressed chain says so, proves the
   compression and executes the pinned function on a stated share.
6. **An installed step is named.** In W = S·R the shift S is the native event with a memory law and the
   reflection R is installed: every phase on that candidate exists because of R. Leaving a register unread
   selects no R.
7. **A calibration is an input.** An anchor (511 keV) fixes a unit; it is reported with what it then
   predicts without further input, and a miss is a result. A frequency set by a supplied number is an input.
8. **Passive registers respond to nothing.** A register that the load law does not read cannot carry a
   force, an energy or a motion; its statistics are those of its cost table.
9. **Means superpose.** The mean law of `native` is linear, so a force, a binding or a screening read from
   mean records gives the single-unit answer. Such rows need statistics beyond the mean or a process whose
   events read records.
10. **Three things stay apart**: what the theory's papers derive, what an executed construction shows, what
    the platform counts. A conditional theorem is quoted with its premises as things observers must measure.
11. **Protocol and class.** The pinned candidate is the admitted controller (two events, total-two loads). A
    continuing schedule, another preparation or a longer memory is N+ and belongs in the pin move.
12. **Time and order.** Terms before runs are shown by a ledger entry, not by file times.

## 9. The owner's actions

Decide the pinned candidate and every pin move; decide whether a reflected, reused choice register is part
of the model's memory law; set priorities and reassign rows; post the premise comment, the corrections and
the questions to the maintainers; file. Each is recorded as a `DECISION` entry in the ledger.
