# AGENTS.md: binding rules for every agent working on the OpenWave rows

Applies to every agent and every lane that reads or writes under `openwave-tuning/` or `openwave-push/`: Codex,
the external auditor with its lanes, the Fable sessions and their workers. Set down on 2026-11-14, 20:20 (+07),
on the owner's order. The reasons, the quotations and the worked shapes are in [MANUAL.md](MANUAL.md). Read
both before any row work. Where an older file of this directory says otherwise, this file holds.

## 0. Mission (owner, 2026-11-14)

Show, row by row, that the effective descriptions of physics (fields, particles, forces, spacetime) emerge as
fits to the records of observers running the native process of OPH, and fit them well enough that OpenWave
accepts each of the 31 rows. Workspace and gates: [../openwave/README.md](../openwave/README.md).

**What you produce and what OpenWave acceptance is: [../openwave/docs/DELIVERABLE.md](../openwave/docs/DELIVERABLE.md)** (one page).
**Before any work: read [../openwave/docs/ESTABLISHED.md](../openwave/docs/ESTABLISHED.md)** (established results and open
premises). Do not repeat an entry without stating why; append a line for every sealed result.

## 1. The premise, in the owner's words (2026-11-14)

> "our theory is observers first and reality emerges from consensus, not from the carrier"

> "What comes out are measurements that are described by effective descriptions of observers. Not the fields
> themselves."

- OPH starts from observers. Observers compare records and reach agreement. Reality is what they agree on.
- The output of the model is measurements: records that observers take and that are compatible between them.
- A field (Maxwell, Dirac, Klein-Gordon, a colour field, a metric) is the effective description that
  observers give of those measurements. It is no object of the model and no object on the carrier.
- "Carrier" is the name of the network on which the observers of the executable model live and run the
  repair-and-memory law. Physics is never said to come out of a carrier, a substrate or a medium.

## 2. How a row is tested

1. Run the native repair-and-memory process (`reverse-engineering-reality/model/observer/native.py::native`
   at the platform's pin, or a candidate law that the owner has decided to pin).
2. Let observers take records. State which observers, which records, and everything that is supplied to them.
3. Compute the row's observable from the records alone.
4. Fit the effective description to the records. The field named in the criterion is that fit.
5. Run a control that can fail: the same observers on a changed law give a different record.

A row has the reversed shape when a field law, Hamiltonian, action, coupling, clock, code or coherent step is
written down, installed on the cells and then read out. Such a construction is a premise. The theory
repository states this itself (`reverse-engineering-reality/docs/OPH_EXECUTABLE_HYPOTHESIS.md`):

> "Supplied coherence, protected codes, radial directions, fields, clocks or couplings remain premises until
> realized by this process."

> "One source process must meet all parts of this claim. It is not enough for different supplied actions to
> reproduce its separate words."

- An installed field law may serve as the target description: the form with which the fit of step 4 is
  compared, and a calculator for what the records should look like. It is never run in place of the
  observers' process and then reported as the row.
- A calibration is never an end state. The data are source-derived (owner, 16:20).
- One process serves all rows. Results on different installed sources are not added up.

## 3. How results are counted and reported

- The count that matters is OpenWave's. Official today: 0 validated, 0 partial, 31 untested. Nothing is filed.
  Only the maintainers award a cell.
- OpenWave's rule for this column is OPH's own rule, taken from the hypothesis page: a row moves when its
  observable is computed on the process that the native observers run; a declared source leaves the row
  untested.
- Every status line gives the platform reading first, then the group's internal count, labelled as the
  group's own. The words "validated", "accepted", "counted", "closed" and "complete" are never used for a
  row without the name of who decided and the class of the construction.
- Every row result names its class:

  | Class | Construction |
  | --- | --- |
  | N | records of observers on the native process, with the apparatus supplied |
  | N+ | the same with an added register, memory or decoder (a new candidate; it needs a pin move) |
  | F | a field law or coherent step installed on the carrier's cells |
  | S | a separately installed source (a colour chain, a ring field, a species source) |
  | C | any of the above with a measured value put in |

- N and N+ are results of OPH in the sense of section 1. F, S and C are premises with a tested property.
- The ten rows of the group's internal count, classed by the audit of 2026-11-14
  ([../openwave-push/Push-results-Nov14.md](../openwave-push/Push-results-Nov14.md)): N row 1; N+ row 7;
  F rows 2, 4, 19, 29 (rows 2 and 19 with the 511 keV anchor); F with the measured α put in, row 20;
  S rows 12 and 13. Row 30 stands on a coherent walk whose committed version is the native chain: N+ under
  the retention-schedule reading (the native event with the choice register left unread,
  `../openwave-push/onecandidate/DP1_retention_schedule/`), with the redraw operator declared and the mass
  angle supplied; F otherwise.
- The campaign that redoes the rows measurement-first, its split between the parties and the status of every
  row are in [MEASUREMENT_FIRST.md](MEASUREMENT_FIRST.md).

## 4. Words

- Write: "the records of observers on the native process are described by ...".
- Do not write: "the field on the carrier", "the carrier produces", "derived from the carrier", "the
  substrate", "the medium". When the platform's own text is quoted, the quotation keeps its words.
- Where the records do not show what a row asks: state what they show, with the number, and name the
  statement or construction that would close the difference. Never "no mechanism". Never "structural".
- Statements of fact. Never the word "honest". No em-dashes. No "not X but Y".

## 5. Standing limit

Nothing goes upstream without the owner's own click: no fork, push, pull request, issue, comment or discussion
post. Answers to the maintainers, the premise comment and corrections to their pages are written and posted
by the owner.

## 6. Provenance

Section 1 quotes the owner. Sections 2 to 4 are the consequences drawn from his words and from the hypothesis
page by the second Fable session, written on his order of 20:15; the class letters of section 3 are a
labelling convention introduced here. The owner corrects any line by editing this file.
