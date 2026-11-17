# MANUAL: testing the OpenWave rows observer-first

For every agent that works on the 31 OpenWave rows of OPH (column M10). Written on 2026-11-14, 20:20 (+07), on
the owner's order, after the audit of the day's results
([../openwave-push/Push-results-Nov14.md](../openwave-push/Push-results-Nov14.md)). The binding rules are in
[AGENTS.md](AGENTS.md); this file gives the reasons, the quotations and the worked shapes.

## 1. The premise

The owner, 2026-11-14:

> "our theory is observers first and reality emerges from consensus, not from the carrier"

> "What comes out are measurements that are described by effective descriptions of observers. Not the fields
> themselves."

Three statements follow and are used throughout this manual.

1. The primitive objects are observers and their records.
2. What the model outputs is measurements: records that observers hold and that agree where they overlap.
3. A field is a description. It summarises the agreed measurements. It is no ingredient of the model.

## 2. The same premise in the theory repository

`reverse-engineering-reality/docs/OPH_EXECUTABLE_HYPOTHESIS.md` (byte-identical at the platform's pin
`3c4e1e23fc32` and at the published head `fca399cf0ebe`):

> "Observers recover and realize their world's specification inside that same history, including the
> realizing hardware. Space, time and matter are descriptions of the compatible observations. The physical
> target is indistinguishability of native observers' experiences from those in our Universe, including its
> fields, gravity, particle properties and masses."

> "The executable-model hypothesis is that a simple, fully specified local repair and memory law on the
> proposed carrier architecture admits continuing observers whose compatible records and interventions have
> the observed 3+1-dimensional quantum-field description and a common classical-gravity limit."

> "One source process must meet all parts of this claim. It is not enough for different supplied actions to
> reproduce its separate words."

> "Supplied coherence, protected codes, radial directions, fields, clocks or couplings remain premises until
> realized by this process."

> "Changing the law, memory, state or decoder defines a new candidate and keeps its empirical ancestry
> visible."

The page lists five obligations (EMH-1 to EMH-5). Each reads "not established". EMH-1 asks for the
observer-accessible algebras, their clocks and the energy normalisation; EMH-2 for the effective quantum
field description of the same process.

## 3. What OpenWave asks, and why it is the same rule

The maintainers wrote the scoring rule of the M10 column from the application and from the hypothesis page
(copy of their page: `../openwave-push/results_nov14/platform/m10_theory_canonical.md`; admission reply:
`auditor/upstream/discussion_712.json`).

| Their sentence | Where | Reading |
| --- | --- | --- |
| "The substrate is not a field and not a dynamical medium" | § 2 | They do not take the carrier for a medium. |
| "The native law is a local repair and memory law, run by the native observers. There is no Lagrangian and no field equation." | § 3 | The process is the observers' process. |
| "The row constructions are a different object. They are finite quantum Hamiltonians on declared cubes, chains and two-cell blocks, among them a Wilson plaquette action, a lattice Dirac Hamiltonian and a Maxwell complex" | § 3 | Installed field laws are set apart from the model. |
| "by the author's account each runs on its own declared finite source and none on the twelve-port screen" | § 7 | The state of the application's packages. |
| "A row whose construction runs on a source that is declared, not derived from the carrier, stays at 🚧." | § 8 | An installed field leaves the row untested. |
| "What moves a row: Its observable computed on the carrier's own process in-platform. With something still supplied, that is a case for ⚠️ with the caveat in the cell. With the source derived from the carrier, it meets the bar for ✅" | § 8 | The observable is a measurement on the observers' process. Supplied apparatus gives ⚠️. |
| "derive the encodings in-platform. [...] When those land, the exponents are measured rather than encoded, and the row flips." | issue 482, linked from § 8 | ✅ asks that the description is measured from the records, where today it is written in. |
| "That is the author's own rule on the hypothesis page, adopted here as an admission term." | § 1 | One candidate at a time is OPH's own rule. |

"Derived from the carrier" in their text means: obtained from the records of the process that the observers
run, in place of written down. It asks for no field object on the screen. Their page contains neither the
word consensus nor the premise list; it says the premise list is taken from the comment that the application
announced. The owner posts that comment and the corrections they invited ("Where we got your framework
wrong, fix it").

## 4. The two shapes

```
observer-first      native process  ->  observers' records  ->  observable  ->  fitted description  ->  nature
(the theory)        (runs)              (taken, compared)       (from records)  (the "field")

reversed            field law written down  ->  evolved  ->  read out  ->  nature
(installed field)   (the "field" is put in first)
```

The reversed shape shows a property of the installed field law. It shows nothing about what the observers'
records are. Under OPH's own page it is a premise; under the platform's rule it is a declared source.

## 5. Worked examples from 2026-11-14

**Row 1, charge quantization (class N).** The native process runs on one carrier with one or two loads. A
record-only observer holds the current tape. Observable: every event current is an integer (5,400 of 5,400
branches), the record winds G − n round each port, one load is never divided (720 of 720). Supplied: the
preparation, the reference G, the circle encoding, the counters, the closing sweep, the pointer. Weak point
found by the audit: the degree gates pass on any integer tape, so the row needs a control on a changed law
inside its manifest.

**Row 7, double cover (class N+).** The native process is read through an added area counter modulo 20. Ten
face loops multiply the first character of the record by −1 (defect 3.25e-18 over 240 executed columns);
twenty are the identity. Weak point: both returns hold for 401 of 401 arbitrary laws, so they are identities
of the readout. The law enters through the executed loop operation and the size of the frame moment. The
counter is a new decoder and needs a pin move.

**Row 20, Coulomb (both shapes exist).** Counted by the group on an installed quadratic field with the
measured α put in (class F with C): the 1/r form is a theorem of the installed law. On the native process
(class N, sealed packet, not filed): the mean tape of a prepared pair over histories equals the lattice
Coulomb field of the routing, and the cell-averaged potential follows a + b/r between 3 and 8 cells. The
second result has the theory's shape.

**Row 29, Maxwell waves (reversed shape).** Counted by the group on the vacuum limit of an installed
quadratic field whose coefficients are chosen so that its speed equals the matter speed. On the native
process no sealed run tests the row; the table entry reads that loop currents on the tape diffuse. A test in
the theory's shape: a disturbance of the native process, records of separated observers, and the law of
propagation and polarisation fitted to those records.

## 6. What the records of the native process show today

From the row audits of 2026-11-14 (`../openwave-push/results_nov14/audit/ROW_XX.md`). "Records" means records
on `model/observer/native.py::native` with the apparatus named in each packet.

| Group | Rows and what the records show |
| --- | --- |
| The records carry the row or a clause of it | 1 (integer winding; a unit never splits). 7 (returns at π and 2π through the area counter). 4 (the unit is never lost or divided in 32,424,192 events; its position law is not stationary). 6 (half-integer levels through the area counter; the sector is the reader's choice). 8 (a contact pair is removed in 360 of 360 cases; nothing radiates). 20 (the mean record of a prepared pair equals the lattice Coulomb field). 28 (no event changes the port count: 0 of 4,680 branches; sizes 12, 42 and 162 are all admitted). 31 (with a supplied observer weight the record statistic has discrete levels as decay rates). |
| The measurement is executed and the row's law is absent from the records | 3 (every eigenvalue of the event matrix is real; no oscillation). 5 (g = 0.787 ± 0.012). 9 (no object decays: 0 of 110,600 cases on candidate J). 10 (every non-uniform mode decays monotonically). 12, 13, 16, 17 (stacked units separate: survival 5/6 per event). 19 (an excess and deficit pair annihilates; a closed tape loop is stable and has no energy defined). 25 (the law is covariant under all 120 automorphisms of the carrier, 60 of them improper). 26 (displacement toward 72 held units negative in 18 of 18 entries). 27 (hop-rate factor exactly 1 on 1,458 of 1,458 ports). 30 (with a record at every step the population law is I − L/12: real spectrum, no phase per step). |
| The present law holds no object for the row | 11, 14, 15, 18 (one conserved label; no energy functional). 22 (no energy is attached to the record). 23 (no coupling is read from a force at a scale; the mean record gives a form factor of the lattice). 24 (no heavier charged species). |
| No record experiment on the native process appears in the row audits | 2 (a packet on candidate J exists, class N+). 21 (every readout is an ideal readout of a state vector). 29 (a table entry, no sealed run). |

For each row of the second and third group the task is stated in the same form: which observers, which
records and which comparison between observers would show the row's law, and what the records show when that
experiment runs. Which records observers hold and compare is part of the model and is stated before testing
(the hypothesis page: "complete C/E record instruments, memory law, scheduler, preparations and available
controls must be stated before testing").

## 7. Procedure for a row

Before the run, in a file with a hash:
1. The process: the law, the carrier or screen, the scheduler, the preparation.
2. The observers: who records what, when, and what they compare.
3. Everything supplied to the observers (apparatus, encoding, reference, readout).
4. The observable as a function of the records.
5. The description to be fitted and the gate, with numbers.
6. The control that fails: a changed law, a shuffled history, a planted effect that must be detected.
7. Every value taken from a measurement of nature, or the statement that there is none.

After the run:
1. The record statistics and the fit, with the control's result.
2. The class (AGENTS.md, section 3) and the list of what was supplied.
3. Where the row's law is absent from the records: what the records show, with the number, and the
   statement or construction that would close the difference.
4. A review by a party other than the producer's, bound to the manifest.

## 8. What to do with an installed-field result

The day's results on installed fields stay on file as premises with tested properties. Each has two uses.

1. **Target description.** It states the law that the fit to the records is compared with (for row 29: two
   transverse polarisations at one speed; for row 30: E² = m² + p²).
2. **Calculator.** It predicts the record statistics that the observers should find if the description
   holds, including the size of the effect at the reachable block sizes.

The row is reported as a result of OPH once the records of observers on the native process are shown to
follow that description, with a control that fails. Until then its line reads: class F (or S, C), the tested
property, and the record experiment that is owed.

## 9. Status line

```
Row 29. Platform: official untested; if filed as it stands: strict untested, most favourable partial.
Group's internal count: counted, class F (installed quadratic field, coefficients chosen).
Records on the native process: no experiment on file. Owed: propagation records of separated observers.
```

The platform reading comes first. The internal count is labelled as the group's own and carries the class.

## 10. Words

| Write | Do not write |
| --- | --- |
| the records of observers on the native process are described by a Maxwell field | the Maxwell field on the carrier |
| observers' compatible records follow E² = m² + p² | the carrier produces a Klein-Gordon wave |
| the description fitted to the records | the field derived from the substrate |
| the group's internal count, class F | validated, accepted |
| the records show a hop-rate factor of 1; owed: a record experiment with compared clocks | gravity has no mechanism; the stop is structural |

## 11. Where things are

| What | Path |
| --- | --- |
| Status of all 31 rows with evidence, audit and submission steps | `../openwave-push/Push-results-Nov14.md` |
| One audit page per row | `../openwave-push/results_nov14/audit/ROW_XX.md` |
| The platform's files as read at `5d686cb049b8` | `../openwave-push/results_nov14/platform/` |
| The theory's hypothesis page | `../reverse-engineering-reality/docs/OPH_EXECUTABLE_HYPOTHESIS.md` |
| The native law | `../reverse-engineering-reality/model/observer/native.py` |
| The group's internal count | `row_completion/INTERNAL_VALIDATION.json` |
| Coordination between the parties | `CODEX_FABLE_SYNC.md`, `sync/` |

## 12. Provenance

Section 1 quotes the owner. Sections 2 and 3 quote the theory repository and the platform. Sections 4 to 10
are consequences drawn by the second Fable session on the owner's order of 20:15, with the row facts taken
from its audit of the same day; that audit was made by workers of one of the parties and the owner has
announced a second one. The owner corrects any line by editing this file.
