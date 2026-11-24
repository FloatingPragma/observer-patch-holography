# OpenWave Submission Guidelines

These guidelines bind every contributor and agent who works toward a row of the OpenWave comparison table from this
repository: research aimed at a row, a row packet, an audit, and any pull request to [OpenWave][ow], where OPH is
admitted as column M10. Read the whole page before such work; the repository [agent guidance](../AGENTS.md) makes it
binding. Where this page asks for more than the platform, it holds for our filings. Where the platform's published
rules ask for more, the platform holds. Upstream links point at one revision so that every quotation stays
checkable; the upstream default branch decides.

The central failure these rules guard against is treating reproducible measurements through the native process as
evidence that the process produced the claimed physical property. Arithmetic and execution can be verified while
that inference goes unchallenged, and a disclosed assumption does not make an unsupported property eligible for
partial credit.

**Do not claim that a law produces a property when the evidence only shows a consequence of something supplied by
construction.**

Every claimed scientific test answers one question: **What possible run outcome would falsify each claim?** A run is
a test of a claim only with a concrete answer: an outcome representable in the stated test domain, a way to observe
it, and a rejection condition fixed before the run. The claimed law may predict that this outcome never occurs; the
representation must not exclude it by definition. If every record the apparatus can represent passes by
construction, the run supplies no test of the claimed emergence. A corrupted file or an out-of-domain signal is no
substitute. Distinguish an exact counterexample from statistical rejection and from an inconclusive timeout.

## 1. Authority

| Party | Decides |
| --- | --- |
| The owner, author of the M10 column | the nominated candidate, any proposed pin move, and the release of each audited version |
| Codex | the audit verdict on every outgoing version, for each row and for the package |
| OpenWave maintainers | the adoption of a pin and the mark of every cell, one cell at a time, at review |

- The platform's claims are the cells of [MODELS.md][models] and the scripts behind them. In the words of the
  [column's consumption rules][canon8]: "Statuses on the author's own site or papers are the author's estimates."
  A grade in this repository's [scorecard](../model/evaluation/SCORECARD.md), such as Demonstrated or Partial, is a
  local assessment of a local candidate law. It is no OpenWave mark and does not make a row submittable.
- Codex audits every upstream submission against the checklist of section 10, on the exact outgoing version. When
  Codex produced the work, a separate Codex reviewer within the project performs the independent audit. Codex marks
  a row submittable only when it expects full credit (✅) on that candidate. An unchecked, failed or unresolved
  applicable item blocks the submission.
- Nothing in this repository authorizes an agent to post upstream, open or update an upstream pull request, propose
  or move the pin, or record a cell award. Research authorization covers local work, and the owner's authorization
  of one version does not extend to another.
- Contributors bring M10 work to this repository, and the owner files it upstream. A collaborator's draft pull
  request to OpenWave is merged here and closed upstream until the owner files the result.
- The upstream [M10 agent orientation][orientation3] states the working rules of the column: "These bind the
  author's agents exactly as they bind the maintainers'." Its first rule is to read the open pull request's
  conversation and review state at the start of every session, since maintainer decisions land there.
- The platform's [AI hygiene stance][hygiene] governs every claim: "A model's output is a draft or a hypothesis,
  never a result. It becomes a result only when it is verified by something that is not a language model: a
  derivation checked by hand, a script anyone can run, a measurement on the lattice, or the explicit confirmation of
  the human who holds the relevant authority."

## 2. Start from the evidence

Before work on a row:

1. Record the default-branch heads of [OpenWave][ow] and of this repository, and read the open M10 threads at that
   head. A verdict standard is written from the head of the day, never from notes taken at an older one.
2. Read the row's literal clause in the [simplest-test table][simplest], its regime, the M10 cell, and every
   maintainer decision that applies to it. The local copy in
   [openwave_criteria.json](../model/evaluation/openwave_criteria.json) names the upstream revision it was taken
   from; the upstream text at the recorded head decides.
3. Read the candidate's folder under [model/candidate_laws](../model/candidate_laws/CANDIDATE.md): its
   specification, its `candidate.py`, its `EVALUATION.json` and its OpenWave evidence in the `openwave/` subfolder.
   Identify the existing results, their hypotheses and the precise unresolved obligation before designing work.
4. Read the [M10 specification of record][canonical], with the pinned candidate in its § 1 and the consumption
   rules in its § 8, and the [M10.2 maintainer review][review]. Section 8 of this page carries that review's row
   decisions.
5. Search this repository's papers, Lean library, [premise register](../tracking/premise_register.json) and code
   for the same question (caution 7.15). Record the paths, the commit and whether each confirms, extends or
   contradicts the planned work.

Keep sealed evidence, adopted pins and other workers' uncommitted files intact. Compare them with the sources in use
and document each difference instead of replacing the tested inputs.

## 3. The platform's bar

The bar is one for every column. [MODELS.md][bar] states it: "A cell reports what the model's own process returns.
A property fixed by what the model declares before any run (an assigned integer, a calibrated value, a symmetry of
the declared field or cell) is recorded in the cell's text and earns no icon: a declared or installed source is
inherited, and the row moves when its observable is computed on the model's own process in-platform."

| Mark | When the platform awards it |
| --- | --- |
| ✅ validated | "the property comes out of that process and meets the row's test" |
| ⚠️ partial | the property comes out of the process with a piece supplied beside it, or with part of the row open; a value the model is calibrated to earns nothing by itself |
| 🚧 | nothing runnable returns the property, whatever the model declares; a declared or installed source is inherited |
| ❌ | a documented negative on the model's own process, with the script that produced it |

The [consumption rules][canon8] of the M10 column, taken from the [admission reply][admission]:

- A row whose construction runs on a source that is declared, and not derived from the carrier, stays at 🚧.
- A row moves when its observable is computed on the carrier's own process in-platform. With a piece supplied
  beside it, that is a case for ⚠️ with the caveat in the cell. With the source derived from the carrier, it meets
  the bar for ✅.
- A calibrated source is a calibration. It earns no credit for the values it was fitted to, and its held-out
  failures go in the cell as results.
- The Coulomb, Newton and cosmology rows (20, 26 and 28) credit a strength only when it follows from the model's own
  mechanism, never fitted or supplied.
- The maintainers decide each cell at review, one at a time, from the artifact ([Gate E][gatee]).

The maintainers' reply on the M10.2 thread, which the [maintainer review][review] links, reads this rule for the
pieces of a pinned candidate:

- "Cells, law, memory, schedule and routing, as pinned, are the carrier's own process." A piece of the pinned
  carrier is no declared source by itself.
- "A row stays 🚧 when what it names comes from a declared piece and not from the run." The verdict of section 6 is
  then `declared`.
- A row computed on the process with a piece supplied beside it (a decoder with an anchor, a field, an action, a
  clock, a coupling, a trap, a declared centre) is the case for ⚠️. The bar for ✅ is met when nothing is supplied.
- "Where the row names a measured number, the number must also agree." A number clause is scored apart, with its
  tolerance and with what was known beforehand.

The same reply makes prediction or post-fit the deciding rule, row by row ([onboarding § 1.1][onboarding]). A law
developed until the rows pass has had the targets in view. For every row, record what was known beforehand: whether
the row's target was consulted while the law was being fixed, and which choices of the law answer to it. Keep the
list of laws set aside, because the [parameter-count test][paramcount] counts a choice as free when the framework
admits an inequivalent alternative. A particle mass, a coupling or another row target placed in the law is a
calibration of that row. A claim that the law derives such a value needs the law frozen with the value as an output
and a held-out comparison.

## 4. What may be filed

- A row submission requests full credit (✅) on the complete official criterion. Codex must expect ✅ for the entire
  criterion and support that judgment with evidence for every clause. A scoped, partial or registration-only result
  stays on **HOLD** in the research record, together with prospective designs and historical reproductions.
- One focused pull request may nominate a corrected candidate and present related rows on it, with each row's scope,
  evidence and requested decision stated separately. In the M10.2 filing the author brought a new pinned candidate
  and three rows together; the maintainers adopted the pin and decided each row on its own evidence.
- The owner may file any subset of rows that individually meet the full-credit rule. Anticipated acceptance of
  further rows is described as a research aim. A claim covering all 31 rows needs evidence meeting all 31 criteria
  on the named candidate.
- The column scores one pinned candidate at a time ([specification of record, § 1][canon1]). A row computed on a
  candidate law of this repository is a row on a different candidate unless that law is the pin. Its filing carries
  a pin move: the current and proposed commits, the earlier failure or limitation that motivated the correction, the
  change in law, state, memory or reader, and which results were reproduced on the corrected source. The owner
  decides whether to propose it.
- A pin move is a pull request that states the change of candidate ([§ 1][canon1]). It adds a declaration of the new candidate in
  the form of § 1.1 of the specification of record: carrier and addressing, update law, memory, preparation, each
  observable with its reading, calibration, and the inferences excluded. It states the change from the ancestor and
  keeps the ancestor's failures on record. After the move the cells decided on the old pin describe it as ancestry
  until a task brings each row onto the new candidate. The owner told the maintainers on the M10.2 thread that the
  law is frozen before it is filed and that the pin moves once, so the pinned commit is a frozen law, never a law in
  development.
- A pinned commit stays reachable. The maintainers asked on the M10.2 review for the pinned commit to be tagged,
  because a commit that only a side branch holds can disappear when the branch is deleted. Tag every pinned
  candidate commit in this repository before filing, never delete a branch that holds a pin, and point every source
  link of a filing at a commit reachable from the default branch or a tag.
- Each law under `model/candidate_laws/` is its own candidate. A changed domain, schedule, controller, memory,
  preparation, decoder or outcome alphabet makes a new candidate with its ancestry, even with unchanged endpoint
  code. No result transfers between candidates by ancestry.

| Request | What the submission contains |
| --- | --- |
| Correct the law and move the pin | The current and proposed executable commits, the earlier failure or limitation that motivated the correction, the change in law, state, memory or reader, and the results reproduced on the corrected source. The original evidence and failures stay accessible. |
| Full credit on a row | Each clause of the official criterion matched to its preparation, observable, units, controls, uncertainty and outcome. Every supplied or calibrated ingredient disclosed. The method note, the complete record or immutable manifest, and a small runnable check. The reason Codex expects ✅ for the whole criterion. |
| A link to the public simulator | The accessible URL, the view or mode, how to find the reading, and which law and observable that view implements. A shared link may serve several rows when its relation to each is explained. A demonstration on a related implementation states its scope and cites the separate row receipts. An exact continuum claim needs a derivation on the named law with its hypotheses and a verified connection to the simulator. |
| A prospective run | Inputs, gates and controls frozen before any confirmatory acquisition. Prior knowledge disclosed and historical dates preserved. The registration is kept distinct from approval of a completed row. |
| An administrative change | The exact edit, the affected references verified, and the requested action. |

## 5. Observers and records

The world of OPH is what observers read from their records. An observer is a local part of the structure with its
own history of interactions and admitted records. Measurement, consensus and the synchronization of clocks come
after the records and are built from them. A field, a particle, a distance, a common time and a wave function are
effective descriptions that observers agree to use. The state of a simulation is a mathematical encoding that
computes what observers read, and its evolution parameter is no time that anyone measures. What is compared with
nature is the readout of observers from their records.

1. **State a claim on records.** Name the observer, the record it holds, the reading computed from that record and
   the probability the law gives. A sentence about a state, a field or a spectrum of the simulator is a step of
   derivation or audit and is labelled as one.
2. **Keep the simulator out of the reader.** State vectors, global expectation values, spectral projectors, the
   evolution parameter and the schedule serve to derive and to audit predictions. No reader, decoder or observer
   operation uses them.
3. **A clock is a reading of one observer's own history.** It orders that observer's records and measures intervals
   in them. A time in a claim is the reading of a named observer's clock, and an energy is a reading with a named
   observer's clock. No test takes a global or common clock as an input. An oscillator, a counter and a pulse
   schedule are possible implementations of an observer's clock and no prerequisites of the theory.
4. **Synchronization is an outcome.** Two observers compare and calibrate their clocks through the records on their
   overlap. A claim that uses a common time coordinate names the overlap records behind it.
5. **A field is an effective description.** A gauge response, a potential, a metric or a wave is tested through the
   probabilities of observers' records. Two descriptions that give every observer the same record probabilities are
   the same physics, and a difference that no record shows is a convention.
6. **Agreement is an outcome of a test.** A claim that observers agree names the overlap, the records compared and
   the outcome that would show disagreement.
7. **A statement at every simulator time is a statement at every reading.** A bound that holds for every value of
   the evolution parameter holds whenever an observer reads, and it is written in that form. A result that depends
   on the simulator's parameter, on a global order of events or on a count of the events of the whole carrier says
   so and is labelled as audit.

When a row uses an effective quantity, give the inference an observer would actually make. For OW-28,
identify the emitter and receiver, the clock or spectral standard, the distance measurement and the records
compared across directions and epochs. Explain how the proposed dynamics yields those readings and under
which geometric and matter assumptions they fix Λ. An alternative OPH mechanism can be compared directly
on the same records. State the domain and accuracy of observational equivalence, a discriminating experiment
where possible, and the public criterion’s remaining clauses.

Every protocol declares the preparations, settings, retained quantum references and outcome records each observer
can access, and uses complete instruments with all complementary outcomes. The decoder maps the available records
to the claimed readings and is specified before the run. Unread quantum state is never presumed recoverable from a
classical tape; a quantum reference or a phase instrument is an explicit part of the protocol.

## 6. The form of every test

Every property is asked with three items, the form the maintainers applied in the M10.2 review:

1. a state the candidate can hold in which the property is absent;
2. a changed law that keeps the alphabet, the routing and the symmetry of the candidate and changes the rule;
3. the outcome that counts against the claim.

A test also runs planted laws that the property has to miss. They test the instrument and earn nothing. A number
fixed by the anchor is compared with nothing; its job is a second number that can miss.

Before a claimed experiment, write into its design and its prospective record:

- the exact row clause and the candidate, with law, state domain, preparation, scheduler or controller, memory, read
  and reset operations, decoder and outcomes;
- which observers hold which records, the observable computed from those records, and the map to the proposed
  physical interpretation, with hidden state confined to a separately labelled audit;
- every supplied ingredient, measured input and calibration, and what is left to predict once they are fixed;
- what the law contributes to the target property beyond identities of the alphabet, geometry, reference and reader;
- a representable falsifying outcome and a control that tests that property with the same observer and readout,
  preserving the symmetry that geometry guarantees when the question is whether the law contributes more;
- separate planted-effect and integrity checks, and the list of identities that earn no emergence credit.

If the declaration guarantees the target property, or the property survives the relevant alternative law, narrow
the claim or redesign the test. The research may continue; the property is withheld from any row award. Calling a
row static leaves the own-process requirement in force, and another column's award cannot replace evidence for this
candidate's mechanism.

| Verdict | Meaning |
| --- | --- |
| pass | the property comes out under the law, and a changed law gives another answer |
| open | nothing fails, and a part of the test is unwritten |
| declared | the property holds in every state the candidate can hold, or under every changed rule that keeps the declared parts: the declaration carries it |
| fail | the property is absent under the law |
| n/a | the candidate provides no state or reading that the test needs |
| no test | the test is unwritten |

A pass or an open result lists the parts of the candidate that are put in by hand: routing, schedule, clock,
memory, origin. The rules of every test:

- Observables come from records.
- A test on carriers runs on the smallest and on the largest carrier of the candidate.
- Seeds and generator are fixed in the experiment, with exact arithmetic where the state space allows.
- A control that lets nothing move is excluded where the observable vanishes with the motion.
- A clause that no finite run decides ("from any start", "at every size") names a theorem on the law or a stated
  range, and its check reads unwritten until it does.
- Every claimed test has a sensitivity and integrity control. A planted effect tests the instrument and earns no
  emergence credit.
- A physical number test names its criterion, tolerance, calibration inputs and reserved comparisons before the
  run. A finite algebraic check is never silently a numerical physical comparison.
- This section and the criteria files of the [scorecard](../model/evaluation/README.md) are ours. Before a row reads
  fail against them, apply caution 7.29.

## 7. Research cautions

These cautions bind research and submissions alike. Each one records a mistake that cost a filing or a day of work.

### 7.1 Do not build in the answer and call it emergence

Do not count a property guaranteed by the state alphabet, geometry, preparation, reference, normalization or decoder
as a dynamical discovery. An integer-valued register does not establish physical quantization. A return period
forced by symmetry does not establish a dynamical rotation mechanism. Rewriting the same bookkeeping as a
topological statement supplies no independent explanation.

Declared architecture and apparatus are permitted model choices. Their presence is no failure, and choosing them
requires no uniqueness proof. The error is crediting the target property to the dynamics when its evidence comes
entirely from those choices.

### 7.2 Do not substitute execution for an explanation of the property

Calling the native event function, accumulating many records and obtaining a nontrivial number do not establish
every physical claim attached to the run. Do not use evidence for one behavior to support a different headline. A
law can affect a magnitude while the claimed symmetry or period is fixed elsewhere. More events, exact arithmetic
and a passing replay do not repair that mismatch.

Do not claim the process selected a state because a preparer supplied it or a spectral certificate recognizes it. A
prepared state, an eigenvector, a stationary distribution and dynamical selection are distinct claims.

### 7.3 Do not let a calibration validate itself

Do not present a fitted value as a prediction or count reproduction of its own anchor as an independent success. Do
not hide choices of scale, reference, state, clock or normalization that were informed by the target value.

Do not quietly recalibrate, change a convention or select another state after an unfavorable comparison and
describe the result as a successful test of the original claim. A disclosed calibration supplies no evidence for
the value it was chosen to reproduce.

### 7.4 Do not assign physical meaning by naming a mathematical quantity

Calling a statistic energy, a label charge or a representation spin does not establish the corresponding physical
interpretation. Do not equate a decay rate, transfer spectrum, record cost or modular quantity with physical energy
without the operational argument the claim requires.

Do not present a separately programmed reference signal as a measurement of the source when the source's actual
state or accepted measurement record never controls that signal. Agreement with a supplied parameter is
insufficient. Do not resolve an ambiguous phase or another aliased reading by supplying the answer that the
measurement is meant to establish.

Do not count reading back a programmed coefficient as an independent derivation of that coefficient, even when the
source controls the gate. The derivation of the source operator, its faithful implementation and an independent use
of its output are different obligations.

Do not mix event counts, reconstructed time and physical clocks without an explicit justified relation. Do not
dismiss unexplained composition effects, such as excess energy assigned to separated noninteracting copies, because
the formula producing them is mathematically correct. Do not use hidden simulator state as if it were a measurement
available in the observer's records.

### 7.5 Do not use controls that miss the claimed mechanism

A control that breaks geometry or changes the reader cannot by itself isolate the contribution of the native law. A
changed law that differs only outside the tested state domain does not challenge the tested property. Do not count
a control's failure on one quantity as evidence for another.

Do not substitute corrupted files, impossible source states or deliberately inserted signals for a scientific
alternative. They may test integrity or sensitivity; they do not establish that the mechanism generates the signal.
Do not interpret a null from a detector that misses a planted effect.

A result shared by another law can be universal. Do not claim that such sharing alone disproves the result, or that
the original experiment isolated the proposed mechanism when the relevant alternative was never distinguished.

Do not credit a selected symmetry to the dynamics when the replacement state carries that symmetry beforehand and
the only control disables replacement.

### 7.6 Do not manufacture agreement by discarding the distinction

Do not infer a physical identification from a reader that removes the very information at issue. Squaring a vector,
taking an absolute value, dropping a sign, averaging away orientation or collapsing every state to the same outcome
can force agreement without the claimed dynamics.

Do not credit an imposed counter period, representation sign or planted phase as a property generated by the
source. A convenient projection cannot establish what the discarded records would have shown.

### 7.7 Do not transfer results between different candidates

Unchanged endpoint code does not make two experiments the same model. State domain, preparation, scheduler,
controller, memory, read and reset operations, decoder and outcome alphabet also belong to the candidate.

Do not combine successful clauses from different candidates into one accepted model without an implemented and
tested common-source connection. Do not silently change the candidate to repair an objection, substitute another
tree for the reviewed commit, or overwrite sealed evidence. Ancestry does not transfer validation to a descendant.

### 7.8 Do not strip a theorem or result of its hypotheses

Do not promote a bounded-density result to arbitrary loads, a finite sample to an all-input theorem, or a
conditional construction to an established physical fact. Do not generalize a negative beyond its source, apparatus
and observable. A finite timeout does not automatically refute an asymptotic claim.

Do not credit classical branch averages with coherent-memory dependence because the package contains a coherent
component. Do not borrow a certificate, proof or physical interpretation from another domain without establishing
that its hypotheses apply to this experiment.

Do not substitute a small population outside a subspace for a small distance between full coherent states. Cross
terms and the actual input domain matter. Do not silently replace a recorded source state with a convenient
prepared state when joining two stages, or mistake an encoding identity for an executed preparation and complete
acceptance instrument.

Do not infer agreement of a joint state or joint observable from agreement of its marginals. Correlations, a common
physical transformation and the full outcome law must be checked. Independently chosen coordinate relabelings do not
establish one spatial rotation of the combined object.

Do not call a preparation assumption removed when it enters through the census, reference, controller or decoder.
Distinguish the units stabilized by the source law from the number of units a controller chooses to retain. Hidden
input labels can also enter through random seeds, lookup tables, file names or run order. Do not call a reader blind
while one of those channels encodes the quantity it is meant to infer.

### 7.9 Do not conceal unfavorable or inconclusive outcomes

Do not drop failed trials, timeouts, invalid or complement outcomes, lower modes, background ambiguities or
unsuccessful alternatives. Do not count inconclusive measurements as successes or change the denominator after
seeing the results. Do not remove a known counterexample from later summaries.

A valid negative result is evidence. An unreliable detector, an incomplete outcome record or an unverifiable run is
a defect in the evidence. Do not confuse these cases to claim success or to suppress a legitimate negative result.

Do not omit failures incurred before the final detector when composing instruments. Do not report conditional
readings without their acceptance probability and propagated error. Machine-precision agreement is not the physical
uncertainty when the complete implementation error is larger.

Do not label signed numerical comparison values as physical probabilities. Do not hide trace drift by clipping or
renormalizing it, assume an approximate map has exact unit gain, or turn a small numerical residual into a stronger
claim than the certified complete instrument supports.

### 7.10 Do not register discoveries retrospectively

Gates written before replaying a known result do not make that result a prospective prediction. Do not hide
exploratory runs, target-informed choices or prior knowledge. File timestamps alone do not establish the claimed
order of design, registration and execution.

Do not revise a failed gate, erase the earlier version and report the amended run as if the revised terms had
always applied. Reproduction and prospective testing support different claims.

### 7.11 Do not substitute a large archive for reviewable evidence

A small wrapper that fetches and executes a large package is no small execution payload. Do not assume a reviewer
has inspected or executed every file because the archive is available or its hashes match. The maintainers declined
to execute such a package in the M10.2 review ([safety gate][safety]).

Do not count events represented by enumeration, matrix powers or compressed sampling as native calls unless the
native function executed them. Do not inflate actual call counts with represented possibilities. An independently
reproduced headline does not verify every script, proof, execution count or gate in its source package.

Do not let a failed scientific gate produce a successful process exit. A verification command compares with the
saved evidence; silently replacing that evidence with its own output is regeneration and verifies nothing. Keep
first-record creation explicit and test the actual failure verdict and exit.

### 7.12 Do not inflate independence, verification or acceptance

Codex audits each submission. A separate Codex reviewer within the project can independently challenge a result.
Name the auditor, identify their checks and record their verdict on each claim. A full-packet review covers every
claim in the submitted version.

Do not equate arithmetic agreement with physical validation, or a merged pull request, an adopted pin or an
approved record with a criterion award. Do not claim a proof, payload or gate was checked when the reviewer
explicitly excluded it. Do not report local expectations or stale summaries as the official decision.

### 7.13 Do not reinterpret the criterion to excuse missing evidence

Do not use the word static, another column's award or a narrow reading of a row label to bypass the column's
own-process requirement. Do not treat a disclosed assumption as automatic grounds for partial credit when it
supplies the very property being claimed.

Do not turn an unresolved objection into a harmless caveat by renaming it a physical-dictionary question. Equally,
do not invent obligations beyond the claimed criterion or require every physical property for a narrower claim. The
problem is a claim exceeding its evidence. Limited results that are useful are welcome at their scope.

### 7.14 Do not let polished packaging override a scientific blocker

Passing formatting, custody, replay and arithmetic checks does not establish submission readiness. Do not approve a
packet without challenging the causal and physical meaning of its headline. Do not treat a detailed response, a new
implementation or a larger evidence bundle as proof that a criticism was resolved.

Do not submit an unsupported full or partial claim hoping the maintainers will supply the missing argument. Do not
bypass a failed checklist item with softer wording, a review-request label or an unexplained N/A. Documentation of a
problem does not satisfy the scientific obligation it describes.

### 7.15 Do not compute before reading

This repository holds the finite theorems, their premise register, the executable packages and the earlier OpenWave
evidence of each candidate. Before a study, search them for the same question and record the paths with what they
say: confirms, extends or contradicts. A repeated result earns nothing, and a lead that a filed theorem excludes
costs its time twice.

### 7.16 Do not read the dimension of space from a count of cells

The number of cells within a given number of steps, a growth exponent and a spectral dimension are properties of a
routing put in by hand. The theory's statement about three directions concerns one cell: the slowest nonconstant
band of its response is the band of the three port coordinates, and a position is a reading of the cumulative
record of port currents. Quote that statement with its limits: the schedule it is proved for, the premises its
receipt leaves open, and whether the band survives on the federation the claim uses.

### 7.17 Do not treat the schedule as an implementation detail

Which seams act, in which order and on whose clock is part of the law. The schedule decides the symmetry of a tick,
the speed limit, the energy reading that adds for separated objects and the counts that a register uses. Quote a
number with the schedule it was computed under. A count of the events of a whole carrier is a counter of the solver
and no clock of an observer, and a result that depends on it says so.

Do not take the attractor of one schedule for a property of the law. A bound state, an orbit or a period found under
one order of events and one value of a free number is tested under the other orders the law allows and across that
number before it is named.

### 7.18 Do not credit a force to a law whose records never act back

A force between objects needs records that enter later events. A law under which the records of a unit are the same
in every background has no force to credit. A claim of a force names the event that reads the records and the
control that scrambles them.

Do not read a record law from a mean over runs. The mean record of a source can look like a lattice field while every
single run holds integer currents. A claim about what observers read is tested on the records of single runs, and a
mean is labelled as a mean.

### 7.19 Do not pass on a summary as the document

A statement from a reading report, a worker or an earlier summary enters a result after its quotation was checked at
the cited path and commit. Record the commit. A report that says what it did not open is complete; one that fills
the gap from memory is a defect.

### 7.20 Do not write that an observer reads a field

An observer holds records. The field, the particle, the distance and the clock are descriptions fitted to the
records that observers agree on. A sentence in which an observer reads a field off a state reverses the premise and
hides the step that has to be tested, the step from records to the description.

### 7.21 Do not join results that stand on different constructions

A row result names its construction, for example the literal law with committed records, the same events kept
unread with a declared reflection, or a declared field source on which no repair event occurs. The word native
covers only the first. Results of different constructions are not added up, and one physical thing has one object
and one reading across all rows.

### 7.22 Do not hide a match with nature and do not inflate it

Record every agreement with a measured or standard number as a sanity check, with the date it was seen, the choices
behind it and the alternatives that were tried. Whether it was seen before or after the number was known decides
what it can support, and it is recorded in both cases.

### 7.23 Do not mark a setting of the test against the law

A start state, a carrier, a joining of cells and a schedule are parts of a test. Before a line reads "does not hold",
check that the state is one the candidate can hold, that the carrier is the one the property is defined on, and that
the row's simplest test asks for the thing. Typical errors of this kind are a declared background above the capacity
of the carrier, cells joined by a few bridges taken for a carrier in three dimensions, and an ingredient that no row
asks for.

### 7.24 Do not compute a premise on a state that the candidate does not hold

A premise about content in a region speaks of states that the process holds. A unit fixed on one port is no such
state, and fixing it removes the entropy of its places. An entropy computed on fixed content can carry the opposite
sign to the same premise computed with the content free to move. Write the ensemble beside every regional number:
which law over places, conditioned on what.

### 7.25 Do not word a result beyond its script

A sentence in a result page is tested by a script that fails under a changed rule. An audit with its own
implementation can confirm every number and find sentences that claim more than the scripts test: a
symmetry named for another, a stationary law stated without its set of states, a factor quoted as a range it did not
span. Run a mutation test: change the rule in each way that matters and record which checks notice. A check that
fails under no change tests something else.

### 7.26 Do not take the first object for the object of the row

A row names a physical thing. The record that stands for it is found in the theory's documents before a test is
written, and every record of the candidate that could be read as that thing is listed and tested. The theory's
particle programme reads a particle from the class of a protected defect
([defect classes](../code/defect_particle_classes/README.md)); a count of units of load does not replace it.

### 7.27 Do not judge a law on one member of its memory family

A repair law has two parts: where the load goes, and what the event writes down for an observer. A verdict on records
that are all read at every event is a verdict on one member. Name the member with every row verdict (all records
committed, the drawn direction unread, the hidden state of a blank seam unread, partial events), and classify the
family before a row that needs a phase, an energy or a sign is closed.

### 7.28 Do not close a verdict before the causes are worked in order

A clause that does not pass has one of five causes: the calculation, the observable or the way it was derived, an
ambiguity of the requirement, a requirement at variance with a measured fact, or the law. They are worked in this
order, with this repository and published work searched for the object that the theory means, and the record says
what was checked for each. The fifth cause is stated in one plain sentence when the first four are excluded. It is
never softened, and a sibling of the law that closes the gap is tested, registered and named as a sibling in every
sentence.

A failed candidate assignment belongs in its source record. Before carrying that failure into a row forecast,
identify the physical observable on the same preparation through its observer records. Check the energy reference,
density convention, finite-size dependence or other calibration the comparison uses. An alternative that gives the
same allowed records leaves the assignment open.

### 7.29 Do not hold our own contract for the platform's requirement

The form of test in section 6 and the criteria files of the scorecard are ours, and our own lines can ask more than
the platform does: an energy conserved in committed records, a count of cells for the dimension of space, an object
that keeps its place. The simplest tests of the platform are one model's vocabulary, and the maintainers have set
mechanism words aside on record. Before a row reads fail, quote the platform's wording, the precedent of another
column and the maintainers' text to this column, and say which words are a fact of nature.

### 7.30 Do not quote the platform's bar from a stale reading

Record the upstream head before new work (section 2). Maintainer answers land in the open threads of the column, and
an evaluation run on notes taken at an older head can miss an answer that stood upstream for a day. Read the open
threads at the head before a verdict standard is written down, and record the head.

### 7.31 Do not take the state of the simulator for the universe

What is compared with nature is the readout of observers from their own records. The state vector, the evolution
parameter and the schedule of a simulation compute the probabilities of those records. A result worded on the
simulator's state, on its time or on a field in it has skipped the step to the records. Word the result on records,
with the observer named, and label what was computed on the state as audit.

### 7.32 Do not require a common clock

A clock is a reading of one observer's own history of records. A construction that needs a clock shared by all
observers before any record exists reverses the order: records first, the comparison of clocks through overlap
records afterwards, a common time as an agreed description last. An oscillator, a counter and a pulse schedule are
implementations and no foundations. A claim about synchronized clocks is a claim about overlap records and is tested
as one.

## 8. Rows 1, 2 and 7

The [M10.2 maintainer review][review] decided rows 1, 2 and 7 on the pinned signed-source candidate and named what
would move each. The reasons are general, and a filing on any of these rows answers them at its scope.

| Row | Why it stayed 🚧 | What would move it |
| --- | --- | --- |
| 1, charge quantization | The loads are integers by the carrier's definition and the current is an integer by its formula, so no run could return half a unit. The topological statement and the integer bookkeeping were one statement. | The unit as an outcome of the law: larger loads allowed, with single units the only settled records against a changed-law control; or fractional loads native, with the law confining settled records to whole units; or an integer class the law cannot change, computed where other values are representable. |
| 2, electron rest energy | One anchor fixed one number and the row offered one number. The anchored mode was a chosen one with lower modes below it, two separated copies came out 2.45 keV above twice the anchor, and no physical clock was established. | A state the dynamics holds or returns to, and the anchored quantity doing a second job that can fail: exact additivity for separated copies, or a second unanchored number that agrees with measurement. |
| 7, spin-½ statistics | For any rule equivariant under the cell's rotations, the cell's symmetry fixes the half-turn result; the law contributed magnitudes and a sign. | A factor of two that the process produces and that the cell's symmetry does not fix: a quantity that changes sign under a full turn and returns after two, read by interference; or a field whose return under a half turn is an outcome, against a control law that keeps a direction and its opposite apart. |

The checks that follow from the review:

- **Charge.** Integer states or currents and equivalent homology bookkeeping do not derive quantization. Test
  law-selected units, or another route the review supports, in a domain that allows meaningful alternatives. State
  density limits, background ambiguity and counterexamples. A unit result on a bounded range of loads is never
  generalized to every integer load.
- **Rest energy.** Require dynamical selection, holding or return of the localized state, and a second falsifiable
  use of the energy: separated-copy additivity or an unanchored prediction compared with measurement. Retain lower
  modes and any separated-pair excess. Declare the clock and the composition convention. A transfer eigenray, a
  decay rate or an additive modular statistic is no physical rest energy without the supporting operational
  dictionary. A changed decoder under which additivity holds identically is a new candidate, and it cannot supply
  the job that could fail.
- **Rotation.** Test a symmetry-preserving competing law; a covariance that returns under a half turn is
  insufficient alone. For dynamical antipodal identification, retain directional information and an oriented reader;
  an even projection or vacuum erasure earns no credit. For the coherent route, measure the full-turn sign and the
  two-turn return through branch interference with finer turns. A planted phase only tests sensitivity. Do not
  impose a full exchange theorem on every narrower claim.

Every follow-up to the M10.2 review answers three questions, each linked to evidence or to a named open obligation:
what could falsify each result; which results need coherent memory; and how the separated-pair energy excess is
interpreted or left unresolved.

## 9. Filing gates

Codex audits each outgoing submission at its declared scope. A row submission needs a full scientific audit of every
clause of the official criterion and a Codex finding that full credit is expected.

| Gate | Evidence | Record | Confirmation |
| --- | --- | --- | --- |
| G0 design | Every clause of the row's simplest test is listed. For each: the observers, the records they hold, the observable as a function of the records, the description to be fitted, the control. The candidate is registered, everything supplied is listed, and the sentence of the theory that says the records should show the clause is quoted with file and line. | `DESIGN.md` | independent design reviewer |
| G1 terms before the run | Gates with numbers, inputs, seeds (development apart from held-out) and controls are written and hashed, and the hashes are pushed to this repository before the run. A result that will claim a cell files its upstream task document with its gates and inputs before the run. | `PREREG.json`, the upstream task document | reviewer checks the order from the pushed history |
| G2 execution | The run uses the filed inputs, and the law's implementation is byte-identical to the registered source. The record states whether execution uses literal events, full matrix or column actions, enumeration or justified compression, and proves any reduction. Every outcome is kept, failures and timeouts included. | `RECEIPT.json`, logs | reviewer replays a subset |
| G3 controls | (a) A changed-law control fails on the row's own object for the stated reason; (b) the instrument detects a planted effect; (c) no gate is an identity of the alphabet, the graph or the decoder, and none can only pass; (d) every statistical statement names its interval, its selection and its power. | `RECEIPT.json`, `CONTROLS.md` | reviewer runs own changed laws |
| G4 reproduction | A script regenerates the headline number in under 60 seconds, a replay is byte-identical, and the manifest covers every file. | `headline.py`, `MANIFEST.sha256` | reviewer runs both |
| G5 method note and cell | A method note in the platform's [form][methodnote]: equations first, a map from every term to code, a scope line, what was not computed. Cell text of at most 65 words of prose as the platform's [cell format][results] counts them, with the expected mark and every supplied piece. | `METHOD_NOTE.md`, `CELL.md` | reviewer checks word by word |
| G6 independent scientific review | Codex independently challenges every claim in the full packet, binds the verdicts to the manifest and final wording, and verifies the resolution of each blocking finding. Component checks are recorded at their tested scope. | `REVIEW.md`, `REVIEW.json` | Codex auditor |
| G7 one candidate | The packet runs on the pinned candidate, or the pin-move statement (the change in law, memory, state or decoder) is written and the owner has decided. No clause is borrowed from another candidate. | `PIN.md` | integrator and owner |
| G8 full criterion | Every clause of the simplest test is carried in the full sense on that candidate. Nothing that the criterion asks to follow from the model is supplied or fitted (rows 20, 26 and 28: the strength). The evidence for each clause is recorded, and an open clause keeps the row on HOLD. | `CLAUSES.md` | independent reviewer |
| G9 filed | The pin move, the registered task with its gates, the in-platform run, the method note and the script are upstream, and the maintainers have decided. | the maintainers' decision, recorded with the packet | the owner only |

Prospective designs, negative results and component measurements belong in the research record. Codex marks a row
**SUBMITTABLE** only when the evidence on the nominated pin supports an expected full mark (✅) and the exact outgoing
package passes the checklist. The owner then decides whether to release it.

## 10. Mandatory submission checklist

Keep this master unchecked. Copy it into the row's packet as `CHECKLIST.md`, complete it against the exact outgoing
material, and retain it with the audit. Use the scientific checks at the scope of the requested decision.

### Identify the submission

| Field | Required entry |
| --- | --- |
| Requested decisions | For a row, full credit (✅) on its complete official criterion. Any accompanying pin change, and each requested decision of a combined pull request. |
| Rows and claims | The literal criterion, the dated upstream mark, the precise claim, the requested mark, the supplied inputs and the open clauses of each row. |
| Law and evidence | The pin, the proposed immutable pin if any, the reason for the correction, the packet paths and the manifest hashes. |
| Public source of a row filing | The commit of this repository that carries the sources, the row's directory under `openwave/`, the manifest of the exact publication files, and the commands that reproduce the row from a clean copy of that commit. |
| Simulator | The public URL, the mode or view, the implementation and the observable it demonstrates for each row. |
| Outgoing version | Destination, base and head or exact diff hash, path and hash of the pull-request text, and attachment digests. |
| Participants | Producer, Codex auditor, integrator and owner. |
| Review record | Date, time zone, audit path and the completed checklist. |

A single pull request may present a corrected candidate and several related rows. Use one entry per item and per
affected claim, so that reviewers can decide each row separately. If a demonstration uses an implementation other
than the row run, identify both sources and what each establishes.

### Apply the checklist

- **All submissions.** Codex audits the final diff, the pull-request text, the links and the attachments in a
  dedicated pass. A separate Codex reviewer can audit work that Codex produced within the project. Codex owns the row
  and package verdicts. The owner records the authorization to send under F4.
- **Scientific decisions.** Apply A to E to every scientific claim at its stated scope, with the evidence for each
  clause and the limitations of the instrument or model. Approval requires the law to meet the complete criterion,
  an explicit Codex expectation of full credit (✅), and a full-packet G6 audit. A scoped or partial result stays on
  HOLD.
- **Prospective work.** Check the proposed inputs, gates and controls before the run that is to claim a cell. A
  registration-only row proposal stays on HOLD until the full criterion is demonstrated and the packet requests full
  credit. Record any earlier draw as historical and disclose what was known. An evidence-only or negative report
  needs reliable measurement, complete outcomes and a reproducible result for the claim presented.
- **Administrative work.** Audit the exact edit and every factual or scope claim it contains. Name the scientific
  checks outside that scope, with a reason.
- **Auditor verdict.** Record `PASS`, `FAIL`, `OPEN` or `N/A` for each item, with an evidence path or commit, the
  observed result, the auditor and the date. Codex confirms each scope exemption. It marks the exact row and version
  **SUBMITTABLE** only when every criterion clause warrants full credit and every applicable item from A1 to F3
  passes. A failed or open applicable item, or missing evidence, means **HOLD**, recorded with its reason and what
  would clear it. A revised or unaudited version needs its own verdict before release.
- **Combined submissions.** Check the shared pin, the cover and the other common material of the actual pull request.
  The package is submittable when every included row and every shared part clears its applicable checks. Recheck the
  affected items when a held row is removed and the package is revised. Record the owner's authorization under F4
  separately from the auditor verdict.

### A. The claim meets the actual criterion

- [ ] **A1. Terms read together.** The exact row clauses, the M10 admission terms and every applicable maintainer
  decision are cited with their versions. "Static" and another column's precedent do not waive the own-process
  rule.
- [ ] **A2. Evidence per clause.** Each claimed clause has an evidence-linked assessment: full, weaker, absent or
  negative. The requested mark follows that assessment, and for approval every criterion clause warrants full
  credit. The record says which properties follow from the process and lists every supplied ingredient. A weaker or
  absent clause leaves the row on HOLD.
- [ ] **A3. Law contribution established.** For each claimed property, the record states what the dynamics
  contributes and shows that the headline is not guaranteed by the alphabet, preparation, geometry, calibration or
  decoder. Native execution alone does not satisfy this item.
- [ ] **A4. Operational interpretation supported.** The chain from the process to observer records, observable and
  physical interpretation is explicit and checked. Hidden simulator state is confined to labelled audits. Calling a
  quantity charge, energy or spin is supported by the claimed tests. A reference signal is coupled to the actual
  source or its retained measured outcome; source-independent calibration alone is insufficient. Where a compiler
  consumes a derived source operator or coefficient, both its derivation and its faithful implementation are
  checked, and reading back the programmed value counts as no independent physical test.
- [ ] **A5. Existing work and limits retained.** Relevant established results, theorem hypotheses, prior failures
  and dependencies are read and cited. Neither positive nor negative conclusions exceed the tested domain. The search
  of this repository and of earlier studies for the same question is recorded with paths and commit, and every
  quotation was checked at its source.

### B. The headline can fail for the right reason

- [ ] **B1. Falsifying outcome specified for every claim.** The record answers **"What possible run outcome would
  falsify each claim?"** Before the claimed run, it records for each scientific claim the counterexample or outcome
  representable in the test domain, the measured quantity and detector, and the exact rejection condition or
  statistical threshold, with evidence that the detector can register the relevant effect. Proof counterexamples,
  statistical rejection and inconclusive timeouts are distinguished. A corrupted-input failure or a gate guaranteed
  to pass by construction does not fulfil this item for a claim of dynamical emergence.
- [ ] **B2. Relevant alternative law tested.** The same observer and readout test an appropriate changed law on the
  claimed object. The result separates the target property from the competing law, with the relevant geometry and
  symmetry preserved. A changed magnitude alone does not establish a return-period claim.
- [ ] **B3. Detector and identities audited.** Planted effects are detected; integrity mutations and identities are
  listed separately. An impossible half charge, a planted phase or a corrupted file receives no emergence credit.
  Aliased readings are resolved without inserting the claimed answer.
- [ ] **B4. Coherence dependence stated.** Each claim says whether it uses classical branch statistics, coherent
  amplitudes or a conditional reconstruction, with the supporting calculation or control. Coherence is never
  credited because another component of the package uses it.
- [ ] **B5. Records of a run.** A claim about what observers read is shown on the records of single runs, and a mean
  over runs is labelled as a mean. Where the two differ, as for integer currents against a mean field, the claim
  states which one it uses. A claimed force names the event that reads the records and the control that scrambles
  them.

### C. Calibration and the M10.2 objections are resolved at the claimed scope

- [ ] **C1. Inputs do not validate themselves.** All measured constants, fitted choices, scales and selection
  decisions are listed. Every calibrated claim identifies the independent test left after fitting and reports its
  outcome, misses included. Clocks, units and composition conventions agree. Conditional acceptance and the complete
  implementation error are propagated into the reported interval of the observable; floating-point residuals do not
  replace that interval.
- [ ] **C2. Row 1, if claimed.** The evidence follows a route the review supports to law-selected units or a genuine
  configuration invariant. Supplied integrality and equivalent load-vector topology earn no quantization credit.
  Density and background limits are explicit, and known counterexamples are kept.
- [ ] **C3. Row 2, if claimed.** State selection is attributable to the dynamics beyond the preparer's choice, with
  holding or return demonstrated. The anchored energy passes either separated-copy additivity or an unanchored
  prediction compared with measurement. Lower modes and any separated-pair excess have evidence-linked dispositions.
  A transfer eigenray or an additive record cost alone is insufficient.
- [ ] **C4. Row 7, if claimed.** The evidence defeats the symmetry-preserving alternative. The supported route is
  either dynamical antipodal identification that retains directional information, or a full-turn sign and two-turn
  return measured by interference with finer turns. Even readers, vacuum erasure and planted phases are excluded as
  explanations of the claimed success.
- [ ] **C5. Review questions answered.** Every applicable prior criticism has a disposition and evidence. A
  follow-up to the M10.2 review answers its three questions: possible falsifiers, actual coherence dependence, and
  the separated-pair excess. Any remaining issue lies outside the requested claim and is explicitly scoped; a
  blocker is never relabelled as a caveat to obtain partial credit.
- [ ] **C6. Space, if claimed.** Three directions and a distance rest on the candidate's position record and its
  far-field reading, with the limits of the rank-three theorem stated. A count of cells, a growth exponent, a
  spectral dimension or a declared three-dimensional routing earns no credit for dimension.

### D. One candidate and a traceable experiment

- [ ] **D1. Complete candidate bound.** Law, state, preparation, scheduler, controller, memory, read and reset
  operations, decoder and outcome alphabet are specified, and changed protocols are registered as new candidates
  with ancestry. The schedule is stated as a partial order of events or as local clocks, every schedule-dependent
  number of the claim was computed under it, and any dependence on a count of the events of the whole carrier is
  stated.
- [ ] **D2. No borrowed successes.** All claimed clauses run on the same candidate, with any common-source join
  implemented and tested. The source commit is retrievable and preserved, and a different tree of the main branch
  cannot replace it. Any pin proposal has the owner's recorded decision. A whole-object claim checks the joint state
  and observable, including correlations and one common physical transformation; matching marginals or separate
  relabelings are insufficient. Removed preparation assumptions do not reappear in the observer's references, inputs
  or controller. Random seeds, lookup tables, file names and run order leak no hidden input label into a reader
  described as blind or record-only.
- [ ] **D3. Timing is genuine.** Design, code, inputs, seeds, controls and numeric gates were frozen before the
  claimed confirmatory run, with pushed-history evidence and the applicable upstream task registration. Earlier
  pilots and discoveries are disclosed. Historical reproduction is labelled as such and never backdated into a
  prospective test; amendments preserve the original record.
- [ ] **D4. Execution matches the claim.** Receipts bind the frozen source and inputs. Literal native calls,
  enumeration, matrix powers and sampled histories are distinguished, and any compression has a checked equivalence
  argument.
- [ ] **D5. Full outcomes retained.** Failures, timeouts, invalid and complement outcomes and failed alternatives stay
  available. Statistical claims state uncertainty, selection and power; finite samples establish no universal claim
  without the required proof. Composed instruments retain upstream as well as terminal failures. Numerical comparison
  centres and trace drift are distinguished from positive physical outcome laws, with gain and conditioning bounds
  where required.

### E. A reviewer can inspect and reproduce the submitted evidence

- [ ] **E1. Small offline check.** A reviewer can read the complete primary execution path and its bounded inputs.
  The headline reproduces in under 60 seconds, with command, environment, runtime and exit status recorded. A
  download-and-execute wrapper does not satisfy this item. Failed scientific gates produce a nonzero exit, verified
  by mutations of the actual checks. Verification preserves saved evidence, and an explicit first-record mode is
  distinguished from read-only replay.
- [ ] **E2. Clean replay and custody.** A clean copy of the submitted commit reproduces the reported results, and
  manifest verification passes. Portability accommodations and provenance limitations are explicit and reviewed.
- [ ] **E3. Scientific adversarial review.** For each scientific claim, Codex attempts to refute the interpretation,
  reproduces or independently checks the headline, and tests a relevant competing law or counterexample. The record
  holds its verdict, the evidence and the producer components the audit used.
- [ ] **E4. Mandatory Codex audit.** Codex records a verdict on the exact outgoing version of every submission. The
  audit binds the applicable candidates and manifests, the claim wording, the attachments, the diff and the
  pull-request text. An award request also records the complete G6 audit and the applicable G0 to G8 confirmations.
  A second Codex worker within the project can perform the independent audit. Each blocking finding is resolved and
  its retest recorded.
- [ ] **E5. Repository checks pass.** The packet files, manifests, the candidate's `EVALUATION.json`, its entry in
  [openwave_criteria.json](../model/evaluation/openwave_criteria.json) and the review records agree. The packet's
  replay and tests pass, and these checks pass from the repository root, with their output attached:

  ```sh
  python -m model evaluate --check
  python tools/check_private_links.py
  python tools/check_docs_policy.py
  python -B -m pytest -q model/candidate_laws/<law>/openwave openwave/<law>
  ```

  Submission-related violations are resolved; unrelated workspace diagnostics are identified separately. A clean
  checker replaces none of the scientific checks above. G9 waits for the actual upstream decision.

### F. Final submission and release

- [ ] **F1. Text matches the evidence.** Task, method note, receipt, proposed cell and pull-request description agree
  on candidate, result, supplied ingredients, limitations and requested mark. The cell meets the 65-word limit.
  Official decisions are attributed to the maintainers. The simulator URL was opened and the view, law and
  observable named for each row were verified. Any claimed continuum limit carries its derivation and hypotheses.
  Remaining criteria and anticipated acceptance are identified as research goals.
- [ ] **F2. Outgoing diff and source reviewed.** The final OpenWave diff, pull-request text, paths, existing links
  and contribution requirements are checked against the intended base. The source files at the commit of this
  repository named in the pull request match the bound manifest, trace to the candidate, and pass the documented
  checks from a clean copy using only reviewer-accessible inputs and dependencies. The planned source link names
  that commit. Concurrent work is preserved.
- [ ] **F3. Completed audit record signed.** The Codex auditor records the disposition of A1 through F2, with
  evidence for every applicable check and a reason for each `N/A`. All applicable checks through F2 pass. Producer and
  auditor identify the exact reviewed packet and outgoing diff, and the completed checklist is included or linked in
  the draft.
- [ ] **F4. Owner authorization covers this version.** The owner's final authorization for the audited version and
  the actions it covers are recorded: pushing the source commit if it is unpublished, and opening the OpenWave pull
  request in the same release. The public commit is verified anonymously against the audited manifest before the
  pull request opens, and both immutable links are recorded after filing.

### Completion record

Fill one row per item and per affected claim; the first row is a template.

| Item / row / claim | Status | Evidence and observed result | Reviewer and date |
| --- | --- | --- | --- |
| A1 / … | OPEN | Required source and comparison pending | Pending |

| Decision | Required record |
| --- | --- |
| Mandatory Codex audit | Auditor, scope, pin and evidence hashes if applicable, checks performed, verdicts on each claim or requested change, open findings |
| Final diff review | Codex auditor, exact diff or commit, hashes of the pull-request text and attachments, review date |
| Codex row verdict | **SUBMITTABLE** requires the scope "full criterion" and the expected credit "full" for each row, with clause-by-clause evidence in the audit. Otherwise **HOLD** with its blocker. |
| Codex package verdict | **SUBMITTABLE** or **HOLD** for the exact combined diff, cover and included rows. A single-row filing may use the row verdict for its packet. |
| Owner release (F4) | Pending, or the owner's authorization reference, the exact version and the date |
| Sending status | **READY TO SEND** when Codex has marked the package **SUBMITTABLE** and F4 records the owner's authorization for that version |
| Filed row | Public source commit URL, verification of the published files, pull-request URL and filing date |

The maintainers' decision and any official cell award are recorded apart from this pre-submission record, and only
after the maintainers decide.

## 11. The packet

A candidate's research evidence for a row lives in its folder, `model/candidate_laws/<law>/openwave/`, one folder
per row study; the scorecard's OpenWave entries cite it. The public source of a filing is staged under
`openwave/<law>/<row>/` at the repository root: a short README with the row summary, the reproduction commands and
the simulator page, the candidate source, input data and dependencies a reviewer needs for the claimed replay, and
their provenance back to the pinned candidate. The README of `openwave/` lists filed items with their OpenWave pull
requests and names prepared rows separately. A scientific packet contains:

| Artifact | Required contents |
| --- | --- |
| Design | Exact target clauses, candidate and domain, observers, preparation, records, interpretation, supplied inputs, falsifying outcomes and controls (G0) |
| Prospective manifest | Hashes of source, inputs and gates, prior knowledge, and the pushed commit before the confirmatory run (G1) |
| Code and receipt | The actual execution method, all outcomes and failures, environment, numerical bounds and results (G2, G3) |
| Proof or method note | Equations, hypotheses, source attachment, scope and unresolved obligations (G5) |
| Reproduction entry point | A bounded offline headline check and any explicitly larger replays (G4) |
| Codex audit | Auditor, version and manifest binding, scope, independent checks, claim-level verdicts, responses and outstanding findings (G6) |
| Completed checklist | `CHECKLIST.md`, the copy of section 10 completed for the exact outgoing version, with the row and package verdicts |

A prospective design is a valid deliverable before its execution. A negative or inconclusive result is a valid
result when its method and scope are sound. A reproduced known result is never relabelled as a blind prediction.
Evidence is relocated with a path and digest manifest, and a frozen receipt is never rewritten to describe a
different candidate. When the evidence changes an assessment, update the candidate's `PREMISES` or OpenWave entry and
regenerate the evaluation with `python -m model evaluate --write`.

## 12. The upstream pull request

The upstream column runs headless first, one concern per pull request ([M10 roadmap][roadmap],
[roadmap standards][roadmapstd]):

1. A task gets its roadmap row and its `tasks/m10_<n>_task_details.md` document, with the next unused `M10.<n>`
   identifier, before its scripts take the identifier. A result that will claim a cell files its gates and inputs in
   that document before the run.
2. Every script, data file and finding carries the identifier of a registered task. Earlier pilots, choices informed
   by the target and laws set aside are recorded. When the result is known beforehand, the run is labelled a
   reproduction, with a statement of what an independent recomputation can establish.
3. The run uses the author's compute. The full record stays in this repository at an immutable commit; the pull
   request brings the method note and a small script that regenerates the headline number.
4. The pull request comes from a focused branch of a fork, with a DCO sign-off on every commit (`git commit -s`) and
   maintainer edits enabled ([contribution guide][contributing]). It stays a draft while commits land. Review
   findings are resolved in the same pull request, and affected evidence is rechecked when the outgoing version
   changes.

| Item | What a reviewer needs |
| --- | --- |
| Roadmap and task record | A short roadmap entry linked to the task document, which holds the exact criterion, gates and inputs, the run record, deviations, the review and the disposition |
| Method note | Under `research/findings/`, in the platform's [method-note form][methodnote]: the local update law, the observer instrument and the observable in mathematics first; terms mapped to functions with immutable line links; the result beside its gate with its convergence or uncertainty; a scope line and the adversarial audit |
| Headline check | A small headless script under `research/scripts/` with inspectable inputs, explicit dependencies, one offline replay command, the recorded environment, runtime and exit status, a nonzero exit on a failing gate, and enough raw inputs for a reviewer to recompute the metric independently ([reproduction guide][reproduce]) |
| Pinned record and receipt | Links to the public source and evidence at the immutable commit of this repository, naming the candidate. Ship only what the upstream recomputation needs; a reproduction or verification run adds its JSON receipt under `research/data/` |
| Proposed cell | A local draft with the requested mark, at most 65 words of prose, the result and the supplied or open pieces, followed by evidence links. External contributors bring the result and leave the matrix change to the maintainers unless the agreed scope includes it; a matrix edit aligns the cell, the matrix icon and the score board |
| Simulation or visualization | The owner's public link in the pull request and the method note, with the mode or the steps to the relevant view, its implementation and whether its content is measured or illustrative. The link is opened and its description verified before release |

The linked demonstration helps a reviewer see the run or its records. The script, the method note and the evidence
carry the claim.

What the maintainers do with a filing shapes it:

- **They recompute by their own route.** In the M10.2 review the maintainers recomputed every headline number with
  their own script from one data file of the author package and the formulas of the method note, and imported
  nothing from the package. Write the method note so that the headline follows from its equations and a small
  pinned data file.
- **They read a payload in full before running it.** The M10.2 package of 329 files was too large to read and was
  not executed, so its counts of calls, columns and branches stayed unconfirmed ([safety gate][safety]). A claim
  that only a large package supports stays unchecked.
- **They keep the change inside the column.** An M10 pull request touches only files under
  `openwave/xperiments/m10_oph/`. A change to `MODELS.md`, shared code or a platform document needs the
  maintainers' agreement on scope.
- **They check every link and the shape.** Every source link is pinned to a commit and resolves, and the pinned
  commit is reachable from a default branch or a tag. The files are UTF-8 with one trailing newline and no hidden
  characters, and `check_models_md.py` and `check_roadmaps.py` pass.
- **They read the words as platform claims.** In the review of the column's onboarding the maintainers flagged local
  review vocabulary that a reader can take for a row status, titles stronger than their own boundary lines, and
  results of an open pull request stated in the present tense. Upstream text uses the platform's words for marks
  and criteria, keeps every title within its boundary, and states a result under review as under review.

```markdown
## M10.<n> · Row <number> · <criterion>

<One sentence naming the observer reading, the result and the tested domain.>

- M10 mark at the recorded head, with its source: <link>
- Requested decision: <full credit (✅) and the evidence for each criterion clause>
- Candidate, source and record: <name, commit, manifest and receipt links>
- Falsifier and changed-law outcome: <prediction, gate and observed outcomes>
- Calibration and limitations: <supplied pieces, held-out misses, open clauses>
- Method and replay: <method note, command, expected output and runtime>
- Review and scoped findings: <audit and completed submission checklist>
- Simulation or visualization: <owner-supplied URL and what it shows>
```

## 13. Audit, verdict and release

1. **Bind the version.** Record the upstream repository and base, the exact patch or head, the pull-request wording,
   the applicable pins, the evidence manifest digests and the simulator links. The row's public sources are the
   files of its packet at one commit of this repository. Use a separate checklist entry for each row or materially
   different claim.
2. **Audit the actual request.** Codex performs a dedicated adversarial pass on the complete outgoing material. For
   a scientific claim it reproduces the result, challenges the physical interpretation and exercises a control that
   can fail. For a registration it audits the prospective terms and the disclosure of prior knowledge. For an
   administrative filing it audits the exact change and its declared scope.
3. **Work the findings.** Codex records each actionable finding with an identifier in the packet's `REVIEW.md`. A
   worker who fixes one adds the tested revision, the replay command and the result under its identifier and marks
   it **Ready for re-audit**. The owner requests the re-audit. Only Codex marks a finding **Verified** and strikes it
   through after checking the fix, and it records any new finding.
4. **Mark each row.** The Codex auditor marks a row **SUBMITTABLE** only when the packet requests ✅, the evidence
   warrants an expected full mark on the entire official criterion on the named candidate, all applicable checks from
   A1 to F3 pass, and the complete G6 audit binds the exact version. A partial result, a registration-only proposal,
   an unresolved clause or another blocking finding is **HOLD**, with the evidence that full credit would need. The
   verdict, filing type, claim, candidate pin, version and audit go into the packet's `CHECKLIST.md`.
5. **Clear the package.** A combined pull request needs a submittable verdict for each included row and clearance of
   its shared material, including any pin proposal and cover text. When the claim, the source files, the evidence or
   the outgoing text changes, Codex rechecks the affected items and records a new version-bound verdict.
6. **Release the source and the pull request together.** The owner's final authorization of the cleared version
   makes it **READY TO SEND**. Push the audited commit to this repository, download that commit without credentials,
   verify the packet files against the manifest and run the replay, then open the matching OpenWave pull request
   with source links naming that commit. If public verification fails, the package waits. Record the source and
   pull-request links with the packet, and later the maintainers' decision and the cell revision.
7. **Run the upstream linters.** In the upstream checkout, run `python3 dev_docs/utils/check_roadmaps.py` for
   roadmap changes and `python3 dev_docs/utils/check_models_md.py` for matrix or count changes. Read the outputs,
   check the links and inspect the exact diff against its base ([review standards][prstd]).

## 14. Reporting

- Report the official status first, with its dated decision and source, followed by local findings, supplied
  ingredients and open obligations.
- Report finite proofs, source consistency, operational identification and empirical comparison separately.
- Attribute acceptance only to the maintainers' actual cell decision. Keep task closure, independent recomputation
  and the cell verdict distinct.
- Scope positive and negative findings to the candidate, observable and hypotheses tested, and state which theorem or
  experiment would close the remaining gap.
- Keep every outcome, including failures, timeouts, invalid and complement outcomes and failed alternatives.
  Distinguish inconclusive observations from falsifications and proof claims from finite samples.
- Where an independent check is claimed, state its actual independence, including the producer code it imports.

## 15. Precedents

| Upstream record | What it shows |
| --- | --- |
| [T11][t11] | After the M10.2 review, every ✅ cell in every column was re-read under one written bar. Cells whose runs lacked the property on the model's own dynamics moved, among them M5's Klein–Gordon cell from ⚠️ to 🚧. A cell that moved keeps its result in its text. |
| [M5 particle stability][results] and [M5.21][m5_21] | The ✅ stability cell rests on audited free evolution and names its box, its time window, the energy it absorbs and its non-compact core. A ✅ cell states its limits. |
| Charge and electron cells of [M4 and M8][results] | The M10.2 review cited them as the bar applied elsewhere: M4's charge row is ❌ with its sign "imposed", M8's charge row is 🚧 on an assignment, and M8's electron row is 🚧 with the mass as a "calibration anchor (benchmark input, not a prediction)". |
| [M10.1][m10_1] | A clean-room recomputation confirmed seven frozen coefficients on the declared twelve-direction operator. Arithmetic on a declared object closed the task with no cell moved. |
| [M10.2][m10_2] and its [maintainer review][review] | The author filed a method note, scripts and receipts. The maintainers recomputed the headline numbers from one data file and the method note's formulas without executing the author's package, adopted the pin, and kept rows 1, 2 and 7 at 🚧 because the integer, the anchored value and the half-turn return came from declarations. The review names what would move each row. |

[ow]: https://github.com/openwave-labs/openwave
[models]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md
[bar]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md#summary-status
[simplest]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md#simplest-test-per-criterion
[results]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md#results-by-model
[admission]: https://github.com/openwave-labs/openwave/discussions/712#discussioncomment-19595538
[canonical]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_theory_canonical.md
[canon1]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#1-specification-of-record-and-the-pin
[canon8]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#8-consumption-rules-how-this-column-is-scored
[orientation3]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_agent_orientation.md#3-the-working-rules-this-column-runs-under
[roadmap]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/m10_roadmap.md
[review]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/findings/m10_2_maintainer_review.md
[m10_1]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/tasks/m10_1_task_details.md
[m10_2]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m10_oph/research/tasks/m10_2_task_details.md
[m5_21]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/openwave/xperiments/m5_liquid_crystal/research/tasks/m5_21_task_details.md
[t11]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/tasks/t11_task_details.md
[hygiene]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/AI_HYGIENE.md#1-the-stance
[onboarding]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/ONBOARDING_MODELS.md#11-the-one-question-that-matters-prediction-or-post-fit
[paramcount]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/ONBOARDING_MODELS.md#the-parameter-count-test
[prstd]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/PR_REVIEW_STANDARDS.md
[gatee]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/PR_REVIEW_STANDARDS.md#7-gate-e-modelsmd-cell-changes
[safety]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/PR_REVIEW_STANDARDS.md#31-the-safety-gate-before-anything-runs
[methodnote]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/METHOD_NOTE.md
[roadmapstd]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/dev_docs/ROADMAP_STANDARDS.md
[reproduce]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/REPRODUCE.md
[contributing]: https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/CONTRIBUTING.md
