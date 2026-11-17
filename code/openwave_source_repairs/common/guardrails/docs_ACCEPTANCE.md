# ACCEPTANCE: what OpenWave accepts, and the evidence a row has to produce

Two questions, answered exactly. The quotations are from the platform's own files as read at commit
`5d686cb049b8` (copies in `openwave-push/results_nov14/platform/`). Where this file and a platform file
disagree, the platform file holds.

## 1. What counts as OpenWave acceptance

**A row is accepted when the maintainers have set its cell to ✅.** They decide "at review, one cell at a
time, from the artifact" (`m10_theory_canonical.md` § 8). Nothing that an agent of ours writes is an
acceptance. The platform's four marks:

| Mark | The platform's words | What earns it in column M10 (§ 8) |
| :---: | --- | --- |
| ✅ | validated in-platform | "With the source derived from the carrier, it meets the bar for ✅". The bar is the maintainer's answer to issue 482: "derive the encodings in-platform. [...] When those land, the exponents are measured rather than encoded, and the row flips." |
| ⚠️ | partial / with caveats | "Its observable computed on the carrier's own process in-platform. With something still supplied, that is a case for ⚠️ with the caveat in the cell." |
| ❌ | negative | A runnable script that shows the criterion is not met, with the reason; it "merges on the same evidence bar as a positive". |
| 🚧 | planned / not tested | "A row whose construction runs on a source that is declared, not derived from the carrier, stays at 🚧." Every row of OPH today. |

Conditions that hold for every mark:

1. **In-platform.** The observable is computed by a script that the maintainers can run, registered as a task
   of the column. "The platform's claims are the cells in `MODELS.md` and the scripts behind them."
2. **One candidate at a pinned commit.** "The column scores the executable model at a pinned commit. A changed
   law, memory, state or decoder is a new candidate, with its empirical ancestry kept visible. [...] Moving
   the pin is a pull request that says what changed."
3. **Terms before runs.** "A result that will claim a cell files its gates and inputs first."
4. **What a pull request brings.** "The method note and a small script that regenerates the headline number.
   The run record stays in the author's repository at a pinned commit, and nothing is imported in bulk."
5. **Calibrations.** A calibrated source "is labelled as one and earns no credit for the values it was fitted
   to. [...] Its held-out failures go in the cell as results."
6. **Strengths.** "The Coulomb, Newton and cosmology rows credit a strength only when it follows from the
   model's own mechanism, never fitted or supplied" (rows 20, 26, 28).
7. **The author's own counts are estimates.** "Statuses on the author's own site or papers are the author's
   estimates."

In the terms of OPH (observers first): "the carrier's own process" is the repair-and-memory law that the
native observers run; "computed on it" means computed from the records those observers hold; "the source
derived" means that the description named in the criterion is measured from those records, where today it
is written in.

**The bar for ✅, row by row, therefore has four parts.**

| Part | Statement | Fails when |
| --- | --- | --- |
| A. The process | The run is the pinned candidate: event law, memory law, schedule, preparation, decoder, all registered and pinned | a field law, a coherent step, a register or a schedule is used that the pinned candidate does not contain |
| B. The observable | It is a function of observers' records alone and depends on the law | it is an identity of the alphabet, the graph or the reader; or it is read from the simulator's state |
| C. The criterion | Every clause of the row's simplest test holds in the full sense, on one object that also carries the clauses of the rows it depends on | a clause holds in a weaker sense, on another candidate, or with a value put in |
| D. The procedure | Terms filed before the run; method note; regenerating script; adversarial audit recorded; cell text of at most 65 words | any is missing |

A row with A, B and D met and C met for some clauses, with every supplied piece named, is a case for ⚠️.
A row measured under A, B and D whose records do not show the criterion is a ❌. A row whose result stands on
an installed source is 🚧 whatever the result.

## 2. The evidence we produce

One packet per row, in the row's `work/` directory, with these files and nothing claimed beyond them
(gates: [FRAMEWORK.md](FRAMEWORK.md) section 4).

| File | Content | Gate |
| --- | --- | --- |
| `DESIGN.md` | The clause table of section 3 below, filled for this row | G0 |
| `TASK.md`, `PREREG.txt` | Candidate (four words), observers, records, observable, description to fit, gates with numbers, controls, inputs, development and held-out seeds, every value of nature that enters (or: none); sha256 and time, posted to the ledger before the run | G1 |
| scripts, `RECEIPT.json` | The run on the filed inputs; every gate with its number; every outcome kept | G2 |
| `CONTROLS.md` | The changed-law control on the row's own object, the planted-effect control, the list of statements that are identities of a reader | G3 |
| `headline.py`, `MANIFEST.sha256` | The headline number in under 60 seconds; the seal | G4 |
| `METHOD_NOTE.md`, `CELL.md` | The platform's method note; the cell text of at most 65 words with the expected mark and every supplied piece | G5 |
| `CLAUSES.md`, `PIN.md` | Each clause full, weaker or absent on the candidate; what the pin move changes | G7, G8 |

### The evidence ladder

Every statement in a packet sits on exactly one rung. Only rungs 3 and 4 move a cell.

| Rung | Kind of statement | Example from 2026-11-14 | What it is worth at the platform |
| :---: | --- | --- | --- |
| 0 | A property of an installed description (class F, S, C) | Maxwell waves on an installed quadratic field; a compact state of the reflection walk | nothing: 🚧. Kept as the target that the records should show |
| 1 | An identity of a reader | the integer degree of an integer tape; the period of the area counter | nothing. Marked as an identity |
| 2 | An observation on a candidate that is not pinned, or with a control that does not touch the object | the lattice Coulomb form in first-hit records from a compressed sampler | a component: the basis of the next packet |
| 3 | An observation on the pinned candidate, dependent on the law, with supplied apparatus named | (none yet) | ⚠️ for the clauses it carries |
| 4 | The same with nothing supplied beyond the pinned candidate and every clause in the full sense | (none yet) | ✅ |

### What is never evidence

A name given to a quantity ("energy", "charge", "spin", "mass", "clock"); an anchor used to set a unit; a
count of rows across different candidates; a review by the producer's own group; a gate that cannot fail; a
negative from a detector that misses a planted effect; a statement in a paper without an execution.

## 3. The clause table (what `DESIGN.md` holds for a row)

A row's simplest test is cut into clauses, and for each clause the design states the six items below before
any run. A row is ready when every line is at rung 4; partial-ready when every line is at rung 3 or 4 or is
named as absent in the cell text.

| Column | Content |
| --- | --- |
| Clause | the words of the simplest test, verbatim |
| Observers and records | who holds which registers; what is read, kept, reset; what two observers compare |
| Observable | the function of the records, with its unit in records (events, letters, counts) |
| Description to fit | the effective description the criterion names (a potential, a dispersion, a level ladder), with the gate as a number |
| Control | the changed law that must fail on this object, and the planted effect that must be seen |
| Supplied | every piece outside the pinned candidate, and the sentence of the theory that says where it should come from |

Example, row 20, first clause, as the evidence of 2026-11-14 stands:

| Clause | Observers and records | Observable | Description to fit | Control | Supplied | Rung today |
| --- | --- | --- | --- | --- | --- | :---: |
| "far-field 1/r potential" | a stop observer at a target and the walker's current letters | first-hit frequency at distance d; residence counts by shell | the exact lattice potential of the routing, h(x) = Σ G_B(x,y) e(y), with no parameter; χ² gate | a law whose jump chain is not uniform; a planted exponent | the cubic routing (outside the pin), the continuing schedule, the stop decoder | 2 |
| "at charge e" | (none) | (none: the amplitude is the capacity of the target) | | | | absent |
| "the strength (α) [...] from the model's own mechanism, never fitted" | (none) | (none) | | | | absent |

## 4. The count that is reported

Every status line and every table gives, in this order: the platform's mark today; the mark the evidence on
file would support if filed as it stands; the ceiling of the evidence on file after a pin move and an
in-platform run; the stage in the gates. The README's table is generated from the rows' `STATUS.json` by
`tools/gate.py render`. A count of rows is always a count on one candidate.

## 5. How the column was admitted (2026-11-13)

Source: discussion 712 of the OpenWave repository, saved copy `openwave-tuning/auditor/upstream/discussion_712.json`.

**What the author sent.** One post, "Pre-application question: is OPH admissible, and how are rows on declared
sources scored?" It described OPH observer-first ("Fields, particles and spacetime are the effective
descriptions these observers arrive at when they compare their records"), gave the evidence, and went through
all 31 rows in two groups by the author's own reading of the simplest tests: 24 rows "the constructions
address in full" and 7 "addressed in part" (3, 12, 17, 18, 21, 26, 28). It asked for no cell. It stated: "All
31 rows are 🚧 in-platform, and under your ruling on issue #482 I expect rows on declared or calibrated
sources to earn less than my reading", and that "the row constructions are finite quantum Hamiltonians on
declared lattices" with "one process that carries all rows" as work in progress. It put two questions: does a
model without a field equation qualify; and, for rows on declared sources, "I read your rule as leaving such
rows partial at best, or untested as inherited, until the source is derived from the carrier in-platform. Is
that right, and which of the two applies?"

**What the maintainers answered, the same day.** "OPH is admitted as M10: the column opens at 31 🚧."
Question 1: yes; update rules with no action behind them are tier 2. Question 2: "your reading is right, and
'untested as inherited' is the one that applies. A row on a declared source stays 🚧 until its observable is
computed on the carrier's own process in-platform." Their reasons:

| Reason | Their words |
| --- | --- |
| Precedent | 🚧 inherited: "the source is declared or installed" (M9's Standard Model rows, M8's Einstein equations). ⚠️: "the model's own mechanism runs in-platform" (M4's gravity row). ✅: "the derivation is in-platform" (the answer to issue 482) |
| The installed pieces | "So the Maxwell cochain complex, the Wilson-type Higgs term and the Gauss-law colour string each leave their row at 🚧. The cell says what you built and that its source is declared." |
| The calibrated source | "it earns no credit for the values it was fitted to, and its held-out failures go in the cell as results" |
| The author's own condition | "This is also your own acceptance condition. Your hypothesis page says: 'It is not enough for different supplied actions to reproduce its separate words.' We read that join, one process carrying the rows, as the column's defining open problem." |
| The lepton register | "is a fit, as you count it": two observables against at least five choices; "At N = 302.28 rounds per clock unit the muon ratio equals the measured one" |
| Predictions | the scorecard row "Predictions exist" fails, "as you say, and that does not block the column": OPH is admitted as a runnable system |

No row was marked down: a new column opens with every row at 🚧, and a cell moves only through a registered
in-platform task. Of the two outcomes the author had asked about, the stricter one applies. The first task,
M10.1, a blind recomputation of the twelve-port dispersion coefficients by the maintainers, reproduced the
frozen values: "arithmetic on a declared object; no cell moves". From there "everything after that is yours
to drive". The onboarding follow-up assigned the premise comment, page corrections,
hero image, paper author lists and the first task with its gates to the author.
The [premise/status comment](https://github.com/openwave-labs/openwave/discussions/712#discussioncomment-19624288)
is posted; [PR #724](https://github.com/openwave-labs/openwave/pull/724) submits the
image, bylines and canonical premise summary, and [PR #723](https://github.com/openwave-labs/openwave/pull/723)
submits the first author-run task. Their review state is maintained in the
[upstream tracker](../notes/UPSTREAM_TRACKER.md).
