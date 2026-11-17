# STRATEGY: from finite theorems to accepted rows

Written on 2026-11-14 after the collection of the finite theorems of all 31 rows (`row_NN_*/THEOREMS.md`,
machine form `premises/TA.json` to `TD.json`). Read with [ESTABLISHED.md](ESTABLISHED.md) and
[DELIVERABLE.md](DELIVERABLE.md).

## 1. What the collection shows

Every row has a chain of statements from the observers' process to the row's effective description. The
machine-checked and written theorems carry most of each chain. Every chain stops at a small set of open
premises, and the premises are shared: 31 rows reduce to about a dozen premise groups.

The open premises are of three kinds:

| Kind | Meaning | Who closes it |
| --- | --- | --- |
| **Decision** | a design choice of the model's law that no theorem about committed records can supply | the owner |
| **Theorem** | a finite statement on a declared object, provable once the decision is made | Lean or a written proof, reviewed |
| **Execution** | a fit of the effective description to observers' records on the decided candidate | a row packet ([DELIVERABLE.md](DELIVERABLE.md)) |

The central finding: **the decisions come first, and they cannot be replaced by theorems about the current
law.** `ContinuousPublicStarFlow.toAut_eq_refl` (`Lean/Dynamics/PublicAutomorphism.lean`) proves that a
continuous flow on the committed record algebra is the identity, and `COHERENT_CONTINUATION.md` shows two
continuations that agree on every committed branch. So no theorem derives the coherent continuation (PR-21,
PR-43) from committed records. The uniform-scheduler paired-unit theorem fixes I − L/60 in every background, but does not
exclude every native response: the pinned feedback controller has a nonadditive second-current contrast.
The native convex-averaging rule bounds signed-load pair creation; spin, binding and parity tests must
retain their precise apparatus, observable and protocol hypotheses. Each
of these needs an ingredient of the law that observers' records then reflect. Choosing those ingredients is
the specification of OPH's observer, which the application itself describes: "an observer is a feedback loop
that reads data ... its next action depends on that record."

## 2. The open premises, ranked by rows blocked

| Premise group | Register ids | Rows blocked | Kind |
| --- | --- | ---: | --- |
| Schedule and coin: is the next interaction drawn from records (feedback), or uniform | PR-05, PR-20 | 12 | decision, then theorem |
| Energy and clock: one operational energy and clock; the maintenance-to-energy map | PR-15, PR-22, unregistered map | 14 | theorem and execution (PR-15 stays the admitted 511 keV anchor) |
| Gauge attachment: currents and fields as records | PR-54 | 14 | theorem and execution |
| Species: the event reads colour class, isospin, generation | U-TC-READ, PR-59, PR-55, PR-67 | 7 to 13 | decision |
| Coherent continuation on that specific route: unread register reuse and redraw | PR-21, PR-43, reuse, blind redraw | route-dependent | optional construction; native positive-history reconstruction supplies a separate free-wave route |
| Phase preparation and readout | PR-64, PR-65 | 6 | execution |
| Spin: a register that carries the ½ sign is written by the law | PR-47, area counter | 5 | decision, then the theorem that the doublet is the law's slowest mode |
| Third direction | PR-52, PR-16, PR-53, routing | 3 to 8 | decision, then a theorem bounding the error against 1/(4πr) |
| Einstein hypotheses: one source, stationarity with one η, area at fixed volume, Ward identity | PR-23 to PR-28 | 3 | theorem and execution on the decided law |
| Pairs and chirality: events that create pairs; a handed ingredient | PR-47, PR-67, U-TC-PAIR, U-TC-CHIRAL | 6 | decision |
| Capacity: a law that selects N | PR-17, PR-51, PR-42 | 1 | decision (local data cannot fix N: Lean `CapacityNonidentifiability`) |
| Calibrations to be replaced | U-TC-HADCAL, PR-37, PR-73 | 5 | they stay out of every cell; replaced by executions |

## 3. The plan, in three stages

**Stage A: one candidate specification.** Local candidate construction is authorized by the resumed research instruction. Selecting the column pin remains the owner’s act. Write the candidate "OPH observer law v1" in the platform's four
words, with each ingredient chosen and the rows it opens:

| # | Choice | Options on file | Opens |
| --- | --- | --- | --- |
| A1 | Optional coherent-memory route: choice register retained and blind redraw | candidate D′, if chosen | one route for phases/compact states; not a necessary primitive for reconstructed waves |
| A2 | Schedule: the next seam is drawn from the observer's committed letters (the pinned controller, continued to a screen) | yes (candidate A+c) / uniform | any response: rows 4, 8, 12 to 20, 22, 23, 26, 27 |
| A3 | Spin register: face events write the area counter | yes (candidate J) / no | rows 5, 6, 7 |
| A4 | Species: the event reads the Z6 class (and isospin, generation) | which registers | rows 9 to 18, 22, 24, 25 |
| A5 | Pairs and handedness of the event | which ingredient | rows 8, 9, 23, 24, 25 |
| A6 | Third direction: the cubic routing joins the candidate | yes / another routing | rows 20, 21, 22, 23, 26, 27, 29, 31 |
| A7 | Native layered-read implementation and operational interval-count clock | existing finite compiler plus integer-event adapter; declared clock comparison | row 21 and geometric clock comparisons |
| A8 | Native full-moment field reconstruction | registered MFR-common-native-field-v1, occupation-only, periodic cubic refinement and same-tape first-hit clock | mathematical Maxwell/Klein-Gordon clauses; possible source orbital extension |

Each choice is a pin move at OpenWave, stated once for the whole column. A choice can be "not yet": the
rows it opens wait.

**Stage B: the theorems on the decided law** (each a Lean file or a reviewed written proof, as a shared
study in `openwave/shared/`):

1. The doublet of the area counter is the unique slowest mode of the law (rows 5, 6, 7).
2. The coherent walk on the decided memory law: dispersion, compact states, the mass angle from the movement
   excess, 1 − cos ω₀ = ε (rows 2, 3, 30), and the front speed 1/√48 against the field speed (rows 29, 30).
3. Under A2: the response kernel depends on the background (the converse of E-21), with its sign (rows 20, 26).
4. Faraday kinematics: a face register written by events, with an incidence theorem (row 29).
5. Attach the already-proved golden-family diamond-count limit to the native layered-read program and operational count clock (row 21). Poisson scheduling is not a prerequisite of that theorem.
6. The error of the cubic routing's potential against 1/(4πr), with a bound (rows 20, 31).
7. On the decided species law: a neutral triple lowers record cost; two charge states differ (rows 13, 14).
8. The Einstein hypotheses PR-24 to PR-27 on the decided law, one at a time (rows 26 to 28).
9. Keep every negative within its proved source and observable domain. Primitive one-unit or signed-load
   statements do not exclude decoded excitations. Uniform stationary weights do not exclude source spectral
   binding. Parity covariance must include the complete observer/controller. Positive exchange, relative-source
   and reconstructed field results in ESTABLISHED supersede broader versions of those claims.

**Stage C: one execution per row** on the decided candidate: the fit of the row's description to observers'
records, as a row packet. With the theorems of stage B the fit is predicted in advance, so each packet's gates
are numbers fixed by a theorem, and a failed fit is an informative negative.

## 4. Order

1. Specify and test local architectures concretely (A1, A2 first); present an evidenced pin proposal when one exists. Do not wait for a uniqueness proof of architectural choices.
2. While the owner decides: the exact negatives of B9, which need no decision and close dead ends; and the
   row-1 and row-20 packets on the current law, which need only A6.
3. After each decision: its stage B theorems, then the stage C packets of the rows it opens, in the order of
   [DAG.md](DAG.md).

## Source-interface correction from the current row loop

The finite theorem stack must be consumed before proposing a missing theorem.
The golden-family Lorentz count limit, source-native stored-program compiler,
full reconstructed matrix algebra and controlled modular histories already
exist in RER. New work should prove the specific implementation/observation
transfer they still need. The RER construction contract permits declared
architecture, scheduler and apparatus; proving their uniqueness is not an
extra prerequisite. A new native program remains an unpinned candidate until
the owner selects it. None of these finite components alone establishes the
full physical rows or supplies missing empirical strengths.

## Construction and theorem application in the resumed work

The actual RER `AGENTS.md`, construction contract and hypothesis page allow the
repair law, code, memory, scheduler and observer apparatus to be declared design
choices. Proving their uniqueness is not a prerequisite for construction. The
physical descriptions and target numbers still have to follow from the resulting
observer records; declaring a separate target Hamiltonian is not an implementation
proof. Research authorization permits concrete candidate experiments before a
column pin move. An unselected candidate remains unpinned.

Before proposing a missing theorem, read its row register and the named Lean
statement, hypotheses and source binding. The finite record-energy and recurrence
clock constructions are now linked from rows 2 and 3. In geometry,
`GoldenSourceCausalLimit.generated_eventually_iff` and
`generated_interval_count_tendsto` already prove order/count limits for the
declared deterministic golden population and complete-neighbour read law.
`FlatLorentzVolume.tilted_diamond_volume` and
`SourceCountClockEnclosure.finite_clock_enclosure_counts` supply the Lorentz
diamond formula and finite clock bounds. A Poisson schedule is not a general
necessity theorem for every refinement limit. Applying this stack requires its
actual population/read-law/count-error hypotheses and clock comparison, rather
than a blanket statement either that randomness proves Lorentz symmetry or that
the repository has no relevant derivation. Row 21 remains with its assigned worker.

## Common source free-wave route

E-61, E-64 and E-75–76 establish a constructive route from native currents to full Gaussian moments, two transverse Maxwell states, native massive scalar bands, locality and a common limiting speed. The field time is reconstructed; primitive endpoint attempts need not be unitary or future-directed in that effective spacetime. The explicit common member uses occupation-only targets L+−Pi/(4W), charges the contact term separately, and takes cubic equal-period refinements. Free-wave mathematical clauses are proved; owner pin, other-group review and in-platform execution remain pending. Species, interactions and the Einstein common-source interface remain separate tasks.
