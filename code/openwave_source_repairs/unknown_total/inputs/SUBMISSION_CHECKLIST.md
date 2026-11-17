# Mandatory submission checklist

**No further upstream submission may be sent until its applicable checks below
are fulfilled and recorded against the exact material being submitted.** This
is the owner's local submission requirement, added 2026-11-15. It supplements
[CAUTIONS.md](CAUTIONS.md) and the [framework](docs/FRAMEWORK.md); passing it
does not guarantee a maintainer award.

## How to complete it

Keep this master unchecked. Make a completed copy in the submission's review
directory, outside the sealed scientific packet. For every item ID, record
`PASS`, `FAIL`, `OPEN` or `N/A`, the evidence path/commit and result, the
reviewer's name/group, and the review date. A checked box without that evidence
is incomplete. Use a separate entry for each row or materially different claim.

`FAIL`, `OPEN`, missing evidence or an unresolved blocking review finding means
**DO NOT SUBMIT**. Fix the problem or remove the unsupported claim and review
the revised submission. Disclosing a failure does not turn it into a pass.
Any change to the bound source, data, gates, interpretation or submission diff
requires a new version of this record and renewed affected checks.

For a full or partial row-award request, every scientific check applicable to
the claimed clauses is mandatory. `N/A` requires a specific scope reason and
reviewer agreement; lack of evidence is not a scope reason. Pure task
registration, administrative corrections or explicitly negative/evidence-only
reports may omit inapplicable award checks, but must request no unsupported
cell movement and must retain the exact status and limits. In particular,
registration before execution must not require future results to exist.

A scientifically valid negative result may be reported as a negative. It does
not excuse an unreliable detector, missing outcomes or an unverifiable run.
Mark checks about demonstrating a positive property `N/A` only when that
property is explicitly not claimed. Calling a submission a review request
does not exempt an implicit full/partial-award or readiness claim.

Record the submission identity before checking anything. Fields outside the
declared scope may be `N/A` with a reason: an administrative correction does
not need a scientific packet, candidate or scientific reviewer.

| Field | Required entry |
| --- | --- |
| Task, title and type | Full award / partial award / evidence-only / registration / administrative |
| Rows and exact claims | Each requested clause, proposed mark and current dated upstream mark |
| Candidate | ID, immutable source commit and any proposed pin change |
| Evidence | Packet paths and manifest SHA-256 values |
| Submission | Destination, base/head commits or exact diff hash, draft PR text path/hash |
| People | Producer, scientific reviewer and group, integrator, owner |
| Review record | Date/time with timezone and link to this completed checklist |

## A. The claim meets the actual criterion

- [ ] **A1 — Terms read together.** The exact row clauses, M10 admission terms
  and latest applicable maintainer decision are cited with their versions.
  “Static” and another column's precedent do not waive the own-process rule.
- [ ] **A2 — Evidence per clause.** Each claimed clause has an evidence-linked
  `full / weaker / absent / negative` assessment. The requested mark follows
  that assessment. A partial request identifies the property that actually
  follows from the process and the separately supplied ingredients.
- [ ] **A3 — Law contribution established.** For each claimed property, the
  record states what the native dynamics contributes and demonstrates that
  the headline is not guaranteed by the alphabet, preparation, geometry,
  calibration or decoder. Native execution alone does not satisfy this item.
- [ ] **A4 — Operational interpretation supported.** The chain from native
  process to observer records, observable and physical interpretation is
  explicit and checked. Hidden simulator state is confined to labelled audits.
  Calling a quantity charge, energy or spin is supported by the claimed tests.
  A reference signal is coupled to the actual source or its retained measured
  outcome; source-independent calibration alone is insufficient.
- [ ] **A5 — Existing work and limits retained.** Relevant established results,
  theorem hypotheses, prior failures and dependencies have been read and cited.
  Neither positive nor negative conclusions exceed the tested domain.

## B. The headline can fail for the right reason

- [ ] **B1 — Falsifying outcome specified for every claim.** Answer explicitly:
  **“What possible run outcome would falsify each claim?”** Before the claimed
  run, record each scientific claim's counterexample or outcome representable in the test
  domain, the measured quantity and detector, and the exact rejection condition
  or statistical threshold. Link evidence that the detector can register the
  relevant effect. Distinguish proof counterexamples, statistical rejection and
  inconclusive timeouts. A corrupted-input failure or a gate guaranteed to pass
  by construction does not fulfill this item for a claim of dynamical emergence.
- [ ] **B2 — Relevant alternative law tested.** The same observer and readout
  test an appropriate changed law on the claimed object. The result separates
  the target property from the competing law, with relevant geometry/symmetry
  preserved. A changed magnitude alone does not establish a return-period claim.
- [ ] **B3 — Detector and identities audited.** Planted effects are detected;
  integrity mutations and identities are listed separately. An impossible
  half-charge, planted phase or corrupted file receives no emergence credit.
  Aliased readings are resolved without inserting the claimed answer.
- [ ] **B4 — Coherence dependence stated.** Each claim identifies whether it
  uses classical branch statistics, coherent amplitudes or a conditional
  reconstruction, with supporting calculation/control. Coherence is not
  credited merely because another component in the package uses it.

## C. Calibration and the PR #723 objections are resolved at the claimed scope

- [ ] **C1 — Inputs do not validate themselves.** All measured constants,
  fitted choices, scales and selection decisions are listed. Every calibrated
  claim identifies the independent test left after fitting and reports its
  outcome, including misses. Clocks, units and composition conventions agree.
  Conditional acceptance and complete implementation error are propagated into
  the reported observable interval; floating-point residuals are not substituted
  for that interval.
- [ ] **C2 — Row 1, if claimed.** Evidence follows a review-supported route to
  law-selected units or a genuine configuration invariant. Supplied integrality
  and equivalent load-vector topology earn no quantization credit. Density and
  background limits are explicit; counterexamples such as uniform-two remain.
- [ ] **C3 — Row 2, if claimed.** State selection is attributable to the
  dynamics, beyond the preparer's choice, with holding or return demonstrated.
  The anchored energy passes either separated-copy additivity or an unanchored
  prediction compared with measurement. Lower
  modes and the old approximately 2.45085 keV pair excess have evidence-linked
  dispositions. A transfer eigenray or additive record cost alone is insufficient.
- [ ] **C4 — Row 7, if claimed.** Evidence defeats the symmetry-preserving
  alternative. The supported route is either dynamical antipodal identification
  retaining directional information, or a full-turn sign/two-turn return
  measured by interference with finer turns. Even readers, vacuum erasure and
  planted phases are excluded as explanations of the claimed success.
- [ ] **C5 — Review questions answered.** Every applicable prior criticism has
  a disposition and evidence. PR #723 follow-ups answer all three questions:
  possible falsifiers, actual coherence dependence, and the separated-pair
  excess. Any remaining issue is outside the requested claim and explicitly
  scoped; a blocker cannot be relabelled as a caveat to obtain partial credit.

## D. One candidate and a traceable experiment

- [ ] **D1 — Complete candidate bound.** Law, state, preparation, scheduler,
  controller, memory, read/reset operations, decoder and outcome alphabet are
  specified. Changed protocols are registered as new candidates with ancestry.
- [ ] **D2 — No borrowed successes.** All claimed clauses run on the same
  candidate, with any common-source join implemented and tested. The source
  commit is retrievable and preserved; a different main-branch tree cannot
  replace it. Any pin proposal has the owner's recorded decision.
- [ ] **D3 — Timing is genuine.** Design, code, inputs, seeds, controls and
  numeric gates were frozen before the claimed confirmatory run, with ledger
  evidence and the applicable task registration. Earlier pilots and discoveries
  are disclosed. Historical reproduction is labelled as such, never backdated
  into a prospective test; amendments preserve the original record.
- [ ] **D4 — Execution matches the claim.** Receipts bind the frozen source and
  inputs. Literal native calls, enumeration, matrix powers and sampled histories
  are distinguished; any compression has a checked equivalence argument.
- [ ] **D5 — Full outcomes retained.** Failures, timeouts, invalid/complement
  outcomes and failed alternatives remain available. Statistical claims state
  uncertainty, selection and power; finite samples do not establish universal
  claims without the required proof.
  Composed instruments retain upstream as well as terminal failures. Numerical
  comparison centers and trace drift are distinguished from positive physical
  outcome laws, with gain and conditioning bounds where required.

## E. A reviewer can inspect and reproduce the submitted evidence

- [ ] **E1 — Small offline check.** A reviewer can read the complete primary
  execution path and bounded inputs. The headline reproduces in under 60
  seconds under the local framework, with command, environment, runtime and
  exit status recorded. A download-and-execute wrapper does not satisfy this.
- [ ] **E2 — Clean replay and custody.** A clean copy at the submitted pin
  reproduces the reported results; manifest verification passes. Portability
  accommodations and provenance limitations are explicit and reviewed.
- [ ] **E3 — Scientific adversarial review.** A reviewer challenges the meaning
  of the headline and its competing-law test, beyond checking arithmetic and
  packaging. The review binds the packet manifest and final claim wording.
- [ ] **E4 — Required independent review.** For an award request, G0–G8 meet
  the framework's confirmation requirements, including a full other-group G6
  review with no blocking findings. Same-group Codex checks and external checks
  of individual components do not satisfy full G6. Partial-ready uses G8's
  clause-level rule; it still requires a law-derived property.
- [ ] **E5 — Gate state checked.** Packet files, manifests, `STATUS.json` and
  review records agree. Run `python3 tools/gate.py check` from `openwave/` and
  attach the result. Resolve submission-related violations; separately identify
  any unrelated workspace diagnostics. A clean checker does not replace the
  scientific checks above. G9 remains pending the actual upstream decision.

## F. Final submission and release

- [ ] **F1 — Text matches the evidence.** Task, method note, receipt, proposed
  cell and PR description agree on candidate, result, supplied ingredients,
  limitations and requested mark. The cell meets the local 65-word limit.
  Acceptance and independent verification are attributed only at their actual
  scope; no merge, replay or candidate adoption is presented as a row award.
- [ ] **F2 — Exact outgoing diff reviewed.** The final files and links are
  checked against the intended upstream base, whitespace/path checks pass and
  the current contribution requirements are met. Supporting data and source
  pins are accessible to the reviewer. Concurrent work is preserved.
- [ ] **F3 — Completed record signed.** Every applicable item has evidence and
  reviewer confirmation, all `N/A` reasons are approved, and no `FAIL`/`OPEN`
  item remains. Producer and reviewer identify the exact reviewed packet and
  outgoing diff. The completed checklist is included or linked in the draft.
- [ ] **F4 — Owner authorization covers this version.** Before release, record
  the existing owner authorization covering the reviewed submission and action
  under the upstream workflow. Do not infer it from permission to conduct
  research. Request authorization only if the reviewed action is not already
  covered, after all preceding checks and reviewable material are complete.

## Completion record

Fill one row per item and per affected claim; the first row below is a template.

| Item / row / claim | Status | Evidence and observed result | Reviewer / group / date |
| --- | --- | --- | --- |
| A1 / … | OPEN | Required source and comparison pending | Pending |

| Decision | Required record |
| --- | --- |
| Scientific review | Reviewer, group, bound manifest, verdict, unresolved findings |
| Final diff review | Reviewer, exact diff/commit and draft-text hash, date |
| Owner release | Authorization reference, exact version, date |
| Submission verdict | **DO NOT SUBMIT** until every applicable check above is fulfilled |

Keep G9 and official cell awards separate from this pre-submission record.
Record them only after the maintainers' actual decision.
