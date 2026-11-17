# Cautions for OpenWave research and submissions

Recorded 2026-11-15 from the review of PR #723 and our submission mistakes.
These cautions apply across all rows. Agent obligations are in
[AGENTS.md](AGENTS.md); every submission must fulfill the
[mandatory submission checklist](SUBMISSION_CHECKLIST.md).

The central failure was treating reproducible measurements through the native
process as sufficient evidence that the process produced the claimed physical
property. Our local review established arithmetic and execution while leaving
that inference inadequately challenged. Disclosing the assumptions did not
make the unsupported property eligible for partial credit.

**Do not claim that source dynamics produces a property when the evidence only
shows a consequence of something supplied by construction.**

Every claimed scientific test must answer:
**What possible run outcome would falsify each claim?**
Do not call a run a test of that claim without a concrete answer: an outcome
representable in the stated test domain, a way to observe it and a rejection
condition fixed before the run. The claimed law may predict that this outcome
never occurs; the representation must not exclude it by definition.
If every record the apparatus can represent
passes by construction, the run supplies no test of the claimed emergence.
A corrupted file or an out-of-domain signal is not a substitute. Distinguish
an exact counterexample from statistical rejection and an inconclusive timeout.

## 1. Do not build in the answer and call it emergence

Do not count a property guaranteed by the state alphabet, geometry, preparation,
reference, normalization or decoder as a dynamical discovery. An integer-valued
register does not establish physical quantization. A return period forced by
symmetry does not establish a dynamical rotation mechanism. Rewriting the same
bookkeeping as a topological statement supplies no independent explanation.

Declared architecture and apparatus are permitted model choices. Their presence
is not itself a failure, and choosing them does not require a uniqueness proof.
The error is crediting the target property to the dynamics when its evidence
comes entirely from those choices.

## 2. Do not substitute execution for an explanation of the property

Calling the native event function, accumulating many records and obtaining a
nontrivial number do not establish every physical claim attached to the run.
Do not use evidence for one behavior to support a different headline. A law
can affect a magnitude while the claimed symmetry or period is fixed elsewhere.
More events, exact arithmetic and a passing replay do not repair that mismatch.

Do not claim the process selected a state merely because a preparer supplied
it or a spectral certificate recognizes it. A prepared state, an eigenvector,
a stationary distribution and dynamical selection are distinct claims.

## 3. Do not let a calibration validate itself

Do not present a fitted value as a prediction or count reproduction of its own
anchor as an independent success. Do not hide choices of scale, reference,
state, clock or normalization that were informed by the target value.

Do not quietly recalibrate, change a convention or select another state after
an unfavorable comparison and describe the result as a successful test of the
original claim. A disclosed calibration still supplies no evidence for the
value it was chosen to reproduce.

## 4. Do not assign physical meaning by naming a mathematical quantity

Calling a statistic energy, a label charge or a representation spin does not
establish the corresponding physical interpretation. Do not equate a decay
rate, transfer spectrum, record cost or modular quantity with physical energy
without the operational argument required by the claim.

Do not present a separately programmed reference signal as a measurement of
the source when the source's actual state or accepted measurement record never
controls that signal. Agreement with a supplied parameter is insufficient.
Do not resolve an ambiguous phase or other aliased reading by supplying the
answer that the measurement is meant to establish.

Do not count reading back a programmed coefficient as an independent derivation
of that coefficient, even when the source controls the gate. The derivation of
the source operator, its faithful implementation and an independent use of its
output are different obligations.

Do not mix event counts, reconstructed time and physical clocks without an
explicit justified relation. Do not dismiss unexplained composition effects,
such as excess energy assigned to separated noninteracting copies, because the
formula producing them is mathematically correct. Do not use hidden simulator
state as if it were a measurement available in the observer's records.

## 5. Do not use controls that miss the claimed mechanism

A control that breaks geometry or changes the reader cannot by itself isolate
the contribution of the native law. A changed law that differs only outside the
tested state domain does not challenge the tested property. Do not count a
control's failure on one quantity as evidence for another.

Do not substitute corrupted files, impossible source states or deliberately
inserted signals for a scientific alternative. They may test integrity or
sensitivity; they do not establish that the mechanism generates the signal.
Do not interpret a null from a detector that misses a planted effect.

A result shared by another law can be universal. Do not claim that such sharing
alone disproves the result, or that the original experiment has isolated the
proposed mechanism when the relevant alternative was never distinguished.

Do not credit a selected symmetry to the dynamics when the replacement state
already has that symmetry and the only control merely disables replacement.

## 6. Do not manufacture agreement by discarding the distinction

Do not infer a physical identification from a reader that removes the very
information at issue. Squaring a vector, taking an absolute value, dropping a
sign, averaging away orientation or collapsing every state to the same outcome
can force agreement without the claimed dynamics.

Do not credit an imposed counter period, representation sign or planted phase
as a property generated by the source. A convenient projection cannot establish
what the discarded records would have shown.

## 7. Do not transfer results between different candidates

Unchanged endpoint code does not make two experiments the same model. State
domain, preparation, scheduler, controller, memory, read/reset operations,
decoder and outcome alphabet also belong to the candidate.

Do not combine successful clauses from different candidates into one accepted
model without an implemented and tested common-source connection. Do not
silently change the candidate to repair an objection, substitute another tree
for the reviewed commit, or overwrite sealed evidence. Ancestry does not
transfer validation to a descendant.

## 8. Do not strip a theorem or result of its hypotheses

Do not promote a bounded-density result to arbitrary loads, a finite sample to
an all-input theorem, or a conditional construction to an established physical
fact. Do not generalize a negative beyond its source, apparatus and observable.
A finite timeout does not automatically refute an asymptotic claim.

Do not credit classical branch averages with coherent-memory dependence merely
because the package contains a coherent component. Do not borrow a certificate,
proof or physical interpretation from another domain without establishing that
its hypotheses apply to this experiment.

Do not substitute a small population outside a subspace for a small distance
between full coherent states. Cross terms and the actual input domain matter.
Do not silently replace a recorded source state with a convenient prepared
state when joining two stages, or mistake an encoding identity for an executed
preparation and complete acceptance instrument.

Do not infer agreement of a joint state or joint observable from agreement of
its marginals. Correlations, a common physical transformation and the full
outcome law must be checked. Independently chosen coordinate relabelings do
not establish one spatial rotation of the combined object.

Do not call a preparation assumption removed when it still enters through
the census, reference, controller or decoder. Distinguish the units stabilized
by the source law from the number of units a controller chooses to retain.

## 9. Do not conceal unfavorable or inconclusive outcomes

Do not drop failed trials, timeouts, invalid or complement outcomes, lower
modes, background ambiguities or unsuccessful alternatives. Do not count
inconclusive measurements as successes or change the denominator after seeing
the results. Do not remove a known counterexample from later summaries.

A valid negative result is evidence. An unreliable detector, incomplete outcome
record or unverifiable run is a defect in the evidence. Do not confuse these
cases to claim success or to suppress a legitimate negative result.

Do not omit failures incurred before the final detector when composing
instruments. Do not report conditional readings without their acceptance
probability and propagated error. Machine-precision agreement is not the
physical uncertainty when the complete implementation error is larger.

Do not label signed numerical comparison values as physical probabilities.
Do not hide trace drift by clipping or renormalizing it, assume an approximate
map has exact unit gain, or turn a small numerical residual into a stronger
claim than the certified complete instrument supports.

## 10. Do not register discoveries retrospectively

Gates written before replaying an already known result do not make that result
a prospective prediction. Do not hide exploratory runs, target-informed
choices or prior knowledge. File timestamps alone do not establish the claimed
order of design, registration and execution.

Do not revise a failed gate, erase the earlier version and report the amended
run as if the revised terms had always applied. Reproduction and prospective
testing support different claims.

## 11. Do not substitute a large archive for reviewable evidence

A small wrapper that fetches and executes a large package is not a small
execution payload. Do not assume a reviewer has inspected or executed every
file merely because the archive is available or its hashes match.

Do not count events represented by enumeration, matrix powers or compressed
sampling as native calls unless the native function actually executed them.
Do not inflate actual call counts with represented possibilities. An
independently reproduced headline does not verify every script, proof,
execution count or gate in its source package.

## 12. Do not inflate independence, verification or acceptance

A second worker from the producing group remains a same-group reviewer, even
when its implementation is independent. That check is useful, but do not
report it as other-group confirmation. A component review does not pass the
full-packet review gate.

Do not equate arithmetic agreement with physical validation, or a merged PR,
adopted pin or approved record with a criterion award. Do not claim a proof,
payload or gate was checked when the reviewer explicitly excluded it. Do not
report local expectations or stale summaries as the current official decision.

## 13. Do not reinterpret the criterion to excuse missing evidence

Do not use the word static, another column's award or a narrow reading of a row
label to bypass the column's own-process requirement. Do not treat a disclosed
assumption as automatic grounds for partial credit when it supplies the very
property being claimed.

Do not turn an unresolved objection into a harmless caveat by renaming it a
physical-dictionary question. Equally, do not invent extra obligations beyond
the claimed criterion or require every physical property for a narrower claim.
The problem is a claim exceeding its evidence, not the existence of limited
but useful results.

## 14. Do not let polished packaging override a scientific blocker

Passing formatting, custody, replay and arithmetic checks does not establish
submission readiness. Do not approve a packet without challenging the causal
and physical meaning of its headline. Do not treat a detailed response, new
implementation or larger evidence bundle as proof that a criticism was resolved.

Do not submit an unsupported full or partial claim hoping the maintainer will
supply the missing argument. Do not bypass a failed checklist item with softer
wording, a review-request label or an unexplained N/A. Documentation of a
problem does not satisfy the scientific obligation it describes.

## Sources

These general cautions draw on the
[maintainer review of PR #723](row_01_charge_quantization/reviews/2026-11-15_codex_pr723_integration/UPSTREAM_REVIEW.md),
the [local static adjudication](shared/compact_static_adjudication/REVIEW.md),
the [final submission review](notes/static_submission_audit_20261115/FINAL_REVIEW.md)
and the [follow-up integration review design](row_01_charge_quantization/reviews/2026-11-15_codex_pr723_integration/DESIGN.md).
The original evidence and issue-specific research remain in those records and
the row folders. The [framework](docs/FRAMEWORK.md) and
[submission checklist](SUBMISSION_CHECKLIST.md) govern the required workflow.
