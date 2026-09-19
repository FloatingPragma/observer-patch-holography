# Issue 777 contract audit — 2026-09-19

This is a self-audit of PR #910, starting from commit `5556050`, with the
corrections below included in the follow-up. It is not an external review.
The issue body, its audited frontier and the open receiving issues #779 and
#740 were read again. The README acceptance matrix records the delivered
artifacts; this report records the weaknesses found and the closure limits.

## Findings and corrections

1. **The production workflow was not exercised in CI.** CI compiled the
   independent verifier and replayed the committed tapes, but would not catch
   a broken producer or incompatibility between `run.py`, `account.py` and
   `pack.py`. `check_reproduction.py` now executes all four q=3 variants from
   fresh inputs, performs a second byte-identical execution for accounting,
   packs the tapes, runs both independent verifiers, and compares the entire
   decoded stream with the committed control hash. The Linux controls job
   runs this workflow with a freshly compiled producer and verifier. All four
   variants passed locally; no receipt or event archive needed changing.

2. **The M1 guard could accept a declaration with its law removed.** It checked
   the retained status, owners and exit tag but not the presence of the actual
   read/feedback/control declaration or outstanding derivation. The guard now
   requires all five substantive M1 fields to be nonempty. Negative controls
   remove the read law and blank the derivation obligation. Both are rejected.
   This is a structural guard; it cannot decide whether arbitrary replacement
   prose is scientifically adequate.

3. **The new Lean history result did not itself discharge its transport
   hypothesis.** The analytic argument and inherited real-valued hop theorem
   supplied that bridge, so this was a formal-coverage limitation rather than
   a counterexample. `six_event_hop_exact` now evaluates the concrete word on
   integer half-units for every integer payload and arbitrary incoming scratch
   values. It proves source preservation, capture and cleanup.
   `six_event_hop_cost` proves the 6/7/7 event/read/write count.
   `six_event_compiledHistory_eq` uses that word to instantiate the complete
   history theorem, without leaving an arbitrary exact-transport hypothesis.

4. **The implementation scope needed sharper wording.** Arbitrary integer
   interventions belong to the unbounded mathematical theorem. The native
   binary and retained receipts use bounded signed 64-bit half-units; arbitrary
   inputs are not a supported machine-level contract. Native peak RSS measures
   the producer, not the complete Python preparation/compression pipeline, and
   a 4 MiB raw codec block is not a total codec-memory bound. The README now
   states these limits explicitly. The semantic operation/storage census is
   unchanged.

## Contract review

**Objective and exit.** The stronger objective of deriving the read/feedback
law from A1--A3 is not achieved. Neither is a universal routing impossibility
proved. Closure uses the issue's explicit declared-M1 alternative, supported
by the conditional construction, complete executions and explicit transfer.
The prior audited frontier's instruction to keep the issue open was not
satisfied merely by changing a status tag: it named missing full-family
executions, refinement and accounting, which are now separately delivered.

**Theorem and interventions.** The native schedule implements the original
`build_causal_poset.py` payload recurrence: initial value `site+1`, followed by
`1+sum(previous-layer metric neighbours)`. The exact menu includes self-reads.
The finite six-event word works independently of payload and scratch values;
tree induction then supplies every receiver, and layer induction preserves
arbitrary initial and local-commit interventions in exact arithmetic. The
full production receipts execute baseline and centre-source interventions;
later-commit and scratch interventions are additionally executed at q=3.
The universal theorem is not inferred from those sampled interventions.

**Order and versions.** The order theorem requires every consumed-writer edge
to respect its logical label and every logical read to have a writer path.
In `verify.cpp`, initial/commit versions acquire their logical label;
export/mean/capture propagate the actual input label; each add checks the
preceding-layer source version and accumulator label. Every consumed writer
must be the actual current writer of the named register. The fixed schedule
checks every add before its target commit. Protected versions cannot be
overwritten. Thus the checked local certificates give both order inclusions.
The separate Python oracle constructs ancestry directly and checks all 6,561
logical pairs per control, without producer labels. The native implementation
and its connection to the abstract certificates remain reviewed executable
code, not a kernel-verified C++ interpreter.

Constant resets consume their own protected zero and do not consume the old
scratch payload. Their events, writes, resource hazards and serial positions
remain retained. The proved projection is semantic value-consumption order;
it does not remove physical interactions or identify a physical causal order.
The older 4,928-event transport test covers old-version rereads after newer
commits; no such experiment is newly asserted for q=21.

**Full receipts and custody.** Both production populations have complete
baseline and source-intervention histories: 19,113,548 events per q=13 run and
292,722,053 per q=21 run, totaling 623,671,202 events. Independent verification
reconstructs all metric decisions, support incidence, host assignment, BFS
paths, primitive arithmetic, writers and completed values. Differential
compression stores literal field differences and restores every row, rather
than generating absent events from a route recipe. The original family
receipts are unchanged. The pinned public simulator supplies the captured
wiring; the native program explicitly simulates the additional M1 law.

**Resources and boundary.** Preparation, every six-event hop, resets, local
accumulations and commits enter the census. All immutable relay versions,
old layer versions, zeros, addresses, scratch ports and accumulators enter
storage. Scalar reads/writes, archive metadata, requested path depths and
BFS/pruning controller work are separately accounted. Shared tree prefixes
are charged once, so the per-path cost is not multiplied into a false
unicast total. Serial unit-event duration counts every primitive. Native
allocation/RSS observations are host measurements; no physical clock, patch
capacity, header channel, quantum instrument, cross-regulator refinement or
volume measure for the added operations is inferred. Metric-menu preparation
and the selected control law remain supplied setup, recorded under M1.

**Explicit transfer.** `specification.json` assigns derivation of the whole
record/read/memory/feedback/placement/control law to `[740, 779]` under PR-52.
The assumption dictionary and PR-52 repeat this scientific obligation. Both
receiving issues were confirmed open. PR-52's `consuming_lanes` is the fixed
V3 umbrella-lane inventory (including #740), not a list of every V4 child
issue; the explicit V4 owners are in the specification and register prose.
PR-52 keeps its `remove` disposition, so closure does not mark the source
derivation completed or silently turn it into a derived fact.

## Validation evidence and remaining limits

- Audit regression: **138 passed** — 96 routing/inherited-transport controls
  and 42 premise/register/surface checks. An initial attempt hit a permission
  error in the existing Windows pytest temp directory; the successful run
  used a fresh workspace-local base directory. No tests were skipped to pass.
- Fresh production workflow: all four q=3 variants reproduced byte-identical
  decoded histories, passed resource accounting and both independent oracles.
- Lean 4.29.1 kernel check with both implicit-variable options disabled:
  concrete word, costs, instantiated history and abstract order passed. The
  order theorem has no axioms; the other printed dependencies are only
  `propext` and, for history equality, `Quot.sound`.
- Hosted checks on the reviewed pre-audit commit `5556050`: full q=13 and q=21
  replays, both platform control jobs, the whole Lean build, both Lean
  certificate jobs, all twelve mandatory shards and both collection gates,
  claim validation and paper preview passed. This is
  evidence for that commit, not a claim that the follow-up's CI had already
  run. Final statuses are attached to the PR's current commit.
- The earlier full local standard mandatory run covered all 143 steps; the
  audit follow-up changes only the proof, guard, controls, workflow and scope
  documentation. The native producer, native replay semantics, inputs and all
  retained tapes are byte-unchanged. The broader optional heavy suite was not
  rerun for this audit.

The audit found no incorrect retained value, induced-order certificate or
semantic resource census. With the corrections above, the package supports
closing #777 **under its declared-M1 exit**. It does not support marking M1
derived, closing #779/#740, or claiming the stronger canonical-only objective.
The PR remains subject to its normal current-commit checks and review.
