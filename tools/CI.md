# Observer CI

`observer-affected.yml` is the operational replacement for `observer-model.yml`.
The old workflow is a pinned input to a finite custody contract, so its bytes are
retained. The replacement must pass its infrastructure validation before the
old workflow is explicitly disabled in GitHub Actions. Adding the new file alone
does not disable the old workflow or prove a faster end-to-end run.

The replacement preserves the old path filters and adds its own workflow,
infrastructure tools, tests and this document. Superseded routine branch/PR
runs may be cancelled. Scheduled, manual and release runs have unique concurrency
groups and are never cancelled by that policy.

## Routine changes

Three independent jobs run in parallel:

* Infrastructure mutation tests check dependency selection, custody preflight
  and the guarded metadata runner.
* A targeted preflight checks the complete live dependency closure of
  `model/particle_identification/completion_results.json`, followed by the
  unchanged `python -m model.verify_current --hashes-only` global guard.
  The preflight does not claim support for every legacy metadata schema.
* The planner resolves old and current Python imports and declared receipt input
  edges, including child JSON dependencies. It computes one immutable selection
  plan and publishes it as a same-run artifact.

Regression workers start only after all three jobs succeed. They download the
plan by artifact ID and verify its hash, Git head/base and selector source hash.
The dependency graph is not rebuilt on every worker. Whole test modules are
balanced using the existing cost estimates and assigned exactly once.
Each affected worker collects only its assigned modules. It retains the
unchanged `model.pytest_shards` plugin as a local `0/1` partition, including its
`tools.completion_contract_inputs` documentary fixture hooks and original test
comparisons.

Selection is deliberately conservative. Directory proximity is not ownership:
a test may be omitted only when its nearest enclosing package receipt explicitly
lists it in a source inventory. Unowned tests, unresolved metadata and dynamic
loaders stay selected. Unknown changed inputs, unavailable Git bases, module
additions/deletions, input-inventory changes and shared configuration changes
fall back to the original full command. This is selection under the declared
source contracts, not an inference of every possible effect of arbitrary Python.

Three source-byte-pinned loader contracts cover the fixed P-response paper
loader and the carrier/federation archive loaders. Their literal input pins and
the complete 259-file archive name and byte inventory are checked. Archive
imports back into current model/tools code are rejected. A changed loader,
archive root, filename, source byte or input pin forces full fallback. Ordinary
record serialization is not treated as code loading; callable loader indirection
remains uncertain. These are finite reviewed boundaries, not a general loader
exemption.

The plan reports selected modules, reasons, partitions, largest recorded module
costs and estimated shard load. Those estimates balance workers; they are not
measured current runtimes. The default ceiling remains sixteen workers because
more queued jobs do not necessarily reduce elapsed time. Routine sizing targets
a fixed ten minutes of retained work per worker; it does not shrink parallelism
in proportion to selected work and erase the expected latency improvement.

## Full coverage

Scheduled, manual and release-tag runs retain the original six scientific replay
shards and sixteen full model-regression shards. The full regression command,
extra independent receipt commands, immutable custody checks and numerical
comparison rules are unchanged. Scope planning skips dependency discovery when
full coverage is explicitly requested.

The finite input proofs, scientific producers and retained result comparisons
are not weakened by CI selection. Hash-only custody is labelled as custody; it
is not reported as a numerical cold replay.

## Metadata reuse

`immutable_custody_replay.py` accelerates the particle-interface owner by reusing
only the discovered path keys from four reviewed metadata functions. Every hit
rehashes the real files into a fresh map. First and final original discovery,
source bytes and filenames, current owner imports and the expected receipt are
checked. Every scientific calculation and original strict comparison still runs.
Writer and cold replay use separate processes and separate empty caches. Use a
private immutable checkout and an external audit report.

`immutable_custody_pytest.py` is a separate, narrowly pinned route for the twelve
unchanged particle-interface tests. It requires the exact module/configuration,
all twelve collected IDs, all thirty-six passed test phases, unchanged inputs
and the expected receipt. It is explicitly an adapted pytest execution. It may
replace that module in routine CI only after a genuine Linux execution is
accepted; its unit tests alone are not sufficient. Full scheduled coverage can
retain the original plain pytest route.

The routine fast route uses `immutable_walsh_replay.py` with the exact
integer/Fraction Walsh butterfly from `exact_walsh.py`. This is explicitly an
arithmetic substitution as well as metadata reuse. The Sylvester identity
`H_N² = N I` preserves the original coefficients and inverse checks; every
native pulse, including zero coefficients, and the native resource counts
remain. The plan moves the selected twelve-test module to a separate job only
when all reviewed helper, test and configuration hashes match. Ordinary
partitions plus that exact module must cover the selected set without overlap.
A changed helper leaves the original module in ordinary pytest. Full scheduled,
manual and release runs retain the plain original route. Activation requires
actual Linux success of all twelve unchanged tests against a genuine receipt.

## Validation and operational switch

1. Run the infrastructure mutation tests and authenticate the workflow coverage
   comparison against the unchanged original workflow.
2. In an immutable Linux Actions snapshot, run the plan, targeted preflight and
   original global source guard. Verify any proposed exact-test fast route with
   its actual scientific tests before routing a module away from ordinary pytest.
3. After the immutable Linux validation and review pass, disable the superseded
   observer workflow immediately before publishing the reviewed replacement.
   Confirm the main push succeeds and the replacement creates its expected runs.
   Re-enable the old workflow if publication or activation fails. This explicit
   transaction avoids duplicate full suites on the same final head; the old
   workflow is never disabled merely because a new YAML file exists.
4. Compare actual job times and selected costs. Do not equate removing repeated
   collection with removing all scientific test time: a measured existing shard
   spent about three minutes collecting and another thirty-five minutes testing.

GitHub documents [job prerequisites](https://docs.github.com/en/actions/how-tos/write-workflows/choose-what-workflows-do/use-jobs)
and [workflow events](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows).
The workflow uses commit-pinned checkout, setup and artifact actions.

`completion-program.yml` also partitions its existing 869-test program-contract
suite into four complete-module groups. The fourteen other existing commands
run exactly once on shard zero, and all four jobs are required for success.
The original suite was measured at about 1,396 seconds; a fourfold speedup has
not been measured or promised. Pip caching stores dependency downloads only.
