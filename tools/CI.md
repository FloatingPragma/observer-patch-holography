# Continuous integration

The repository runs three workflows.

* `.github/workflows/ci.yml` runs on pushes to `main`, pull requests and release tags. The required `CI result` collects every job below.
  * `ledgers` checks the claim and observable registers, the QFT/WZ theorem oracle, the Lean trust inventory and CI configuration.
  * `model` runs the complete model test suite and the executable replays in parallel. Both use the same pinned environment. The replays cover the carrier, candidate laws, evaluations, experiments and model certificate lanes.
  * `einstein` runs six suites of matched finite Einstein tests with the model dependencies.
  * `certificates` distributes the mandatory checks and root certificate lanes across six shards. Shard 0 also runs the additional scientific evidence tests.
  * `papers` rebuilds every review PDF and the book and checks the committed artifacts.
* `.github/workflows/lean.yml` runs the change-aware Lean build when `Lean/**` or a Lean-facing tool changes.
  `tools/check_lean_ci_budget.py` keeps every job under 30 minutes.
* `.github/workflows/nightly.yml` runs the full Lean build on a schedule and on dispatch. A failed scheduled run
  opens or updates one issue labelled `nightly-failure`. Dispatch it on a branch with
  `gh workflow run nightly.yml --ref <branch>`. Dispatches also compare the latest public release assets with their tagged commit.

The exhaustive certificate suites run on demand: `python tools/run_mandatory_suite.py --certificates-only`.

The ledger job logs each register and candidate grade change with its ID and old and new status. It uses the pull request's base SHA or the previous commit on a push to `main`. After regenerating candidate assessments with `python -m model evaluate --write`, run `python tools/report_grade_changes.py --base HEAD` to compare your working tree with the local commit. A different commit, tag or branch can be supplied. Register rows, executable checks, full criterion grades, candidate premises and joint requirements are tracked separately. Added and removed grades appear with `[absent]` on the corresponding side. Editing evidence or wording alone leaves the grade report empty. The existing target validator independently checks upgrades against the earlier statement.

## Adding a check

Add a lane to `tools/ci_lanes.py` with its commands, environment, install profile and measured minutes.
`tools/test_ci_lanes.py` checks that every lane runs exactly once and that the paths it names exist. A check that
needs more than a few minutes does not belong in CI.

## Custody pins

Receipts bind the code that produces or verifies a result and the data it reads. `code/pin_scope.py` (loaded
from tools/ and model/ through `tools/pin_scope.py`) defines what a pin may bind: prose, papers, documentation,
READMEs, workflow files, Git attribute files and test modules are never pinned. Verifiers that use it drop such
entries from recorded pin maps before comparing, so editing prose does not require renewing their receipts.
