# Continuous integration

Three workflows run the repository checks, and none runs longer than thirty minutes.

* `.github/workflows/ci.yml` runs on every push to `main`, every pull request and every release tag. It is the one
  required check (`CI result`). Its four jobs run in parallel, each with a hard timeout of at most 25 minutes, so a
  run takes about fifteen minutes.
  * `ledgers`: the claim registry, the generated registers and ledgers, the QFT/WZ theorem oracle, the Lean
    `native_decide` inventory, the Lean CI budget wiring and the tests of this CI layout.
  * `model`: the fixed carrier, every candidate law (`python -m model`, `python -m model run --law ...`), their
    evaluation (`python -m model evaluate --check`), the experiments, every test under `model/`, and the theory
    certificates that use the model pins (lanes of profile `model`).
  * `certificates`: six shards, each running its share of the steps of `tools/run_mandatory_suite.py` and of the
    replay lanes of `tools/ci_lanes.py` (profile `root`). Shard 0 also runs the additional evidence tests.
  * `papers`: rebuilds every review PDF and the book and requires the committed files to match.
* `.github/workflows/lean.yml` runs the change-aware Lean build when `Lean/**` or a Lean-facing tool changes.
  `tools/check_lean_ci_budget.py` keeps every job under 30 minutes.
* `.github/workflows/nightly.yml` runs the full Lean build on a schedule and on dispatch. A failed scheduled run
  opens or updates one issue labelled `nightly-failure`. Dispatch it on a branch with
  `gh workflow run nightly.yml --ref <branch>`.

The exhaustive certificate suites run on demand: `python tools/run_mandatory_suite.py --certificates-only`.

## Adding a check

Add a lane to `tools/ci_lanes.py` with its commands, environment, install profile and measured minutes.
`tools/test_ci_lanes.py` checks that every lane runs exactly once and that the paths it names exist. A check that
needs more than a few minutes does not belong in CI.

## Custody pins

Receipts bind the code that produces or verifies a result and the data it reads. `code/pin_scope.py` (loaded
from tools/ and model/ through `tools/pin_scope.py`) defines what a pin may bind: prose, papers, documentation,
READMEs, workflow files, Git attribute files and test modules are never pinned. Verifiers that use it drop such
entries from recorded pin maps before comparing, so editing prose does not require renewing their receipts.
