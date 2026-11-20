# Continuous integration

Five workflows run the repository checks.

* `.github/workflows/ci.yml` runs on every push to `main`, every pull request
  and every release tag. It is the one required check (`CI result`). Every job
  has a hard timeout, so a run finishes in under 45 minutes.
  * `gates`: claim registry, generated registers and ledgers, the QFT/WZ
    theorem oracle, the Lean `native_decide` inventory and the Lean CI budget
    wiring.
  * `official-model`: the A-contrast-memory-spherical-v1 model (`python -m model` and its
    tests under `model/contrast_memory_spherical_v1`).
  * `experimental-scalar`: the scalar channel of the model, its extended
    preparation and its cage (`model/scalar_channel`), with the replay of
    law 1 from the construction manual.
  * `suite`: the standard steps of `tools/run_mandatory_suite.py` in ten
    shards. Shard 0 also runs the additional evidence tests.
  * `lanes`: the replay lanes of `tools/ci_lanes.py` (tier `push`), packed by
    their measured minutes into fixed buckets per install profile.
  * `paper-preview`: rebuilds every review PDF and the book and requires the
    committed PDFs to match.
* `.github/workflows/lean.yml` runs the change-aware Lean build when `Lean/**`
  or a Lean-facing tool changes. `tools/check_lean_ci_budget.py` keeps every
  job under 30 minutes.
* `.github/workflows/nightly.yml` runs the heavy replays on a schedule and on
  dispatch: the full mandatory suite, exhaustive certificates, the nightly
  replay lanes and the full Lean build. It never blocks a push. A failed scheduled run opens or updates one issue labelled
  `nightly-failure`. Dispatch it on a branch with
  `gh workflow run nightly.yml --ref <branch>`.

## Adding a check

Add a lane to `tools/ci_lanes.py` with its commands, environment, install
profile, tier and measured minutes. `tools/test_ci_lanes.py` checks that every
lane runs exactly once and that the paths it names exist. A lane that needs
more than about 15 minutes or an external checkout belongs to the `nightly`
tier or to its own job in `nightly.yml`.

## Custody pins

Receipts bind the code that produces or verifies a result and the data it
reads. `code/pin_scope.py` (loaded from tools/ and model/ through
`tools/pin_scope.py`) defines what a pin may bind: prose, papers,
documentation, READMEs, workflow files, Git attribute files and test modules
are never pinned. Verifiers that use it drop such entries from recorded pin
maps before comparing, so editing prose does not require renewing their
receipts.

Two workflow files stay byte-identical because the receipts of their lanes
pin them: `m1-source-realization.yml` and `source-gravity-admissibility.yml`.
Their commands also run in ci.yml.
