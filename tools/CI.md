# Continuous integration

Three workflows run the repository checks.

* `.github/workflows/ci.yml` runs on every push to `main`, every pull request
  and every release tag. It is the one required check (`CI result`). Every job
  has a hard timeout, so a run finishes in under 45 minutes.
  * `gates`: claim registry, generated registers and ledgers, the QFT/WZ
    theorem oracle, the Lean `native_decide` inventory and the Lean CI budget
    wiring.
  * `suite`: the standard steps of `tools/run_mandatory_suite.py` in ten
    shards (`--skip-long`). Shard 0 also runs the additional evidence tests.
  * `suite-long`: the long standard steps (`--long-only`, the three-axiom
    certificate campaign).
  * `lanes`: the replay lanes of `tools/ci_lanes.py` (tier `push`), packed by
    their measured minutes into fixed buckets per install profile.
  * `paper-preview`: rebuilds every review PDF and the book and requires the
    committed PDFs to match.
* `.github/workflows/lean.yml` runs the change-aware Lean build when `Lean/**`
  or a Lean-facing tool changes. `tools/check_lean_ci_budget.py` keeps every
  job under 30 minutes.
* `.github/workflows/nightly.yml` runs the heavy replays on a schedule and on
  dispatch: the full mandatory suite, exhaustive certificates, observer
  construction and regression, the completion contract shards, OpenWave
  packages, external-source replays and the full Lean build. It never blocks a
  push. A failed scheduled run opens or updates one issue labelled
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
are never pinned. Builders leave such paths out of the maps they write, and
verifiers drop them from recorded maps before comparing, so editing prose or
a test never requires renewing a receipt.

Only receipts carry pins. Review, campaign, class and work-queue records name
the files and receipts they rely on, and their checkers keep their structural
and scientific assertions, but they do not compare recorded digests: each
receipt is checked by its own source map and replay, and git history records
earlier bytes. Package manifests under `model/` (`model/verify.py` PACKETS)
bind data artifacts only, not code or JSON receipts.

## Regenerating receipts

After a change to code that receipts pin, `python tools/replay_model_receipts.py
--output <new dir>` lists the stale receipts and their dependency order
without writing. To regenerate them on Linux, dispatch the nightly workflow in
regenerate mode:

    gh workflow run nightly.yml --ref <branch> -f regenerate=true

It runs only `tools/replay_model_receipts.py --execute --discover --keep-going`,
which replays every stale registered or discovered receipt with its own
producer in dependency order, and uploads the changed files and the report as
the `regenerated-receipts` artifact. A producer whose payload changes beyond
its pins stops for review and keeps its original bytes.

Two old workflow files stay byte-identical because the receipts of their
lanes pin them, and several openwave package manifests and model receipts pin
those receipts in turn: `m1-source-realization.yml` and
`source-gravity-admissibility.yml`. `observer-model.yml` and
`completion-program.yml` are read only by older workflow-structure tests.
Their commands also run in ci.yml or nightly.yml; they are disabled in
Actions.
