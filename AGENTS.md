# Repository guidance for agents

Research results belong in the papers, the Lean library and the code. The
[writing style guide](docs/policies/STYLE_GUIDE.md) binds all prose.

## The `docs/` directory

The full rules are in the [docs index](docs/README.md#rules-for-this-directory).

- The `docs/` root holds only the selected high-level programs, references
  and ledgers, such as `COMPLETION_PROGRAM.md`, `OPH_FALSIFICATION_PROGRAM.md`
  and `POSTDICTION_LEDGER.md`. Do not add a root document unless it is a new
  program or ledger of that kind.
- Never add a daily note, progress report, audit result, integration record,
  PR handoff or status summary anywhere under `docs/`, and never put a date in
  a file name or heading there. Record integration and review history in the
  GitHub pull request or issue.
- A note in `docs/research/` is allowed only when it states a current result
  that its code package verifies, without history, dated findings or progress
  narrative. Prefer extending the owning paper or package README, because
  separate notes go stale.
- Release papers do not cite files in `docs/`. Cite papers, Lean modules, code
  packages and public preprints.
- Code, receipts, reviews and registers bind code, Lean and receipts as
  evidence, never a note in `docs/`.
