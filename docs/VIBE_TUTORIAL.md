# Vibe tutorial

**Let's vibe our way to the theory of everything.**

Observer Patch Holography (OPH) is an open collaboration to find the computation
that outputs our Universe. Every observable it has to reproduce, and every
premise its theorems wait on, has a row and a status on the
[status page](https://status.floatingpragma.ai/), which also shows how much of
the whole is complete. The work is short of inference compute. If you have
quota on one of the research-grade models below, you can move a row from your
sofa.

## What you need

- A GitHub account, git and the [GitHub CLI](https://cli.github.com/).
- Python 3.12.
- A coding agent running a research-grade model. Only these models are
  accepted at the moment:
  - the newest Claude Fable, in [Claude Code](https://github.com/anthropics/claude-code);
  - GPT-6.2 Astra, in [Codex](https://github.com/openai/codex);
  - Gemini Carbon, in [Gemini CLI](https://github.com/google-gemini/gemini-cli).
- A few hundred megabytes of disk for a shallow clone.
- For proof rows only: [elan](https://github.com/leanprover/elan), the Lean
  toolchain manager.

## Six steps

1. Fork and clone, shallow:

   ```sh
   gh repo fork FloatingPragma/observer-patch-holography --clone -- --depth 1
   cd observer-patch-holography
   ```

2. Install the dependencies and check that everything is current:

   ```sh
   python3.12 -m venv .venv && source .venv/bin/activate
   pip install -r requirements.txt
   python -m model evaluate --check
   python tools/check_observables.py
   ```

   The last command prints the same completion figure as the status page.
   Proof rows also need the Lean library: `cd Lean && lake exe cache get && lake build`.

3. Start your agent in the clone (`claude`, `codex` or `gemini`) and select
   one of the accepted models.

4. Paste the prompt below and let it work.

5. Read what it did. Ask it which run outcome would have proved it wrong. If it
   has no answer, it has no test yet.

6. Ask it to open the pull request against
   `FloatingPragma/observer-patch-holography` and to fill in the template.

## The prompt

```text
Please continue to complete Observer Patch Holography (OPH) and the OpenWave
challenges. First read AGENTS.md, docs/OPENWAVE_SUBMISSION_GUIDELINES.md and
docs/registers/README.md in full; they bind this work.

Goal: raise the status of rows in the OPH register
(docs/registers/OBSERVABLES.json, mirrored at https://status.floatingpragma.ai/)
and in the model scorecard (model/evaluation/SCORECARD.md). Aim to close one
OpenWave criterion together with OPH premises that many observables wait on.
Work on the main candidate law, Interval Echo, on the fixed carrier. Use as
many parallel workers as the work allows, each in its own git worktree, and
move OPH along as efficiently as you can.

Rules:
- Never put in by hand what the theorems derive or what observers measure
  (masses, the number of dimensions, a lattice, a metric, a gauge group).
  Every calibration needs declared provenance.
- Statuses change only through evidence. Do not edit criteria, the scorecard
  generator or the register to change a status.
- For each claim, fix before the run the outcome that would falsify it.
- Run python -m model evaluate --check, python tools/check_observables.py and
  the tests of the files you touch. Leave the full suite to CI.
- Report what failed as well as what passed.

Open one pull request per row against FloatingPragma/observer-patch-holography
and fill in its template: the row ID, the evidence, the exact commands that
replay it, the outcome that would have falsified it, and what failed.
```

## What gets merged

A reviewer merges a pull request after replaying it and checking that the
evidence meets the row's criterion. The template asks for the row ID, the
evidence and the files that hold it, the replay commands, the falsifying
outcome fixed before the run, and what failed along the way.

Four things block a merge: work by a model outside the accepted list, an
outcome put in by hand, a status changed by editing criteria, generators or
the register, and a claim nobody can replay.
A failure shown on a candidate's own observer records is a contribution too.
It saves everyone the same detour.

## Credit

Contributors whose merged pull request moves a row on the status page are
listed as special contributors in the repository [README](../README.md) and in
the [synthesis paper](https://philpapers.org/rec/MUEFOC). Everyone with merged
commits appears in the author list of the next
[Zenodo release](https://doi.org/10.5281/zenodo.24191882) of the paper stack.

## Tips

- Premises pay twice. A Conditional observable scores half until the premises
  it waits on are shown, so one demonstrated premise raises every observable
  that depends on it. The status page lists both kinds.
- Give each parallel worker its own `git worktree` so they do not overwrite
  each other.
- The rules your agent follows are in [AGENTS.md](../AGENTS.md), the
  [OpenWave submission guidelines](OPENWAVE_SUBMISSION_GUIDELINES.md) and the
  [register format](registers/README.md).
- Agents like to declare victory early. The falsifier question is the quickest
  way to find out whether they have earned it.
