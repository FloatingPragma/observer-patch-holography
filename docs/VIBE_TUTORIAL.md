# Vibe physics tutorial

**Let's vibe our way to the theory of everything.**

Observer Patch Holography (OPH) is an open collaboration to find the computation
that outputs our Universe. Every observable it has to reproduce, and every
premise its theorems wait on, has a row and a status on the
[status page](https://status.floatingpragma.ai/), which also shows how much of
the whole is complete. The work is short of inference compute. If you have
quota on one of the research-grade models below, you can move a row from your
sofa.

This tutorial assumes no experience with GitHub or coding agents. It takes you
from installing the tools to running one of the prepared prompts, checking
what the agent did and opening a pull request.

## How it works

You give a coding agent a prepared prompt. The agent reads the repository,
does the research in code, Lean proofs or paper text, runs the checks and
writes up what it found. You read the result, ask it which outcome would have
proved it wrong, and send the work to the maintainers as a pull request. A
reviewer replays it and merges it when the evidence holds.

Three kinds of work are open, each with its own prompt:

- **Rows.** Raise the status of an observable or a premise on the status page.
- **Audit.** Find and repair defects in existing evidence, under
  [issue #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
- **Obstruction search.** Look for results that would justify rejecting OPH,
  under [issue #1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078).

## Words you will meet

- **Repository (repo):** the project's files together with their full history,
  hosted on GitHub. This one is
  [FloatingPragma/observer-patch-holography](https://github.com/FloatingPragma/observer-patch-holography).
- **Fork:** your own copy of the repository on GitHub. You can change it freely
  without touching the original.
- **Clone:** a copy of a repository on your computer.
- **Branch:** a named line of changes inside a repository. Each piece of work
  gets its own branch, so it stays apart from `main`, the branch everyone
  shares.
- **Commit:** a saved set of changes with a short message describing them.
- **Push:** uploading your commits to your fork on GitHub.
- **Pull request (PR):** a request to the maintainers to pull the commits of
  your branch into the main repository. GitHub shows the changed lines, runs
  the automated checks and holds the review conversation on the pull request.
- **CI (continuous integration):** the automated checks GitHub runs on every
  pull request. A green tick means they passed, a red cross means one failed.
- **Merge:** the maintainers accept the pull request, and its commits join
  `main`.
- **Issue:** a numbered discussion thread on GitHub for a task, a bug or a
  research question.
- **Coding agent:** a program in your terminal that lets an AI model read
  files, run commands and edit code in a folder you choose. It asks for your
  approval as it goes.
- **Row:** one line of the status page, an observable or a premise, with one of
  four statuses: Open, Partial, Conditional or Demonstrated.
  [docs/registers/OBSERVABLES.json](registers/OBSERVABLES.json) defines them.

## What you need

- A computer running macOS, Linux, or Windows with WSL (step 1 covers it).
- A free [GitHub account](https://github.com/signup).
- Quota on one research-grade model. Only these models are accepted at the
  moment:
  - the newest Claude Fable, in Claude Code;
  - GPT-6.2 Astra, in Codex;
  - Gemini Carbon, in Gemini CLI.
- About 1 GB of free disk for the clone and its Python packages.
  Proof rows need several gigabytes more for the Lean library.

## Step 1: Install the basic tools

You need git (version control), the GitHub CLI `gh` (GitHub from the
terminal), Python 3.12 (the model and its checks) and Node.js 22 or newer (the
installer of the coding agents). Open a terminal and follow the part for your
system.

**macOS.** Install [Homebrew](https://brew.sh), then:

```sh
brew install git gh python@3.12 node
```

**Linux (Ubuntu 24.04 or newer, Debian).**

```sh
sudo apt update
sudo apt install git gh python3.12 python3.12-venv
```

The Node.js in the system package list is often too old. Install version 22 or
newer from [nodejs.org](https://nodejs.org/en/download).

**Windows.** Install the
[Windows Subsystem for Linux](https://learn.microsoft.com/windows/wsl/install):
open PowerShell as administrator, run `wsl --install`, restart, and open the
Ubuntu app. Follow the Linux part inside Ubuntu. Every command in this tutorial
runs in that Ubuntu terminal.

Check that everything is in place; each command prints a version number:

```sh
git --version
gh --version
python3.12 --version
node --version
```

Then sign in to GitHub from the terminal:

```sh
gh auth login
gh auth setup-git
```

Choose GitHub.com, HTTPS and "Login with a web browser". The second command
lets git push with the same login.

## Step 2: Install a coding agent

Install the agent that runs your model:

| Model | Agent | Install | Start |
| --- | --- | --- | --- |
| the newest Claude Fable | [Claude Code](https://github.com/anthropics/claude-code) | `npm install -g @anthropic-ai/claude-code` | `claude` |
| GPT-6.2 Astra | [Codex](https://github.com/openai/codex) | `npm install -g @openai/codex` | `codex` |
| Gemini Carbon | [Gemini CLI](https://github.com/google-gemini/gemini-cli) | `npm install -g @google/gemini-cli` | `gemini` |

The linked pages list other installers and the sign-in options. On first start
the agent asks you to sign in; use the account that holds your quota. Inside
the agent, the `/model` command lists the models your account can use: select
the accepted one. You can also pass it at start, for example
`claude --model <name>`.

The agent asks before it edits a file or runs a command. Approve edits and
commands inside the clone. Decline anything that reaches outside it, such as
your home folder, your SSH keys or other repositories.

## Step 3: Fork and clone the repository

```sh
gh repo fork FloatingPragma/observer-patch-holography --clone -- --depth 1
cd observer-patch-holography
```

The first command creates your fork on GitHub and a clone of it on your
computer. `--depth 1` downloads the newest commit without the full history,
which saves several gigabytes. In the clone, the remote `origin` points at your
fork and `upstream` at the main repository; `git remote -v` shows both.

## Step 4: Install the Python packages and run the checks

```sh
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m model evaluate --check
python tools/check_observables.py
```

The first two lines create and enter a private Python environment for the
project; run `source .venv/bin/activate` again in every new terminal.
`python -m model evaluate --check` confirms that the
[scorecard](../model/evaluation/SCORECARD.md) matches the executable model.
`python tools/check_observables.py` validates the register and prints the same
completion figure as the status page. If either check fails on a fresh clone,
open an issue with the output. That is a contribution too.

**Proof rows only.** Rows that need Lean proofs need the Lean toolchain. Install
its manager [elan](https://github.com/leanprover/elan) and build the library:

```sh
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
source ~/.elan/env
cd Lean
lake exe cache get
lake build
cd ..
```

`lake exe cache get` downloads the prebuilt Mathlib library, which saves hours
of compiling.

## Step 5: Choose your work and run its prompt

Start each piece of work on its own branch, named after the row or the defect:

```sh
git switch -c row-born-rule
```

Then start your agent in the clone (`claude`, `codex` or `gemini`), select the
accepted model and paste one of the three prompts below. The agent can work for
an hour or longer and spend much of your quota. You can watch, answer its
questions, and interrupt it with Esc or Ctrl+C. A prompt that opens several
pull requests creates one branch per pull request itself.

### Rows

The status page lists every row with its status, its premises and the files
behind it. The [scorecard](../model/evaluation/SCORECARD.md) lists the
criteria of the OpenWave comparison for the executable model. Premises pay
twice: a Conditional observable scores half until the premises it waits on are
shown, so one demonstrated premise raises every observable that depends on it.

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

### Audit: issue #1033

[Issue #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
is the standing audit of existing evidence. Its scope:

- numerical and logical defects in code: wrong results, silent fallbacks or
  floors, invalid inputs accepted, missing validation;
- statements in papers, Lean or docs that are wrong or claim more than their
  proof supports;
- receipts, registries and custody chains that fail to reproduce;
- public overclaims, corrected to the scientific status of the evidence;
- simplifications that merge duplicated machinery, with any change in results
  documented.

Each audit pull request repairs one defect class, reproduces the defect before
the fix, confirms the fix with an independent control (exact arithmetic, a
separate implementation, a countermodel, or a replay without the original
producer), states the impact on other code, papers and public text, and
references #1033. Frozen registrations and pinned evidence stay unchanged
unless the pull request states and justifies the change. A finding you do not
fix goes into a comment on the issue, with a reproduction. A corrected theorem
or receipt improves the evidence; it supplies no missing physical premise.

```text
Please audit existing Observer Patch Holography (OPH) evidence for defects,
under issue #1033 of FloatingPragma/observer-patch-holography. First read
AGENTS.md, docs/OPENWAVE_SUBMISSION_GUIDELINES.md and docs/registers/README.md
in full, then read the issue with
gh issue view 1033 --repo FloatingPragma/observer-patch-holography --comments

Goal: find wrong results, silent fallbacks, invalid inputs that are accepted,
proofs or paper statements that claim more than they show, and receipts that
fail to reproduce. Pick one defect class and repair it.

Rules:
- Reproduce the defect before you fix it, and keep the reproduction as a test.
- Confirm the fix with an independent control: exact arithmetic, a separate
  implementation, a countermodel, or a replay without the original producer.
- State the impact on other code, receipts, papers and public text.
- Leave frozen registrations and pinned evidence unchanged unless the pull
  request states and justifies the change.
- Run python -m model evaluate --check, python tools/check_observables.py and
  the tests of the files you touch. Leave the full suite to CI.
- Report what failed as well as what passed.

Open one pull request per defect class that references #1033, and fill in its
template. For a defect you find and do not fix, draft a comment for #1033 with
a reproduction and show it to me before posting it.
```

### Obstruction search: issue #1078

[Issue #1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
searches for results that would justify rejecting OPH or its claim to explain
physics from one process. OPH's physical correctness is the hypothesis under
test, and a decisive negative result counts as a successful research outcome.
The issue names four lines of attack:

| Line of attack | What a result has to show |
| --- | --- |
| Core inconsistency | A contradiction derived from the fixed basis of the theory. |
| Joint-realization impossibility | No admissible source can supply a necessary combination of geometry, reversible response, clocks, quantum histories or physical-limit premises, with every model class and extra assumption stated. |
| Failure of explanatory unification | Parts the theory claims to derive (the kinetic law, clocks, matter, couplings) can vary independently while every premise holds, shown by fully checked countermodels. |
| Empirical exclusion | An unavoidable physical consequence that every remaining admissible realization gets wrong against an eligible observation. |

A failed candidate, or a no-go for a restricted class, rejects that candidate
or class and leaves the programme standing; the result says exactly what
survives. A search that finds nothing is inconclusive and supplies no
confirmation of OPH.

```text
Please work on the programme-wide obstruction search in issue #1078 of
FloatingPragma/observer-patch-holography. First read AGENTS.md,
docs/OPENWAVE_SUBMISSION_GUIDELINES.md, docs/registers/README.md and
docs/BASIS.md in full, then read the issue with
gh issue view 1078 --repo FloatingPragma/observer-patch-holography --comments

Goal: a result that would justify rejecting OPH or a narrower OPH claim. A
decisive negative result counts as success. Choose one of the four lines of
attack in the issue: core inconsistency, joint-realization impossibility,
failure of explanatory unification, or empirical exclusion.

Rules:
- Fix the target before testing: the commit you start from, the exact theorem
  statements and premises, the candidate laws and the calibrations allowed.
- Start with analytic reductions and small exact counterexamples. A finite
  search that finds no witness proves no impossibility.
- A failed candidate or a restricted no-go excludes that candidate or class.
  State exactly what survives.
- Record for each result: the target commit, the premises, the class it
  covers, the evidence, the escape routes that survive and the consequence.
- Run python -m model evaluate --check and the tests or Lean builds of the
  files you touch. Leave the full suite to CI.

Put code and proofs in a pull request that references #1078, and fill in its
template. Draft a summary comment for #1078 and show it to me before posting
it.
```

Further research issues carry the label
[program:physics](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics).

## Step 6: Check the work

Read what the agent did before anything goes to GitHub. Ask it for a summary,
then ask the two questions below. If it has no answer to the first, it has no
test.

```text
Which run outcome, fixed before the run, would have proved this result wrong?
Show me where the run checks for it.
```

```text
List every approach you tried that failed, with the command that shows the
failure.
```

Then look at the changes yourself and repeat the checks:

```sh
git status
git diff --stat
python -m model evaluate --check
python tools/check_observables.py
```

If the agent changed a status by editing the register, the scorecard or a
criterion, ask it to undo that edit: statuses change only through evidence.

## Step 7: Open the pull request

A pull request asks the maintainers to take your branch into the main
repository. It carries the changed files, a description filled in from the
repository's [template](../.github/pull_request_template.md), the CI results
and the review. The simplest route is to let the agent open it:

```text
Commit the work on a branch named after the row or the defect, push it to my
fork (the remote origin) and open a pull request against main of
FloatingPragma/observer-patch-holography with gh pr create. Fill in every
section of .github/pull_request_template.md. Tick a check only if you ran it,
and name the coding agent and model you used.
```

By hand, add the files that `git status` lists as changed, commit, push and
open the form in your browser:

```sh
git add path/to/changed/file another/changed/file
git commit -m "One line saying what the change shows"
git push -u origin HEAD
gh pr create --repo FloatingPragma/observer-patch-holography --base main --web
```

The form opens with the template. Fill in each section: the row or issue, the
evidence, the commands that replay it, the falsifying outcome, what failed and
the checks. If you are unsure whether the work is ready, choose "Create draft
pull request"; a draft can be marked ready later.

## Step 8: After you open it

- **CI.** The checks run within minutes of opening. A red cross links to the
  log of the failed check. Ask the agent to fix it with the prompt below and
  push to the same branch; the pull request picks up new commits by itself.
- **Review.** A reviewer replays the work and comments. Answer each comment,
  with new commits where the comment asks for a change.
- **Outcome.** The pull request is merged, or closed with the reason. A closed
  pull request that records a failure clearly saves the next person the same
  detour.

```text
Read the checks and the review comments on my open pull request with
gh pr checks and gh pr view --comments. Fix what they point out, run the
checks again and push to the same branch. Draft a reply to each review comment
and show it to me before posting it.
```

If your quota runs out midway, ask the agent to commit what it has to the
branch with a message saying where it stopped. Later, start it in the same
clone and paste:

```text
Continue the work on this branch. Read git log and the diff against main
first, then carry on from where the last commit says it stopped.
```

## What gets merged

A reviewer merges a pull request after replaying it and checking that the
evidence meets the row's criterion, or the rules of the issue it answers. The
template asks for the row or issue, the evidence and the files that hold it,
the replay commands, the falsifying outcome fixed before the run, and what
failed along the way.

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

- One pull request per row or defect. Small pull requests are reviewed faster.
- Give each parallel worker its own `git worktree` so they do not overwrite
  each other.
- The rules your agent follows are in [AGENTS.md](../AGENTS.md), the
  [OpenWave submission guidelines](OPENWAVE_SUBMISSION_GUIDELINES.md) and the
  [register format](registers/README.md).
- Agents like to declare victory early. The falsifier question is the quickest
  way to find out whether they have earned it.
- Before starting a new piece of work, bring your clone up to date with the
  main repository: `git switch main && git pull upstream main`.
- Questions about this tutorial or a contribution go to
  [info@floatingpragma.ai](mailto:info@floatingpragma.ai).
