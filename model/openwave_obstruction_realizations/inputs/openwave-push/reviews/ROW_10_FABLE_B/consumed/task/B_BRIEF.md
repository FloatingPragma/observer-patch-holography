# Binding brief for workers of the second Fable session (prefix B), 2026-11-14

Owner, 13:50 and 13:57 (+07): list the open steps to 31 of 31 validated on OpenWave, execute the priority steps
(those that unblock several rows first) with maximum parallelism, coordinate with the other parties, and repeat
until every row is validated.

Parties: Codex (`openwave-tuning/sync/CODEX.md`), an external auditor with three lanes (`sync/AUDITOR.md`,
`sync/CHARGE_SPIN_LANE.md`, `sync/INTERACTING_LANE.md`), the first Fable session (`sync/FABLE.md`; its open steps
S1 to S13 are in entry F037, its state in F039), and this session (`sync/FABLE_B.md`). Index of all entries:
`/Users/muellerberndt/Projects/oph-meta/openwave-tuning/CODEX_FABLE_SYNC.md`. The group counts a row as
internally complete when the complete literal criterion holds on a named source with reproducible evidence and
two producer-independent full-criterion reviews (`openwave-tuning/row_completion/README.md`). Official platform
status is separate: 0 of 31. Cells move only by a filed task and the maintainer's review.

All relative paths below are under `/Users/muellerberndt/Projects/oph-meta/openwave-push/` unless they start with
`openwave-tuning/` or name the theory repository.

## Rules

- Write only in your own directory (named in your task file). Never write under `openwave-tuning/`, under the
  theory repository `/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/` (read only), or in
  another worker's directory. `steps/`, `blockers/`, `packets/`, `kernel/`, `chain/`, `issues/`, `bridgework/`,
  `rows_native/`, `reviews/`, `rowfull/ROW_31/` belong to the first Fable session: read only.
- Reuse first (`REUSE_BRIEF.md`). Import existing code read only with `PYTHONDONTWRITEBYTECODE=1`, or copy it
  into `consumed/` in your directory and record its SHA256. Never run a case that a producer marks reserved,
  held out or unrun.
- Custody: every file of another worker or party that you cite as evidence and that is not frozen under a
  manifest (a verdict, a status file, a README, a ledger entry) is copied into `consumed/` in your directory
  with its SHA256 and the time of the copy. A hash without the bytes is not custody: such files are rewritten
  by their owners.
- Interpreter `/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/.venv/bin/python` with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`.
- The machine is in memory thrash (10 cores, 24 GB, swap nearly full; wall time is 10 to 30 times CPU time).
  One python process at a time per worker, arrays under 200 MB, each script under five minutes of CPU. Exact
  linear algebra, Bloch matrices and transfer matrices come before Monte Carlo; seeded runs stay small. Never
  `pkill` or `killall`; stop only processes you started, by pid. Use `date` for timestamps.
- Every claim gets an executed gate that can fail and at least one executed control that fails it (a mutated
  law or a wrong mechanism). Thresholds are fixed and written with a timestamp before the run they judge.
- What may be declared (the group's rule): carrier, routing, coherent continuation, compiler, preparation and
  schedule are architecture. A number set from the target it is compared with is calibration and is named with
  its ancestry. Strengths (alpha, G, Lambda) are never put in. Say for every number which source it is on;
  results from different sources combine only with a transfer proof.
- The stricter note on OPH's column at the platform: a result that follows from installed standard physics
  does not earn the cell; an observable computed on the carrier's own process with something still supplied is
  a case for partial; a source derived from the carrier in-platform is the bar for validated. State for each
  result what a declaration puts in and what is found.
- Writing: statements of fact in the present tense, no narration of progress, never the word "honest", no
  em-dashes. No status is written as awarded.
- Your final message (under 400 words) goes to the coordinator: the decisive finding first, then the file
  list. If the task cannot be done, deliver the exact remaining obligation and the smallest construction that
  would supply it.
