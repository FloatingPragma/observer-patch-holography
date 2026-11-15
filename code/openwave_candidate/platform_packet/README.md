# OpenWave submission handoff

Three issues addressed here: a complete 31-row literal-gate inventory tied to one
candidate, a runnable platform-layout native-charge subgate, and verified upstream
submission/access requirements. Physics completion is not asserted.

## Actual route and verified access

Upstream is `openwave-labs/openwave`, main at
`5d686cb049b82cc55f8822210a254d5833be23a7` when inspected 2026-11-14.
The authenticated GitHub account can read but cannot push upstream (`pull=true`,
`push=false`). No open M10 pull request was returned by the GitHub query.
The current roadmap has no active/backlog task and reserves M10.1 as completed.
M10.2 here is a proposal, not an allocated or filed task.

The fastest route is author-run headless computation plus one concern per fork
branch and signed-off pull request. Research tasks are roadmap rows, not platform
issues. Read any newly opened PR conversation before revising it. The platform
runs no general remote physics service for this work; the author supplies compute.
The roadmap says the run record stays at a pinned author commit, while the PR
brings a method note and small headline reproducer. Terms precede any cell-claiming
run. The canonical theory scores one candidate at a time. Declared-source results
must not be presented as carrier-derived physical validation.

`research/` mirrors the real destination
`openwave/xperiments/m10_oph/research/`: it contains a proposed task and standalone
standard-library reproducer. Copy only the concerned task artifacts into a fork,
add its roadmap row, and supply the pinned author commit and method note before
submission. This is a reviewable local handoff; no fork, PR, message or upstream
write was made. The handoff is ready for a task PR once its proposed task ID and pinned author
source are settled. No PR was sent. Upstream merge needs maintainer review;
GitHub authentication itself is available.

## Preparation checks

Run from the bundle root (`code/openwave_candidate` in RER or `openwave-tuning` in the development workspace).

```sh
python platform_packet/validate.py
python platform_packet/research/scripts/m10_2_native_charge.py --run
```

The packet validator verifies 31 distinct literal gates and a shared candidate ID;
it reports zero physical gates ready. The native-charge subgate enumerates 8450
signed local events with exact arithmetic and rejects a corrupted-record control.
It exposes a supplied integer register; it does not prove electron topology.
These already-run controls are not held-out physical tests.

`packet.json` carries every literal criterion, supplied inputs and per-row missing
prerequisites. No invented numeric tolerances are substituted for missing physical
observables. Those entries are deliberately blocked until actual instruments,
state identifications, error budgets and fresh acceptance configurations exist.

## Sources checked

Read directly through GitHub's API: upstream repository permissions,
`openwave/xperiments/m10_oph/research/m10_roadmap.md` and
`tasks/m10_1_task_details.md`. Submission details are cross-checked against
`../../openwave-push/platform/PRECEDENTS.md`, especially sections 5f and 5g.
Literal gates come from `../criteria_snapshot.json`, pinned to the same commit.
No CI run was requested or awaited.

## Portable integrated implementation export

```sh
python platform_packet/export.py /path/to/new-bundle
```

The destination must be new. The exporter preserves source layout and writes a
SHA-256 manifest plus exact replay commands. It bundles the actual candidate,
joined integer-link implementation, ring geometry, 31-row assessment runner and
frozen packet. It excludes prior results, dirty assessments, RER and worker code.
The separate native arithmetic script supplies an independent local control.
All 31 rows name their specific missing physical observable as well as execution
prerequisites; blocked rows cannot become accepted by passing the packet checker.

The preparation validator rejects altered row IDs (including booleans), changed
literal tests, provenance/candidate mismatches and status promotions, including
under `python -O`. `criteria_frozen.json` pins the exact literals, upstream commit
and MODELS.md digest; this checks recorded provenance, not an online refetch.
Run `python platform_packet/test_validate.py` for mutation checks.
