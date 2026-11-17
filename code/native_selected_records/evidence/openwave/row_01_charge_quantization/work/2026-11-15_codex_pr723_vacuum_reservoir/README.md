# Observed unit selection with a native vacuum reservoir

Official row1 is 🚧 under the 2026-11-15 PR #723 decision. This packet supplies
new local evidence for the review's proposed multi-load unit-selection route.
It is an unpinned N+ candidate; no full-row, G6 or G9 award is claimed.

With a public capacity bound B fixed before preparation, any twelve-port input
in[-B,B]^12 is coupled to B−1 empty native cells. The resulting N=12B source has
|Q|≤N. Native repair therefore settles almost surely to absolute0/+1 or0/−1
loads. This applies to every admitted input, including signed and uniformly
double/triple FIRST-cell preparations. A known INITIAL vacuum datum and complete
native current records provide the observer's absolute measurement; the
reference is allowed to fill during the process.

The frozen B=2,3 execution reports:

- **36/36 native histories select units**, first entry at most4641 attempts;
  fixed burn8192. **0/36 swap/hold controls select units**.
- **72/72 current-only readings** are conclusive and recover the complete
  initial audit vector. The independent reader imports no producer code.
- **207084 exact native branches**,1028 states on the new graphs, closed unit
  classes, full reachability and exact stationary current-rate checks.
- One earlier record calibration gives u=1. All native held-out currents have
  magnitude0/u; every control transports a larger magnitude. **4211 nonzero
  events in73728 attempts**, versus conditional prediction241106/57, pass the
  frozen alpha10^-6 martingale test without a stationarity approximation.
- Five feasible detector states, a retained timeout, and the full-graph
  uniform-two boundary are checked. Eight integrity mutations reject.

Main literal execution counts are591418 native and377856 changed-law events;
diagnostics are separate. Initial execution took15.607seconds and independent
verification7.150seconds. Exact timing for the composite replay and mutation
suite is in FINAL_RUNS.json.

## Reproduce and inspect the whole primary path

```sh
PYTHONDONTWRITEBYTECODE=1 python3 headline_strict.py
```

Python standard library only, no downloads or external package imports. The
entire primary source path is **282 lines** including comments and blanks:
native_function.py(8), model.py(60), reader.py(69), experiment.py(145), plus the
explicit thirty-edge seams.json and TASK.json. The reader takes only settings,
current tapes and the declared vacuum datum. The 232-line independent verify.py
uses floor/ceiling endpoints and a separately implemented reference-current
anchor. Wrappers only check hashes and replay these local files.

PREREG.json predates the first new-protocol run. Prior adjacent-level theorems
and tomography were already known and are disclosed; this is an extension and
new execution, not a claim of blind mathematical discovery. Later verification
and diagnostic registrations explicitly identify the already known main result.
The complete source/body/probe records are in EVIDENCE.json.gz. The exact
scientific claims and prior-review dispositions are in DESIGN.md, PROOF.md,
CLAUSES.md and RESPONSE_TO_REVIEW.md; the proposed cell has62words.

The supplied vacuum/capacity architecture does not force unit selection: the
same prepared source preserves doubles under the control. It does set the
admitted preparation domain and absolute reference. Uniform two over the ENTIRE
graph remains fixed and lacks that vacuum preparation. Integer current output,
continuity and equivalent load-vector topology receive no quantization credit.
No fractional-state, energy/spin join or new physical topological-charge
invariant is claimed. The corresponding release checklist is maintained in a
separate row review directory, with required other-group review still open.
