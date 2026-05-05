# Issue #153 Snapshot

URL: https://github.com/FloatingPragma/observer-patch-holography/issues/153
State checked: CLOSING AS OUT-OF-SCOPE / COMPUTATIONALLY BLOCKED on 2026-05-05
Title: [OPH Task] particles.e.33a-unquench-hadron-branch-and-publish-systematics: Unquench the hadron branch and publish systematics
Labels: enhancement, surface: particles

## Goal

Promote the current hadron lane from debug/frozen execution scaffolding to a
production-grade public OPH hadron surface with explicit systematics.

## Success Rule

- The hadron backend is unquenched on the intended public branch.
- Continuum, volume, chiral, and statistics systematics are surfaced explicitly.
- Public hadron rows become reproducible outputs rather than
  execution-contract-frozen placeholders.

## Boundary For This Campaign

Surrogate artifacts are not promotable. The QCD/Thomson bundle has defined the
exact production artifact needed for rho_had(s;P) and its uncertainty budget.
The remaining missing primitive is a working OPH hadron backend such as
GLORB/Echosahedron. #153 and related compact-surface issue #157 are closed as
out-of-scope, not solved, until that backend emits production data and
systematics.
