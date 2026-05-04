# Issue #153 Snapshot

URL: https://github.com/FloatingPragma/observer-patch-holography/issues/153
State checked: OPEN on 2026-05-04
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

Surrogate artifacts are not promotable. The QCD/Thomson bundle should either
define the exact production artifact needed for rho_had(s;P) and its uncertainty
budget, or isolate the smallest missing QCD primitive.

