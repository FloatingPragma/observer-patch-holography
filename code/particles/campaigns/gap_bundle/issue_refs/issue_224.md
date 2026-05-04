# Issue #224 Snapshot

URL: https://github.com/FloatingPragma/observer-patch-holography/issues/224
State checked: OPEN on 2026-05-04
Title: [OPH Task] Update all particle-derivation code paths to use the derived P closure root
Labels: enhancement, surface: particles

## Scope

Every live particle-derivation code path should consume the live derived P root
instead of legacy reference values, frozen endpoints, or stale public defaults.

## Current Boundary

The immediate circularity cleanup has landed locally:

- code/P_derivation no longer carries built-in CODATA/NIST inverse-alpha or
  Thomson endpoint constants;
- external inverse-alpha enters only through explicit compare metadata;
- the particle source-transport readout has no default Thomson reference;
- OPH Lab no longer computes the displayed fine-structure row from a hidden
  measured endpoint/factor.

This issue remains open because the broader particle stack still needs one
canonical, certificate-backed live P_* source after the interval/uniqueness
certificate and zero-momentum endpoint theorem are closed.

## Acceptance Criteria

- There is one canonical programmatic source of the derived P value for the live
  path.
- All live particle-derivation scripts read from that source instead of
  hard-coded P inserts.
- Legacy reproduction, compare-only, and historical helpers are clearly named.
- The emitted report prints the derived P, the corresponding alpha object, and
  downstream key quantities from the same run.
- README/code docs state which paths are live, compare-only, or historical.

