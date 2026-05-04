# Issue #223 Snapshot

URL: https://github.com/FloatingPragma/observer-patch-holography/issues/223
State checked: OPEN on 2026-05-04
Title: [OPH Task] Close the Ward-projected Thomson endpoint theorem gap on the P-closure branch
Labels: enhancement, surface: particles

## Summary

The synthesis-paper P derivation fixes the pixel ratio by outer/inner
self-consistency, but the remaining exact-closing gap is the low-energy Thomson
endpoint.

## Current Boundary

The code and paper cleanup removed the default imported Thomson/CODATA endpoint
from the live D10 electroweak source-transport readout. The default artifact
leaves Thomson endpoint comparison fields null unless an external value is
supplied explicitly as compare-only metadata.

The required closure is still the theorem-grade Ward-projected zero-momentum
transport law from the D10 source branch, including endpoint normalization,
quadrature/remainder control, and a clear statement of theorem-path versus
continuation/readout model.

## Acceptance Criteria

- Define the exact theorem-grade object replacing the imported Thomson endpoint
  on the Ward-projected U(1)_Q lane.
- Derive the low-energy transport law from OPH data, or isolate the first
  scalar/object that still cannot be internalized.
- State whether the charged-spectrum package, screening law, or endpoint
  normalization is the actual blocker.
- Remove theorem-path dependence on imported Thomson endpoint constants from
  the live derivation path.

