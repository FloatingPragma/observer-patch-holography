# Frozen task

Candidate: CAR-current-loop-clock-v1. Source, observation domain and complete
clock are fixed in DESIGN.md. No seeds, optimization, held-out selection or
empirical calibration enter: enumerate all finite columns, all stated times
and all diagnostic bands. q runs over0,1 in that order. L=1. Loops run forward
and reverse. Stationary phase signs are + then -, except endpoints have one.

Freeze DESIGN.md, TASK.md, probe.py and source pins before evidence execution.
No existing evidence may be overwritten. Any failed execution is retained
with its exact code, freeze and log before a separately registered correction.

Numerical gates: complete permutation/inverse and integer Gauss identities
are exact; source matrix identities, reference normalization, Chebyshev
recurrence, stationary residuals and electric balance each <=1e-9; bounds on
lazy and common amplitude spectra are checked with tolerance1e-9. Changed
law full-reference difference must exceed1e-6 for omitted-loop and zero-current
controls. The circulation reader must detect a planted ± oriented cycle at
error<=1e-12 and its zero-reader mutant must fail by>1e-3.

Report every stationary state, spin, realized current, returned-loop and
nonreturned-loop electric increment, coherent J quadrature and axial
transport moment. Nonzero diagnostic threshold is1e-8. A positive threshold
crossing is a finite stationary circulation result only. Zero or adverse
results remain in the receipt and are not a reason to change the recipe.

The source loop state and reference outcomes are retained unconditionally.
Basis-word intermediate audits are separate diagnostic preparations; no
intermediate pointer measurement is inserted into the coherent clock.
The random-unitary loop channel retains all eight branch addresses; its
amplitude kernel is not asserted to generate that channel's time evolution.

The primary replay is read-only, compares every retained array and canonical
receipt, and must run in under60seconds. Machine-dependent elapsed time is
printed, not part of the immutable scientific receipt.
