# Patch notes

## Post-integration fail-closed hardening

The upstream pole algebra was retained, but its v1 promotion path was removed.
The original checker trusted fixture evidence booleans and recomputed the
physics from receipt-controlled values without comparing them to the fixture.
This allowed an unrelated self-consistent pole to pass while claiming physical
promotion.

The integrated v1 now:

- emits the explicit sign translation `Delta=-Pi`;
- treats caller evidence booleans only as unverified claims;
- fixes every external evidence gate to false and physical promotion to false;
- binds the complete numerical subject and provenance to the fixture;
- derives status and blockers;
- recomputes every redundant pole/readout and neutral diagnostic field;
- applies checker/schema tolerance caps; and
- rejects promotion flips, unrelated receipts, inflated tolerances, corrupted
  derived fields, and corrupted neutral diagnostics.

Physical promotion is delegated to a future aggregate verifier that must load
and validate independent, hash-bound evidence receipts.

## Corrected architectural error

The previous `survival-proof-3/src/pole_map.py` is a valid target-convention
map from mass-dependent-width Breit--Wigner parameters to an energy-pole
coordinate. It is not a theory pole map.

The previous `run_pipeline.py` is a conditional adapter from a split OPH/SMDR
input model to SMDR outputs. It does not expose or independently check the
finite-order inverse-propagator equation.

This patch adds that missing layer and keeps all three artifacts distinct:

```text
experimental convention map
conditional backend adapter
strict finite-order theory pole map
```

## Corrected finite-order reporting

The square-root `M,Gamma` obtained from a one-loop-truncated `s` are useful
physical coordinates, but their nonlinear relation contains kinematic powers
beyond one loop. Gauge/BRST comparisons therefore use the emitted strict
coefficients, not the square-root coordinates.

## Fail-closed neutral mixing

At strict one loop the Z root uses `Delta_ZZ^(1)(z)`. The package can store the
full neutral matrix, but it refuses to insert the one-loop-squared AZ/ZA term
into the one-loop root. A future two-loop package must include it explicitly.
