# Repeated centered scalar readout

This package constructs a binary quantum instrument for the same 64-site
scalar action, original-vacuum preparation and recovered reference-clock
intervals used by the existing scalar packets. The observer-like detector
has local pointer state, field couplings, local readback, a retained record
algebra and an explicit backaction map.

Repeated centered readouts multiply the original detector response by an
exact product of cross-time commutator cosines. The certificate retains all
21 unconditional marginals, includes all 210 prior-readout pairs, and bounds
their difference from the **continuously evolving field with the same
repeated instrument**. All 15 resolved time slots of the parent readout are resolved here.

The packet does not contain sampled quantum outcomes. It defines their joint
instrument mathematically; the reported marginals average over earlier
outcomes and are not independent trials or branch-conditioned probabilities.
Controlled quantum gates, pointer resets, Born readout, the field vacuum
and its physical attachment are supplied. The clock intervals are inherited
conditional calibration data, not a clock reconstructed from pointer bits.

A GHZ pointer compilation uses only the detector half of the **scalar action
graph**, with 31 declared entangling edges. The verifier checks its stabilizers,
local coupling coefficients, readout parity and ideal operation counts. This
graph is not identified with W12 seam wiring. The circuit is instantaneous in
the comparison model; its finite-duration and noise errors are not certified.

From the repository root:

```bash
python3 code/source_scalar_instruments/sequential_instrument.py --check
python3 code/source_scalar_instruments/verify_sequential_instrument.py
python3 -m pytest -q code/source_scalar_instruments/test_sequential_instrument.py
```

The independent verifier rebuilds the full action using the existing
independent edge-energy assembler, propagates first-order phase-space kicks
(the producer uses a second-order recurrence), and reruns the complete parent
quantum/clock verification. `verify_arithmetic` is an internal test helper;
its result explicitly does not claim that parent replay occurred. No parent
receipt is rewritten. The mathematical proof is in
`paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex`.
