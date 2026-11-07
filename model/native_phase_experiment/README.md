# Source preparation and readout of a one-port phase contrast

The twelve-bit occupation register represents the binary total-two native
source. Complete source words prepare the count-derived state
`rho=P/Tr(P)` and its port-zero phase conjugate. A local inverse beam and
the source computational readout measure the single-occupation exchange
operator on ports zero and one. Its full uniform-seam C and E records give
the contrast factors `5/6` and `7/10` in one native attempt.
Uniform independent seam selection is the supplied native scheduling
condition. The experiment contracts every seam's executed local instrument
with that classical weight.

The exact ideal readouts are `+/-1/168` before repair, `+/-5/1008` after C
and `+/-1/240` after E. The source word execution uses complex128 arithmetic;
its admission tolerance is `2e-10` in absolute probability and readout units
(`2e-20` for selector timeout). Admission against rational targets uses exact
subtraction of the represented floating-point values. Full density admission
compares the exact sum of squared real and imaginary deviations with the
squared Frobenius norm bound; it does not first round the target or take a
floating-point square root. These are
bounded numerical checks, with no certified outward rounding interval.

The preparation selects a seam and one of its 56 symmetric-basis vectors.
A source fair-bit selector accepts 1680 of 2048 words per draw, for at most
16 draws. Its timeout has probability `(368/2048)^16`. Successful outputs
retain their unconditional weight. Source export failures, raw aborts,
zero-weight records and all 36 terminal measurement labels are retained.
The declared proper-code controls, classical storage and computational
readout are apparatus services. The phase acts on a logical occupation bit.

```sh
python -m model.native_phase_experiment --check model/native_phase_experiment/results.json
python -m pytest -q -W error model/native_phase_experiment
```

`--write PATH` executes and independently verifies a fresh receipt, refusing
an existing path. `--archive PATH` retains every preparation factor and each
local source run's input, success, failure Gram and first-failure weights.
The emitted words and bound source reconstruct individual first-failure
outputs. Regional QR contractions keep every row and report their marginal
error. They trace the untouched ports and classical preparation archive
only for the declared marginal observations.
An archive name without `.npz` receives that suffix before destination
checks. Receipt and archive creation are exclusive.

The independent implementation reads the pinned source literals and counts
every endpoint, coin and raw record using exact rational arithmetic. It
imports no producer module. It compares both full prepared densities and
all 720 raw slots per preparation, including their pointer probability and
readout contribution. The exact phase-pair loss quotient invokes the
protected-obstruction classifier in
`Lean/ObservableNormalForms/tools/verify_protected_obstruction_models.py`.
The mathematical transfer context uses the SourceLaw results described in
`model/native_effective_fields/README.md`.
The receipt binds the experiment's local Python import closure, package
initializers and declared source data. It replays without Lean installed.
The standalone verifier checks the complete receipt profile, including
preparation slot labels and weights, selector bits, loss-quotient metadata
and source-word resource counts. The density and outcome targets come from
the separate rational implementation; source-word metadata is reconstructed
from the emitted services. Component weights also reconstruct a canonical
point/pair mixture whose total mass and full density must satisfy the same
global bounds and agree with the recorded preparation. Primitive selector
bits have a complete probability budget. The package tests check the retained receipt
against a fresh source execution.

The measured quantity is a finite transfer contrast in native attempts.
Physical particle identification requires the common protection, transport,
clock and energy conditions in the owning particle construction.
