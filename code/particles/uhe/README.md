# UHE Coefficient Emission Receipts

This folder mirrors the high-energy messenger coefficient-emission simulator
contract inside the paper-stack code tree. It is a synthetic algebra scaffold,
not a source receipt, particle-spectrum prediction, or UHE source detection
claim. The default coefficients are planted examples. No source artifact,
baseline, source moments, minimality/interiority calculation, refinement
comparison, or coefficient solve is supplied by this command, so every such
readiness gate remains unverified and public promotion is always false.

Run:

```bash
python3 particles/uhe/build_uhe_coefficient_emission_receipts.py \
  --output particles/runs/uhe/coefficient_emission
```

The generated bundle illustrates the required shapes for:

- source-release quotient
- source law
- compact-engine source loads
- baseline measure
- feature map
- source moment targets
- coefficient solver
- emitted coefficient examples
- source dependency DAG and no-UHE-data-use receipt
- claim ladder

Event coordinates, event energies, residual maps, likelihood values, diagnostic
overlays, or post-event catalog choices invalidate the source-only label.
Passing the no-target-data scan is necessary but not sufficient for a source
receipt. A future producer must bind a real source artifact, compute every
readiness receipt, run the declared solver, and add independent verification
before it may emit `SOURCE_ONLY_COEFFICIENT_EMITTED`.

## Conditional universal pair threshold

`universal_pair_threshold.py` evaluates the leading head-on photon--photon
pair-production threshold when the hard photon, electron, and positron share
the same negative dimension-six coefficient. It retains the
Lorentz-invariant-lepton photon-only branch as a separate comparator. The
calculation is a conditional kinematic attachment: it assumes a preferred
frame, additive energy--momentum conservation, physical particle
identifications, and a negligible soft-photon correction. It derives no
interaction vertex, optical depth, propagated flux, carrier scale, or
empirical verdict.

Run the producer and focused tests from the repository root:

```bash
python3 code/particles/uhe/universal_pair_threshold.py \
  --output code/particles/uhe/runtime/universal_pair_threshold.json
python3 -m pytest -q code/particles/uhe/test_universal_pair_threshold.py
```

The JSON uses analytic or bracketed roots rather than grid onsets,
distinguishes a true no-window result from a window outside a scan range, and
labels every external limit by the assumptions needed to translate it into
the OPH scale coordinate. It also distinguishes the unconstrained AM--GM
lower bound from the minimum over the physical share domain: below the
transition the latter is the equal-share endpoint, while above the transition
the two coincide. Every receipt binds the producer and the corresponding Lean
source by SHA-256.
