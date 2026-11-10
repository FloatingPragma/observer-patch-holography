# Certified numerical replay of the protected readout

The complete protected-cube experiment uses 681,283 source pulses per setting.
Small differences in six-level matrix exponentials accumulate in its floating
state. The retained raw probabilities and complete NumPy archive are preserved.
This package checks both that archive and a fresh execution against an
independent interval calculation of the same finite circuit.

The calculation proves an exact protected-subspace identity for the ideal
Strang circuit. Literal finite controls have a separate source-error bound;
their terminal repair probabilities are not set to zero. All failure labels,
padding checks, references, complex phases and source-clock conventions are
retained. The source, experiment and original numerical checker are unchanged.
The [proof](PROOF.md) specifies the numerical admission and its limits.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -m model.protected_readout_numerics --check
.venv/bin/python -B -m pytest -q -Werror model/protected_readout_numerics
```

The CLI executes the original cold producer. `--replay-output /path/run.json`
exports its complete fresh JSON/NPZ pair. The original command
`python -m model.protected_particle_readout.packet --check` retains its stricter
backend-sensitive comparison and can reject the same replay. No physical gate,
particle identification or hardware accuracy is inferred from this validator.

`--check --renew-live` renews only the two enumerated parent hashes in the
legacy live receipt after the complete producer and numerical checks pass. Its
original scientific fields and NPZ remain exact. `retained.json` preserves the
original JSON bytes. The current full-suite entry point is
`python -m model.verify_current`; the original verifier remains available.
