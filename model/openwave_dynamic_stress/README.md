# Dynamic source stress and native clocks

The autonomous wave source supplies its own field force through a
positive-square-root Sylvester equation. The same graph action fixes
stiffness and inertia. Exact finite Ehrenfest equations retain the
oscillator boundary correction, and geometry/native/field energy rates
sum to zero. Complete stress and momentum records include every outcome.

All64 native vacuum sectors have interval-certified spectra. The seven
dressed native clocks have a strictly positive frequency-ratio spread of
about87.89913ppm; an exact single redshift factor does not apply at finite
probe strength. A spectral-gap theorem bounds the weak-probe deviation
from the background lapse expectation. Preparations and this auxiliary
probe domain are explicit in `PROOF.md`.

The result is a finite action/backreaction and clock-response theorem,
without physical G, a continuum Einstein stress tensor or laboratory
clock calibration. Force records are numerical; the native clock-gap
intervals have exact residual certificates.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_dynamic_stress --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_dynamic_stress --check
python3 -m unittest model.openwave_dynamic_stress.test_response
```
