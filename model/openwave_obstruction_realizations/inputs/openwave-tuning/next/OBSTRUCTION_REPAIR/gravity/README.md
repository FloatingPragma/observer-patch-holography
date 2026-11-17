# Einstein realization: scoped obstruction and responding-regulator candidate

The reported Dirac construction still does not realize the Einstein tower.
The conclusion that all fixed-lattice routes are excluded is too broad.
This additive packet contains an independently reproduced negative, a precise
correction to an exact-zero claim, and a new explicit local candidate.

| Problem | Work completed | Status |
|---|---|---|
| F087 says no q² response on the finite-mass Dirac source | Independent exact spectral derivative, continuous-q check, reproduction of the producer's L32 mode, and an exact positive witness in the same source family | Exact-zero claim corrected; correct gravity coefficient not obtained |
| Fixed-site entropy substituted for fixed-proper-volume variation | Explicit same-topology metric calculation with both curvature signs and a failing fixed-coordinate control | Geometric bookkeeping solved; microscopic entropy attachment open |
| Original conformal handles leave the ultraviolet state unresponsive | New finite-range quartic proper-regulator/coframe candidate; actual 3D slab entropy responds with the geometric sign and exponent approaching 2 | Constructive candidate, declared dynamics |
| Area response might destroy the modular bulk relation | Same candidate's transverse-zero low-energy modular response approaches the CFT coefficient; coframe clock and numerical controls included | Partial compatibility; full 3D/ten-observer premise open |
| Newton field and strength | No matched stationarity mechanism or same-source G normalization supplied | Unsolved |
| Full continuum Einstein tower and Lambda | Proper-cutoff/refinement issue exposed; finite residual theorem reusable in a stated mesoscopic window | Unsolved |

`PROOF.md` gives the exact statements and remaining obligations. The quadratic
Wilson attempt is retained as a failed area-scaling control. Existing Fable
connectivity, anisotropic-handle, and typed-coin producers are not duplicated.

Reproduction, from the workspace root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/OBSTRUCTION_REPAIR/gravity/audit.py
reverse-engineering-reality/.venv/bin/python openwave-tuning/next/OBSTRUCTION_REPAIR/gravity/exact_witness.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/OBSTRUCTION_REPAIR/gravity/regulator_candidate.py 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/OBSTRUCTION_REPAIR/gravity/regulator_candidate.py 4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/OBSTRUCTION_REPAIR/gravity/candidate_checks.py
```

These are development/source diagnostics. No reserved acceptance run is used.
No official or internal gravity row is promoted. `RESULT.json`,
`EXACT_WITNESS.json`, `REGULATOR_ORDER4_RESULT.json`, and
`CANDIDATE_CHECKS.json` are the receipts. `REGULATOR_RESULT.json` records the
initial order-2 attempt; the script was then generalized with an order argument.
The first candidate-check run reached serialization and failed on a NumPy
boolean; its log is preserved, the serialization was repaired, and the complete
unchanged numerical calculation was rerun.
