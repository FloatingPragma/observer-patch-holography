# Clocks and waves on the reciprocal geometry

The source-generated finite geometry changes complete scalar-wave detector
records by several parts per million. Its transverse metric ray time has
a positive delay, while its radial null speed stays one. Exact native
clock calculations also expose a missing physical condition: different
local transitions fail to share the same lapse ratio.

The [proof](PROOF.md) distinguishes the original eight-level clocks from
the declared spatial scalar-wave extension. It retains the finite wave
cutoff, every detector outcome and the scalar geometry's failure of the
tested vacuum Einstein balance. These are model-coordinate observations,
without a physical gravitational-redshift or light-delay identification.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_metric_observables --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_metric_observables --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/openwave_metric_observables
```
