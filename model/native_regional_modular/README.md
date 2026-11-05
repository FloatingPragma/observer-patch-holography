# A native clock for the interacting regional reference

The corrected modular generator of the prepared Weyl/contact plaquette
has a native source operation on its complete sixteen-state regional
algebra. At dimensionless modular parameter `s=1/2`, the operation uses
37,568 source pulses and has complete channel error below `1e-8` against
the actual prepared regional reference. Every processor export failure
is retained. The bound applies to arbitrary inputs and untouched quantum
references, including the exterior field and the preparation record.

The actuator uses all 39 exact Pauli coefficients of the
[regional modular polynomial](../regional_modular_expansion/README.md),
including its noncommuting boundary corrections. Its error combines
literal port-control rounding, a proved symmetric product bound and the
[certified preparation](../certified_reference/README.md) comparison to
the actual regional logarithm. The [proof](PROOF.md) fixes the modular
sign, scalar gauge, complete output domain and counted services.

The same operation continues all 256 original preparation records
without a source reset or success conditioning. The default replay
executes all sixteen regional input columns and checks their complete
output. A separate full-preparation check executes the original 256
conditional words, then continues their recorded columns and retains
both preparation and actuator aborts.

This is an operational regional clock for a specified finite source.
It does not identify regional flow with the global modular restriction,
calibrate physical time, or supply a geometric stress charge, area law
or Einstein stationarity.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_regional_modular --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_regional_modular
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_regional_modular --prepared-check
```

`--write` regenerates the exact receipt in a fresh interpreter. Tests
include the complete native preparation and therefore take longer than
the all-input actuator replay. Raw matrix residuals are numerical
diagnostics; the rational control certificate supplies the source-level
error bound.
