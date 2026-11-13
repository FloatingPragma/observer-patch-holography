# Active six-state charge memory

A coherent memory coupled to actual native load currents has six exactly
protected charge sectors. Two native loops act on the memory with coprime
order-two and order-three shifts. All 396 source/memory coordinates and all
25 raw outcomes are retained, with a tensor proof and complete column controls
on a 240-port source.

The memory voltage is declared architecture. The calculation identifies its
exact conservation law and the conditions under which memory errors destroy
it. Its modular no-equal-halves statement depends on a chosen class. Finite
centre labels alone cannot define a nonzero additive real electric charge.
Physical electric units, a localized quark string and rest masses are not
outputs of this package.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_active_native_charge --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_active_native_charge --check
python3 -m pytest -q model/openwave_active_native_charge
```

The [proof](PROOF.md) states the complete source domain and normalization
boundary. Actual coherent service and source-quadrature execution belong to
`openwave_native_charge_readout`.
