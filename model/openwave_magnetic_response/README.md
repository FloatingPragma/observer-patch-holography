# Magnetic response from the native kinetic matrices

Squaring the complete gauged native kinetic operator produces its magnetic
spin-curvature term. The tree-level coefficient is `g=2` in the smooth,
nonrelativistic limit, without an inserted Pauli term or fitted magnetic
moment. Finite wave-packet energy and coherent phase measurements converge
toward that value with explicit matrix-function error bounds. All native
spin and valley rows are retained through exact internal reduction.

This is the response of a declared gauge extension to an external field;
a self-sourced dipole and physical electron identification are separate
requirements.

[Source, equations and scope](PROOF.md) · [Executed results](results.json)

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_magnetic_response --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_magnetic_response --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_magnetic_response
```
