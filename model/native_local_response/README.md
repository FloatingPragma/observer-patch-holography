# Local repair records and boundary effects

The native repair program has an exact local observation rule: only its
backward dependency cone can affect selected local records and the final
local state. For an actual connected cell, a 15-port halo gives a four-state
calculation of the boundary error, including the cell's own repair clock.
Arbitrary quantum references and measuring back-action are retained.
Adaptive FIFO signals require their full controller dependencies; the
result does not identify physical spacetime or particle masses.

[Proof](PROOF.md), [replay](packet.py), [checked result](results.json).

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_local_response.packet --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_local_response
```
