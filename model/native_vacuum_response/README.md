# A positive-energy response from native records

The fair repair source fixes a transfer ground and a positive generator.
Its complete four-record experiments recover finite-time ground response
functions through an explicit complex linear combination of native lag
statistics. Blocking keeps the approximation and coefficient bounds uniform
in capacity at fixed correlation degree and counted time. The
[proof](PROOF.md) includes the actual preparation budgets, coefficient
amplification and finite-sample costs.

The recovered response is not an executed unitary or a prepared physical
vacuum. An exact phase-intervention witness shows why the transfer ground
cannot be identified with the native tracial working state under the same
field and intervention dictionary. A physical identification must supply
that operational map in addition to the recovered correlations.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_vacuum_response.packet --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/native_vacuum_response
```
