# Three complete finite observer readouts

The [proof](PROOF.md) attaches full-source operator witnesses to the periodic
modular Y readout, joint event-clock readings and the seven-mode temporal
record recurrence. The [receipt](results.json) retains original preparation
failures, timeouts, quantum references and source bindings.

These are proved and executed finite model laws. The coherent apparatus,
preparation and scheduler are declared parts of each construction. The
modular phase is an active target service; it is not native attempt time or
physical duration. Exact correlation laws establish the seven-mode rank;
finite sampled correlations do not.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/finite_observer_readout_closure
OPENBLAS_NUM_THREADS=1 python -m model.finite_observer_readout_closure --write
OPENBLAS_NUM_THREADS=1 python -m model.finite_observer_readout_closure --check
```
