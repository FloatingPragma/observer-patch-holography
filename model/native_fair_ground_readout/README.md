# Preparing and measuring a native source excitation

The original fair source has an exact exchange ground. An admitted native
word prepares it, and a local occupation phase creates an excitation with
exchange energy `50/33`. The complete energy readout distinguishes it from
the matched control by a bounded score contrast above `3/50`, retaining
failures and residual outcomes. The [proof](PROOF.md) gives the preparation,
error bounds and one common counted clock. This is a finite excitation test;
a particle mass requires a stable sector and a relativistic dispersion law.

Exact regional ground states and local-history estimates extend the
mathematics to large capacities without allocating a global state. The
same ground differs from the faithful transfer reference: neither its
regional algebra nor the thermal reference can silently supply the vacuum
premises of a relativistic field theorem. Growing coherent preparation and
routed spatial refinement require separate source constructions.

The owner replays the ground word and both phase controls, certifies their
stored primitives, checks every terminal route and retains complete readout
errors. It authenticates the audited 66-word energy receipt by its exact
SHA and recomputes its live source map, spectral certificate and consumed
bounds. It does not rerun that full energy producer. Floating successful
vectors are diagnostics; the instrument bounds use exact rational errors.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_fair_ground_readout.packet --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_fair_ground_readout.packet --check
python -W error -m pytest -q model/native_fair_ground_readout
```
