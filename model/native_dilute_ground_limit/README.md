# Dilute ground observations

A fixed number of source occupations becomes invisible to bounded observers
on a fixed number of ports as capacity grows. The exact ground theorem bounds
both connected spectral weight and the complete source-plus-record output of
finite local measurement protocols. For norm at most K and connected weight
at least w, visibility requires s*q/V >= w/(4*K^2). Every failure and
coherence is retained.
The density modes use the source Laplacian frequency. Finite-density states,
growing observables and admitted charged fields require distinct continuum
constructions rather than an identification from these dilute observations.

[PROOF.md](PROOF.md) gives the assumptions and bounds;
[results.json](results.json) binds the literal source and exact endpoint checks.

```sh
python -m model.native_dilute_ground_limit.packet --check
python -m pytest -q model/native_dilute_ground_limit
```
