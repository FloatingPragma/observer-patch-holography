# Local energy of the shared-link source

A specified positive Hamiltonian combines the actual onsite repair unitary
with the shared-link spatial hopping. Its constrained empty-wave states are
exact grounds, and its finite source has populated energy bands and a
nonadditive two-wave response. The [proof](PROOF.md) states the clock choice,
the exact spatial calculation and its finite-time detector error bound.
The [receipt](results.json) binds the source; no GeV mass is assigned.

An exact six-state invariant sector supplies a finite collision calculation:
two waves prepared one per cell acquire a joint double-occupation response.
All six source columns, both local constraints, the full occupation readout
and preparation failures are retained. A countercontrol distinguishes the
positive-clock contact interaction from independent hopping. This is a
finite-time detector result, without an asymptotic scattering interpretation.

```sh
python -m model.coherent_defect_wave.energy --check
python -m pytest -q model/coherent_defect_wave/energy
```

The continuous clock differs from the discrete circuit. Its physical
relativistic dispersion and particle interpretation require their own checks.
