# Particles and holes of a graded flux field

A specified fermionic occupation law carries the literal paired-flux code
through a local Wilson-Dirac Hamiltonian. Its filled-sea reference has positive
particle and hole energies, a source-unit gap of two, complete preparation
records and a controlled massive free-field limit. The grading, regulator
and reference are construction choices. Interactions and identification with
measured particle species require a source-compatible attachment.

The [proof](PROOF.md) states the actual link encoding, finite ground,
antiparticle meaning and continuum bounds. The
[receipt](results.json) records complete detector outcomes and source hashes.

```sh
python -m pytest -q model/coherent_defect_wave/many_body_route
python -m model.coherent_defect_wave.many_body_route --check
```

`spectrum.symbol(a,p)` and `energy(a,p)` use physical momentum p;
`symbol_lattice(a,k)` uses dimensionless lattice momentum k=a*p.
`projectors(a,p)` returns positive and negative bands. The same negative
projector is the filled-sea covariance. All five literal face channels have
the same four-component symbol.
