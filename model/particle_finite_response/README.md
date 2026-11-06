# Finite particle response

A finite repair process determines recorded correlations after its state and
detector are specified. This package checks the native neutral response,
conserved passive-repair zero, and the distinction between a scalar mass and a
family potential. It also gives an exact comparison criterion for complete
observation histories at different capacities. These finite results constrain
particle identification; they do not supply the physical vacuum, clock,
matching convention or measured spectrum.

The [proof](PROOF.md) gives the domains and finite-stack connections. The
[receipt](results.json) binds the source code and records exact algebra alongside
explicit numerical controls.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_finite_response.certificate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_finite_response.certificate --check
python -m pytest -q model/particle_finite_response
```

The package imports the native family memory and passive Schur calculation
directly. It uses neither measured masses nor the calibrated mass producer.
The neutral hierarchy bound applies to one specified amplitude reference and
common mass scale, with `0 <= p <= 1/12` and `0 < h <= 1/2`. Floating residuals
are numerical checks, not measurement uncertainties or interval certificates.

## Coverage

| Result | Evidence in this package |
|---|---|
| Complete-record comparison across capacities | Exact outcome-wise criterion, induction proof, error bound, mutation controls |
| Neutral finite response | Native projector/transition join, exact quartic, proved finite-step hierarchy bound, retained detector backaction |
| Passive family metric | Direct reuse and replay of native harmonic extension and conservation certificate |
| Scalar/potential and finite translation obstruction | Proofs with exact matrix controls |
| Quantum-link occupation boundary | Complete exact CAR counterexample |
| Full compact gauge band, polar transport, coordinate-link completion | Complete public replay in [particle_gauge_band](../particle_gauge_band/README.md); separate from this compact replay |
| Full all-lag memory correction | Inherited continuation result; scalar projection applies conditionally to its coherent potential |
| Physical spectrum identity | Requires selected preparation, field/reference/clock dictionary, and complete-record comparison; no such identification is asserted |
