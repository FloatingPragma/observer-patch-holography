# Measuring a filtered physical response

A coherent finite time filter applies to the selected source's authenticated
charge projections and Wilson plaquette contractions, as well as its
boundary-charge parity probe. Every selector, neutral probe and Ramsey outcome is retained,
including filter failures and the untouched quantum reference. Their
original scores read the filtered Gram correlation without normalizing by
success. The [proof](PROOF.md) gives source-bound quadrature, clock and
compiler errors, and an exact small matrix execution checks the apparatus.

A separate complete Lüders construction implements the common field
algebra's binary effects. Its effects differ from coherent-filter success
effects. The proof bounds every finite measurement history, including all
outcomes and quantum references. Two-level exact controls check this
distinction; the three-level control supplies a positive original spectral
window using the inherited contrast bound.

The finite construction uses declared coherent apparatus services. Its
spatial matrix and compiled tape are not executed, and no surviving
primary continuum excitation or particle mass is established.

```sh
python -m model.particle_filtered_probe --write
python -m model.particle_filtered_probe --check
python -m pytest -q model/particle_filtered_probe
```
