# Exact source energies and their dispersion

The declared local energy law has exactly calculable three-dimensional
triplet branches. Their first positive source energies have ratios
`1:2.170820393…:2.618033989…`. The [proof](PROOF.md) identifies their
preparation, spectral weights and invariant reference, and shows why their
dispersion does not supply a common relativistic particle spectrum.

These are source energies, without a lepton or GeV-mass interpretation.

```sh
python -m model.coherent_defect_wave.dispersion --check
python -m pytest -q model/coherent_defect_wave/dispersion
```
