# Finite primary response and its survival boundary

The actual selected Wilson source has a strictly positive finite
plaquette contrast, derived from a native loop-character matrix element
and the prepared state's faithful floor. A complete same-clock apparatus
can resolve its sign at explicitly chosen, level-dependent times and
precision. The [proof](PROOF.md) keeps every original outcome weight.

The exact empty and fully occupied Fock sectors also exist, but their
combined selected probability is at most `7/(627*2^(8n+7))`. Conditioning
their spectra cannot establish surviving primary particle weights.

The finite contrast bound does not prove uniform weight in a fixed
positive-energy band. It supplies no mass pole or GeV assignment.

```sh
python -W error -m pytest -q model/particle_primary_survival
python -W error -m model.particle_primary_survival --write
python -W error -m model.particle_primary_survival --check
```
