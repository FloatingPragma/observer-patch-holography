# A physical odd local energy spectrum

The complete native quantum cell has an exactly resolved neutral odd
response with 49 energy records. Its lowest visible state is a gauge
singlet of spin one-half and lies below the vacuum-connected minimum.
Charge conjugation pairs the two spin directions at centered charge zero.
This comparison constrains the choice of a many-cell reference; it does
not identify a measured particle or establish a global ground across
every local sector.

```sh
python -m model.coherent_defect_wave.chiral_joint.neutral_local_spectrum --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/neutral_local_spectrum
```

See the [proof](PROOF.md), [exact operator calculation](model.py) and
[complete energy receipt](results.json).
