# Certified bath mixing and cooling lag

The unchanged declared bath has rigorously enclosed relaxation gaps on
all 26 components at four supplied cost temperatures. Complete cooling
laws give certified distances to their componentwise Gibbs references.
The [proof](PROOF.md) separates component weights, the bath attempt clock,
matrix arithmetic and retained history.

```sh
python3 -m pytest -q model/native_bath_mixing
python3 -m model.native_bath_mixing --write
python3 -m model.native_bath_mixing --check
```
