# Momentum and quantum references

The full interacting three-dimensional source has exact momentum readouts
around a physical atomic preparation. On eight sites its source momentum has
weight 23/128 at zero and 15/128 in each other bin, while the retained quantum
reference makes their joint momentum sharp. A normalized neutral pump has a
complete success/failure readout and places the source on the odd-parity
momentum grid. These records do not determine a particle dispersion relation.

The [proof](PROOF.md) states the source and apparatus. The
[receipt](results.json) contains exact records at sides 2, 4 and 6.

Run from the repository root:

```sh
python -m model.coherent_defect_wave.chiral_joint.atomic_momentum_response --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/atomic_momentum_response/test_momentum.py
```
