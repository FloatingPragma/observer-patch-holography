# Charge position after a finite energy measurement

The finite charged-energy detector is followed by a joint measurement of
charge position and whether the transverse field is excited. The
[proof](PROOF.md) bounds all twenty-six outcomes, including the finite
clock's ground-state error and original preparation probability. Every
rejected energy label remains recorded. These are observations of the
specified finite source, without a measured particle assignment.

```sh
python -m model.affine_gauss_clock_observables --check
python -m pytest -q -W error model/affine_gauss_clock_observables
```
