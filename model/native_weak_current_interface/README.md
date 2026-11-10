# Weak-current interface and scalar background

The exact full-CAR map uses actual L and conjugated e^c modes and gives
weak-current action only on its designated left components, independently
of scalar state. The quantum five-state scalar cannot be replaced by a
fixed nonzero eigenbackground without a strictly bounded residual. This
does not exclude collective expectations or the proved same-H dressed
elementary/Higgs energy channel.

The [proof](PROOF.md) also supplies an actual gauge-invariant Higgs-dressed
weak-current density with full operator norm two, and an exact copy-symmetry
control for generation claims. The [receipt](results.json) replays the
native current map and scalar matrices. Physical Lorentz chirality and
weak-current form factors require their same-source dynamical readouts.

```sh
python -B -m model.native_weak_current_interface --write
python -B -m model.native_weak_current_interface --check
python -B -Werror -m pytest -q model/native_weak_current_interface
```

Reusable exact operators are returned by
`model.native_weak_current_interface.model.dressed_weak_polynomial()`
and `dressed_weak_matrix()`. The first uses the original native CAR/scalar
polynomial convention; the second is the complete eighty-state active matrix.
