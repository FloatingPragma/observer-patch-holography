# Maxwell behavior at weak fields

Smooth local field actions and smooth conservative dynamics can have
Maxwell theory as their leading response without being exactly quadratic
or linear. On the committed finite carrier, symmetry and the stated vacuum
susceptibilities select that response. Energy conservation supplies a
uniform finite-time error bound. An exact quartic example retains gauge
and Gauss constraints while allowing nonlinear corrections.

The [proof](PROOF.md) states the assumptions and
[receipt](results.json) records the exact checks. Native selection of the
field energy, its coefficients and the physical clock is a separate task.

```sh
python3 -m model.maxwell_weak_field.certificate --check
python3 -m pytest -q model/maxwell_weak_field
```
