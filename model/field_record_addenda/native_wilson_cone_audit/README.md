# Wilson dispersion audit and an explicit repair candidate

The actual fixed-spacing source has matter/photon squared-dispersion
ratio two along each coordinate axis. Merely enlarging its box does not
remove it. `PROOF.md` distinguishes this energy curvature from a signal
velocity, and gives two explicit repair routes.

`cone.py` checks the literal Clifford matrices, the proposed orthogonal
Wilson regulator, charge conjugation, spin rotations, full bond norms
and a finite gauge-background gap-square identity. The changed regulator
is a new source member, not a recalibration or a proof of interacting
Standard Model masses. The old source and its receipts are unchanged.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_wilson_cone_audit --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_wilson_cone_audit -W error
```
