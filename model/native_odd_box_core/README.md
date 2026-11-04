# A charged core without an imposed texture

The unmodified centered Weyl action on an odd open cube has two exact
zero modes. Its original repulsive contact interaction gives an actual
twofold ground state and a rigorously positive finite charged gap at the
declared weak coupling, with every Fock sector included. A complete native
circuit reads the four zero-mode occupations; the three-site version runs
on the full carrier with a retained quantum reference. The mode is an
open-boundary standing wave, whose gap bound vanishes in the controlled
growing-box family. Thermodynamic particle identification requires a
different limit.

The [proof](PROOF.md) gives the filled-sea interaction calculation,
actual-source ground and charged bounds, and native detector resources.
The [receipt](results.json) binds the source laws and complete replay.

```sh
python -m model.native_odd_box_core --check
python -m pytest -q model/native_odd_box_core
```
