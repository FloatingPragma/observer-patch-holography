# Regional quantum energy and entropy

Restricting the declared interacting scalar source to one fermion family
produces an exact modular pair interaction. Its spectrum, complete local
probabilities and a controlled entropy response follow from the same
prepared source state. The calculation retains every local occupation and
the source's conjugate reference. These internal regions supply explicit
modular data; identifying a spatial diamond kernel and geometric area
requires a spatial region and metric construction.

See the [proof](PROOF.md), [executable](regional.py) and
[same-source preparation](../particle_observer_outputs/THERMAL_CLOCK.md).

```sh
python3 -m pytest -q model/particle_regional_geometry/test_regional.py
python3 -m model.particle_regional_geometry.regional
```
