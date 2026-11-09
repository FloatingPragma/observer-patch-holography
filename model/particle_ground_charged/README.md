# Charged readouts from the complete finite ground

Charged probes act on the certified ground of the entire finite scalar
source. Their complete spectra contain multiple lines with calculated
operator weights; they differ from the response of the prepared 307-state
sector. Integer checks certify the energy intervals and probabilities,
including dark directions and the original addition/removal weights.
These are source-energy readouts. Identifying a particle mass requires
the same excitations' spatial propagation and relativistic transformation.

The [proof](PROOF.md) gives the full probe and error accounting. The
[executable](spectrum.py) also shows that a common number-energy shift
cannot turn these finite addition/removal spectra into identical conjugate
laws about this reference.

```sh
python -m model.particle_ground_charged
python -m pytest -q model/particle_ground_charged/test_spectrum.py
```
