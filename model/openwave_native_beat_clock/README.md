# Free native charge beats

An actual native charged pole and its retained reference define fixed
quadrature observables that oscillate under the unchanged interacting
Hamiltonian after one finite coherent preparation. No periodic drive or
reference rotation is applied during the clock evolution. Complete
full-H energy bins and source norm bounds control its frequency window,
waveform and sensitivity to initial-state or generator disturbances.

The ground-to-e^c beat reads `4+Delta_ec`. A second actual e^c-to-L
transition has one reference particle on both arms and reads
`Delta_ec-Delta_L`, cancelling their common reference energy. The L
transition has no measured-particle assignment. Preparation work, the
coherent services and the distinction from a physical Compton clock
remain explicit in the [proof](PROOF.md).

```sh
python -m model.openwave_native_beat_clock --write
python -m model.openwave_native_beat_clock --check
python -m pytest -q model/openwave_native_beat_clock
```
