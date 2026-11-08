# Reading particle-shell candidates without discarding failures

Complete energy-clock and spatial records bound how much of the original
source response lies near a proposed dispersion surface. The bound includes
clock resolution, preparation and simulation error, and finite-shot
uncertainty. On the selected eight-site spin source, an exact full-spectrum
certificate gives a positive finite excitation weight in one spatial-symmetry
sector. The [proof](PROOF.md) states the uniform growing-volume and growing-time
conditions needed to turn such readouts into a sharp continuum component;
the finite observation does not prove those conditions or a measured mass.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.spin_shell_readout --check
python -W error -m pytest -q model/spin_shell_readout
```
