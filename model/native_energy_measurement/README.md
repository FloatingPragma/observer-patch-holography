# Measuring finite transfer energy

Existing coherent source services implement a destructive energy measurement
for the fair twelve-port transfer. Sixty-six complete filters retain every
energy label, processor failure and final residual. Exact source residuals
certify the energy dictionary and the complete measurement error. The
[proof](PROOF.md) distinguishes this controlled experiment from autonomous
repair evolution and laboratory particle masses.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_energy_measurement.packet --check
python -m pytest -q -W error model/native_energy_measurement
```
