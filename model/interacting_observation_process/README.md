# Interacting dynamics with measurements between observations

The compact 180-mode memory source admits a controlled comparison for
complete experiments with measurements and feedback between evolution
intervals. Its actual memory continues across each detector operation;
every native and detector record is retained. The [proof](PROOF.md) gives
an explicit accuracy schedule, valid with an arbitrary quantum reference
and uniformly in the exterior volume. It extends the finite interacting
source's observation law, without constructing a spatial continuum theory.

A separate complete 16-row nonlinear probe executes two noncommuting
detectors with feedback using native and detector records. Its exact
12-state memory retains all 512 fine labels per input through a lossless
record factorization. The probe is an independent check of the comparison,
not a reducing sector of every compact family pulse. Detector compilation,
physical clocks and an interacting chiral continuum require additional
constructions.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.interacting_observation_process --check
python -W error -m pytest -q model/interacting_observation_process
```
