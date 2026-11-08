# Source representations and particle identification

The [particle map](PARTICLE_MAP.md) exports every Standard Model target channel
beside its computed internal source representation. Exact occupation records
recover all fifteen matter weights, including color conjugation, and the
charged weak and neutral Higgs directions. The [machine export](particle_map.json)
keeps unproved physical sources, masses and defect assignments empty. A common
interacting field, source-selected families and complete particle detectors are
required to identify these representations with laboratory particles.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.standard_model_source_map --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.standard_model_source_map --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/standard_model_source_map
```
