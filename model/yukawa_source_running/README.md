# Higgs thresholds and Yukawa running

The complete finite Higgs source has an exact Pauli-blocking mechanism:
charged-addition energies depend approximately on squared Yukawa singular
values. Both spectral lines and every blocked creation outcome reconstruct
the original coupling, including an unknown common reference offset. The
[proof](PROOF.md) and [calculation](results.json) retain the full 438- and
1,596-state carriers. Hopping controls show that the tiny addition energies
depend on the declared interaction and preparation.

The same source matrices also define a conditional one-loop running
calculation. Minimal Standard Model running contracts charged-lepton
hierarchies toward the infrared. The [separate comparison](comparison.json)
finds that neither this running nor the finite addition energies supplies
the observed charged-lepton hierarchy. The source-to-field matching and
the scalar, coupling and reference choices require a physical construction.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/yukawa_source_running
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.yukawa_source_running.certificate --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.yukawa_source_running.certificate --check
.venv/bin/python -m model.yukawa_source_running.compare --check
```

`spectral_records.npz` retains complete spectra, amplitudes, residues,
creation vectors and blocked effects for every declared finite control.
No empirical particle file enters the producer. `compare.py` reads the
explicit empirical common-scale packet after the source calculation.
