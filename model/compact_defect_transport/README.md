# Protected compact defects and their transport

A specified local texture-repair law preserves both compact winding charges on
every retained record. Its four localized fermion modes have uniform gap and
leakage bounds. An explicit band transport preserves the same native currents
and pairing operator and has a finite compilation into local quadratic gates.
The controls supply work and preserve the band at record endpoints; autonomous
defect motion and physical particle identification require additional dynamics.

The [proof](PROOF.md) states the complete law and its bounds. The
[receipt](results.json) contains every proposal branch and complete gate-replay
controls on the 686-coordinate source.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.compact_defect_transport --check
python -m pytest -q model/compact_defect_transport
```
