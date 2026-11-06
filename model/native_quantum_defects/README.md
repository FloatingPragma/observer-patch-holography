# Native modes and the defect search

Native repair counts give a one-port excitation of the reconstructed state
with unchanged exterior, positive energy and an exact contrast decay law.
They also determine all 66 transfer modes: an exhaustive search finds no
amplitude mode confined to nine or fewer ports. The [proof](PROOF.md)
distinguishes these count-derived constructions from preparing a protected
physical particle; no native phase operation or mass identification is claimed.

```sh
python -m model.native_quantum_defects.certificate --check
python -m pytest -q -W error model/native_quantum_defects
```
