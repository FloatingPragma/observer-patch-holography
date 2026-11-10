# Spin and mass conditions for particle candidates

The declared routing field has an explicit spin-one-half action and an exact
finite charge current. Its free relativistic window has a controlled rotation
error. The source's fifteen-state internal representation has cancelling
anomaly traces and permits the three familiar Higgs mass channels. It forbids
an unbroken constant chiral mass and leaves the Higgs state and Yukawa
coefficients unspecified. The [proof](PROOF.md) separates these constructions
from native particle identification; the [receipt](results.json) checks their
matrices and complete finite current records.

If the map retains the full lifted `U(3) x U(2)` symmetry, all three Yukawa
channels carry its determinant character. A declared even charged-memory
transition cancels that character. Three complete sixteen-state controls verify the resulting
charge exchange and physical parity. This provides a concrete permitted
interaction architecture; its coupling, Higgs state and memory preparation
do not follow from the native neutral repair rule. Retaining only the
trace-free gauge group is a different permissible extension; this memory
is not a universal Higgs requirement.

```sh
python -m model.coherent_defect_wave.physical_attachment --check
python -m pytest -q model/coherent_defect_wave/physical_attachment
```
