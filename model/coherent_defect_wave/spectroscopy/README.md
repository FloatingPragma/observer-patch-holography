# Reading a frequency from the interacting source

An exact finite clock reads the same unitary that couples the defect and wave.
A complete filter has probability `17/144` on its gauge-constrained input and
removes both exact flat phases. Every other filter outcome, clock bin and
preparation failure is retained. The [proof](PROOF.md) and
[receipt](results.json) give a positive non-flat phase weight without a dense
matrix calculation. Physical momentum and an energy-branch convention are
required before this readout can support a rest-mass claim.

A second complete experiment disturbs the explicit stationary background by
reading its five defect faces. Its forty filter/clock records retain total
active weight `283/7200`; every active output is orthogonal to the stationary
reference. This gives a populated source-relative excitation response with
no particle or ground-energy assignment.

```sh
python -m model.coherent_defect_wave.spectroscopy --check
python -m pytest -q model/coherent_defect_wave/spectroscopy
```
