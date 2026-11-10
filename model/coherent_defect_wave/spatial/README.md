# Spatial propagation between repair cells

Neighboring cells exchange coherent wave records while retaining their local
link defects. A received wave changes the next local repair and its complete
defect-position readout. The exact circuit preserves each interior gauge
constraint and has a finite propagation cone on cubic boxes of any size.
Its uniform triplet sector has a quadratic low-momentum band; a relativistic
particle interpretation requires a different justified propagation law.

The [proof](PROOF.md) defines the connected family, conserved current,
stationary product background and exact finite-volume comparison.
The [receipt](results.json) checks the two-cell experiment and the finite
translation symbol against the literal circuit.
The separate [Clifford routing construction](ROUTE.md) executes a finite
position-space walk with a controlled free Dirac limit and specifies a
local coupling to native defects. Its mass coefficient is a declared input.

```sh
python -m model.coherent_defect_wave.spatial --check
python -m pytest -q model/coherent_defect_wave/spatial/test_spatial.py
```
