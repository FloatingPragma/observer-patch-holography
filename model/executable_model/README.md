# The executable model

Everything here is fixed for every candidate law.

| File | What it fixes |
| --- | --- |
| [carrier.py](carrier.py) | one cell, the echosahedron: 12 ports, 30 seams and 20 faces with the incidence of the icosahedron; the federations of 20, 80 and 320 cells, in which every cell is complete and external seams join three ports of each cell; the emitter seam of each carrier; the checks of these facts |
| [validation.py](validation.py) | the live federation check counts external ports on every owner. The carrier command, candidate evaluator and contrast-memory verification use it. Frozen campaigns pin the carrier source, whose gateway summary counts owners appearing on external seams. |
| [carriers/](carriers/) | the federation arrays `FEDERATION_LEVEL{0,1,2}.npz` and `CELL.npz`, one cell with its seams and its 120 wedges |
| [refinement.py](refinement.py) | optional reader for a supplied Openwave L0–L7 refinement at its deposited SHA-256. It returns read-only labelled external seams and the shared cell incidence. |
| [refinement.py](refinement.py) | loads an explicitly supplied deposited spherical board at levels 0 through 7, verifies its identity and native external matching and returns read-only seam addresses |
| [symmetry.py](symmetry.py) | the sixty rotations of the cell as port permutations, and the incidence of ports, seams and wedges that decides which registers a local term may couple |
| [laws.py](laws.py) | the registry of candidate laws: every package in [candidate_laws/](../candidate_laws/CANDIDATE.md) with a `candidate.py` |

A port `12 k + j` of a federation is port `j` of cell `k`, and every cell keeps its 30 seams with the labelling of
`CELL_SEAMS`. The carrier is computational architecture: ports, seams and cells are no physical space, and every
test compares what observers record.

```sh
python -m model carrier                       # the checks of the fixed carrier
python -m model run --law NAME                # a candidate law on it
python -m model evaluate                      # every candidate against the requirements
```

Large runs can supply a file from Openwave revision `b189afa2420912c02e2ab5cee9808490ee6dd4d9`.
The reader checks the file's SHA-256 for the requested level, the three gateways per cell and, at
levels 0 through 2, every labelled external seam against the bundled federation.

```python
from model.executable_model.refinement import load_refinement

board = load_refinement("simulation/carriers/FEDERATION_LEVEL7.npz", level=7)
assert board.cells == 327_680 and board.external.shape == (491_520, 2)
```
