# The executable model

Everything here is fixed for every candidate law.

| File | What it fixes |
| --- | --- |
| [carrier.py](carrier.py) | one cell, the echosahedron: 12 ports, 30 seams and 20 faces with the incidence of the icosahedron; the federations of 20, 80 and 320 cells, in which every cell is complete and external seams join three ports of each cell; the emitter seam of each carrier; the checks of these facts |
| [carriers/](carriers/) | the federation arrays `FEDERATION_LEVEL{0,1,2}.npz` and `CELL.npz`, one cell with its seams and its 120 wedges |
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
