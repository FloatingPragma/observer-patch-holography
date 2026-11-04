# Higher-winding defects on the supplied source

The full sixteen-row continuum action has chiral index **2d** for a
prescribed degree-d mass texture. A degree-three texture therefore has
at least six zero modes, including all eight lattice valleys. This is a
net chiral count; it does not identify three Standard Model generations.
The [proof](PROOF.md) gives the exact matrix map to the Callias operator
and states the analytic theorem it uses.

On the original lattice, every unit, inversion-odd texture at
`m>sqrt(3)*v/a` instead has exactly two zero modes. This extends the bounded
radial theorem to every odd degree. Its inversion grading differs from the
continuum internal grading, and its strong-mass condition fails in the
fixed-mass continuum limit.

The finite search retains all eigenvalues of 28 complete matrices. It
starts with degrees 1–4, radii 1–3 and masses 1 and 4. Four explicitly
identified follow-up calculations extend degrees 1 and 3 to radii 4 and 5
at mass 1. With `a=v=m=1` and the fixed window `|E|<1/4`:

| Degree | Radius | Window rank | Nonzero window energies | Boundary projector trace |
| --- | ---: | ---: | --- | ---: |
| 1 | 4 | 2 | none | 0.0235444 |
| 1 | 5 | 2 | none | 0.00455718 |
| 3 | 4 | 6 | ±0.1657820223, each twice | 0.385777 |
| 3 | 5 | 6 | ±0.1451903447, each twice | 0.154790 |

Each row includes two numerical zeros. The degree-three low-energy
subspace becomes less concentrated on the boundary, but two volumes do
not establish its limiting energies or interacting particles. At mass 4,
the odd-degree cases satisfy the exact two-zero strong-mass theorem.

Degree, mass, texture dynamics and interactions remain inputs. This
calculation supplies neither an action-selection rule nor a current
intertwiner assigning these modes to quarks or leptons. A protected
single-particle zero mode is not a calculated massive-defect rest energy.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.winding_sector_search --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m model.winding_sector_search --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -W error -m pytest -q model/winding_sector_search
```

The [receipt](results.json) binds the source and proof files, preserves
every trial and distinguishes exact algebra from floating diagnostics.
