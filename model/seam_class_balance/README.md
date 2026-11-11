# Class balance of the seam link

The electroweak step of the working chain has the form `v = E_cell b <A>`.
`A` is the flip of one seam by the class-3 central element. On a
continuous weak link it is the shift of the link variable by the centre
element `-1` of `SU(2)`, diagonal in the electric sectors with eigenvalue
`(-1)^(2j)`. The class balance is

    <A> = sum_j p_j (-1)^(2j),

with `p_j` the sector law of the seam link. This package computes `<A>`
exactly for single registers and on the Cartan circle, and measures it by
Monte Carlo for `SU(2)` on three complexes with two plaquette weights and
with a doublet.

Every result is written `<A> = M exp(-c S)` with `S = 2 pi/(4 alpha)`. The
working law of the chain is `c = 1` with `b M = 1`.

| Ensemble | c | M at `1/alpha = 24.3166` | Status |
| --- | --- | --- | --- |
| one free link, one tick, Cartan circle | 1 | 2 | exact |
| one free link, one tick, full link (weights `d^2`) | 1 | `2 (2S - 1) = 150.79` | exact |
| one-sided register (weights `d`), labels (weights 1) | no exponential | 0.0338, 0.167 | exact |
| Cartan circle, cubic lattice in `d` dimensions | `1 - 2/d` | 2 | exact |
| Cartan circle, icosahedral surface as space, many ticks | 0.398231 | 2 | exact |
| `SU(2)` heat kernel, `6^4`, `8^3`, icosahedral space | 0.500, 0.335, 0.399 | 180, 147, 137 | measured |
| `SU(2)` Wilson, `6^4`, `8^3`, icosahedral space | 0.400, 0.320, 0.369 | 48, 173, 126 | measured |

The measured coefficients carry errors of 0.0005 to 0.0014 with the prefactor power one. On the Wilson `6^4`
lattice that fit gives 0.420 and is excluded by its chi-square. The table shows the fit with a free power there,
0.400(5).

The identity behind the measurement: with `A` read at the midpoint of the
electric step, `<A>` equals the ratio of two partition functions that
differ in the weight of one plaquette, the time-like plaquette of the seam,
whose argument is multiplied by the centre element. On the Cartan circle
the Villain form of the ensemble is Gaussian with integer windings and the
leading weak-coupling ratio is `2 exp(-S (1 - Pi))`, where `Pi g^2` is the flux variance of that
plaquette. A pinned plaquette has `Pi = 0` and gives the one-link law. A
plaquette inside a lattice shares the half turn with its surroundings and
has `c = 1 - Pi` below one.

The counts `8/d` and the pair weight `exp(-2S)` on the icosahedral surface
are densities of the plaquette flux at a half turn. They are exact
statements about the same ensembles and they are a different object from
`<A>`.

No lattice ensemble measured here has `c = 1`. The coefficient one belongs
to a link whose plaquette cannot move: one free link in its one-tick
thermal law, or a link held by a stiff doublet. The measured sector law of
a lattice link has weight power 2 and a heat time longer than `g^2/2`.

From the repository root:

```sh
python -m model.seam_class_balance --write
python -m model.seam_class_balance --check
python -m pytest -q -W error model/seam_class_balance
```

The exact tables, the fits and the 73 named checks are recomputed by these
commands from [measured_points.json](measured_points.json), which holds the
reduced Monte Carlo outputs. The check and the tests run without `numba`
in under a minute. The Monte Carlo itself needs `numba`:

```sh
python -m model.seam_class_balance.tasks
python -m model.seam_class_balance.staged tasks_prod.json 0 1 prod_out
python -m model.seam_class_balance.staged_face_step tasks_fm.json 0 1 fm_out
```

[PROOF.md](PROOF.md) states each result with its label and separates exact
computation, measurement and declared premises. The comparison with the
measured electroweak scale is in its last section and is labelled as not
blind.

- [exact.py](exact.py), [complexes.py](complexes.py): the exact laws and the complexes.
- [fits.py](fits.py), [mini.py](mini.py): the fits and a staged measurement in plain numpy.
- [engine.py](engine.py), [engine_face_step.py](engine_face_step.py): the compiled Monte Carlo kernels.
- [tasks.py](tasks.py), [staged.py](staged.py), [staged_face_step.py](staged_face_step.py): the task lists and
  the drivers of the production.
- [certificate.py](certificate.py) and [results.json](results.json): the derived numbers and the named checks.
