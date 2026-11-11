# Colour deconfinement

The colour factor of the link system is pure `SU(3)` gauge theory with the Wilson plaquette action
`S = beta sum_P (1 - Re tr U_P / 3)`, `beta = 6/g^2`. At finite temperature its links deconfine. The order parameter
is the Polyakov loop, and the symmetry is the colour centre `Z3`, the flux classes 2 and 4 of the defect register.
This package states the centre symmetry in exact form, records a Monte Carlo study of the transition, and gives the
transition temperature in units of the string tension from the same lattice.

| Quantity | Value | Label |
| --- | --- | --- |
| Colour factor of the seam cost `1 - Re tr_5(U)/5` | `(3/5)(1 - Re tr W / 3)` | EXACT |
| Centre transformation on one time slice | action unchanged, `L -> exp(2 pi i/3) L` | EXACT |
| Classes charged under the colour centre | 1, 2, 4, 5 | EXACT |
| Transition coupling at `N_t = 4`, three sizes | `5.69105(39)` | MEASURED |
| Same, two largest sizes | `5.69248(71)` | MEASURED |
| Published value, same theory | `5.6925(2)` | IMPORTED |
| Growth of the susceptibility maximum, `N_s^x` | `x = 2.555(29)` | MEASURED |
| `a sqrt(sigma)` at `beta = 5.6925` on `12^4`, 980 measurements | `0.4071(120)` | MEASURED |
| `T_c / sqrt(sigma)` at `N_t = 4` | `0.614(18)`, statistical 0.011, systematic 0.014 | MEASURED |
| Published value, same theory, same spacing | `0.5982(34)` | IMPORTED |
| Published value, continuum limit | `0.629(3)` | IMPORTED |
| Register ensemble, `Z3` | two-peak energy histogram, `C_max ~ L^1.35(5)`: first order | MEASURED |
| Register ensemble, `Z2` | no resolved second peak, `C_max ~ L^0.30(4)`: continuous | MEASURED |

The simulation holds the links alone, at one lattice spacing, as an ensemble scan in the coupling. The published
values belong to the same theory, so the comparison is not blind. With light quarks the transition is a crossover.
A temperature in MeV exists only through an imported ratio and is listed in [PROOF.md](PROOF.md) section 8 with
that label.

From the repository root:

```sh
python -m model.colour_deconfinement --write
python -m model.colour_deconfinement --check
python -m model.colour_deconfinement --demo
python -m pytest -q -W error model/colour_deconfinement
```

`--check` replays the exact statements, the kernel gates on small lattices and the consistency of the recorded
tables, and compares with [results.json](results.json). `--write` rebuilds `results.json` from the modules and the
recorded tables. `--demo` scans five couplings on a `6^3 x 2` lattice. With the numpy kernel it takes a few
minutes.

- [PROOF.md](PROOF.md) states each result with its label and citation.
- [su3.py](su3.py): numpy kernel and observables.
- [su3_numba.py](su3_numba.py): the compiled kernel that produced the tables. It needs numba. The check and the
  tests run without it.
- [centre.py](centre.py): exact centre statements.
- [register.py](register.py): the `Z_N` register ensemble.
- [scale.py](scale.py): scale arithmetic.
- [lattice_tables.json](lattice_tables.json): the recorded Monte Carlo tables. The raw series are not stored.
- [certificate.py](certificate.py) and [results.json](results.json): the derived numbers and 42 named checks.
- `test_colour_deconfinement.py`: gates. They run in under a minute.
