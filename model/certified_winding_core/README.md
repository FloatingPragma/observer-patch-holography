# Certified winding spectra and finite interacting poles

All 28 matrices in the [winding search](../winding_sector_search/README.md)
have exact rational spectral enclosures with every eigenvalue and its
multiplicity retained. The [proof](PROOF.md) controls source-matrix rounding
and the complete eigenbasis using exact integer residuals. Numerical
eigenvectors propose a certificate; they do not supply its error bound.
Binary64 limb products stay in a proved exact-integer range.

For every retained odd-degree cube, exactly two enclosed levels can be
zero. The independent inversion index forces both to be exactly zero and
exhausts the kernel, including the weak-mass cases outside the strong-mass
theorem. The six-level degree-three windows contain this pair and four
nonzero levels. Their index identifies no generations. Even-degree cases
retain certified invertibility or an unresolved kernel count.

The same symmetries and full-Fock contact theorem give interacting ground
and charged poles whenever the displayed inequalities hold. The supplied
contact grid is `1/65536, 2^-24, 2^-40`; inconclusive cases are retained.
The four-state comparison controls the full Fock complement. No matrix
with `2^(2V)` rows is claimed to be executed.

The [receipt](results.json) contains all intervals and contact cases.
Mass coefficients, textures, couplings and boundary conditions are source
inputs. An interacting relativistic limit, physical particle identification
and a defect formation law are not supplied by these finite certificates.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.certified_winding_core --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.certified_winding_core --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/certified_winding_core
```
