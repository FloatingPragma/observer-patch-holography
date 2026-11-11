# Quark binding in quenched SU(3)

The colour links of the model carry the Wilson plaquette action (`model/colour_deconfinement`, section 1). This
package puts Wilson quarks on those links without quark loops and measures three things at one coupling: the
potential between static colour charges, the masses of the pion, the rho and the nucleon at four quark masses, and
the free energy of a static pair below and above the deconfinement transition.

Lattice `12^3 x 24` at `beta = 5.7`, 128 configurations in 32 independent streams. All masses are in lattice units.
A single bracket is a statistical error. With two brackets the first is statistical and the second systematic.

| Quantity | Value | Label |
| --- | --- | --- |
| String tension `a^2 sigma` | 0.1667, statistical 0.0017, systematic 0.0135 | MEASURED |
| Sommer scale `r0 / a` | 2.946 | MEASURED |
| Critical hopping parameter `kappa_c` | 0.16919(6)(3) | MEASURED |
| `a m_rho` in the chiral limit | 0.5440(70)(29) | MEASURED |
| `a m_N` in the chiral limit | 0.7886(162)(703) | MEASURED |
| `m_N / m_rho` in the chiral limit | 1.4495 +- 0.132 (statistical 0.034, systematic 0.127) | MEASURED |
| `m_N / m_rho` at the physical `m_pi / m_rho` | 1.4544 +- 0.124 (0.034, 0.120) | MEASURED, input `m_pi / m_rho` IMPORTED |
| Measured `m_N / m_rho` (Particle Data Group 2024) | 1.211 | IMPORTED, comparison NOT BLIND |
| `m_rho / sqrt(sigma)` in the chiral limit | 1.3326 +- 0.0183 (0.0169, 0.0071) | MEASURED |
| `m_N / sqrt(sigma)` in the chiral limit | 1.9317 +- 0.177 (0.041, 0.172) | MEASURED |
| Slope of `F/T` for `r >= 2`, `beta = 5.6`, `N_t = 4` | 0.9896(51) per lattice unit | MEASURED |
| Same at `beta = 5.85`, `N_t = 4` | 0.00107(1) per lattice unit | MEASURED |

The static potential rises linearly, the pion mass squared falls linearly to zero at `kappa_c`, the rho and the
nucleon stay massive there, and above the transition the free energy of a static pair is flat. The nucleon mass at
the lightest of the four quark masses is not resolved ([PROOF.md](PROOF.md), section 7.3).

Not established: masses in GeV. The light-quark mass and the lattice spacing in GeV are open. Quark loops are
dropped. The Wilson fermion is a discretisation chosen here. There is one lattice spacing.

From the repository root:

```sh
python -m model.quenched_hadron_binding --write
python -m model.quenched_hadron_binding --check
python -m pytest -q -W error model/quenched_hadron_binding
```

`--check` runs the fast gates of the code on small lattices, repeats the analysis of the recorded tables and
compares with [results.json](results.json). It does not rerun the production.

- [PROOF.md](PROOF.md) states each result with its label and citation.
- [lattice.py](lattice.py): the lattice theory in numpy. It is the reference of the compiled kernel.
- [gates.py](gates.py): fast gates of the code on small lattices.
- [analysis.py](analysis.py): the fit chain from jackknife samples to the potential, the masses and the ratios.
- [literature.py](literature.py): imported numbers with their sources.
- [recorded_tables.json](recorded_tables.json): the reduced outputs of the production. Jackknife samples of the
  Wilson loops and of the hadron correlators, the Polyakov-loop correlators and the gates of the compiled kernels.
  The configurations and the raw measurement files are not stored.
- [certificate.py](certificate.py) and [results.json](results.json): the derived numbers and 68 named checks.
- [production/](production/): the compiled kernel and the run drivers. They need numba. The check and the tests
  run without it.
- `test_hadron_binding.py`: gates. They run in under a minute.
