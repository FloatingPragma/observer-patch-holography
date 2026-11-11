# Quark binding in quenched SU(3)

Labels: EXACT, MEASURED (Monte Carlo, with error), ARGUED, IMPORTED (with citation), DECLARED. The measured numbers
below are stored in [results.json](results.json), which is recomputed from
[recorded_tables.json](recorded_tables.json). Prefixes of its 68 named checks: `1_gamma`, `2_gauge`, `3_dirac` and
`4_hadrons` for the fast gates of section 6, `5_potential`, `5_gauge`, `6_solver`, `7_hadrons`, `8_chiral` and
`9_volume` for section 7, `10_melting` for section 8, `11_recorded_gates` for the gates of the compiled kernels.
Errors: a single bracket is statistical. With two brackets the first is statistical and the second systematic. A
number after `+-` is the two added in quadrature.

## 1. The system

- The model's. The gauge group of colour is `SU(3)`, the colour factor of the group of port phases. The colour
  factor of the seam cost is the Wilson plaquette action `S = beta sum_P (1 - Re tr V_P / 3)`, `beta = 6/g^2`
  (EXACT, `model/colour_deconfinement`, section 1). Quanta sit in the triplet of the three colour levels, so a quark
  carries one colour index.
- DECLARED. The links on the seams are continuous. A law that writes the port phases is open.
- Chosen here. The Wilson fermion is a standard discretisation of one Dirac species per flavour. It is not an
  operator of the model.
- Chosen here. `beta = 5.7` is a working point of the literature. No value of `g` at the cell is used.
- Open in the model. The quark masses. The runs use four hopping parameters, the chiral limit and one measured
  ratio, `m_pi / m_rho`, to fix the light quark mass at the physical point.
- Approximation. Quenched: the quark determinant is set to one, so quark loops are dropped. The weak and
  hypercharge factors and the `Z6` identification are left out.

## 2. Gauge configurations

The kernel is `model/colour_deconfinement/su3_numba.py`: Cabibbo-Marinari heat bath over three `SU(2)` subgroups
and subgroup overrelaxation. One compound sweep is one heat-bath sweep and four overrelaxation sweeps. Lattices
`8^3 x 16` (A) and `12^3 x 24` (B) at `beta = 5.7`, periodic. 32 independent streams per lattice, even streams from
a cold start and odd streams from a hot start, 500 compound sweeps of thermalisation, 50 between saved
configurations. Lattice A has 5 configurations per stream, 160 in total. Lattice B has 4 per stream, 128 in total.
Errors are jackknife errors over the 32 streams.

## 3. Static potential

- Smearing (Albanese et al. 1987): APE smearing of the spatial links, `V = (1 - alpha) U + (alpha/4)` times the
  sum of the four spatial staples, projected to `SU(3)`, `alpha = 0.5`, 12 steps. Time-like links stay thin, so the
  transfer matrix in time is the original one.
- Wilson loop `W(r, t)`: a spatial transporter of smeared links, two time-like lines of `t` thin links and the
  transporter back. Off-axis transporters average the shortest lattice paths closest to the straight line. All
  orientations are averaged.
- `W(r, t) = sum_n c_n(r) exp(-V_n(r) t)` with `c_n >= 0`, because the same operator creates and annihilates the
  string. Estimator `V(r) = ln[W(r, 2) / W(r, 3)]`.
- Primary fit `V(r) = V0 + sigma r - e/r`, three parameters, all separations with `r >= sqrt 2`, uncorrelated,
  weights from the jackknife errors. Alternative fits: all separations; `r >= 2` with `e = pi/12`; the estimator
  `ln[W(r, 3) / W(r, 4)]` with the points whose error is below 10 per cent. The systematic error of `sigma` is half
  the spread of the four fits. Sommer scale `r0 = sqrt((1.65 - e)/sigma)` (Sommer 1994).

## 4. Quarks and correlators

- Wilson-Dirac operator `M = 1 - kappa H`,
  `H psi(x) = sum_mu [(1 - gamma_mu) U_mu(x) psi(x + mu) + (1 + gamma_mu) U_mu(x - mu)^dag psi(x - mu)]`
  (Wilson 1975). Hermitian gamma matrices in the Dirac basis, `gamma_4 = diag(1, 1, -1, -1)`. `gamma_5 M gamma_5 =
  M^dag`. Antiperiodic boundary in time, periodic in space.
- Solver: BiCGStab on the even-odd reduced system to a relative residual `1e-10`, with conjugate gradient on the
  normal equations as fallback. The true residual of the full system is computed with the plain operator.
- Source: the origin, Gaussian smearing `F = (1 + s^2/(4n) Lap)^n` with the covariant Laplacian of the APE links,
  `s = 3.0`, `n = 20`. Twelve solves per hopping parameter. Primary correlator: smeared source, point sink. Check:
  smeared source, smeared sink.
- Pion `psibar gamma_5 psi`. Rho `psibar gamma_i psi`, averaged over `i`. Nucleon
  `eps_abc (u_a^T C gamma_5 d_b) u_c` with the projectors `(1 +- gamma_4)/2` and degenerate quark masses. The two
  Wick pairings of the two up quarks are both kept. Zero momentum by the sum over each time slice.
- Folding: mesons `C(t)` and `C(T - t)` averaged. Nucleon: the `(1 + gamma_4)/2` correlator at `t` averaged with
  the `(1 - gamma_4)/2` correlator at `T - t`.
- Hopping parameters 0.1600, 0.1650, 0.1663, 0.1675, four of the values of Butler et al. 1994 at this coupling.

## 5. Analysis

Every derived number is recomputed on every jackknife sample: potential, fits, masses, chiral fits, ratios.

- Mass fits: mesons `A [exp(-m t) + exp(-m (T - t))]`, nucleon `A exp(-m t)`, uncorrelated least squares.

| Lattice | pion | rho | nucleon | alternative windows |
| --- | --- | --- | --- | --- |
| A, `T = 16` | 4 to 8 | 4 to 8 | 4 to 7 | 5 to 8, 5 to 8, 5 to 7 |
| B, `T = 24` | 5 to 12 | 5 to 10 | 5 to 9 | 6 to 12, 6 to 10, 6 to 9 |

- Systematic error of a mass: the larger of the differences primary window against alternative window and point
  sink against smeared sink.
- Chiral fits, linear, four points: `m_pi^2 = A (1/kappa - 1/kappa_c)`, `m_rho = r_0 + r_1 m_pi^2`,
  `m_N = n_0 + n_1 m_pi^2`. Variants: the alternative windows, the smeared sink, the three lightest hopping
  parameters. The systematic error of a chiral number is its largest deviation among the variants.
- Physical point: the quark mass where `m_pi / m_rho = 0.1780`, the measured ratio (IMPORTED, Particle Data Group
  2024).
- `sigma` and `r0` come from the primary potential fit of the same lattice inside the same jackknife.
- The windows and the forms were fixed before the runs. No configuration is removed. One rule was added after
  exceptional configurations were seen in a test on part of lattice A and before the analysis of either lattice: a
  hopping parameter is dropped from the chiral fits of a lattice when the jackknife error of its pion mass exceeds
  25 per cent of the mass, or when the folded pion correlator is not positive in the fit window.

## 6. Code and gates

[lattice.py](lattice.py) is the reference in numpy. [production/kernel.py](production/kernel.py) holds the compiled
kernels of the production. The fast gates run in `--check` on small lattices with the reference, on a fixed random
stream:

| Fast gate | Result | Reference |
| --- | --- | --- |
| Clifford algebra, `gamma_5`, charge conjugation | below 1e-14 | 0 |
| plaquette at `beta = 1`, `4^4`, against `u + 4 u^5` | 0.05922(72) | 0.06013 |
| links in `SU(3)`; plaquette and smeared loops under a gauge transformation | below 1e-12 | 0 |
| `gamma_5` hermiticity of the operator, relative | 1e-15 | 0 |
| true residual of the solves | 9.8e-11 | bound 1e-9 |
| gauge covariance of the propagator | 1.4e-11 | 0 |
| free field, `4^3 x 16`, `kappa = 0.120`, against the momentum-space form: propagator, pion, rho | 5e-14, 1.7e-12, 1.8e-12 | 0 |
| pion, rho and nucleon correlators under a gauge transformation: point, smeared source | 1.7e-10, 4.8e-11 | 0 |
| nucleon contraction against the brute-force Wick sum | 7e-15 | 0 |
| free nucleon, `6^3 x 32`, `kappa = 0.110`: effective mass at `t = 10` | 1.3062 | three quark masses 1.3060 |
| free pion, same lattice: effective mass at `t = 10` | 0.87064 | two quark masses 0.87064 |

The residual-size entries are one run. They vary at their own size between platforms and are held to the bounds of
the named checks.

Gates of the compiled kernels, recorded with the production (key `production_gates`):

| Gate | Result |
| --- | --- |
| compiled hopping term against the reference with explicit gamma matrices | 2.4e-16 |
| `gamma_5` hermiticity: operator; propagator from two sources | 4.9e-16; 8.8e-10 |
| even-odd solve: true residual; difference to conjugate gradient on the normal equations | 7.9e-11; 2.9e-10 |
| free field, `4^3 x 16`: pion and rho against the momentum-space sum; propagator | 3.1e-12; 8.3e-14 |
| gauge covariance at `kappa = 0.15`: propagator, mesons, nucleon | 1.7e-11, 4.5e-10, 1.4e-12 |
| smeared correlators under a gauge transformation: point sink, smeared sink | 3.9e-11, 8.9e-11 |
| nucleon contraction against the brute-force Wick sum | 7.7e-16 |
| free nucleon, `8^3 x 32`, `kappa = 0.110`: effective mass against three quark pole masses | 1.3075 against 1.3060 |
| APE links in `SU(3)`; smeared loops gauge invariant; loops equal the routine of the gauge kernel | 3e-15; 1.7e-15; 3e-15 |

`test_hadron_binding.py` compares the compiled hopping term with the reference when numba is installed.

## 7. Results at zero temperature

MEASURED throughout.

- **7.1 Configurations.** Plaquette 0.54928(15) on A and 0.54911(8) on B. IMPORTED: 0.5492 (Lepage and Mackenzie
  1993). Integrated autocorrelation time of the plaquette: 3.5(3) and 3.2(2) compound sweeps. Hot and cold starts
  agree after 23 and 25 sweeps. Lag-one correlation of the plaquette between saved configurations: -0.11(9) and
  0.05(10). Solver on B: 91, 193, 263 and 381 iterations on average at the four hopping parameters, true residual
  at most 1.4e-10, no fallback. On A one solve at `kappa = 0.1675` used the fallback. No configuration is excluded.

- **7.2 Static potential, lattice B.**

| Separation | r | V(r) |
| --- | --- | --- |
| (1,0,0) | 1.000 | 0.5054(4) |
| (1,1,0) | 1.414 | 0.6828(7) |
| (1,1,1) | 1.732 | 0.7786(9) |
| (2,0,0) | 2.000 | 0.8073(10) |
| (2,1,0) | 2.236 | 0.8771(13) |
| (2,2,0) | 2.828 | 1.0004(20) |
| (3,0,0) | 3.000 | 1.0211(21) |
| (2,2,2) | 3.464 | 1.1235(31) |
| (4,0,0) | 4.000 | 1.2041(29) |
| (3,3,0) | 4.243 | 1.2573(34) |
| (5,0,0) | 5.000 | 1.3863(58) |
| (3,3,3) | 5.196 | 1.4238(68) |
| (4,4,0) | 5.657 | 1.4960(69) |
| (6,0,0) | 6.000 | 1.5533(107) |

| Fit | `a^2 sigma` | e | chi2 / dof | `r0 / a` |
| --- | --- | --- | --- | --- |
| primary, `r >= sqrt 2` | 0.1667(17) | 0.204(6) | 598 / 10 | 2.946(10) |
| all separations | 0.1407(13) | 0.369(2) | 1145 / 11 | 3.018(12) |
| `r >= 2`, `e = pi/12` | 0.1664(12) | 0.262 | 152 / 9 | 2.889(10) |
| estimator `t = 3` | 0.1677(35) | 0.206(12) | 290 / 10 | 2.934(20) |

  The potential rises at every on-axis step. Result: `a^2 sigma = 0.1667(17)(135)`, 12 combined errors from zero,
  and `r0 / a = 2.946(10)(64)`. No fit has an acceptable chi2: the points are precise enough to show the breaking
  of rotational symmetry at this lattice spacing, and the three-parameter form does not describe it. Lattice A
  gives `a^2 sigma = 0.1744(36)(184)` and `r0 / a = 2.913(20)(95)`.

- **7.3 Hadron masses, lattice B.** Smeared source, point sink.

| kappa | `a m_pi` | `a m_rho` | `a m_N` | `m_N / m_rho` |
| --- | --- | --- | --- | --- |
| 0.1600 | 0.6897(18)(8) | 0.8109(26)(9) | 1.3065(91)(79) | 1.611(9)(8) |
| 0.1650 | 0.4584(21)(4) | 0.6626(40)(15) | 1.027(12)(16) | 1.550(17)(26) |
| 0.1663 | 0.3809(34)(4) | 0.6231(64)(24) | 0.936(16)(25) | 1.502(26)(44) |
| 0.1675 | 0.232(32)(1) | 0.582(19)(16) | 0.52(22)(13) | 0.89(38)(21) |

  The pion lies below the rho at every hopping parameter. The nucleon lies above the rho at the three heavier
  quark masses. At `kappa = 0.1675` the nucleon mass is not resolved. There the largest pion correlator at
  mid-window is 110 times the median and two configurations lie above ten medians. At the other three hopping
  parameters the largest is 1.6, 1.8 and 3.3 medians. Such exceptional configurations are a property of the
  quenched Wilson operator near `kappa_c` (Bardeen et al. 1998). The point enters the fits with its error.

- **7.4 Chiral limit, lattice B.** `m_pi^2` is linear in `1/kappa` with slope 1.401(10)(11) and vanishes at
  `kappa_c = 0.16919(6)(3)`, chi2 3.8 for two degrees of freedom. `a m_rho = 0.5440(70)(29) + 0.561 (a m_pi)^2`,
  chi2 0.3. `a m_N = 0.789(16)(70) + 1.091 (a m_pi)^2`, chi2 3.4. At the physical `m_pi / m_rho`:
  `kappa = 0.16899(6)(3)`, `a m_rho = 0.5494(70)(29)`, `a m_N = 0.799(16)(67)`.

| Ratio | Chiral limit | Physical `m_pi / m_rho` |
| --- | --- | --- |
| `m_N / m_rho` | 1.4495(344)(1274) | 1.4544(335)(1196) |
| `m_rho / sqrt(sigma)` | 1.3326(169)(71) | 1.3458(169)(71) |
| `m_N / sqrt(sigma)` | 1.932(41)(172) | 1.957(40)(163) |
| `m_rho r0`, `m_N r0` | 1.603(20)(9), 2.32(5)(21) | |

  Variants of `m_N / m_rho` in the chiral limit: 1.439(38) with the alternative windows, 1.408(43) with the smeared
  sink, 1.322(82) with the three lightest hopping parameters. The last one sets the systematic error.

- **7.5 Lattice A.** The rule of section 5 drops `kappa = 0.1650` and `0.1675`: at 0.1650 one configuration lies at
  100 medians, at 0.1675 eleven lie above ten medians. At 0.1600: 0.6849(54)(65), 0.7998(69)(44), 1.334(13)(31).
  At 0.1663: 0.357(19)(13), 0.627(17)(9), 1.041(54)(12). Two points give `kappa_c = 0.16879(34)(19)` and
  `m_N / m_rho = 1.65(15)(1)` with no test of the linear form.

## 8. Above the transition

MEASURED. Lattice `16^3 x 4`, 16 streams, 300 compound sweeps of thermalisation and 1500 measured sweeps per
stream, 24000 measurements, at `beta = 5.60` and `beta = 5.85`. The transition at `N_t = 4` is at
`beta_c = 5.6925(2)` (IMPORTED, Boyd et al. 1996; `model/colour_deconfinement` measures 5.69248(71) from its two
largest volumes). `F(r)/T = -ln <L(x) L(x + r)^*>` with `L` the Polyakov loop, volume averaged, all separations
(McLerran and Svetitsky 1981).

| | `beta = 5.60` | `beta = 5.85` |
| --- | --- | --- |
| `<abs L>` | 0.01070(5) | 0.21615(11) |
| `F/T` on axis, `r = 1, ..., 6` | 4.432, 5.923, 7.120, 8.001, 8.931, 9.765 | 2.883, 3.025, 3.054, 3.060, 3.062, 3.063 |
| slope of `F/T` over `r >= 2`, weighted line | 0.9896(51) | 0.00107(1) |
| separations in the line, largest `r` | 73, 10.2 | 115, 13.9 |

Below the transition the free energy of a static pair rises linearly. Above it `F/T` is constant within 0.005
from `r = 4` to `r = 12`. The slope falls by a factor of about 930.

## 9. Comparisons (NOT BLIND)

- `m_N / m_rho`. Measured: `938.92 / 775.26 = 1.211` (IMPORTED, Particle Data Group 2024, isospin-averaged
  nucleon). This lattice: 1.4495 +- 0.132 in the chiral limit, 1.8 errors above. The measured value was known
  before the run.
- String tension. With the conventional `sqrt(sigma) = 440 MeV` (IMPORTED convention from Regge slopes and
  quarkonia, about 5 per cent wide) the measured masses give `m_rho / sqrt(sigma) = 1.762` and
  `m_N / sqrt(sigma) = 2.134`. This lattice: 1.3326 +- 0.0183 and 1.9317 +- 0.177.
- ARGUED. The coarse lattice, the finite volume and the dropped quark loops each move these ratios. Their sizes
  are not computed here.

## 10. What is not established

- Masses in GeV. The light-quark mass is open in the model: the physical point is set by the measured
  `m_pi / m_rho`. The lattice spacing in GeV is open: no value of `g` at the cell is used.
- Quark loops are dropped.
- The Wilson fermion is a discretisation chosen here.
- One lattice spacing. No continuum limit is taken.
- The nucleon at the lightest quark mass is not resolved, and the fit over the three lightest quark masses lowers
  `m_N / m_rho` by 0.13.
- The potential fits have no acceptable chi2. The error of `sigma` is the spread of the fits.
- The mass fits are uncorrelated. Their chi2 is no test of the fit form.

## 11. References

- K. G. Wilson, Phys. Rev. D 10, 2445 (1974); in New Phenomena in Subnuclear Physics, ed. A. Zichichi (1977).
- M. Creutz, Phys. Rev. D 21, 2308 (1980). N. Cabibbo and E. Marinari, Phys. Lett. B 119, 387 (1982).
- A. D. Kennedy and B. J. Pendleton, Phys. Lett. B 156, 393 (1985).
- M. Albanese et al. (APE), Phys. Lett. B 192, 163 (1987).
- S. Guesken et al., Phys. Lett. B 227, 266 (1989).
- H. A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 (1992).
- R. Sommer, Nucl. Phys. B 411, 839 (1994), arXiv:hep-lat/9310022.
- F. Butler, H. Chen, J. Sexton, A. Vaccarino, D. Weingarten, Nucl. Phys. B 430, 179 (1994), arXiv:hep-lat/9405003.
- W. Bardeen, A. Duncan, E. Eichten, G. Hockney, H. Thacker, Phys. Rev. D 57, 1633 (1998).
- G. P. Lepage and P. B. Mackenzie, Phys. Rev. D 48, 2250 (1993), arXiv:hep-lat/9209022.
- G. Boyd et al., Nucl. Phys. B 469, 419 (1996), arXiv:hep-lat/9602007.
- L. D. McLerran and B. Svetitsky, Phys. Rev. D 24, 450 (1981).
- Particle Data Group, Phys. Rev. D 110, 030001 (2024), https://pdg.lbl.gov/.
