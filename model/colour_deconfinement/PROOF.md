# Colour deconfinement of the link system

Labels: EXACT (algebra or exact arithmetic, replayed by the named checks), MEASURED (Monte Carlo with a
statistical error), ARGUED (reasoning without a check), IMPORTED (published value, cited), DECLARED (premise).
The measured numbers below are stored in [results.json](results.json). Prefixes of its 42 named checks:
`1_action`, `2_centre` and `3_lattice` for sections 1 to 3, `4_strong_coupling` and `5_kernel` for section 4,
`6_tables` for sections 5 to 7, `7_scale` for section 8.

## 1. The colour factor of the plaquette action

The seam cost of a code state carried around a plaquette with holonomy `U` in `S(U(3) x U(2))`, averaged over the
five-level code `V = C^3 + C^2`, is `s(U) = 1 - Re tr_5(U)/5`. For a colour holonomy `U = diag(W, 1_2)` with `W` in
`SU(3)`,

\[
s = 1 - \frac{\mathrm{Re}\,\mathrm{tr}\,W + 2}{5} = \frac{3}{5}\Bigl(1 - \frac{\mathrm{Re}\,\mathrm{tr}\,W}{3}\Bigr).
\]

With the total action `(10/g^2) sum_P s(U_P)` the colour factor is the Wilson plaquette action
`S = beta sum_P (1 - Re tr W_P / 3)` with `beta = 6/g^2` (EXACT). The weak factor has `4/g^2`. A plaquette whose
holonomy is the central element of class `f` costs `s = 0, 17/10, 9/10, 4/5, 9/10, 17/10` for `f = 0, ..., 5`
(EXACT). For the colour classes 2 and 4 this is `3 beta / 2` in the Wilson normalisation. It is the cost of a thin
flux on one plaquette. The cost of the smooth flux quantum in `model/defect_particle_classes` is a different quantity.

Premises (DECLARED): the links on the seams are continuous elements of `S(U(3) x U(2))`, and their weight is the
exponential of the summed seam cost. A law that writes the port phases is open. The value of `g` at the cell is
not used. Every curve of this package is an ensemble scan in the coupling.

## 2. Centre classes

The hypercharge on the code is `Y = diag(-1/3, -1/3, -1/3, 1/2, 1/2)`. The central element of class `f` is
`exp(2 pi i f Y)`.

| Flux class `f` | Phase on colour levels | Phase on weak levels | Classes `k` whose loop turns |
| --- | --- | --- | --- |
| 2 | 1/3 turn | 0 | 1, 2, 4, 5 |
| 4 | 2/3 turn | 0 | 1, 2, 4, 5 |
| 3 | 0 | 1/2 turn | 1, 3, 5 |
| 1, 5 | 2/3, 1/3 turn | 1/2 turn | 1, 2, 3, 4, 5 |

Classes 2 and 4 are the colour centre `Z3`. Class 3 is the weak centre `Z2` (EXACT). A static quantum of class `k`
carried around a flux of class `f` acquires `exp(2 pi i f k / 6)`. For `f = 2` the phase is `k/3` of a turn, the
colour triality of the class. The quark classes 1, 2, 4, 5 are charged under the colour centre. Classes 0 and 3 are
neutral (EXACT).

## 3. Centre transformation, Polyakov loop and centre flux

Lattice `N_s^3 x N_t`, periodic, links `U_mu(x)`, time direction 4. The Polyakov loop is
`L = (1/N_s^3) sum_x tr prod_t U_4(x, t) / 3`. Let `z = exp(2 pi i / 3)`.

1. Global transformation (EXACT). Multiply every time-like link leaving one time slice by `z`. A plaquette in a
   plane `(i, 4)` based on that slice contains `U_4(x + i)` and `U_4(x)^dagger`, so the factors cancel. No other
   plaquette contains a link of the slice. The action is unchanged, and `L -> z L`.
   Checks: plaquette change and loop mismatch below `1e-12` on a thermalised lattice.
2. Order parameter. The transformation is a symmetry of action and measure. In a finite periodic volume the
   distribution of `L` is invariant under a rotation by a third of a turn (EXACT). The expectation of the loop of a
   static colour charge is `exp(-F/T)` with `F` its free energy (IMPORTED: Polyakov 1978, Susskind 1979, McLerran and
   Svetitsky 1981). `L = 0` means an infinite free energy: confinement. A nonzero `L` breaks the centre symmetry:
   deconfinement. The loops of the quark classes are the charged ones. Classes 0 and 3 have no such order parameter.
3. Transformation on a region (EXACT). Multiply the time-like links of one slice by `z` on a spatial region `R`
   only. The loops inside `R` turn by `z`. The loops outside stay. The plaquettes that change are the time-like
   plaquettes of the slice that straddle the boundary of `R`. Their holonomy is multiplied by `z` or `z^2`: a
   class-2 or class-4 flux is added there. On the cold lattice with a `3 x 2 x 4` block the flux sits on 52 plaquettes,
   one for each of the 52 links that cross the boundary. On every cube the six plaquette classes add to zero
   modulo 3: the flux surface is closed.
4. Reading (ARGUED). A sheet of class-2 or class-4 flux is the boundary between two centre sectors of the loop. A
   loop that links the sheet is turned by `z`. Percolating flux scrambles the phase of the loop and gives `L = 0`.
   Dilute closed flux leaves a net loop. The register ensemble of section 7 holds exactly this flux.
5. A transformation of the spatial links of a slice leaves the time-like loop unchanged (EXACT).

With matter in a charged class the action is not invariant. With light quarks the transition of quantum
chromodynamics is a crossover (IMPORTED: Aoki et al. 2006). The simulation of sections 4 to 6 holds the links alone.

## 4. Kernel and gates

Cabibbo-Marinari heat bath over the three `SU(2)` subgroups with the Kennedy-Pendleton or Creutz generator, and
subgroup overrelaxation. One compound sweep is one heat-bath sweep and four overrelaxation sweeps. The numpy kernel
in [su3.py](su3.py) updates links of one direction and one site parity together. The compiled kernel in
[su3_numba.py](su3_numba.py) runs in sequential order and produced the tables.

| Gate | Result | Reference | Label |
| --- | --- | --- | --- |
| gauge invariance of plaquette, Polyakov loop and Wilson loops | 2e-16 | 0 | EXACT |
| centre transformation, loop mismatch | 1e-16 | 0 | EXACT |
| overrelaxation, change of the action | 2e-16 | 0 | EXACT |
| unitarity drift after 200 compound sweeps | 8e-14 | 0 | MEASURED |
| plaquette, `8^4`, `beta = 5.7` | 0.54886(29) | 0.5492 (IMPORTED) | MEASURED |
| same, heat bath alone | 0.54934(22) | 0.5492 (IMPORTED) | MEASURED |
| plaquette, `12^4`, `beta = 6.0` | 0.59371(6) | 0.5937 (IMPORTED) | MEASURED |
| same on `8^4` (finite-size shift) | 0.59411(12) | 0.5937 (IMPORTED) | MEASURED |
| strong coupling, `6^4`, `beta = 1.0` | 0.06014(9) | 0.06013 (EXACT) | MEASURED |
| strong coupling, `6^4`, `beta = 2.0` | 0.12890(14) | 0.12877 (EXACT) | MEASURED |

Strong coupling (EXACT): with the Haar moments `<(Re tr U)^n> = 1/2, 1/4, 3/4, 15/16` for `n = 2, 3, 4, 5`, the
one-plaquette average is `u = beta/18 + beta^2/216 - 5 beta^4/93312 + O(beta^5)`. On the four-dimensional lattice the
first correction is the cube term, `<P> = u + 4 u^5` (IMPORTED: Drouffe and Zuber 1983).

The checks `5_kernel.*` and the tests repeat the exact gates and the strong-coupling gate on `4^4` with the numpy
kernel.

## 5. The transition at `N_t = 4`

Runs (MEASURED): 28 runs with 263300 compound sweeps in total on `8^3 x 4`, `12^3 x 4` and `16^3 x 4`. The first
15 per cent of each run is discarded.

| `N_s` | Compound sweeps in the reweighting | `beta_c(N_s)` | Susceptibility maximum |
| --- | --- | --- | --- |
| 8 | 78200 in 7 runs | 5.6881(10) | 1.314(14) |
| 12 | 64600 in 5 runs | 5.68952(26) | 3.60(5) |
| 16 | 63412 in 6 runs | 5.69123(39) | 7.90(15) |

`beta_c(N_s)` is the maximum of `chi = N_s^3 (<|L|^2> - <|L|>^2)`, found by multi-histogram reweighting of all runs of
one size. Errors are jackknife errors over ten blocks.

- Infinite volume (MEASURED). The first-order form `beta_c(N_s) = beta_c - h (N_t/N_s)^3` through the three sizes
  gives `beta_c = 5.69105(39)` with `chi^2 = 5.9` for one degree of freedom. Through the two largest sizes it gives
  `beta_c = 5.69248(71)`. Published values (IMPORTED): `5.6925(2)` in infinite volume and `5.6908(2)` on `16^3 x 4`
  (Boyd et al. 1996); `5.69226(41)` (Iwasaki et al. 1992).
- Order (MEASURED). The susceptibility maximum grows as `N_s^x` with `x = 2.555(29)` (pairs: 2.48, 2.74). A
  first-order transition has `x = 3` asymptotically. At `N_s = 16` the histogram of `|L|` at the maximum has two
  peaks, with a dip of 0.48 of the lower peak. The long runs pass between the two phases 30 to 45 times each.
- Sectors (MEASURED). At `beta = 5.80` on `12^3 x 4`, 6 independent runs (5 from random starts) settle at `|L|` near
  0.20: 2 in the sector of phase 0, 2 at a third of a turn and 2 at two thirds.

## 6. The transition temperature in units of the string tension

`T_c / sqrt(sigma) = 1 / (N_t a sqrt(sigma))` at the transition coupling is a pure number of the simulated theory.

- String tension (MEASURED). Lattice `12^4` at `beta = 5.6925`, 980 measurements. Spatial links smeared by
  12 APE steps with coefficient 0.4 and a gauge-covariant projection; time-like links thin. `V(R)` from
  `ln[W(R, 3)/W(R, 4)]`, fit `V = V0 + sigma R - e/R` over `R = 2, ..., 6`: `a sqrt(sigma) = 0.4071(74)`
  (statistical). Checks: `T = 2` gives 0.3981(29); 6 smearing steps give 0.4103(81); `e = pi/12` fixed over
  `R = 3, ..., 6` gives 0.4113(42); `R = 2, ..., 5` gives 0.4055(91). The largest difference is the systematic
  error.
- Result (MEASURED): `a sqrt(sigma) = 0.4071(120)` and

  `T_c / sqrt(sigma) = 0.614(18)` at `N_t = 4`, with 0.011 statistical and 0.014 systematic.

- Published lattice values for the same theory (IMPORTED). At this lattice spacing `a sqrt(sigma) = 0.4179(24)`
  (Boyd et al. 1996, Table 2), which gives `T_c / sqrt(sigma) = 0.5982(34)`, and `a sqrt(sigma) = 0.395(4)`
  (Teper 1998), which gives `0.633(6)`. In the continuum limit `0.629(3)` (Boyd et al. 1996) and `0.646(5)`
  (Lucini, Teper, Wenger 2004).
- The measured value lies 0.9 combined standard errors above the first same-spacing value. The comparison is not
  blind.
- Unsmeared loops (MEASURED). Creutz ratios from 500 measurements on the same lattice: `chi(2,2) = 0.3782(5)`,
  `chi(3,3) = 0.2555(14)`, `chi(4,4) = 0.200(13)`. They fall with the loop size towards the fitted
  `sigma a^2 = 0.1658(60)`. Key `string_tension` of `results.json`.

## 7. The register ensemble

The centre projection of the links is a `Z6` register. Its colour part is a three-dimensional `Z3` link field and
its weak part a `Z2` link field, each with weight `exp(kappa sum_p Re W_p)`.

| Group | `L` | Specific-heat maximum | Coupling of the maximum |
| --- | --- | --- | --- |
| `Z3` | 8 | 7.20(13) | 1.0668(8) |
| `Z3` | 12 | 12.1(4) | 1.0781(6) |
| `Z3` | 16 | 19.6(12) | 1.08085(17) |
| `Z3` | 24 | 31.3(23) | 1.08315(16) |
| `Z3` | 32 | runs stay on one branch (0 passages) | not determined |
| `Z2` | 8 | 3.566(28) | 0.7323(7) |
| `Z2` | 12 | 4.01(13) | 0.7455(8) |
| `Z2` | 16 | 4.40(14) | 0.7505(5) |
| `Z2` | 24 | 4.3(4) | 0.7583(7) |
| `Z2` | 32 | 5.9(4) | 0.75751(28) |

- `Z3` (MEASURED). The maximum grows as `L^x` with `x` = 1.27(10), 1.69(26), 1.16(24) between successive sizes
  and `x = 1.35(5)` in one fit over `L` = 8 to 24. The energy histogram at the maximum has two peaks. The dip
  between them, relative to the second peak, is 0.57, 0.41, 0.31, 0.20 for `L` = 8, 12, 16, 24. The extrapolated
  coupling is `kappa_c = 1.08387(20)`. At `L = 32` the runs of 42500 sweeps pass between the branches 0 times. The
  mean plaquette is 0.886 at `kappa = 1.083` and 0.897 at 1.084 from random starts, and 0.950 at 1.085 from an
  ordered start.
- `Z2` (MEASURED). The maximum grows as `L^x` with `x` = 0.29(8), 0.32(16), -0.07(24), 1.1(4) between successive
  sizes and `x = 0.30(4)` in one fit over `L` = 8 to 32, with `chi^2 = 4.0` for three degrees of freedom. The
  histogram has no resolved second peak: the same ratio is 0.93, 1.00, 1.00, 1.00, 0.97 for `L` = 8, 12, 16, 24, 32.
  The extrapolated coupling is `kappa_c = 0.7612(5)`.
- Binder minimum (MEASURED). The distance of the minimum of the energy cumulant from 2/3 is 8.5e-03, 3.9e-03,
  2.6e-03, 1.2e-03, 1.0e-03 for `Z3` and 4.0e-03, 1.3e-03, 5.7e-04, 1.6e-04, 9.2e-05 for `Z2`, at `L` = 8, 12, 16,
  24, 32.
- The `Z2` maximum at `L = 32` lies 2.8 standard errors above the value at `L = 24`. The exponent of the last pair
  alone does not separate the two ensembles. The order assignment rests on the histograms, the Binder minimum and
  the fit over all sizes. Checks `6_tables.register_*`.
- Duality values (IMPORTED spin-model couplings): `kappa_c = 1.0843` for `Z3` from the three-state Potts model
  (first order; Janke and Villanova 1997) and `kappa_c = 0.761413` for `Z2` from the Ising model (continuous,
  `alpha/nu = 0.17`).

The colour centre has a first-order transition and the weak centre a continuous one, in the register's own ensemble.

## 8. Scale

The simulation gives pure numbers. A temperature in MeV needs a scale.

- From the coupling chain (DECLARED forms, derivation open): `alpha_s(m_Z) = 0.11834`, which runs to
  `Lambda3 = 0.3400 GeV` (four-loop running and three-loop decoupling, IMPORTED).
- From lattice quantum chromodynamics with physical quarks (IMPORTED): `T_pc = 156.5(1.5) MeV` (HotQCD 2019) and
  `Lambda3 = 338(12) MeV` (FLAG 2021), ratio `0.4630(170)`.
- Product: `157.4 MeV`, with `5.8 MeV` from the imported ratio and `6.9 MeV` from the chain's threshold corners.
  Radiation-era time `1.24e-5 s` (`g* = 61.75`) to `2.35e-5 s` (`g* = 17.25`), which is `1.8e38` to `3.4e38` ticks.

This number exists only through the imported ratio. It is no computed result of the source, and it is not blind:
the product returns the imported temperature times the ratio of the two `Lambda3` values, `1.006`.

## 9. What is not established

- The simulation holds the colour links alone. Quarks, the weak and hypercharge factors and the `Z6`
  identification are absent.
- `N_t = 4` is a coarse lattice. No continuum limit is taken. The published continuum value of
  `T_c / sqrt(sigma)` lies 5 per cent above the published `N_t = 4` value.
- The volume exponent is measured on `N_s <= 16`. It is below 3. The three-size fit of `beta_c` has `chi^2 = 5.9`
  for one degree of freedom.
- The string tension comes from on-axis loops up to `R = 6` on `12^4`.
- The register tables are an ensemble scan in `kappa`. No law fixes `kappa`.
- The agreement with the published values tests the code and the identification of the colour factor with the
  Wilson action. It is no prediction.

## 10. References

- K. G. Wilson, Phys. Rev. D 10, 2445 (1974).
- A. M. Polyakov, Phys. Lett. B 72, 477 (1978). L. Susskind, Phys. Rev. D 20, 2610 (1979).
- L. D. McLerran and B. Svetitsky, Phys. Rev. D 24, 450 (1981).
- G. 't Hooft, Nucl. Phys. B 138, 1 (1978).
- B. Svetitsky and L. G. Yaffe, Nucl. Phys. B 210, 423 (1982).
- M. Creutz, Phys. Rev. D 21, 2308 (1980). N. Cabibbo and E. Marinari, Phys. Lett. B 119, 387 (1982).
- A. D. Kennedy and B. J. Pendleton, Phys. Lett. B 156, 393 (1985).
- F. R. Brown and T. J. Woch, Phys. Rev. Lett. 58, 2394 (1987).
- M. Albanese et al., Phys. Lett. B 192, 163 (1987).
- J.-M. Drouffe and J.-B. Zuber, Phys. Rep. 102, 1 (1983).
- A. M. Ferrenberg and R. H. Swendsen, Phys. Rev. Lett. 63, 1195 (1989).
- M. Fukugita, M. Okawa, A. Ukawa, Phys. Rev. Lett. 63, 1768 (1989).
- Y. Iwasaki et al., Phys. Rev. D 46, 4657 (1992).
- G. P. Lepage and P. B. Mackenzie, Phys. Rev. D 48, 2250 (1993), arXiv:hep-lat/9209022.
- G. Boyd et al., Nucl. Phys. B 469, 419 (1996), arXiv:hep-lat/9602007.
- M. Teper, arXiv:hep-th/9812187 (1998).
- B. Lucini, M. Teper, U. Wenger, JHEP 01 (2004) 061, arXiv:hep-lat/0307017.
- Y. Aoki et al., Nature 443, 675 (2006), arXiv:hep-lat/0611014.
- W. Janke and R. Villanova, Nucl. Phys. B 489, 679 (1997).
- A. Bazavov et al. (HotQCD), Phys. Lett. B 795, 15 (2019), arXiv:1812.08235.
- Y. Aoki et al. (FLAG 2021), Eur. Phys. J. C 82, 869 (2022), arXiv:2111.09849.
- E. W. Kolb and M. S. Turner, The Early Universe (1990).
