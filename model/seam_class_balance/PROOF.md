# The class balance of one seam link

Labels: EXACT (closed form or exact arithmetic, a named check in [results.json](results.json)), MEASURED (Monte
Carlo, with one standard error), ARGUED, IMPORTED, DECLARED. Conventions: weak factor `SU(2)`, coupling `g^2`,
`alpha = g^2/(4 pi)`, `S = 2 pi/(4 alpha) = 2 pi^2/g^2`, unit heat time `tau = g^2/2`, Wilson coupling
`beta = 4/g^2`. Chain coupling `1/alpha_U = 24.3166`, `S_U = 38.1964`. Every result has the form
`<A> = M exp(-c S)`. The working law of the chain is `c = 1` with `b M = 1`.

## 1. The operator and the identity (EXACT)

The link space is `L^2(SU(2))` with electric sectors `j`, dimension `d_j = 2j + 1`, Casimir `C_j = j (j + 1)`.
The centre shift `(A psi)(V) = psi(-V)` is diagonal in the sectors with eigenvalue `(-1)^(2j)`, which is
`(-1)^(6Y)` on the link of the full group. With the sector law `p_j`,

\[ \langle A\rangle = \sum_j p_j\,(-1)^{2j}. \]

One tick is the transfer matrix `T = K^(1/2) exp(-W) K^(1/2) P`. `K` is the product over seams of the one-tick
kernel `k(V' V^dag)` of the electric term, `W` the face term and `P` the Gauss projector. `A` commutes with `K`
and with `P`, so `K^(1/2) A K^(1/2) = A K`, and the kernel of `A K` on the seam is `k(-V' V^dag)`. Hence

\[ \langle A\rangle = \frac{Z[\text{one time-like plaquette of the seam with weight } k(-V_P)]}{Z}. \]

The heat-kernel weight `k = K_tau = sum_j d_j chi_j exp(-tau C_j)` is `exp(-(g^2/2) C)` exactly. The Wilson weight is
`k = exp((beta/2) tr V)`. `A` is read at the midpoint of the electric step (DECLARED). Section 10 gives the
reading at the midpoint of the face step.

## 2. One register in the one-tick thermal law (EXACT)

A sector law `p_j` proportional to `d_j^w exp(-t C_j)` has weight power `w` and heat time `t`. At `t = tau` and the
chain coupling:

| Register | w | `<A>` | c | M |
| --- | --- | --- | --- | --- |
| Cartan circle | | 5.15861e-17 | 1 | 2.000000000000 |
| full link | 2 | 3.88923e-15 | 1 | 150.7857 = 2 (2S - 1) |
| one-sided register | 1 | 0.03378 | no exponential | leading form `tau/8` = 0.03230 |
| labels alone | 0 | 0.16740 | no exponential | `1/(2 sqrt(pi/tau) - 1)` |

The full link in this law is one free link with a pinned plaquette: `<A> = K_tau(-1)/K_tau(1)`. The exponential
needs the two-sided sum over sectors: weight power 2, or the Cartan circle. Weight powers 1 and 0 give power laws
in `alpha` at every heat time. The dual form of the heat kernel used by the Monte Carlo agrees with the character
sum to 5.6e-13 in the logarithm.

## 3. The Cartan circle on a complex (EXACT)

Link angle `eps` of period `4 pi`, plaquette flux `F = d eps`, Villain weight
`sum_n exp(-(F - 4 pi n)^2/(2 g^2))`, twist `F -> F + 2 pi` on one plaquette. Integer shifts of the links extend
the angles to the real line and the Gaussian integral is exact:

\[ Z[v] = \sum_{[n]} e^{-8\pi^2\,(n+v)^{T}(1-\Pi)(n+v)/g^2}, \]

summed over integer plaquette fields modulo coboundaries, with `Pi` the orthogonal projector on the range of the
coboundary, `v = 0` for `Z` and `v = -e_P/2` for the twisted `Z`. At weak coupling the two orientations of the
half turn dominate:

\[ \langle A\rangle = 2\,e^{-S\,(1-\Pi_{PP})}, \qquad c = 1-\Pi_{PP}, \qquad M = 2 . \]

`Pi_PP g^2` is the Gaussian flux variance of the twisted plaquette. The electric variance of the link is
`<E^2> = (1 - Pi_PP)/g^2`.

| Complex | c | count `4/c` | `<A>` at the chain coupling |
| --- | --- | --- | --- |
| one free link, one tick, no Gauss projection | 1 | 4.000 | 5.159e-17 |
| two faces on one seam, outer links pinned (Euclidean) | 1/2 = 0.500000 | 8.000 | 1.016e-08 |
| icosahedral surface as a Euclidean two-complex | 1/20 = 0.050000 | 80.000 | 0.2959 |
| icosahedral surface as space, one tick with the Gauss law | 19/30 = 0.633333 | 6.316 | 6.237e-11 |
| icosahedral surface as space, 8 ticks | 0.398328 | 10.042 | 4.936e-07 |
| icosahedral surface as space, many ticks | 0.398231 | 10.044 | 4.954e-07 |
| face pair as space, many ticks | 0.353553 | 11.314 | 2.73e-06 |
| cubic lattice, three dimensions | 1/3 = 0.333333 | 12.000 | 5.909e-06 |
| cubic lattice, four dimensions | 1/2 = 0.500000 | 8.000 | 1.016e-08 |

- Cubic lattice of side `L` in `d` dimensions: `Pi_PP = (2/d)(1 - L^-d)`, checked in rational arithmetic for
  `(d, L) = (2, 3), (3, 3), (4, 2)`.
- Icosahedral space with many ticks: `c = (1/30) sum sqrt(lam/(lam + 4))` over the nonzero eigenvalues of the face
  Laplacian, `3 - sqrt 5` (3), `2` (5), `3` (4), `5` (4), `3 + sqrt 5` (3). In continuous time the same sum with
  `sqrt(lam)/2` gives `0.5405`.
- A closed two-dimensional complex has classes `n = N e_P` and the sum is a theta function. The pinned face pair
  agrees with direct quadrature to ten digits.
- A free two-dimensional boundary absorbs the twist: `<A> = 1`.

## 4. Flux densities (EXACT)

The density of a plaquette flux at a half turn relative to zero is `2 exp(-S/Pi_PP)` in the Gaussian ensemble.

- Icosahedral surface: 12 corners, 30 seams, 20 faces, rank 19, every seam in two faces with opposite orientation,
  face variance `(19/20) g^2`. One seam at `2 pi` gives `(+2 pi, -2 pi)` on its two faces with Gaussian action
  `4 pi^2/g^2`. The joint density of the pair is the square of the one-face density: `2.661e-33` at the chain coupling.
- Cubic lattices: variance `2 g^2/d`, density exponent `(d/2) S`, count `8/d = 4, 8/3, 2`.

These are statements about the density of the flux. The class balance of the same ensembles has `c = 1 - 2/d`,
count `12` for `d = 3` and `8` for `d = 4`. The two agree only for a free face (density) against a pinned
plaquette (balance), both with count 4.

## 5. Monte Carlo for SU(2): estimator and gates (MEASURED)

The twisted plaquette carries `ln W = (1 - lam) ln k(V) + lam ln k(-V)`, `lam = 0, 1/N, ..., 1`, and
`D = ln k(-V) - ln k(V)`. Each step uses the exact identity
`ln Z_(k+1)/Z_k = ln <exp(D/2N)>_k - ln <exp(-D/2N)>_(k+1)`. The full-step forward and backward ratios and the
Simpson integral of `<D>` are checks. Errors: jackknife over 20 blocks per stage. Updates: multi-hit Metropolis
on unit quaternions, one sweep of the lattice followed by 20 local sweeps around the twisted plaquette.
The same estimator in plain numpy on the pinned face pair gives `-4.2469 (0.0052)` against the exact `-4.2417`.

| Gate | measured `ln <A>` | exact | pull |
| --- | --- | --- | --- |
| heat kernel, icosahedral two-complex, `g^2 = 0.5168` (Migdal) | -0.1847 (0.0127) | -0.1662 | -1.46 |
| heat kernel, icosahedral two-complex, `g^2 = 0.1` (Migdal) | -6.2046 (0.0300) | -6.2458 | 1.37 |
| Wilson, icosahedral two-complex, `beta = 7.7402` (Migdal) | -0.1419 (0.0130) | -0.1289 | -1.00 |
| Wilson, icosahedral two-complex, `beta = 30` (Migdal) | -3.9269 (0.0329) | -3.9747 | 1.46 |
| Cartan circle, `4^3`, `g^2 = 1` | -6.0888 (0.0084) | -6.0922 | 0.41 |
| Cartan circle, `4^4`, `g^2 = 1` | -9.2329 (0.0098) | -9.2150 | -1.83 |
| Cartan circle, icosahedral space, 4 ticks, `g^2 = 1` | -7.2507 (0.0072) | -7.2518 | 0.15 |
| Cartan circle, `4^4`, chain coupling | -18.4656 (0.0130) | -18.4797 | 1.08 |
| Cartan circle, `6^4`, chain coupling | -18.4253 (0.0173) | -18.4198 | -0.32 |
| Cartan circle, `8^3`, chain coupling | -12.1024 (0.0207) | -12.0887 | -0.66 |
| Cartan circle, icosahedral space, 8 ticks, chain coupling | -14.5243 (0.0124) | -14.5215 | -0.22 |

Further gates of the engine: Wilson plaquette 0.6032 (0.0006) at `beta = 2.3` and 0.6540 (0.0004) at `2.5` on `6^4`
against the references 0.6024 and 0.6520 (IMPORTED). The second lies 0.3 per cent, five standard errors, above its
reference on this small lattice. Three gates were run with the production and are not stored in
`measured_points.json`: compact doublet at `beta_H = 0` with `sin^2 theta = 0.7996 (0.0003)` against `4/5`; local
weights against differences of the total weight to `1e-11`; gauge invariance to `1e-12`.

## 6. SU(2): the measured laws (MEASURED)

Couplings `beta = 2.5, 3, 4, 5, 6, 7.7402, 9` (`1/alpha = pi beta`). Fit `ln <A> = ln M0 + p ln S - c S` over
`beta >= 4`, five couplings. `M_U = <A>(S_U) exp(c S_U)` is the multiplicity at the chain coupling, quoted with c at
its prediction.

| Complex | Weight | c (p = 1) | chi2/dof | c (p free) | p | prediction | `ln <A>` at `S_U` | `M_U` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cubic `8^3` | heat kernel | 0.3347 (0.0012) | 5.2/3 | 0.3447 (0.0112) | 1.30 (0.34) | `1 - Pi` = 0.3346 | -7.795 (0.021) | 146.5 (3.1) |
| cubic `8^3` | Wilson | 0.3197 (0.0013) | 6.9/3 | 0.3381 (0.0081) | 1.55 (0.24) | classical minimum 0.3193 | -7.042 (0.021) | 172.9 (3.6) |
| cubic `6^4` | heat kernel | 0.5001 (0.0008) | 1.7/3 | 0.4978 (0.0067) | 0.93 (0.20) | `1 - Pi` = 0.5004 | -13.921 (0.018) | 179.7 (3.3) |
| cubic `6^4` | Wilson | 0.4196 (0.0014) | 18.0/3 | 0.3999 (0.0049) | 0.40 (0.15) | classical minimum 0.4053 | -11.612 (0.013) | 47.9 (0.6) |
| icosahedral space, 8 ticks | heat kernel | 0.3988 (0.0005) | 0.2/3 | 0.3968 (0.0044) | 0.94 (0.13) | `1 - Pi` = 0.3983 | -10.292 (0.012) | 137.4 (1.7) |
| icosahedral space, 8 ticks | Wilson | 0.3693 (0.0005) | 3.2/3 | 0.3655 (0.0045) | 0.89 (0.13) | classical minimum 0.3664 | -9.162 (0.012) | 125.7 (1.5) |

- Heat kernel: the three coefficients equal `1 - Pi_PP` of the Cartan circle on the same complex, within one
  standard error and within 0.2 percent.
- Wilson: the coefficients equal the classical minima of section 8 within 1.4 percent. On the `6^4` lattice the fit
  with `p = 1` is excluded by its chi2 and the free fit gives `p = 0.40 (0.15)`. The classical minimum there is the
  pinned configuration, which keeps the global rotation. On the icosahedral space the fit with `p = 1` lies 0.8
  percent and 5.7 standard errors above the minimum and the free fit agrees with it.
- The half turn breaks the global rotation to a circle. This gives the prefactor one power of `S`.
- No coefficient is 1. The measured values lie between 0.32 and 0.50.

## 7. The sector law of the lattice link (MEASURED)

`p_n = f_n <chi_n(V_P)/k(V_P)>` over the time-like plaquettes of the untwisted ensemble, `f_n` the character
coefficients of `k` (`n exp(-tau (n^2-1)/4)` for the heat kernel, `2 n I_n(beta)/beta` for Wilson). Fit
`ln p_n = const + w ln n - t (n^2 - 1)/4` at the chain coupling:

| Complex | Weight | weight power w | `tau/t` (w free) | `tau/t` (w = 2) | c of section 6 |
| --- | --- | --- | --- | --- | --- |
| cubic `8^3` | heat kernel | 1.999 (0.062) | 0.3039 | 0.3038 (0.0028) | 0.3347 |
| cubic `8^3` | Wilson | 1.950 (0.035) | 0.2997 | 0.2943 (0.0022) | 0.3197 |
| cubic `6^4` | heat kernel | 2.005 (0.002) | 0.4755 | 0.4763 (0.0003) | 0.5001 |
| cubic `6^4` | Wilson | 1.961 (0.012) | 0.4766 | 0.4700 (0.0019) | 0.4196 |
| icosahedral space, 8 ticks | heat kernel | 1.929 (0.053) | 0.3868 | 0.3780 (0.0035) | 0.3988 |
| icosahedral space, 8 ticks | Wilson | 1.957 (0.041) | 0.3780 | 0.3721 (0.0027) | 0.3693 |

- The lattice link carries weight power 2 for both weights.
- Heat kernel: the heat time is `tau/c'` with `c'` five to nine percent below `c` at the chain coupling. On the
  `6^4` lattice `c'` rises with `beta`: 0.4505, 0.4621, 0.4691, 0.4763, 0.4807 for `beta = 4` to `9`. The link is
  colder than the one-tick law by the factor `1/c`.
- Wilson: the bulk of the sector law has nearly the heat-kernel width, while the coefficient of the class balance
  is the classical minimum. The balance is an alternating sum and is set by the tail of the law.

## 8. Wilson weight (EXACT minimum, MEASURED law)

One pinned link: `<A> = exp(-2 beta) = exp(-8/g^2)` exactly, so `c = 4/pi^2 = 0.4053`, `M = 1`, count `pi^2`.
On a lattice the exponent is the minimum of the twisted action, `s_W/g^2`, `c = s_W/(2 pi^2)`:

| Complex | `s_W` | c | flux left on the twisted plaquette, in units of `2 pi` |
| --- | --- | --- | --- |
| cubic `6^4` | 8.0000 | 0.4053 | 0.000 |
| cubic `8^3` | 6.3021 | 0.3193 | 0.609 |
| icosahedral space, 8 ticks | 7.2327 | 0.3664 | 0.490 |

In four dimensions the minimum is the pinned configuration, `s_W = 8`. The measured coefficients of section 6
agree with these minima within the spread of the two fits.

## 9. With the doublet (MEASURED)

Compact gauge-Higgs form: `n = (cos theta, sin theta q)` on the unit four-sphere per site, hop term
`beta_H [1 - n0 n0' - Re q^dag V q']`, on-site term `-h sin^2 theta`, `beta_H = 2 beta_G` (DECLARED). The hop term
has unit kinetic coefficient on `Phi = (2/g) theta`. Lattice `4^4`, `beta_G` at the chain coupling. At `h = 1.5`
both sides are stable for the length of the run.

| Weight | `beta_H/beta_G` | h | start | `sin^2 theta` | side | `ln <A>` |
| --- | --- | --- | --- | --- | --- | --- |
| Wilson | 2 | 0.0 | pole | 0.040 | symmetric | -11.640 (0.011) |
| Wilson | 2 | 0.5 | pole | 0.045 | symmetric | -11.633 (0.012) |
| Wilson | 2 | 1.5 | equator | 0.989 | Higgs | -13.328 (0.006) |
| Wilson | 2 | 1.5 | pole | 0.053 | symmetric | -11.635 (0.012) |
| Wilson | 2 | 3.0 | equator | 0.989 | Higgs | -13.335 (0.005) |
| Wilson | 2 | 4.0 | equator | 0.991 | Higgs | -13.345 (0.005) |
| Wilson | 8 | 4.0 | equator | 0.995 | Higgs | -14.505 (0.002) |
| Wilson | 32 | 4.0 | equator | 0.996 | Higgs | -15.152 (0.001) |
| heat kernel | 2 | 0.0 | pole | 0.040 | symmetric | -13.985 (0.017) |
| heat kernel | 2 | 0.5 | pole | 0.041 | symmetric | -13.983 (0.017) |
| heat kernel | 2 | 1.5 | equator | 0.989 | Higgs | -18.234 (0.013) |
| heat kernel | 2 | 1.5 | pole | 0.045 | symmetric | -14.111 (0.017) |
| heat kernel | 2 | 3.0 | equator | 0.990 | Higgs | -18.262 (0.013) |
| heat kernel | 2 | 4.0 | equator | 0.990 | Higgs | -18.274 (0.013) |
| heat kernel | 8 | 4.0 | equator | 0.997 | Higgs | -24.044 (0.007) |
| heat kernel | 32 | 4.0 | equator | 0.998 | Higgs | -29.476 (0.004) |

- Heat kernel: `ln <A>` changes by -4.230 (0.012) from the symmetric side to the Higgs side, a factor 0.015. Gaussian estimate with the link mass `m^2 = 2 sin^2 theta`: -3.91.
- Wilson: `ln <A>` changes by -1.700 (0.007) from the symmetric side to the Higgs side, a factor 0.18.
- A stiffer doublet pins the link. The heat-kernel value moves toward the one-link law
  `ln(2 (2 S_U - 1)) - S_U = -33.181`, the Wilson value toward `-2 beta = -15.480`.

## 10. Second reading: the midpoint of the face step

With `T = exp(-W/2) K exp(-W/2)` the insertion of `A` between two ticks twists the time-like plaquette above the
seam and half-twists the faces of the seam in that time slice. Gaussian exponent on the Cartan circle (EXACT):

| Complex | c, face step | c, electric step |
| --- | --- | --- |
| cubic `6^3` | 0.7008 | 0.3364 |
| cubic `4^4` | 1.2874 | 0.5020 |
| icosahedral space, 8 ticks | 0.7424 | 0.3983 |

Cartan circle through the Monte Carlo at the chain coupling (MEASURED against the Gaussian form with `M = 2`):

- cubic `6^3`, face step: -26.093 (0.020) against -26.074.
- cubic `4^4`, face step: -48.485 (0.020) against -48.480.
- cubic `4^4`, electric step through the second engine: -18.484 (0.014) against -18.480.
- icosahedral space, 8 ticks, face step: -27.668 (0.014) against -27.662.

SU(2) (MEASURED), fit with `p = 1` over `beta = 4, 6, 7.7402, 9`:

| Complex | Weight | c | Gaussian c | `ln <A>` at `S_U` |
| --- | --- | --- | --- | --- |
| cubic `4^4` | heat kernel | 1.2919 (0.0016) | 1.2874 | -39.922 (0.024) |
| cubic `4^4` | Wilson | 1.0830 (0.0015) |  | -32.385 (0.034) |
| cubic `6^3` | heat kernel | 0.7046 (0.0022) | 0.7008 | -19.399 (0.023) |
| cubic `6^3` | Wilson | 0.5976 (0.0013) |  | -15.325 (0.026) |
| icosahedral space, 8 ticks | heat kernel | 0.7443 (0.0007) | 0.7424 | -21.372 (0.015) |
| icosahedral space, 8 ticks | Wilson | 0.6397 (0.0007) |  | -17.383 (0.016) |

With `n` electric steps per tick a transverse mode of spatial frequency `w` contributes `w/sqrt(w^2/n^2 + 4)` at the
midpoint of the electric step, `w sqrt(w^2/n^2 + 4)/4` at the midpoint of the face step and `w/2` in continuous time
(EXACT, Cartan circle):

| Steps per tick | cubic space, 3 dimensions: electric, face | cubic space, 2 dimensions | icosahedral space |
| --- | --- | --- | --- |
| 1 | 0.5000, 1.2844 | 0.3333, 0.6975 | 0.3982, 0.7422 |
| 2 | 0.6731, 0.9436 | 0.4236, 0.5428 | 0.4892, 0.5982 |
| 3 | 0.7329, 0.8649 | 0.4515, 0.5085 | 0.5155, 0.5670 |
| 4 | 0.7583, 0.8355 | 0.4629, 0.4959 | 0.5259, 0.5556 |
| 8 | 0.7859, 0.8060 | 0.4748, 0.4833 | 0.5367, 0.5443 |
| continuous | 0.7959, 0.7959 | 0.4790, 0.4790 | 0.5405, 0.5405 |

The two readings bracket the continuous-time value. In three space dimensions the face-step coefficient passes
through 1 at 1.67 steps per tick, between 1.284 for one step and 0.944 for two. No whole number of steps gives the
coefficient 1 in either reading on any of the three spaces. The measured coefficient closest to 1 is 1.083 (Wilson,
`4^4`, face step), 8 percent above.

## 11. Which ensemble has the working law

- `c = 1` (EXACT): one free link in its one-tick thermal law, with the heat-kernel weight. `M = 2` on the Cartan
  circle and `M = 150.79` on the full link. `b M = 1` then needs `b = 1/2` or `b = 1/150.79`.
- `c = 1` as a limit (MEASURED): a link held by a stiff doublet.
- Not `c = 1` (EXACT on the Cartan circle, MEASURED for `SU(2)`): every lattice ensemble with the Gauss law and a
  moving link. Icosahedral space 0.398, cubic three dimensions 1/3, cubic four dimensions 1/2 for the heat
  kernel; 0.37, 0.32, 0.40 for Wilson.
- Weight power 1 or 0 at any heat time: no exponential.

## 12. Declared and imported items

- DECLARED: the tick equals the cell (isotropic lattices); `A` read at the midpoint of the electric step; the
  compact doublet with `beta_H = 2 beta_G` and the on-site term `h sin^2 theta`; the fit window `beta >= 4`.
- ARGUED: the power `p = 1` of the prefactor; the Gaussian estimate with a link mass.
- IMPORTED: the reading of the class-3 seam flip as the centre shift of the weak link, and `v = E_cell b <A>`,
  from the electroweak step of the working chain; `E_cell = 9.5599e18 GeV` and `1/alpha_U = 24.3166`, the chain's
  declared forms; Migdal's solution of two-dimensional gauge theory; Wilson plaquette references at `beta = 2.3, 2.5`.
- Open: the sector law of the seam link at the smallest collar that has `c = 1` and `b M = 1`.

## 13. Comparison with the measured scale (NOT BLIND)

Run once, after the fits were fixed. `v = 246.22 GeV` (Particle Data Group 2024, from `G_F`,
https://pdg.lbl.gov/). Two readings: `b M = 1`, `v = E_cell exp(-c S_U)`; and `b = 1`, `v = E_cell <A>(S_U)`.
The measured value needs `c = 1.00004` at `b M = 1`.

| Row | c | v at `b M = 1` (GeV) | v at `b = 1` (GeV) |
| --- | --- | --- | --- |
| one tick, Cartan circle (exact) | 1.0000 | 246.6 | 493.2 |
| one tick, full link, weight power 2 (exact) | 1.0000 | 246.6 | 3.718e+04 |
| one tick, one-sided register, weight power 1 (exact) |  |  | 3.23e+17 |
| one tick, labels, weight power 0 (exact) |  |  | 1.6e+18 |
| Cartan circle, cub4_infinite (exact) | 0.5000 | 4.855e+10 | 9.71e+10 |
| Cartan circle, cub3_infinite (exact) | 0.3333 | 2.825e+13 | 5.649e+13 |
| Cartan circle, icoT_infinite (exact) | 0.3982 | 2.368e+12 | 4.736e+12 |
| Cartan circle, icoT_NT1_Gauss (exact) | 0.6333 | 2.981e+08 | 5.963e+08 |
| Cartan circle, pair2_pinned_boundary (exact) | 0.5000 | 4.855e+10 | 9.71e+10 |
| Cartan circle, ico2 (exact) | 0.0500 | 1.416e+18 | 2.829e+18 |
| SU(2) cub4_L6_hk, electric step | 0.5001 | 4.831e+10 | 8.599e+12 |
| SU(2) cub4_L6_w, electric step | 0.4196 | 1.046e+12 | 8.658e+13 |
| SU(2) cub3_L8_hk, electric step | 0.3347 | 2.68e+13 | 3.938e+15 |
| SU(2) cub3_L8_w, electric step | 0.3197 | 4.75e+13 | 8.359e+15 |
| SU(2) icoT_NT8_hk, electric step | 0.3988 | 2.318e+12 | 3.242e+14 |
| SU(2) icoT_NT8_w, electric step | 0.3693 | 7.152e+12 | 1.004e+15 |
| SU(2) cub4_L4_hk, face step | 1.2919 | 0.003552 | 43.89 |
| SU(2) cub4_L4_w, face step | 1.0830 | 10.37 | 8.24e+04 |
| SU(2) cub3_L6_hk, face step | 0.7046 | 1.963e+07 | 3.593e+10 |
| SU(2) cub3_L6_w, face step | 0.5976 | 1.168e+09 | 2.113e+12 |
| SU(2) icoT_NT8_hk, face step | 0.7443 | 4.31e+06 | 4.995e+09 |
| SU(2) icoT_NT8_w, face step | 0.6397 | 2.34e+08 | 2.698e+11 |

- The one-tick law of a single link gives 246.6 GeV at `b M = 1`. It is the working law itself.
- Every lattice ensemble in the electric-step reading gives between `4.8e10` and `4.8e13 GeV` at `b M = 1` and
  between `8.6e12` and `8.4e15 GeV` at `b = 1`.
- One row of the second reading lies within a factor of ten of the measured value: heat kernel, `4^4`, face step,
  `b = 1`, 43.9 GeV, 0.18 of the measured value. Its coefficient is 1.29 and its multiplicity `1.2e4`. It is one of
  22 rows compared and counts as a look.
