# Electroweak crossover of the declared coupling chain

Labels: EXACT, MEASURED (numerical, with error), ARGUED, IMPORTED (with citation), DECLARED. Every number below is
stored in [results.json](results.json) or [production.json](production.json). The 47 named checks in
`results.json` carry the section number as prefix.

## 1. The three-dimensional theory

IMPORTED: Kajantie, Laine, Rummukainen, Shaposhnikov, Nucl. Phys. B458 (1996) 90 and B466 (1996) 189.

At temperatures near 100 GeV the static modes of the electroweak theory at distances above `1/(gT)` are described
by the three-dimensional SU(2) gauge theory with one fundamental Higgs doublet,

\[
L = \tfrac14 F^a_{ij}F^a_{ij} + (D_i\phi)^\dagger D_i\phi + m_3^2\,\phi^\dagger\phi + \lambda_3(\phi^\dagger\phi)^2 .
\]

The coupling `g3^2` sets the scale. The dynamics depends on `x = lambda3/g3^2` and `y = m3^2(g3^2)/g3^4`. The
temperature enters almost only through `y`. Hypercharge adds a U(1) gauge field with `z = g3'^2/g3^2`. This package
simulates `z = 0` and keeps the hypercharge terms of the map of section 3.

## 2. Lattice action and its relation to the continuum

IMPORTED: hep-lat/9510020, section 2; Laine, Nucl. Phys. B451 (1995) 484; Laine and Rajantie, Nucl. Phys. B513 (1998) 471.

\[
S = \beta_G\sum_{x,\,i<j}\bigl(1-\tfrac12{\rm Tr}P_{ij}\bigr)
 - \beta_H\sum_{x,i}\tfrac12{\rm Tr}\,\Phi^\dagger(x)U_i(x)\Phi(x+i)
 + \sum_x R^2 + \beta_R\sum_x (R^2-1)^2,\qquad R^2=\tfrac12{\rm Tr}\,\Phi^\dagger\Phi .
\]

\[
\beta_G=\frac{4}{g_3^2a},\qquad x=\frac{\beta_R\beta_G}{\beta_H^2},\qquad
y=\frac{\beta_G^2}{8}\Bigl(\frac1{\beta_H}-3-\frac{2x\beta_H}{\beta_G}\Bigr)+\frac{3\Sigma\beta_G}{32\pi}(1+4x)
+\frac{1}{16\pi^2}\Bigl[\bigl(\tfrac{51}{16}+9x-12x^2\bigr)\bigl(\ln\tfrac{3\beta_G}{2}+\zeta\bigr)+4.9941+5.2153\,x\Bigr],
\]

with `Sigma = 3.175911535625` and `zeta = 0.08849`. The theory is super-renormalizable, the divergences end at two
loops, and both are included. The relations are exact for `a -> 0`. At finite spacing the error is of order
`1/beta_G`.

- **2.1** EXACT. The substitution `phi^dag phi = (beta_H/2a) R^2` in the naive lattice action gives the hopping
  coefficient `beta_H`, the quartic `beta_R = x beta_H^2/beta_G` and the bare mass term of the relation for `y`.
- **2.2** The relation reproduces the six critical temperatures that hep-lat/9510020, table 4, quotes for its
  critical couplings, within 0.02 GeV. Checks `2_relations.*`.
- **2.3** Condensate and susceptibility. `<phi^dag phi>/g3^2 = (beta_G beta_H/8)(<R^2> - Sigma/(pi beta_H)) -
  (3/16 pi^2)(ln(3 beta_G/2) + 0.66796)` in the MS-bar scheme at the scale `g3^2`. `chi_3 = N^3 (beta_H^2/beta_G)
  Var(R2bar)` for the volume average `R2bar` on `N^3` sites. The susceptibility of D'Onofrio and Rummukainen is
  `chi = V T <(phi^dag phi - <phi^dag phi>)^2> = chi_3 T/g3^2`.

## 3. The map from four-dimensional parameters and temperature

IMPORTED: hep-ph/9508379, sections 4.2, 4.3 and 5. The formulas are transcribed in [reduction.py](reduction.py).

Input set: `G_mu`, `m_W`, `m_Z`, `m_H`, `m_t` (pole masses) and `alpha_s`. The vacuum relations of section 5 give
the MS-bar parameters at one loop. Section 4.2 integrates out the modes with non-zero Matsubara frequency and
section 4.3 the temporal gauge fields. The result is `g3^2`, `x` and `y` to relative order `g^4` in the counting
`lambda ~ g^2`, `g_Y ~ g`, `g'^2 ~ g^3`, with the two-loop logarithm of the mass parameter. The cited papers
estimate the accuracy of the three-dimensional description at about one percent.

- **3.1** Form. Each parameter is the tree value times one plus the vacuum and the thermal one-loop terms. In this
  form the renormalization scale cancels exactly. Check `3_map.scale_independent`.
- **3.2** EXACT. The scale dependence of the transcribed vacuum relations equals the one-loop running of `g^2`,
  `nu^2`, `lambda` and `g_Y^2` printed in the source, and the QCD part of the top relation at `mu = m_t` is
  `-(8/3) alpha_s/pi`. Checks `3_map.running_*`, `3_map.top_qcd_term` and `test_reduction.py`.
- **3.3** Against the literature. For Standard Model inputs (`m_H = 125.1`, `m_t = 173.3`, `m_W = 80.385`,
  `m_Z = 91.1876` GeV, `alpha_s = 0.118`) the map gives `x = 0.2913`, `y = -0.0644`,
  `g3^2/T = 0.39246` at 160 GeV and `y = 0` at 162.54 GeV. The published curves of
  D'Onofrio and Rummukainen, Phys. Rev. D93 (2016) 025003, figure 1, read from the figure source: 0.2899, -0.0536,
  0.39282 and 162.1 GeV. That paper does not list its inputs. Tree-level reduction gives `y = 0` at
  139.9 GeV.

## 4. Code and gates

MEASURED. [lattice.py](lattice.py) is a numpy implementation with checkerboard updates. [kernels.py](kernels.py)
holds the compiled kernels of the production runs. Updates: Kennedy-Pendleton heat bath and overrelaxation for the
links; for the Higgs field a Gaussian heat bath with an accept step for the quartic term, the overrelaxation of
hep-lat/9510020, section 3, and its global radial update.

| Gate (production kernels) | Measured | Reference |
| --- | --- | --- |
| gauge invariance of H, sum R^2, sum R^4, plaquette (largest relative change) | 1.4e-16 | 0 |
| free field, unit links, beta_H = 0.30: <R^2> | 2.46799(85) | 2.46774 |
| free field: hopping term per site | 1.5604(16) | 1.55912 |
| global radial update: acceptance | 0.73 | between 0.2 and 0.95 |
| plaquette at zero hopping, beta_G = 0.5, against u + 2 u^5 | 0.12389(13) | 0.123776 |
| plaquette at zero hopping, beta_G = 1, against u + 2 u^5 | 0.24190(18) | 0.241793 |
| plaquette at zero hopping, beta_G = 8 (reference: leading weak coupling 1 - 1/beta_G) | 0.870699(43) | 0.875 |
| Schwinger-Dyson identity, x = 0.036 | 2.00064(50) | 2 |
| Schwinger-Dyson identity, x = 0.29, heat bath alone | 1.99981(55) | 2 |
| Schwinger-Dyson identity, x = 0.29 | 2.00021(70) | 2 |
| <R^2>, heat bath alone minus heat bath with overrelaxation, x = 0.036 | difference 0.01(12) | 0 |
| <R^2>, heat bath alone minus heat bath with overrelaxation, x = 0.29 | difference 0.006(12) | 0 |

The same gates run on small lattices in [test_lattice.py](test_lattice.py) with the numpy implementation. The two
implementations agree on the plaquette at zero hopping, `beta_G = 8`, `8^3`: 0.87068(16) and 0.870699(43).

## 5. Validation against published lattice results

- **5.1 Condensate in the broken phase.** MEASURED against IMPORTED. Point `m_H* = 80 GeV`, `T* = 145 GeV` of
  hep-lat/9510020 (`x = 0.11885`, `y = -0.23614`), `beta_G = 6, 8, 12, 16, 20` at fixed physical volume. Linear
  extrapolation in `1/beta_G` from `beta_G >= 8`: `<phi^dag phi(T*)>/T* = 0.7072(22)`.
  Published continuum value: 0.7061(7) (ibid., table 1).
- **5.2 Susceptibility for Standard Model inputs.** MEASURED against IMPORTED, NOT BLIND. The pipeline of section 7
  with Standard Model inputs gives the maximum at 159.042(42) GeV on `24^3` at `beta_G = 8`.
  The published maxima are 158.73 GeV at `beta_G = 6` and 159.01 GeV at `beta_G = 9` (arXiv:1508.07161, figure 5),
  158.94 GeV interpolated to `beta_G = 8`. Curve against the published points at `beta_G = 9` with the U(1) field
  (figure 6, square symbols, read from the figure source):

| T (GeV) | This pipeline, `beta_G = 8`, `24^3`, `z = 0` | Published, `beta_G = 9`, `z = 0.31` |
| --- | --- | --- |
| 156 | 6.152(43) | 6.22 |
| 157 | 6.680(49) | 6.84 |
| 158 | 7.868(57) | 8.22 |
| 158.5 | 8.92(12) | 9.24 |
| 158.8 | 9.54(15) | 9.94 |
| 159 | 9.73(14) | 10.02 |
| 159.2 | 9.64(12) | 9.96 |
| 159.5 | 8.85(10) | 8.75 |
| 160 | 6.59(10) | 6.06 |
| 160.5 | 4.556(66) | 4.30 |
| 161 | 3.293(40) | 3.12 |

- **5.3 Critical coupling at x = 0.0644, beta_G = 8.** Susceptibility maximum on 24^3 at beta_H = 0.347982, on 16^3 at 0.347996. IMPORTED infinite-volume value: 0.3479735(6) (hep-lat/9510020, table 4). MEASURED against IMPORTED.

## 6. Phase structure

MEASURED, `beta_G = 8`. The endpoint of the first-order line is IMPORTED: `x_c = 0.0983(15)` in the continuum and
0.1133(25) at `beta_G = 5` (Rummukainen, Tsypin, Kajantie, Laine, Shaposhnikov, Nucl. Phys. B532 (1998) 283).

- **6.1 First-order point `x = 0.036`.**

| N | y at the maximum | chi_3 maximum | histogram valley / peak at equal heights |
| --- | --- | --- | --- |
| 8 | 0.01178(22) | 62.09(49) | 0.328 |
| 12 | 0.02283(18) | 149.3(15) | 0.0743 |
| 16 | 0.02406(40) | 341.9(80) | 0.00596 |

  The maximum grows with the exponent 2.163(32) between 8^3 and 12^3 and 2.881(88) between 12^3 and 16^3; volume scaling is 3. The two-peak histogram of the volume average of
  `phi^dag phi` deepens with the volume. Hysteresis with 300 sweeps per step of 0.004 in `y`: on 16^3 the downward sweep enters the broken phase at y = 0.012 and the upward sweep leaves it at y = 0.032; on 24^3 the downward sweep enters the broken phase at y = -0.004 and the upward sweep leaves it at y = 0.032.
- **6.2 Crossover point `x = 0.13`.**

| N | y at the maximum | chi_3 maximum |
| --- | --- | --- |
| 16 | -0.04516(34) | 13.23(15) |
| 24 | -0.04090(40) | 20.17(33) |
| 32 | -0.03904(44) | 29.0(13) |

  Growth exponent: 1.040(49) between 16^3 and 24^3; 1.26(16) between 24^3 and 32^3.
- **6.3 Growth of the maximum against x.** Exponent `d ln chi_max / d ln N` between `16^3` and `24^3`. The value 3
  is first-order scaling and 1.96 is the Ising value at the endpoint. A crossover gives 0 at large volume. At
  `x = 0.036` the exponent between `12^3` and `16^3` is given in 6.1.

| x | chi_3 max, 16^3 | chi_3 max, 24^3 | exponent |
| --- | --- | --- | --- |
| 0.0644 | 64.79(77) | 203.2(34) | 2.819(50) |
| 0.0897 | 29.50(27) | 67.1(12) | 2.028(49) |
| 0.1150 | 17.01(27) | 29.6(11) | 1.367(96) |
| 0.1300 | 13.23(15) | 20.17(33) | 1.040(49) |
| 0.2960 | 3.407(23) | 3.736(52) | 0.227(38) |

## 7. The chain's point

Inputs DECLARED (the chain's declared forms): `v = 246.60 GeV`, `1/alpha_2(m_Z) = 29.605`,
`sin^2 theta_W = 0.23074`, `alpha_s = 0.11834`, tree masses `m_W = 80.330`, `m_Z = 91.588 GeV`, and the declared
midpoint `m_H = 125.84`, `m_t = 171.99 GeV`. They enter the map as its physical input set with
`G_mu = 1/(sqrt2 v^2)`.

- **7.1** The map gives `x = 0.2961`, `g3^2/T = 0.39083` at the transition and
  `y = 0` at 163.98 GeV. The point lies at three times the endpoint value of `x`: the transition is a
  crossover. Check `5_rows.declared_is_crossover`.
- **7.2** MEASURED. Maximum of the susceptibility at `x = 0.296`, reweighted along the chain's curve `(x(T), y(T))`:

| beta_G | N | size in 1/g3^2 | runs | y at the maximum | chi_3 maximum | T at the maximum (GeV) |
| --- | --- | --- | --- | --- | --- | --- |
| 6 | 12 | 8 | 7 | -0.10293(41) | 3.578(24) | 159.945(17) |
| 6 | 18 | 12 | 7 | -0.09700(89) | 3.841(32) | 160.168(33) |
| 8 | 16 | 8 | 7 | -0.09726(68) | 3.407(23) | 160.161(24) |
| 8 | 24 | 12 | 16 | -0.0908(11) | 3.736(52) | 160.406(43) |
| 8 | 32 | 16 | 5 | -0.0900(27) | 3.833(82) | 160.44(10) |
| 12 | 24 | 8 | 6 | -0.0906(11) | 3.320(41) | 160.414(44) |

- **7.3** Continuum. Fit `T_max = T_c + C/beta_G + D(size)` over these sets, the largest size as reference:
  `T_c = 161.12(12) GeV`, slope C = -5.50(45) GeV, `chi^2 = 0.31` for 2 degrees
  of freedom; volume offsets: size 12: -0.04(10) GeV; size 8: -0.27(10) GeV. The same fit to the position in `y` gives
  `y_c = -0.0719(31)`.
- **7.4** Result.

\[
T_c = 161.13\ \text{GeV},\qquad
\sigma_{T,\mathrm{linear}}\in[0.11810,0.13680]\ \text{GeV}.
\]

  This statistical range assumes independent run sets and permits every covariance of the shared
  reference peaks compatible with their retained marginal errors. It is conditional on those variance
  estimates and linearization. The separate reduction allowance is one percent, `1.6113 GeV`; it has
  no calibrated confidence coverage. The covariance calculation is given below.

  The susceptibility in the normalization of 2.3 peaks at 9.54(13) on `24^3`. It exceeds 80 percent of
  its maximum over 1.9 GeV. It
  falls to half its maximum 1.42 GeV above the peak. Below the peak it levels near 5,
  about half the maximum.
- **7.5** Variants of the map at fixed `y_c` (GeV): strong coupling run to `m_W`: 161.35;
  product form at `mu = pi T`: 161.82; hypercharge dropped from the map:
  162.86; tree-level reduction: 139.2.
- **7.6** NOT BLIND. Standard Model lattice result: `159.5 +- 1.5 GeV` (D'Onofrio and Rummukainen 2016, abstract).
  Difference: +1.63 GeV. No combined significance is assigned: the lattice standard deviation is
  unidentified and the reduction allowance is not a calibrated standard deviation. The same pipeline with
  Standard Model inputs gives 159.76 GeV; the published continuum position of the maximum is 159.58 GeV.

## 8. The source as built with one doublet

Inputs of the row with one Higgs doublet: `v = 222.9 GeV`, `alpha_s = 0.1123`, `sin^2 theta_W = 0.2056`,
`1/alpha_em = 142.25`, hence `g^2 = 4 pi alpha_em/sin^2 = 0.42967` and `g'^2 = 0.11120`. Its Higgs
and top masses are not derived. DECLARED assumptions:

| Assumption | m_H, m_t (GeV) | x at the transition | T_c (GeV) |
| --- | --- | --- | --- |
| A: quartic and top couplings of the declared row | 113.75, 155.46 | 0.2930 | 145.11; linearized SD 0.10755 to 0.12407 |
| B: masses of the declared row | 125.84, 171.99 | 0.3600 | 149.33; linearized SD 0.10421 |

Both are crossovers. Row A is reweighted from the runs at `x = 0.296`. Row B uses runs at `x = 0.36` on `24^3`,
`beta_G = 8`, with the continuum shift of 7.3.
The statistical entries assume independent run sets. Their separate one-percent reduction allowances
are `1.4511 GeV` and `1.4933 GeV`. Row B's temperature peak and the reference peaks belong to different
run sets, so that assumption fixes its linearized variance from the retained marginal errors.

## 9. The endpoint as a Higgs mass

With the other inputs of the declared chain fixed, the map reaches `x_c = 0.0983(15)` on the transition line at
`m_H = 72.5 +- 0.7 GeV`. The error is that of `x_c`. The
reduction adds about 1 GeV. Below this mass the same chain gives a first-order transition. The same procedure with
Standard Model inputs gives 72.3 GeV; the published value is 72(2) GeV.
IMPORTED endpoint, map of section 3.

## 10. Cosmic time

EXACT relation with DECLARED inputs. Radiation era, `t = 1/(2H)`, `H^2 = (8 pi/3) rho/M_Pl^2`,
`rho = (pi^2/30) g* T^4`, `g* = 106.75`, `M_Pl = 1.220890e19 GeV`. At `T_c` of the declared chain:
`t = 9.023e-12 s`, that is 1.310e+32 ticks of `hbar/E_cell` with
`E_cell = 9.5599e18 GeV`.

## 11. What is not established

- The chain's couplings are declared forms. This package does not derive them.
- The Higgs mass term is an input through `m_H`. The package says nothing about its origin.
- The reduction is perturbative. Its one-percent allowance is a declared estimate, without a remainder or coverage theorem.
- The three-dimensional U(1) field is not simulated. The published simulation with that field has its maximum at y = -0.082 at `beta_G = 8`; this simulation has -0.090. The difference corresponds to +0.3 GeV. It is not applied and is smaller than the reduction allowance; that comparison does not bound the omitted field's effect.
- The endpoint value `x_c` is imported. The runs here show the first-order signal at `x = 0.036` and its weakening
  with `x` at one lattice spacing; they do not determine `x_c`.
- No real-time quantity is computed: no sphaleron rate, no nucleation.

## Matched volume and cutoff sensitivity

A cutoff comparison holds the dimensionless physical side `L g3^2 = 4L/beta_G` fixed. The sets `(beta_G,L)`
`(6,12)`, `(8,16)` and `(12,24)` all have side eight. The package fits their peak coordinates with a constant
plus `1/beta_G`, a constant plus `1/beta_G^2`, and both powers. It also repeats the linear fit omitting each
spacing. These six declared forms are all retained. No form is selected using agreement with a reference
temperature. Statistical fit covariance treats the separately simulated sets as independent.

For the declared-chain temperature at side eight, the linear intercept is `160.85458 +- 0.07203 GeV`.
The largest central-value displacement among these forms is `0.34896 GeV`. The fit with both powers and all
three leave-one-out fits have zero residual degrees of freedom. Their exact interpolation supplies no
goodness-of-fit test. This spread is finite design sensitivity, not a bound on omitted cutoff terms or an
infinite-volume extrapolation. It is distinct from the joint size/cutoff fit whose reference side is sixteen.

At `beta_G=8`, the adjacent-volume effective exponent is
`log(chi_max(L2)/chi_max(L1))/log(L2/L1)`. It is `0.22693 +- 0.03804` for sides 16 to 24 and
`0.08886 +- 0.08895` for sides 24 to 32 at the declared chain's point. The first-order control at `x=0.036`
gives `2.16336 +- 0.03170` and `2.88122 +- 0.08832` for sides 8 to 12 and 12 to 16. These are measured
finite-volume responses, not fitted asymptotic critical exponents. The peak temperature shifts by
`0.24438 GeV` from side 16 to 24 and `0.03071 GeV` from side 24 to 32 at the common cutoff.

## Covariance of the shared reference peaks

The retained producer `inputs/EW/code/ew_oph.py.raw` evaluates the peak in `y` and each temperature peak
on the same delete-one-block samples. Its `ew_wham.py.raw` helper returns the full estimate and each
coordinate's jackknife error, discarding the paired deletion vectors. The eleven retained numerical EW
inputs contain JSON summaries and curves; they contain no raw simulation series, paired peak replicas or
cross-covariance matrix. `comparison_budget.py` verifies their manifest and records that input boundary.

Let `Y_i` denote the six peak estimates at `x=0.296`, and let `r` identify the `(beta_G,L)=(8,24)` reference.
The fixed weighted least-squares design gives an intercept `Y_c=sum_i w_i Y_i`. The weights are computed
as exact rationals from the serialized design and marginal errors; they preserve a constant exactly and
annihilate the cutoff and volume-offset columns. The temperature computation is the implicit equation

```
y(T_c) = y(T_lat) + sum_i w_i Y_i - Y_r.
```

With the reduction parameters held fixed, its first variation is

```
delta T_c = [y'(T_lat) delta T_lat + sum_i (w_i-delta_ir) delta Y_i] / y'(T_c).
```

In particular, `T_lat`, `Y_r` and the fitted intercept are not independent measurements. If the control
`y'(T_lat) delta T_lat = delta Y_r` holds exactly, the two reference terms cancel before propagation.
The physical temperature and `y` scans use different susceptibility normalizations and peak fits, so
this control is not imposed on the retained production. The receipt retains `0.125474 GeV` as the
uncorrected independence control, without assigning it to the unidentified statistical error.

For coefficients `a_i` and supplied marginal standard deviations `s_i`, any positive-semidefinite
covariance matrix with diagonal `s_i^2` is a Gram matrix of vectors of lengths `s_i`. The variance of
`sum_i a_i X_i` is the squared length of their weighted sum. Writing `l_i=abs(a_i)s_i`, the sharp bounds are

```
max(0, 2 max_i l_i - sum_i l_i)^2 <= variance <= (sum_i l_i)^2.
```

The upper bound aligns every weighted vector. The lower bound follows from the reverse triangle
inequality and is attained by opposing the longest vector to the others, or by a closed planar polygon
when no length exceeds the sum of the others. These Gram constructions prove attainability. If distinct
run sets are independent, apply this argument separately to each set and add the variances. Only the
two peaks from the reference set then share an unspecified covariance. This is a statistical sampling
assumption, not a consequence of the marginal error table.

For the declared row the resulting linearized standard deviation lies in `[0.1180963,0.1367957] GeV`.
Allowing unspecified cross-set covariance widens the sharp range to `[0,0.2763059] GeV`. The corresponding
ranges for Standard Model inputs are `[0.1169101,0.1352495]` and `[0,0.2732030] GeV`.
The coefficient and variance arithmetic uses exact rationals; square-root endpoints have integer-square
enclosures. The reduction derivatives use a five-point numerical stencil with a reported step-halving
check, not a rigorous derivative remainder bound. The variance estimates themselves have sampling error.
These results therefore identify covariance-compatible delta-method ranges, without nonlinear confidence
coverage, a distributional law or a rigorous total uncertainty for the physical comparison.

`comparison_budget.json` retains the exact influence weights, both covariance domains, numerical derivative
diagnostics, input hashes and cancellation controls. The one-percent reduction estimate, `0.3489639 GeV`
finite-design displacement, finite-volume limit and omitted hypercharge simulation keep separate scopes.
The null primary error and pull fields prevent an unidentified standard deviation or a declared allowance
from becoming a combined significance through a scalar fallback.
