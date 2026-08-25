# Joint SPARC rotation-curve and BTFR profile likelihood: declared conventions

Producer: `joint_rar_likelihood.py`. Receipt:
`runtime/joint_likelihood_receipt.json`. Independent replay:
`verify_joint_likelihood_independent.py` (byte-identical recomputation from
the snapshot). Tests: `test_joint_rar_likelihood.py`. This document lists
every choice the numbers depend on, with its status, and closes with the
exact list of what this likelihood does not do.

Scope of the word "joint" in the title: the two channels share one set of
profiled nuisance values and are resampled jointly in the paired bootstrap;
the likelihood object itself is the per-point rotation-curve likelihood,
and no combined two-channel likelihood over both data channels is formed.
The BTFR channel enters as a matched estimator compared through the paired
bootstrap.

Status vocabulary: **declared** (a fixed analysis choice, made before
reading the results and carried in the receipt), **catalogue-supplied**
(taken per galaxy or per point from the committed SPARC snapshot),
**open** (a needed construction that does not exist here; the receipt says
so).

## Data

| Choice | Value | Status |
| --- | --- | --- |
| Snapshot | `../data/table1.dat`, `../data/table2.dat`, CDS J/AJ/152/157, pinned by sha256 in the receipt | catalogue-supplied |
| Point cuts | quality <= 2, catalogue inclination >= 30 deg, v_obs > 0, e_Vobs > 0, e_Vobs/v_obs <= 0.10, g_bar_cat > 0 | declared (the committed standard cuts) |
| Retained sample | 2700 points in 149 galaxies; 0 points excluded for nonpositive e_Vobs | consequence of the cuts |
| Deep subsets | g_bar_cat < f * 1.2e-10 m/s^2, f in {0.3, 0.1}, fixed before fitting at catalogue nuisance values; nuisance shifts do not re-select points | declared; the reference scale is a comparison-only import, never an OPH output |
| Outermost point | largest radius among the retained points of each galaxy | declared |

## Model

| Choice | Value | Status |
| --- | --- | --- |
| Acceleration law | g_model = g_bar + sqrt(g_bar a0), the deep-regime candidate profile of the dark-matter paper | declared candidate; the source does not fix a0 |
| Velocity form | v_model^2(r) = d [v_bar_cat^2(r) + sqrt(v_bar_cat^2(r) a0 r_cat)] at r = d r_cat | derived from the law plus the distance scaling below |
| Baryonic speed | v_bar_cat^2 = v_gas\|v_gas\| + Upsilon_d v_disk\|v_disk\| + 0.7 v_bul\|v_bul\| (signed gas, catalogue components at unit mass-to-light) | declared, the committed convention |
| Distance scaling | r = d r_cat (fixed angular sizes); M proportional to D^2 (photometry), so v_bar^2 and the anomalous term both scale as d; table-1 masses scale as d^2 | re-derived in the module docstring |
| Inclination scaling | v_obs(i) = v_obs_cat sin(i_cat)/sin(i), e_Vobs scaled by the same factor; differential form dv/v = -cot(i) di | re-derived in the module docstring |
| Full-sample row | the same law applied to all retained points is the additive all-gradient extension; its misfit bears on that extension submodel only | declared scope |

## Error model and covariance

| Choice | Value | Status |
| --- | --- | --- |
| Per-point error | Gaussian in velocity, sigma = e_Vobs sin(i_cat)/sin(i); the 2 n ln s(i) variance-normalization term is kept in the objective | declared; e_Vobs is catalogue-supplied |
| Intra-galaxy covariance | equal-correlation block R = (1-rho) I + rho J per galaxy, closed-form quadratic form | declared family |
| rho treatment | scanned over the fixed grid {0.0, 0.2, 0.4, 0.6}; the a0 interval is reported at every rho; no single rho is selected | declared scan; the calibration of rho is **open** |
| Not propagated | luminosity and HI-mass measurement errors, stellar-population systematics beyond the mass-to-light prior, baryonic decomposition errors, inter-galaxy correlations | open |

## Nuisance parameters (profiled on declared grids)

| Nuisance | Prior | Grid | Status |
| --- | --- | --- | --- |
| Disk mass-to-light Upsilon_d | log-normal, center 0.5, width 0.1 dex | 9 points over +-2.5 prior sigma | declared |
| Distance factor d | Gaussian at 1, sigma = e_Dist/Dist per galaxy | 7 points over +-2.5 sigma, clipped below at 0.05 | catalogue-supplied errors |
| Inclination i | Gaussian at i_cat, sigma = e_i per galaxy | 7 points over +-2.5 sigma, clipped to [10, 89.9] deg | catalogue-supplied errors |
| Bulge mass-to-light | fixed at 0.7 | none | declared fixed fallback; the snapshot carries no per-galaxy bulge datum |
| Absent error datum | nuisance fixed at the catalogue value, counted in the receipt | none | declared fallback; count is 0 on this snapshot for both distance and inclination |

## Fitting, intervals, goodness of fit

| Choice | Value | Status |
| --- | --- | --- |
| a0 grid | log10 a0 in [-11.0, -9.5], 301 points, step 0.005 dex | declared |
| Profiling | per galaxy and per a0, minimize penalized -2 ln L over the nuisance grid; curves are summed over galaxies | declared |
| ML point | grid argmin (inside its own intervals by construction); the three-point parabola vertex is a labeled refinement row | declared |
| Intervals | delta(-2 ln L) thresholds 1.0 (68.3 percent) and 3.84 (95 percent), one-parameter asymptotic chi-square, linear crossing on the grid | declared; approximate under grid profiling with penalties; intervals narrower than one grid step carry a resolution flag |
| Goodness of fit | data quadratic form at the ML grid point, reduced by two dof conventions (with and without counting the profiled nuisances) that bound the effective count | declared |
| Misfit statement | likelihood-only: reduced chi-square, boundary-pinned endpoints, nuisance grid-edge stacking; no posterior quantity | declared |
| Seed | 20260825, fresh generator per bootstrap block | declared |
| Bootstrap | 1000 replicates, common-set galaxies resampled with replacement; channel A re-minimized per replicate; channel-B nuisances held at point estimates (declared simplification); zero tallies as counts with plus-one fractions | declared |
| Channel-B estimator | a0_B = (v_A^2)^2 / (G M_b d^2), v_A^2 = (v_obs(r_out) s(i))^2 - d v_bar_cat^2(r_out, U), at the profiled (U, d, i) of the subset ML for the same rho | declared; matched to the channel-A observable family |
| Verdict rule | direction-neutral tension/consistent labels on the 95 percent paired interval; consistency is shared with the standard null and is not evidence for OPH; a tension counts only against this declared submodel | declared |

## What the receipt shows (headline rows, committed snapshot)

Reduced chi-square sits well above one on every row (3.6 to 16.4 across
subsets and rho values under the point-count dof convention; 5.7 to 19.6
counting the profiled nuisances as free), so the declared error family understates the
observed scatter and the quoted intervals are conditional on that family;
many 68 percent intervals carry the grid-resolution flag for the same
reason. At f = 0.1 the a0 maximum-likelihood point moves from 7.94e-11
(rho = 0) to 5.82e-11 m/s^2 (rho = 0.6): the correlation choice, which is
open, dominates the systematic budget. The paired channel comparison is
CONSISTENT_INTERVAL_CONTAINS_ZERO at every (f, rho) combination, for
example [-0.027, +0.178] dex at f = 0.1, rho = 0. Consistency selects
nothing: every model reproducing both relations, including the standard
null, shares it.

## What this likelihood does NOT do

1. It does not derive a source value of a0. The OPH source does not fix
   a0; a0 enters as a fitted comparison parameter, and the source-value
   question is open. Order item 4's "derive or falsify a source value"
   therefore resolves as: no source value exists to test, and this receipt
   supplies the measurement-side interval that any future source derivation
   must meet.
2. It does not calibrate the intra-galaxy covariance; the equal-correlation
   family is a declared stand-in and the rho scan bounds the sensitivity.
3. It does not model intrinsic scatter; the reduced chi-square rows above
   one say the declared family understates the observed scatter.
4. It does not treat relativistic lensing, clusters, the CMB, structure
   growth, or any cosmological observable.
5. It does not resolve the bulge mass-to-light convention beyond the
   declared fixed 0.7, and it does not propagate luminosity or HI-mass
   measurement errors.
6. It does not arm, freeze, or discharge any prediction, and it is not a
   preregistered contract: every number is a labeled postdiction on the
   seen committed snapshot.
7. It does not authenticate the snapshot; the sha256 pins bytes, never
   provenance or custody.
