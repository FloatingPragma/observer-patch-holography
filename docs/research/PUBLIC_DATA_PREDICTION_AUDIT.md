# OPH predictions and public-data comparisons

The target is a prediction for our observed Universe that adds information
beyond the Standard Model and General Relativity, followed by an executable
comparison with existing observations. Building a new laboratory realization
of the finite source is not a prerequisite. The observable dictionary, physical
premises and uncertainty model must nevertheless be explicit.

Coverage consists of the selected observations listed below. The target is
one fixed construction with shared parameters and declared calibrations
across experiments, without a separate fitted law for each case. No
comparison here establishes agreement with all public observations or at
every physical scale.

The executable model of OPH runs the [OPH law](../../model/README.md)
([specification](../../simulator/SPEC.md)). Its tests compare
distributions of admissible observer records with what observers in our world
record. The state of the simulation is no physical space: the flagship's
proposition on smooth effective field descriptions proves that observers who
inhabit the simulation can describe their experiences as relativistic 3+1D
continuum spacetimes. The observables compared below come from the
theory's conditional maps; deriving them from the model's own records is work
in progress.

Two kinds of claim need separate interpretation. Fixing a parameter that a
standard model leaves free is a testable restriction of that model. Predicting
an intrinsic vacuum dispersion or another incompatible response changes its
dynamics. A successful parameter comparison does not establish incompatible
dynamics, and agreement alone does not uniquely identify OPH.

This audit is a retrospective assessment. It does not arm or change any
[frozen prediction](../FROZEN_PREDICTION_LADDER.md), replace an eligible future
dataset, or assign a prospective verdict to already exposed observations.

The coupling/tilt calculation has no instrument input from the executable
model. Its physical audit must address the
[hadronic current and scheme matching](../../code/P_derivation/README.md),
and the [primordial construction's](../../cosmology/oph_cosmology_finite_source_cmb_program.tex)
round-sphere generator, infinitesimal reserve density, release law and
screen-to-volume map. The ACT result below uses the combined conditional
Lambda-CDM posterior, including Planck, lensing and DESI DR1; its tail is
neither an ACT-only likelihood nor a model-independent rejection probability.
The fixed-background likelihood slice is a separate overlapping-data
diagnostic.

| Candidate | Additional prediction and physical premise | Public-data disposition |
| --- | --- | --- |
| Primordial scalar tilt | `n_s = 1 - P/48 = 0.9660214956374175...`, with the inherited alpha-calibrated `P = 1.6309682094039593`. The physical primordial generator, mode identification and standard CMB transfer remain conditional. | The official combined posterior places the fixed value at its lower 0.612% position. The ACT-only conditional likelihood slice below holds the background and absolute calibration fixed while fitting amplitude and polarization calibration. These overlapping-data analyses are retrospective; neither issues a prospective verdict. |
| Intrinsic photon dispersion | Equal-weight primitive-port or seam-edge propagation gives linked fourth- and sixth-order coefficients below, if these directions are physical translations of the photon field. | A conditional departure from minimal locally Lorentz-invariant vacuum propagation. Published timing limits are far too weak at the proposed scale; the Auger transfer is source- and interaction-dependent. |
| Fine-structure constant | Several declared endpoint maps exist, with incomplete physical transport and map-selection obligations: the bare source map's Thomson-endpoint root `136.994835177` (300.388218 ppm below CODATA), whose hadronic part is the declared quark-screening ansatz `S = 1 - N_c alpha_3/pi` (about `0.887`, with no perturbative derivation); the gauge-width map's root `137.035660137` (2.474 ppm below); the separately evaluated spectral candidate `137.035999189566` (`+0.60` reference standard uncertainties); and, replacing the ansatz by the measured hadronic payload, the endpoint `136.3827548` on `[136.3670, 136.3985]` (0.48 % low; CODATA lies outside the interval). | Withhold a new prediction claim. The bare source map has the 300.388218 ppm gap; the [Thomson endpoint](../../code/P_derivation/runtime/empirical_thomson_endpoint_current.json) would need a hadronic payload of `0.032244` against the measured `0.027609 ± 0.000112`. No map carries a justified theory-error distribution, and P is located by the measured coupling. |

## Executed CMB comparison

The [executable reanalysis](../../code/cosmology/tilt_posterior/README.md)
uses the collaboration's actual ordered posterior samples, not a Gaussian
distribution manufactured from its published mean and standard deviation.
We downloaded the official `p-actbase-l-b_lcdm_camb` archive: full
multifrequency ACT DR6 plus Planck primary CMB, ACT/Planck lensing and DESI
2024 (DR1) BAO. The model is flat six-parameter power-law LambdaCDM with the
collaboration's nuisance treatment and priors, including `n_s` uniform on
`[0.85,1.1]`, at pivot `0.05 Mpc^-1`. These are conditional cosmological
inferences, not model-independent direct measurements of a repair event.
[ACT release](https://act.princeton.edu/act-dr6-data-products),
[likelihood and parameter paper](https://arxiv.org/html/2503.14452v1).

The [official burn-in instructions](https://portal.nersc.gov/project/act/dr6.02/chains/README_chains.html)
specify 50% of stored rows for this CAMB product. Applying that fraction
separately to all four chains leaves 379,175 rows with total multiplicity
508,985. The result is:

| Quantity | Actual weighted-chain result |
| --- | ---: |
| Posterior mean `n_s` | `0.974258242` |
| Posterior standard deviation | `0.003377503` |
| Fixed edge-center branch | `0.9660214956374175...` |
| Posterior mass at or below that branch | `3114/508985 = 0.00611805849...` |
| Alternative clock branch `1-e(P-phi)` | `0.96484114303...` |
| Posterior mass at or below that alternative | `1020/508985 = 0.00200398833...` |

The primary branch lies in the lower 0.612% of the empirical posterior.
This indicates tension with this particular dataset combination and its
cosmological assumptions. It is not a frequentist p-value, the probability
that OPH is true, a likelihood ratio or a Bayes factor. The alternative is
retained as an exposed diagnostic, not substituted for the fixed branch.
This posterior reduction fits no parameters or likelihoods to improve either
agreement.

The mean and standard deviation reproduce the published rounded
`0.9743 +/- 0.0034`, while the measured tail differs from a Gaussian shortcut.
The four individual edge-branch tail positions range from about 0.00373 to
0.01031. Only 85 contiguous tail excursions occur in the retained ordering;
the many repeated `n_s` values arise during nuisance-parameter updates.
Rows and multiplicities are not independent observations. Per-chain and
leave-one-chain-out sensitivity is retained; no rare-tail effective-sample-size
or confidence certificate is claimed. The release's ordinary convergence
checkpoint is not such a certificate.

The 163.5 MB source archive remains external. A 567 kB lossless projection
retains every original weight/`n_s` token, all chain ordering and burn-in rows,
with archive and original-member hashes. The independent archive parser and
exact-rational replay agree on both integer tail totals. Reproduce locally:

```sh
python code/cosmology/tilt_posterior/audit.py --check
python code/cosmology/tilt_posterior/likelihood_slice.py --check
python -m pytest -q code/cosmology/tilt_posterior
```

These observations were released in 2025 and were already exposed before
FZ-15. Its registered comparison concerns a specified future Simons
Observatory product, with additional convergence and eligibility conditions.
This analysis changes none of those terms and issues no FZ-15 verdict.
It supplies an actual-data constraint on the declared primordial branch;
deriving that branch from the executable model's agreed records is work in
progress.

The separate [conditional likelihood slice](../../code/cosmology/tilt_posterior/README.md#conditional-act-only-likelihood-slice)
uses the [official foreground-marginalized ACT likelihood](https://github.com/ACTCollaboration/DR6-ACT-lite)
with 135 measured TT/TE/EE bandpowers, their exact windows and full covariance.
It fixes `omega_b=0.0226`, `omega_c=0.118`, `theta_MC=0.010417`,
`tau=0.0544` and `A_act=1`, and fits primordial amplitude and polarization
calibration under each tilt choice. CAMB supplies flat Lambda-CDM transfer,
neutrino, recombination and nonlinear-lensing assumptions; both fitted points
receive full recalculations. The retained data projection is independently
hash-pinned to the [measured SACC release](https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/dr6_data_cmbonly.tar.gz),
and the offline check recomputes the covariance scores. It does not certify
a global optimum or rerun the radiation-transfer calculation.

The retained full-recalculation scores are `chi^2=171.929744332` for the
fixed tilt and `167.860973954` for the free tilt, whose fitted value is
`n_s=0.977274822`. The conditional difference is `4.068770378` in favor
of freeing the tilt. This is an actual-data likelihood comparison under
the stated fixed background, with no assigned significance or OPH-wide
verdict.

The [finite-source CMB construction](../../cosmology/oph_cosmology_finite_source_cmb_program.tex),
[inflation construction](../../cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.tex)
and [vacuum/structure construction](../../cosmology/oph_cosmological_vacuum_and_structure_formation.tex)
therefore support conditional early-Universe tests. A full cosmological and
nuisance profile or marginalization, with the same parameter treatment for
fixed and free tilt, is unperformed. The conditional slice fixes several of
those parameters and cannot replace it. It also overlaps the ACT data used in
the combined posterior; the two calculations are not independent evidence.
Neither supplies a p-value or changes FZ-15. Running and tensor initial
conditions retain their own premises;
vanishing source forcing alone does not imply zero primordial tensors.
The executable model supplies no curvature perturbation, release state,
reserve generator, expanding background or photon-transfer map for these
papers; that attachment is work in progress.

## Conditional dispersion: exact branch relation and actual sensitivity

In natural units, write the candidate physical dispersion as

```text
omega^2 = k^2 + C4 k^4 + [B0 + B6 I6(direction)] k^6 + O(a^6 k^8).
```

`I6` is the fixed normalized icosahedral angular invariant, not particle spin.
The [primitive-port producer](../../code/a5_fingerprint/spin_six_primitive_port_prediction.py)
and [seam-edge producer](../../code/a5_fingerprint/seam_current_edge_prediction.py)
give distinct alternatives:

| Coefficient | FZ-11 primitive ports | FZ-12 seam edges |
| --- | --- | --- |
| `C4` | `-a^2/20` | `-a^2/20` |
| `B0` | `a^4/840` | `a^4/840` |
| `B6` | `2 a^4/7875` | `-a^4/12600` |
| `B0/C4^2` | `10/21` | `10/21` |
| `B6/C4^2` | `32/315` | `-2/63` |

The finite geometry does not select the physical translation action, sole
propagation support, photon identification, kinetic term, physical clock or
preferred-frame transport. Those are premises of these branches. In
particular, taking the square root of a positive spatial operator as a wave
frequency is additional dynamics; it is not a consequence of the finite
occupation process relaxing by that operator.

Under the existing additional scale assumption `a = sqrt(P) l_P`, with
`l_P = 1.616255e-35 m` and `hbar c = 1.973269804e-7 eV m`, direct substitution
gives

```text
d = a^2/20 = 5.470944672629916e-58 eV^-2,
E_QG,2 = 1/sqrt(d) = 4.275322088547057e19 GeV,
v_group = 1 - (3/2) d E^2 + O(a^4 E^4).
```

The LHAASO collaboration's GRB 221009A analysis reports the subluminal
quadratic bound `E_QG,2 > 6.9e11 GeV` at 95% confidence using its ML/MINOS
procedure; its calibrated variant gives `7.2e11 GeV`. These are published
results under the paper's light-curve, spectral, propagation and detector
assumptions, not a new fit to the events here. [LHAASO, Table 1](https://arxiv.org/html/2402.06009v2)

For the quoted `6.9e11 GeV` bound, the proposed coefficient divided by the
bound's coefficient is `(6.9e11 / 4.275322088547057e19)^2 = 2.60472e-16`.
Thus this timing comparison is compatible but has negligible power to test
the assumed scale. It supplies no empirical preference over zero dispersion
and no test of the linked sixth-order angular coefficients. No claim that
this is the strongest currently available timing analysis is needed.

The [Auger electromagnetic analysis](https://arxiv.org/abs/2112.06773)
instead reports `delta_gamma,2 > -1e-58 eV^-2` under a high-energy
proton-containing source scenario, and no electromagnetic bound for its
reference composition without protons beyond `1e19 eV`. In the photon-only
branch `delta_gamma,2 = -d` to leading order, so the candidate magnitude is
5.47094 times that quoted limit. This is not an OPH exclusion: source
composition, lepton dispersion, interactions and shower response determine
whether that bound applies. The [existing feasibility calculation](DISPERSION_FEASIBILITY.md)
shows that shared photon/lepton dispersion substantially changes the pair
threshold; free shells alone do not supply the interaction rates.

A decisive test needs an isolated intrinsic propagation response and the
joint coefficient uncertainties. A resolved positive `C4`, or incompatible
linked sixth-order coefficients, can reject the specified cosine branch.
A null can reject a fixed-scale branch only with sufficient sensitivity.
If arbitrarily small positive `a` remains allowed, finite-precision null
measurements cannot exclude the entire family, whose closure contains zero.
Even a matching nonzero pattern would require comparison with other
Lorentz-violating theories and systematic explanations before supporting a
distinctive OPH interpretation.
