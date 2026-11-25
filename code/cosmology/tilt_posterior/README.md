# ACT DR6 retrospective fixed-tilt comparison

The conditional edge-center branch gives
`n_s = 1 - P_C/48 = 463690317905960407/480000000000000000`, or
approximately **0.966021496**. In the official ACT DR6 P-ACT-LB baseline
posterior, approximately **0.612%** of the retained multiplicity is at or below
this value. The historical clock alternative gives approximately **0.964841143**,
with **0.200%** at or below it. These are posterior positions conditional on the
published cosmological model and priors; they are not p-values, Bayes factors,
probabilities that OPH is true, or evidence that the native observer model has
derived our Universe.

This is a reproducible reanalysis of measured-Universe data already used in the
[cosmology postdiction ledger](../postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md).
It includes exact weighted marginal calculations from the released ordered
chains, with per-chain and leave-one-out sensitivity, and a separate conditional
ACT-only likelihood slice below. These analyses use overlapping ACT data and
are not independent evidence. Both exposed formulas remain visible in the
posterior comparison; the clock alternative does not replace the edge-center
commitment.

## Posterior data and provenance

The [ACT collaboration release page](https://act.princeton.edu/act-dr6-data-products)
links the official [NASA LAMBDA archive listing](https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_chains_lcdm_get.html)
and [NERSC chain documentation](https://portal.nersc.gov/project/act/dr6.02/chains/README_chains.html).
The selected product is exactly `p-actbase-l-b_lcdm_camb`, not the smaller CLASS
or foreground-marginalized ACT-lite alternatives. It combines Planck and the
full multifrequency ACT DR6 primary CMB likelihood, ACT/Planck CMB lensing,
and DESI 2024 (DR1) BAO. No Pantheon+ supernova likelihood is included.

The [ACT DR6 parameter paper](https://arxiv.org/html/2503.14452v1), section VIII.3
and Appendix H, supplies the physical interpretation. The dimensionless scalar
spectral index describes the primordial scalar power law with amplitude pivot
`k_* = 0.05 Mpc^-1`. The chain assumes baseline six-parameter flat Lambda-CDM
and the collaboration's standard radiation transfer, recombination, foreground,
calibration, and nuisance models. Its configuration gives a uniform `n_s` prior
on `[0.85, 1.1]`. The normal reference distribution in that configuration only
initializes the sampler; it is not the prior. Marginalizing the released joint
posterior retains the effect of its fitted parameter correlations, without
inventing a diagonal observational covariance.

The archive contains four ordered Cobaya chains, configuration, and convergence
records. This package retains only the original ASCII `weight` and `ns` tokens
from every numbered-chain row, together with member sizes/hashes and the
checkpoint. There is no thinning, rounding, or removal of initial rows in the
stored projection. The full 163.5 MB archive remains outside the repository.

| Artifact | Size in bytes | SHA256 |
| --- | ---: | --- |
| Official archive | 163494761 | `8f2801a69663ef1142bb05660e0118ef7b1e0bdf3d05afba434a0569023b98ab` |
| Retained `act_dr6_ns_weight.json.gz` | 567477 | `13dab540b652300bd000cf1ac3b197af55c7c411de215682482ea55993f40485` |

The archive is available from
[NASA LAMBDA](https://lambda.gsfc.nasa.gov/data/act/chains/lcdm/p-actbase-l-b_lcdm_camb.tar.gz)
and its [NERSC mirror](https://portal.nersc.gov/project/act/dr6.02/chains/compressed/lcdm/p-actbase-l-b_lcdm_camb.tar.gz).
The listing dates delivery to 18 March 2025, before the FZ-15 stamp of
4 October 2026. The data are retrospective relative to that commitment. The projection is a
lossless projection onto two columns, not a copy of the complete posterior
state or a replacement for the full likelihood.

## Posterior calculation and its limits

The official README prescribes **50% burn-in for this exact CAMB product**.
The calculation discards `floor(0.5 * stored_row_count)` initial rows separately
from each chain, then applies the remaining original integer multiplicity
weights. It does not remove half the total weight and does not apply the 10%
burn-in recommended for the separate CLASS product. Of 758350 source rows,
379175 remain, with total multiplicity 508985. Comparisons with the fixed
rational edge prediction use the exact decimal tokens, without a float
comparison at the decision boundary. Quantiles invert the weighted empirical
CDF without interpolation.

| Retrospective result | Value |
| --- | ---: |
| Pooled posterior mean of `n_s` | 0.9742582419612366 |
| Pooled posterior standard deviation | 0.0033775031447940 |
| Edge-center mass at or below prediction | `3114/508985 = 0.00611805849` |
| Historical clock-alternative mass at or below prediction | `1020/508985 = 0.00200398833` |
| Edge-center per-chain lower-tail mass range | 0.00372945 to 0.01031207 |
| Clock-alternative per-chain lower-tail mass range | 0.00077226 to 0.00377020 |

The mean and standard deviation reproduce the published rounded
`0.9743 +/- 0.0034`. The fixed branch lies in the lower posterior tail, so this
comparison supplies tension rather than a numerical confirmation. A tail mass
alone does not quantify a model's relative evidence.

The released checkpoint reports `converged: true` with
`Rminus1_last = 0.009799812902692016` and four MPI chains. This ordinary
convergence record does not certify rare-tail accuracy. Many consecutive rows
have identical `n_s` during fast nuisance-parameter updates. The 2319 retained
rows below the edge prediction occur in only 85 contiguous excursions, so
neither stored rows nor summed multiplicities count independent observations.
Per-chain and leave-one-out sensitivity is reported without treating it as a
confidence interval. This package supplies no rank/folded split-Rhat or tail-ESS
certificate and cannot satisfy the complete
[FZ-15 adjudication requirements](../../../docs/FROZEN_PREDICTION_LADDER.md).
No prospective FZ-15 verdict is issued or changed.

The model-side coordinate `P_C = 1.6309682094039593` is inherited from a
CODATA-alpha calibration. A physical primordial generator, mode embedding,
and standard CMB transfer remain additional hypotheses. No joint uncertainty
propagation over that calibration and these missing theoretical premises is
claimed. Standard Lambda-CDM treats `n_s` as a fitted primordial input; the
conditional OPH branch fixes that input, but this comparison does not establish
its derivation from the runnable native observer construction.

## Conditional ACT-only likelihood slice

[likelihood_slice.py](likelihood_slice.py) uses the official
[foreground-marginalized ACT DR6 likelihood](https://github.com/ACTCollaboration/DR6-ACT-lite/tree/627aeafb88ae5ad1aa66b406bea2d65cfa66a27d)
and [measured SACC archive](https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/dr6_data_cmbonly.tar.gz).
It compares the fixed edge-center tilt with a free tilt using all 135 TT/TE/EE
bandpowers, their windows and full covariance. The data are measured ACT
observations; the repository's bundled mock data are inadmissible.

The calculation fixes `omega_b=0.0226`, `omega_c=0.118`,
`theta_MC=0.010417`, `tau=0.0544` and absolute calibration `A_act=1`.
Each branch optimizes `log(10^10 A_s)` and polarization calibration `P_act`;
the second branch also optimizes `n_s`. Their bounds are `[2.6,3.5]`,
`[0.9,1.1]` and `[0.9,1.1]`, respectively. These optimization bounds are
not a posterior marginalization. The pivot is `0.05 Mpc^-1`. The supplied
flat Lambda-CDM transfer model has total neutrino mass `0.06 eV`,
`N_eff=3.044`, one massive neutrino, BBN-consistent helium, zero scalar
running and tensors, `A_lens=1`, Recfast recombination and HMcode-2020
nonlinear lensing. The full settings are retained in
[likelihood_slice.json](likelihood_slice.json).

CAMB transfer reuse sets `only_time_sources=True`, so changing amplitude
or tilt recalculates nonlinear-lensing corrections as required by the
[CAMB implementation contract](https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.power_spectra_from_transfer).
Both fitted points also receive independent full CAMB recalculations.
The retained scores use those full calculations and independent Cholesky
whitening; agreement with the official ACT quadratic likelihood is checked.

| Conditional fit | `n_s` | `log(10^10 A_s)` | `P_act` | Full-recalculation `chi^2` |
| --- | ---: | ---: | ---: | ---: |
| Fixed edge-center tilt | 0.966021496 | 3.046192764 | 1.003886089 | 171.929744332 |
| Free tilt | 0.977274822 | 3.037913732 | 1.002521246 | 167.860973954 |

Thus `chi^2_fixed - chi^2_free = 4.068770378`: freeing the tilt improves
this conditional fit. Full recalculation changes either cached-transfer
score by less than `4e-13`; this is a numerical consistency check, not a
bound on cosmological modeling uncertainty.

This is a numerical likelihood comparison conditional on the specified
background and calibration. It is not the collaboration's joint cosmological
fit, a fully profiled or marginalized cosmological/nuisance likelihood, or a
certified global optimum. No p-value or FZ-15 verdict is assigned. Its ACT
observations overlap the posterior analysis above; the two results cannot be
combined as independent evidence. The supplied primordial-to-CMB map also
does not identify the native observer's field with primordial curvature.

The receipt retains every observed and predicted bandpower and the full
covariance triangle. Independent extraction from the measured FITS pins the
ordered bin indices, ell values, tracers and window hashes together with
those observations. The offline checker verifies that projection, metadata
and both quadratic scores; it does not rerun CAMB or establish optimality.

## Reproduce

From the repository root, obtain the pinned archive outside the repository:

```sh
curl --fail --location \
  https://portal.nersc.gov/project/act/dr6.02/chains/compressed/lcdm/p-actbase-l-b_lcdm_camb.tar.gz \
  --output /tmp/p-actbase-l-b_lcdm_camb.tar.gz
python3 code/cosmology/tilt_posterior/project_archive.py \
  /tmp/p-actbase-l-b_lcdm_camb.tar.gz \
  /tmp/act_dr6_ns_weight.json.gz
```

The generator checks the archive size and pinned SHA256 before processing. It
streams regular tar members without extracting paths or executing configuration
code, hashes the original member bytes, preserves original token strings,
and writes sorted compact JSON with gzip timestamp zero and no filename header.
The rebuilt projection reproduces the retained bytes and digest. Only Python's
standard library is required for this projection step.

Replay the retained comparison using the repository's Python dependencies:

```sh
python3 code/cosmology/tilt_posterior/audit.py --check
python3 code/cosmology/tilt_posterior/likelihood_slice.py --check
```

To analyze the independently rebuilt projection without changing the retained
receipt:

```sh
python3 code/cosmology/tilt_posterior/audit.py \
  --projection /tmp/act_dr6_ns_weight.json.gz \
  --output /tmp/act_dr6_tilt_results.json
```

`results.json` records the input and analysis hashes, exact comparisons,
weighted posterior summaries, and sensitivity checks. CI can replay it from
the small retained projection without a network download.

The optional likelihood reproduction requires CAMB `1.6.6`, SACC `0.16`,
NumPy `2.2.6`, SciPy `1.15.3`, Cobaya `3.6.2` and the other official
likelihood dependencies, together with the pinned ACT source linked above. Install the measured archive under
`<data-root>/data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits`, then run:

```sh
python3 code/cosmology/tilt_posterior/likelihood_slice.py \
  --data-root /path/to/cobaya-packages \
  --official-code /path/to/DR6-ACT-lite \
  --output /tmp/act_dr6_likelihood_slice.json
```

The producer verifies measured-data and official-source hashes before
optimizing. Ordinary offline CI needs no download or optional cosmology
packages. The independent covariance/calibration tests and mutation gates run
with `python3 -m pytest -q code/cosmology/tilt_posterior`.

## Planck acoustic peak comparison

The [peak table](planck_peaks.json) gives the measured and calculated multipole and signed bandpower for each extremum, with the published marginal errors on the measured band. The [calculation](planck_peaks.py) propagates the conditional repair tilt and a control tilt through identical scalar adiabatic CAMB 1.6.6 transfer. Both use the Planck-calibrated expansion and matter parameters and the supplied primordial amplitude $A_s=2.1\times10^{-9}$. All three spectra are lensed, in thermodynamic microkelvin squared. The theory is sampled at the effective bin centres. The first five positive maxima of the smooth control and the first five negative TE minima each define a window of 60 multipoles on either side. Each curve locates its own extremum within that window.

| Spectrum | Bins matching Planck, repair tilt | Bins matching Planck, control | Repair-tilt height ratios relative to the first extremum, fractional error in percent |
| --- | ---: | ---: | --- |
| TT maxima | 5/5 | 5/5 | +1.51, +1.53, +2.03, +2.38 |
| TE maxima | 5/5 | 5/5 | +10.57, −2.88, +75.26, −28.35 |
| TE minima | 5/5 | 5/5 | +4.14, +9.78, +5.31, +2.59 |
| EE maxima | 4/5 | 4/5 | −1.19, −3.63, −3.71, −6.93 |

The third EE maximum is one bin apart. The large fractional difference for the fourth positive TE crest compares a calculated 6.49 with a measured $3.70\pm2.49$ in the same bin. Peak picking in noisy data, finite bin widths and inter-bin covariance limit the interpretation of these ratios. Inter-observer records would have to fix the primordial normalization, the acoustic initial conditions and the photon transfer for the same comparison to test a repair-law sky.

Download the three small binned text files identified by URL and SHA-256 in the peak table to one directory. From the repository root, run `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python code/cosmology/tilt_posterior/planck_peaks.py DATA_DIRECTORY`. [Planck public release](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/)
