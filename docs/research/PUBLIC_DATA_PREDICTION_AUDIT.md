# OPH predictions and public-data comparisons

The target is a prediction for our observed Universe that adds information
beyond the Standard Model and General Relativity, followed by an executable
comparison with existing observations. Building a new laboratory realization
of the finite source is not a prerequisite. The observable dictionary, physical
premises and uncertainty model must nevertheless be explicit.

Current coverage consists of the selected observations listed below. The
target is one fixed construction with shared parameters and declared
calibrations across experiments, not a separate fitted law for each case.
No comparison here establishes agreement with all public observations or
at every physical scale. The [common local-clock observer](../../model/native_modular_reference/README.md#a-common-observer-from-two-local-stopping-clocks)
derives one finite preparation and continuing readout contract, with native
scheduler bits and explicit failure budgets. Its shared gaps `log(25/8)` and
`log(8/5)` are internal dimensionless frequencies. Neither that construction
nor its finite theorem connections identify the atomic transitions, detector
response, laboratory seconds or primordial modes used in the comparisons
below. Those interfaces cannot be selected separately to force agreement.

The [same-parameter replay](../../model/predictions/shared_parameter.py)
uses the inherited structured-running branch's forward numerical enclosure
`P ≈ 1.630972095858897377` in both its coupling endpoint and the conditional
map `n_s=1-P/48`, with no recalibration. This gives
`n_s ≈ 0.9660214146696063`, a shift of `−8.09678112e−8` from the
alpha-calibrated tilt below. Both enclosure endpoints give exactly the same
ACT posterior tail weight, `3114/508985`; harmonizing P does not remove that
tension. The same branch retains the inverse-alpha deficit of `300.388218`
ppm. The prediction payload is sealed before observations are loaded, and
the [receipt](../../model/predictions/shared_parameter_results.json) records
the source and data hashes. This retrospective diagnostic shares a parameter;
the low-energy transport and primordial/transfer maps remain conditional.
Its numerical enclosure supplies no physical theory-error budget, and no
joint likelihood or all-scale agreement verdict is assigned.

The coupling/tilt calculation has no C/E instrument input. The
[P-response observer](../../model/observer/p_response/experiment.py) imports
the structured-running closure value into its detector; it does not derive
P from either repair instrument. The [C/E readout test](../../model/predictions/observer_instrument.py)
compares their distinct coherent laws, but supplies no repair-dependent
coupling or primordial generator. Thus the alpha/tilt comparison is not a
separate empirical test of C and E. Its physical audit must address the
[hadronic current and scheme matching](../../code/P_derivation/README.md),
and the [primordial construction's](../../cosmology/oph_cosmology_finite_source_cmb_program.tex)
round-sphere generator, infinitesimal reserve density, release law and
screen-to-volume map. None is repaired by switching an instrument label.
The ACT result below uses the combined conditional Lambda-CDM posterior,
including Planck, lensing and DESI DR1; its tail is neither an ACT-only
likelihood nor a model-independent rejection probability. The retained
fixed-background likelihood slice is a separate overlapping-data diagnostic.

The [earlier-state readout test](../../model/interacting_band_extension/README.md#a-held-out-record-law-after-central-weight-calibration)
fixes one Hamiltonian, temperature and instrument family, infers its remaining
sector weight from complete records, and predicts a separate model-time
readout without another fit. The source execution passes and four changed
validation contracts fail. Its illustrative counts are synthetic; the
sampling and service-error bounds expose weak sensitivity of this short
protocol. This is an executed model prediction, pending a physical dataset
with the same preparation and readout contract.

Two kinds of claim need separate interpretation. Fixing a parameter that a
standard model leaves free is a testable restriction of that model. Predicting
an intrinsic vacuum dispersion or another incompatible response changes its
dynamics. A successful parameter comparison does not establish incompatible
dynamics, and agreement alone does not uniquely identify OPH.

This audit is a retrospective assessment. It does not arm or change any
[frozen prediction](../FROZEN_PREDICTION_LADDER.md), replace an eligible future
dataset, or assign a prospective verdict to already exposed observations.

| Candidate | Additional prediction and physical premise | Public-data disposition |
| --- | --- | --- |
| Primordial scalar tilt | `n_s = 1 - P/48 = 0.9660214956374175...`, with the inherited alpha-calibrated `P = 1.6309682094039593`. The physical primordial generator, mode identification and standard CMB transfer remain conditional. | The official combined posterior places the fixed value at its lower 0.612% position. The ACT-only conditional likelihood slice below holds the background and absolute calibration fixed while fitting amplitude and polarization calibration. These overlapping-data analyses are retrospective; neither issues a prospective verdict. |
| Charged-tau mass | Balanced-circulant mass interpretation and measured electron/muon masses give approximately `1776.969027 MeV`. This relation is algebraically the Koide relation; the balance premise has observed-mass ancestry. | A conditional mass restriction, not a change to SM dynamics or a unique OPH signature. The historical Belle II 2023 comparison below is compatible at its reported precision. |
| Executed chiral particle candidate | One protected 630-mode field has complete finite source/clock records. A distinct local two-hop completion gives the same restricted seven-momentum histories; a second-shell joint read separates the two fixed laws. A separate Euclidean overlap operator implements a gauge-covariant chiral kinetic projector. Reference-preserving passive repair schedules leave a rank-two family seed in the tested Higgs route. | The [mass-response audit](OPH_PARTICLE_MASS_RESPONSE.md) rejects the chosen charged hierarchy. The second-shell contrast is a model test without a public particle dataset. The kinetic regulator does not supply a chiral measure, interacting continuum, physical clock, source-selected Yukawas or observed-particle attachment; no cross-scale agreement verdict follows. |
| Intrinsic photon dispersion | Equal-weight primitive-port or seam-edge propagation gives linked fourth- and sixth-order coefficients below, if these directions are physical translations of the photon field. | A conditional departure from minimal locally Lorentz-invariant vacuum propagation. Published timing limits are far too weak at the proposed scale; the Auger transfer is source- and interaction-dependent. |
| Calibrated matter control | One supplied common-scale Standard Model parameter point passes exact matter/current identities. External SMDR inverse/forward matching reproduces the fitted charged mass inputs; the unfitted W value misses the selected one-standard-uncertainty interval. | [Calibration evidence](../../model/common_source/calibrated_matching_results.json), not native observer records or a mass derivation. Fitted-input reuse, omitted CKM matching effects and the absent joint theory/input uncertainty budget prevent an independent physical verdict. |
| Calibrated observer control | A decoder without the mass packet encloses all nine supplied running gaps from 84 million simulated complete rest-field/clock records, with timeouts retained. | [Source and record replay](../../model/common_source/calibrated_spectrum_observer_results.json). Conditional sampling and admitted service/frame errors; several intervals exceed measured-mass precision. This is a calibrated readout sanity check, not independent physical data, a parameter derivation or an interacting spectrum. |
| Fixed repair parameter selection | The complete corrected logical transfer is exactly identity; native one-port and edge observables become scalar. A two-hop correlation distinguishes some code levels but leaks from the code. | [Exact source audit](../../model/common_source/parameter_selection_results.json), using no calibrated packet. Neither its physical coherent coefficient nor a three-family mass dictionary is derived. |
| Specified coherent continuation | A supplied native two-hop phase and continuing schedule have a certified effective encoded Hamiltonian under the actual compiled correction. The finite-proposal bound includes arbitrary references and a separate accumulated service-error term. | [Continuation certificate](../../model/common_source/repair_continuation_results.json). Terminal stopping and within-cycle phase ambiguity are tested separately. This is a conditional dynamics derivation, not a physical mass or clock identification. |
| Native correlation/exchange field | A fixed orbit average and whole-carrier exchange give independent fermion modes, charge continuity and a density interaction. Twelve million reduced-processor detector records enclose the native gaps without measured-mass inputs. | [Source and record control](../../model/common_source/repair_field_attachment_results.json). The gap ratios 1:3:12 fail the charged-lepton hierarchy. Local native attachments are separate from a complete service history; encoded charge, one-dimensional locality and model time are not Standard Model identifications. |
| Persistent native geometric object | Two record-correctable equivariant matching codes supply actual geometric 3 + 3prime, with source preparation on the same opposite-code pair used in the interaction test. | [Object](../../model/common_source/native_gauge_object_results.json), [source execution](../../model/common_source/native_code_transport_results.json) and [response](../../model/common_source/native_persistent_response_results.json). The tested preserving response has dimension 17 and fails the complete 12-current premise. No physical particle, mass or public-data agreement follows. |
| Specified same-object gauge/matter law | A chosen coherent port law gives twelve closed currents on the geometric object; its declared spin/exterior lift supplies fifteen internal matter weights and three Higgs lines. | [Current law](../../model/common_source/native_sm_current.py) and [matter attachment](../../model/common_source/native_matter_attachment.py). Completeness is over the declared response menu; universal apparatus controls remain explicit. No measured masses enter this algebraic construction, and no additional empirical agreement is claimed. |
| Native dimer spectroscopy | Fixed occupation preparation and transfer quadratures identify both frequencies using four training and two held-out model-coordinate samples, retaining 158 terminal records. | [Spectroscopy receipt](../../model/common_source/native_charged_spectroscopy_results.json). The inverse receives no target pole; the action supplies its coefficients. Numerical complete-record contrasts do not constitute experimental mass precision or an additional public dataset. |
| Native free cubic field and composed history | Specified inverse/phase/parity controls derive a staggered field, two Dirac copies per sector and a fixed-band continuum contract. A finite source cube carries all 256 Fock columns. A distinct two-carrier protocol retains 319 complete native detector/repair outcomes. | [Free field](../../model/common_source/native_staggered_field_results.json) and [native history](../../model/common_source/native_field_history_results.json). The free ratios remain 1:3:12. The two-carrier number law probes exchange/backaction but cannot identify its uniform onsite gap. Neither supplies chiral Standard Model selection, physical service execution or an additional public dataset. |
| Interaction-dressed native poles | The full alternating-sign exchange dimer has exact charged poles sqrt(J²+m²)−J and nonzero CAR residues. An exclusion/Poisson calculation controls finite-volume response through a stated perturbative domain. | [Charged-spectrum receipt](../../model/common_source/native_charged_spectrum_results.json). Ratios 1:8.999872:143.963411 fail both broad lepton windows. The volume-dependent coefficient cannot be extrapolated to a bulk mass outside its domain; no empirical parameter was fitted. |
| Moving calibrated field | One 1,170-CAR field covers rest and two momentum shells in one clock. Fixed preparation/readouts and classical phase tracking remove mass-dependent probe controls. | [Complete-record replay](../../model/common_source/calibrated_transport_results.json). The evolution retains calibrated parameters; conditional energy and invariant-mass intervals do not establish a physical spectrum or measured-mass precision. |
| Same-action fermion entropy | One regulated calibrated free action supplies the regional vacuum, modular response and a separate 32-CAR source preparation. | [Regional receipt](../../model/common_source/fermion_ball_results.json). No independent physical dataset, interacting gauge-region algebra, native vacuum selection or joint entropy--geometry stationarity is supplied. |
| Reciprocal thermal geometry | Fixed-gap thermal collisions derive a joint Gibbs fixed point and nonuniform geometric probabilities for the finite scalar control and a fixed-charge native fermion-chain sector. | [Reciprocal source control](../../model/native_equilibration/reciprocal_vacuum_results.json). Hamiltonians, temperature, reservoir preparation and geometry coupling are supplied. Probability-space free-energy stationarity is not metric entropy--area stationarity, and no additional public dataset is reproduced. |
| Zero-temperature free-field wedge | The same uncalibrated staggered action supplies a filled vacuum, half-space modular/stress relation, continuum normalization bound and local preparation-size schedule. A small source preparation has an independent exact binary fidelity certificate; transverse modes give a bare planar entropy density. | [Wedge receipt](../../model/native_equilibration/wedge_modular_results.json). The half-space theorem and prescribed free vacuum do not supply interacting matter, a renormalized Newton coefficient or joint small-ball metric stationarity. These are mathematical/source controls without an additional public measurement comparison. |
| Native clock ratios and noise | Same-history movement and incident counts have definite waiting laws and spatial correlations. Under simultaneous source refinement, bounded-incidence clock jitter has Gaussian finite-dimensional covariance fixed by shared seams and is asymptotically independent of density snapshots. | The [joint clock/field prediction](../../model/relative_observer_geometry/README.md#incident-clocks-on-the-simultaneous-density-continuum-scale) uses explicit stationary preparation, IID seams and readout access. No identified atomic-clock, pulsar or photon observable has this joint record law. Ordinary stochastic realizations also satisfy it; matching it alone would not distinguish OPH. |
| Native record-to-phase relation | Four movement-bit covariances on the complete native attempt tape determine a two-mode transfer generator and predict its coherent quadratures under the fixed readout protocol. | A conditional cross-readout prediction with uncertainty and failed-identification outcomes. Public optical servo ratios are not this bit tape; the [inverse-clock construction](../../model/predictions/README.md) does not assign them a new experimental verdict or establish a distinction from standard quantum instruments. |
| Prepared-state clock and response | One finite prepared density fixes eleven fluctuation pairs, their noise, ordered response and positive-log clock; three band representatives have complete source pointer reads. The native growing angular family has derived partition/gap bounds and a constructive finite-word/cap/precision schedule. A same-product quartic projection has complete interacting Gibbs readouts. | The [joint-limit contract](../../model/predictions/README.md#growing-shared-state-field) retains preparation, control, precision-service and spectral-tail premises; growing execution is unperformed. The [product field](../../model/observer_time_fields/REFINEMENT.md#one-product-field-state-and-active-clock) requires a radial design and an anchor for its same-density star construction. The [finite interacting projection](../../model/interacting_band_extension/README.md#interacting-readouts-on-the-same-radial-product) adds a declared coupling and model clock. None is an atomic-spectrum identification or a discriminator against standard quantum theory. The fixed three-gap Ar13+/Yb+ attachment is excluded only under the stated error budgets below. |
| Source-coded neutral coherence | A complete finite native C/E trial predicts total variation `29/39600`, using two current heralds, source mixers and nine pointer copies. The full error bound includes finite warm-up and binary source arrays. | The [executable instrument](../../model/native_effective_fields/README.md) retains failures, timeouts and every raw address. Existing coherent/reset/read services and global event access remain; no matching public particle or clock dataset has been identified. This discriminates two model instruments, without changing the coupling/tilt comparison or establishing a physical mass. |
| Native repair response | A disclosed-bit-triggered repair predicts complete future clock/field events from passive correlations. A controlled execution bias has continuum mean derivative `2(I-exp(-tau*A))f` and a positive fixed-amplitude signal under a uniform likelihood bound. | The [complete native forecast](../../model/predictions/README.md#native-repair-response-from-passive-clock-and-field-records) retains timeouts and distinguishes executed repairs from additional coin-service costs. Its controller and readout dictionary are declared. Public optical frequency-ratio traces contain no matching intervention/current records; no empirical verdict or distinction from ordinary stochastic or quantum realizations follows. |
| Shared oscillator preparation | One real-displaced thermal state, selected from the 21 complex axis points, predicts both quadratures at all 100 reserved off-axis points. Exact outer interval cells propagate axis uncertainty without per-point refits. | The [shared-calibration receipt](../../model/scalar_observations/experimental_readout/shared_calibration_results.json) contains a simultaneous common-state witness for each unresolved published coordinate extent, also linked to the same retained source parameters. The 5% total sampling allowance is conditional on 250 independent accepted shots per component, the declared reported-mean rounding bound and zero additional physical bias. Raw accepted counts are unavailable; this is neither authenticated uncertainty agreement nor an OPH-specific or dynamical test. |
| Executed oscillator dynamics | Native source words implement the supplied spin–oscillator apparatus, preparation, noise, detector and published duration convention. Later observations are excluded from calibration. | The [707-mean comparison](../../model/anharmonic_observations/README.md) has binary-probability RMSE `0.055582`; [125 changed-potential estimates](../../model/anharmonic_transfer/README.md) have dimensionless estimator RMSE `0.287643`. Both improve on their frozen and harmonic controls, with substantial retained discrepancies and incomplete physical error budgets. These retrospective tests of a supplied quantum apparatus do not identify the separate spatial field or native incident clock. |
| Diagonal local response | Independent settings and factorized local readouts restrict the recorded classical branch to CHSH `S<=2`, allowing shared randomness and device memory under the stated trial contract. | The [Hensen replay](../../model/predictions/BELL_RECORDS.md) reproduces 196 wins in 245 eligible trials and fixed-horizon tail bound `0.03917464`. This excludes that conditional classical extension at 5%; it neither excludes the quantum instrument nor confirms a physical OPH implementation. |
| C/E coherent record law | One fixed coherent preparation gives plus probabilities `5/8` and `3/8` under the two instruments. | These are ordinary quantum instruments. The distinction selects a record law, not OPH over standard quantum mechanics. |
| Joint timed coherent response | The charged total-two source history fixes a joint cost/clock/detector law. Under one explicitly synthetic calibration, its two settings have identical `0.52` deadline mass and distinct `0.2451`/`0.2715` joint detector masses. | An executable conditional forecast with shared timing coefficients and all timeouts retained. No measured system or physical duration is identified; these finite quantum instruments are not unique to OPH. |
| Fine-structure constant | Several declared endpoint maps exist, with incomplete physical transport and map-selection obligations: the bare source map's Thomson-endpoint root `136.994835177` (300.388218 ppm below CODATA), whose hadronic part is the declared quark-screening ansatz `S = 1 - N_c alpha_3/pi` (about `0.887`, with no perturbative derivation); the gauge-width map's root `137.035660137` (2.474 ppm below); the separately evaluated spectral candidate `137.035999189566` (`+0.60` reference standard uncertainties); and, replacing the ansatz by the measured hadronic payload, the endpoint `136.3827548` on `[136.3670, 136.3985]` (0.48 % low; CODATA lies outside the interval). | Withhold a new prediction claim. The [existing endpoint audit](../../model/predictions/README.md#inherited-coupling-comparison) finds the 300.388218 ppm gap for its named map; the [Thomson endpoint](../../code/P_derivation/runtime/empirical_thomson_endpoint_current.json) would need a hadronic payload of `0.032244` against the measured `0.027609 ± 0.000112`. No map carries a justified theory-error distribution, and P is located by the measured coupling. |
| Stipulated mass branch | A declared circulant and six-power texture, with imported running and criticality, give twelve nonzero outputs. | [Comparison receipt](../../model/particle_spectrum_audit/README.md). All twelve lie outside printed PDG intervals and the quark texture is rejected by scale-free ratios. The architecture is target-informed, so this is no prediction. |
| Higgs mass | The historical near-match used a target-conditioned calibration. The [replacement electroweak producer](../../code/particles/calibration/source_ew_vev_matching.py) uses additional criticality, matching and energy-scale premises; it is not derived from the finite observer. | Withhold a distinctive prediction claim: the physical matching and theory-error budget remain open. The replacement's [retrospective comparison](../../code/particles/runs/calibration/source_ew_vev_matching_comparison.json) should retain those qualifications. |

The [certified transfer comparison](../../model/occupation_control/README.md#complete-changed-potential-comparison)
joins the unchanged apparatus to **all 125 deposited positive-duration means**
across five potential/preparation cases. An independent 64-level evaluator
bounds the unrestricted apparatus, including occupation truncation, integration,
roundoff and detector evaluation. A trace-aware conditioning inequality retains
signed approximate states and tightens the selected-readout bounds without
changing any apparatus parameter.

| Case | Certified endpoints | Central values outside numerical intervals | Largest required additional allowance |
| --- | ---: | ---: | ---: |
| a | 25 | 14 | > 0.60059383 at 8.4 ms |
| b | 25 | 11 | > 0.42598315 at 1.6 ms |
| c | 25 | 14 | > 0.46868933 at 1.2 ms |
| d | 25 | 12 | > 0.36578125 at 2.8 ms |
| e | 25 | 9 | > 0.33039699 at 3.6 ms |

Allowances use the published dimensionless position-estimator units. Each is
computed exactly as `max(0, abs(observed-centre)-radius)` from that endpoint's
own binary path and outward bound. Thus numerical and occupation error alone
cannot cover 60 of the published central values. Experimental, calibration,
model and physical-service errors remain outside these intervals; this is
**not a statistical rejection** or a full-uncertainty agreement test.

The [shared-calibration diagnostic](../../model/occupation_control/shared_calibration.py)
allows a common nonnegative output gain and unrestricted offset on all 125
certified intervals. Even this family requires a minimum uniform additional
allowance of **0.4916177969622167** estimator units (rounded). The
[exact receipt](../../model/occupation_control/shared_calibration_results.json)
contains a feasible rational point and nonnegative dual witness; verification
checks every endpoint and both nonnegative boundaries. Allowing independent
values inside each interval makes the result conservative for the actual
shared dynamics. The same exact minimum holds with separate nonnegative gains
and unrestricted offsets for each of the five cases: the dual witness lies
entirely in case a, while the common feasible point remains feasible in the
enlarged family. Its three weighted endpoint-error allowances, at 1.2, 3.2 and
8.4 ms, must sum to at least that minimum. Thus requiring one calibration
across different potentials is not the cause of this obstruction.
The optimizing gain and offset are never adopted. This
excludes an exact calibration-only explanation within this affine family;
changing the physical probe coordinate, preparation, timing or dynamics is
outside that family. No experimental uncertainty or rejection probability is
inferred from this deterministic residual.

The same existing observations are used throughout; predictions are unchanged
and no datum is counted twice as independent evidence. Five zero-time diagnostics
remain separate. The [paper](https://arxiv.org/html/2603.04744v1#S3.SS2) reports
200 nominal measurements per point and spin-noise error bars in Figure 5.
The [deposit](https://zenodo.org/records/18752287) contains no numerical
per-point uncertainties, accepted counts, raw trials or calibration covariance. The
[retained receipt](../../model/occupation_control/transfer_results.json.gz)
binds all endpoints and source/data hashes; automatic replay recomputes the
complete numerical paths. This tests a supplied quantum apparatus, not a unique
OPH prediction or a native clock/particle identification.

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
| Historical clock alternative `1-e(P-phi)` | `0.96484114303...` |
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
deriving that branch from the runnable observer's joint records remains open.

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
The native model needs the specific curvature perturbation, release
state, reserve generator, expanding background and photon-transfer map used
by these papers. Its static two-mode product interaction supplies none of
those by identifying its scalar field with primordial curvature.

## Conditional tau comparison

For `p = sqrt(m_e) + sqrt(m_mu)` and `s = m_e + m_mu`, the ordered Koide
branch is `m_tau = [2p + sqrt(6p^2 - 3s)]^2`. The
[independent executable reproduction](../../code/independent_postdictions/README.md)
gives `1776.9690272931573 MeV` from the declared electron/muon calibrations;
its first-order propagated uncertainty is at most `0.000035343 MeV` for any
input covariance with the stated marginal errors.

[Belle II's 2023 measurement](https://arxiv.org/abs/2305.19116) reports
`1777.09 +/- 0.08(stat) +/- 0.11(syst) MeV`. Combining the two reported
uncertainties in quadrature gives `0.1360 MeV`; the prediction is lower by
`0.12097 MeV`, about `0.89` of that experimental uncertainty. This is a
retrospective conditional postdiction, not a new likelihood analysis or an
eligible post-freeze measurement. The balanced OPH and Koide formulations
predict exactly the same mass for every allowed input; with the same
measurement model their likelihood ratio is one. The
[existing comparison](INDEPENDENT_POSTDICTION_COMPARISON.md) details this
equivalence and the premise's observed-mass ancestry.

## Why the native clocks and C/E test are not yet distinguishing predictions

The [single-carrier clock comparison](../../model/native_clock/COMPARISON.md)
already gives a nontrivial internal spacetime correspondence. If `n` counts
global movements and `K_v` counts movements incident to port `v`, then
`sum_v K_v = 2n` on every history and `K_v/n -> 1/6` in probability. The
asymptotic count covariance distinguishes three graph-distance classes and
recovers the finite seam geometry. These are dimensionless event counts;
neither external physical timestamps nor observer federations are needed to
derive them. They do not yet identify a measured physical clock or establish
Lorentzian spacetime.

On the refining source, uniform incident-count and inverse-clock bounds
transport the actual density continuum paths to preselected local event
clocks. For fixed bounded incidence patterns, centered count fluctuations
have limiting covariance `min(t,u)*c_ij`, where `c_ij` is the shared-seam
count. The normalized correlations are `1/5` for adjacent degree-five ports
and `1/6` across an external seam. Their finite-vector asymptotic independence
from density snapshots follows by conditioning on the shared stirring history;
finite nonlinear clock–field dependence is retained. This joint source
prediction supplies no atomic transition or mapping from published servo
records to native density. The optical comparisons below test their stated
direct phase-count attachment and gain no broader verdict from this limit.

The classical histories are finite Markov processes. More generally, for any
finite instrument with Kraus maps `K_r` and `sum_r K_r^* K_r = I`, the map

```text
V |psi> = sum_r K_r |psi> tensor |r>
```

is an isometry, extends to a unitary with an ancillary system, and reproduces
the instrument when that system is measured. This applies to both executed
[C/E laws](../../model/predictions/README.md#physical-scope). Their differing
interference probabilities therefore admit ordinary quantum realizations.
This is a statement about these finite laws, not a theorem excluding every
possible OPH distinction or every ultraviolet completion.

The [charged clock-history attachment](../../model/laboratory_clock/README.md#compiled-operations-require-separate-timing-coefficients)
strengthens the internal test: one shared two-resource calibration predicts
the joint timed response without choosing separate clocks for the two
coherent settings. It preserves the complete stopping distribution and
proves that calibration workloads with fixed pulse-to-repair ratio cannot
identify the two operation durations separately. Its additional attempt
counter is a declared controller readout, and its timestamps are synthetic.
An elapsed forecast can instead sum over latent costs without disclosing
that counter. Equal-record, equal-quantum-state source histories prove that
the original readout cannot reconstruct every hidden attempt count; the
duration hypothesis is not supplied by that inference procedure.
The [protocol](NATIVE_CLOCK_EXPERIMENT_PROTOCOL.md#a-joint-timed-coherent-prediction-with-charged-operations)
states the exact conditional forecast and the physical correspondence needed
to use it. It changes no frozen prediction or empirical verdict.

### Direct event-clock attachments tested on public optical records

An event series defines a clock without external validation. Its unit can
be conventional. Identifying its readings with a particular published
measurement is a hypothesis that can be stated and tested directly.
The [optical comparison](../../model/predictions/optical_clock.py) does this
for King et al.'s [corrected Ar13+/Yb+ traces](https://zenodo.org/records/6901524).
It does not require an experimental carrier or observed microscopic repairs.

Let A and B count incident repairs at two distinct native ports, with
independent uniform seam scheduling and coincident terminal ticks included.
The [visible-clock derivation](../../model/laboratory_clock/README.md#geometry-and-fields-from-visible-incident-clocks)
gives independent interval increments with generating function

```text
g_delta(z) = (5-delta+delta*z)/(10-delta-(5-delta)*z),
E[B(N)] = N,       Var[B(N)] = N*(2-2delta/5),
```

where `delta=1` for adjacent ports and `0` otherwise. These are predictions
in A's clock units. The tested attachment identifies a reported 300-second
ratio block with `gain*B(N)/N + error`. The gain is a constant at least one;
`N` is fixed within 300 times a declared ±1 Hz bracket around
642121496772645 Hz. That bracket is based on the independent
[Yb+ frequency measurement](https://arxiv.org/abs/2010.06620),
642121496772645.10(8) Hz, whose quoted uncertainty is not a hard bound. The allowed absolute ratio error is `10^-12` per block.
Each candidate calendar block comprises exactly N consecutive A ticks,
including blocks absent from the published file. The hypothesis includes
direct phase-count readback without an additional temporal filter. These are
explicit tested assumptions, not derived atomic physics.

The [regional comparison](../../model/predictions/regional_optical.py) also
covers every nonempty port region on the same carrier. A regional clock ticks
once per event incident to any of its ports. The 4,095 regions expose 3,433
distinct event sets. If the two sets contain `a,b` seams and share `c`, then

```text
E[z^X] = (a-c+cz)/(a+b-c-(b-c)z),
E[a B(N)/(bN)] = 1,  Var[a B(N)/(bN)] = (a+b-2c)/(bN).
```

The regional attachment uses `gain*a*B(N)/(bN)+error`, retaining the same
frequency and error hypotheses. Distinct event sets have variance at least
`1/(30N)`; at equal mean rate the sharp floor is `2/(29N)`. Identical event
sets instead have `B(N)=N` on every history: they define the same clock and
are excluded from this distinct-clock test. No independent replicas are added.

Both complete public files are SHA-256 pinned and independently reduced with
exact rational arithmetic. The reducer uses complete 300-second UTC-grid
blocks, applies no value-based filter and leaves published removals as gaps.
Only compact sufficient summaries are retained; the raw-file option rechecks
the complete extraction. The measured block ranges are:

| Isotope | Published one-second rows | Complete blocks | Ratio range across blocks |
| --- | --- | --- | --- |
| 36Ar | 53,682 | 120 | `1.425293886457181e-14` |
| 40Ar | 102,888 | 210 | `1.0714598669982e-14` |

The test allows arbitrary discarded blocks and post hoc choice of the gain,
calendar subwindow and retained count within the fixed year 2021. It uses
an exact concentration bound rather than a Gaussian approximation: a
geometric-parity decomposition gives `max_k P[B(N)=k] <= 2/sqrt(N)`. If a
count interval has `w+1` integer values, its mass is at most
`p=(w+1)*2/floor(sqrt(N))`. For `K` retained draws among `T` possible blocks,
anchoring at the smallest retained value gives the range-tail bound
`T*choose(T-1,K-1)*p^(K-1)`. A further union factor `105120^3` charges every
2021 subwindow and retained-count choice. Exact integer inequalities give
adjusted conditional bounds below `10^-159` and `10^-104`, respectively.
For the regional family, each interval is an independent Bernoulli plus a
geometric count. A Bernoulli factor with variance v gives the exact Fourier
bound `max_k P[B(N)=k] <= sqrt(pi/(8Nv))`. Evaluating the range bound over
all 2,848 incidence classes and charging all 11,782,056 ordered distinct
event-set pairs gives conservative bounds `10^-121` and `10^-22`. These
rejections use the full count law, not variance alone. Searching other years,
block sizes, nonuniform schedulers or other readout mechanisms remains outside
these tests.

This rejects the **direct attachment together with its readout assumptions**.
A separate eight-replica extension assumes independent native sources,
equally spaced A ticks in the declared calendar coordinate, and a fixed
causal, nonnegative, unit-integral linear filter with compact support
`[0,30 s]`. Conditioning on all noninterior increments leaves independent
`[0,270 s]` innovations, giving bounds `10^-123`, `10^-38` for primitive
ports and `10^-123`, `10^-25` for the regional family (36Ar, 40Ar). The
[experiment](https://arxiv.org/abs/2205.13053) reports eight servo outputs
and a roughly 30-second response time constant. Neither establishes
independent native sources or compact filter support. The extension tests
its explicit readout hypothesis, not the actual atomic servo transfer.

Equating the unfiltered, zero-error population Allan variance to the sample
value gives moment targets of about `2.54e12`, `2.41e12` independent primitive
replicas, or `5.30e10`, `5.01e10` at the thirty-seam regional floor. These are
not empirical confidence bounds or inferred atom counts. No continuum level
is identified: a fixed degree-d port in the refining source has expected
`2*d*4^l` incidents per density-continuum unit, and no seconds calibration
has been fixed. The finite regional floor also depends on its thirty-seam
source. These retrospective comparisons establish neither an atomic
transition nor an OPH distinction from SM + GR.

The same-source clock/field constraint is `P=Id-diag(p)+(Sigma+p*p^T)/2`,
where p and Sigma are incident-indicator means and covariances. This reuses
the finite pair-incidence reconstruction and is checked against every native
branch. A fixed N-event IID window supplies operator uncertainty `6*delta`
with failure at most `78/(4*N*delta^2)`. It uses complete coincident clock
tapes and predicts the fair-repair density response without a second fitted
rate. These public optical ratios do not supply that joint measurement.

```sh
python -m model.predictions.optical_clock --check model/predictions/optical_clock_results.json
python -m model.predictions.regional_optical --check model/predictions/regional_optical_results.json
# After downloading both original files into a local directory:
python -m model.predictions.optical_clock --raw-dir /tmp/oph-clock-public --check model/predictions/optical_clock_results.json
```

The separate [minimal modular-transition audit](../../model/predictions/README.md#minimal-modular-transition-optical-attachment)
tests every assignment of two optical transitions to the three positive gaps
of the native three-mode generator. A shared frequency scale cancels. All
nine ordered assignments miss each published Ar13+/Yb+ absolute ratio by more
than `0.057`, under the declared total absolute ratio-error budget `1e-12`.
The retained values come from King et al., accepted manuscript, Extended Data
Table 1. This deterministic retrospective comparison treats neither quoted
standard uncertainties as hard bounds nor mismatch as a p-value. It excludes
that minimal direct-transition attachment; no atomic spectrum or universal
modular-clock exclusion follows.

The [prepared-register audit](../../model/predictions/README.md#shared-state-fluctuation-field)
separately tests the three positive gaps of its fixed `H=-log(M/m0)`, with
`p=R=1/2` and the existing four-filter cap. Exact rational log enclosures
are computed before loading observations. All nine ordered assignments
remain more than `0.057` from each published ratio, including independent
absolute per-gap attachment allowances `1e-6` and a total ratio allowance
`1e-12`. One common positive frequency scale cancels. There is no retuning
of preparation, offset, harmonic or species-dependent scale. These error
allowances are hypotheses of the restricted test; compiler guards and
quoted standard uncertainties do not establish physical hard bounds.
The [primary table](https://arxiv.org/pdf/2205.13053v2#page=20) therefore
excludes this direct constant-to-band attachment under those hypotheses,
not other transitions, composite clocks, OPH or standard quantum theory.

The growing native clock bounds and the same-product quartic execution do
not inherit the optical verdicts or the oscillator comparisons above. The
quartic model supplies its coupling, inverse temperature, radial profiles
and observer control/readout. Its occupation limit removes the finite
oscillator cutoff at fixed two-mode spatial projection, and its native
incidents drive the same interacting Hamiltonian. Neither result identifies
a measured transition or seconds per incident. The executed four-level
example has substantial occupation error; its finite detector probabilities
cannot be assigned unrestricted accuracy. A public comparison requires a
fixed apparatus/readout map and physical uncertainty model. The
[McGarry apparatus and data](https://arxiv.org/html/2603.04744v1) concern a
supplied single-mode spin–oscillator pulse experiment, not this two-mode
quartic Gibbs preparation and adaptive history.

A candidate benchmark is Aimet et al.'s published massive-to-massless
Klein–Gordon quench, with Neumann boundaries and an inferred initial
temperature. The [experiment](https://www.nature.com/articles/s41567-025-02930-9)
and its [relative-phase and density data](https://zenodo.org/records/15205139)
provide spatial snapshots at declared laboratory times. Detection is
destructive, so the time points come from repeated preparations, not
successive disturbed reads of one system. Their Neumann zero mode,
mass quench and inferred initial temperature differ from the declared
Dirichlet quartic `g=1`, `b=1/4` Gibbs menu. A comparison needs their field,
boundary, imaging and energy/time dictionary; those snapshots cannot test
the current adaptive feedback history directly. No fit or empirical verdict
for that menu is retained.

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
