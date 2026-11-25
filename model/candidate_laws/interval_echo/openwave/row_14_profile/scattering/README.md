# Scattering from a neutral material core

A balanced material probe resolves the two signs inside the guarded four-defect source. When the four defects form a connected cluster, the probe scatters elastically at their intrinsic center with probability **0.51571**. Across all admitted sources that probability is **0.0044960**. A phase-sensitive reading on a separate cohort finds nearly equal contributions from excesses and vacancies in the connected core. A small mean signed signal can therefore accompany a large scattering probability. [Method and proof](#how-the-probe-couples-to-the-source) · [Complete simulated banks](results.json) · [Count-only decoder](readout.py)

| Neutral source and aperture | Elastic match per admitted source | Mean signed material load in aperture |
| --- | ---: | ---: |
| Entire cell, coherent paths | 0 | 0 |
| Intrinsic center | 0.0044960 | −0.0015864 |
| Intrinsic outer six sites | 0.037636 | +0.0040534 |
| Gateway 1 center | 0.33515 | +0.33256 |
| Gateway 1 outer six sites | 0.060132 | −0.097729 |
| Intrinsic center, connected sources only | 0.51571 | −0.000382 |

The table reports the computed law. Every detector bank enrolled 2⁴⁰ fresh simulated starts. The neutral preparation admits 0.0588454% through the two source filters. On independent cohorts, the intrinsic-center signed interval is [−0.001786, −0.001396] and the intrinsic outer interval is [+0.002816, +0.005153]. Gateways 3 and 10 have the same computed profile as gateway 1. Closing the private references gives zero signed mean on every neutral aperture by sign symmetry, with intrinsic-center elastic probability 0.0071093. Thus these two readings can distinguish a sign bias due to the landmarks from an unpolarized population of material defects. [Preparation and direct port evidence](../README.md)

## How the probe couples to the source

Owner zero filters the original neutral or positive preparation. After success it holds the source and makes a binary host/nonhost read at each port for the intrinsic-center rule and the connectivity flag. A uniform tie is retained. The scattering receiver uses the enrolled aperture and obtains the sign from its own two-rail probe. Its controller supplies coherent routing and path analysis. Each address branch swaps the addressed material along a path inside the owner to gateway 1, contacts a neighboring exterior probe rail, and reverses the swaps. The second probe rail is an adjacent exterior spectator. Full-duration blank swaps at strengths 8 and 0 have unit transfer amplitude for each of the nine load words. The enumerated paths therefore restore every transited rail. The receipt checks the shared cell and feasible probe contacts on 20, 80 and 320 owners. The spectral source calculation uses the twenty-owner guarded preparation. [Program](measure.py) · [Frozen choices and comparison](PROTOCOL.md)

Write the active and spectator stock in the order of their two loads. The probe starts in the normalized sum of (2,0) and (0,2). It is accepted in their normalized difference. A sharp test load identical to the active probe load is unchanged. For distinct loads at difference d, the amplitude for both to stay is exp(−iκtd²) exp(−it/2) cos(t/2). This follows by diagonalizing the swap on the two orders of the blank native seam. An exchanged unequal load misses the balanced output because the spectator is untouched. The difference of the two retained amplitudes is therefore gq, with q equal to the tested load minus one and

$$
g=\tfrac12[1-e^{-i(4\kappa+1/2)t}\cos(t/2)],\qquad t=\pi/32 .
$$

Uniformly addressing A sites and projecting back to the same address state gives the diagonal amplitude K = gQ/A, with Q their total signed load. All material words survive in this derivation. Their off-diagonal entries contribute zero to these effects even if the filtered source is coherent. The binary support read commutes with them. An elastic detector clicks with probability |K|². In the independent interference bank an equally weighted bypass prepares the same address and the odd analyzer state directly. The two path ports have probabilities |1±K|²/4 and the complementary outcome has probability (1−|K|²)/2. Single-run records contain a detector port, a source decision and the connectivity flag. The assigned graph lengths are counts of cell links. The contact duration is an apparatus setting relative to exchange. No inferred energy or common spacetime clock enters the decoder.

The difference of the two path probabilities gives Re(g) times the mean Q/A. Sharp standards of load 0, 1 and 2 calibrate its sign and gain. The summed path probabilities on the positive standard also give |g|². At a single site the first moment C and second moment S recover the two populations as (S+C)/2 and (S−C)/2. These formulas use only local probe counts and calibration counts. On the connected neutral subsample they put the intrinsic excess fraction in **[25.51%, 26.14%]** and the vacancy fraction in **[25.55%, 26.17%]** with simultaneous coverage of at least 99% (outward rounding). On all admitted neutral sources the fractions are [0.131%, 0.161%] and [0.290%, 0.320%]. The separate full 99-state seam check evolves a coherent two-site source through both address branches and verifies the elastic and interference reductions. [Tests](test_scattering.py)

## Forward cancellation and its controls

The neutral forward channel gives zero of 647,020,015 admitted simulated starts. Recording the probed address before returning to the uniform mode removes its interference with the other addresses. The matched probability becomes |g|² times the sum of squared site charges divided by 144. The same neutral inventory gives **0.0277276**, with 17,950,459 simulated matches. The coherent positive source gives 0.00693190 and 13,188,719 matches. Whole-cell cancellation is an inventory identity for the coherent signed instrument. The relative weights at specified apertures and the large increase at a connected center depend on the filtered source. At the same connected center the exchange-only probe gives 0.00031097 at the same duration. The unit standard and the undefected unit-filled host both have zero elastic amplitude.

The source shutters, prepared inventory, contact standard, absolute phase reference and spatial addressing are apparatus inputs. This measurement resolves material charge in a guarded cell. Electric form factors of a freely propagating neutron and the proton–neutron rest-energy ordering require electromagnetic and inertial calibrations on identifiable baryon states. The local full OW-14 criterion is Open. [Official wording and M10 cell](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L449)

## Reproduction

The replay regenerates the spectral source, scattering response and every multinomial count. The count decoder reads the archived detector outcomes. It evaluates bounded-score Hoeffding intervals conditional on the admitted count. A multiplicative binomial bound supplies a simultaneous lower acceptance for propagation of the fixed-grid probability error. The budget covers up to 1,000 tail or two-sided interval statements with total failure probability at most 1%. These are sampling intervals for the simulated apparatus. The source diagonalization uses floating-point linear algebra. Finite sampling bounds do not enclose that separate numerical error.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_14_profile.scattering
nice -n 15 .venv/bin/python -B -O -W error -m "$p.measure"
.venv/bin/python -B -m "$p.readout"
.venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_14_profile/scattering
~~~
