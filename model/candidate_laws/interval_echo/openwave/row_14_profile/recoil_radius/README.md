# Does the heavier local formation cost predict a heavier recoil?

On the native shuttered cell the neutral four-defect stock costs more to prepare locally than the positive three-defect stock. Giving both the same mechanical impulse, the neutral centroid moves **more** readily. A body-centred material probe gives negative charge on its nearest occupied shell and positive charge on its farthest. All three signs concern the same source preparation, with fresh copies for incompatible instruments. The full OW-14 assignment of neutron and proton remains open. [Frozen measurements](PROTOCOL.md) · [Single-cell records](counts.json.gz) · [Source and earlier formation costs](../README.md)

| Material during preparation and flight | Neutral inverse inertia | Positive inverse inertia |
| --- | ---: | ---: |
| Native with private references | 0.04518 ± 0.00091 | 0.02679 ± 0.00088 |
| Native with private references shuttered | 0.06980 ± 0.00118 | 0.13922 ± 0.00094 |
| Exchange alone | 0.16752 ± 0.00117 | 0.25110 ± 0.00093 |

These are local zero-impulse current responses in exchange-energy times squared cell-radius units, bounded with impulses of magnitude one quarter. The neutral-to-positive effective inertia ratio in the native cohort is below 0.626. The same comparison after closing the private references gives a neutral-to-positive ratio above 1.94. Mechanical inference here uses the trace-average curvature in three axes, with a common action unit. The later motion of a freely propagating body requires further readings.

## Two independent ways to read the response

Each of the 12 target ports carries load 0, 1 or 2. Against a filled host those give signed defects −1, 0 or +1. Shutter the same data cell used for the formation-energy measurements. Prepare as before from the uniform mixture of all words at fixed inventory and accept the two physical clock reads near the sector's registered ground. The impulse uses the centroid X of all unsigned defects along a calibrated Cartesian axis, dividing by the number D of defects. A local coherent ternary read controls a phase proportional to the unsigned defect indicator, followed by its inverse read. Multiplying the phases gives exp(iuX), with total impulse u independent of D. The coordinates, shuttered contacts, source clock, reversible detector controls and choice of equal defect weights are supplied.

One observer reads the odd exchange parity on a uniformly chosen seam. Its score is 30 on the odd exit or zero on the even exit. Subtracting the no-impulse bank from the mean energy of the two opposite kicks cancels every diagonal term. At half a unit of impulse the resulting curvatures are 0.04408 ± 0.01252 and 0.02734 ± 0.00768 for the neutral and positive stocks. This meter agrees with the current readings within its wider intervals. The published intervals include the finite-step allowance.

On independent copies another observer measures the centroid current i[H,X]. If a seam exchanges words i and j, their positions differ by d. Read the unordered two-site occupation and the other sites, preserving the coherence between the two orders. For a nonzero d, add a relative quarter phase with the sign set by which order has larger X, then read the same swap parity. The two quadrature exits give opposite scores of magnitude 15|d|. The factor 15 comes from 30 seams and the exchange amplitude −1/2. At d=0 the score is zero on either exit. Opposite kicks cancel a constant offset. On a single-axis source with real density matrix, the velocity equals one half of the sum over all words i and seams of rho(i,Si) d sin(ud), where S is the chosen swap. A direct two-word commutator test checks the sign and normalization. At initial release the same inner-cell current is obtained if the exterior seams are opened, since the filled exterior has no coherence with a word in which a defect has already left. The measurements do not demonstrate a later interval of free rigid motion.

For the near-zero curvature, one half of the same sum with d² replaces d sin(ud)/u. With B the largest row sum of d⁴ over the 30 seams, the current difference divided by twice the impulse has curvature bias at most u²B/12. The symmetric energy difference divided by u² has bias at most u²B/24. These bounds use twice the absolute off-diagonal density entry being at most the sum of its two diagonal entries. Each reported statistical radius applies Hoeffding's inequality to the full score range, conditional on the observed number of accepted starts, with a family of 4,096 elementary comparisons at 99% coverage. Both original refusals and every quadrature or parity exit enter the archive. Tracing over three Cartesian axes removes the choice of their orientation at zero impulse. [Acquisition and decoder](measure.py) · [Independent checks](test_radius.py)

The exchange-only ground provides an independent normalization. Its equal-amplitude word on the connected permutation sector puts precisely one defect on a specified seam with probability 2D(12−D)/(12·11). That swap changes the centroid by one edge vector divided by D. Every edge has the same squared length. The ratio of positive to neutral inverse inertias is therefore 3/2, as the exchange-only meters confirm. Counting three versus four equally weighted defects without the exclusion of occupied destinations would give 4/3.

## What the sign probe sees from the centre of the material

First read the unsigned occupied sites. Locate their Euclidean centroid using the unit-radius icosahedral chart. A uniformly selected port enters the balanced material interferometer with probability proportional to its squared distance, or with probability one if it belongs to the enrolled occupied shell. All equal-distance ties are included. The matched bypass, two opposite charge calibrations and occupied unit-host control give the signed probe gain. The angular interference measurement uses its ordinary and inverted analysis on separate copies. In a further neutral control, the terminal read is performed on the sign-conjugate word. This can be implemented after a sharp word read by swapping each recorded excess with a recorded vacancy, preserving the support and centroid. The shape probe is diagonal in the recorded words. [Balanced physical probe](../scattering/README.md)

| Neutral read on native admissions | Forward signed score | Inverted signed score |
| --- | ---: | ---: |
| Squared material radius | +0.00290 ± 0.00163 | −0.00315 ± 0.00163 |
| Nearest occupied shell | −0.00401 ± 0.00062 | +0.00395 ± 0.00062 |
| Farthest occupied shell | +0.00407 ± 0.00062 | −0.00404 ± 0.00062 |

The first score has material-charge times squared cell-radius units, the two shell scores material-charge units. Sign conjugation changes the nearest and farthest forward means to +0.00392 and −0.00405. The filled host gives zero within the quoted precision. Shuttering the private references gives a radius consistent with zero. On the connected-only native cohort the near and far uncertainties both exceed 0.0127 and neither sign is resolved. The exchange-only analyzer has very small sign gain, as its independent standards record. A fixed gateway and a centre obtained from the material ask different experimental questions.

There is a geometry check on every neutral word, also valid for any four positions in Euclidean space. Twice the signed second moment about the four-defect centroid equals the squared distance between its two excesses minus the squared distance between its two vacancies. The local pointer consequently tests which pair is wider. The identity is proved coordinate by coordinate and checked for all 2,970 neutral words. [Lean proof](RecoilRadius.lean)

## Where does the compensating vacancy go?

There are exactly two excesses and two vacancies in each two-cell trial. One endpoint holds all four defects in the data cell, with its store filled. The other has a positive three-defect data cell and a single vacancy in its store. They occupy the same pair of independently shuttered cells as the earlier work comparison. Both stocks are clock-selected. The separate storage-only current calibration gives 0.75109 ± 0.00759 on the native material. For each held joint launch the receiver keeps both admission bits, chooses a cell uniformly and reads one of its 30 seam quadratures.

| Native joint endpoint and kick | Quarter-impulse inverse inertia |
| --- | ---: |
| Neutral data, filled store | 0.04515 ± 0.00051 |
| Positive data alone, single vacancy left at rest | 0.02681 ± 0.00150 |
| Positive data and storage vacancy, equal impulse per defect | 0.06237 ± 0.00105 |

For the last row the data centroid receives 3/4 of the impulse, the storage vacancy 1/4. Their independent calibrations forecast the joint curvature by adding 9/16 of the data response and 1/16 of the vacancy response. The homogeneous kick to all four defects yields a heavier neutral endpoint. Leaving the storage vacancy at rest reverses the inference. On all 18 held size-of-kick, endpoint and generator comparisons, the independent one-cell forecasts agree with the two-cell records within the joint coverage bounds. On a positive data cell and a store of inverse responses p and v, the least quadratic work at fixed total impulse gives inverse response pv/(p+v), attained by sending the fraction v/(p+v) to the data. This optimized response is at most p and thus also below the native neutral response. [Frozen comparison](STORAGE.md) · [Joint score and admission records](storage.json.gz) · [Local readers](storage.py) · [Square completion](RecoilRadius.lean)

The body-centred charge and controlled kicks specify why the native prepared stocks fail as an OW-14 assignment. A free particle family with an electric probe, a common physical length scale and a storage-independent rest-energy or inertial comparison remains outstanding. OpenWave lists M10 as untested on this row. [Official comparison](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#oph-m10)

From the RER root, replay the recorded banks with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_14_profile.recoil_radius
.venv/bin/python -B -m "$p.measure"
.venv/bin/python -B -m "$p.storage"
.venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_14_profile/recoil_radius
~~~
