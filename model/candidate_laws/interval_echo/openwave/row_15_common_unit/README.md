# Can a common unit reproduce both nucleon masses?

On the clock-selected neutral and positive material stocks, a single conversion from the two local formation energies to MeV misses at least one of the published nucleon masses by **8.72%**. Including the positive stock's balancing vacancy lowers the minimum discrepancy to **3.29%**. A different zero-impulse detector reads the local mechanical response of each preparation. Its two inverse curvatures also fail the common-unit test. **OW-15 remains open.** OpenWave lists M10 as untested on this row. [Criterion and M10 table](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#oph-m10) · [Frozen assignment and control](PROTOCOL.md) · [Original simulated count banks](counts.json.gz)

The positive stock contains two excesses and one vacancy in the filled cell. The neutral stock has two of each. The same previous double-clock preparation accepts the lowest resolved line of each inventory. Resetting the energy zero to the ordinary filled host and retaining every clock refusal makes the two data energies comparable on this fixture. A second matched cell is either filled or carries the balancing vacancy. [Source preparation and direct-kick measurements](../row_14_profile/recoil_radius/README.md)

| Native assignment | Neutral / positive, simultaneous 99% interval | Lower bound on the larger fractional mass error |
| --- | ---: | ---: |
| Local formation energies | 1.19292–1.19397 | 8.72% |
| Formation energies including the balancing store | 0.93675–0.93755 | 3.29% |
| Local impulse inertias | 0.59331–0.59446 | 25.49% |
| Fixed inventory with the same impulse per defect | 1.38017–1.38163 | 15.90% |
| Fixed inventory with minimum quadratic work | 0.57300–0.57409 | 27.12% |

The upstream values are 938.3 and 939.6 MeV, a ratio of 1.00138548. Anchoring the positive local formation reading to 938.3 MeV would give a neutral reading between 1,119.31 and 1,120.30 MeV. Re-anchoring the complete positive endpoint gives 878.95–879.70 MeV. These illustrative conversions use the proton mass as an input. A rest-energy prediction needs independently identified baryons and a physical energy scale.

## A zero-impulse measurement

Choose one of the 30 seams uniformly. Read whether exactly one endpoint contains a defect relative to the filled host. This request is invariant under exchanging the endpoints. Its outcome can be copied by a coherent ternary-occupation comparison followed by uncomputation of the ordered inputs. On a mobile pair, read the even or odd swap parity. A blocked pair gives zero, even gives +1 and odd −1. In a preparation with D defects, let X be the Cartesian projection of their unsigned centroid in the unit-cell-radius chart. For a swap S the second derivative of the energy of exp(iuX) at zero impulse is one half the sum of rho(i,Si) times the squared change of X. Tracing over three axes gives the squared seam length divided by 3D² on a mobile pair. All other pairs carry zero. Thus the averaged inverse inertia equals the mean detector score times 30 times the squared seam length divided by 6D². The three outcomes cover the original preparation on every admitted launch. The sampled native curvatures are in [0.045141, 0.045173] for the neutral stock and [0.026801, 0.026835] for the positive stock, in exchange energy times squared cell-radius units. They agree with both finite-impulse strengths of the independent quadrant-current experiment.

On fresh copies a different observer reads the ordinary energy of a random internal seam or a private reference. A seam allows simultaneous contrast and parity measurements, while a private reference reads squared occupation relative to the filled host. Multiplication by 39 restores the total energy. Each detector bank uses 2^44 original starts, with the source's clock calibrated independently of the held terminal parity and energy panels. Hoeffding bounds for the full score range cover the 256 budgeted elementary comparisons simultaneously with at least 99% probability conditional on the admitted counts. The complete receipt retains all three stocks, both endpoint comparisons and the filled-host zero standard. [Readout, mean and interval calculation](measure.py) · [Direct commutator and accounting checks](test_ratio.py)

For a measured mass ratio r and target ratio t, let a and b be the two scaled outputs divided by their respective targets. Their quotient is r/t. If each lies within fractional error e of one, then |b−a| ≤ e(a+b). Equality is attained by choosing the common scale so that a+b=2. Hence even an optimally chosen unit has error at least |r−t|/(r+t) on one of the two masses. The table takes the point closest to t in each confidence interval and rounds its lower bound down. A separate additive offset or a charge chemical potential would introduce another fitted constant and change the filled-host reference. [Lean scale inequality](CommonScale.lean)

## Does an adjustable impulse share rescue the assignment?

If the fraction f of total impulse goes to the positive data, the remaining fraction goes to the stored vacancy. For independent cells the inverse response is p f² + v(1−f)², with p and v determined by the two local current readers. There is one tuning on the native material capable of fitting the mass ratio. All possible solutions lie in **0.80858 ≤ f ≤ 0.80881**. At the central positive estimates the curvature crosses the target between f=3/4 and f=1. This is a continuous calibration using the requested masses. The material that receives each part of the impulse is observable.

When the contrast is halved, the minimum over *all* f gives a mass ratio above **1.18864**, excluding the target even if the impulse split is retuned. The fixed-inventory formation ratio on this control is 0.93432–0.93540. To bound every possible split, cover [0,1] with 100,000 closed bins and bound both positive terms of the quadratic separately on the whole bin using the simultaneous detector intervals. The two reference hypotheses, separately prepared baryon rest energies and free inertial masses, remain open beyond these specified stocks. [Prespecified split check](MOMENTUM.md) · [Full-interval decoder](momentum.py)

For a single-thread replay from the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_15_common_unit
.venv/bin/python -B -m "$p.measure"
.venv/bin/python -B -m "$p.momentum"
.venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_15_common_unit
~~~
