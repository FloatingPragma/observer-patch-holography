# Does the correlated native collar change the release focusing of four tracers?

An ordinary load changes the initial focusing of a native four-tracer tetrahedron by **−27.65519** per squared exposure. At fixed ordinary energy, switching from fresh to shared preparation tags gives **−0.00395 ± 0.06467** in the same instrument. The release response to the shared history is exactly zero under the registered generator. At exposure two the joint volume does resolve which collar ensemble was prepared. [Frozen protocol](PROTOCOL.md) · [Complete observer histograms](results.json)

| Change of preparation | Exact change in release focusing | Measured change in release focusing | Change in measured signed volume at exposure two, divided by release volume |
|---|---:|---:|---:|
| Vacant to ordinary load with fresh tags | −27.65519 | −27.65797 | +0.72389 |
| Fresh to reversed tags, matched ordinary energies | 0 | −0.00100 | +0.000005776 |
| Fresh to shared tags, matched ordinary energies | 0 | −0.00395 | −0.000005833 |

The three focusing contrasts have joint 99% half-width 0.06467 conditional on the independently measured clock. The finite-stencil bias has been removed using the archived word dynamics. This changes the shared-to-fresh estimate by less than $10^{-12}$. The absolute shared-history focusing response is bounded by **0.249%** of the magnitude of the ordinary response on this tetrahedron under the same simultaneous intervals. The delayed volume contrasts have joint 99% half-width 0.000001756 across the preregistered times and arms. The two-exposure shared and reversed forecasts are −0.000005875 and +0.000005875.

## What the volume reader measures

The service and source laws come from the [native two-collar experiment](../shared_variance.py). A sharp word survives all four flights. Four load tracers are launched sequentially with resetting between them. Conditional on that word their terminal measurements are independent. Each tracer receives a balanced, zero-current launch on a blank edge. The tetrahedron uses rails (4,5), (0,7), (4,11), and (8,9). Its shape was chosen by maximal release volume among eligible layouts including the archived ordinary-sensitive launch. Every final home port has coordinates given by the first nonconstant native triplet; an outside event has coordinate zero. The same rule assigns a signed volume to every four-port record, including outside events. No probability is conditioned on return. The estimated populations use independent simulated acquisitions at each appointment and 2^48 four-port draws per arm and time.

Define the focusing as minus the second material-clock derivative of the mean signed tetrahedral volume divided by its release value, 0.11470 in native angular-coordinate units cubed. The angular coordinates are the calibrated finite home-cell triplet. A spacetime Ricci interpretation would additionally require proper separation vectors and freely falling geodesic deviation. The Bianchi I scalar radar comparator has such an explicitly supplied geometric model; the native angular tetrahedron uses the material service specified above.

## Why delayed distortion and initial focusing differ

For a protected word $w$ the flying load evolves under $H_w=H_0+D_w$, with a diagonal potential linear in its four contact energies. For each diagonal coordinate $X$, the first derivative vanishes on a real launch. Its second derivative is the expectation of $-[H_0,[H_0,X]]-[D_w,[H_0,X]]$, which is affine in the four contact energies. Their means agree under fresh, shared and reversed tags.

The oriented volume is affine in each of its four vertices. Conditional independence makes its expectation equal to the volume of the four conditional mean positions. At release these positions are common to all words and all their first derivatives vanish. Differentiating twice leaves a sum of four linear acceleration contributions with common geometric coefficients. Averaging gives the same initial volume acceleration in all three correlated-source arms. At finite time the four conditional means have evolved, and products of their word-dependent displacements can distinguish the histories. The held terminal observations detect that effect with the source's complete two-word ordinary-energy histogram held fixed.

The [dark-sector source definition](../../../../cosmology/oph_dark_matter_paper.tex) assigns an anomalous local stress to a fourth-order modular energy variation on independently calibrated shrinking proper balls. This tetrahedron measures one angular frame and one source footprint. The repeated exposures change duration. The fresh, shared and reversed arms change history covariance. Neither operation calibrates a decreasing proper radius. The native result constrains the proposed correction to this focusing instrument while leaving the physical fourth-order collar coefficient undetermined.

Run the following commands from the repository root.

```sh
.venv/bin/python -m model.experiments.interval_echo_einstein_station.dark_curvature.measure --output /tmp/native-dark-curvature.json
.venv/bin/python -m pytest -q model/experiments/interval_echo_einstein_station/dark_curvature
```
