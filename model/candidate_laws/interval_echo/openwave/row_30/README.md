# OpenWave 30 · Klein–Gordon dispersion of a guided matter wave

**Demonstrated on a driven one-dimensional guide.** A doubled load propagates in both directions along a closed path on Interval Echo’s spherical federation. The controller alternates two sets of contacts, and the difference in their durations supplies the gap. The frequencies inferred from local detector counts agree with the low-energy Klein–Gordon prediction within **1.05%** at all 48 held-out signed momenta across three screens and two gap settings. The largest error using exact probabilities is **0.58%**. Under unrestricted free evolution, the same doubled load gives a different dispersion for squared graph wave number. [Guided protocol](PROTOCOL.md) · [Counts and probabilities](results.json) · [Free control](../row_30_free/README.md)

## What the observer measures

The 20-, 80- and 320-cell screens supplied closed paths of 110, 214 and 246 actual ports, respectively. They cross 15, 36 and 41 real intercell seams. One spatial unit is the translation between consecutive two-rail units along the path. A sharp load on the gapless guide moved six spatial units in six periods with unit probability. This measures the velocity used in the predictions. The spatial scale is distance counted along the guide; independent physical ranging of the individual seams remains open. The choice of a waveguide, the two alternating sets of gates, their durations, the launch state and the readout windows belong to the apparatus.

For the independent momentum reading a patch mixes the A-rail amplitudes from neighbouring units on a native half contact. The two local click probabilities are proportional to 1−sin k and 1+sin k. The reader retains both outcomes and quiet. We use the prespecified principal branch |k|<π/2. On separate trials a single B-rail reader measures the beating of the two bands. The fitted frequency ω is half the interband beat in the low-quasienergy branch. A zero-winding source measures the rest frequency m for each installed detuning. None of the held frequencies enters the speed or gap calibration. Finite tables use 33,554,432 simulated preparations per read setting. The high exposures are aggregated binomial or multinomial draws. [Native trajectory and reader implementation](experiment.py)

| Cells | Gap setting δ | Largest held frequency error from exact probabilities | Largest held error from counts | Largest error of the tangent parabola |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.1 | 0.162% | 0.708% | 145.0% |
| 20 | 0.2 | 0.574% | 0.520% | 45.6% |
| 80 | 0.1 | 0.162% | 0.980% | 151.0% |
| 80 | 0.2 | 0.579% | 0.657% | 48.1% |
| 320 | 0.1 | 0.162% | 1.046% | 146.3% |
| 320 | 0.2 | 0.575% | 0.805% | 46.1% |

The KG forecast is ω = √(m²+c²k²). The competing parabola ω = m+c²k²/(2m) uses the same independently measured rest frequency and speed. At the larger detuning δ=0.8 and |k| near 1.15 the KG forecast misses by **9.20–9.44%** on exact local probabilities. On the 320-cell guide two finite-width packets moved at **−0.8481 and +0.8481** spatial units per period. The KG group-velocity forecasts at their central momenta were **−0.8728 and +0.8728**. Every position contributes to those packet centroids.

## Timing by local positive repairs

A five-cell observer on the 320-cell screen uses its ordinary stationary repair source and four off-guide clock cells. It opens its B-rail read when its own physical count first reaches each requested threshold. A material witness reached the 32nd tick at period 59. For the simulated ensemble, the waiting-time distribution is determined by the native one-half positive-repair probability. The fit incorporates this clock noise and preserves quiet outcomes. At the two held signed momenta it gives frequency residuals of **0.89% and 0.48%**, with rest frequency **0.40029 radians per mean home tick**. Frequency in these units refers to the band splitting inferred through the waiting-time response. Individual clock histories fluctuate. [Clock-triggered protocol](CLOCK_PROTOCOL.md) · [Detector and physical hand records](patch_time.json) · [Clock-triggered replay](patch_time.py)

## Why this branch passes

On a blank seam the pair of basis states with one mobile double has generator (4κ+½)I−X/2, where X exchanges the rails. The program reduces the registered 99-state native operator to this invariant sector and evolves real port amplitudes with its matrix exponential. Empty contacts preserve the vacuum. After removing a common phase, a period consisting of contacts of durations π and π−2δ has two reciprocal eigenphases. Their half-splitting in the low branch satisfies

$$\cos\omega=\cos\delta\cos k.$$

The determinant-one two-band step therefore gives the following real-space identity for either rail's envelope after the same common phase has been removed:

$$\Phi_{t+2,j}+\Phi_{t,j}
   =\cos\delta\,(\Phi_{t+1,j+1}+\Phi_{t+1,j-1}).$$

The native trajectory of a **localized broadband excitation** verifies this recurrence at fifteen consecutive triples of time steps on every guide and detuning. In the joint scaling δ=μh, k=ph, ω=Ωh the identity gives Ω²=p²+μ²+O(h²) at fixed p and μ. This is the (1+1)-dimensional Klein–Gordon envelope. At fixed δ and arbitrary microscopic k the cosine law is the exact prediction.

The full A/B band gap is 2δ in radians per service period. The controller creates it by detuning the exchange, and the local beat measures it. The κ-dependent native contrast term supplies a common dynamical phase in this sector on each perfect matching. The κ=0 exchange-only control gives the same propagation modulo that phase on all nine fixtures. These controls identify the mechanism and the apparatus dependence. Rest frequency here is a Floquet half-gap in the supplied drive frame; an absolute inertial mass and an independently relaxed, unguided vacuum field have yet to be identified. [Underlying repair interaction](../../../contrast_memory_spherical_v1/SPEC.md)

For comparison the unfiltered one-load and mobile-double generators on all three complete screens give ω=b+λ/2, with λ the independently specified graph-Laplacian eigenvalue. At this energy zero b is 48 for the single load and 192 for the double. Squaring leaves λ²/4 after the constant and linear terms. No choice of a constant energy reference removes its quadratic coefficient. This rules out exact KG for those autonomous branches when k²=λ. Other excitations and spatial momentum observables require their own tests. [Registered generator](../../candidate.py) · [Unfiltered controls](results.json)

## Reproduction

From the RER repository root, with its model dependencies installed, run:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B -m model.candidate_laws.interval_echo.openwave.row_30.experiment
python -B -m model.candidate_laws.interval_echo.openwave.row_30.patch_time
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30
```
