# What survives an angularly uniform launch?

An observer can predict the return of a freely moving load after partially changing its phase at a measured radius. The same observer loses most of the fine-board return contrast by assuming that the load equilibrates over all addresses at that radius between flights. These results apply to a supplied stationary seam mark and a freely moving load under the native Interval Echo generator.

## Resetting the angle changes the flight

Prepare the marked gateway by a sharp odd-coin word with recorded energy 32 on a guarded seam. Calibration in the native 99-state local model gives zero probe-dependent contact shift for unit and double mobile loads. In the matched guarded preparation both mobile sectors propagate with the same Laplacian divided by two after subtraction of a charge-dependent scalar. The observable is active-seam hop distance from this marker. Choose one of two radial windows A and sample uniformly over all its port addresses. Read a fresh one-leg ensemble at time t and an independent uninterrupted ensemble at 2t. Separate copies receive a sharp window observation at t, or that observation followed by uniform re-preparation in the measured window. [Acquisition protocol](PROTOCOL.md) · [Native propagator and instruments](measure.py)

| Cells | Width in hops | t | Window | Launch addresses | Free return | Observed at midpoint | Uniformly re-prepared | Reset error / free departure |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 80 | 1 | 0.5 | 3 | 9 | 0.540893 | 0.742225 | 0.717100 | 38.38% |
| 80 | 1 | 0.5 | 6 | 17 | 0.516547 | 0.728903 | 0.701369 | 38.23% |
| 1,280 | 2 | 1 | 3 | 35 | 0.648467 | 0.724895 | 0.649849 | 0.39% |
| 1,280 | 2 | 1 | 6 | 87 | 0.673924 | 0.735307 | 0.638443 | 10.88% |
| 20,480 | 4 | 2 | 3 | 146 | 0.725619 | 0.741816 | 0.619027 | 38.85% |
| 20,480 | 4 | 2 | 6 | 326 | 0.737614 | 0.751548 | 0.635651 | 38.86% |
| 327,680 | 8 | 4 | 3 | 614 | 0.832385 | 0.835820 | 0.766509 | 39.30% |
| 327,680 | 8 | 4 | 6 | 1,334 | 0.828294 | 0.832279 | 0.761535 | 38.88% |

The values are native ensemble probabilities. The archive also contains independent synthetic readout counts with 536,870,912 trials per cohort. For window multiplicities d and measured one-leg transition probabilities K, microscopic reciprocity gives the uniformly re-prepared return as dₐ Σᵦ Kₐᵦ²/dᵦ. Direct propagation on the full 880-port mobile sector of the 80-cell host reproduces this prediction by summing transitions measured from each possible intermediate launch window. An independent reversed launch from window 6 to 3 gives 12,682 returns. The forward counts and surveyed multiplicities predict 12,740. [Probabilities, sampled counts and full-board validation](results.json)

Both finest windows fail the preregistered accuracy requirements for the freely propagating load. Sharp midpoint measurement changes their return by 0.003435 and 0.003985. Uniform re-preparation changes the measured-arm return by a further 0.069311 and 0.070744 in the opposite direction. Thus the fine-board discrepancy mainly comes from the conditional state inside the reported window, which the uniform reset discards. The table gives finite-scale evidence and does not establish an asymptotic law.

## A radial phase that the observer can read later

Apply a phase exp(iθ) at the midpoint only to the launch window. This externally supplied intervention leaves the address populations at the midpoint identical for all θ. The subsequent native flight reads their coherence. If X and Y are the two return amplitudes through the launch window and its complement, the probability is the squared norm of exp(iθ)X+Y divided by the number of launch addresses. Three independent settings at 0, π/2 and π determine its constant, cosine and sine terms. Those measured coefficients predict fresh counts at −π/2 and π/3 without a fitted radius, energy or exposure. [Phase protocol and enrollment order](PHASE_PROTOCOL.md)

All 32 held predictions pass the simultaneous 99% bounds over both loads, both windows and four hosts. The largest held frequency error is 0.00006341. The bounds are 0.00037980 and 0.00022465 for the two held settings. The largest change among the five phase settings is 0.006870 in the inner window and 0.007970 in the outer window on the finest host. The full-board check applies the phase and second propagation independently in the address basis. The finer-board computation uses the symmetric one-leg propagator and agrees with separately propagated uninterrupted two-leg amplitudes at zero phase.

These are held predictions for phase-controlled return at one pair of exposures on each refinement. They identify a radial interference observable. Selecting a proper radial measure, predicting general unperturbed flights and reconstructing causal geometry require further evidence. The phases, guard pattern and choice of window belong to the preparation and readout. [Finest-board and control records](results.json)

## Reproducing the acquisition

Run from the checkout using the four pinned federation deposits consumed by the shared refinement loader. The output path must be fresh.

```sh
python -m model.experiments.interval_echo_einstein_station.relational_ruler.angular_equilibrium.measure \
  --directory /path/to/deposited/carriers --output /tmp/radial-flight.json
```

Both flight times use Chebyshev polynomials of the native bounded Hamiltonian. A source window's maximum occupied radius plus the polynomial degree bounds every visited port. The stencil retains the full degree at its boundary. The Bessel tail bounds each exponential in operator norm below 10⁻¹³, and the second-leg return control tests agreement to 5×10⁻¹². Independent random seeds are recorded in the two protocols.
