# What wavelength do observers see on the selected charged line?

Four groups of sharp counters read the vacancy along the bare interior of the three arms. The counters pool the three ports at equal arm distance and discard all six terminal memory bits. The finite coherent clock selects the previously measured 16.43120649 line. Its control is supplied. The spatial read uses only the four charge counts and the number of original starts. The stations were chosen on the simulated profile before the fresh count bank. [Physical ports and raw count bank](wavelength.json) · [Clock and pole tests](results.json)

After the one-third input aperture, the four groups counted 424,947, 171,706, 1,260,333 and 745,309 detections in 274,877,906,944 simulated launches. All other outcomes enter the fifth bin. The inferred standing-wave phase per step is 0.9627–0.9688 radians or 2.1728–2.1789 radians at simultaneous 99% confidence. Position intensity alone has this spatial alias. Both intervals exclude zero. This establishes a finite wavelength inside the material arms. Center-of-mass momentum in a reconstructed three-dimensional space and the rest energy of the composite require other measurements.

## Derivation

At the measured stations every nonzero native column has diagonal 17 and hopping coefficients −1/2 to its two neighbors with identical head records. For an exact energy E the amplitude of each word obeys

    u[j+1] + u[j−1] = 2 cos(k) u[j],    cos(k) = 17 − E.

The density at station j sums the squared magnitudes over the arms and all memory records. Each summand has the form A + B cos(2kj) + C sin(2kj). Its first differences d[j] = n[j+1] − n[j] therefore satisfy d[j+1] + d[j−1] = 2 cos(2k) d[j]. This also holds for mixtures at one spatial frequency. At a band endpoint the amplitude is affine up to an alternating sign. The density is then quadratic and its third difference is zero.

The readout divides the four-counter contrast with coefficients (−1,1,−1,1) by the contrast (0,−2,2,0). It propagates simultaneous Clopper–Pearson intervals for all four binomial marginals of one multinomial draw. A separate band-edge test bounds the third difference with coefficients (−1,3,−3,1) between −1.079×10⁻⁵ and −1.064×10⁻⁵ per launch. These intervals describe the finite-filter count measurement. Independently, inserting the isolated spectral interval into the exact native recurrence selects the lower phase interval. [Spatial reader and native column checks](arm_wavelength.py) · [Spectral interval](results.json)

Run from the repository root with one numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.arm_wavelength
```

## Predicting a smooth profile from four counters

Can the observers use a smooth wave to predict intensities at places they have never measured? Four existing counter banks at ordinals 3 through 6 determine an interval for the spatial frequency. No energy, eigenvector or simulation amplitude enters this inference. At a single energy each head component of the interior has a smooth interpolation solving the one-dimensional Helmholtz equation. Its pooled intensity at distance x has the form A + B cos(2kx) + C sin(2kx). Positions use one physical guide step as the length unit. Different branches of the interpolation can give exactly the same integer-site readings.

Using only those four measured probabilities, the wave predicts a deep minimum at the next ordinal and much larger returns at the two following ordinals. At ordinals 7, 8 and 9 the fresh native simulation recorded 21,046, 1,043,814 and 1,046,606 detections in a shared bank of 274,877,906,944 original launches. All eight held stations through ordinal 14 lie within the predicted intervals, including both subsequent maxima and minima. The fit and test jointly use 99% simultaneous intervals over all twelve counters. A uniform density and a quadratic band-edge density already disagree with the training readings. Every refused launch, missed charge or detection elsewhere is in the other-outcome bin. [Count-only forecast](smooth_profile.py) · [Predictions and original counts](smooth_profile.json)

For the forecast, c = cos(2k) is inferred from the four consecutive counts. The next intensity obeys n[j+1] = (1+2c)(n[j]−n[j−1])+n[j−2]. Keeping the three interpolation coefficients as polynomials in c allows interval arithmetic to propagate the correlated frequency uncertainty. At this confidence the two admissible k intervals are 0.9623–0.9692 and 2.1724–2.1793 radians per guide step. Integers sample both profiles identically after reversing the sine coefficient. The spatial inference alone cannot distinguish the interpolations between ports. The selected finite clock also admits small amplitudes at other energies, so the eight-site comparison is an empirical test of the single-frequency description. It is a stationary wave on the interior of these guides. A field description in three spatial dimensions would additionally have to predict records on different physical paths, momentum directions, and responding bodies.

```sh
python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.smooth_profile --write
```
