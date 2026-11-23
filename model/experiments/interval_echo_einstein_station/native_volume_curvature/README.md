# Can a volume deficit in the native route survey determine curvature?

On the finest blank board, doubling the clocked route radius increases the enrolled population from 1,158 to 4,630. If these counts measure volume in a homogeneous three-dimensional geodesic ball, the apparent contraction relative to cubic growth fixes its sectional curvature and predicts a **21.65% area deficit at equal volume**. Independently clocked distances between five stations in that very ball contradict the inferred geometry. Even allowing each distance to change independently, some range must change by more than **6.4%**. The inferred curvature also fails to predict the two held finer volume ratios. These deterministic tests use the registered blank addressed-route surveys. [Original census and clocked ranges](../proper_ball_ruler/coarse.json) · [Analysis protocol](PROTOCOL.md) · [Numerical results](results.json)

## What a volume can determine

In a convex three-dimensional region of constant positive sectional curvature K, define x = √K r and F(x) = x − sin(2x)/2. Integrating the area of geodesic spheres gives V(r) = 2πF(x)/K^(3/2) and A(r) = 4π sin²(x)/K. The volume doubling ratio is F(2x)/F(x). For 0 < x ≤ π/2, this ratio equals eight times the sin²(t)-weighted mean of cos²(t) over 0 ≤ t ≤ x, and hence decreases strictly from eight to two. Flat space gives eight and the negatively curved space form gives a larger ratio. Treating the counts as volume fixes the positive solution without adjusting any volume unit.

The flat ball of the same volume has radius (3F(x)/2)^(1/3)/√K. Its curved-to-flat boundary area ratio is sin²(x)/(3F(x)/2)^(2/3). At sufficiently small x, this fixed-volume deficit measures curvature to leading order. In a general three-dimensional spatial slice the first isotropic volume correction measures scalar curvature. Recovering a spacetime field equation additionally requires observer frames, an independently identified matter source, and controlled variations about a reference state. [Jacobson’s fixed-volume derivation](https://arxiv.org/abs/1505.04753)

## Predictions on other records

All four rows use the existing board census. The volume estimator sees only the radius and the two ball populations on its row. The twelve landmark ranges are reserved for the separate distance test. The labels in the last column refer to changes expressed as a percentage of the recorded length.

| Host cells | Ball radius / home range | Volume ratio | Inferred √K r | Predicted area / flat at equal volume | Lower bound on held range distortion, rounded |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 80 | 0.093750 | 2.9583 | 1.29051 | 0.69339 | 13.70% |
| 1,280 | 0.046875 | 4.0282 | 1.06409 | 0.78577 | 6.45% |
| 20,480 | 0.023438 | 4.0490 | 1.06003 | 0.78730 | 6.59% |
| 327,680 | 0.011719 | 3.9983 | 1.06993 | 0.78354 | 6.42% |

The one-swap-thick outer shell on the finest board accounts for 10.02% of the enrolled ports, compared with the volume model forecast of 10.31%. The distance measurements discriminate more strongly here. Every matrix of geodesic distances in a common convex spherical 3-space patch must give a positive semidefinite cosine Gram matrix of rank at most four, since its entries are inner products of normalized four-vectors. The initial spectral test gives required range distortions of at least 9.11%, 5.80%, 5.49%, and 5.02%. A supplemental interval calculation on exactly the same five-landmark witnesses gives the tighter limits tabulated above. For each pair of stations the allowed interval of its angular distance gives the maximum contribution it could make to the Gram quadratic form of the observed negative eigenvector. At 5% tolerance the sum of even these pairwise maxima on the finest board remains −0.06756. Any genuine spherical Gram must give a nonnegative value. [Calculation and numerical certificates](measure.py)

If curvature in units of the independent home-to-farthest range stays fixed at the 1,280-cell calibration, the predicted volume doubling ratios on the next two boards are **6.7436 and 7.6669**, compared with **4.0490 and 3.9983** observed. Fitting each separately drives K in squared inverse home-range units from 515 to 2,046 to 8,336, even as the balls occupy progressively smaller fractions of the surveyed host. The approximately quadratic count growth documented in the route census accounts for the persistence of a ratio near four. [Previous native scaling result](../proper_ball_ruler/README.md)

All four panels describe the same blank preparation class. They support a test of whether registered port count and clocked route length have a common constant-curvature three-dimensional interpretation. The 0.78354 area ratio above comes from the rejected joint model. A source-induced proper-area deficit requires a comparison of perturbed and reference preparations with independently equalized proper volume and one common validated ruler. Static route ballots and arbitrary enclosed shapes do not supply this comparison. The finite-density contact apparatus separately exhibits scale-dependent cubed-aperture to squared-volume ratios on the same hosts. The freely stationary source family has independent entropy slopes at fixed mean routed hull area and volume. [Finite-density contact surveys](../density_area/README.md) · [Free sources at matched geometric expectations](../proper_ball_ruler/README.md)

Reanalyze the archived observations from the repository root with the following commands.

```sh
python -m model.experiments.interval_echo_einstein_station.native_volume_curvature.measure --output /tmp/native-volume-curvature.json
python -m pytest -q model/experiments/interval_echo_einstein_station/native_volume_curvature
```
