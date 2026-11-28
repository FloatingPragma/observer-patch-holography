# Acoustic polarization from three paid and retained record trails

The 80-owner board gives a conserved excess-to-initial energy ratio of **1037/202 = 5.134**. If that excess acts as cosmological cold matter, the model gives a conditional polarization residual of **307.30 on 310 measured bands**, compared with **293.16** for the standard benchmark. Its ratio lies in the central 68% intervals of both Planck polarization fits. The other two commissioned boards give ratios of **2.104** and **10.520** and large polarization residuals. The three geometries were supplied, and the Planck data and earlier cold-ratio tests were already known. This is a retrospective comparison of all three available levels. [Protocol](PROTOCOL.md) · [Polarization spectra](trail_polarization.pdf) · [Results for every level](trail_results.json)

## What the source fixes

For each geometry, a controller routes one load down the shortest path from the origin to a farthest port. At each seam the controller repeats the addressed echo and sharp read until the load exits with a written head. The completed stock is released into its own free dynamics. The total energy of both the written board and the initial blank control is then conserved. Their whole-board energy difference equals the energy paid during preparation. The three provided levels have the same initial blank energy of 50.5 native units, also equal to the matched blank at the final port. [Native experiment](../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/stock/README.md) · [Native work receipt](../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/stock/trail.json)

| Owners | Written seams | Mean echo/read pairs | Echo work | Sharp-read work | Final whole-board energy | Exact excess / initial |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 14 | 35 | 46.75 | 59.5 | 156.75 | 425/202 |
| 80 | 32 | 80 | 123.25 | 136 | 309.75 | 1037/202 |
| 320 | 64 | 160 | 259.25 | 272 | 581.75 | 2125/202 |

The 320-owner experiment also certifies a lower bound of 478 on the local path excess above a one-load blank at every free time. The cosmological input here is the exact *whole-board* excess, 531.25 for this level. The local theorem does not set the global abundance.

## What the sky tests

The common-mass identification sets the original one-load energy equal to the baryon energy and gives all of the retained excess standard cold, pressureless, adiabatic perturbations by recombination. Each level uses the existing CAMB calculation and the same TT-only calibration of amplitude, angular dilation and tilt. The predictions use every preselected TE/EE band. Their joint score accounts for covariance between TE and EE and conditions on their covariance with the measured TT. The standard control uses the published Planck ratio of 0.12/0.02237. The two posterior columns use original Planck multiparameter chains and give the percentage of posterior weight at or below each fixed candidate ratio. [Measurement and likelihood](README.md) · [Planck chains](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/)

| Owners | Cold/baryon | Equality redshift | TT (155 bins) | TE given TT | EE given TT | Joint given TT (310 bins) | TE+lowE percentile | EE+lowE percentile |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 2.104 | 1,659 | 2,005.78 | 2,540.08 | 2,455.89 | 4,611.86 | 0 | 0 |
| 80 | 5.134 | 3,280 | 159.11 | 163.08 | 145.15 | 307.30 | 21.81 | 75.85 |
| 320 | 10.520 | 6,160 | 1,758.47 | 2,481.85 | 2,265.74 | 4,271.88 | 100 | 100 |
| Standard control | 5.364 | 3,403 | 153.92 | 153.76 | 138.29 | 293.16 | — | — |

Zero and 100 percent indicate positions beyond these finite sampled ranges. The Planck TE+lowE central 95% interval is 4.986–5.484 and the EE+lowE interval is 4.071–5.737. Only the middle of the three geometries lies within both. Their placement is also visible in two of the fixed windows, where the low- and high-ratio models miss in opposite directions. Entries are means of \(D_\ell\) in µK² with correlated errors for the observed averages.

| Channel and multipole | Planck | 20 owners | 80 owners | 320 owners | Standard |
|---|---:|---:|---:|---:|---:|
| TE 400–480 | −57.58 ± 1.20 | −81.71 | −56.76 | −31.44 | −55.29 |
| EE 270–330 | 7.70 ± 0.13 | 10.51 | 7.90 | 5.25 | 7.74 |

At the nominal baryon density the 80-owner and standard conditional means have squared separation 13.84 in the Planck covariance norm. The observed score favors the standard by 14.14. The ordering changes with the baryon density. At 0.0215, 0.02237 and 0.0235 the 80-owner scores are 427.39, 307.30 and 333.42, while the corresponding standard-control scores are 365.41, 293.16 and 384.13. Propagating the temperature calibration errors gives 306.02 for the nominal 80-owner model. Allowing the wide tilt range leaves its nominal score unchanged. All three baryon densities and both tilt ranges leave large residuals for the 20- and 320-owner models. [Scores, nuisance sensitivities and signed window predictions](trail_results.json)

The successful middle level establishes that a commissioned stock can reach the needed amount and reproduce the acoustic pattern under this transfer. Its 259.25 units of excess required 259.25 units of external work. An actual cosmic supplier must furnish that work early enough, and its own stress and any depletion belong in the energy budget. The calculation also requires a physical law determining the abundance and lengths of paths and their pressure, velocities and gravitational perturbations. The data can distinguish these three provided levels. The native experiment supplies no independent physical scale selecting 32 seams.

## Replay

Run in the same environment as the acoustic likelihood comparison.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -m cosmology.cold_acoustic_polarization.trails
python -m cosmology.cold_acoustic_polarization.figure --trails
python -m pytest -q code/cosmology/cold_acoustic_polarization
~~~
