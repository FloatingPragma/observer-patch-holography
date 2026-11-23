# Transfer of the deep acceleration scale to gas-rich ultra-diffuse galaxies

A scale fitted on 584 radii in 74 SPARC galaxies predicts limiting speeds of 58.0–69.3 km/s for six external gas-rich ultra-diffuse galaxies. Their published resolved H I speeds are 19–37 km/s. The median predicted-to-measured ratio is 2.27. A higher-resolution observation of AGC 114905 independently measures its inclination and resolves five rings. Applying the fitted radial acceleration law to the measured stars-and-gas distribution predicts 62.91 km/s at 9.95 kpc. The corrected circular speed there is 26.94 ± 1.83 km/s. [Comparative figure](comparison.pdf) · [Acquisition and test conventions](PROTOCOL.md) · [Numerical results](results.json)

## Transfer between catalogues

Fix a₀ to 7.603 × 10⁻¹¹ m/s² with the 74-training-galaxy cluster bootstrap range 6.192–9.221 × 10⁻¹¹ m/s². The constant comes from fitting gobs = gbar + √(gbar a₀) in the SPARC sample. The OPH deep profile yields the asymptotic anomalous speed (G Mbar a₀)¹⁄⁴, with mass from a separate optical and H I catalogue. All six external objects from the authors' first table enter this test. The IAU nominal solar mass parameter is 1.3271244 × 10²⁰ m³/s². No observational scale or nuisance parameter is refitted to these galaxies. [Archived SPARC calibration](../receipts/held_galaxies.json) · [Six-galaxy source](https://doi.org/10.3847/2041-8213/ab40c7) · [Transcribed measurements](observations.json)

| AGC | Measured circular speed (km/s) | Predicted limiting speed (km/s) | Published inclination | Inclination required by the central prediction |
|---:|---:|---:|---:|---:|
| 114905 | 19 (+6/−4) | 63.60 | 33° | 9.36° |
| 122966 | 37 (+6/−5) | 63.60 | 34° | 18.98° |
| 219533 | 37 (+5/−6) | 69.34 | 42° | 20.92° |
| 248945 | 27 ± 3 | 58.01 | 66° | 25.17° |
| 334315 | 26 (+4/−3) | 67.76 | 52° | 17.60° |
| 749290 | 26 ± 6 | 62.15 | 39° | 15.26° |

Even at the lower training-scale endpoint, lowering each log mass by twice its tabulated uncertainty and raising each measured speed by twice its upper uncertainty leaves positive speed discrepancies for all six. The smallest is 1.27 km/s. These boxes are a sensitivity calculation with no joint coverage assignment. Inclinations use the projected-velocity approximation and the central mass. The published uncertainties in speed already include inclination uncertainty. The authors' gas-pressure corrections are small in the six-object reduction.

The observed radii are approximately 8–18 kpc. Applying the compact limiting speed to those finite radii assumes they have reached the asymptotic regime. The observations constrain that assignment. A slower approach outside the measured H I discs requires independent dynamics. These six objects were chosen from a published sample known for slow rotation. They test external transfer to that population and cannot establish its frequency among unselected galaxies.

## Resolved outer rotation at a measured inclination

Deep optical imaging gives inclination 31 ± 2° for AGC 114905 at distance 78.7 ± 1.5 Mpc. The updated stars-and-gas mass is (1.47 ± 0.14) × 10⁹ solar masses. In the five fiducial rings, baryonic accelerations span 1.63–1.83 × 10⁻¹² m/s², all inside the SPARC comparison cut of 1.2 × 10⁻¹¹ m/s². The applied local-force version of the model adds √(R Vbar² a₀) to the baryonic squared speed. It extends the deep exterior profile to the force of the measured disc. [Optical and kinematic measurements](https://doi.org/10.1051/0004-6361/202450230) · [Fiducial and sensitivity ring inputs](rings.json)

| Radius (kpc) | Measured circular speed (km/s) | Baryonic speed (km/s) | Predicted circular speed (km/s) |
|---:|---:|---:|---:|
| 1.11 | 10.81 ± 1.80 | 7.45 | 20.85 |
| 3.32 | 18.13 ± 1.83 | 13.67 | 37.31 |
| 5.53 | 24.82 ± 1.77 | 17.58 | 48.08 |
| 7.74 | 27.19 ± 2.17 | 20.86 | 56.97 |
| 9.95 | 26.94 ± 1.83 | 22.57 | 62.91 |

The last ring misses by 19.65 plotted speed-error widths, or 18.15 at the lower endpoint of the training scale. A compact limiting-speed interpretation of the same updated mass gives 62.06 km/s and would need an inclination near 12.92° at the observed outer speed. A second reduction at a distance of 80.27 Mpc and an inclination of 35.86°, chosen by the observational authors for their baryons-only fit, gives a direct stars-and-gas speed of 20.75 km/s and a predicted 60.25 km/s versus 24.25 ± 1.83 at 10.13 kpc. Both reductions use the same underlying H I observations. The plotted speed errors condition on a geometric model. The error widths are discrepancy units and do not supply a joint significance over geometry, intrinsic disc ellipticity, baryonic force, noncircular motions and inter-ring covariance.

## Reproduction

The six-object snapshot transcribes numerical columns from the publisher's Table 1. The five-ring extraction uses the authors' vector figures 7 and 8 from arXiv v2. The reader checks the source archive and individual figure hashes, maps marker and uncertainty endpoints through labeled axis ticks, and interpolates the encoded baryonic curve at the marker locations. In figure 8 the two displayed forward curves encode the same separately estimated baryonic force. Their inversion uses the authors' equation, Vbar² = Vmodel⁴/(Vmodel² + aR). The inferred baryonic speeds differ by at most 0.051 km/s from plotting and interpolation. Figure 7 supplies a direct photometric force curve at its separate geometry. The uncertainty band in the right panel of the comparative figure varies only the SPARC training scale. [Published six-galaxy table](https://openresearch-repository.anu.edu.au/server/api/core/bitstreams/fa77458d-b5d2-41a4-bcc3-e72e8bcf24a2/content) · [Author figure archive](https://arxiv.org/src/2404.06537v2)

From the repository root, run the following commands to check the published data extraction and the two scientific regressions.

```sh
curl -fLsS https://arxiv.org/src/2404.06537v2 -o /tmp/agc114905-source.tar
python code/cosmology/rar_deep_regime/udg_transfer/read_figures.py /tmp/agc114905-source.tar
python code/cosmology/rar_deep_regime/udg_transfer/measure.py
python -m pytest -q code/cosmology/rar_deep_regime/udg_transfer
```
