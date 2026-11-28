# Does a gas-only gravitational source transfer beyond the motivating pair?

NGC 3521 / NGC 3198 gave a 1.56 ratio of measured gas boundaries and a 1.54 outer excess ratio when the ordinary force used the same full-phase maps. A fixed gas-only law performs substantially worse on other galaxies. With NGC 3198 as its **only calibrator**, two-dimensional H I and CO maps give a median absolute total-speed error of 18.6 km/s over the other eight galaxies. In NGC 2841 the predicted speed falls 99.0 km/s below the observation. NGC 5055 and NGC 6946 overshoot by 78.6 and 75.8 km/s. Four of eight galaxies require gains incompatible with the calibrator when both speeds move independently by two quoted errors at fixed geometry and baryonic force. [All image predictions and error boxes](gas_image_transfer.json) · [Measured source strengths, FITS hashes and gas coverage](gas_image_inputs.json)

## Selection and evidence reuse

The motivating pair was first selected on baryonic radial profiles. Its gas-only agreement was noticed later, and the other eight orbital measurements have already been used in this project. Before computing their gas boundaries we fixed a single calibrator, NGC 3198, and used the same estimator for that galaxy and each of the eight. NGC 3521 supplies neither a gain nor a score. Every input to each held source is a baryonic measurement. No speed from the eight galaxies enters calibration. All original outer orbits and molecular non-detections remain in the sample. The second, image-based design was fixed after the radial results were known. These are retrospective held applications with respect to the motivating pair. [Radial transfer protocol](GAS_TRANSFER_PROTOCOL.md) · [Image transfer protocol](GAS_IMAGE_TRANSFER_PROTOCOL.md)

| Input to Y=kS | Median absolute total-speed error (km/s) | RMS total-speed error (km/s) |
|---|---:|---:|
| Optical radial profiles, all measured baryons | 22.1 | 89.8 |
| Optical radial profiles, H I and detected H₂ | 34.9 | 68.7 |
| Optical radial profiles, H I alone | 13.7 | 57.3 |
| Full H I image and detected molecular image in the optical aperture | 18.6 | 53.1 |
| Same measured gases, finite-source response at the actual orbit | 15.0 | 56.7 |
| Same gases with H I normalized to the catalogue’s total atomic mass | 17.3 | 51.8 |
| Full H I image alone | 18.2 | 42.5 |
| Full H I and molecular images plus the measured radial stars | 16.8 | 75.8 |

The all-baryon optical predictor has smaller median speed error than the optical gas-only predictor, alongside much larger errors for specific galaxies. Their respective RMS errors in excess speed squared are 48,324 and 28,852 (km/s)². The finite optical two-gas model gives 26.6 km/s median speed error and 68.6 km/s RMS. Using actual extended gas images reduces the compact two-gas median and RMS errors. These aggregate comparisons do not select a universal material law. [Signed radial predictions and excess errors](gas_transfer.json) · [All map variants and sources](gas_image_transfer.json)

| Galaxy | Measured speed (km/s) | Optical H I+H₂ residual | Mapped H I+H₂ residual | Mapped gases + radial stars residual |
|---|---:|---:|---:|---:|
| DDO 154 | 45.5 | −7.3 | −4.0 | −15.9 |
| IC 2574 | 67.5 | +42.8 | +22.4 | −2.4 |
| NGC 2403 | 134.0 | −27.1 | −14.7 | −17.6 |
| NGC 2841 | 294.0 | −141.0 | −99.0 | +4.0 |
| NGC 5055 | 172.0 | +63.6 | +78.6 | +147.4 |
| NGC 6946 | 154.0 | +103.9 | +75.8 | +96.1 |
| NGC 7331 | 238.0 | +17.9 | +8.5 | +119.6 |
| NGC 7793 | 90.8 | +9.3 | −7.0 | +7.9 |

Residuals are predicted minus observed **total speed** in km/s. The fitted quantity on NGC 3198 is exclusively its outer **excess speed squared** after subtracting the published ordinary gravity and measured molecular correction. The same subtraction defines each held target. The radial two-gas source predicts only 0.155 times the observed excess for NGC 2841 and 4.03 times for NGC 6946. At 0.8 kpc resolution the image-based two-gas law fails the two-speed-error common-gain check in NGC 2403, 2841, 5055 and 6946. Those with compatible error intervals include both good central predictions and the relatively uncertain IC 2574 orbital speed. The test measures conditional fit under nominal geometry. It is neither a survey likelihood nor a statement of formal exclusion for every stellar mass, inclination or projection model. [Speed uncertainties and raw endpoint masses](results.json) · [Fixed-geometry image comparisons](gas_image_transfer.json)

## Flux, aperture and noise

The atomic source uses the full optically thin THINGS footprint. The independent HERACLES molecular maps use detections ≥3 local errors within each galaxy's optical source aperture. The archival CO images include structure that the earlier Leroy radial table did not report, including small amounts of emission in IC 2574 and NGC 2403. DDO 154 has a CO image but no material detections at this threshold to the displayed 10⁶ M☉ precision. NGC 7793 has no HERACLES release, so the two-gas estimator uses atomic gas alone there. No galaxy is dropped. Reported H I mass coverage versus the SPARC global inventory spans 0.64–1.35 including the calibrator. NGC 3198’s downloaded map has a lower integrated flux than the survey table. The other eight image totals reproduce their THINGS catalogue fluxes to about 1%, before comparison to the separately compiled SPARC inventories. Ordinary gas phases beyond the molecular survey footprint, CO-dark gas and unresolved compact structures can alter the true source.

Galaxies remain at their original catalogue distances in this transfer. A common 0.8 kpc cell size does not erase differences in the physical telescope beams or the disk inclinations. The 1.6 kpc combined-gas control gives a 22.4 km/s held median error and 48.6 km/s RMS. The finite-source control assigns the measured boundary at each location to a spherical shell before reading gravity. The new data test a conversion from visible gas columns to proposed source strength. They supply no microscopic probability law selecting the gas or a force normalization supplied by OPH. [THINGS survey](https://arxiv.org/abs/0810.2125) · [Original atomic maps](https://www2.mpia-hd.mpg.de/THINGS/Data.html) · [HERACLES survey](https://doi.org/10.1088/0004-6256/137/6/4670) · [Original IRAM maps](https://www.iram.fr/ILPA/LP001/)

The two offline runs are `python -m cosmology.rar_deep_regime.boundary_source.gas_transfer` and `python -m cosmology.rar_deep_regime.boundary_source.gas_image_transfer` with `PYTHONPATH=code:.` from the repository root. The image run also accepts `--rebuild-from-fits /tmp/oph-boundary-fits` when Astropy is available. This reprojects the original intensity and uncertainty images, reconstructs each source and writes a fresh source summary before scoring.
