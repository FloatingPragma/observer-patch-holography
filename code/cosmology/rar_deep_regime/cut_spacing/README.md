# Where the compact-source law can be tested

At 9.95 kpc in AGC 114905 the fitted local disk formula predicts at least 62.97 km/s. The published, pressure-corrected circular speed is 26.94 ± 1.83 km/s. The authors measured the gas distribution as well as the rotation. Integrating their fiducial surface density, independently of the speed points, gives **66.9% of the modeled gas mass inside this cylindrical radius**. The gas supplies 94% of the baryonic mass. Even putting every star inside gives an upper bound of **68.9%** on the central baryonic enclosure for the fitted profile; a sphere of that radius contains at most that much. The model gas profile reaches 90% enclosure in projection at 14.30 kpc. These are integrations of the authors’ fitted radial profile, including its outer extrapolation. Their fit parameters and distance specify the geometry; the enclosure values have no propagated observational uncertainty. [Original mass model and kinematics](https://doi.org/10.1051/0004-6361/202450230) · [Calculation for both published reductions](disk_domain.json)

## Which inference fails on this disk?

The compact deep exterior uses a single enclosed baryonic mass, an asymptotically linear count of active cuts and one normalization. The step to an interior disk ring replaces the spherical acceleration by the independently estimated midplane force from stars and gas. The OPH exterior theorem does not specify how those spatially distributed baryons excite the cuts. A gas tracer can measure an exterior orbit around a compact mass. A high gas fraction alone neither excludes a system nor establishes the exterior approximation.

There is an elementary ambiguity even in spherical symmetry. Let the exterior enclosed excess be C r beyond source radius Rₛ. Inside Rₛ, either C Rₛ (r/Rₛ)³ or C Rₛ (r/Rₛ)⁵ is continuous at the match, has nonnegative bounded density and gives exactly the same exterior. At half Rₛ the first interior acceleration is four times the second. These examples show why the exterior profile leaves the interior force undetermined. Neither function was fitted to the galaxy.

| Outer-ring quantity | Fiducial geometry | Published alternative geometry |
| --- | ---: | ---: |
| Radius (kpc) | 9.95 | 10.13 |
| Published circular speed (km/s) | 26.94 ± 1.83 | 24.25 ± 1.83 |
| Baryonic circular speed (km/s) | 22.57 | 20.75 |
| Lowest phase-conditioned mean speed (km/s) | 62.97 | 61.20 |
| Observed excess squared speed / lowest modeled excess | 0.0626 | 0.0475 |

The fiducial central excess is only one sixteenth of the least excess allowed by the fitted disk model. The equivalent spherical excess masses are 0.50 versus 7.99 billion solar masses. These restate a radial acceleration using G; they do not measure a spherical halo directly. The velocity-only inclination needed at unchanged projected speed would be about 12.7°, against the photometric 31 ± 2°. The inclination estimate assumes the published morphological interpretation. The 12.7° diagnostic holds the other inputs fixed and neglects the change in pressure correction; reanalysing the data cube would require recomputing both. Both columns use the same H I observations. Their difference measures sensitivity to the reduction. For the fiducial column the authors’ two plotted acceleration curves recover the same outer baryonic speed to 0.004 km/s. [Published figure extraction](../udg_transfer/README.md) · [Reproducible diagnostics](disk_domain.py)

Using the whole baryonic mass in the compact asymptote gives 63.58 km/s for the anomalous speed alone. It answers a question about a distant orbit. Determining the onset radius requires a source law. An exterior test can select on baryonic enclosure and low acceleration without consulting the measured speed. In a separate sixteen-galaxy early-type sample, NGC 3838 alone passes the adopted 90% enclosure, five effective radii and deep-acceleration cuts. Its outer speed is 159 ± 14 km/s, compared with 153–160 km/s using the published mean-shear normalization intervals and nominal geometry. This is one conditional cross-check. Other radii and the shapes of the shear profiles constrain the same interpretation. [Outer-orbit and shear comparison](../outer_orbit_transfer/README.md)

With the SPARC normalization and spacing held fixed, the point-mass exterior itself gives the following ranges. They vary only the shared-phase count for each orbit. The ranges condition on the adopted baryonic masses and distances. None of these galaxies enters the SPARC fit. [Enclosure and speed figure](comparison.pdf) · [All transferred orbits](disk_domain.json)

| Galaxy | Lower bound on enclosure | Predicted speed (km/s) | Published speed (km/s) |
| --- | ---: | ---: | ---: |
| NGC 3838 | 95.5% | 136.8–138.3 | 159 ± 14 |
| NGC 2685 | 82.9% | 146.8–148.0 | 144 ± 10 |
| NGC 6798 | 89.9% | 186.2–187.5 | 190 ± 8 |

## Rotation fit and external shear

In 121 SPARC galaxies, 1,240 radii meet the selection on quality, inclination and baryonic acceleration. One unknown cut phase is integrated jointly over all radii in each galaxy. All five training folds select a common spacing of 0.5 kpc. The held predictive gain over a smooth mean is 0.470 nats per galaxy, with a paired fixed-prediction bootstrap interval [0.092, 0.895]. A smooth comparison allowing extra small-radius scatter gives a difference of 0.052 nats per galaxy, interval [−0.481, 0.735]. These retrospective data discriminate poorly between a cut signature and radial scatter. The fitted full-sample scale is 8.375 × 10⁻¹¹ m/s². Averaging held predictions over the catalogue distance error gives a gain of 0.486 nats per galaxy against the original smooth comparison. This last sensitivity can overlap with the fitted common scatter. [Source catalogue](https://arxiv.org/abs/1606.09251) · [Protocol](PROTOCOL.md) · [Fit and external predictions](results.json) · [Smooth-scatter and distance controls](controls.json)

At the fitted scale, even putting the entire one-cut deficit at the last dwarf ring requires a spacing of at least 8.77 kpc to reach the measured central speed plus two plotted errors. This bound concerns the phase-conditioned mean. With the fitted scatter, the outer-ring lower-tail probability is 0.00036 under the conditional velocity model. For compact spherical lenses and selection independent of phase, the mean excess mass at each radius is exactly the smooth mass. Linearity of weak shear therefore preserves the mean prediction at fixed normalization for every spacing. With nominal baryonic masses, negligible slip and bin-centre evaluation, the six bins wholly above the published isolation threshold give covariance-weighted squared residuals of 44.37 (KiDS-bright) and 17.07 (GAMA). The two samples use shared imaging. Restricting the orbit claim to compact exteriors leaves the lens-shape test in force. [Published lensing measurements](https://doi.org/10.1051/0004-6361/202040108) · [Effects of mass and bin weighting](../outer_orbit_transfer/README.md)

## Reproduction

From the repository root, run the following commands. The fit verifies the SPARC source snapshot. The figure reader in the linked data package checks both author figure hashes and ten extracted ring records.

```sh
python code/cosmology/rar_deep_regime/cut_spacing/fit.py
python code/cosmology/rar_deep_regime/cut_spacing/controls.py
python code/cosmology/rar_deep_regime/cut_spacing/disk_domain.py
python code/cosmology/rar_deep_regime/cut_spacing/figure.py
python -m pytest -q code/cosmology/rar_deep_regime/cut_spacing
```

To check the shear inputs against the original survey files, download the [released KiDS and GAMA archive](https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar) and run `python code/cosmology/rar_deep_regime/lensing_profile.py /path/to/brouwer2021_rar.tar` from the repository root. The reader verifies the archive digest, reconstructs the shear calibration and covariance, and compares them with the inputs used by the fit.
