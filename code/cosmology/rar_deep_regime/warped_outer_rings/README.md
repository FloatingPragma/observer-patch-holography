# Resolved exterior orbits of NGC 2685

The smooth exterior profile calibrated on SPARC predicts 149.2, 147.5, 146.1 and 145.5 km/s at four successive radii of NGC 2685. The measured circular speeds are 143.9, 139.6, 146.7 and 168.8 km/s, each with an assigned error of 15 km/s. The last ring sits 1.56 quoted errors above the forecast. Spectral-line fits across neighboring rings share distance, inclination systematics and the same observations. These four radii test one held galaxy. [Published H I data, Appendix C](https://arxiv.org/abs/0810.3025) · [Reconstructed numbers](results.json)

## Locating the baryons before comparing speeds

The sixteen-galaxy exterior comparison assigned every unlocated gas atom to the outside when testing enclosure. That gave a lower bound of 82.89% for NGC 2685. The resolved 21-cm observations locate much of the gas. Use the published stellar enclosure at 320 arcsec as a lower bound at all larger angles, the same conservative stellar mass-to-light ratio 0.6, a one-percent stellar allowance, and the original upper atomic inventory including helium and a 0.10 dex allowance. The prior speed entry at 320 arcsec was based on this same H I study. The 360, 400 and 420 arcsec rings extend the tested radial range. Selection uses stellar enclosure, gas columns, five effective radii and compact-source acceleration below 1.2 × 10⁻¹¹ m/s². NGC 2685 supplied no SPARC training radii. The selection and data reuse are retrospective. [Original photometric and gas-mass calculation](../outer_orbit_transfer/PROTOCOL.md) · [Provenance of the outer H I velocity](https://arxiv.org/abs/1509.05236)

At angular radius θ integrate the face-on H I columns over tilted rings with radii at most 0.9θ. Lower each published knot by one quoted error, clip negative lower columns at zero, and linearly interpolate between knots. The integrated mass, with the same helium convention as the catalogue, exceeds 1.688 billion solar masses under these allowances. The tilt of a circular ring leaves its area and its three-dimensional distance to the center unchanged. Gas with vertical distance at most θ√0.19 from that ring plane is inside the spherical orbit. At 320 arcsec at least 90.61% of the conservative inner gas column must satisfy this condition to pass the 90% overall mass cut. The corresponding height is 139.5 arcsec, and an rms height below 42.76 arcsec suffices by the second-moment bound. The published high- and low-resolution scale heights are 5.4 and 3.0 arcsec. Converting a scale height to an rms height requires a vertical profile.

Enrolling 95% of the inner column gives an enclosed baryon fraction of 90.345% at all four radii. Holding the certified interior fixed would admit another 82.5 million solar masses wholly outside before failing the cut. A twice-error downward allowance on all columns gives 88.10%. This is a sensitivity calculation on the declared inventory. The errors estimated from the two sides of the galaxy, with manual error assignments at some extreme rings, have no joint coverage specification. Additional interior molecular gas adds to both the numerator and denominator. Observations have detected molecular gas in the inner ring. A complete inventory of diffuse material outside the measured disk would require another observation. [H I geometry and error construction](https://arxiv.org/abs/0810.3025) · [Interferometric CO observations](https://arxiv.org/abs/astro-ph/0209004)

## An orbit bound on a finite halo

For a compact baryonic monopole of mass M and a spherical exterior excess, write the squared speed as GM/r + √(GM a₀) min(1,L/r) when the halo ends sharply at L. This specifies an effective pressure and a classical orbit response. The held smooth SPARC scale is 8.7504 × 10⁻¹¹ m/s². The adopted stellar mass-to-light ratio is 0.8, catalogue distance 16.1 Mpc, and modelled baryon mass 2.681 × 10¹⁰ solar masses. Enclosure alone does not fix the force of the real warped disk. [SPARC radial scale](../receipts/lensing_profile.json) · [Finite halo model](../finite_halo/README.md)

| H I radius, kpc | Measured km/s | Smooth forecast, km/s | Finite-spacing forecast range, km/s |
| ---: | ---: | ---: | ---: |
| 24.98 | 143.9 ± 15 | 149.2 | 146.8–148.0 |
| 28.10 | 139.6 ± 15 | 147.5 | 146.0–147.0 |
| 31.23 | 146.7 ± 15 | 146.1 | 144.3–145.3 |
| 32.79 | 168.8 ± 15 | 145.5 | 143.6–144.5 |

At fixed nominal geometry, mass and smooth calibration, every sharp edge misses at least one ring by one quoted speed error. At two errors the complete feasible range is **L ≥ 29.26 kpc**, with no upper constraint. Increasing stellar mass-to-light to 1.0, atomic gas by 0.10 dex and the SPARC scale to its upper galaxy-bootstrap endpoint gives a necessary lower bound of 18.67 kpc after also raising the distance and lowering all speeds by two quoted errors. This bound applies at the stated inclinations. A common distance adjustment to the smooth model can reduce its largest deviation to 1.044 errors at 17.84 Mpc. That adjustment uses the measured speeds and lies outside the held forecast. The separate frozen 0.5-kpc spacing and its own SPARC normalization give a largest residual of 1.62–1.68 errors when all four rings share one phase. The errors do not resolve that spacing. [Ring intersections and common-phase calculation](measure.py) · [Cut-spacing calibration](../cut_spacing/README.md)

## Can one mass profile also describe the diffuse dwarf?

If the same nondecreasing excess mass fraction f applies to galaxies as a function of x = r/√M, ordering x orders f. This test needs no acceleration scale or specified outer edge. The outermost NGC 2685 ring is at smaller x than the outer AGC 114905 ring. With the compact baryonic force, NGC 2685 needs an excess of 7.10 times its baryonic mass there. At the dwarf ring, assigning all measured squared speed to the excess gives an upper ratio of 1.142. The discrepancy persists if NGC 2685 is allowed higher stellar and atomic mass, one higher distance error and one lower speed error. Its required ratio exceeds 3.14. Even the generous dwarf mass, distance and tilt case admits only 2.43. Allowing two speed and distance errors on NGC 2685 reduces its lower bound to 1.74, or 1.34 with an additional tilt variation. Those broad boxes overlap.

This excludes a universal monotone extension on the central inputs. AGC 114905 is a diffuse disk whose outer measured orbit fails the compact-source assumption. Describing that disk would require a source law beyond this exterior test. Neither galaxy measures the microscopic collar charge or its metric coupling. The astronomical result is independent of the static OpenWave criterion requiring a stable neutral massive state. [Dwarf velocity data](https://arxiv.org/abs/2404.06537) · [Dwarf mass and inclination allowances](../finite_halo/transfer.json) · [OpenWave criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L121)

## Reproduction

The released source archive and the typeset PDF independently reproduce all 21 input rows. Rebuilding the comparison and tests requires only the archived columns. To check their extraction separately, download v1 of the source and PDF and pass the local paths to the extraction command. The PDF check uses Poppler.

~~~sh
PYTHONPATH=code:. python -m cosmology.rar_deep_regime.warped_outer_rings.acquire source.tar.gz --pdf paper.pdf
PYTHONPATH=code:. python -m cosmology.rar_deep_regime.warped_outer_rings.measure --write
PYTHONPATH=code:. python -m pytest -q -W error code/cosmology/rar_deep_regime/warped_outer_rings
~~~

[Authors’ source archive](https://arxiv.org/src/0810.3025v1) · [Typeset Table C.1](https://arxiv.org/pdf/0810.3025v1)
