# Radiation changes the area of a ball read by the same Weyl field

A finite Weyl occupation word drives one additional homogeneous clock scale. Calorimetry, regional fermion correlations, clock phases, interior volumes and boundary areas are distinct readings on this material. In its radiation sector a single calibrated coefficient gives the field stress and the contraction of the proper boundary at fixed volume in every enrolled timelike frame. For the specified family of interior packet detectors, their entropy change cancels that area response to leading order. [Construction, hypotheses and proof](PROOF.md)

The additional law gives a cost proportional to proper volume times the square of its fractional expansion rate, and couples a clock phase to the scale. Varying its time lapse relates the squared expansion rate to the energy of the occupied word. Free propagation of the radiation then fixes the pressure to one third of the proper energy density. These choices introduce gravity into the auxiliary homogeneous law. The field's measured occupations, stresses and modular entropy supply independent predictions for its rod readings. Native Interval Echo does not presently select the added scale degree of freedom, its sign, its action or its clock coupling.

## What the held source predicts

Eight occupied particle/hole pairs calibrate the coupling at 2 from independent energy and pulse-spacing readings. The held preparation flips 56 pairs on a \(48^3\) antiperiodic box. Its comoving energy density is 0.302182 and its measured pulse-spacing ratio is 1.056105. The observer locates the boundary by launching spacelike rods in its rest directions and stopping at prescribed enclosed proper volume. One quadrature adds the interior, another measures the boundary. The table reports the largest discrepancy across the ten frames between the resulting curvature coefficient and the independent stress prediction. [Complete records](results.json)

| Target radius | Maximum error in the area-derived Einstein contraction |
| ---: | ---: |
| 0.12 | 0.813355 |
| 0.06 | 0.184117 |
| 0.03 | 0.044984 |
| 0.015 | 0.011183 |
| 0.0075 | 0.002792 |

At the same rest-frame rod radius 0.12, the enclosed proper volume changes from 0.00723823 in the vacuum to 0.00723402 in the occupied word. To enclose the original volume the observer extends the rod to 0.12002328. The corresponding boundary area changes from 0.18095574 to 0.18085029. These are direct comparisons of the reconstructed space at one appointment; the homogeneous clock scale is increasing with cosmological time.

Keep the held word fixed and change the added clock dynamics. At radius 0.015 in the fastest frame the area coefficient is 7.06209 in the declared law and 3.52825 after doubling its geometric kinetic coefficient. Flat clocks and rods give zero within \(1.1\times10^{-11}\). This control measures how much of the result depends on the new law.

Two coarser rest balls also admit direct finite covariance measurements. At radius 0.55 a common interior aperture reads an entropy increase of 0.118246 against the leading prediction 0.145552. The measured entropy plus the calibrated area change is −0.042095. A different reader takes every point inside the changing rim. It enrolls 147 sites in the vacuum and 171 in the excited run, and reports an entropy increase of 2.462998. The instrument choice affects the answer markedly. The theorem states the controlled refinement for an increasingly complete interior packet reader; these coarse results do not establish that limit numerically.

On a separate \(24^3\) realization all ten observers acquire the same eight normalized packet channels. The source is one pure word, with the same cutoff and gravity coefficient. At target radius 0.51 the rest-frame entropy difference is 0.033353; one axial boost at rapidity 0.32 gives 0.037442. Their independent leading stress forecasts are 0.105859 and 0.120812. The finite mode count systematically under-resolves this response. The packet support extends to 0.4386 in their optical chart, inside the sampled minimum clearances 0.513624 and 0.511972 to the actual retimed rims. [Finite moving detector and vacuum transport](core_modes.json)

In moving frames a round conformal diamond and the exponential ball can have different rims, even at the same volume. The proof locates each enrolled packet inside the genuine exponential domain before applying the flat-space conformal vacuum theorem. It then takes the packet count, detector occupation floor, target radii and finite hopping resolution along a common diagonal. Every recorded matter mode enters the entropy. The result applies to isotropic massless radiation; it supplies neither a native realization nor an inhomogeneous or dark-matter stress calculation. [Earlier matter-only experiment](../README.md)

## Reproduction

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.radiation.measure
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.radiation.core_modes
python -m pytest -q code/matched_einstein_limit/source_implication/modular_orbit/weyl_balls/radiation
~~~
