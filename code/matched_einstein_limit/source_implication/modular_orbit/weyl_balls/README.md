# Correlated Weyl matter on overlapping observer balls

A massless Weyl hopping law supplies regional entropy and a separately measured energy on a single preparation across radii and moving frames. For fixed nonzero excitations the entropy difference approaches \(8\pi^2 R^4 T_{uu}/15\), the coefficient required by the finite Einstein tower. The result includes a spatially homogeneous source with directed energy flux and localized charged packets under free evolution. The analytic construction uses a common stationary ground state, one observer readout at each resolution and a single finite approximation for all the named balls and frames. [Proof and order of limits](PROOF.md)

The physical reason is visible in the records. An observer can measure correlations between separated modes. Changing the radius changes which correlated partners enter the same ball, so the energy of a given packet has a different modular weight on the larger region. The clock comes from the hopping action. The modular weights come from independent measurements of the regional covariance. In the Weyl vacuum the continuum theorem identifies these two readings, including the factor \(2\pi\). As a control, the identical finite homogeneous change placed on a tracial product reference has zero linear entropy response and a negative finite entropy change. [Finite comparison](results.json) · [Vacuum modular flow of the massless Weyl field](https://doi.org/10.1007/s00220-025-05253-1)

## Finite records

The common vacuum fills the negative band of a nearest-neighbour spinor hopping on an antiperiodic box. From any exactly half-filled input, a work reservoir and uniformly sampled mode exchanges prepare this ground state with probability one. The stopping certificate counts occupied negative-energy modes. A calibrated local readout mixes a fraction \(2\delta\) of independent tracial signal into every channel. This makes the logarithm of every measured regional density well defined. Both the action energy and the modular entropy use these same output modes. In the homogeneous experiment the particle occupation is
\[
 f(p)=\tfrac14(1-|p|^2/P^2)_+^2[1+d\cdot p/P],
 \qquad P=2.5,\quad d=(0.6,-0.3,0.2);
\]
the hole occupation at \(p\) equals \(f(-p)\). This fixed spectrum gives nonzero momentum. The signal's continuum rest density is \(P^4/(96\pi^2)\), its pressure in each direction is \(\rho/3\), and its energy flux is \(64\rho\,d/315\). These follow by integrating the independently specified dispersion and occupations. No field equation for a metric enters the calculations.

Two balls share a single acquisition at each of three resolutions. The smaller ball has radius \(0.8R\), with center displaced by \(0.1R\). The last two columns use the exact measured rest-energy density and the geometric continuum coefficient. The numerical ratios include lattice boundary and resolution errors.

| Box side in cells | Outer radius | Outer \(\Delta S/[(8\pi^2/15)R^4\rho]\) | Smaller ball ratio |
| ---: | ---: | ---: | ---: |
| 32 | 0.75 | 0.639 | 0.992 |
| 48 | 0.60 | 0.795 | 0.882 |
| 64 | 0.50 | 0.834 | 0.859 |

On the finest outer ball the entropy difference is 0.011236, the independently summed weighted hopping is 0.013536, and the relative entropy is \(7.18\times10^{-5}\). The product-reference control gives \(-1.12\times10^{-5}\) for that source. The small grids illustrate the finite discrepancy. The convergence assertion comes from the proof with its stated ordering of the limits. [Six regional measurements and controls](results.json)

Ten observers read the same source on one further \(12^3\) box. Applying an action-derived unitary to a packet of mode columns on the common Cauchy record implements each tilt. The ten timelike contractions span the ten-dimensional space of symmetric tensors. At rapidity 0.32 the two opposite \(x\)-direction readings give a signed energy-flux test. Across the whole set the signal is finite and the model's high-resolution tensor prediction uses only its rest spectrum and the chosen rapidities. [Joint moving-frame record](boosts.json)

| Frame generator direction | \(\Delta S\) | Independently summed weighted hopping |
| --- | ---: | ---: |
| Rest | 0.065490 | 0.070458 |
| \(+x\) | 0.074258 | 0.080789 |
| \(-x\) | 0.067120 | 0.072310 |
| \(+y\) | 0.068909 | 0.074430 |
| \(-y\) | 0.072477 | 0.078669 |

One other experiment prepares a localized particle with occupation 0.65 and lets the hopping act freely. In the first unit of elapsed time its occupation excess drops from 0.03796 to 0.02918 in the launch ball and grows from 0.02692 to 0.03662 in the receiving ball. There the entropy change grows from 0.042985 to 0.056161, and the independent weighted energy from 0.053324 to 0.070917. The packet has a nonzero local occupation excess. The proof controls its relative entropy by \(O(R^5)\); charge-balanced homogeneous sources have the sharper \(O(R^7)\) bound. [Free packet readings](packet.json)

## Application to Interval Echo

This is an auxiliary Weyl sector with an explicit propagator, prepared vacuum, readout algebra and boost. An independent copy can be appended to an already published Interval Echo process while preserving the marginal distribution of every original observer protocol. Entropy differences on the added sector then have the properties above. For an Interval Echo gravity result, native observers still need to identify those fields and clocks in the published records, or an explicit joint law must connect this sector to their calibrated proper volume and area. These flat matter measurements supply no area susceptibility, ordinary universal gravitational coupling or independently measured collar stress. The homogeneous source's positive pressures also differ from those of a static set of parallel strings. [Einstein theorem and matching coefficients](../../../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) · [Earlier native modular-orbit controls](../README.md)

Two opposite directions already distinguish the existing parallel-sheet rods from the response required by this Weyl source. Calibrate to equal rest energy. In the sheet geometry both axial null curvatures vanish. The independently calculated Weyl null stresses are 0.051620 and 0.058321. At rapidity 0.32 the mean of the two axial observers' stress readings exceeds the rest reading by 14.13%, while boosts along the sheets preserve their measured proper-ball area exactly. The same test can use three fixed-volume area readings with momentum-sign effects removed by the average. [Derivation](PROOF.md#an-axial-area-test-distinguishes-the-existing-rod-geometry) · [Directional and free-signal records](axial_response.json)

## Reproduction

From the repository root, use NumPy and SciPy in the project Python environment.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.measure
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.boosts
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.packet
python -m matched_einstein_limit.source_implication.modular_orbit.weyl_balls.axial_response
python -m pytest -q code/matched_einstein_limit/source_implication/modular_orbit/weyl_balls
~~~
