# Which boundary measurement cancels the regional entropy variation?

Two conserved load-two excitations propagated under the registered contact Hamiltonian inside a monitored patch. An observer at its exterior contacts could measure how many sharp boundary words were accessible. A coefficient trained on one matter variation was **0.694068 nats per accessible cut bit**. It canceled the interior entropy slope to within **0.02583 nats per unit probability transfer** on 54 held variations (protocol margin 0.05). An occupied-contact count calibrated on the same variation missed by up to 1.10218 nats. The positive result identifies a finite, operational boundary-information term. Relating it to independently measured proper area and establishing autonomous stationarity remain open G4 premises.

## Comparison with an independently surveyed boundary

For the two-appointment ball on the one-cell source, every single-station marker reading stays fixed under a held change of preparation. Joint rim flags resolve a change of −0.45815 nats in the exterior completion information. A timing-selected aperture with overlapping panels detects the sign and predicts −0.26824 nats after calibration on a different transfer. On the two-cell source, all twelve directional launch and return probability tables can be held fixed while the interior entropy has a strictly positive first variation. These experiments identify the completion response in shared records and test two specific area reconstructions. [Timed ball and angular apertures](CLOCK_BALL.md) · [Four-station comparison](area_survey/README.md)

## Experiment

The nine enrolled regions consist of one, two and three nested cells on each 20-, 80- and 320-cell federation. They have 12, 24 and 36 rails, 3, 4 and 5 exterior-facing contacts, and 66, 276 and 630 possible sharp two-load words. The full 99-state local Hamiltonian leaves the binary load-zero/load-two sector invariant. On unequal occupied and vacant endpoints it has diagonal energy 32.5 and swap amplitude −0.5 at the registered contrast coefficient. The regional generator sums every internal contact and retains the full federation's on-site degree correction. The external interaction gates are parked during each free interval. Monitoring the complete word, closing the gate and preparing each sharp input are apparatus operations.

A cut endpoint is read directly against an exterior vacancy. The bulk observer routes an interior rail to the cut using native contact swaps, reads the pointer and reverses the route. The code exhaustively verifies recovery of every sharp pair word and restoration of its preparation on every region. The exterior entropy meter uses only the directly read cut word and its own occupation count. Another fresh-copy meter conditions its boundary entropy separately on every interior word, providing a check on the independence of the exterior-only estimator.

The sharp reads are interleaved with 0.5 native exposure units. The probability of word j at the next read from word i is the squared magnitude of the corresponding native propagator element. The measured kernel is connected, symmetric and aperiodic. These properties select the uniform two-load ensemble for this monitoring schedule. The experiment computes the finite endpoint distribution from the reversible spectrum after 268,435,456 read intervals, then draws 4,194,304 terminal outcomes per independent meter. It starts the interior acquisition with both loads inside and the exterior acquisition with both loads on the cut. These are simulated terminal counts after spectral propagation. The archived spectra bound total-variation distance from equilibrium below 10⁻¹⁰ from any sharp input within 63,329 intervals, even for the slower biased controller.

## Held readings at fixed mean energy and cut load

The six directions per region include uniform transfers between each pair of cut-load classes, a localized redistribution between two interior words, a mean-energy-matched transfer, and a transfer preserving both mean native energy and mean occupied cut load. For the last direction, a localized interior-only pair transfer cancels the mean energy change. The native generator fixes this correction before evaluating either entropy reading. The geometry, total load and preparation normalization also stay fixed. Slopes are per unit of probability moved, evaluated about the prepared equilibrium. They are derivatives of the recorded interior distribution along physically preparable mixture variations.

| Region | Interior entropy slope across hosts, nats | Exterior-only capacity slope, bits | Largest absolute residual, nats |
| --- | ---: | ---: | ---: |
| One cell | +0.409 to +0.413 | −0.5944 | 0.00373 |
| Two cells | +0.432 to +0.476 | −0.6881 to −0.6215 | 0.00222 |
| Three cells | +0.420 to +0.448 | −0.6456 to −0.6082 | 0.00265 |

The occupied-cut meter has zero slope in every entry. The conditional exterior record still responds because at the same mean load the preparation can redistribute event weight between zero-, one- and two-occupied boundary configurations. The two independent entropy meter coefficients are 0.694068 for the exterior-only meter and 0.694071 for the joint-access check. Both fit only the zero-to-one transfer on fresh calibration counts from the single-cell region of the smallest board.

## Why the cancellation holds, and how it breaks

Let m count the exterior-facing rails, k their occupied count, and i denote the full interior word. At equilibrium the sharp word pairs are equiprobable. Each interior word admits exactly dₖ = binom(m,k) distinct cut words. Thus p(i) = dₖ / binom(M,2), where M is the number of enrolled rails, and the externally accessible conditional capacity is Cₖ = log₂ dₖ. The exact bulk entropy gradient is −ln p(i), up to a constant that drops out of normalized variations. For every such interior variation about this equilibrium,

$$
\delta H(I)+(\ln 2)\,\delta\langle C_k\rangle=0.
$$

The coefficient approaches the bits-to-nats conversion. The finite cancellation follows from uniform sharp words and the entropy chain rule. Physical content comes from the actual support, propagation, native contacts, addressable exterior records, and selection of that uniform law by the supplied monitoring. Outside this preparation family, calculating a cut entropy from the marginal alone can lose correlations with the interior.

Changing just the monitored acceptance rule gives a counterexample. When a proposed sharp word moves more loads outward, the controller accepts with probability 2^(−Δk). A rejected read returns to the previously resolved word by sharp swaps. A round trip through a common two-load reference certifies native re-preparation for every pair word, and the 99-state local gate certifies the fair bits for the decision. All inward moves are accepted. Detailed balance gives a unique equilibrium weight proportional to 2^(−k) for each sharp pair, with the same local Hamiltonian and exterior words as before. The frozen unbiased coefficient now misses by **1.39151 nats**, close to 2 ln 2 for a transfer of two loads. The signed residual predicted by the controller is ln 2 times the change in expected k. That prediction agrees within 0.02288 nats across all held directions. In particular the controller preserves first-order cancellation along the stricter directions that keep mean cut load fixed. With monitoring removed altogether, a pure preparation under the finite Hamiltonian retains unit purity. The sharp measurement entropy can grow even during a single unmonitored interval.

These nested contact regions test an operational finite equilibrium across more states and host resolutions. Their measured information capacity is nonlinear in occupied-contact count and sensitive to the preparation. Closing G4 still calls for an independently reconstructed proper volume and area, an area coupling in physical units, and material dynamics that select stationarity under the relevant variations.

## Reproduction

Run from the repository root. The JSON contains every fresh-copy count vector, the exterior-only word labels and pooled counts, source indices, matched energy variations, spectra-based convergence bounds and held residuals.

```bash
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_einstein_station.equilibrium_refinement.measure --output model/experiments/interval_echo_einstein_station/equilibrium_refinement/results.json
.venv/bin/python -m unittest model.experiments.interval_echo_einstein_station.equilibrium_refinement.test_measure
```
