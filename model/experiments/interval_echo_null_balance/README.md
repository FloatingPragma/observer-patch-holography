# Two physical gateways account for a patch’s energy change

On the registered Interval Echo sphere a patch measured the energy carried through each of two open contacts by the same two-load source used for the gauge-current experiment. Across one held interval it lost **0.4384 ± 0.0057** native energy units. Independently interrogating both gateways accounted for the loss to **0.0017 ± 0.0084**. Omitting one left **−0.1582 ± 0.0077** unaccounted for. Uncertainties are one standard error over simulated preparations. This supplies a finite matter-energy continuity check toward the Einstein Ward premise G3. [Held protocol](INTERVAL_PROTOCOL.md) · [Measurement and controls](interval.json) · [All count histograms](counts.npz)

## Source, service and energy meter

On the 20-cell screen, cell zero opens its ordinary port 1 to port 82 of cell six and port 3 to port 48 of cell four. Each of the three cells runs all its internal seam generators. The other services are gated closed during flight, including the degree-reference services. No site, edge, port or internal register is added to the carrier. At port one of cell zero the source prepares the normalized difference of two antipodal stars of adjacent two-load configurations. The registered rotation-selected pulse acts before release. The full 630-state binary two-load space of these three serviced cells propagates, including the part of the source outside the two selected triplets. The guide, gate windows, source state, phases, duration and readout preparations are supplied. [Current and source specification](../../candidate_laws/interval_echo/CURRENTS.md) · [Native Hamiltonian and measured releases](experiment.py)

Define the regional energy as the sum of the 30 internal seam generators. The two bridge terms store energy separately. On binary blank-memory inputs the native 99-state operator reduces on each seam to $h=8P_{\rm mixed}+P_-$, where $P_-$ selects the antisymmetric one-end-occupied state. One terminal reading of the patch's complete occupation word gives the number of mixed seams. In each of 30 other settings the patch phases one endpoint, applies a native half exchange, and counts the antisymmetric outcome on the chosen seam. The phase convention and the contraction onto the physical detector are checked directly against the native operator. Different settings and interrupted versus free flights use independent preparations. The sampled experiment uses the resulting local outcome distributions and records every word and click count; the routed implementation of the same word and interference readout was tested previously. [Meter checks](test_energy_flux.py) · [Routed four-cell meter](../../../../LAW/current-interval-clock/theorem-calibration/search/stress/README.md)

## Measuring the flux through one gate

Let $E$ denote the central regional operator. Its rate of change depends on each open bridge generator $h_b$ through $i[h_b,E]$. Every other cell's internal energy commutes with $E$, giving the exact scalar identity on all 630 pair states. Checking only a matrix identity would leave open which gateway a patch can actually distinguish. To read $b$ on one released state, interrupt the flight and apply only that gate for a known exposure $s$. Measure $E$ on new preparations at 19 exposures $s_j=2\pi j/19$. Every angle is physically forward. The binary gate has eigenvalues 0, 8 and 9 in native units. For $f_b(s)=\langle e^{ish_b}Ee^{-ish_b}\rangle$, its slope at zero follows from physical readouts as

$$ f_b'(0)=\sum_{j=0}^{18} f_b(s_j)\,{2\over19}\sum_{m=1}^{9}m\sin(ms_j). $$

There are no fitted coefficients. The two gateway measurements use disjoint preparations. Compare their sum with the central energy readings advanced or delayed on the original unperturbed flight. Both signs of the current pulses at ports 0 and 9 were tested at native times 1, 4 and 16. Across these 12 comparisons the RMS energy slope was **1.24189**, the RMS two-gateway discrepancy was **0.000000032**, and dropping gateway two raised the RMS discrepancy to **0.57094**. Halving the free difference step from 0.0001 to 0.00005 changed the estimated slope by at most 0.000000129. With both gates closed the central energy stayed constant to the receipt's precision. With gates open, neglecting energy on the links themselves produced a change in total cell energy as large as 0.90048. Adding the stored link energies restored conservation. [Initial protocol](PROTOCOL.md) · [Individual directions, times and detector responses](results.json)

## Resolving the energy transfer with finite records

The fresh interval ran from 1.125 to 1.375 after a positive port-nine pulse. Use six Gauss–Legendre appointments to integrate each separately measured flux. Doubling the quadrature order changed either noiseless integral by at most $2.1\times10^{-8}$ energy units. The simulation samples the exact count distributions for 524,288 preparations per setting, 230 energy readouts with 31 settings each, totalling 3,738,173,440 equivalent independent preparations. The deposited counts include zero-occupancy words and quiet seams. The two-gateway residual is consistent with zero at the reported resolution. The one-gateway control differs from zero by 20.5 standard errors. [Sampling code and uncertainty propagation](interval.py) · [Frozen sample and results](interval.json)

| Observable in native energy units | Expected | Sampled ± one SE |
| --- | ---: | ---: |
| Change in central energy | −0.432427 | −0.438417 ± 0.005721 |
| Inflow through port 1 | −0.277681 | −0.280202 ± 0.005194 |
| Inflow through port 3 | −0.154747 | −0.159911 ± 0.003372 |
| Change minus both measured inflows | 0 within integration error | +0.001695 ± 0.008431 |
| Change minus port 1 alone | −0.154747 | −0.158216 ± 0.007727 |

The two directions here are physical intercell contacts. They have no independently measured null vectors or pressures attached. A covariant matter Ward identity on the reconstructed observer geometry, common cap normalization, and curvature measured on the same state remain open. The earlier angular counterexample to identifying gauge currents directly with null stress still applies. [Einstein hypotheses G1–G6](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) · [Earlier angular test](../../../../LAW/current-interval-clock/theorem-calibration/search/ledger/README.md)

From the RER root the following commands replay the two measurements and test the operator, meter and archive controls.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_null_balance
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_null_balance.interval
python -m pytest -q model/experiments/interval_echo_null_balance
```
