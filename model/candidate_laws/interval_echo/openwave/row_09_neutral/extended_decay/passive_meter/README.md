# A remote material dimer stores energy from the decaying core

A cold calorimeter in owner six gained one exchange quantum in 2,571,093 of 4,194,304 simulated starts after a pulse in owner zero. The pulse and the read were 32 exchange times apart. On **1,930,863 starts (46.04%)** the calorimeter gained energy and the source observer read one of the two surviving lower heavy levels. The corresponding energy gain in the idle control was 265,442 starts (6.33%). The largest heating count in any control at this appointment was 693,560. All counts include every launched trial. [Held protocol](PROTOCOL.md) · [Recorded outcomes](records.json)

The contrast between pulsed, idle and severed preparations tests energy transported through this commissioned system. Neutrality of the balanced receiving dimer and its isolated gap are fixed before flight. A charge measurement on the outgoing disturbance in transit, joined to a light-energy record before absorption, would distinguish a neutral carrier from energy passing through the material. A positive terminal record alone does not identify what carried the energy.

The charged stocks 00 and 11 calibrate the two ports of the detector. The balanced orientations 01 and 10 give the identical neutral read, including their phases. The empty calorimeter is the lower exchange level of this independently prepared balanced dimer. Its excited level also has the neutral charge reading. An isolated two-quadrature clock on the native blank contact measures a gap of **1.000005 ± 0.000049** exchange units. The source costs 33.60603 units to prepare in its top isolated level. The two lower heavy readings cost 31.78359 and 32.67193 relative to its ground. [Physical circuit and calibrations](calorimeter.py) · [Five-dimer material and source preparation](../README.md)

## Heat at a single appointment

The first five dimers form the crossing material. The calorimeter occupies ports 77 and 81, respectively four and five ordinary seams from the nearest source rail. Its incident triangle touches the 76 rail of the end dimer. Fifteen active written contacts and twenty-three idle writes commission the six-dimer material from 234 conserved stock units. The first five dimers start in their connected cold stock, independently of the sixth. The pulse acts on the two source dimers. The coupled generator stays fixed from the start of free flight to the simultaneous level reads. The aperture has one terminal appointment. The independent stock preparation can cause heating on its own. That heating is measured in the idle and cut cohorts.

With the D–E link cut, the pulsed stocks gave 693,006 hot readings and idle stocks gave 693,560. Of the former, 451,609 coincided with a lower core. Contact between the initial E stock and the cold detector produces a substantial background even with the upstream route cut. Exact generator factorization across D–E makes the remote marginal independent of the source pulse on these same inputs, including their initial correlations. The final detector-link cut left the meter cold on every trial at every wait. Removing the extreme-contrast reversal on the same prepared stocks eliminated the heating. At wait 32 the pulsed joint frequency exceeded the greatest total control frequency by **29.50 ± 0.22 percentage points** (simultaneous 99% counting bounds). This bounds a difference between cohorts. Each separate material click has the same readout as a background click. [Six appointments, both source conditions and every control](records.json)

## Energy flow into the receiver

A second experiment sampled the work done on the meter through the E–F contact. In the orientation basis let e and f be the first-rail loads at E and F. The stored energy above the cold reading is D=(1−X_F)/2. Only the weak contact V=−(1−e)f contributes to i[H,D], which equals (1−e)Y_F/2. A simultaneous E rail read and F quadrature read therefore return −1/2, 0 or +1/2 units of power. Fresh cohorts measured all 1,024 midpoint appointments up to time 32 with 65,536 simulated starts each. Their integrated current gave **0.61177 ± 0.00918** detector energy units with the pulse and **0.06271 ± 0.00918** with the idle sham. The directly measured endpoint energies were 0.61300 and 0.06329. Under the upstream cut, the integrated currents were 0.16660 and 0.16655, compared with endpoint energies 0.16523 and 0.16536. The four current intervals cover the discrepancies even before allowing for endpoint counting uncertainty.

The common current half-width includes a simultaneous 99% bound on sampling at all four arms and the composite midpoint error T Δt² ‖[H,[H,i[H,D]]]‖/24. The endpoint readouts have their own simultaneous 99% counting intervals. Local contact power supplies an independent measurement of how much energy entered the meter through ordinary evolution. The declared static Hamiltonian and the exact 64-word sector bound this conclusion. The apparatus selects its gap and geometry. The result establishes passive reception in a finite material. A free daughter’s mass, weak interaction and identity between different materials require additional records. [Power protocol](POWER.md) · [Power circuit](power.py) · [All power counts](power_counts.json)

## How much of the signal belongs to the calorimeter?

Two spectator contacts on the receiving dimer set its energy gap. Both are native written memories coupled to fixed exterior units. Equal contact costs reproduce the original meter. Raising either cost by one exchange unit gives the same isolated Ramsey gap, $\sqrt{2}$, for opposite bias directions. All mounts use the same five-dimer input, parent pulse, source-to-receiver contact and read time. The meter starts in its own calibrated ground state. The table gives computed probabilities over **all launches** at 32 exchange times. The settings, controls, and comparisons were specified locally before computing these flights. The performance of the original meter was known. [Design](DETUNING.md) · [Complete joint laws and Ramsey quadratures](detuning.json)

| Receiver contact costs | Pulse excess in hot probability over idle | Lower heavy core and hot meter after pulse | Largest total hot probability in any control |
| --- | ---: | ---: | ---: |
| 32, 32 | 54.9657 percentage points | 46.0134% | 16.5330% |
| 33, 32 | 0.2587 percentage points | 1.2012% | 4.1444% |
| 32, 33 | 6.1521 percentage points | 5.7112% | 14.2204% |

The upstream cut erased the pulse excess for every meter. The last bias direction clears the one-point marginal contrast. Only the equal-cost meter clears the stronger comparison between its joint yield and every control's total heating. The detuned meters consequently give no second calibrated gap at which that comparison establishes excess joint events.

The difference between the two equal-gap meters tests an interpretation of the measured excess. The capture interaction is $-P_E\otimes n_F$, where $P_E$ is the occupation of the incident source rail and $n_F$ the meter's first-rail occupation. In weak absorption, the leading excitation probability from a cold meter is proportional to
$|\langle 1|n_F|0\rangle|^2\,\|\int_0^t e^{i\Delta s} P_E(s)|\psi\rangle\,ds\|^2$.
Here $P_E(s)$ evolves under the five-dimer Hamiltonian with the meter disconnected. Both bias directions have the same gap $\Delta$ and the same squared matrix element $1/8$. The leading weak-absorption prediction is therefore identical for each of the two receivers on any common input. Their pulse excesses differ by 5.8934 points at the full native contact. Detector orientation and its interaction with the source affect the reading beyond this weak-absorption description. Both meter levels stay in the same one-load sector. That charge accounting belongs to the receiver. The energy and electric charge of a freely emitted daughter remain unresolved. [Native generator and all flight settings](detuning.py) · [Direct propagation and matrix-element checks](test_detuning.py)

Run from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_09_neutral.extended_decay.passive_meter
python -Bm "$p.calorimeter"
python -Bm "$p.power"
python -Bm "$p.detuning"
python -Bm pytest -q --pyargs "$p.test_meter" "$p.test_detuning"
~~~

Each of the three measurement modules accepts `--acquire PATH` to generate its complete simulation record.
