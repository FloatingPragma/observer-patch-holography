# Formation energy of the freely arriving neutral pair

Two observers record **the same two zero-ring-phase energy peaks, 49.0625 and 50.0625 exchange units**, when an arriving pair produces opposite signed loads on their neighboring ports. A single vacancy and a single excess on the same material both peak at 16, and a filled guide peaks at zero. Across all three original winding preparations and both input signs, more than 99.98% of the paired, zero-phase arrivals lie between 44 and 58 units (simultaneous 99% lower bounds). This particular mobile neutral pair is energetically costly compared with either charged constituent. [Original triggered release](../../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/neutral_ladder/release/README.md) · [Full clock and signed-charge counts](counts.npz) · [Decoded results](results.json)

## What reaches the detector

The charged signal first travels from its original source to port 160. Only 7.26–8.17 in ten million original starts reach the local donor operation. A full native pulse supplies the opposite charge on port 167. We reuse the five-owner release configuration in which the observer retires two original winding heads and releases both material loads on the resulting blank ring. The incoming source winding was independently prepared; the present detector does not read a transmitted winding qubit.

At the established free appointment of 256, owners 17 and 8 respectively read the ordinary blank edges 207–206 and 106–96. Their local echo records the difference of the two actual ternary loads. A signed difference of +2 or −2 locates the full vacancy–excess pair and gives its order. A separate direct acquisition reads the magnitude, undoes the write, and measures signed contrast on each of the three incident ring bonds on fresh trials. The two outer bonds yield squared differences of one and the middle bond yields four. For incident vacancy and winding zero, owner 17 infers a **post-arrival contrast energy of 47.9933–48.0068**, with simultaneous 99% coverage. All twelve sign/winding/receiver results lie within 47.9559–48.0439. This measurement uses local seam operations alone. [Native meter and charge propagation](apparatus.py)

The clock cohort measures formation energy **before** the signed coincidence. One setting also keeps zero phase under a one-port translation of the ring. Both zero-phase and all-phase settings use the same freely arriving material. Uniform delay labels at intervals of 2π/128 interfere at a Fourier output, so each clock produces a finite-width spectrum. The reference is a matched unit-filled ring of zero formation energy. A ring translation acts on the actual carried loads through thirteen native whole-load swaps. Coherent control and recombination of the translated paths and of the delay histories are supplied apparatus. The spatial projection can spread the charge around the ring, so the direct, controller-free joint-arrival experiment remains the evidence for free physical transport. The clock run determines the spectrum associated with this source and this material. [Acquisition and decoder](measure.py) · [Protocol](PROTOCOL.md)

| Original vacancy, winding zero, owner 17 | Signed pair arrivals | In the soft clock window | In 44–58 units |
|---|---:|---:|---:|
| Registered flight, all ring phases | 10,233,064 | 60 | 10,226,240 |
| Registered flight, zero ring phase | 291,412 | 0 | 291,394 |
| Exchange-only flight, zero ring phase | 61,504 | 39,551 | 0 |
| Compensator removed, registered single charge | 0 | — | — |

Each row starts with 2⁴⁶ simulated original launches. The third row changes only the ring's flight and clock-delay contrast after commissioning. Its neutral arrivals put 63.39–65.21% in the soft window, showing that this receiver can detect a soft neutral response. The last row records 587,468 ordinary charged arrivals at energy 16. A filled guide gives zero. For the registered zero-phase pair, the upper confidence limit below 1.6 units is 0.00468% of joint arrivals. On the shorter 1,024-path clock, every fine zero-phase peak is resolved only to a 0.125-unit bin, and soft readings sometimes occur. Even the finer 4,096-path all-phase clock has sidelobes. The 60 low readings do not imply a low eigenstate. An independent exact lower bound is given below. The declared one-tenth benchmark is relative to the measured single-charge guide standard. The original two-wall decay source belongs to a separate experiment. [Wall-source controls](../README.md)

The compressed raw bank contains every signed difference at every clock bin, spatial rejections, undelivered charged starts and original source refusals. Across all twelve fine zero-phase banks the reported peaks agree to one 0.03125-unit bin. Count intervals use exact binomial tails with 4,096 simultaneous comparisons at 99%. Banks are seeded multinomial draws from the numerical material amplitudes, with all original starts retained. The local three-contact cohorts each use 2³⁸ original starts. [Complete bank layout](results.json)

## Exact finite energy and decay checks

On each blank edge the free Hamiltonian is the positive exchange operator (I−swap)/2 plus 8 times the squared load difference. The reference material has zero energy. For a single vacancy and excess on this 14-port ring there are 28 neighboring ordered placements with contrast energy 48 and 154 separated placements with contrast energy 32. **Every state in this entire two-charge sector has formation energy at least 32**, including traveling states, states with zero spatial phase, and states having no definite relative separation. The projective energy of an arriving pair and a postcoincidence contrast read are different physical measurements.

Grouping the translations by the separation of the two charges gives a thirteen-dimensional integer matrix for twice the zero-phase generator. Exact real-root counts place one energy in each of [49.066660, 49.066670] and [50.062495, 50.062505]. Its eleven other energies are below 36. They are other relative-motion levels; two close peaks do not identify multiple particle species. Zero lattice phase and a rested particle in observer spacetime require a further spatial interpretation. [Exact certificate](certificate.py)

There is also a same-guide test of a proposed autonomous decay. Every freely prepared **bare single charge**, at any ring momentum, has energy between 16 and 18. A charged daughter accompanied by a freshly emitted vacancy–excess pair would have two like signed loads and one opposite load. Enumerating all 1,092 placements gives contrast energies 32, 48, 64 and 80 with multiplicities 140, 658, 280 and 14. Exchange positivity makes these source and product energies disjoint. Passive blank contacts also preserve load multiplicities. The experimental donor supplies the missing opposite stock through an addressed operation. A dressed, sufficiently energetic parent or a different neutral degree of freedom remains available to future tests. [Earlier inventory result](../../row_01_charge/README.md)

The official M10 OW-09 cell remains untested. Its criterion is “light charge-0 states released in decays” in both regimes. The upstream M10 cell separately describes the author’s Dirac construction with oscillation inputs and a declared lightest mass; it has no in-platform run. The remote observer has identified a mobile neutral material pair with substantial matched-host formation energy produced using a donor. Neither this source nor its measured spectrum qualifies for the row. [Upstream criterion and cell, head c531265eb](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#simplest-test-per-criterion)

Run from the repository root at base 3edcfa2af, or a descendant containing this package.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_09_neutral.flight_calorimetry
python -B -m "$p.measure" --write
python -B -m pytest --pyargs "$p.test_calorimetry" -q -W error
~~~
