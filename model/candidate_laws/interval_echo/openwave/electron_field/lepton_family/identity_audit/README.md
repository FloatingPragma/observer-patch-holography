# Can an observer identify the three guided charges?

A receiver on the seventh cell reads the signed charge and reconstructs the same outcome by counting inward written heads. An undisturbed uniform rest state clicks there with probability **1/12** for each of the three register characters, for either sign, at every native wait. Repeated measurements can distinguish their motion. The same records admit a realization by three spinless orbital amplitudes. The physical electron, muon and tau assignments require additional observables. [Complete pointer records](results.json) · [Source and charge-to-flux map](../dirac_read/README.md)

The apparatus uses the twelve-port route and the same vertex and auxiliary charge addresses as the first-order incidence comparison. Each trial prepares a chosen register character. The receiver receives no character label. It projects the local charge at native times 1, 2, 4, 8 and 12, copying each binary outcome to its own record. Every miss projects onto the complementary state before the next interval. The receipt enumerates all 32 record words for each input and contains 64 sampled five-event trials per input and source setting. The Born probabilities, supplied native schedule and seeded sampling distinguish the exact instrument distribution from an individual simulated execution.

| Independently prepared source | Best single-word accuracy for three equiprobable characters | Light versus middle | Light versus heavy |
|---|---:|---:|---:|
| Uniform rest state | 34.4588% | 50.0088% | 51.6882% |
| Sharp charge at the first port | 41.0709% | 50.0553% | 61.6063% |

For three equal priors the optimum is one third of the sum, over all words, of the largest conditional word probability. The two-input accuracy is one half plus one half of the total-variation distance. The undisturbed rest wave is uniform at all waits. The first projective read can disturb it, so later clicks carry information about its motion. Giving all three characters the mean mass effect keeps the initial source, contact geometry, charge detector and cyclic symmetry, and reduces the three-way optimum to precisely one third.

## Spin and weak transitions on this path

For each character \(v_j\), the isometry \(L_j x=x\otimes v_j\) intertwines the full generator with the one-component incidence generator of mass \(f_j\), and intertwines every local charge projector with the same spatial projector. Products of evolution and pointer copying therefore reproduce every observation sequence in a direct sum of three spinless orbital branches, even outside the five sampled times. A two-state spin carried as a spectator also gives these same charge and flux probabilities.

Each rested energy is a simple level for a fixed sign. The path incidence has one uniform zero mode and its first nonzero singular value is \(2\sqrt{4}\sin(\pi/24)=0.522105\), above even the largest resting energy of 0.209694. The three distinct positive mass effects supply exactly one uniform rest state apiece. An energy-preserving intrinsic spin-½ doublet cannot fit inside any of these three one-dimensional rest eigenspaces. Likewise any antiunitary reversal preserving this Hamiltonian and fixed charge sign squares to +1 on a simple eigenline, because the second application conjugates its phase. Conjugation in the real orbital blocks exhibits such a +1 symmetry. A fermionic reversal with square −1 would need a partner or an action outside this prepared sector. An unconfined medium or a Hamiltonian with additional physical spin components requires a separate test.

The generator and both charge outcomes commute with every family-character projector. All heavy-to-light transitions under the audited free evolution and charge reads have zero amplitude, at arbitrary times. The distal detector registers a supplied local address-phase kick converting the heavy input into the light character with joint probability 1/12. The kick changes the mean native energy by −0.0174694 and supplies the conversion operation. A spontaneous, energy-accounted weak transition with tagged daughters is a different measurement. [Matched charged and neutral daughter experiment](../decay/README.md)

## Absolute formation energy

Assign an additional cost of one native unit to every one-charge state and keep the neutral reference at zero. This comparison preserves the entire charge instrument word by word. The propagator in each such sector gathers the same global phase on every trajectory. The modified rest ratios are **1 : 1.012407 : 1.209621**, while the largest numerical change over the two source settings, three characters and 32 words is below \(10^{-13}\). Measuring an absolute formation energy requires a reference or a work measurement that compares the compensated charged preparation with its vacuum. Conversion to keV additionally requires an independently calibrated physical unit.

Reproduce the receipt in the repository environment with `python -B -W error measure.py --output FILE` in this directory, using an unused output path. [Executable](measure.py) · [Measurement design](DESIGN.md)
