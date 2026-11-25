# A charged excitation gives its energy to neutral matter across a cell boundary

A charged vacancy dimer in owner 0 can relax while a neutral dimer in owner 6 gains the same independently measured excitation energy. One fixed contact acts throughout the flight. At 64 native exchange times the observers jointly read **charged ground and neutral excited in 4,158,050 of 4,194,304 simulated preparations**. Exchanging the initial energy occupations gives 4,157,987 reverse transfers. At time 128 the excitation has largely returned. This apparatus realizes resonant energy delivery into prepared neutral matter. [Physical layout and all measurements](results.json) · [Executable](measure.py) · [Protocol](PROTOCOL.md)

## What the observers measure

A vacancy can move on the 0–1 dimer in owner 0. A second vacancy moves on the 74–82 dimer in owner 6, which also contains one stationary compensating excess at port 73. Each gateway has degree six and each outer dimer port has one native private reference. A symmetric extreme-contrast collar on the 1–82 seam couples the two dimers. Twenty-one surrounding contacts carry native zero-contrast marks. The two active dimer seams stay blank and all material settings stay fixed during the release. The signed terminal reads give charge −1 in the source owner and zero in the receiver. [Stored contact on every endpoint word](../../../../../experiments/interval_echo_consensus/signed_effects/COLLAR_STORAGE.md) · [Silent contact construction](../silent_seams/README.md)

Each regional observer first fills the opposite dimer and measures the local clock. A balanced vacancy is released at 32 delays in steps of 1/32. Return and complement counts give excitation gaps **40.00870** on the charged arm and **40.01214** on the neutral arm, against the same exact gap $\sqrt{1601}=40.012498\ldots$. The three-point cosine recurrence on the first differences gives the frequencies with no target frequency in the estimator. The two local energy-basis readouts use a native dimer exchange and a calibrated private phase. Each read includes all four possible energy pairs. [Calibration counts, frame rotations and decoder](measure.py)

| Free exposure | Charged ground and neutral excited, per launch | Simulated joint reads / 4,194,304 |
| --- | ---: | ---: |
| 0 | 0 | 0 |
| 32 | 0.502521665 | 2,106,523 |
| 64 | 0.991364297 | 4,158,050 |
| 128 | 0.000920122 | 3,701 |

On separate copies the observers read whether both vacancies occupy the central gateways. This measures the part of the interaction energy left after local calibration. At exposure 64 their local counts give a neutral energy gain of **39.99911**. Summing both observed local energies and the separately measured interaction term gives **64.99523**, against the conserved full-system mean 64.99500. The electrical charge of the parked excess is a material load assignment. Removing it changes the receiver charge to −1 and subtracts a constant 24 from its energy while preserving the transfer probabilities. An initially filled receiver without a mobile vacancy has exactly zero absorption in this apparatus. [Additional material controls](results.json)

## Why the energies match

Write $P=|n\rangle\langle n|$ for the local vacancy at the near gateway, and use the basis (far, near) for each vacancy. On the invariant nine-word blank collar the active dimer has exchange amplitude $-1/2$. On the stored crossing the 99-state native seam gives $32(1-\mathbf1_{11})$. Relative to the filled reference this contributes $32(P_A+P_B-P_A P_B)$. When the other dimer is filled, the single-arm Hamiltonians are
$$
H_A=\tfrac12 I+40P-X/2,\qquad H_B=H_A+24I .
$$
The 40 comprises 32 from the calibrated crossing and 8 from the private-reference energy. The 24 is the cost of the parked excess. The full coupled Hamiltonian is
$$
H=H_A\otimes I+I\otimes H_B-32P\otimes P .
$$
The program reconstructs these matrices from all possible native outputs of each physically addressed word. The filled-host energy convention is common to every preparation.

The charged levels are $(41\pm\sqrt{1601})/2$ and the neutral levels $(89\pm\sqrt{1601})/2$. Thus the two observed product channels both have independent energy **65 exactly**. The local projector has $|\langle e|P|g\rangle|=1/(2\sqrt{1601})$. On the equal-energy pair the native coupling is $-8/1601$ in the deposited real phase convention. The magnitude of the neutral energy-current matrix element for $J_B=i[H,I\otimes H_B]$ is $8/\sqrt{1601}$. This identifies a positive coupling on the independent energy shell as well as a finite-time regional energy transfer. [Formal energy identity and coupling](ResonantTransfer.lean)

Exchanging the two mobile vacancy positions commutes with the complete Hamiltonian and carries one product preparation to the other. Their *entire* full-energy distributions therefore agree. The antisymmetric combination is a sharp full-Hamiltonian eigenstate at 65, with probability one half in either preparation. A second eigenstate at 64.95040063 carries another 49.75587%. Interference between the stationary lines transfers energy between the two owners. The localized preparations are finite packets; selecting the 65 line alone gives a stationary state. An independent rational enclosure of the full four-word evolution gives $0.991364296635\le P(64)\le0.991364296636$. [Dyadic propagation bound](enclosure.py) · [Independent spectral replay](test_delivery.py)

## Which alternatives the records distinguish

Replacing the shared crossing with separate identically marked contacts to filled posts at 5 and 78 preserves both local calibration spectra. The joint transfer then has exactly zero probability. Adding a static contact between the receiving gateway and its filled post changes its gap to $\sqrt{5185}$; even with its recalibrated energy reader the transfer probability at 64 is $2.477\times10^{-7}$. Halving contrast during free evolution, with both resulting arms independently calibrated, gives a joint probability of 0.001486 at that same appointment. The native 99-state law produces both the local frequencies and the coupling. The geometry, collar preparation, compensated receiving matter, clock and terminal operations are supplied.

The neutral vacancy and its compensating excess exist before the charged source relaxes. A freely emitted neutral daughter in an initially empty receiver, an independently identified massive parent and a decay width require different observations. The independently moving, three-defect upper band on a blank host has a separate energy budget that excludes its two-body breakup into one vacancy and one neutral pair. [Blank-host daughter spectra](../daughter_shell/README.md)

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.resonant_delivery.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/resonant_delivery
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_11_spectrum/resonant_delivery/ResonantTransfer.lean)
```
