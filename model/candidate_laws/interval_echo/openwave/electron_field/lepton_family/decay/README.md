# Does a higher charged mode decay into a lower one?

**The three modes exchange energy with the written seam field.** At time π/4, a freely released excess prepared in level 2 reaches receiver level 0 with an excited neutral donor field in **19.7201%** of all preparations. Preparing level 1 gives **19.2397%**. The entrance contact supplies energy even to the lowest mode. Comparing the complete neutral field spectrum with every downward charged transition finds zero separated daughter states at the isolated parent's energy in this two-cell apparatus. [Measurement design](DESIGN.md) · [Source](measure.py) · [Numerical record](results.json)

## What was prepared and measured

The two written-field cells share the degree-preserving material seam 3–48. Before joining them, prepare one of three ideal charged energy eigenstates in the donor and the independent zero-energy field ground state in the receiver. The only mobile matter defect is the signed excess or vacancy. All 60 internal material links, 40 triangular field rings, and the seam then evolve with their archived native Hamiltonian. The receiver measures a charged energy level as the donor measures its entire neutral field spectrum. The joint dimensions are 2,211,780 for a vacancy and 734,704 for an excess. A cut seam keeps the charge in the donor and the initially neutral receiver in its ground state at all times. [Archived free Hamiltonian and gateway construction](../../records.zip) · [Independent preparations and transfer reads](../transfer/README.md)

| Charge | Parent → daughter | Charged daughter and positive field energy at π/4, all trials | Parent minus daughter energy | Closest positive donor field energy | Absolute mismatch |
| --- | --- | ---: | ---: | ---: | ---: |
| Vacancy | 1 → 0 | 0.100223% | 0.277087 | 0.288779 | 0.011692 |
| Vacancy | 2 → 0 | 0.129643% | 0.376479 | 0.412130 | 0.035651 |
| Vacancy | 2 → 1 | 0.000024604% | 0.099392 | 0.098528 | 0.000863 |
| Excess | 1 → 0 | 19.2397% | 0.280236 | 0.308444 | 0.028207 |
| Excess | 2 → 0 | 19.7201% | 0.469771 | 0.484752 | 0.014981 |
| Excess | 2 → 1 | 0.302663% | 0.189534 | 0.206407 | 0.016872 |

The energies are native generator units. Each row includes every positive neutral energy in the donor, with degenerate eigenspaces summed. The full record gives charge arrival outside the three selected levels, every neutral spectral weight at π/4, and the corresponding daughter and neutral probabilities at π/8 and π/2. The positive bare bridge matrix element squared at the closest level of each transition ranges from 9.41 × 10⁻⁷ to 0.03052. The shell comparison uses absolute tolerance 10⁻⁸. [Numerical record](results.json)

## Energy supplied by opening the seam

For a charged source with gateway weight *w*, the sudden opening adds mean energy 4*w* and produces open-Hamiltonian energy variance 32*w* − 16*w*². The bridge is the native positive contact with eigenvalues 0 and 8. The work ranges from 0.0580 to 0.1170 across the vacancy sources and from 1.3996 to 1.5139 across the excess sources. The lowest charged source itself produces a level-0 daughter accompanied by positive neutral field energy with probability **0.001954%** for a vacancy and **0.867909%** for an excess. Across two numerical environments, the total-energy drift stays below 8.7 × 10⁻¹⁴ native units and the norm error below 1.6 × 10⁻¹⁴. Ground-field amplitudes for interchanging source and daughter agree in both directions within 2.1 × 10⁻¹⁶. [Generator, controls and readout](measure.py) · [Numerical record](results.json)

The energy shell compares separated cells, the preparations before opening, and the products after isolation. The full open Hamiltonian conserves its energy throughout flight. The receiving charged level includes its confined orbital energy. Recoil into a travelling momentum state and neutral escape into a third cell lie outside this geometry. A strong connection can shift modes and mix many such products. Its finite-time coincidences therefore measure conversion following a quench. Establishing an autonomous decay would require a parent resonance coupled to energy-conserving outgoing products under continuous evolution, with release and escape demonstrated in a geometry large enough to distinguish return from a physical lifetime.

## What the neutral detector can identify

Upon heralding a receiver daughter at π/4, isolate the donor and let its field evolve. Each local ring term is a projector on the antisymmetric combination of its two flippable orientations. Summing the 20 ring expectations recovers the total neutral energy. For the vacancy 2 → 0 channel, the far-face click probability per original release changes from 1.0494 × 10⁻⁵ to 8.8655 × 10⁻⁶ during a further wait π/2. At that last wait, both vertex-disjoint near and far detectors click with joint probability 3.9932 × 10⁻⁸. In the excess sector five of the 20 faces, including the selected far face, have zero flippable pairs across the entire neutral basis. The participation of every face is recorded so this geometric limitation can be checked directly. These ideal clicks locate energy in the neutral field of one donor cell. They provide no neutral-particle multiplicity or species identification. [Detector geometry and time series](results.json)

The existing weak pair service has support on cells containing precisely two absolute loads. The tested cells contain 12 and 11 loads, or 12 and 13, interchanged when the mobile defect crosses. Zero basis states in either complete joint space are eligible. Neither charged energy-level selection nor either ring detector supplies an independent lepton-flavour read. A larger radiative continuum or a different charged composite remains outside this experiment. [Original species service and controls](../../records.zip) · [Support check](measure.py)

The experimental comparison includes μ⁻ → e⁻ + ν̄ₑ + νμ and the leptonic tau channels τ⁻ → ℓ⁻ + ν̄ℓ + ντ for ℓ = e or μ. The daughter energy and polarization distributions and the inverse weak process provide additional tests of the family assignment. [Particle Data Group, Muon Decay Parameters](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-muon-decay-params.pdf)

## Replay

Run from the repository root with `model/requirements.txt` installed and a fresh output path. The script loads the original charged sector, physical gateway and joint propagator directly from the repository's `electron_field/records.zip`.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -B model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/decay/measure.py \
  --output /tmp/written-family-decay.json
```
