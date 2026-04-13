# Observer Patch Holography (OPH)

> Observer Patch Holography starts from one claim: no observer sees the whole world at once. Each observer accesses only a local patch, and neighboring patches must agree on their overlap. OPH asks how much physics can be reconstructed from that starting point once the full axiom and branch ledger is made explicit.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

OPH is a reconstruction program for fundamental physics. Spacetime, gauge structure, particles, records, and observer synchronization are treated as consequences of the OPH package rooted in overlap consistency on a finite holographic screen, together with the explicit branch premises stated in the papers.

## What OPH Delivers

OPH is unusual because it tries to recover the shape of the world before it fits the numbers. At the structural level, it predicts the exact shape of the universe we appear to inhabit: a `3+1D` Lorentzian spacetime, de Sitter static-patch cosmology on the gravity side, and the Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6` with the exact hypercharge lattice and the counting chain `N_g = 3`, `N_c = 3`.

The quantitative side is deliberately small. OPH uses only two calibrated inputs: the screen-pixel scale `P = a_cell / l_P^2` and the total screen capacity `N_scr = log dim H_tot`, inferred from the cosmological constant. From that two-constant surface, OPH makes concrete numerical predictions for couplings, masses, and gravity-facing quantities instead of fitting each sector independently.

- A finite-resolution theorem package for observer patches, collars, overlap repair, higher gauge structure, records, and checkpoint/restoration.
- A conditional route to Lorentzian geometry, modular time, Jacobson-type Einstein dynamics, and de Sitter static-patch cosmology on the extracted prime geometric subnet; the Einstein branch uses fixed-cap stationarity, the null-surface bridge, and the separate bounded-interval projective branch, while the remaining UV/BW scaffold is geometric cap-pair realization on that subnet plus ordered cut-pair rigidity, with the eventual fixed-local-collar modular-transport common floor as the smallest lower blocker.
- A conditional compact gauge route in the bosonic branch to the realized Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6`, under the transportable-sector reconstruction premises and MAR, together with the exact hypercharge lattice and the realized counting chain `N_g = 3`, `N_c = 3`.
- A particle program with exact structural massless carriers, a forward-emitted Phase II electroweak calibration branch with a closed target-free public `W/Z` theorem surface plus a compare-only exact frozen pair, a quantitative Higgs/top stage, an exact selected-class quark closure with explicit exact forward Yukawas, exact non-hadron mass surfaces, and explicit continuation lanes where theorem boundaries remain open.
- A concrete screen-microphysics architecture that puts measurement, records, and observers inside the physics.

### Precise Derivations

This table focuses on OPH outputs with a clean one-number comparison against the latest official
PDG or NIST reference values. Structural results such as the `3+1D` Lorentz branch, the Standard
Model gauge quotient `SU(3) x SU(2) x U(1) / Z_6`, the exact hypercharge lattice, and the counting
chain `N_g = 3`, `N_c = 3` are stated in the papers and are not repeated here.

| Constant or particle | Standard shorthand | OPH predicted value | Latest PDG / NIST value | Agreement | Measurement source |
| --- | --- | --- | --- | --- | --- |
| Newtonian constant of gravitation | `G` | `6.674299995910528 x 10^-11 m^3 kg^-1 s^-2` | `6.67430(15) x 10^-11 m^3 kg^-1 s^-2` | `0.00003σ` | [NIST 2022](https://physics.nist.gov/cuu/pdf/wall_2022.pdf) |
| Speed of light in vacuum | `c` | `299792458 m s^-1` | `299792458 m s^-1 (exact)` | `100% match` | [NIST 2022](https://physics.nist.gov/cuu/pdf/wall_2022.pdf) |
| Inverse fine-structure constant | `alpha^-1(0)` | `137.035999177` | `137.035999177(21)` | `100% match` | [NIST 2022](https://physics.nist.gov/cgi-bin/cuu/Value?eqalphinv=) |
| Photon mass | `m_gamma` | `0 eV` | `< 1 x 10^-18 eV` | `Below PDG bound` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| Gluon mass | `m_gluon` | `0 GeV` | `0 GeV` | `100% match` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| Graviton mass | `m_graviton` | `0 eV` | `< 1.76 x 10^-23 eV` | `Below PDG bound` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| W boson mass | `m_W` | `80.377000015 GeV` | `80.3692 ± 0.0133 GeV` | `0.59σ, not yet fully derived` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| Z boson mass | `m_Z` | `91.187978078 GeV` | `91.1880 ± 0.0020 GeV` | `0.01σ, not yet fully derived` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| Higgs boson mass | `m_H` | `125.218922060 GeV` | `125.20 ± 0.11 GeV` | `0.17σ, not yet fully derived` | [PDG 2025 gauge/Higgs](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf) |
| Top quark mass | `m_t` | `172.388645595 GeV` | `172.56 ± 0.31 GeV` | `0.55σ, not yet fully derived` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |
| Bottom quark mass | `m_b(m_b)` | `4.183 GeV` | `4.183 ± 0.007 GeV` | `100% match` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |
| Charm quark mass | `m_c(m_c)` | `1.273 GeV` | `1.2730 ± 0.0046 GeV` | `100% match` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |
| Strange quark mass | `m_s(2 GeV)` | `93.5 MeV` | `93.5 ± 0.8 MeV` | `100% match` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |
| Down quark mass | `m_d(2 GeV)` | `4.70 MeV` | `4.70 ± 0.07 MeV` | `100% match` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |
| Up quark mass | `m_u(2 GeV)` | `2.16 MeV` | `2.16 ± 0.07 MeV` | `100% match` | [PDG 2025 quarks](https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf) |

Agreement is reported as a sigma distance where PDG or NIST quotes a one-standard-deviation
uncertainty. For exact definitions, exact digit matches, or published upper bounds, the table reports
`100% match` or `Below PDG bound` instead.

`W/Z` sit on the Phase II calibration lane, and Higgs/top sit on the quantitative D11 branch. Those rows are therefore not yet fully end-to-end derived in the same sense as the structural zeros, the exact fine-structure row, or the selected-class exact quark rows.

For the quark rows, PDG uses its standard quark-mass conventions: `u`, `d`, and `s` at `2 GeV`,
`c` and `b` in the `MS` scheme at their own mass scale, and the direct top-mass summary for `t`.
The papers also contain the structural Standard Model derivations listed above and a theorem-grade
neutrino family, which are not included in this table because they do not have a single direct
PDG/NIST one-number comparison row.

## Local Unification Surface

OPH places a local unification surface around the calibrated local UV input. The same `P`-driven scale carries the electroweak boson and Higgs lane together with the gravity-side entropy lane, while the Lorentz branch supplies the invariant causal speed and the local readout package supplies the SI display.
On the public constants surface, `hbar` and `k_B` remain part of that downstream familiar-unit readout rather than standalone OPH-emitted constants.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260407" alt="OPH unification diagram" width="92%">
  </a>
</p>

Constants, theorem chains, and open proof fronts for this surface are tracked in [extra/OPH_PHYSICS_CONSTANTS.md](extra/OPH_PHYSICS_CONSTANTS.md).

**Theorem stack and open fronts**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="OPH theorem stack and open proof fronts" width="92%">
  </a>
</p>

<p align="center"><sub>The OPH stack from axioms to relativity, gauge structure, particles, observers, and the open proof fronts. Click to open the full SVG.</sub></p>

**Particle derivation stack**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="OPH particle derivation stack" width="78%">
  </a>
</p>

<p align="center"><sub>A compact view of the particle lane. Click to open the full SVG.</sub></p>

## Papers

- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: synthesis paper for the whole OPH stack.
- **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: SM/GR derivation paper for the recovered core.
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle derivation, exact-hit surface, and theorem-boundary map.
- **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**: fixed-point, repair, and consensus formulation.
- **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)**: finite screen architecture, records, and observer machinery.

## More

- **Website:** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Theory explainer:** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Simulation-theory explainer:** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Book:** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Guided study app:** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions and detailed explanations:** OPH Sage on [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage), or [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **Lab:** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Common objections:** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **IBM Quantum note:** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Repository Guide

- **[`paper/`](paper):** PDFs, LaTeX sources, and release metadata.
- **[`book/`](book):** OPH Book source.
- **[`code/`](code):** computational material, particle outputs, and experiments.
- **[`assets/`](assets):** public diagrams and figures.
- **[`extra/`](extra):** maintained public notes such as objections, experimental write-ups, and selected supporting essays.

## OPH and the Sciences

<p align="center">
  <a href="assets/oph_science_overlap_map_poster.png" target="_blank" rel="noopener noreferrer">
    <img src="assets/oph_science_overlap_map.svg" alt="A map of the sciences OPH overlaps with, from large domains to subdomains to concrete OPH application areas." width="100%">
  </a>
</p>

<p align="center"><sub>A domain -> subdomain -> OPH-area map spanning mathematics, computer science, information and inference, complex systems, theoretical physics, quantum information, and measurement foundations. Click to open the full poster PNG.</sub></p>
