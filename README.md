# Observer Patch Holography (OPH)

> Observer Patch Holography starts from one claim: no observer sees the whole world at once. Each observer accesses only a local patch, and neighboring patches must agree on their overlap. OPH asks how much physics can be reconstructed from that starting point once the full axiom and branch ledger is made explicit.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

OPH is a reconstruction program for fundamental physics. Spacetime, gauge structure, particles, records, and observer synchronization are treated as consequences of the OPH package rooted in overlap consistency on a finite holographic screen, together with the explicit branch premises stated in the papers.

## What OPH Delivers

- A finite-resolution theorem package for observer patches, collars, overlap repair, higher gauge structure, records, and checkpoint/restoration.
- A conditional route to Lorentzian geometry, modular time, Jacobson-type Einstein dynamics, and de Sitter static-patch cosmology on the extracted prime geometric subnet; the Einstein branch uses fixed-cap stationarity, the null-surface bridge, and the separate bounded-interval projective branch, while the remaining UV/BW scaffold is geometric cap-pair realization on that subnet plus ordered cut-pair rigidity, with the eventual fixed-local-collar modular-transport common floor as the smallest lower blocker.
- A conditional compact gauge route in the bosonic branch to the realized Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6`, under the transportable-sector reconstruction premises and MAR, together with the exact hypercharge lattice and the realized counting chain `N_g = 3`, `N_c = 3`.
- A particle program with exact structural massless carriers, a forward-emitted Phase II electroweak calibration branch with a closed target-free public `W/Z` theorem surface plus a compare-only exact frozen pair, a quantitative Higgs/top stage, exact non-hadron mass sidecars, and explicit continuation lanes for flavor and hadrons.
- A concrete screen-microphysics architecture that puts measurement, records, and observers inside the physics.

## Local Unification Surface

OPH places a local unification surface around the calibrated local UV input. The same `P`-driven scale carries the electroweak boson and Higgs lane together with the gravity-side entropy lane, while the Lorentz branch supplies the invariant causal speed and the local readout package supplies the SI display.
On the current public constants surface, `hbar` and `k_B` remain part of that downstream familiar-unit readout rather than standalone OPH-emitted constants.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260407" alt="OPH unification diagram" width="92%">
  </a>
</p>

Constants, theorem chains, and open proof fronts for this surface are tracked in [extra/OPH_PHYSICS_CONSTANTS.md](extra/OPH_PHYSICS_CONSTANTS.md).

**Current theorem stack and open fronts**

<p align="center">
  <a href="assets/prediction-chain.svg?v=20260412" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg?v=20260412" alt="Current OPH theorem stack and open proof fronts" width="92%">
  </a>
</p>

<p align="center"><sub>The current OPH stack from axioms to relativity, gauge structure, particles, observers, and the remaining open proof fronts. Click to open the full SVG.</sub></p>

## Particle Highlights

### Theorem-grade and structural hits

- Exact structural zeros for the photon, gluons, and graviton.
- Electroweak output on the D10 calibration branch, with closed target-free public `W/Z` rows and an exact frozen compare-only pair
  `W = 80.377 GeV`, `Z = 91.18797809193725 GeV`.
- A quantitative Higgs/top stage downstream of the electroweak core, with a closed one-scalar forward seed carrying the public rows
  `H = 125.218922 GeV`, `t = 172.388646 GeV`.

### Exact non-hadron output surface

| Lane | Exact output(s) | Status note |
| --- | --- | --- |
| Structural carriers | `m_photon = m_gluon = m_graviton = 0` | theorem-grade structural exactness |
| Electroweak sidecar | `W = 80.377 GeV`, `Z = 91.18797809193725 GeV` | exact frozen repair surface |
| Higgs/top exact sidecar | `(H, t) = (125.1995304097179, 172.3523553288311) GeV` | exact compare-only inverse slice on the same Higgs/top Jacobian |
| Charged witness | `(e, mu, tau) = (0.00051099895, 0.1056583755, 1.7769324651340912) GeV` | exact same-family witness |
| Quark witness | `(u, d, s, c, b, t) = (0.00216, 0.00470, 0.0935, 1.273, 4.183, 172.3523553288311) GeV` | exact same-family witness |
| Neutrino theorem branch | `(m1, m2, m3) = (0.017454720257976796, 0.019481987935919015, 0.05307522145074924) eV` with emitted weighted-cycle `Δm21²`, `Δm31²`, `Δm32²` | theorem-grade weighted-cycle bridge-rigid absolute family |

Public Higgs/top rows are carried by the closed one-scalar forward seed. The exact inverse pair above is a compare-only sidecar on the same Jacobian and does not replace the public forward branch.

**Particle derivation stack**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="OPH particle derivation stack" width="78%">
  </a>
</p>

<p align="center"><sub>A compact view of the particle lane. Click to open the full SVG.</sub></p>

### Continuation-grade particle successes

- The quark continuation lane emits public rows for `u`, `d`, `s`, `c`, and `b` on the selected continuation branch.
- The neutrino weighted-cycle branch reaches the observed PMNS and hierarchy regime with
  `theta12 = 34.2259°`, `theta23 = 49.7228°`, `theta13 = 8.68636°`, `delta = 305.581°`,
  and `Δm21² / Δm32² = 0.03072111`.
- The exact non-hadron surface is collected in
  [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

## Papers

- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: synthesis paper for the whole OPH stack.
- **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: SM/GR derivation paper for the recovered core.
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle derivation, exact-hit surface, and continuation map.
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
