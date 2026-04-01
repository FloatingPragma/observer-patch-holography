# Observer Patch Holography (OPH)

> Observer Patch Holography starts from one claim: no observer sees the whole world at once. Each observer accesses only a local patch, and neighboring patches must agree on their overlap. OPH asks how much physics follows from that demand alone.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [website](https://floatingpragma.io/oph/) | [Paper 1: synthesis](paper/observers_are_all_you_need.pdf) | [Paper 2: SM/GR derivation](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) | [Paper 3: particle derivation](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) | [OPH Book](https://oph-book.floatingpragma.io) | [OPH Lab](https://oph-lab.floatingpragma.io) | [challenge](https://challenge.floatingpragma.io)

OPH is a reconstruction program for fundamental physics. Spacetime, gauge structure, particles, records, and observer synchronization are treated as consequences of overlap consistency on a finite holographic screen.

## What OPH Delivers

- A finite-resolution theorem package for observer patches, collars, overlap repair, higher gauge structure, records, and checkpoint/restoration.
- A conditional route to Lorentzian geometry, modular time, Jacobson-type Einstein dynamics, and de Sitter static-patch cosmology.
- A compact gauge route to the realized Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6`, together with the exact hypercharge lattice and the realized counting chain `N_g = 3`, `N_c = 3`.
- A particle program with exact structural massless carriers, theorem-grade electroweak closure, a quantitative Higgs/top stage, exact non-hadron mass sidecars, and explicit continuation lanes for flavor and hadrons.
- A concrete screen-microphysics architecture that puts measurement, records, and observers inside the physics.

**Overall theorem and derivation stack**

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="Overall OPH theorem and derivation stack" width="92%">
  </a>
</p>

<p align="center"><sub>The full OPH stack from axioms to relativity, gauge structure, particles, observers, and the remaining open fronts. Click to open the full SVG.</sub></p>

## Particle Highlights

### Theorem-grade and structural hits

- Exact structural zeros for the photon, gluons, and graviton.
- Electroweak closure on the D10 calibration chain, with target-free public `W/Z` rows and an exact frozen pair
  `W = 80.377 GeV`, `Z = 91.18797809193725 GeV`.
- A quantitative Higgs/top stage downstream of the electroweak core, with strong public forward rows
  `H = 125.218922 GeV`, `t = 172.388646 GeV`.

### Exact non-hadron output surface

| Lane | Exact output(s) | Status note |
| --- | --- | --- |
| Structural carriers | `m_photon = m_gluon = m_graviton = 0` | theorem-grade structural exactness |
| Electroweak sidecar | `W = 80.377 GeV`, `Z = 91.18797809193725 GeV` | exact frozen repair surface |
| Higgs sidecar | `H = 125.1995304097179 GeV` | exact compare-only inverse slice |
| Charged witness | `(e, mu, tau) = (0.00051099895, 0.1056583755, 1.7769324651340912) GeV` | exact same-family witness |
| Quark witness | `(u, d, s, c, b, t) = (0.00216, 0.00470, 0.0935, 1.273, 4.183, 172.3523553288311) GeV` | exact same-family witness |
| Neutrino adapter | `(m1, m2, m3) = (0.01745663295, 0.01948419960, 0.05308139066) eV` with exact representative `Δm21²`, `Δm31²`, `Δm32²` | exact compare-only adapter |

**Particle derivation stack**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="OPH particle derivation stack" width="78%">
  </a>
</p>

<p align="center"><sub>A compact view of the particle lane. Click to open the full SVG.</sub></p>

### Continuation-grade particle successes

- The quark continuation lane emits public rows for `u`, `d`, `s`, `c`, and `b` on the selected D12 sheet.
- The neutrino weighted-cycle branch reaches the observed PMNS and hierarchy regime with
  `theta12 = 34.2259°`, `theta23 = 49.7228°`, `theta13 = 8.68636°`, `delta = 305.581°`,
  and `Δm21² / Δm32² = 0.03072111`.
- The exact non-hadron surface is collected in
  [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

### Brief gap summary

The remaining particle gaps are narrow and explicit: charged leptons wait on `C_hat_e^{cand}` and `mu_phys(Y_e)`, quarks wait on `intrinsic_scale_law_D12` on the selected D12 mass ray, neutrinos wait on `C_nu`, and hadrons wait on the production backend bundle plus full production systematics.

## Papers

- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: synthesis paper for the whole OPH stack.
- **Paper 2. [Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: SM/GR derivation paper for the recovered core.
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle derivation, exact-hit surface, and continuation map.
- **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**: fixed-point, repair, and consensus formulation.
- **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)**: finite screen architecture, records, and observer machinery.

## More

- **Website:** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Theory explainer:** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Simulation-theory explainer:** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Book:** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Lab:** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Common objections:** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **IBM Quantum note:** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Repository Guide

- **[`paper/`](paper):** PDFs, LaTeX sources, and release metadata.
- **[`book/`](book):** OPH Book source.
- **[`code/`](code):** computational material, particle outputs, and experiments.
- **[`assets/`](assets):** public diagrams and figures.
- **[`extra/`](extra):** explainers, objections, and supporting notes.
