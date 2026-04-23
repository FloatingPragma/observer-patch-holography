# Observer Patch Holography (OPH)

> Observer Patch Holography starts from one claim: no observer sees the whole world at once. Each observer accesses only a local patch, and neighboring patches must agree on their overlap. OPH asks how much physics can be reconstructed from that starting point once the full axiom and branch ledger is made explicit.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io)

OPH is a reconstruction program for fundamental physics. Spacetime, gauge structure, particles, records, and observer synchronization are treated as consequences of the OPH package rooted in overlap consistency on a finite holographic screen, together with the explicit branch premises stated in the papers.

## Authority and Reading Rule

For recovered-core theorem status and claim tier, consult **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)** first. Lane-specific claim tier remains with the corresponding companion papers, including **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**, **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**, and **Paper 5. [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)**. This README, Paper 1, and the book are synchronized synthesis surfaces; they summarize and organize results but do not promote claim tier.

## What OPH Delivers

OPH is unusual because it tries to recover the shape of the world before it fits the numbers. At the structural level, it predicts the exact shape of the universe we appear to inhabit: a `3+1D` Lorentzian spacetime, de Sitter static-patch cosmology on the gravity side, and the Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6` with the exact hypercharge lattice and the counting chain `N_g = 3`, `N_c = 3`.

The quantitative side is deliberately small. OPH uses only two calibrated inputs: the screen-pixel scale `P = a_cell / l_P^2` and the total screen capacity `N_scr = log dim H_tot`, inferred from the cosmological constant. From that two-constant surface, OPH makes concrete numerical predictions for couplings, masses, and gravity-facing quantities instead of fitting each sector independently.
The exact numerical value of `P` comes from the calibration surface. Its proximity to the golden ratio is part of the quantitative explanation. `φ = (1 + sqrt(5)) / 2` is the exact self-similar balance point of the total/bulk/edge hierarchy, and a universe with durable records, structure, and dynamics must sit at a small equilibrium-breaking detuning away from that balance point. This belongs to the Phase II calibration layer, not the structural core.

- A finite-resolution theorem package for observer patches, collars, overlap repair, higher gauge structure, records, and checkpoint/restoration. On the genuinely noncentral branch this already includes the fixed-cutoff crossed-module defect hierarchy; the remaining open item there is the refinement-limit transportable-sector lift into the compact-gauge lane.
- A conditional route to Lorentzian geometry, modular time, Jacobson-type Einstein dynamics, and de Sitter static-patch cosmology on the extracted prime geometric subnet; local null data fix the Einstein branch only up to a metric term, so the cosmological-capacity lane fixes `Lambda` through the global screen-capacity identification. Local null data alone do not fix it. The Einstein branch uses fixed-cap stationarity, the null-surface bridge, and the separate bounded-interval projective branch. The UV/BW scaffold consists of geometric cap-pair realization on that subnet plus ordered cut-pair rigidity, with the eventual fixed-local-collar modular-transport common floor as the smallest lower blocker.
- A conditional compact gauge route in the bosonic branch to the realized Standard Model quotient `SU(3) x SU(2) x U(1) / Z_6`, under the transportable-sector reconstruction premises and MAR, together with the exact hypercharge lattice and the realized counting chain `N_g = 3`, `N_c = 3`.
- A particle program with exact structural massless carriers, a forward-emitted Phase II electroweak calibration branch with a closed target-free public `W/Z` theorem surface plus a compare-only exact frozen pair, an exact source-only Higgs/top split theorem on the declared D10/D11 running, matching, and threshold surface together with a compare-only exact inverse Higgs/top slice, an exact selected-class quark closure with explicit exact forward Yukawas, exact non-hadron mass surfaces, and explicit continuation lanes where theorem boundaries remain open.
- A concrete screen-microphysics architecture that puts measurement, records, and observers inside the physics.

### Precise Derivations

This condensed table keeps only OPH rows with either an exact match, a quoted sigma agreement, or
a clean upper-bound success against the PDG/NIST reference values used in the papers. Structural
results such as the `3+1D` Lorentz branch, the Standard Model gauge quotient
`SU(3) x SU(2) x U(1) / Z_6`, the exact hypercharge lattice, and the counting chain `N_g = 3`,
`N_c = 3` are stated in the papers and are not repeated here. The `W/Z/H` boson lane sits on the
Phase II calibration branch, so it is discussed in the papers but omitted from this quick
comparison.

| Quantity | Symbol | OPH | PDG/NIST | Δ |
| --- | --- | --- | --- | --- |
| Gravitational constant | G | 6.6742999959e-11 | 6.67430(15)e-11 | 0.00003σ |
| Speed of light | c | 299792458 | 299792458 (exact) | match |
| Fine-structure (inv) | α⁻¹(0) | 137.035999177 | 137.035999177(21) | match |
| Photon mass | m_γ | 0 eV | <1e-18 eV | below bound |
| Gluon mass | m_g | 0 GeV | 0 GeV | match |
| Graviton mass | m_grav | 0 eV | <1.76e-23 eV | below bound |

**Quark sector**

| Quark | Symbol | OPH | PDG | Δ |
| --- | --- | --- | --- | --- |
| Bottom | m_b(m_b) | 4.183 GeV | 4.183 ± 0.007 | match |
| Charm | m_c(m_c) | 1.273 GeV | 1.2730 ± 0.0046 | match |
| Strange | m_s(2 GeV) | 93.5 MeV | 93.5 ± 0.8 | match |
| Down | m_d(2 GeV) | 4.70 MeV | 4.70 ± 0.07 | match |
| Up | m_u(2 GeV) | 2.16 MeV | 2.16 ± 0.07 | match |

`Δ` reports the sigma distance where PDG or NIST quotes a one-standard-deviation uncertainty.
Otherwise it records `match` or `below bound`.

For the quark rows, PDG uses its standard quark-mass conventions: `u`, `d`, and `s` at `2 GeV`,
and `c` and `b` in the `MS` scheme at their own mass scale.
The papers also contain the structural Standard Model derivations listed above and a theorem-grade
neutrino family, which are not included in this table because they do not have a single direct
PDG/NIST one-number comparison row.
The public neutrino surface also includes theorem-grade physical Majorana phases on the
shared-basis weighted-cycle transport branch; see `code/particles/RESULTS_STATUS.md`.

The declared electroweak calibration surface also carries an exact source-only Higgs theorem with
`m_H = 125.1995304097179 GeV` and a companion top coordinate
`m_t = 172.3523553288312 GeV` on the same Jacobian surface.
At the precision quoted by PDG, the Higgs row lands on the 2025 Higgs average.
The exact public running-top row on the selected quark frame uses the PDG 2025 cross-section entry
`Q007TP4`.
The bridge to the auxiliary direct-top average `Q007TP = 172.56 ± 0.31 GeV` is open and tracked in
[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207).

Charged leptons sit on a split claim tier. The repo carries an exact same-family witness, the same-label `q_e` readback, a source-side determinant character defined for a fixed source multiplicity vector, a conditional determinant-line lift on physical charged data, and an algebraic mass readout from theorem-grade absolute charged scale. The public theorem lane does not emit a theorem-grade sector-isolated charged determinant exponent vector, and it does not identify the source-side determinant character with the physical charged determinant line. Public electron, muon, and tau masses are therefore not emitted from `P`.

## Local Unification Surface

The local unification surface is organized around the calibrated local UV input. The same `P`-driven scale carries the electroweak boson and Higgs lane together with the gravity-side entropy lane, while the Lorentz branch supplies the invariant causal speed and the local familiar-unit package supplies meters, seconds, GeV, and Kelvin. On the stated local extension surface, the lifted product presentation of the realized quotient branch gives `ellbar_shared = ellbar_SU(2) + ellbar_SU(3)`; the same D10 pixel law on that surface fixes `ellbar_shared = P/4`, and the gravity-side emitted row is `G_SI = c^3 a_cell / (hbar P)` relative to the declared microscopic datum `a_cell`.
The familiar-unit package on that surface is explicit: `L_loc = sqrt(a_cell) * Lhat(P)`,
`t_loc = sqrt(a_cell) * That(P) / c`, `E_loc = hbar c * Ehat(P) / sqrt(a_cell)`, and
`Theta_loc = hbar c * Thetahat(P) / (k_B sqrt(a_cell))`, with dimensionless branch outputs
`Lhat`, `That`, `Ehat`, `Thetahat`. So `a_cell` supplies the one local ruler, `c` is the
structural Lorentz output, and `hbar` and `k_B` are downstream familiar-unit display conventions
for GeV and Kelvin rather than standalone OPH-emitted constants.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260415" alt="OPH unification diagram" width="92%">
  </a>
</p>

Particle status surfaces for this repo live in [code/particles/RESULTS_STATUS.md](code/particles/RESULTS_STATUS.md) and [code/particles/EXACT_NONHADRON_MASSES.md](code/particles/EXACT_NONHADRON_MASSES.md).

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

- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: synthesis paper for the whole OPH stack; it organizes the suite on one surface and inherits claim tier from the compact recovered-core paper and the corresponding companion-paper ledgers. It does not upgrade them.
- **Paper 2. [Recovering Relativity and the Standard Model from the OPH Package Rooted in Observer Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: authoritative recovered-core and claim-tier surface for the Lorentz/gravity chain and the realized Standard Model structural branch.
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
- **[`book/`](book):** OPH Book source. Print-PDF build notes live in [`book/README.md`](book/README.md).
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
