# Observer Patch Holography: Simulation Theory, the Theory of Everything, and Observer-Centric Physics

> OPH is a mathematical reconstruction program that starts from a simple requirement: overlapping local observer descriptions on a holographic screen must agree where they overlap. From that starting point, OPH develops routes toward general relativity, toward reconstruction of the Standard Model gauge structure, toward a first-principles particle-mass program, and toward string/worldsheet descriptions emerging from the same framework. Early IBM Quantum Cloud benchmarks also provide preliminary real-device evidence for the corresponding local reduced-sector signatures. OPH also develops an interpretive strange-loop lane supported by the theorem-usable state-and-law habitat discussed in the paper set, while the stronger closure-map and stability questions remain open.

> **Status disclaimer:** OPH is an active research program and not yet fully proven. Several derivations remain incomplete, some proofs currently exist only as sketches, and certain auxiliary assumptions still need to be removed. The framework should therefore be regarded as under active development.

**French version:** [README_FR.md](README_FR.md) -- translated by Martin S.

> **Official OPH Website:** Start at [floatingpragma.io/oph/](https://floatingpragma.io/oph/).

> **OPH Disproval Challenge:** A USD 10,000 challenge to disprove OPH is currently running at [challenge.floatingpragma.io](https://challenge.floatingpragma.io).

OPH should be read directly as a mathematical and practical implementation of simulation theory. It does not just claim that reality is "like" a simulation. It specifies a concrete simulation architecture: quantum information lives on a holographic screen, observer patches act as local computational viewpoints, overlap consistency acts as the synchronization rule, and the laws of physics emerge as the stable rules that let the whole system run coherently.

For the same reason, OPH should also be read as a broad reconstruction program that aims to connect spacetime, gauge structure, particle physics, and the interpretive closure question within one observer-consistency architecture.

On that basis, OPH offers one route toward explaining both why physics has the form it does and how the existence question might be framed internally to the program. In the paper set collected here, that stronger closure story remains an interpretive strange-loop hypothesis and is not part of the recovered core.

The simulator-hardware side of OPH is developed in [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf), which turns the simulation-theory picture into explicit local screen models, overlap observables, record layers, observer criteria, and synchronization moves.

## Core Idea

OPH takes the strong observer-first position that objective reality is not fundamental but emergent from a network of subjective perspectives that must agree where they overlap.

The laws of physics are the consistency rules that make this intersubjective agreement possible.

From this starting point (plus entropy and Markov constraints), OPH treats spacetime, gauge structure, and particle physics as emergent consequences of consistency.

## Current Successes

- **Recovered core / conditional relativity branch:** OPH derives schedule-independent overlap-consistency structure, a conditional Lorentz branch, and a conditional gravity route under the stated scaling-limit premises.
- **Recovered core / gauge structure:** OPH reconstructs the Standard Model gauge structure on the stated branch and develops supplement-backed quantitative branches downstream from that structure.
- **Fixed-cutoff microphysics and observer interfaces:** OPH has an explicit constructive screen-microphysics note with local screen models, overlap observables, record layers, and observer-side fixed-cutoff theorem packages.
- **Interpretive strange-loop lane:** OPH develops an interpretive closure lane supported by the theorem-usable state-and-law habitat, while the stronger closure-map and stability questions remain open.
- **String/worldsheet continuation:** OPH connects the edge-sector heat-kernel branch to a continuation-level worldsheet reorganization rather than treating string structure as a separate starting axiom set.
- **Reduced-sector hardware fingerprints:** Early IBM Quantum Cloud benchmarks reproduce the predicted reduced-sector recoverability ordering and exact-ratio patterns on real hardware.

As a program target, OPH aims to recover the familiar low-energy effective-action form from a compact axiom set, explicit bridge premises, and minimal external inputs:

$$
\mathcal L_{\mathrm{eff}}^{\mathrm{OPH}}
\approx
\sqrt{-g}\left[
\frac{1}{16\pi G}(R-2\Lambda)
+
\mathcal L_{\mathrm{SM}}^{\mathrm{realized\ branch}}
\right]
+
\sum_i \frac{c_i}{M_*^{\Delta_i-4}}\mathcal O_i.
$$

Informally, this equation says what a completed low-energy OPH world should look like once the reconstruction is fully spelled out. The first term is gravity, with spacetime curvature and the cosmological-constant contribution. The second term is the realized Standard Model sector: gauge fields, matter fields, and their interactions. The final sum collects higher-energy corrections that should become small at ordinary scales. In other words, the displayed action is the compact "all known low-energy physics plus controlled corrections" target that the OPH derivation DAG is trying to recover from the observer-overlap framework.

## Program Goals

At full completion, OPH aims to:

- **Fully recover general relativity:** derive the gravitational sector all the way down to the familiar low-energy Einstein description, including the cosmological branch.
- **Fully recover the Standard Model:** derive the realized Standard Model structure, interactions, and low-energy effective sector from the OPH framework itself.
- **Fully recover the precise particle zoo:** derive the observed spectrum of quarks, leptons, neutrinos, hadrons, and related masses and mixings from first principles rather than treating them as inputs.
- **Fully recover string/worldsheet structure:** show how string-theoretic and worldsheet descriptions emerge as a continuation of the same observer-overlap substrate rather than as a separate axiom system.

## Papers

**Observers are all you need** is the primary technical paper. It gives the broadest technical statement of the OPH program in this repo, its main derivation branches, and the IBM Quantum Cloud benchmark surface. Its interpretive epilogue discusses the strange-loop hypothesis and the theorem-usable state-and-law habitat that supports that lane, while keeping the stronger closure-map and stability questions explicit.

- **PDF (main paper):** [Observers are all you need](paper/observers_are_all_you_need.pdf)
- **LaTeX source:** [observers_are_all_you_need.tex](paper/observers_are_all_you_need.tex)

**Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency** is the compact submission paper. It is the authoritative recovered-core and claim-tier surface for the OPH paper set in this repo. It concentrates the falsifiable core: relativity from observer-overlap consistency, the conditional gravity branch, and the gauge reconstruction program.

- **PDF (compact submission paper):** [Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)
- **LaTeX source:** [recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex)

**Reality as a Consensus Protocol** is the CS companion paper. It presents the distributed-computation view of OPH and reformulates objective law as the fixed point of reconciliation across observer patches.

- **PDF:** [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)
- **LaTeX source:** [reality_as_consensus_protocol.tex](paper/reality_as_consensus_protocol.tex)

**Screen Microphysics and Observer Synchronization** is the constructive microphysics note. It develops the simulator-hardware layer of OPH: finite local screen models, overlap observables, observer criteria, record layers, and synchronization moves.

- **PDF:** [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)
- **LaTeX source:** [screen_microphysics_and_observer_synchronization.tex](paper/screen_microphysics_and_observer_synchronization.tex)

**Toward a Particle-Spectrum Derivation from Observer-Overlap Consistency** is the particle-program note in this repo. It records the public dependency boundary for the OPH particle lane, keeps D10/D11 separate from the open charged, neutrino, and hadron continuations, and points at the compact mirrored `code/particles/` surface.

- **PDF:** [Toward a Particle-Spectrum Derivation from Observer-Overlap Consistency](paper/toward_a_particle_spectrum_derivation_from_observer_overlap_consistency.pdf)
- **LaTeX source:** [toward_a_particle_spectrum_derivation_from_observer_overlap_consistency.tex](paper/toward_a_particle_spectrum_derivation_from_observer_overlap_consistency.tex)

Release-tracked PDFs share a visible paper release line sourced from
[`paper/release_info.tex`](paper/release_info.tex). Supplemental notes can be built with the
same TeX helper, but the release manifest continues to track the core release bundle only.

## Resources

Useful entry points for reading and exploring OPH:

- **Official OPH website:** [floatingpragma.io/oph](https://floatingpragma.io/oph/)
- **Simulation theory explainer:** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Theory of everything explainer:** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything/)
- **USD 10,000 OPH disproval challenge:** [challenge.floatingpragma.io](https://challenge.floatingpragma.io)
- **OPH Book (web edition):** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Interactive OPH Lab:** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **NotebookLM notebook:** [Introduction video and guided Q&A](https://notebooklm.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a?artifactId=fb7c0ebd-4375-4997-9cae-6558ff8977b4)
- **Third-party chapter-by-chapter video course:** [Sriharsha Karamchati's OPH playlist on YouTube](https://www.youtube.com/playlist?list=PLff0tYtg64Egc2sTtKgThcPRNRdR6i83O)
- **Practical outlook:** [Potential practical applications of OPH](extra/PRACTICAL_APPLICATIONS.md)
- **OPH Sage on Telegram:** [t.me/HoloObserverBot](https://t.me/HoloObserverBot)
- **OPH Sage on X:** [x.com/OphSage](https://x.com/OphSage)
- **OPH Sage on Bluesky:** [ophsage.bsky.social](https://bsky.app/profile/ophsage.bsky.social)

## IBM Quantum Experiments

A first public IBM Quantum Cloud benchmark bundle is included in this repo. It summarizes the initial reduced-sector recoverability and exact-ratio tests, the measured hardware outputs, and the public code/data bundle used for the runs. These are benchmark checks of local OPH-motivated structure, not direct tests of the full recovered-core theorem package.

- **Experimental note:** [IBM Quantum Cloud Evidence for OPH](extra/IBM_QUANTUM_CLOUD.md)
- **Public code and data:** [code/ibm_quantum_cloud/](code/ibm_quantum_cloud/)

## Common Objections

These are rebuttals to common objections to OPH.

- [Deriving `P` from gauge data and then using `P` downstream is completely circular](extra/COMMON_OBJECTIONS.md#objection-1-circularity)
- [A fixed cell size breaks Lorentz invariance, so OPH can only recover a Newtonian limit](extra/COMMON_OBJECTIONS.md#objection-2-lorentz)
- [OPH has a Type I / Type III discontinuity, so its modular-time story is internally inconsistent](extra/COMMON_OBJECTIONS.md#objection-3-type-i-type-iii)

## Observer Patch Holography

We model reality as a network of subjective perspectives that must agree where they overlap. Concretely, we start with observer patches on a 2D holographic screen. Each patch represents a perspective with its own local data. Where patches overlap, their descriptions must agree. On the OPH interpretation, "objective reality" is the overlap-consistent backbone shared across those perspectives rather than an assumed primitive.

In OPH, Lorentz structure, gauge structure, and conservation laws are treated as consistency requirements among overlapping descriptions rather than as primitives added by hand. On the current claim ladder, the Lorentz/gravity and gauge branches still depend on the stated scaling-limit and reconstruction premises, but the program's organizing idea is that physical law is recovered from agreement constraints rather than imposed externally.

The model rests on four core concepts:

- **Screen**: A horizon-like 2D sphere (like a cosmic horizon surrounding each observer) that carries quantum information. This is where the fundamental data lives.

- **Patch**: A connected region of the screen accessible to a particular observer. Each patch has its own algebra of observables, the questions that observer can ask about reality.

- **Overlap consistency**: Where two patches share a region, their descriptions must agree. This is the central axiom. It makes objectivity a reconstructed intersubjective structure rather than a starting assumption.

- **Observer**: A stable pattern within the screen data that maintains records and participates in consistency relations.

### Reality from Computation

In the current constructive reference architecture, the screen is modeled as a gauge-invariant quantum system on a 2-sphere. It resembles a quantum cellular automaton, but with important structure: finite-dimensional local systems sit on a cellulated screen and are coupled by gauge constraints. Not all configurations are physical; only those satisfying the gauge constraints survive.

Within that reference architecture, **observer patches** are subsystems defined by boundary-gauge-invariant algebras. Each patch is a computational thread: a connected region where an observer can ask questions and get answers. The algebra $\mathcal{A}(R)$ defines what that observer can measure, namely the operators invariant under boundary gauge transformations.

**Overlap consistency** is automatic. Where two patches intersect, they access the same gauge-invariant observables. Both observers are reading the same underlying data, just from different angles. The gauge redundancy at boundaries is what makes gluing non-trivial and gives rise to the "edge modes" that carry geometric information.

**Observers are not external users.** They are emergent computational structures *within* the screen data. They are stable patterns that process information, maintain records, and create correlations. Think of them as programs running on the substrate.

**The 4D bulk is not on the sphere.** It emerges from the entanglement structure between patches. The sphere is the boundary; the interior is reconstructed holographically. When you look around and see three-dimensional space, you are experiencing a compressed encoding of how your patch is entangled with others.

*The screen is the computation. Observer patches are the threads. Reality is what they agree on.*

<a href="assets/screen.svg"><img src="assets/screen.svg" alt="The Holographic Screen" width="800"></a>

### What Drives the Computation?

In local screen models, the dynamics come from gauge-constrained quantum degrees of freedom on the screen. Observer time is not assumed to be the microscopic Hamiltonian time; it is tied to modular flow on the observer patch. On this picture, geometry and time are reconstructed from consistency and entanglement structure rather than imposed as primitives.

### Why This Approach Works

The guiding move of OPH is to treat locality, Lorentz structure, gauge structure, and gravity as consistency requirements among overlapping descriptions, not as separate ingredients added by hand. That is what makes the framework mathematically distinctive and why the papers focus so heavily on overlap algebras, entropy constraints, recoverability, and modular structure.

## The Axioms

The entire framework rests on five core axioms:

| Label | Name | Content |
|-------|------|---------|
| **A1** | Screen net | Observable algebras live on a closed 2D surface $S^2$ |
| **A2** | Overlap consistency | Where patches share a region, their descriptions agree |
| **A3** | Local MaxEnt and Refinement Stability | The realized branch is selected by a finite local MaxEnt family that persists under refinement |
| **A4** | Recoverable Generalized Entropy | $S_{\rm gen} = S_{\rm bulk} + \langle L_C \rangle$ with recoverability/collar control on the stated branch |
| **A5** | Minimal Admissible Realization | Among admissible sectors, Nature realizes the lexicographically minimal one |

Additional theorem-local technical premises and scaling-limit branch assumptions are detailed in the compact note and synchronized main paper.

## The Prediction Chain

The following infographic summarizes the full planned OPH derivation chain, from the five axioms and primary inputs through the bridge premises and major theorem branches to the intended end-state, while marking which lanes are recovered core, implementation-backed, or active completion work:

<a href="assets/prediction-chain.svg"><img src="assets/prediction-chain.svg" alt="OPH full derivation DAG poster" width="1200"></a>

*Click the poster to open the full SVG in the browser. The top band separates axioms, primary constants, and removable bridge assumptions; the lower bands trace the full planned OPH derivation DAG toward completed TOE synthesis and global closure.*

> **Particle Spectrum Derivation**: The derivation from pixel area to the particle-mass program, together with comparisons against PDG data and audit checks, is documented in **[the spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex)**.

## Particle Results

Selected particle-side successes:

| Result | OPH output | Reference |
|------|------|------|
| Photon mass | exactly `0` | experimental upper limit `< 1e-18 eV` |
| Gluon mass | exactly `0` | structural massless color sector |
| `W` boson mass | `80.3863 GeV` | `80.377 GeV` |
| Higgs mass | `126.5 GeV` | `125.20 GeV` |
| Top mass | `171.1 GeV` | `172.4 GeV` |

The fuller particle ledger, including open and incorrect rows, is collected in
[the particle note](paper/toward_a_particle_spectrum_derivation_from_observer_overlap_consistency.pdf),
[the SVG ledger](assets/particle_mass_derivation_graph.svg),
and [code/particles/RESULTS_STATUS.md](code/particles/RESULTS_STATUS.md).

## Repository Contents

This repository is organized around the OPH paper set in this repo and its supporting material.

- **[`paper/`](paper):** release-tracked PDFs, supplemental notes, LaTeX sources, and shared paper metadata. This is the canonical home of the main paper, the compact submission paper, the CS companion paper, the screen-microphysics note, and the particle-program note.
- **[`paper/tex_fragments/`](paper/tex_fragments):** shared derivation fragments used by the longer papers, including the gauge, spectrum, technical-supplement, and string-theory source files.
- **[`book/`](book):** Markdown source for the OPH Book web edition.
- **[`code/particles/`](code/particles):** compact public export of the active particle program: `core/` compatibility modules, closure-lane directories, frozen artifacts, and status ledgers.
- **[`code/ibm_quantum_cloud/`](code/ibm_quantum_cloud):** IBM Quantum Cloud experiments, data, and hardware-facing utilities.
- **[`extra/`](extra):** supporting notes such as common objections, the IBM Quantum writeup, and practical-application notes.
- **[`assets/`](assets):** figures and diagrams used across the papers, README, and public materials.

## Code

The repository code is organized by paper-facing lanes, ledgers, and artifacts. The main entry points are:

| Path | Purpose |
|------|---------|
| [code/particles/README.md](code/particles/README.md) | Map of the compact public particle mirror and its claim boundary |
| [code/particles/RESULTS_STATUS.md](code/particles/RESULTS_STATUS.md) | Human-readable public particle-status ledger |
| [code/particles/results_status.json](code/particles/results_status.json) | Machine-readable mirror of the particle-status ledger |
| [code/particles/core/particle_masses_paper_d10_d11.py](code/particles/core/particle_masses_paper_d10_d11.py) | Paper-synchronized D10/D11 reconstruction and transport layer |
| [code/particles/core/oph_predict_compare.py](code/particles/core/oph_predict_compare.py) | Predictor backbone and comparison surface used by the compact export |
| [code/particles/calibration/](code/particles/calibration) | D10/D11 calibration and exactness-audit lane |
| [code/particles/flavor/](code/particles/flavor) | Shared charged-sector transport, projector, and quark continuation lane |
| [code/particles/leptons/](code/particles/leptons) | Charged-lepton continuation lane |
| [code/particles/neutrino/](code/particles/neutrino) | Neutrino continuation lane with forward artifacts and guards |
| [code/particles/hadron/](code/particles/hadron) | Hadron debug/systematics lane and thin lattice wrapper |

## Book
 
Core explanatory chapters

| Chapter | Title | Topic |
|---------|-------|-------|
| [Prologue](book/prologue.md) | Prologue | Setting the stage |
| [1](book/chapter-01-consistency.md) | Consistency | Observer agreement as the fundamental principle |
| [2](book/chapter-02-lineage.md) | Lineage | Historical roots of holographic ideas |
| [3](book/chapter-03-screen.md) | The Screen | Holographic screens and information bounds |
| [4](book/chapter-04-entropy.md) | Entropy | Thermodynamics and the arrow of time |
| [5](book/chapter-05-algebra.md) | Algebra | The mathematical structure of observables |
| [6](book/chapter-06-overlap.md) | Overlap | Consistency conditions and Bell's theorem |
| [7](book/chapter-07-recovery.md) | Recovery | Information preservation and quantum error correction |
| [8](book/chapter-08-holography.md) | Holography | AdS/CFT and bulk reconstruction |
| [9](book/chapter-09-entanglement.md) | Entanglement | Geometry from quantum correlations |
| [10](book/chapter-10-error-correction.md) | Error Correction | Reality as a quantum code |
| [11](book/chapter-11-maxent.md) | MaxEnt | Entropy, time, and modular flow |
| [12](book/chapter-12-symmetry.md) | Symmetry | Conservation laws from consistency |

Physics branch chapters

| Chapter | Title | Topic |
|---------|-------|-------|
| [13](book/chapter-13-desitter.md) | De Sitter | Our universe's holographic screen |
| [14](book/chapter-14-standard-model.md) | Standard Model | Particles from gluing constraints |
| [15](book/chapter-15-relativity.md) | Relativity | Spacetime from modular time |
| [16](book/chapter-16-matter.md) | Matter | Classical physics as emergent stability |
| [17](book/chapter-17-darwin.md) | Darwin's Laws | Laws as evolutionary survivors |

Interpretive chapters

| Chapter | Title | Topic |
|---------|-------|-------|
| [18](book/chapter-18-synthesis.md) | Synthesis | Putting it all together |
| [19](book/chapter-19-metaphysics.md) | Metaphysics | Qualia and the hard problem |

Speculative epilogue

| Chapter | Title | Topic |
|---------|-------|-------|
| [Epilogue](book/epilogue.md) | Epilogue | One last surprise |

The book source lives in [`book/`](book). It is the Markdown source for the public OPH Book and is organized as a 21-chapter sequence from the prologue and observer-consistency foundations through holography, symmetry, the Standard Model, relativity, synthesis, and metaphysics.

## Release Workflow

Paper releases are managed from the shared files under [`paper/`](paper):

1. bump [`paper/release_info.tex`](paper/release_info.tex)
2. rebuild the release-tracked PDFs with `python3 tools/build_tex_papers.py --release-only`
3. regenerate [`paper/paper_release_manifest.json`](paper/paper_release_manifest.json)

To build every hand-authored TeX paper in the repo, including supplemental notes such as the
particle-program paper, run `python3 tools/build_tex_papers.py`.

The shared release line applies across the current release-tracked paper set. The manifest records the synchronized release state for the tracked paper bundle.

## Contributing

For corrections, suggestions, or additions, please open a pull request.

## License

This repository is licensed under the
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

You are free to share and adapt the material for non-commercial purposes,
provided proper attribution is given and derivative works are licensed
under identical terms.

For commercial licensing inquiries, please contact:
bernhard@floatingpragma.io
