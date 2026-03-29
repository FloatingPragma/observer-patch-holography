# Observer Patch Holography: Simulation Theory, the Theory of Everything, and Observer-Centric Physics

> OPH is a mathematical reconstruction program that starts from a simple requirement: overlapping local observer descriptions on a finite-capacity holographic screen must agree where they overlap. From that starting point, the current paper suite develops an explicit fixed-cutoff collar and higher-gauge package, a conditional route to Lorentzian geometry and Einstein dynamics, a conditional route to compact gauge reconstruction and the realized Standard Model quotient, a current particle-spectrum branch, and a fixed-cutoff measurement/observer package. The stronger strange-loop closure story remains interpretive rather than part of the recovered core.

> **Status disclaimer:** OPH is an active research program and not yet fully proven. Several derivations remain incomplete, some proofs currently exist only as sketches, and certain auxiliary assumptions still need to be removed. The framework should therefore be regarded as under active development.

**French version:** [README_FR.md](README_FR.md) -- translated by Martin S.

> **Official OPH Website:** Start at [floatingpragma.io/oph/](https://floatingpragma.io/oph/).

> **OPH Disproval Challenge:** A USD 10,000 challenge to disprove OPH is currently running at [challenge.floatingpragma.io](https://challenge.floatingpragma.io).

OPH can be read as a route toward a mathematical and practical simulation-theory interpretation. It does not stop at saying that reality is "like" a simulation: it specifies a concrete information-processing architecture in which quantum data live on a holographic screen, observer patches act as local computational viewpoints, overlap consistency acts as the synchronization rule, and the laws of physics emerge as the stable rules that let the whole system run coherently.

For the same reason, OPH should also be read as a broad reconstruction program that aims to connect spacetime, gauge structure, particle physics, and the closure question within one observer-consistency architecture, while keeping the stronger closure story separate from the recovered core.

On that basis, OPH offers one route toward explaining both why physics has the form it does and how the existence question might be framed internally to the program. In the paper set collected here, the stronger global-closure story remains an interpretive strange-loop hypothesis and is not part of the recovered core.

The simulator-hardware side of OPH is developed in [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf), which turns the simulation-theory picture into explicit local screen models, overlap observables, record layers, observer criteria, and synchronization moves.

## Core Idea

OPH takes the strong observer-first position that objective reality is not fundamental but emergent from a network of subjective perspectives that must agree where they overlap.

The laws of physics are the consistency rules that make this intersubjective agreement possible.

From this starting point (plus entropy and Markov constraints), OPH treats spacetime, gauge structure, and particle physics as emergent consequences of consistency.

## Current Achievements

- **OPH now has a real synthesis paper.** *Observers Are All You Need* presents the current program on one common theorem and status surface, so readers can see in one place what OPH is, what is already sharp, and what is still conditional.
- **The fixed-cutoff core is explicit.** The collar package is constructive and interacting at fixed cutoff, and the topological branch is closed on the ordinary, central-defect, and genuinely noncentral higher-gauge branches.
- **There is a disciplined relativity/gravity route.** OPH derives schedule-independent overlap consistency and develops a conditional Lorentz/null-modular/Einstein branch under an explicit scaling-limit premise ledger.
- **There is a disciplined gauge route.** OPH develops a conditional compact-gauge reconstruction route and, on the realized branch, recovers the Standard Model quotient \(\SU(3)\times\SU(2)\times\U(1)/\mathbb Z_6\) together with the realized hypercharge/counting chain.
- **The particle branch has moved beyond pure structure.** Structural massless carriers are fixed; the electroweak and Higgs/top stages emit quantitative rows; the wider quark, lepton, and neutrino sectors remain openly status-split rather than silently promoted, while hadrons are explicitly execution-bound rather than claimed as paper-derived outputs.
- **Measurement and observers are part of the theorem surface.** OPH now has fixed-cutoff record, Born-Luders, checkpoint/restoration, and observer-interface packages rather than treating observers as pure interpretation.
- **The microphysics side is concrete enough to simulate.** The screen-microphysics paper gives a finite reference architecture with explicit registers, overlap observables, records, repair moves, and synchronization tests.
- **The worldsheet/string lane is integrated but not oversold.** OPH relates the edge-sector heat-kernel theorem to a continuation-level worldsheet reorganization while keeping that branch outside the recovered core.
- **Some hardware-facing signatures exist already.** Early IBM Quantum Cloud benchmarks reproduce the expected reduced-sector recoverability ordering and exact-ratio patterns on real devices.

## Particle Status

The particle derivation is hard and still a work in progress. The current particle paper should be
read as a derivation document with real outputs in hand, not as a claim that the full observed
particle spectrum is finished.

- **Exact structural outputs already on the current surface:** photon, gluons, graviton.
- **Current strongest quantitative rows:** \(W\) boson and \(Z\) boson on the electroweak calibration stage, plus the Higgs boson and top quark on the Higgs/top critical stage.
- **Still unfinished:** the broader quark family remains continuation-level; charged leptons and flavor-labeled neutrinos remain open theorem/closure fronts; and hadrons are execution-bound because promotable rows require real nonperturbative production computation and systematics rather than theorem-only completion.

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
- **Fully recover the Standard Model particle zoo:** derive the full measured Standard Model particle spectrum, not just the gauge structure and partial quantitative descendants.
- **Fully recover string/worldsheet structure:** show how string-theoretic and worldsheet descriptions emerge as a continuation of the same observer-overlap substrate rather than as a separate axiom system.

## Papers

**Observers Are All You Need** is the best place to start if you want the big picture. It is the synthesis paper: what OPH is trying to do, what the current strongest achievements are, which parts are theorem-level, and where the real open problems still sit.

- **PDF (main paper):** [Observers are all you need](paper/observers_are_all_you_need.pdf)
- **LaTeX source:** [observers_are_all_you_need.tex](paper/observers_are_all_you_need.tex)

**Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency** is the compact recovered-core paper. If you want the sharpest and most falsifiable OPH claim surface, start here: it isolates the relativity branch, the gravity branch, the gauge-reconstruction route, and the fixed-cutoff higher-gauge closure package.

- **PDF (compact submission paper):** [Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)
- **LaTeX source:** [recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex)

**Deriving the Particle Zoo from Observer Consistency** is the particle-spectrum status paper. It asks how far one shared calibration scale can be pushed through the observed spectrum, records the current bosonic/Higgs-top/quark outputs, keeps the charged-lepton and neutrino theorem frontiers explicit, and marks hadrons as execution-bound rather than paper-derived output rows.

- **PDF:** [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)
- **LaTeX source:** [deriving_the_particle_zoo_from_observer_consistency.tex](paper/deriving_the_particle_zoo_from_observer_consistency.tex)

**Reality as a Consensus Protocol** is the computational spine paper. It reformulates OPH as asynchronous local reconciliation on a patch network and makes objective law precise as a schedule-independent fixed point on the gauge quotient.

- **PDF:** [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)
- **LaTeX source:** [reality_as_consensus_protocol.tex](paper/reality_as_consensus_protocol.tex)

**Screen Microphysics and Observer Synchronization** is the engineering-side paper. It gives a concrete finite screen architecture with local registers, overlap observables, record layers, observer criteria, repair interfaces, and synchronization targets.

- **PDF:** [Screen Microphysics and Observer Synchronization](paper/screen_microphysics_and_observer_synchronization.pdf)
- **LaTeX source:** [screen_microphysics_and_observer_synchronization.tex](paper/screen_microphysics_and_observer_synchronization.tex)

Release-tracked PDFs share a visible paper release line sourced from
[`paper/release_info.tex`](paper/release_info.tex). Supplemental notes can be built with the
same TeX helper, but the release manifest continues to track the current numbered paper bundle only.

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

Additional theorem-local technical premises and scaling-limit branch assumptions are detailed in *Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency* and synchronized into the main synthesis paper.

## The Prediction Chain

The following infographic summarizes the full planned OPH derivation chain, from the five axioms and primary inputs through the bridge premises and major theorem branches to the intended end-state, while marking which lanes are recovered core, implementation-backed, or active completion work:

<a href="assets/prediction-chain.svg"><img src="assets/prediction-chain.svg" alt="OPH full derivation DAG poster" width="1200"></a>

*Click the poster to open the full SVG in the browser. The top band separates axioms, primary constants, and removable bridge assumptions; the lower bands trace the full planned OPH derivation DAG toward completed TOE synthesis and global closure.*

## Repository Contents

This repository is organized around the OPH paper set in this repo and its supporting material.

- **[`paper/`](paper):** release-tracked PDFs, supplemental notes, LaTeX sources, and shared paper metadata. This is the canonical home of the main paper, the compact submission paper, the particle-spectrum paper, the CS companion paper, and the screen-microphysics note.
- **[`paper/tex_fragments/`](paper/tex_fragments):** shared derivation fragments used by the longer papers, including the gauge, technical-supplement, and string-theory source files.
- **[`book/`](book):** Markdown source for the OPH Book web edition.
- **[`code/ibm_quantum_cloud/`](code/ibm_quantum_cloud):** IBM Quantum Cloud experiments, data, and hardware-facing utilities.
- **[`extra/`](extra):** supporting notes such as common objections, the IBM Quantum writeup, and practical-application notes.
- **[`assets/`](assets):** figures and diagrams used across the papers, README, and public materials.

## Code

The repository code is organized by paper-facing lanes, ledgers, and artifacts. The main public entry point is:

| Path | Purpose |
|------|---------|
| [code/ibm_quantum_cloud/](code/ibm_quantum_cloud/) | IBM Quantum Cloud experiments, hardware-facing utilities, and public benchmark data |

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
