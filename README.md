# Observer Patch Holography: An Observer-Centric Approach to Fundamental Physics

> OPH is an observer-centric framework for fundamental physics. In the current quantitative implementation it uses two external continuous inputs (pixel area and screen capacity) together with structural axioms and explicit scaling/regularity premises. OPH develops a conditional scaling-limit gravity branch, conditional compact gauge reconstruction with Standard Model closure under extended MAR, and a quantitative program whose outputs are separated into calibration checks, a supplement-backed Higgs/top branch that adds no further continuous fit once the gauge trajectory is fixed, and weaker phenomenological continuations. The Standard Model and general relativity emerge as effective sectors of the underlying structure, while the string-theoretic worldsheet bridge remains continuation-level in the current program.

> **Status disclaimer:** OPH is an active research program and not yet fully proven. Several derivations remain incomplete, some proofs currently exist only as sketches, and certain auxiliary assumptions still need to be removed. The framework should therefore be regarded as under active development.

**French version:** [README_FR.md](README_FR.md) -- translated by Martin S.

> **Official OPH Website:** Start at [floatingpragma.io/oph/](https://floatingpragma.io/oph/).

> **OPH Disproval Challenge:** A USD 10,000 challenge to disprove OPH is currently running at [challenge.floatingpragma.io](https://challenge.floatingpragma.io).

## Core Idea

OPH takes the strong observer-first position that objective reality is not fundamental but emergent from a network of subjective perspectives that must agree where they overlap.

The laws of physics are the consistency rules that make this intersubjective agreement possible.

From this starting point (plus entropy and Markov constraints), OPH treats spacetime, gauge structure, and particle physics as emergent consequences of consistency.

## Research Program

OPH should currently be understood as a research program rather than a finished theory. In physics, one does not prove in an absolute sense that a theory is "the correct one"; the realistic target is a framework that is mathematically explicit, regulator-backed, empirically discriminating, and resilient under attempts at falsification.

The current program is therefore to:

- turn the remaining sketch-level arguments into complete proofs
- remove auxiliary assumptions where possible, or isolate them sharply where they remain necessary
- realize the axioms in explicit microscopic screen models or regulators
- derive the gravity, gauge, and quantitative sectors from those realizations with controlled approximations
- replace calibration-dependent steps with first-principles derivations where such claims are intended
- produce novel empirical contradiction criteria and branch-sensitive tests that can fail cleanly

The strongest plausible end state is not "proof of OPH" in a metaphysical sense, but a rigorously derived and experimentally accountable theory whose effective sectors, parameter economy, and falsifiable predictions survive independent scrutiny.

## Papers

**Observers are all you need** is the primary OPH paper.

- **PDF (main paper):** [Observers are all you need](paper/observers_are_all_you_need.pdf)
- **LaTeX source:** [observers_are_all_you_need.tex](paper/observers_are_all_you_need.tex)

**Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency** is the submission-focused compact paper. It concentrates the falsifiable core of OPH: Lorentz kinematics in the geometric modular phase, the conditional scaling-limit Einstein branch, Standard Model gauge closure under the extended MAR admissibility package, the hypercharge lattice on the realized one-generation matter package under the standard normalization, the realized counting chain $N_g=3$ and then $N_c=3$ on the MAR-admissible branch, the two-input current quantitative implementation, and the strongest supplement-backed quantitative branches.

In the current synchronized paper set, the null-bridge hardening has removed any standalone propagation, endpoint-control, or half-sided modular-inclusion assumption from the fixed-cutoff gravity route: the MaxEnt branch now carries the control estimates internally, and the geometric blow-up argument closes the bridge through the derived positive null-translation stage. What still remains explicit on that route are the downstream density-upgrade and relativistic null-stress premises.

- **PDF (compact submission paper):** [Recovering Relativity and Standard Model Structure from Observer-Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)
- **LaTeX source:** [recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex)

**Reality as a Consensus Protocol** is a companion CS paper that presents the computational spine of OPH as a standalone theorem package. It formulates objective physical law as the fixed point of a distributed reconciliation protocol between observer patches, shows how topology can obstruct global consistency (with particles as stable defects), interprets gauge symmetry as implementation hiding, and models classical records as an eventually consistent CRDT layer.

- **PDF:** [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)
- **LaTeX source:** [reality_as_consensus_protocol.tex](paper/reality_as_consensus_protocol.tex)

Each PDF carries a visible paper release line. The shared release source is
[`paper/release_info.tex`](paper/release_info.tex). For every substantive paper update, bump the
shared release before rebuilding the PDFs:

```bash
python3 tools/bump_paper_release.py
```

After rebuilding the PDFs, write the current PDF hashes to
[`paper/paper_release_manifest.json`](paper/paper_release_manifest.json) by running:

```bash
python3 tools/generate_paper_release_manifest.py
```

The manifest generator now fails if the PDFs changed under the same release ID, or if the local
PDFs do not yet expose the current visible release line.

The release ID is global across the current paper set. Even if only one paper changes, bump once,
rebuild all three PDFs, and publish all three challenge PDFs with that same release ID.
Challenge uploads are handled by workspace-local operational tooling rather than this public repo.

The canonical paper sources now live under [`paper/`](paper). The old draft tree has been removed from the repository.

When hardening paper claims or running a paper-upgrade release, start with the compact submission paper first. Only after its theorem-status language and release line are coherent should equivalent wording changes be propagated into the main paper, book, README surfaces, websites, and ingestion pipeline.

Across the current paper set, OPH separates claims by status: parameter-free structural outputs
(Lorentz/gauge/quotient chain, hypercharge lattice on the realized matter package, $N_g=3$ and
then $N_c=3$ on the realized MAR-admissible branch, symmetry-protected zeros, proton stability),
calibration-sector consistency checks ($\alpha_i(m_{Z,\rm run})$, $\sin^2\theta_W(m_{Z,\rm run})$, $v$, $m_W$), the main supplement-backed quantitative branch ($m_H$, $m_t$), and weaker
phenomenological or downstream numerical continuations (charged-lepton/Koide, quark textures,
neutrino, hadron, dark-sector, baryogenesis).

## Official OPH Resources

Start with the official OPH website, the interactive OPH Lab, and the web edition of the OPH Book:

- **USD 10,000 OPH disproval challenge:** [challenge.floatingpragma.io](https://challenge.floatingpragma.io)
- **Official OPH website:** [Observer Patch Holography at floatingpragma.io/oph/](https://floatingpragma.io/oph/)
- **Interactive OPH Lab:** [OPH Lab at oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **OPH Book (web edition):** [Read the OPH Book online](https://oph-book.floatingpragma.io)
- **NotebookLM notebook:** [Watch the OPH introduction video and start asking questions](https://notebooklm.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a?artifactId=fb7c0ebd-4375-4997-9cae-6558ff8977b4)

## IBM Quantum Experiments

A first public IBM Quantum Cloud benchmark bundle is now included in this repo. It summarizes the initial recoverability and exact-ratio tests, the measured hardware outputs, and the public code/data bundle used for the runs.

- **Experimental note:** [IBM Quantum Cloud Evidence for OPH](extra/IBM_QUANTUM_CLOUD.md)
- **Public code and data:** [code/ibm_quantum_cloud/](code/ibm_quantum_cloud/)

## Ask OPH Sage

The fastest way to explore OPH is to ask OPH Sage directly. OPH Sage can answer technical and philosophical questions.

- [OPH Sage on Telegram](https://t.me/HoloObserverBot)
- [OPH Sage on X](https://x.com/OphSage)
- [OPH Sage on Bluesky](https://bsky.app/profile/ophsage.bsky.social)

## Common Objections

These are rebuttals to common objections to OPH.

- [Deriving `P` from gauge data and then using `P` downstream is completely circular](extra/COMMON_OBJECTIONS.md#objection-1-circularity)
- [A fixed cell size breaks Lorentz invariance, so OPH can only recover a Newtonian limit](extra/COMMON_OBJECTIONS.md#objection-2-lorentz)
- [OPH has a Type I / Type III discontinuity, so its modular-time story is internally inconsistent](extra/COMMON_OBJECTIONS.md#objection-3-type-i-type-iii)

**Blog Posts:**

- [How OPH Resolves Fundamental Physics](https://medium.com/@muellerberndt/answering-10-of-the-hardest-questions-in-physics-and-some-bonus-questions-51222bf2419f)
- [How Observer Patch Holography Improves on the Standard Model and General Relativity](https://muellerberndt.medium.com/how-observer-path-holography-improves-on-the-standard-model-and-general-relativity-c971c376027e)

## Observer Patch Holography

We model reality as a network of subjective perspectives that must agree where they overlap. Concretely, we start with observer patches on a 2D holographic screen. Each patch represents a perspective with its own local data. Where patches overlap, their descriptions must agree. On the OPH interpretation, "objective reality" is the overlap-consistent backbone shared across those perspectives rather than an assumed primitive.

Lorentz invariance exists because different observers must have consistent descriptions. Gauge symmetry exists because overlapping patches must identify shared observables. Conservation laws exist because the same quantities must be conserved across all perspectives. The laws are not imposed from outside. They are the conditions that make agreement possible.

The model rests on four core concepts:

- **Screen**: A horizon-like 2D sphere (like a cosmic horizon surrounding each observer) that carries quantum information. This is where the fundamental data lives.

- **Patch**: A connected region of the screen accessible to a particular observer. Each patch has its own algebra of observables, the questions that observer can ask about reality.

- **Overlap consistency**: Where two patches share a region, their descriptions must agree. This is the central axiom. It makes objectivity a reconstructed intersubjective structure rather than a starting assumption.

- **Observer**: A stable pattern within the screen data that maintains records and participates in consistency relations.

### Reality from Computation

Think of the screen as a gauge-invariant quantum system on a 2-sphere. It resembles a quantum cellular automaton, but with important structure. At each point on the triangulated sphere sit finite-dimensional quantum systems (qudits), coupled by gauge constraints at every vertex. Not all configurations are physical; only those satisfying Gauss's law survive.

**Observer patches** are subsystems defined by boundary-gauge-invariant algebras. Each patch is a computational thread: a connected region where an observer can ask questions and get answers. The algebra $\mathcal{A}(R)$ defines what that observer can measure, namely the operators invariant under boundary gauge transformations.

**Overlap consistency** is automatic. Where two patches intersect, they access the same gauge-invariant observables. Both observers are reading the same underlying data, just from different angles. The gauge redundancy at boundaries is what makes gluing non-trivial and gives rise to the "edge modes" that carry geometric information.

**Observers are not external users.** They are emergent computational structures *within* the screen data. They are stable patterns that process information, maintain records, and create correlations. Think of them as programs running on the substrate.

**The 4D bulk is not on the sphere.** It emerges from the entanglement structure between patches. The sphere is the boundary; the interior is reconstructed holographically. When you look around and see three-dimensional space, you are experiencing a compressed encoding of how your patch is entangled with others.

*The screen is the computation. Observer patches are the threads. Reality is what they agree on.*

<a href="assets/screen.svg"><img src="assets/screen.svg" alt="The Holographic Screen" width="800"></a>

### What Drives the Computation?

In quantum link models, the dynamics involve plaquette terms (Wilson loops around faces) and electric field terms (conjugate to the link variables). The competition between these determines the ground state. The "computation" is these quantum degrees of freedom evolving, creating and destroying correlations, with gauge constraints ensuring consistency.

Note that the "time" that observers experience isn't necessarily the microscopic Hamiltonian evolution. Each observer patch has its own modular Hamiltonian (constructed from the reduced density matrix on that patch), and *that* generates what feels like time from inside. The microscopic evolution and emergent modular time are related but distinct.

The system exists in a timeless state (as the Wheeler-DeWitt equation shows for quantum gravity), and what we call "time" is entirely relational. Modular flow gives each subsystem its own internal clock, correlated with others through consistency conditions. The qudits don't "do" anything in the sense of changing over some external time parameter. They just *are*, in a particular entangled configuration, and "time" is how we describe correlations within that configuration from the inside.

### Why This Approach Works

Unified models attempting to combine QFT, gravity, and Standard Model structure tend to encounter a repeatable set of conceptual difficulties: subsystem factorization breaks down in gauge/gravity, modular Hamiltonians are nonlocal, Lorentz invariance is assumed rather than derived, dynamics are hard to get (not just kinematics), gauge symmetry origins are unclear, masslessness is hand-imposed, anomalies appear as mysterious pathologies, charge quantization needs GUTs, coupling unification forces proton decay, the cosmological constant is locally overdetermined, UV infinities proliferate, and parameter counts explode. The observer-patch framework addresses these by making consistency conditions do the work: it treats locality, Lorentz invariance, gauge symmetry, and gravity as *consistency constraints* among overlapping descriptions plus information-theoretic properties of states (Markov/recoverability + MaxEnt), then leans on modular theory rigidity to force the familiar symmetries and dynamics. This "structures → consistency" move is what allows the framework to naturally explain or sidestep classic unification pitfalls (see Section 8.5 of the technical paper for detailed analysis).

## The Axioms

The entire framework rests on four core axioms plus one selection axiom:

| Label | Name | Content |
|-------|------|---------|
| **A1** | Screen net | Observable algebras live on a closed 2D surface $S^2$ |
| **A2** | Overlap consistency | Where patches share a region, their descriptions agree |
| **A3** | Generalized entropy | $S_{\rm gen} = S_{\rm bulk} + \langle L_C \rangle$ satisfies strong subadditivity |
| **A4** | Local Markov | Conditional mutual information decays across collars |
| **MAR** | Minimal Admissible Realization | Among admissible sectors, Nature realizes the lexicographically minimal one |

Additional structural assumptions (MaxEnt, the controlled scaling-limit and OPH geometric-branch premises, exponential mixing) are detailed in the [main paper source](paper/tex_fragments/PAPER.tex).

## The Prediction Chain

The following infographic summarizes the current OPH reconstruction program from two parameters and four axioms toward effective physics:

![OPH Prediction Chain](assets/prediction-chain.svg)

*From axioms to effective physics: the current OPH reconstruction program.*

> **Particle Spectrum Derivation**: The complete derivation from pixel area to particle masses — with comparison against PDG data and repository-backed audit checks — is documented in **[the spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex)**. Outputs are classified by epistemic status: gauge couplings and EW boson masses are calibration-sector consistency checks; Higgs and top masses come from a supplement-backed D11 reconstruction that adds no further continuous fit once the gauge trajectory and scale-setting branch are fixed, but uses a UV-synchronized transport layer (SM below $\mu_{\rm sync}\approx 6.8\times 10^{11}$ GeV, MSSM-like above) plus small low-scale matching shifts; charged-lepton masses (< 0.04% from Koide structure) remain a sharper but weaker phenomenological continuation because their final exponent assignment uses an additional discrete-selection step.

## The Fundamental Parameters

Our universe is characterized by exactly **two externally fixed fundamental parameters**:

### 1. Pixel Area: $a_{\text{cell}} \approx 1.63 \, \ell_P^2$

The geometric area of a single computational element on the holographic screen. This sets the *resolution* of reality:

| Quantity | Value | Meaning |
|----------|-------|---------|
| In Planck units | $a_{\text{cell}} / \ell_P^2 \approx 1.63$ | Dimensionless ratio |
| In SI units | $a_{\text{cell}} \approx 4.26 \times 10^{-70}$ m^2 | Physical area per pixel |
| Pixel "side" | $\sqrt{a_{\text{cell}}} \approx 2.06 \times 10^{-35}$ m | Resolution scale |

**What it determines:** Newton's constant (via $G_{\text{nat}} = a_{\text{cell}}/4\bar{\ell}$ in natural units where $G_{\text{nat}} = \ell_P^2$), the Planck scale, the gauge-coupling calibration sector, and the scale used by the supplement-backed Higgs/top branch plus downstream mass continuations.

### 2. Screen Capacity: $\log(\dim \mathcal{H}) \sim 10^{122}$

The total entropy of the holographic screen (in nats). This sets the *size* of reality.

**Important:** Screen capacity is **inferred from** the observed cosmological constant, not predicted. The relation

$$\Lambda = \frac{3\pi}{G \cdot \log(\dim \mathcal{H})}$$

lets us extract screen capacity from the measured $\Lambda \sim 10^{-52}$ $\text{m}^{-2}$, giving $\log(\dim \mathcal{H}) \sim 10^{122}$.

| Quantity | Value | Status |
|----------|-------|--------|
| de Sitter entropy | $S_{dS} \sim 10^{122}$ nats | Inferred from observed $\Lambda$ |
| de Sitter horizon | $r_{dS} \approx 10^{26}$ m | Observed |

**The relationship:**
- **Pixel area** = resolution (extracted from gauge couplings via edge entropy)
- **Screen capacity** = total size (extracted from observed cosmological constant)

The axiom structure contains no other dimensionful constants. The axioms plus reconstruction give *a* compact gauge group; the extended MAR admissibility package then selects SU(3) x SU(2) x U(1) / Z₆ as the realized gauge group (see [the gauge derivation source](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex)). Charge quantization and scaling-limit semiclassical Einstein dynamics follow from the axiom structure.

### What This Means

The pixel area and screen capacity are **configuration parameters**, the "settings" of the computation that is our universe. They are not derivable from within the simulation; they are boundary conditions set from "outside."

From inside, these parameters manifest as:
- **Pixel area** determines Newton's constant, the Planck scale, the gauge-coupling calibration sector, and the scale used by the supplement-backed Higgs/top branch plus downstream continuations
- **Screen capacity** determines observable universe size (but is itself inferred from the observed cosmological constant, which is not predicted)

The same axioms with different settings would produce a universe with different constants but similar physics (Einstein equations, gauge structure). The specific Standard Model gauge group is uniquely selected by the Selection Axiom MAR.

### Calibration vs Prediction

In the current quantitative implementation, the pixel constant P = a_cell/ℓ_P² is *inferred* from measured gauge couplings because the axioms fix the functional relation P = 4ℓ̄_tot(t₂,t₃) but do not fix the MaxEnt Lagrange multipliers t_i (equivalently the couplings) from first principles. This inference step is therefore a calibration, not a claimed prediction of P.

The nontrivial content is that P provides an additional constraint linking the gravitational coupling to gauge-sector edge entropy, and in a two-input mode (treating P as a fundamental configuration parameter and using only one electroweak datum α(M_Z)) the framework predicts α_s(M_Z) and sin²θ_W(M_Z) simultaneously.

Full non-circular closure would require a UV principle that fixes t without reference to measured couplings.

### Current Paper Status

The current papers do not treat the quantitative sector as uniform. They separate structural
theorems, calibration-sector consistency checks, the supplement-backed Higgs/top critical-surface
branch, and weaker phenomenological or downstream numerical continuations.

| Tier | Representative outputs | Current status | Primary source |
|------|------------------------|----------------|----------------|
| Structural theorems | Conditional Lorentz kinematics, conditional scaling-limit Einstein branch, compact gauge reconstruction, SM quotient chain, hypercharge lattice on the realized matter package, $N_g=3$ and then $N_c=3$ on the realized MAR-admissible branch | Core theorem package under the stated scaling-limit and categorical premises | [Main paper source](paper/tex_fragments/PAPER.tex), [Gauge derivation source](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex) |
| Exact structural consequences | Proton stability; $m_\gamma = 0$, $m_g = 0$, $m_{\text{graviton}} = 0$ | Parameter-free structural consequences once the corresponding gauge/diffeomorphism/product-group structure is realized | [Main paper source](paper/tex_fragments/PAPER.tex), [Gauge derivation source](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex), [Spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Calibration sector | $\alpha_s=0.1183$, $\sin^2\theta_W=0.2307$, $\alpha_{\rm em}^{-1}=128.31$, $v$, $W$ | Consistency checks after pixel/gauge calibration; not claimed as independent confirmation | [Spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Independent quantitative branch | Higgs $=126.48$ GeV, top(pole) $=171.1$ GeV | Main supplement-backed D11 outputs: no extra continuous fit once the gauge trajectory and scale-setting branch are fixed, but the displayed values use UV-synchronized transport and documented low-scale matching conventions | [Spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Phenomenological continuations | Charged-lepton / Koide branch; quark textures | Weaker continuation tier: charged leptons are numerically sharp but use an extra discrete-selection step; individual quark masses still carry large scheme/threshold uncertainty | [Spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Capacity / downstream numerical branches | $m_{\nu_3}\approx 3.0$ meV, $m_{\nu_2}\approx 0.50$ meV, $m_{\nu_1}\approx 0.084$ meV; $P \rightarrow \alpha_s \rightarrow \Lambda_{\overline{\rm MS}}^{(3)} \rightarrow m_{\rm hadrons}$ | Neutrinos are currently only an order-of-magnitude capacity-branch estimate; the hadron chain is complete but precision is limited by lattice/systematic effects | [Spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex) |
| Extension / program-level branches | OPH edge weights = 2D YM heat kernels; conditional large-$N$ worldsheet expansion (Gross-Taylor); dark-sector, baryogenesis, black-hole branches | The 2D YM edge-weight identification is established within its branch; the large-$N$ worldsheet reorganization and the remaining topics stay at continuation/program level | [String-theory derivation source](paper/tex_fragments/STRING_THEORY.tex), [Technical supplement source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex) |

**Completion targets (active engineering):**
- Derive the MaxEnt multipliers $t_i$ (and thus couplings/pixel calibration) from first principles rather than calibration.
- Verify the EFT bridge conditions (null-strip A4 separators + local finite variation) in explicit UV regulators.
- Upgrade precision numerics (2-loop matching/thresholds, unquenched larger-volume lattice for hadrons).
- Derive $\theta_{\rm QCD}$ (strong-CP) non-perturbatively from the gluing/cohomology structure.

## What the Model Derives

### Direct Consequences of Patch Consistency

These results follow from the absence of a privileged "third-party" frame. There are only observer perspectives that must agree where they overlap:

**No preferred reference frame**: There is no objective vantage point from which to define absolute rest. Every observer's perspective is equally valid. This follows directly from there being no "view from nowhere." Crucially, there is no observer *outside* the sphere looking in. Observers are patterns *within* the screen data, each with access only to their own patch. The global state exists mathematically, but no entity within the model can observe it.

**Lorentz invariance**: Different observers through the same point must have consistent descriptions. The transformations relating their perspectives form a group. On S^2, this group is Conf(S^2) ~ SO(3,1), the Lorentz group. Special relativity emerges as a consistency condition. Note that the qubits on the screen do not "move." What we call motion in the emergent 4D spacetime is a pattern in how correlations change. A Lorentz boost is a transformation relating how two observers describe the same correlation pattern. The substrate is not in spacetime; spacetime emerges from how patches relate.

**Gauge symmetry as overlap redundancy**: When patches overlap, there is redundancy in how they identify shared observables. This redundancy has the mathematical structure of a gauge group. The photon and graviton are massless because mass terms would break this redundancy, making consistent gluing impossible.

**4D from 2D**: How does a 2D screen give rise to 4D spacetime? The conformal group of S^2 is isomorphic to the Lorentz group of 3+1 dimensions. The screen's symmetries *are* spacetime symmetries. The extra dimensions emerge from how screen regions relate through entanglement.

### Einstein's Equations

**Einstein's equations** are derived from the following chain:

1. *MaxEnt selection* (Assumption B): The global state maximizes entropy subject to local constraints.

2. *Local Gibbs form* (Lemma 2.6): Derived from MaxEnt + local constraints.

3. *EFT bridge* (N1-N3): The null-surface modular bridge is derived from core axioms A1-A4 under two testable conditions: (i) null strips qualify as A4 separators, (ii) local finite variation holds.

4. *QNEC/Focusing*: Derived internally via relative entropy monotonicity.

5. *Entanglement equilibrium* leads to a conditional scaling-limit Einstein branch via the Jacobson mechanism (1995, 2016).

The remaining explicit assumptions are MaxEnt, the OPH Bisognano--Wichmann geometric-branch/scaling-limit premises behind the BW step, and the EFT-bridge conditions listed above. In this formulation, the BW branch premise supplies the standard cap-preserving conformal dilation, while the BW theorem fixes its modular normalization and controlled collar remainder on that branch.

### Particles and Forces

The gauge-as-gluing principle (Section above) reconstructs compact gauge groups from edge-sector fusion rules. Specific consequences:

**Massless photon**: The electromagnetic U(1) symmetry emerges as the redundancy structure for charged observables. A photon mass would break this redundancy, making consistent gluing impossible. Therefore the photon mass is exactly zero. Current experimental upper bounds are consistent with this prediction, constraining any deviation from zero at roughly the 10⁻¹⁸ eV level (often summarized as ~27 orders of magnitude).

**Massless graviton**: Diffeomorphism invariance (the symmetry underlying general relativity) emerges from the fact that bulk spacetime is a compression of screen data. A graviton mass would break this invariance. The graviton mass is exactly zero, with current observational upper bounds consistent with this prediction (e.g. m_g ≲ 1.8×10⁻²³ eV from GW dispersion; often summarized as ~22 orders of magnitude).

**Three colors**: The Witten anomaly requires an odd number of quark colors. Under the Selection Axiom MAR, $N_c = 1$ fails admissibility (trivial color), and MAR selects the minimum: $N_c = 3$.

**Three generations**: CP violation requires at least three generations; asymptotic freedom of SU(2)_L allows at most five. MAR selects the minimum: $N_g = 3$.

**Proton stability**: This is where our model diverges from traditional Grand Unified Theories. Standard GUTs embed SU(3) x SU(2) x U(1) into a larger simple group, introducing leptoquark bosons that mediate proton decay. Super-Kamiokande's limit $\tau_p > 10^{34}$ years excludes minimal SU(5).

Our model takes a different path. Under the extended MAR admissibility package, the minimal coupled carrier $\mathbb{C}^3 \otimes \mathbb{C}^2$ enforces a *product* gauge structure. There is no embedding in a larger simple group and no leptoquark generators. See [the gauge derivation source](paper/tex_fragments/GAUGE_GROUP_DERIVATION.tex) for the complete proof.

## How the Model Explains Known Physics

To be clear: most of what follows are *postdictions*. The photon mass, the number of quark colors, the structure of spacetime, etc. were known long before this model existed. The question is whether the model *explains* them from deeper principles, rather than taking them as inputs.

### Derived from Framework Assumptions

These results follow from the core axioms (A1-A4) plus additional assumptions (MaxEnt, gauge-as-gluing, and the OPH geometric-branch/scaling-limit premises). The full required assumption set is detailed in the technical paper.

| Result | Derivation Chain |
|--------|------------------|
| Photon mass = 0 | Given emergent gauge redundancy + unbroken U(1)_em, mass term forbidden |
| Graviton mass = 0 | Given entanglement equilibrium + diffeo invariance, mass term forbidden |
| Gluon mass = 0 | Given gauge-as-gluing for SU(3), mass term forbidden |
| Lorentz group | A1-A4 + BW theorem on the OPH Bisognano--Wichmann geometric branch gives Conf(S^2) ~ SO(3,1) |
| CPT invariance | Lorentz kinematics + locality via CPT theorem |
| Charge conservation | Unbroken U(1)_em gauge symmetry |
| Newton's constant formula | $G = a_{\rm cell}/4\bar{\ell}(t)$ from edge entropy density |

### Derived Given Assumed Matter Content

These follow once the Standard Model matter content is assumed:

| Result | Status |
|--------|--------|
| Hypercharges (exact rationals) | SM matter content assumed |
| Charge quantization | Z6 quotient from realized spectrum |
| Z6 congruence rule | SM global group structure |
| Conditional center-label entropy = log2 6 bits | Uniform sixfold center-label ensemble on the realized Z6 quotient sectors |
| Yukawa hierarchy $y_f \propto 6^{-n_f}$ | Z6 defect suppression + integer charges |

### Quantitative Outputs vs Existing Data

The complete derivation chain from pixel area $P = 1.63094$ to particle masses is documented in **[the spectrum derivation source](paper/tex_fragments/SPECTRUM_DERIVATION.tex)**. Outputs are classified by epistemic tier there. Additional results (dark matter, Koide formula, baryogenesis, proton spin) are in the **[technical supplement source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex)**.

**Calibration-sector consistency checks** (gauge couplings used to fix $P$; agreement expected by construction):

| Quantity | OPH | PDG | Rel. Error | Status |
|----------|----:|----:|-----------:|--------|
| $\alpha_s(m_{Z,\rm run})$ | 0.1183 | 0.1179 ± 0.0009 | +0.37% | Consistency check |
| $\sin^2\theta_W(m_{Z,\rm run})$ | 0.2307 | 0.23122 ± 0.00004 | −0.21% | Consistency check |
| $\alpha_{\text{em}}^{-1}(m_{Z,\rm run})$ | 128.31 | 127.952 ± 0.009 | +0.28% | Consistency check |
| W boson | 80.386 GeV | 80.377 ± 0.012 GeV | +0.012% | Consistency check |
| Higgs VEV | 246.77 GeV | 246.22 GeV | +0.22% | Consistency check |

**Independent quantitative branch** (main non-calibration mass outputs):

The displayed Higgs/top values are the synchronized D11 outputs. The literal low-scale all-SM appendix flow lands much lower, at about $m_H \approx 115.18$ GeV and $m_t^{\rm pole} \approx 164.25$ GeV. The published $126.48 / 171.1$ branch uses UV-synchronized transport with SM running below $\mu_{\rm sync}\approx 6.8\times 10^{11}$ GeV, MSSM-like running above, and only order-one-GeV residual matching corrections.

| Quantity | OPH | PDG | Rel. Error | Origin |
|----------|----:|----:|-----------:|--------|
| Higgs mass | 126.48 GeV | 125.20 ± 0.11 GeV | +1.02% | D11 supplement reconstruction after $\lambda = \beta_\lambda = 0$ |
| Top quark (crit. surf.) | 171.1 GeV | 172.57 ± 0.29 GeV | −0.87% | D11 supplement reconstruction after $\lambda = \beta_\lambda = 0$ |

**Phenomenological continuations** (reported for completeness, but weaker than the theorem package and the Higgs/top branch):

| Quantity | OPH | PDG | Rel. Error | Origin |
|----------|----:|----:|-----------:|--------|
| Electron | 5.109 × 10⁻⁴ GeV | 5.110 × 10⁻⁴ GeV | −0.023% | Koide continuation ($Q=2/3$, $\delta=2/9$) |
| Muon | 0.10564 GeV | 0.10566 GeV | −0.022% | Koide continuation |
| Tau | 1.7766 GeV | 1.7769 GeV | −0.020% | Koide continuation |
| Top quark (Z₆ texture) | 174.5 GeV | 172.57 ± 0.29 GeV | +1.1% | $\mathbb{Z}_6$ texture ($n_t = 0$) |

**Program-level / supplementary continuations** (from the [technical supplement source](paper/tex_fragments/TECHNICAL_SUPPLEMENT.tex); not all are at the same status as the core paper results):

| Quantity | OPH | Observed | Agreement |
|----------|-----|----------|-----------|
| MOND acceleration scale $a_0$ | $1.03 \times 10^{-10}$ m/s² | $1.2 \times 10^{-10}$ m/s² | 15% |
| Koide ratio $Q$ | 2/3 (continuation ansatz) | 0.666664 | 10⁻⁵ |
| Koide phase $\delta$ | 2/9 (continuation ansatz) | 0.222225 | 10⁻⁵ |
| Baryon asymmetry $\eta_B$ | $4.6 \times 10^{-10}$ | $6.1 \times 10^{-10}$ | Factor 1.3 |
| Proton spin fraction $\Delta\Sigma$ | 0.308 | 0.29 ± 0.03 | Within 1σ |
| MSSM $\beta$-shift ratio $\Delta b_3/\Delta b_2$ | 0.91 | 0.96 (MSSM) | 5% |

**Structural predictions:**

| Prediction | Experimental Status | Notes |
|------------|---------------------|-------|
| Z6 charge quantization | PDG: $\|q_p + q_e\|/e < 10^{-21}$ | Structural prediction |
| Casimir log-gap ratios | Lattice SU(3) confirms 9/4, 5/2, 4, 9/2, 6 | Parameter-free ratios |
| Photon mass = 0 | PDG: $m_\gamma < 10^{-18}$ eV | Exact (gauge symmetry) |
| Graviton mass = 0 | PDG: $m_g < 1.76 \times 10^{-23}$ eV | Exact (diffeomorphism) |

### Derived Under Extended Theory ($T_{\text{ext}}$: A1–A4 + R0 + R1 + [z]=0 + MAR)

These are derived from the Selection Axiom MAR:

| Result | How Derived |
|--------|-------------|
| Product gauge group | Minimal coupled carrier $\mathbb{C}^3 \otimes \mathbb{C}^2$ (MAR) |
| SM global group SU(3) × SU(2) × U(1)/Z₆ | MAR + admissibility conditions |
| $N_c = 3$ (three colors) | Witten anomaly + MAR minimality |
| $N_g = 3$ (three generations) | CP + asymptotic freedom + MAR minimality |
| Proton stability | Product structure from MAR (no leptoquark generators) |
| No magnetic monopoles | Product structure (no GUT-scale symmetry breaking) |

### Consistency Checks (Not Novel Predictions)

These match known results but were derived elsewhere first:

| Result | Original Source |
|--------|-----------------|
| $(b_2, b_3) = (1, -3)$ from Dynkin indices at $t^* \approx 4/3$ | New result (reduces MSSM beta to collar geometry) |
| $\alpha_s(M_Z) \approx 0.118$ with MSSM-like spectrum | MSSM GUT analyses (1990s) |
| $\sin^2\theta_W(M_U) = 3/8$ | Georgi and Glashow (1974) |
| Witten anomaly constraint | Witten (1982) |
| GIM mechanism (no tree-level FCNC) | Glashow, Iliopoulos, Maiani (1970) |

### Structural Features

| Observation | How the Model Explains It |
|-------------|---------------------------|
| Holographic entropy bounds | Information lives on boundaries where observer patches meet |
| Bell violations <= Tsirelson bound | Quantum correlations are maximal consistent correlations |
| No preferred reference frame | No privileged observer exists to define one |
| Measurement affects outcomes | "Measurement" is observer synchronization |
| 4D spacetime from 2D screen | Conf(S^2) ~ SO(3,1); screen symmetries *are* spacetime symmetries |

### Novel Testable Predictions

| Prediction | How to Test | Status |
|------------|-------------|--------|
| Casimir log-gap ratios (9/4, 5/2, 4, 9/2, 6) | Future lattice measurements of SU(3) edge-sector probabilities | Parameter-free from group theory |
| Conditional Z₆ center-label entropy = log₂ 6 ≈ 2.585 bits | Edge-sector entropy measurement under the uniform center-label ensemble ansatz | Global-structure observable |
| BH emission line ratios $E_k/E_2 = \ln k / \ln 2$ | PBH gamma-ray analysis | Parameter-free arithmetic pattern |
| Mass-independent fractional linewidth ~3-5% | PBH emission line profiles | Shape prediction, not just positions |
| Yukawa exponents $-\ln y_f / \ln 6$ near integers | Extract from fermion masses | Phenomenological continuation tied to Z₆ structure |
| Coupling unification without proton decay | Proton lifetime + precision coupling data | Geometric vs algebraic unification |
| GR deviation bound from Markov defect | Precision gravity tests in low-CMI regime | Exponential suppression with collar width |
| Neutrino masses: $m_{\nu_3} \approx 3.0$ meV, $m_{\nu_2} \approx 0.50$ meV, $m_{\nu_1} \approx 0.084$ meV | JUNO, DUNE, KATRIN, cosmological $\sum m_\nu$ | Current order-of-magnitude capacity-branch estimate; normal ordering with $\sum m_\nu \approx 3.6$ meV |
| String theory from edge sectors | Theoretical verification | Mathematical bridge to 2D Yang-Mills / Gross-Taylor is established; critical-superstring completion remains extension work |

### Dynamical Predictions

These predictions follow from the discrete area spectrum under the assumption that integer-multiplication transitions dominate in Hawking emission. Measurement contradiction would rule out the specific dynamical selection rule, not necessarily the framework.

| Prediction | How to Test | Required Assumption |
|------------|-------------|---------------------|
| GW horizon spectroscopy comb | Stack LIGO/Virgo events at $x_k = \ln k / 8\pi$; absence of coherent peaks contradicts | Integer-k transitions dominate |
| Discrete Hawking comb | PBH gamma-ray bursts should show $E_k/E_2 = \ln k / \ln 2$ structure | Integer-k transitions dominate |

## Empirical Contradiction Criteria

OPH makes a mix of exact structural claims, calibrated consistency checks, and weaker branch-level continuations. The following observations provide direct contradiction criteria against the core framework or the specific branches noted:

### Physical Falsifiers

- **Nonzero photon or graviton mass**: Any confirmed mass, however small, would break the symmetry structure the model requires.

- **Proton decay via gauge bosons**: Detection of proton decay through X or Y boson exchange would indicate the gauge group is embedded in a larger simple group, contradicting the model's prediction of a product structure.

- **Information exceeding Bekenstein bound**: If more information is stored in a region than the holographic bound permits, the screen picture would be wrong.

- **Bell violations exceeding Tsirelson bound**: Quantum mechanics saturates the Tsirelson bound. Stronger-than-quantum correlations would indicate reality has structure beyond what the model describes.

- **Unitarity violation in black hole evaporation**: If information is genuinely lost in black holes, the consistency structure would fail.

- **Absence of GW horizon spectroscopy comb**: After rescaling by remnant mass and spin, LIGO/Virgo events should show coherent spectral features at universal coordinates $x_k = \ln k / 8\pi$. Absence of this pattern in a sufficiently large dataset identifies a measurement contradiction with the log-integer area spectrum.

- **Wrong Casimir ratios in lattice data**: If future precision lattice measurements show $\Delta_8/\Delta_3 \neq 9/4$ (e.g., closer to 2.67 or 5.06), the heat-kernel edge mechanism would be contradicted. The full prediction set is $\Delta_8/\Delta_3 = 9/4$, $\Delta_6/\Delta_3 = 5/2$, $\Delta_{10}/\Delta_3 = 9/2$, $\Delta_{15}/\Delta_3 = 4$, $\Delta_{27}/\Delta_3 = 6$.

- **Wrong Z₆ center-label entropy**: If edge-sector entropy measurements disagree with the conditional $\log_2 6$ center-label value under the same ensemble ansatz, the Z₆ continuation would be contradicted.

- **Non-integer Yukawa exponents**: If the extracted exponents $-\ln y_f / \ln 6$ show no clustering near integers across fermion generations, the Z₆ defect mechanism for mass hierarchy would be contradicted.

### Conceptual Falsifier

- **Proof that objective reality is necessary**: The model's central claim is that observer consistency is sufficient. No objective reality independent of observers is needed. If it is demonstrated that some physical phenomenon *requires* a world-in-itself that exists independently of all observation, the model would be contradicted.

This is worth elaborating. The question "does objective reality exist?" has been debated for millennia, but modern physics has sharpened it. Bell's theorem (1964) showed that quantum correlations cannot be explained by local hidden variables (pre-existing definite values that observers passively discover). The Kochen-Specker theorem (1967) showed that quantum observables cannot all have simultaneous definite values independent of measurement context. The PBR theorem (2012) constrained the space of "epistemic" interpretations where the quantum state merely represents knowledge.

None of these theorems *prove* that objective reality does not exist. But they progressively narrow the space for it. Each shows that another intuitive feature of objective reality (locality, non-contextuality, epistemic hiddenness) is incompatible with quantum mechanics.

The model takes the next step: rather than defending an increasingly constrained objective reality, it abandons the assumption entirely. Subjective perspectives and their consistency conditions are all there is. The "objective world" is a useful approximation (the overlap-consistent backbone that all observers share) but it is not a fundamental feature of reality.

A measurement contradiction would require showing that some phenomenon cannot be explained by observer consistency alone. No such phenomenon is currently known.

## Current Status

As summarized in the status matrix above, the framework derives Lorentz kinematics, a conditional semiclassical Einstein branch, compact gauge symmetry, the SM gauge group (via MAR), three generations, three colors, massless force carriers, a two-input quantitative particle-physics program, proton stability, and a continuation-level OPH-to-string route through 2D Yang-Mills edge sectors. The main engineering extensions are:

1. **Screen microphysics**: Quantum link models realize the regulator premises and give edge-center completion automatically. Remaining: force or verify selection of the OPH geometric modular branch in the continuum limit.

2. **EFT bridge verification**: The null-surface modular bridge is derived from core axioms under two testable conditions. These need verification in explicit UV regulators.

3. **Cosmological constant**: Structurally explained as $\Lambda = 3\pi/(G \cdot \log \dim \mathcal{H}_{\rm tot})$ from screen capacity. The numerical value is inferred from observation.

4. **Strong CP problem**: a gluing-obstruction route is outlined, but $\theta_{QCD}$ remains a continuation-level proposal rather than a recovered-core prediction.

## Code

| Script | Description |
|--------|-------------|
| [oph_predict_compare.py](code/particles/oph_predict_compare.py) | Full particle spectrum prediction + PDG comparison (main entry point) |
| [particle_masses_stage5.py](code/particles/particle_masses_stage5.py) | Main spectrum entry point: gauge closure, delegated paper-driven D10/D11 critical-surface branch, Z₆ texture |
| [particle_masses_paper_d10_d11.py](code/particles/particle_masses_paper_d10_d11.py) | Paper-driven D10/D11 reconstruction: literal appendix flow, UV-synchronized transport, and supplement matching layer |
| [oph_qcd.py](code/particles/oph_qcd.py) | 4-loop MSbar QCD running and Λ extraction |
| [oph_lattice_su3_quenched_v5.py](code/particles/oph_lattice_su3_quenched_v5.py) | Quenched Wilson SU(3) lattice for hadron mass ratios |
| [oph_no_cheat_audit.py](code/particles/oph_no_cheat_audit.py) | Static + runtime anti-leak audit |

### Book Chapters

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
| [13](book/chapter-13-desitter.md) | De Sitter | Our universe's holographic screen |
| [14](book/chapter-14-standard-model.md) | Standard Model | Particles from gluing constraints |
| [15](book/chapter-15-relativity.md) | Relativity | Spacetime from modular time |
| [16](book/chapter-16-matter.md) | Matter | Classical physics as emergent stability |
| [17](book/chapter-17-darwin.md) | Darwin's Laws | Laws as evolutionary survivors |
| [18](book/chapter-18-synthesis.md) | Synthesis | Putting it all together |
| [19](book/chapter-19-metaphysics.md) | Metaphysics | Qualia and the hard problem |
| [Epilogue](book/epilogue.md) | Epilogue | One last surprise |

## Building Your Own Reality Simulator

Here's a thought that brings the project full circle: once you've reverse engineered how something works, you can build it yourself.

The screen model is a specification for a particular kind of quantum system: a gauge-invariant lattice on a closed 2D surface with the right entanglement structure. In principle, we could build one.

Current quantum computing architectures are approaching the required ingredients. Topological quantum computing, particularly Microsoft's approach using Majorana fermions, naturally implements non-abelian anyons that could realize the edge sectors our model requires. Fractional quantum Hall states already exhibit the kind of topological order and boundary physics that the screen demands. Quantum link models on triangulated surfaces have been simulated on trapped-ion and superconducting platforms.

A small-scale implementation wouldn't recreate our universe (the screen has $\sim 10^{122}$ degrees of freedom). But it could test whether the claimed emergent properties actually emerge: whether gauge symmetry appears at patch overlaps, whether something like spacetime geometry reconstructs from entanglement, whether observer-like subsystems can form and synchronize.

Where would patches be? They would emerge from the dynamics. The system evolves under a local Hamiltonian (plaquette terms creating and destroying flux loops, electric field terms penalizing large fluxes, gauge constraints at every vertex). Under this evolution, stable correlation patterns spontaneously form. Regions with high internal mutual information but lower correlation to the outside become natural patches. We wouldn't draw them beforehand; we'd watch them self-organize.

Note that Hamiltonian evolution is not the same as the time experienced "inside". The microscopic dynamics tick away according to our external laboratory clock. But any observer-like pattern that emerges within the simulation has its own modular flow, its own internal sense of time constructed from the density matrix on its patch. The lab clock and the simulated observer's clock are related but distinct. We'd be watching from outside as the qudits evolve, while any emergent observers inside would experience their own temporal flow, potentially at a completely different rate or structure.

This is established, but it's not science fiction. The physics is well-defined enough to simulate. The ultimate validation wouldn't be a mathematical proof. It would be switching on a small holographic screen and watching spacetime flicker into existence.

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
