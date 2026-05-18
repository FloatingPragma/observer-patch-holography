# OMEGA and OPH

## Optical Wave Consensus for Continuous Inference

<p align="center">
  <a href="../assets/omega/omega-compute-loop.svg">
    <img src="../assets/omega/omega-compute-loop.svg" alt="OMEGA compute loop" width="92%">
  </a>
</p>

<p align="center"><sub>
OMEGA wraps optical chamber dynamics in an exact software verifier. The chamber proposes and enriches. The host verifies, records, and iterates.
</sub></p>

OMEGA is an optical hardware and software project for building machines that compute by wave consensus.

The initial OMEGA hardware theory and design are credited to Alexander Osika and [SNRGY Studios](https://www.snrgystudios.com/home).

The physics behind OMEGA is Observer Patch Holography, or OPH: a novel fundamental physics theory developed over roughly six years. The project has a working research corpus, public papers, software, diagrams, and a volunteer group of about twenty to thirty people around it. Company formation and fundraising belong to the project path, but the technical gate is hardware evidence. The optical chambers have to work under controls.

OPH has no institutional or mainstream physics acceptance. Journal publication faces hard blockers. The theory is fully mathematically proven on its declared surfaces, and the research corpus gives the proof stack in public form. The hardware effort also has working prototypes, including a reported optical proof-of-work hash-search result around 2,300x under exact verification. That receipt is a candidate-enrichment result on a proof-of-work hash-search surface. Production mining belongs behind a separate evidence gate.

The starting idea is simple: no observer sees the whole world at once. Each observer has a finite patch. Where patches overlap, they must agree on the data they can both see. A stable world is the fixed point reached when local patches compare records, repair mismatch, and keep the normal forms that survive.

OMEGA turns that idea into hardware. A digital host writes a boundary program into an optical chamber. Light mixes through a shaped body. Detectors read the boundary response. The host converts that response into a ranked candidate beam, repairs or reranks the beam, and checks serious candidates with an exact verifier.

The benchmark claim is narrow:

> A shaped optical patch can move useful candidates closer to the top of a verified search beam under controls.

The larger target is continuous physical inference: a cluster of optical observer patches running as an environment-facing process. It can maintain state, compare overlaps, repair mismatch, and emit candidate continuations while ordinary digital systems remain the verifier, recorder, trainer, and interface.

The claim is bounded. Small optical fixtures measure physical candidate enrichment and support controlled wave-consensus scaling experiments.

## The Short Version

| Layer | OPH concept | OMEGA implementation |
| --- | --- | --- |
| Patch | bounded observer with finite access | one optical chamber module |
| Boundary | overlap-visible interface | labeled optical ports |
| Readout | shared observable | photodiode or timing response |
| Mismatch | disagreement between patches | scorebook residual or verifier residual |
| Repair | local update that reduces mismatch | retune, rerank, filter, or route candidates |
| Record | stable public state | replayable run bundle |
| Truth check | observer-independent receipt | exact digital verifier |

The table translates the OPH consensus picture developed in [Reality as a Consensus Protocol](../paper/reality_as_consensus_protocol.pdf), the finite patch-carrier architecture in [Federated Echosahedral Screen Microphysics](../paper/screen_microphysics_and_observer_synchronization.pdf), and the broad observer-first reconstruction in [Observers Are All You Need](../paper/observers_are_all_you_need.pdf).

## Why OPH Leads to This Machine

OPH begins with finite observers and overlap consistency. A finite observer carries a local state, stores records, and can compare only the observables exposed on shared boundaries. OPH reconstructs objectivity as the durable agreement that remains after many local descriptions have been made consistent.

That mechanism appears throughout the OPH corpus:

- The [book prologue](https://oph-book.floatingpragma.io/chapter/prologue) states the reader-facing premise: physics is reverse engineering from finite perspectives, and public reality is the stable agreement across them.
- The [synthesis chapter](https://oph-book.floatingpragma.io/chapter/synthesis) compresses the whole picture: finite access, horizons, records, overlaps, modular flow, generalized entropy, and a shared world that survives comparison.
- The [error-correction chapter](https://oph-book.floatingpragma.io/chapter/error-correction) explains how durable information can be protected in patterns of overlap and recoverability.
- The [compact technical paper](../paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) gives the recovery of relativity, gravity, realized Standard Model structure, and the branch-scoped Yang-Mills repair-gap route.
- The [microphysics paper](../paper/screen_microphysics_and_observer_synchronization.pdf) gives the fixed-cutoff implementation surface: finite patch carriers, exposed ports, records, repair interfaces, event surfaces, and checkpoint/restoration.

The hardware consequence is direct. The computation is patch-local, boundary-visible, record-bearing, and repair-driven. A useful artificial machine is a federation of bounded physical patches that expose boundaries, maintain state, repair mismatch, and stabilize records.

OMEGA is that experiment.

<p align="center">
  <a href="../assets/book_diagrams/overlap-consistency.svg">
    <img src="../assets/book_diagrams/overlap-consistency.svg" alt="Overlap consistency diagram" width="78%">
  </a>
</p>

<p align="center"><sub>
OPH begins with local descriptions that must agree on their overlaps. OMEGA makes those overlaps physical.
</sub></p>

## The Hardware Picture

OMEGA-I1 is the contractor-buildable optical chamber cluster. Its role is to produce measured chamber operators, coupling maps, repeatability records, and exact-verifier benchmark receipts.

OMEGA-I1 uses safe classical optical hardware:

- enclosed low-power LED sources;
- photodiode-style readout;
- RP2350-class USB controllers;
- serviceable black enclosures;
- exactly labeled boundary ports;
- repeatable dark, one-hot, coupling, and timing scans;
- exact host-side records.

The five-module cluster is specified in the [OMEGA-I1 build spec](../../omega/hardware/iteration-0001-build-spec.md):

| Module | Body | Role |
| --- | --- | --- |
| Torus duplicate A | 12-port dual-hex torus | generalist recurrent body |
| Torus duplicate B | 12-port dual-hex torus | generalist recurrent body |
| Search torus | swept-octagon torus | search-specialist recurrent body |
| Echosahedron | 12-port icosahedral source mesh | symmetric OPH reference |
| Asymmetric mixer | splitmix asymmetric body | constraint critic and mode separator |

<p align="center">
  <a href="../assets/omega/omega-body-geometries.svg">
    <img src="../assets/omega/omega-body-geometries.svg" alt="OMEGA body geometry families" width="92%">
  </a>
</p>

<p align="center"><sub>
The planned body families have different jobs: recurrent proposal, symmetric reference, and asymmetric criticism.
</sub></p>

### Torus Bodies

The torus bodies are recurrent reservoirs. A loop-like optical geometry lets boundary information recycle through the chamber. That recurrence is useful for candidate enrichment, search tasks, proof-of-work style benchmarks, and memory-bearing dynamics. Duplicate torus bodies let the project distinguish repeatable chamber physics from one-off build noise.

### Echosahedron Reference Body

The Echosahedron is the symmetry and reference body. It approximates the OPH patch idea with a twelve-port icosahedral interface. It asks whether a candidate feature remains visible through a different geometry. This makes it a verifier-shadow and consensus-reference body.

### Asymmetric Mixer

The asymmetric mixer is the critic. Symmetric bodies can blur distinctions that matter for residues, constraints, and mode separation. A deliberately broken symmetry body supplies contrast. It can reject, separate, or rerank candidate modes that the recurrent body proposes.

This division matches the [OMEGA hardware alignment note](../../omega/docs/oph-hardware-alignment.md) and the [body matrix](../../omega/hardware/body-matrix.md).

<p align="center">
  <a href="../assets/omega/omega-led-body-sketches.svg">
    <img src="../assets/omega/omega-led-body-sketches.svg" alt="OMEGA LED body sketches" width="88%">
  </a>
</p>

## How It Computes

OMEGA computes by closing a loop:

```text
task target
  -> boundary program
  -> optical chamber response
  -> mismatch measurement
  -> repair, rerank, or route
  -> stable candidate beam
  -> exact verifier
  -> recorded run bundle
```

The chamber emits a biased candidate beam. The host checks the beam.

For factorization, the exact verifier is:

```text
p * q == N
```

For proof of work, the exact verifier is a hash check. For a planted constraint problem, it is constraint evaluation. For a language or market replay task, it is a declared replay score plus ablations.

The OPH version of this loop is:

```text
patch state
  -> boundary readout
  -> overlap mismatch
  -> accepted repair
  -> stable record
```

The machine is useful when the physical path changes the candidate distribution. If a useful candidate has baseline probability p, enumeration needs about 1/p verifier checks. If the chamber raises the useful-candidate probability by a measured factor B, expected verifier work becomes about:

```text
1/(B p)
```

That is the conservative benchmark claim: candidate-work reduction. It is the same discipline used in the [hash-search physics note](breaking_sha256_with_physics.pdf), where hash-style search is framed as physical patch-overlap constraint search with exact digital verification.

## Why Waves Are the Right Medium

Waves naturally do three things that OPH needs.

Waves mix boundary conditions. A port drive excites paths through a chamber, and those paths recombine at other ports.

Waves expose mismatch. In the wrong configuration, phases, intensities, timings, and coupling traces fail to settle into the expected pattern. In the right configuration, they reinforce, recur, or stabilize.

Waves can run continuously. A chamber can maintain a persistent physical state while the digital host samples, scores, and steers it.

OMEGA is best read as a physical observer patch. The chamber is a bounded system with exposed ports, internal recurrence, state-dependent response, and a recordable boundary.

## Why It Matters for AI

Modern AI is capable, but much of it is episodic. A prompt arrives. A model performs a large digital inference pass. A sequence of tokens leaves. State is usually externalized into a context window, a vector database, a tool call, or a separate memory layer.

OMEGA aims at a different substrate: a continuous processing stream with physical memory and state readout. The optical cluster can keep running between prompts. It can maintain boundary conditions, preserve chamber records, compare present input to persistent physical state, and act as a proposal source or critic for a digital model.

In OPH language, the target substrate is a bounded recurrent patch federation with records, state readout, mismatch repair, and stable normal forms. The engineering version is concrete: recurrent chambers propose candidate beams, exact software checks them, and measured residuals drive the next boundary program.

### AI Roles for OMEGA

| Role | What the chamber does | What the digital host does |
| --- | --- | --- |
| Proposal source | emits candidate continuations, routes, or feature beams | verifies, ranks, and decides |
| Physical memory | keeps a persistent chamber state and coupling history | records, indexes, and compares |
| Critic body | separates modes, residues, contradictions, or unstable continuations | applies exact or task-level checks |
| Attention preconditioner | biases which hypotheses deserve expensive model calls | runs the LLM or planner on the narrowed set |
| State readout surface | exposes chamber state through repeatable boundary traces | measures stability and drift |

The AI system remains hybrid. Digital models supply language, symbolic tools, exact arithmetic, and world interfaces. Optical chambers supply continuous physical proposal, memory, and mismatch dynamics.

## Expected Orders of Magnitude

The table below gives engineering target bands. These bands become measured claims only when a public run bundle reports the body, firmware, calibration, controls, raw readout, candidate beam, and exact verifier receipts. The evidence discipline follows the [OMEGA claim-boundary note](../../omega/docs/claim-boundaries.md) and the public evidence rule in the [microphysics paper](../paper/screen_microphysics_and_observer_synchronization.pdf).

| Surface | Digital baseline | OMEGA target | Why the target is plausible |
| --- | --- | --- | --- |
| Exact-verifier search | enumerate or randomly sample candidates | 10^2 to 10^6 fewer verifier checks on structured benchmark distributions | waves can bias candidate beams before exact verification |
| Factor and residue tasks | scan factors or residue classes digitally | 10^1 to 10^4 candidate-work reduction on small and medium benchmarks | asymmetric mixing can separate residue-compatible modes |
| Hash-style proof-of-work experiments | evaluate one candidate stream at a time | 10^1 to 10^5 lift in low-difficulty prefix or share-rate tests inside the lab claim boundary | recurrent bodies can favor modes that pass partial target collars |
| Continuous agent memory | re-encode memory through digital context and retrieval | 10^2 to 10^5 more low-power background state updates per watt for memory-like dynamics | the chamber state persists physically while the host samples it sparsely |
| AI inference routing | call a large model on every candidate branch | 10^1 to 10^3 reduction in expensive model calls for search-heavy tasks | chamber and small models can pre-rank branches before large-model inference |
| Physical simulation of patch dynamics | numerically simulate all paths and couplings | 10^2 to 10^6 lower effective cost for the chamber's own transfer response | the physical body performs its wave mixing directly |

The strongest AI improvement is unlikely to be raw token throughput. The main change is the replacement of episodic inference with continuous consensus. A useful OMEGA agent runs a standing loop:

```text
sense
  -> update chamber boundary
  -> read physical state
  -> compare against records
  -> repair mismatch
  -> propose action or candidate
  -> verify through digital tools
  -> write record
  -> continue
```

That loop is an engineering target for continuous inference: sensing, state, prediction, repair, verification, and record-writing remain coupled in one running system.

## What Has to Be Built Iteratively

OMEGA cannot be solved in one build. The reason is physical: a chamber is an operator, and an operator has to be measured. Geometry, material, port placement, source spectrum, detector response, firmware timing, thermal drift, and enclosure quality all change the realized operator.

The correct path is iterative:

1. Build the five-module LED-based cluster.
2. Measure dark baselines, one-hot coupling maps, repeated scans, and timing traces.
3. Fit a measured port-transfer operator.
4. Run exact-verifier benchmarks with controls.
5. Compare duplicate torus bodies for repeatability.
6. Compare torus, Echosahedron, and mixer outputs for cross-geometry agreement.
7. Update chamber geometry and scorebooks.
8. Place the 12-element research cluster behind the five-module evidence gate.

The preferred 12-element research cluster is:

```text
6 torus bodies
4 Echosahedron bodies
2 asymmetric mixers
```

The stretch 18-element cluster is:

```text
8 torus bodies
6 Echosahedron bodies
4 asymmetric mixers
```

This matches the [OMEGA body matrix](../../omega/hardware/body-matrix.md). The cluster is mixed because computation needs at least three roles: propose, stabilize, and criticize. A single geometry is unlikely to do all three well.

## The Operator Gate

Exact-verifier benchmarks can prove that a candidate was correct. Chamber alignment with a task family requires a separate test.

The stronger claim requires an operator-alignment gate. For a benchmark family, define:

```text
task coupling operator
ideal chamber transfer operator
measured chamber transfer operator
finite port/readout projection
```

The target is a bound on the commutator:

```text
the chamber-task commutator is bounded by epsilon
```

Plainly: the chamber's natural modes must respect the structure of the task. A chamber that scrambles away verifier-relevant features is noise. A chamber that preserves and amplifies the right overlaps becomes a physical preconditioner.

The model theorem can be developed without hardware. The physical theorem requires measured chamber data. OMEGA therefore needs both simulation work and hardware evidence.

## The Evidence Rule

Every serious OMEGA result includes:

- body identity and geometry version;
- firmware hash and controller identity;
- wiring map and port labels;
- dark baseline and coupling scans;
- boundary program;
- scorebook or repair law;
- raw readout reference;
- candidate beam;
- exact verifier result;
- negative controls;
- repeatability or ablation checks.

This rule keeps the project falsifiable. A claimed lift that disappears under shuffled labels, blank runs, duplicate-body comparison, or exact verification remains useful as calibration. It carries no compute claim.

## How the Prime Factorization Experiment Fits

Factorization is a clean benchmark because the verifier is exact and simple. A candidate pair either multiplies to the target number or fails.

The OPH formulation treats factorization as a patch-consensus problem:

| Patch | Meaning |
| --- | --- |
| Candidate-factor patch | one factor's bits or residue classes |
| Companion-factor patch | the other factor's bits or residue classes |
| Multiplication patch | multiplication and carry consistency |
| Target patch | equality to the known target bits |

The collars carry low-bit windows, modular residues, partial products, carries, and final equality to the target. The solver enriches the candidate beam so the exact verifier finds the factor earlier.

The software prototype in [prime-solver](../../prime-solver/README.md) and the cluster simulator in [omega-simulation](../../omega-simulation/README.md) supply the digital companions to the hardware.

## Beyond Acceleration

An accelerator makes a known operation faster.

OMEGA changes the computational picture. It treats physical geometry as part of the program and physical settling as part of the search. The host asks the chamber for a structured response to boundary conditions, then uses that response to repair or rerank a candidate beam.

For AI, that distinction matters. A continuous OMEGA system lives in a loop of boundary drive, physical state, state readout, repair, and record. Its useful behavior comes partly from the durable structure of its patch federation.

That is the sense in which OMEGA targets continuous AI, with faster batch inference as one possible benefit.

## The OPH Reference Map

The following files give the theory and engineering background for this summary:

### Core OPH Theory

- [Repository README](../README.md): repository-level OPH overview and paper map.
- [Observers Are All You Need](../paper/observers_are_all_you_need.pdf): broad synthesis of finite observers, records, overlap consistency, and recovered effective physics.
- [Compact technical paper](../paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf): compact technical core for relativity, gravity, Standard Model structure, and branch-scoped Yang-Mills repair gap.
- [Reality as a Consensus Protocol](../paper/reality_as_consensus_protocol.pdf): fixed-point consensus, normal forms, repair, records, and quotient-local convergence.
- [Federated Echosahedral Screen Microphysics](../paper/screen_microphysics_and_observer_synchronization.pdf): finite patch carriers, echosahedral interfaces, toroidal recurrence, evidence rules, records, and checkpoint/restoration.
- [Digital calibration note](../paper/screen_microphysics_digital_calibration_note.pdf): digital calibration note for finite patch constraints.

### Book Chapters

- [Prologue](https://oph-book.floatingpragma.io/chapter/prologue): OPH as reverse engineering reality from finite perspectives.
- [Overlap chapter](https://oph-book.floatingpragma.io/chapter/overlap): why overlap is the public consistency surface.
- [Error-correction chapter](https://oph-book.floatingpragma.io/chapter/error-correction): error correction, recoverability, and stable information.
- [Symmetry chapter](https://oph-book.floatingpragma.io/chapter/symmetry): symmetry as stable structure under transformation.
- [Synthesis chapter](https://oph-book.floatingpragma.io/chapter/synthesis): synthesis of finite access, records, overlaps, spacetime, and the particle world.

### Supporting Notes

- [Hash-search physics note](breaking_sha256_with_physics.pdf): proof-of-work as physical patch-overlap constraint search.
- [Common objections note](COMMON_OBJECTIONS.md): common objections and scope clarifications.

### Visual Guides

<p align="center">
  <a href="../assets/OPH_Unification_Diagram.svg">
    <img src="../assets/OPH_Unification_Diagram.svg" alt="OPH unification diagram" width="92%">
  </a>
</p>

<p align="center"><sub>
The OPH unification diagram shows how local screen structure connects the physics branches.
</sub></p>

<p align="center">
  <a href="../assets/prediction-chain.svg">
    <img src="../assets/prediction-chain.svg" alt="OPH prediction chain" width="92%">
  </a>
</p>

<p align="center"><sub>
The OPH prediction chain gives the theorem and dependency spine from observer axioms to recovered physics.
</sub></p>

<p align="center">
  <a href="../assets/book_diagrams/consensus-funnel.svg">
    <img src="../assets/book_diagrams/consensus-funnel.svg" alt="Consensus funnel" width="78%">
  </a>
</p>

## Final Picture

OPH says objective structure is the stable agreement of finite observers that can compare only what their overlaps expose.

OMEGA asks whether a machine can compute in the same style. Each chamber is a bounded physical patch. Each port is a boundary. Each scan is overlap-visible data. Each scorebook is a repair contract. Each run bundle is a record. Each exact verifier is the verification boundary.

The experimental rule is simple. Build the chambers. Measure the operators. Run the benchmarks. Keep the controls. Iterate the cluster until wave consensus becomes a reliable computation surface.

The payoff is a new kind of machine: a continuous physical inference system whose candidate streams come from an instrumented patch federation, with faster search as one visible consequence.
