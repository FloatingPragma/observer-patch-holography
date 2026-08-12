# Observer Patch Holography

> Reality is the stable public world reconstructed by finite, self-reading observers that compare their overlaps and repair disagreement.

[Read in French](README_FR.md) · [Technical paper](flagship/from_observer_consensus_to_standard_physics.pdf) · [Textbooks](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [Hardware](https://omega.floatingpragma.io/)

Observer Patch Holography (OPH) is a zero-dial theory-of-everything research
program built on one central thesis: **observers are primary, and objective
reality is emergent.** Physics normally begins by supplying spacetime, quantum
fields, a gauge group, and a table of measured constants. OPH begins with
observers: bounded systems that carry local state, read part of themselves and
their neighbors, keep records, and repair disagreement. Reality emerges from
observer overlap repair on a holographic screen. From this architecture OPH
reconstructs an exact finite structural core: conditional quantum-record
identities, a conditional finite four-law package, a three-dimensional
observer-frame carrier, and explicit order/clock interfaces. It also derives
Lorentz kinematics on the stated global-support branch, the Standard Model
gauge Lie type, and a conditional one-generation exterior matter pair.

Three axioms govern the simulator architecture and how observers reach
consensus. Beside them sit two proposed closure programs. The first seeks a
fixed point for the pixel constant $P$, with an open physical attachment
to the fine-structure constant. The second seeks a fixed point for the
capacity $N$, with an open source-capacity bridge to the cosmological
constant. Identifying the simulated and simulating universe motivates those
self-consistency equations; it does not by itself prove that a solution
exists, is unique, or has the observed numerical value.

## Start Here

Physics has revised its idea of what is fundamental before. Space was
absolute until it was relative; matter was continuous until it was quantized.
Each revision looked outrageous from inside the previous picture and obvious
from inside the next one. OPH makes the next revision. The observer, treated
for a century as a nuisance at the edge of quantum mechanics, moves to the
foundation. Spacetime, matter, and the constants become precise reconstruction
problems, with exact finite results and open physical identifications kept
apart. The material below takes you through that shift from a standing start.

- **The technical paper.** [*From Observer Consensus to Standard Physics*](flagship/from_observer_consensus_to_standard_physics.pdf)
  gives the primary technical account of the observer-first reconstruction.
- **The textbooks.** The [OPH textbooks](https://learn.floatingpragma.io/)
  teach the theory the long way. Every basic derivation is worked in full,
  with the required math built up as you go. Volume one covers the
  computational substrate and the consensus machinery; volume two connects
  that machinery to classical physics. Each is readable online or as a PDF.
- **The simulation.** The [interactive visualizations](https://simulation.floatingpragma.io/)
  render real data from the repair dynamics. They expose finite settling,
  signature tests, and candidate carrier structure, with each finite receipt
  available for direct inspection.

The rest of this README is the technical entrance to the repository.

Two ledgers carry the quantitative record. The
[postdiction ledger](docs/POSTDICTION_LEDGER.md) is the compare-only
scoreboard: every certified comparison against a measured value,
with its premises and input ancestry stated on the row. The
[frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) is the forward
instrument: stances registered with cryptographic custody and kill bands
before their comparison data is examined, with fixed rules that permit
refutation by qualifying measurements.

## One Architecture, All Of Physics

Every mainstream theory starts by assuming most of physics: a spacetime,
quantum fields, a gauge group, a list of measured constants. OPH assumes none
of that. It starts with observers, finite systems that read part of
themselves and their neighbors, keep records, and repair disagreement, and
derives the rest as theorems. From that one architecture:

- **Quantum mechanics as theorems.** On the finite observer surface, public
  records form an event algebra with Born probabilities, Lüders
  conditioning, and the Tsirelson bound. Schrödinger dynamics is the unique
  continuous symmetry flow, and the Born weights follow without a continuity
  axiom in every finite dimension, including dimension two.
- **The four laws of thermodynamics from disagreement repair.** One
  conditional theorem package about how observers resample toward consensus
  yields all four laws, with the second law appearing as data processing
  applied to repair and the Landauer bound as a corollary.
- **Relativity on the screen.** Lorentz kinematics is a theorem on its
  stated branch. Deterministic simulations measure a Lorentzian event form
  with one time and three space directions, and an Einstein-equation
  composition stands on the modular and entropy layer.
- **The Standard Model gauge group from twelve ports.** OPH makes an
  architectural choice at the simulation hardware layer: each observer patch
  has twelve boundary ports wired as the corners of an icosahedron. A classification theorem forces
  the complete port response to have the Standard Model's gauge Lie type,
  with no gauge group chosen from a catalogue, and an exhaustive finite
  search returns the fifteen states and charge pattern of one Standard Model
  generation with exact anomaly cancellation.
- **Constants as fixed points.** The core has zero adjustable parameters.
  Koide's charged-lepton relation holds exactly under a stated balance
  premise, interval arithmetic certifies the tau-mass comparison, and a
  fixed-capacity mechanism gives the de Sitter time-advance sign. The
  constants of nature enter as fixed-point problems to be solved.
- **Machine-checked and falsifiable.** Thousands of Lean theorems with no
  admitted proofs, exact rational arithmetic in place of floating-point
  trust, and deterministic simulations with pinned receipts. A
  frozen-prediction ladder registers kill bands under cryptographic custody
  before comparison data is examined, so OPH commits in advance to what
  would refute it.

Exact finite results and open physical identifications stay strictly
separate across the corpus; every result above carries its premises and
boundary in the linked papers and proofs. The condensed version of this
case, with the receipts and their evidence in one table, is the
[compact case for OPH](extra/compact_proof_of_oph.pdf); the full technical
route is the
[flagship paper](flagship/from_observer_consensus_to_standard_physics.pdf).

The rest of this README is the architecture that case comes from.

## The Three Axioms

The whole construction stands on three core axioms. The canonical statements
live in [the axiom reference](docs/AXIOM_REFERENCE.md) and the machine
registry `claims/axiom_registry.yaml`; the papers include the shared formal
basis.

1. **A1: Oriented twelve-port observer screen.** There exists an observer
   patch net on an oriented spherical screen. At every finite resolution,
   each local carrier has twelve primitive boundary ports forming the
   vertices of an oriented triangular boundary with 30 edges and 20 faces,
   combinatorially the boundary of an icosahedron. Carriers join through
   typed seams and coherent triple overlaps, refine to an oriented spherical
   support, and expose local state, readback, records, repair moves, and
   checkpoints. Formally: for every regulator $r$ there is a typed object
   $\mathfrak N_r=(\mathcal P_r,\mathcal A_r,\mathcal R_r,\mathcal I_r,
   \mathcal U_r,\mathcal C_r,N_r,S_r,b_r)$ whose carriers carry twelve
   primitive central port projections and the exact boundary packet
   $K=(P,E,F,o)$, joined by seam algebras into a nerve with a degree-one
   bridge to the oriented spherical support, all commuting with refinement.
   The local carrier, the federation of carriers, and the global $S^2$
   support stay typed and distinct throughout the corpus.
2. **A2: Observer agreement.** Observers operating on the screen agree on
   the meaning of the data they jointly interpret. Formally: the
   interpretation map $\mathcal J_r$ from observer-accessible data to
   operational meanings is natural with respect to every visible overlap
   restriction, recharting, seam translation, higher-overlap map, federation
   map, and refinement map on accepted public data. No patch sees the whole universe;
   a fact becomes public only when it survives comparison across overlaps.
3. **A3: Conditional maximum randomness.** Everything that observer
   agreement leaves unconstrained is maximally random. Formally: the
   realized state is the information projection of an exact reference family
   onto the convex set of compatible local state families satisfying the
   finite observer-visible constraints. The finite A1-generated observer
   cover is state-determining on that feasible set, and its exact weights are
   strictly positive:
   $\rho_r=\arg\min_{\rho\in\mathcal K_r}\sum_P w_{r,P}
   D(\rho_{r,P}\Vert\tau_{r,P})$.

None of the axioms contains a gauge group, a particle list, a recovery law,
or a rule that selects field content or multiplicity; A3 selects one state
inside one fixed feasible space and nothing else. Collar recovery,
generalized-entropy structure, and
sector completions enter as named interfaces and declarations at the results
that consume them, each classified as an exact theorem, an exact result
inside a named finite realization, a discovery-level observation, a declared
open interface, an independence result with countermodels, a physical
identification, or a withdrawn claim.

Everything else in the repository is the working-out of what these three
axioms force, and of exactly how much further structure each physical
conclusion consumes.

## The Idea In Plain Language

OPH asks: **what is the smallest kind of system capable of having a world at
all?**

The answer is an observer patch. It need not be a person. It is any
bounded physical or computational system that has a local state, a boundary,
memory, the ability to read part of itself and its neighbors, and a way to
repair disagreement. No patch sees the whole universe. A fact becomes
objective only when it can be written, compared across overlaps, recovered
after further evolution, and retained as part of the public record.

OPH treats this process as the mechanism that selects a public physical world.
The theory has no external ruler, master clock, preferred observer,
or list of adjustable physical constants. “Zero dials” means zero fitted
continuous theory values. The finite observer contract and each discrete
branch condition remain visible.

“Observer” is a structural role. A human mind, an organism, an instrument, or a
software process can instantiate it when it has the required state, boundary,
records, readback, and repair loop. OPH does not claim that human thoughts
manufacture reality. It claims that a world with no possible local perspective,
record, or self-consistent readback lacks public physics.

## Claim Scope

The [claim scoreboard](tracking/claims_scoreboard.md) states the scope,
premises, and evidence class of each branch. This README concentrates on the
strongest parts of the case.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## The Two Constants: P and N

**$P$ is the local pixel ratio**: the observation cell's size in natural
units, informally the universe's **resolution**. Two declared trial maps ask
the cell to agree with the observation process it supports:

$$
\boxed{P_\star=\varphi+\frac{\sqrt\pi}{A_T(P_\star)}}.
$$

Each map has one exact interval-certified candidate root. A physical
fine-structure prediction needs a map selected without using the measured
constant, proof that its two sides read one quantity, and same-scheme transport
to the Thomson limit. The numerical match has diagnostic status. See the
[claim scoreboard](tracking/claims_scoreboard.md) for the exact values.

**$N$ is the public-record capacity** of the whole observer system: how much
correctable memory the substrate carries. It sits opposite $P$, tied to the
cosmological constant rather than to the fine-structure constant.

The direct route reads $N$ off the universe itself. The self-read condition
$N=\log M_0(\mathfrak U_N)$ asks the capacity handed to a trial universe to
match the record capacity reconstructed inside it, and if both sides are
readings of one quantity, self-reference forces them to agree. The proof that
they are one quantity does not exist, so this route returns no
number. Nothing else in the reconstruction waits on it.

A second route goes through $P$. At the pixel value supplied to that declared
branch the uncorrected capacity is
$N_0=\pi\exp[6\pi/(P\alpha_U(P))]=3.5321315\times10^{122}$. Two ways of
applying the finite survival correction to it give

$$
N_{\rm pres}=N_0\left(1-\frac{P}{24}\right)=3.2920979\times10^{122},
\qquad
N_{\rm Pois}=N_0e^{-P/24}=3.3000722\times10^{122},
$$

about $0.63$ and $0.39$ percent below the Planck base-$\Lambda\mathrm{CDM}$ comparison
value $3.3129271\times10^{122}$. The theory does not select between the two
corrections, and both numbers were computed after the comparison value was
known, so neither is a prediction. The
[claim scoreboard](tracking/claims_scoreboard.md) states what each step
assumes and what is missing.

## Technical status

The case above is the reader-facing summary. Exact premises,
comparison ancestry, and falsification rules live in the
[claim scoreboard](tracking/claims_scoreboard.md), the
[postdiction ledger](docs/POSTDICTION_LEDGER.md), and the
[frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md). The exact finite
and structural results are the strongest part of the stack. Source-to-physical
attachments, physical scales, and prospective data comparisons form the main
research route.

## Why Take The Claim Seriously?

A successful theory of everything should explain why facts that appear
unrelated arrive as one package. OPH starts from a bounded self-reading patch
instead of a spacetime manifold, field content, gauge group, or table of
constants. It returns exact dimensions, compact Lie types, conditional global
quotients, charge assignments, anomaly cancellations, representation
multiplicities, and fixed-point equations. These outputs come from one typed
carrier, overlap, and repair architecture. The local icosahedral theorem
forces the Standard Model Lie type. The separate compact-sector route reaches
that type only on its declared Standard Model packet, and a common physical
source identity is an open test. Their shared dependence is the main case that
OPH describes one physical world rather than a collection of coincidences.

The evidence also comes in different forms: paper proofs, exact arithmetic,
interval certificates, finite receipts, simulations,
and explicit falsifiers. Agreement among those forms is more informative than
another numerical match produced by another adjustable model.

## Evidence You Can Inspect

The evidence comes in several complementary forms:

- hand proofs in the TeX papers;
- interval and uniqueness certificates for declared numerical maps;
- finite carrier and hierarchy receipts;
- particle, geometry, dark-sector, and quantum-hardware code;
- a small-scale simulation harness that supplies receipts where the hand proofs
  and the Lean development do not reach, in the companion
  [oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim) repository;
- a claim registry connecting prose claims to artifacts.

## Audit The Finite Core

The shortest scientific audit checks the claim graph, the exact twelve-port
algebra, public-record capacity, the reversible $N$ packet, and finite
consensus:

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

The [reproduction guide](REPRODUCE.md) gives the clean-clone setup and the
fuller finite-core lane, which adds the two W/Z convention and
survival-boundary calibration tests.

## The Twist: The Universe Is Its Own Simulator

Everything above stands on the three axioms together with the stated
premises and named interfaces of each result; none of it uses the hypothesis
of this section. That hypothesis arrives as a twist rather than a
foundation. It is itself
an indirect consequence of consistency: something that exists with no outside
support must be capable of creating itself. A completely consistent
observer-built reality must therefore evolve observers, and those observers
eventually build the hardware the reality runs on. The simulated universe and
the simulating universe turn out to be the same system. The patches,
computation, records, and resulting world all belong to one closed loop; no
external computer or programmer appears in the formal construction. The
organizing equation of that closure is

$$
T(\mathfrak U_{\mathrm{OPH}})=\mathfrak U_{\mathrm{OPH}}:
$$

the universe as a fixed point of its own observer-accessible readback and
repair process.

If the loop closes, $P$ and $N$ cannot be arbitrary. They must satisfy
self-referential closure conditions: the cell
must agree with the observation process it supports, and the record capacity
must agree with the records the system keeps about itself. Part of that
closure is machine-checked in Lean. The two declared $P$ maps have certified
fixed points, while their comparison with the physical fine-structure
constant has diagnostic status. The evaluation boundaries of the closure conditions and
their missing physical inputs are stated in the
[OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md).

A physical closure of both constants would give a zero-continuous-parameter
branch with both values returned by the architecture. That physical
attachment is open. The fixed-point theorems certify roots of declared maps;
they do not turn an observed basin or target-defined coordinate into a
physical derivation. On the $N$ side the finite counting is exact, but the
capacity source it would close over is incomplete, so the direct condition is
not evaluable and the common-load route stays conditional on its physical
identifications. Reading $N$ from the universe leaves every consequence of the
three axioms intact.

Under full closure, the loop answers the last question a theory of
everything can be asked: why anything exists, and why it is the way it is.
The universe is the unique structure consistent with reading itself into
existence.

## Open Problems And The Falsification Boundary

The reconstruction runs from the three axioms toward a public quantum
theory, event geometry, a macroscopic spacetime description, and, at the
end of the chain, the Standard Model Lagrangian. Some links in that chain
are proved, some are proved in bounded form, and some are open. Each open
step is a tracked
[research question](https://github.com/FloatingPragma/observer-patch-holography/issues)
with its dependencies, and every claim in the papers carries its own scope
note. A mismatch with the Standard Model at any step is an allowed outcome
that the protocol may not tune away.

The [OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md) lists
the mature claims together with the exact observations that would break
them.

## Choose A Reading Path

| If you want... | Start here |
| --- | --- |
| The flagship introduction to OPH | [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf) |
| The shortest persuasive overview | [A Compact Case for OPH](extra/compact_proof_of_oph.pdf) |
| The spacetime and Einstein derivation | [Recovering Observer Spacetime and Einstein Dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Both Standard Model gauge routes | [Deriving Standard Model Gauge Structure](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| The finite consensus mechanism | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| The particle construction | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| The twelve-port screen architecture and finite modular-gearing theorem | [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Supporting evidence | [`code/`](code) and the [reproduction guide](REPRODUCE.md) |
| Observer continuation and interpretation | [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf) |

The [paper index](paper/) gives the curated publication map. Focused research PDFs remain in [`extra/`](extra/) for repository readers and are not part of the publication release.

## Dependency Map

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH reconstruction chain" width="92%">
  </a>
</p>

<p align="center"><sub>The typed OPH dependency map. It separates exact and conditional branches from the open source, support, current, attachment, and scale bridges that would make them one physical realization.</sub></p>

## Repository Guide

- [`flagship/`](flagship): the primary standalone OPH paper, its TeX source, and release PDF.
- [`paper/`](paper): core papers, TeX sources, PDFs, and release metadata.
- [`extra/`](extra): the published compact proof plus repository-only focused research PDFs.
- [`code/`](code): certificates, simulations, particle calculations, and experiments.
- [`book/`](book): legacy book source and downloadable PDF, retained outside the primary reading path.
- [`cosmology/`](cosmology): dark-sector and cosmology research.
- [`physics-problems/`](physics-problems): focused applications and open-problem notes.
- [`docs/`](docs): claim policy, falsification program, and technical audit material.
- [`assets/`](assets): diagrams and public figures.

The simulation source is maintained in the companion
[oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim)
repository, which produces the simulation receipts and evidence artifacts
cited here.

## Explore OPH

- [The technical paper](flagship/from_observer_consensus_to_standard_physics.pdf)
- [Core paper index](paper/)
- [Textbooks](https://learn.floatingpragma.io)
- [Interactive simulation](https://simulation.floatingpragma.io)
- [OMEGA applications and hardware](https://omega.floatingpragma.io)
- [Blog](https://blog.floatingpragma.io/)
- OPH Sage on [Telegram](https://t.me/HoloObserverBot) and [X](https://x.com/OphSage)

## Contribute

OPH welcomes proofs, counterexamples, simulations, audits, and readable
explanations. The [reproduction guide](REPRODUCE.md) rebuilds the certificates
and checks from a clean clone. The
[scoped research questions](https://github.com/FloatingPragma/observer-patch-holography/issues)
identify suitable contributions, while the
[selection ledger](docs/SELECTION_LEDGER.md) states their exact theorem
premises and unresolved mathematical inputs.

## License

The repository uses split licensing. All software, including the Lean library, [`code/`](code), and [`tools/`](tools), is licensed under [Apache-2.0](code/LICENSE). Papers, the book, documentation, figures, and data are licensed under [CC BY-NC-SA 4.0](LICENSE). Hardware design files use CERN-OHL-W 2.0. The [LICENSE](LICENSE) file gives the per-directory map.
