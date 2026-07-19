# Observer Patch Holography

> Reality is the stable public world reconstructed by finite, self-reading observers that compare their overlaps and repair disagreement.

[Read in French](README_FR.md) · [OPH website](https://floatingpragma.io/oph/) · [Book](https://oph-book.floatingpragma.io/) · [Guided study](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

Observer Patch Holography (OPH) is an observer-first reconstruction program for fundamental physics. It begins with finite patches carrying local state, boundaries, records, readback, and repair moves. A public fact is one that survives when overlapping patches compare what they can see and settle into a common normal form.

The organizing equation is simple:

$$
T(\mathfrak U_{\mathrm{OPH}})=\mathfrak U_{\mathrm{OPH}}.
$$

The universe is modeled as a fixed point of its own observer-accessible readback and repair process.

## The Idea In Plain Language

Physics normally begins with a universe already equipped with spacetime,
quantum fields, a gauge group, and a list of measured constants. OPH asks a
more radical question: **what is the smallest kind of system capable of having
a world at all?**

The answer is an observer patch. It is not necessarily a person. It is any
bounded physical or computational system that has a local state, a boundary,
memory, the ability to read part of itself and its neighbors, and a way to
repair disagreement. No patch sees the whole universe. A fact becomes
objective only when it can be written, compared across overlaps, recovered
after further evolution, and retained as part of the public record.

OPH proposes that this is not merely how inhabitants learn about a pre-existing
world. It is the mechanism by which a public physical world is selected. A
candidate universe must support observers that can read it, and their readback
must reproduce the same universe without contradiction. That is the meaning of

$$
T(\mathfrak U)=\mathfrak U.
$$

The theory therefore has no external ruler, master clock, preferred observer,
or list of adjustable physical constants. Its numbers must be returned by the
same self-reading consistency loop that returns its laws.

## Two Closures, No Dials

OPH has two quantitative self-readback problems. They answer complementary
questions.

### P: how fine is one local act of observation?

$P$ is the local pixel ratio: the size of the elementary observation cell in
natural geometric units. OPH does not choose this grain by fitting the
fine-structure constant. It asks a cell to agree with the observation process
that the cell itself supports. The local inside/outside readback closes at

$$
\boxed{P_\star=\varphi+\frac{\sqrt\pi}{A_T(P_\star)}}.
$$

Informally, $P$ is the universe's **resolution**. If it were changed by hand,
the cell geometry, repair spectrum, gauge widths, and particle-side hierarchy
would no longer describe the same observer system. The closure equation makes
$P$ an output of the architecture rather than a dial.

### N: how much of the universe can its observers read back?

$N$ is the global correctable public-record capacity. Its official equation is

$$
\boxed{N=M_0(\mathfrak U_N)}.
$$

Read it literally: construct a trial universe with capacity $N$, put its
observers inside it, and count how much information remains publicly readable
after every allowed observation and continuation. Closure occurs when the
inside readback equals the outside supplied capacity.

At finite size, compatible observer records form a graph in which two records
are connected when some allowed checkpoint can confuse them. The largest set
that never becomes confused is the exact readback code,

$$
M_0(q)=\alpha(G_q).
$$

On the exact reversible branch, nothing is erased and every checkpoint is
injective, so this becomes the especially transparent count

$$
M_0(q)=|X_{\rm reach}(q)|.
$$

In other words, $P$ says how finely the universe resolves itself; $N$ says how
much stable public world it can hold. Neither is selected from a data table.
Both are demanded by self-readback.

The Higgs then supplies a remarkable cross-check rather than an input. The same
load that appears globally on the observer screen appears locally in the
four-copy weak sector. Identifying those two independently constructed carriers
gives

$$
\alpha_U(P)\log\!\left(\frac{N}{\pi}\right)=\frac{6\pi}{P},
$$

while horizon saturation gives

$$
\Lambda\ell_\star^2=\frac{3\pi}{N}.
$$

One pair $(P,N)$ therefore links the local observation grain, the Higgs/weak
hierarchy, and the cosmological horizon.

## A Complete Universe From Consistency

This is the main OPH claim: **one observer-consistency architecture, with no
fitted continuous inputs, reconstructs the architecture of our observed
universe.** The logic is cumulative and overdetermined:

- finite readback and repair turn private states into stable public records;
- the algebra of those records gives quantum probabilities and repeatable
  observation;
- the conformal geometry of a shared spherical screen gives the connected
  Lorentz group and three observer-frame spatial dimensions;
- modular flow, recoverability, null transport, and entropy stationarity give
  the Einstein branch;
- the maximally symmetric twelve-port icosahedral screen decomposes into the
  exact dimensions needed for
  $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$;
- trace balance and global consistency produce
  $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ rather than merely its Lie
  algebra;
- the minimal chiral one-Higgs realization gives the Standard Model
  hypercharges, anomaly cancellation, three colors, three generations, and one
  Higgs doublet;
- the same $P$ and $N$ closures connect the gauge, Higgs, gravitational, and
  cosmological scales.

These are not unrelated curve fits. They are different consequences of the
same demand: every local description must glue into one recoverable world that
can contain observers capable of reading it back. Any individual match could
be accidental. The joint recovery of quantum observation, Lorentzian
spacetime, Einstein gravity, the exact Standard Model global form and matter
counts, and the Higgs/cosmological hierarchy from zero fitted theory values is
why OPH presents itself as a full **theory of everything with no dials**.

The detailed papers state every theorem-local premise. The proof obligations
that turn the remaining conditional bridges into unconditional theorems are
collected near the end of this README, where they belong as a finite completion
program rather than as the introduction to the idea.

## OPH For Absolute Beginners

Imagine a simulator that has no finished universe map to consult. There is no pre-drawn four-dimensional grid and no master clock outside the system. The simulator is made from many small **observer patches**. Each patch has only:

- a local state: what is happening from this patch’s perspective;
- ports or boundaries: what it can exchange with neighboring patches;
- memory: records of earlier reads;
- self-readback: the ability to inspect part of its own state;
- repair moves: ways to respond when overlapping records disagree.

This is OPH’s operational model of a subjective perspective. Every patch sees only its own finite piece of reality. “Subjective” means locally accessible, not arbitrary. Two patches that overlap must agree about the information both can inspect.

The patch-net repeatedly performs a simple computation:

```text
read local state
      ↓
exchange boundary records
      ↓
compare overlapping descriptions
      ↓
repair disagreement
      ↓
write the stable result and repeat
```

The public universe is the part that remains stable after this process. OPH calls that stable result **consensus** or a **normal form**.

The icosahedron supplies the first highly symmetric finite scaffold for this computation. It has twelve vertices, twenty triangular faces, and rotation group $A_5$. In the OPH screen model, the twelve locations become observer-facing ports, the faces organize local overlaps, and repeated refinement turns the finite scaffold into a progressively finer patch-net. The “echosahedral” language used throughout the project refers to this icosahedral observer-screen architecture.

On the certified spherical branch, spacetime kinematics comes out of the
computation instead of being supplied beforehand. Stable relations among
patches define public adjacency, angle, and distance. Record order supplies a
local clock, compatible clocks supply public time, and the conformal symmetry
of the shared spherical screen gives Lorentz symmetry with a
three-dimensional space of observer frames. Populating that kinematic chart
with a physical event manifold requires the separate receipts stated in the
compact paper.

Matter and forces are stable patterns inside the same network. A particle is a reproducible, transportable pattern in the public record structure. A gauge force describes how internal labels must transform when those patterns move between overlapping patches. Gravity describes the smooth large-scale geometry required by the shared information and entropy relations.

This is the **consciousness-first** or **observer-first** thesis in its simplest form. Conventional physics usually starts with matter fields in a four-dimensional spacetime and asks how observers or consciousness arise inside it. OPH starts with the minimal structure required for any perspective to exist at all: bounded access, self-readback, memory, comparison, and stable continuation. The shared physical world is reconstructed from the agreement conditions among those perspectives.

“Observer” is a precise structural role here. A human mind, an organism, a physical instrument, or a software process can instantiate it when it has the required local state, boundary, records, readback, and repair loop. OPH does not say that human thoughts manufacture reality. It says that a reality with no possible local perspective, record, or self-consistent readback is not yet the public physical world described by observers.

In simulation language, the universe is its own simulator. The patches, computation, records, and resulting world all belong to the same closed system. No external computer or programmer is required by the formal construction.

## What OPH Derives From Observer Consistency

OPH uses one mathematical architecture across subjects that are normally introduced separately:

- finite consensus gives stable public records and quotient normal forms;
- central record algebras give quantum event probabilities and conditional update;
- the conformal geometry of an observer screen gives the connected Lorentz group and a three-dimensional observer-frame space;
- modular flow, null transport, entropy stationarity, and small-ball geometry compose into the Einstein relation on a source-derived common-domain tower;
- transportable charges and compact reconstruction give the Standard Model gauge structure;
- a finite twelve-port $A_5$ construction produces the Lie algebra
  $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$;
- trace balance integrates this algebra to
  $S(U(3)\times U(2))\cong(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$;
- the exterior/MAR branch selects the Standard Model charge lattice, three colors, three generations, and one Higgs doublet;
- correctable public-record closure supplies the global capacity equation and its independent horizon and weak/Higgs bridges;
- public Lean, exact-arithmetic certificates, simulations, and executable receipts check the finite mathematical core.

Measurement, spacetime, gravity, and gauge structure are tested against the
same mechanism: finite observers form public records by comparing overlaps and
repairing disagreement. The finite theorems reach from quotient consensus to
the Standard Model gauge algebra and global quotient. The continuum theorem
stack reaches the Lorentz and Einstein branches under its stated geometric,
modular, stress, entropy, and scale premises. This reuse of one mechanism is
the programme's central result.

## The Strongest Finite Result

On the declared twelve-port icosahedral branch, the permutation module splits as

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5.
$$

An explicit equivariant pullback of the block commutator then constructs

$$
(P_{12},[\ ,\ ]_\Theta)
\cong
\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

This is the local Lie algebra of the Standard Model gauge forces. It is obtained from the finite coefficient geometry rather than inserted as the starting symmetry.

The same construction gives two independent appearances of the number $24$:

$$
m_{\mathrm{rep}}=2(8+3+1)=24,
$$

while the twelve screen ports have $24$ oriented slots. One count comes from the recovered gauge algebra; the other comes from the screen geometry.

## One Reconstruction Chain

```text
self-reading patches
        ↓
records, overlap comparison, repair
        ↓
public quotient normal forms
        ↓
screen conformal geometry and modular flow
        ↓
Lorentz kinematics, observer time, Einstein dynamics
        ↓
transportable charges and compact reconstruction
        ↓
SU(3) × SU(2) × U(1) / Z6 and Standard Model matter
        ↓
quantitative fixed-point and physical-readout programs
```

The detailed hypotheses and receipt types are stated in the papers. This front
page intentionally leads with the positive reconstruction; the consolidated
proof obligations and falsification boundary appear near the end.

## Results At A Glance

| Result | What OPH contributes | Main source |
| --- | --- | --- |
| Finite observer consensus | Terminating repair, protected readout, schedule-independent quotient normal forms, and central records | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| Quantum event surface | Born probabilities, Lüders conditioning, and the Tsirelson bound on the finite central record surface | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Relativity | $\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ and $H^3\cong\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$ | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Einstein dynamics | The typed modular/null, stress, entropy, and small-ball chain to $G_{ab}+\Lambda g_{ab}=8\pi G\langle T_{ab}\rangle$ on its declared common-domain tower | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Finite $A_5$ gauge algebra | Exact twelve-port construction of $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$ | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Standard Model global form | $S(U(3)\times U(2))$ and the shared-center $\mathbb Z_6$ quotient | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Matter structure | Hypercharge lattice, three colors, three generations, and one Higgs doublet on the realized Minimal Admissible Realization (MAR) branch | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Global N closure | Official equation $N=M_0(\mathfrak U_N)$; finite implementation $M_0(q)=\alpha(G_q)$, whole-fiber robust saturation, exact reversible count $M_0=|X_{\rm reach}|$, and the unique finite-size slack target | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| N–Higgs bridge | $R_{\rm EW}=\alpha_U(P)\log(N/\pi)-6\pi/P$ from the common screen/weak load carrier | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Exact verification | Lean theorem subset, interval certificates, finite receipts, and reproducible simulations | [`Lean/`](Lean) and [`code/`](code) |

## Choose A Reading Path

| If you want... | Start here |
| --- | --- |
| The shortest persuasive overview | [A Compact Case for OPH](extra/compact_proof_of_oph.pdf) |
| The technical center | [Recovering Relativity and the Standard Model](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| The full observer-first synthesis | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| The finite consensus mechanism | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| The particle construction | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| The twelve-port screen architecture | [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Machine-checkable evidence | [`Lean/`](Lean), [`code/`](code), and the [closure ledger](docs/CLOSURE_LEDGER.md) |
| Observer continuation and interpretation | [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf) |

The [paper index](paper/) and [supplement index](extra/) give the complete curated publication map.

## Evidence You Can Inspect

The repository contains several complementary forms of evidence:

- hand proofs in the TeX papers;
- a sorry-free 111-theorem Lean subset covering finite observer consensus and related algebraic results;
- interval and uniqueness certificates for declared numerical maps;
- finite carrier and hierarchy receipts;
- particle, geometry, dark-sector, and quantum-hardware code;
- simulations reaching more than one million patches;
- a claim registry and closure ledger connecting prose claims to artifacts.

## Remaining Proof Obligations And Falsification Boundary

The positive case above is the organizing narrative. The remaining work is
finite and named rather than hidden in prose. The highest-leverage object is
the source-only packet

```text
PUBLIC_CHECKPOINT_PACKET(r,D)
  = (X_reach, publicness policy, global joint kernels,
     carrier projections, refinement maps).
```

Build it first on the exact reversible branch. Injective checkpoint generators
reduce the capacity theorem to $M_0=|X_{\rm reach}|$, computable by exact CSP or
model counting. The remaining hard N theorem is then the exact finite-size
slack law with one physical zero.

The other named obligations are:

- prove horizon-record saturation on the same refinement tower;
- construct the common screen/EW load carrier without feeding the Higgs target
  into N;
- discharge the physical current, determinant, spin-lift, deck-descent,
  carrier-selection, no-extra-sector, and family-attachment gates that promote
  the exact exterior witness to a forced physical Standard Model;
- instantiate the complete common-domain gravity tower and the source-only
  quantitative particle endpoints;
- continue the active phenomenology programs:

- quantitative particle readout and flavor transport;
- neutrino susceptibility and mixing geometry;
- record-capacity cosmology;
- a conditional source-screen spectrum with a source-functional amplitude and
  edge-center tilt; the radial packet proves one-shell non-identifiability and
  gives physical source dilation and cross-covariance tomography as separate
  uniqueness routes. One finite source evidence bundle satisfying every receipt
  is work in progress;
- dark gravity as a repair-charge condensate with dust-like and deep-galaxy regimes;
- finite Yang–Mills transfer and repair-gap constructions;
- observer-like hardware and software systems with local state, boundaries, readback, records, and repair.

These programs share the same design principle as the core theory: every proposed physical system must be represented as a bounded, self-reading patch with a public evidence bundle.

The [OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md) is deliberately limited to mature mathematical and realized-branch claims. It is a verification index, not the organizing narrative of the repository.

## Dependency Map

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH reconstruction chain" width="92%">
  </a>
</p>

<p align="center"><sub>The OPH line from observer consistency to public records, spacetime, gravity, gauge structure, matter, and quantitative readout.</sub></p>

## Repository Guide

- [`paper/`](paper): core papers, TeX sources, PDFs, and release metadata.
- [`extra/`](extra): compact proof and focused mathematical supplements.
- [`Lean/`](Lean): machine-checked theorem development.
- [`code/`](code): certificates, simulations, particle calculations, and experiments.
- [`book/`](book): the book source and downloadable PDF.
- [`cosmology/`](cosmology): dark-sector and cosmology research.
- [`physics-problems/`](physics-problems): focused applications and open-problem notes.
- [`docs/`](docs): closure ledger, claim policy, and technical audit material.
- [`assets/`](assets): diagrams and public figures.

## Explore OPH

- [Theory explainer](https://floatingpragma.io/oph/theory-of-everything)
- [Interactive simulation](https://simulation.floatingpragma.io)
- [OMEGA applications and hardware](https://omega.floatingpragma.io)
- [Book](https://oph-book.floatingpragma.io)
- [Guided study](https://learn.floatingpragma.io)
- [Blog](https://blog.floatingpragma.io/)
- OPH Sage on [Telegram](https://t.me/HoloObserverBot) and [X](https://x.com/OphSage)

## License

The authored material is licensed under [CC BY-NC-SA 4.0](LICENSE). The repository-wide [OPH Open Use and Anti-Patent Covenant](PATENTS.md) keeps OPH-derived ideas, software, methods, devices, simulations, and hardware open to study, implementation, modification, and sharing without private patent monopolies.
