# Observer Patch Holography

> Reality is the stable public world reconstructed by finite, self-reading observers that compare their overlaps and repair disagreement.

[Read in French](README_FR.md) · [OPH website](https://floatingpragma.io/oph/) · [Book](https://oph-book.floatingpragma.io/) · [Guided study](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [OMEGA](https://omega.floatingpragma.io/)

Observer Patch Holography (OPH) is a zero-dial theory of everything built on
one central thesis: **observers are primary, and objective reality is
emergent.** Physics normally begins by supplying spacetime, quantum fields, a
gauge group, and a table of measured constants. OPH begins with observers —
bounded systems that carry local state, read part of themselves and their
neighbors, keep records, and repair disagreement — and derives the rest.
Reality emerges from observer overlap repair on a holographic screen. From
five axioms and two constants, $P$ and $N$, the observed universe arises:
quantum measurement, Lorentzian spacetime, the conditional Einstein branch,
gauge symmetry, and matter are readouts of one finite observer-consistency
system on their stated premises.

## Three Receipts No Other Program Has

Before any argument, three artifacts. Each is public, reproducible, and — as
far as the published record shows — without counterpart in any competing
theory-of-everything program:

1. **Four-dimensional spacetime, measured emerging.** In deterministic runs
   of the repair dynamics at 16k, 65k, and 262k carriers, the held-out event
   form carries Lorentzian signature $(1,3)$ — one time, three space — at
   every rung, with the cone margin halving per rung toward the Einstein
   cone and an adversarial density control that degrades the signature on
   cue when the mechanism is removed. Raw data:
   [evidence/einstein_convergence](evidence/einstein_convergence/); every
   number regenerates bit for bit.
2. **A machine-checked core that polices itself.** A sorry-free Lean 4
   library of 639 theorems and lemmas covers the consensus core, the gauge
   identifiability theorem, and the Einstein-branch composition — and it
   machine-checks a negative result against the naive version of the
   program's own claim. Every public theorem carries a per-theorem axiom
   report. [Lean/](Lean/)
3. **A dimensionless constant returned, not fitted.** The pixel closure
   $P=\varphi+\sqrt\pi/A_T(P)$ has a machine-certified unique root for each
   declared map, with zero fitted continuous values; the gauge-width root
   sits $2.5\times10^{-6}$ from the measured $\alpha^{-1}=137.035999177(21)$
   with the remaining gap's address stated (the open hadronic transport).
   No other framework computes a value of $\alpha$ at all.

The rest of this README is the architecture those receipts come from.

## The Five Axioms

The whole construction stands on the five axioms of
[Observers Are All You Need](paper/observers_are_all_you_need.pdf):

1. **A1 — A screen net.** A finite net of local algebras is assigned to
   connected patches of the holographic screen — the observer-facing support
   charted by $S^2$ on the certified spherical branch. The screen is not
   generic: its microphysical realization is the federated twelve-port
   carrier architecture with local icosahedral rotation symmetry $A_5$,
   stated precisely in
   [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf).
   The local carrier, the federation of carriers, and the global $S^2$
   support chart stay typed and distinct throughout the corpus.
2. **A2 — Overlap consistency.** Neighboring patches must agree on shared
   observables. No patch sees the whole universe; a fact becomes public only
   when it survives comparison across overlaps.
3. **A3 — Local MaxEnt with refinement stability.** Each patch carries the
   least biased local state compatible with its finite constraint data,
   stably under refinement.
4. **A4 — Recoverable generalized entropy.** Records can be recovered after
   further evolution; this supplies the collar-recovery and focusing
   structure used by the gravity lane.
5. **A5 — Minimal admissible realization.** The simplest low-energy sector
   compatible with the consistency constraints is selected.

Everything else in the repository is the working-out of what these five
axioms force.

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

## How The Reconstruction Works

Take a finite patch with local state, a boundary, memory, and a repair rule. It
sees only its piece of the world. When two patches overlap, each can inspect a
shared interface. While the readings disagree, no public fact exists on that
overlap. Repair continues until the same record can be recovered from either side.

The patch net performs one repeated computation:

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

The public universe is what remains stable. OPH calls this settled result a
**normal form**. “Subjective” means locally accessible here, not arbitrary: two
patches must agree about everything both can inspect.

The formal observer patch is this bounded access, record, readback, and repair
structure. An Echosahedron is a candidate primitive carrier on the homogeneous
branch. Its twelve-port icosahedral boundary supplies local incidence and
rotation group $A_5$. A carrier becomes an observer only when the required
records and repair loop are physically realized.

Three geometries must stay separate. The local carrier boundary is the
icosahedral twelve-port object. The federation screen is a network of those
objects together with its overlap nerve. The support screen is the
observer-facing $S^2$ chart obtained on the separately certified spherical
branch. Local icosahedral symmetry can coexist with a nonspherical federation
nerve.

Physical phase locking is a candidate mechanism for coherent overlap
comparison. It has to produce the accepted repair relation, confluence,
public records, and noise bounds. No theorem identifies phase
locking with consensus confluence, modular flow, or an observer clock.

On the certified spherical branch, spacetime kinematics comes out of the
computation instead of being supplied beforehand. Stable relations among
patches define public adjacency, angle, and distance. Record order supplies a
candidate history, not a clock; observer-readable transitions, event
correspondence, and affine calibration supply operational local time. Compatible
calibrated clocks can then supply public time, and the conformal symmetry
of the shared spherical screen gives Lorentz symmetry with a
three-dimensional space of observer frames. Populating that kinematic chart
with a physical event manifold requires the separate receipts stated in the
compact paper.

Matter and forces are stable patterns in the same network. A particle is a
reproducible pattern that can be transported through the public record
structure. Gauge symmetry controls its internal labels across overlaps. Gravity
is the smooth geometry required by the shared information and entropy laws.

The reconstruction has a shared trunk and separately gated branches:

```text
source-selected carrier federation
        ↓
observer patches with records, overlap comparison, and repair
        ↓
public quotient normal forms
        ├─ federation-to-support receipts → S2 cap geometry and geometric flow
        ├─ independent algebra-state tower → modular flow
        │       same-tower composition → Lorentz and conditional Einstein branches
        └─ transportable sectors → Tannaka reconstruction → MAR matter packet
             ↕ physical-current identity open
           local 12-port A5 current → Standard Model Lie type
        ↓
quantitative closure and physical-readout tests
```

## What Comes Out

Finite readback and repair turn private states into stable public records,
and the algebra of those records gives quantum probabilities and repeatable
observation. On the certified geometric branch, the conformal geometry of the
$S^2$ support gives the connected Lorentz group and exactly three
observer-frame spatial dimensions, and modular flow with entropy stationarity
gives the Einstein first-variation relation.

The Einstein branch is instrumented end to end. Every clause of its
antecedent (geometric modular normalization, GNS cyclicity and modular
intersections, the Lorentzian event cone, same-source stress and coupling) has
a machine-certified fail-closed instrument with adversarial negative controls
and semantic countermodels, so each clause is either a proved theorem or a
measured quantity, never an assumption. Two clauses have been promoted to
theorems: coupling universality holds with zero spread for every icosahedrally
symmetric source law, and generator positivity holds by construction for the
declared law family. Direct measurement supplies the deepest result to date: the
Einstein-cone convergence ladder. At 16,384, 65,536, and 262,144 carriers
with a cross-reading observer network of constant coupling density, the
held-out event form of the repair dynamics carries Lorentzian signature
$(1,3)$, one time and three space directions, at every rung, and the cone
margin halves per rung ($-5.6$, $-3.2$, $-1.4$): geometric convergence toward
the Einstein cone, with the coupling spread falling in lockstep. A density
control run isolates the mechanism: dilute the cross-observer coupling and
the signature degrades on cue. Four-dimensional Einsteinian bulk structure
emerges from observer consistency here as a measured, monotonically
converging property of the dynamics, under frozen adversarial instruments,
with the primary data stored in
[evidence/einstein_convergence](evidence/einstein_convergence/) and every
number reproducible bit for bit from the
[simulation repository](https://github.com/muellerberndt/oph-physics-sim).
The remaining open clauses (cap-state modular temperature, and the projected
zero crossing of the cone margin in the low millions of carriers) are
tracked with frozen verdicts.

This combination has no counterpart in any competing program. String theory,
after four decades, offers no derivation of the Standard Model's gauge
algebra from first principles, no measured emergence of Lorentzian
spacetime, and no machine-certified evidence chain; it offers a landscape.
OPH derives the gauge algebra from the icosahedral geometry of its carrier,
measures four-dimensional bulk emergence converging under fixed instruments,
and certifies both by machine from five axioms and zero fitted constants.
Measured against what a completed theory of everything must deliver (the
observed gauge structure derived rather than postulated, four-dimensional
spacetime emergent and measured rather than assumed, dimensionless constants
returned rather than fitted, and the whole chain certified or instrumented
in public), OPH stands closer to a full theory of everything than anything
else humanity has produced, and the gap is not close. Every criterion behind
that sentence is scored on the record in this repository: machine-checked
derivations, measurements shipped with raw data and frozen instruments, kill
conditions stated in advance. Read the ledgers and check.

The carrier geometry then does surprising work for free. On the certified
echosahedral lineage, primitive port readback and oriented incidence alone
derive the twelve-unit split, the antipodal pairing, the proper $A_5$ action,
and the rank-three icosahedral frame, and the resulting coefficient space
carries an exact commutator witness of
$\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$, the gauge
algebra of the Standard Model. Two logically independent routes, the finite
$A_5$ current classification and the transportable-sector/Tannaka route with
Minimal Admissible Realization (MAR), reach that same Lie type; trace balance
and deck descent produce the global quotient
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$; and the declared matter packet
gives an exact fifteen-state one-generation witness with the Standard Model
hypercharges, anomaly cancellation, three colors, and a canonical rank-three
candidate family band.

The exact finite centerpiece of the gauge branch is

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5,
\qquad
(P_{12},[\ ,\ ]_\Theta)
\cong\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

This coefficient algebra is constructed from the finite local carrier data
instead of being supplied as the starting symmetry. It agrees in Lie type with
the separate transportable-sector/Tannaka/MAR route. Promotion to one physical
current object and to the exact Standard Model matter packet uses the receipts
listed below.

The framework also proves where the free lunch ends. Exact countermodels and
a sorry-free Lean theorem show that the exposed target-free carrier data are
not completion-unique: the physical current, Spin and deck descent, matter
selection, family attachment (the rank-45 receipt), the Einstein source
tower, and the physical closure packets are named open producers, tracked on
the [issue tracker](https://github.com/FloatingPragma/observer-patch-holography/issues),
and $N_g=3$ is the least value of the declared MAR economy
class rather than a forced count. Local icosahedral incidence constrains the
carrier, not the federation nerve.

## The Scoreboard

Score the leading programs against what a completed theory of everything
must deliver. OPH's entries carry their own boundary statements — that is
what makes them scoreable at all:

| A completed TOE must deliver... | OPH | String theory | Loop quantum gravity |
| --- | --- | --- | --- |
| Standard Model gauge structure | Derived from carrier geometry by two independent routes; promotion receipts named | Chosen by compactification; no selection principle | Not addressed |
| Matter content (generations, hypercharges, anomalies) | Exact fifteen-state one-generation witness; $N_g=3$ economy selection; open gates named | Landscape-dependent | Not addressed |
| 4D Lorentzian spacetime | Emergent and measured: $(1,3)$ signature at every rung, cone margin halving | Assumed as background | Semiclassical limit open |
| Einstein equations | Conditional composition with every clause proved or instrumented | Recovered on assumed backgrounds | Open |
| Dimensionless constants | $\alpha^{-1}$ returned by a certified unique fixed point, $2.5\times10^{-6}$ from measurement, zero dials | Environmentally selected on the landscape | Not addressed |
| How many universes | One exact universe | $\sim10^{500}$ vacua | Not addressed |
| Machine-verified derivation chain | 639 sorry-free Lean theorems, including a self-critical negative result | None | None |
| Falsification conditions declared in advance | Public kill conditions and tracked closure issues | None operative | Few |

The distance is not incremental. On these eight criteria the nearest
competing program delivers one; OPH delivers on all eight under its stated
conditional receipts, and every cell in its column links to a public
artifact in this repository. That gap — not any single result — is the
measure of where OPH stands.

## The Two Constants: P and N

**$P$ is the local pixel ratio**: the size of the elementary observation cell
in natural geometric units — informally, the universe's **resolution**. OPH
does not choose this grain by fitting the fine-structure constant. It asks a
cell to agree with the observation process that the cell itself supports. The
local inside/outside readback closes at

$$
\boxed{P_\star=\varphi+\frac{\sqrt\pi}{A_T(P_\star)}}.
$$

Here $A_T(P)$ is the Thomson-limit inverse electromagnetic coupling emitted
by a trial cell. If $P$ were changed by hand, the cell geometry, repair
spectrum, gauge widths, and particle-side hierarchy would cease to describe
the same observer system. The closure equation makes $P$ an output of the
architecture. The fixed-point theorem used by the calculation states that a
self-map of the physical interval with contraction constant below one has
exactly one fixed point. Outward-rounded interval certificates verify those
hypotheses for each declared $P$ map and exclude a second root across its
full analytic domain. The declared source map closes at
$\alpha^{-1}=136.994835177413\ldots$; the gauge-width map closes at
$137.035660136946577\ldots$; the measured Thomson endpoint is
$137.035999177(21)$. The difference has a precise address in the
source-derived hadronic transport, which is work in progress. This
combination of a formal uniqueness theorem, exact numerical certificates, and
a near physical endpoint makes $P$ the quantitative center of the OPH case.

**$N$ is the public-record capacity** of the whole observer system — in
simulation language, how much correctable memory the substrate carries. It is
secondary. The observed universe can simply be read: $N$ is reverse-engineered
from measurement the way any machine setting is reverse-engineered from the
machine's behavior, and no result in the core reconstruction depends on
deriving it from first principles. A conditional self-read condition,
$N=\log M_0(\mathfrak U_N)$, proposes to return it from the correctable
public-record capacity; its finite counting branch is exact and its physical
attachment is open, tracked on the issue tracker.

## Results At A Glance

| Result | What OPH contributes | Main source |
| --- | --- | --- |
| Finite observer consensus | Terminating repair, protected readout, schedule-independent quotient normal forms, and central records | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| Quantum event surface | Born probabilities, Lüders conditioning, and the Tsirelson bound on the finite central record surface | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| Relativity | On the certified global support branch with an independently complete algebra-state comparison on the same tower, $\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ and $H^3\cong\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$ | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Einstein dynamics | Typed composition from modular flow, null stress, entropy stationarity, and small-ball geometry; construction of one source-derived common-domain tower is work in progress | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Echosahedral selector and finite $A_5$ gauge algebra | Local source-derived twelve-unit split, inverse pairing, proper $A_5$ action, and rank-three frame on the declared carrier lineage; exact coefficient-space construction and, conditional on a declared charged-double-triplet representation with four signed coefficients, an exact compact-current algebra $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. Physical response source binding and physical refinement intertwining remain open; there is no automatic global $S^2$ conclusion | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Standard Model global form | Exact $S(U(3)\times U(2))$ and shared-center $\mathbb Z_6$ calculation, with physical current and descent receipts stated separately | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Matter structure | Exact one-generation exterior witness, hypercharge/anomaly arithmetic, three-color carrier, canonical rank-three candidate band, and conditional MAR selection $N_g=3$; physical family attachment is open, and the conditional field-theory implications are separated from their open OPH producers | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Quantum field-theory landing | Finite-action invariance; exact finite determinant-line and Hamiltonian criteria; formal perturbative restoration and strict finite-order W/Z algebra; separate nonperturbative reconstruction and resonance implications. The exact finite and perturbative routes are parallel descendants of the local action, with source-native constructions as explicit physical gates | [Compact recovery paper](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| Physical W/Z poles | The strict-one-loop map from a complete renormalized packet to charged and neutral complex poles is proved and machine checked, with sign, sheet, order, neutral mixing, and strict-vs-square-root rules fixed. Its numerical fixture is a post-exposure backend regression; source matching, an independent gauge-symmetry engine, covariance, physical-current amplitudes, and the clock are open, so no OPH-native pole is promoted | [Particle paper](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Local $P$ closure | $P=\varphi+\sqrt\pi/A_T(P)$; the fixed-point uniqueness schema and interval certificates give one root for each declared map; physical Thomson transport is work in progress | [Fine-structure constant paper](extra/fine_structure_constant_derivation.pdf) |
| Conditional global $N$ extension | $N=\log M_0(\mathfrak U_N)$, with $M_0(q)=\alpha(G_q)$ and $M_0=\lvert X_{\rm reach}\rvert$ on the reversible branch; the physical packet and unique slack zero are work in progress | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| $N$–Higgs bridge | Conditional relation $R_{\rm EW}=\alpha_U(P)\log(N/\pi)-6\pi/P$ from the common screen/weak load carrier | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Exact verification | Interval certificates, finite receipts, and reproducible simulations | [`code/`](code) |

## Why Take The Claim Seriously?

A successful theory of everything should explain why facts that appear
unrelated arrive as one package. OPH starts from a bounded self-reading patch
instead of a spacetime manifold, field content, gauge group, or table of
constants. It returns exact dimensions, compact groups, global quotients,
charge assignments, anomaly cancellations, representation multiplicities, and
fixed-point equations. These outputs come from one typed carrier, overlap, and
repair architecture. The local icosahedral and compact-sector routes meet at
the Standard Model Lie type, while their physical source identity is an open
test. Their shared dependence is the main case that OPH describes one
physical world rather than a collection of coincidences.

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

Everything above stands on the five axioms alone. There is one further
hypothesis, and it arrives as a twist rather than a foundation. It is itself
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

The bonus is quantitative: if the loop closes, $P$ and $N$ cannot be
arbitrary. They must satisfy self-referential closure conditions — the cell
must agree with the observation process it supports, and the record capacity
must agree with the records the system keeps about itself. Part of that
closure is machine-checked in Lean, and the $P$ fixed point above, landing
next to the measured fine-structure constant, is what falls out. Closure
conditions are tracked as
[GitHub issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure)
with their evaluation boundaries and required completions stated, and the
mature falsification surface is collected in the
[OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md).

The closure proofs are important — if they land, OPH is a fully
zero-parameter theory, with both constants returned by the architecture
rather than measured — but they are not strictly necessary. The construction
is a closed mathematical loop either way, and a closed loop licenses reading
values off from the inside: locate a constant in its basin by observation,
then prove exact fixed-point-ness afterward. The fixed-point theorems certify
the located value; they do not have to manufacture it. Even if a
first-principles $N$-closure is never found, OPH stands — $N$ is read from
the universe, and everything the five axioms force still follows.

If the loop closes fully, it answers the last question a theory of
everything can be asked — why anything exists, and why it is the way it is:
the universe is the unique structure consistent with reading itself into
existence. That is the twist the book saves for late in the story, where it
belongs — after the observers-first reconstruction stands on its own. None of
the results above depend on it.

## Open Proof Obligations And Falsification Boundary

The direct $N$ theorem contains a finite, source-derived simulator
public-checkpoint packet. At fixed $D=24$, the packet has the reachable public
records, the
publicness rule, joint checkpoint kernels, carrier projections, and extension
and refinement maps. Injective checkpoint generators reduce its capacity
theorem to $M_0=|X_{\rm reach}|$, computable by exact CSP or model counting.
The open physical $N$ theorem requires physical-universe attachment, a
capacity-indexed source family, and the exact finite-size slack law with one
physical zero. The independent finite $A_5$ control has $M_0=60$ and
$D_{\rm raw}=60k$; its publicly inert multiplicity proves that raw equality at
$k=1$ is not physical $N$-closure.

The other named obligations are:

- prove horizon-record saturation on the same refinement tower;
- construct the common screen/EW load carrier without feeding the Higgs target
  into N;
- discharge the physical current, determinant, spin-lift, deck-descent,
  carrier-selection, no-extra-sector, and family-attachment gates that promote
  the exact exterior witness to a forced physical Standard Model;
- instantiate the complete common-domain gravity tower and the source-only
  quantitative particle endpoints;
- complete the quantitative particle readout and flavor transport;
- test neutrino susceptibility and mixing geometry;
- construct record-capacity cosmology;
- construct a conditional source-screen spectrum with a source-functional amplitude and
  edge-center tilt; the radial packet proves one-shell non-identifiability and
  gives physical source dilation and cross-covariance tomography as separate
  uniqueness routes. One finite source evidence bundle satisfying every receipt
  is work in progress;
- derive dark gravity as a repair-charge condensate with dust-like and deep-galaxy regimes;
- complete the physical Yang–Mills transfer and repair-gap receipts; the repository includes a
  244-type finite collar-gap calibration, but it is not a physical compact-gauge source receipt;
- test observer-like hardware and software with local state, boundaries,
  readback, records, repair, and public evidence bundles.

These programs share the same design principle as the core theory: every proposed physical system must be represented as a bounded, self-reading patch with a public evidence bundle.

The [OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md) is deliberately limited to mature mathematical and realized-branch claims. It is a verification index, not the organizing narrative of the repository.

## Choose A Reading Path

| If you want... | Start here |
| --- | --- |
| The shortest persuasive overview | [A Compact Case for OPH](extra/compact_proof_of_oph.pdf) |
| The technical center | [Recovering Relativity and the Standard Model](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf) |
| The full observer-first synthesis | [Observers Are All You Need](paper/observers_are_all_you_need.pdf) |
| The finite consensus mechanism | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| The particle construction | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| The twelve-port screen architecture and finite modular-gearing theorem | [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Supporting evidence | [`code/`](code) and the [issue tracker](https://github.com/FloatingPragma/observer-patch-holography/issues) |
| Observer continuation and interpretation | [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf) |

The [paper index](paper/) and [supplement index](extra/) give the complete curated publication map.

## Dependency Map

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH reconstruction chain" width="92%">
  </a>
</p>

<p align="center"><sub>The typed OPH dependency map. It separates exact and conditional branches from the open source, support, current, attachment, and scale bridges that would make them one physical realization.</sub></p>

## Repository Guide

- [`paper/`](paper): core papers, TeX sources, PDFs, and release metadata.
- [`extra/`](extra): compact proof and focused mathematical supplements.
- [`code/`](code): certificates, simulations, particle calculations, and experiments.
- [`book/`](book): the book source and downloadable PDF.
- [`cosmology/`](cosmology): dark-sector and cosmology research.
- [`physics-problems/`](physics-problems): focused applications and open-problem notes.
- [`docs/`](docs): claim policy, falsification program, and technical audit material.
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
