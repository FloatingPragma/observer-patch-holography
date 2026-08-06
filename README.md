# Observer Patch Holography

> Reality is the stable public world reconstructed by finite, self-reading observers that compare their overlaps and repair disagreement.

[Read in French](README_FR.md) · [Paper](flagship/from_observer_consensus_to_standard_physics.pdf) · [Book](https://oph-book.floatingpragma.io/) · [Textbooks](https://learn.floatingpragma.io/) · [Simulation](https://simulation.floatingpragma.io/) · [Hardware](https://omega.floatingpragma.io/)

Observer Patch Holography (OPH) is a zero-dial theory-of-everything research
program built on one central thesis: **observers are primary, and objective
reality is emergent.** Physics normally begins by supplying spacetime, quantum
fields, a gauge group, and a table of measured constants. OPH begins with
observers: bounded systems that carry local state, read part of themselves and
their neighbors, keep records, and repair disagreement. Reality emerges from
observer overlap repair on a holographic screen. From this architecture OPH
reconstructs an exact finite structural core: conditional quantum-record
identities, a conditional finite four-law package, a three-dimensional
observer-frame carrier, and explicit order/clock interfaces that do not
instantiate a physical observer-local time. It also derives Lorentz kinematics
on the stated global-support branch,
the Standard Model gauge Lie type, and a conditional one-generation exterior
matter pair.

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

- **The book.** [*Reverse Engineering Reality*](https://oph-book.floatingpragma.io/),
  also available as a [print-quality PDF](https://cfxrbtseaimxxqsxlrku.supabase.co/storage/v1/object/public/books/reverse-engineering-reality.pdf),
  tells the whole story: what the theory says, how it was discovered, and why
  the observer-first turn is the one physics has been circling for a century.
  It is written to entertain and it keeps the science exact.
- **The flagship paper.** [*From Observer Consensus to Standard Physics*](flagship/from_observer_consensus_to_standard_physics.pdf)
  gives the primary technical account of the observer-first reconstruction.
- **The textbooks.** The [OPH textbooks](https://learn.floatingpragma.io/)
  teach the theory the long way. Every basic derivation is worked in full,
  with the required math built up as you go. Volumes cover gravity, the
  Standard Model, and unification, each readable online or as a PDF.
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

## Eight Reproducible Physics Receipts

These public results link directly to their papers, proofs, data, and
certificates:

1. **Three-dimensional space emerges from the algebra of repair records.**
   The declared twelve-port response contains an exact abstract
   three-dimensional Euclidean completion. Adding comparison records and
   completing their distance gives ordinary continuous three-space without an
   assumed coordinate grid. Physical identification of these points and their
   scale is work in progress. See the
   [spacetime and Einstein paper](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf),
   the Lean [completion proof](Lean/Screen/PrimitivePortFrameQuotient.lean),
   and the [finite signature evidence](evidence/einstein_convergence/).
2. **Observer agreement reproduces the quantum rules.** On the declared finite
   model, agreement gives the usual quantum probabilities and measurement
   rules, together with exact limits on correlations and copying. Connecting
   this mathematics to physical measuring devices is work in progress. See the
   [flagship paper](flagship/from_observer_consensus_to_standard_physics.pdf),
   the [public-record proof](Lean/EventAlgebra/PublicRecordAlgebra.lean), and
   the [Born-rule boundary](Lean/EventAlgebra/FiniteEffectClosureBoundary.lean).
3. **Dynamics, quantum weights, and the action are forced, not chosen.** On
   the finite observer algebra, continuous time evolution has no freedom
   beyond a Hamiltonian, so the Schrödinger form is a theorem. The quantum
   weight is the only additive assignment on measurement effects, down to
   the qubit case. The realized dynamics fixes its own action up to gauge:
   least action and most probable history are two readouts of one
   functional, joined to Hamiltonian flow by a Legendre bridge. Attaching
   units, clocks, and outcome frequencies is work in progress. See the
   [observers paper](paper/observers_are_all_you_need.pdf), the
   [Born representation proof](Lean/EventAlgebra/FiniteBuschGleason.lean), and
   the [derived-action proof](Lean/InformationProjection/LogTransitionAction.lean).
4. **The four laws of thermodynamics follow in the finite observer model.**
   Given a common reference state and the repaired information visible to all
   observers, the same finite repair rule gives equilibrium, entropy increase,
   heat and work accounting, and the low-temperature limit. The remaining work
   connects energy and clock units to a source-produced physical system. See the
   [observers paper](paper/observers_are_all_you_need.pdf), the Lean
   [conditional-repair proof](Lean/Thermodynamics/FiniteConditionalRepair.lean),
   and the [exact certificate](code/thermodynamics/conditional_repair_certificate.py).
5. **The Standard Model gauge structure from twelve ports.** The twelve-port
   geometry, complete reversible response, and observer agreement recover the
   symmetries of the strong, weak, and electromagnetic forces. A separately
   specified matter structure gives
   their familiar global form
   $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. The finite result also lacks
   the extra symmetries responsible for proton decay in minimal grand
   unification. Deriving the matter structure and physical gauge fields from
   the source is work in progress. See the
   [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf),
   the Lean [gauge proof](Lean/Screen/A5OPH.lean), and the
   [global-form proof](Lean/Screen/Z6Descent.lean).
6. **One generation of matter out of a finite search.** An exhaustive scan of
   the declared possibilities leaves the fifteen particle states and exact
   charges of one Standard Model generation, with all anomalies cancelled.
   A separate finite selection gives a rank-three family candidate. The
   remaining work attaches these structures to physical particles and excludes
   extra light sectors. See the
   [particle paper](paper/deriving_the_particle_zoo_from_observer_consistency.pdf),
   the Lean [matter-selection proof](Lean/Screen/ExteriorSelection.lean), and
   the [family-band proof](Lean/Screen/A5FamilyBand.lean).
7. **The Koide lepton relation comes out as a theorem.** A Hermitian $C_3$
   response gives the exact positive-chamber relation among the electron,
   muon, and tau masses. With two masses supplied, the formula fixes a
   72-eV-wide interval centered at $1776.969027$ MeV, compatible with the
   measured tau mass. This is a target-informed conditional postdiction because
   its balance premise comes from the known lepton pattern. Source derivation
   of that premise is work in progress. See the
   [Koide paper](extra/koide_identity_from_positive_c3_face_circulants.pdf),
   the [Lean proof](Lean/ObserverPatchHolography/KoideCirculant.lean), and the
   [comparison ledger](docs/POSTDICTION_LEDGER.md).
8. **A frozen fingerprint in how waves travel.** Two declared twelve-port wave
   rules fix distinctive directional patterns whose first anisotropy appears
   at sixth order. Their ratios and rejection rules sit in pre-comparison
   cryptographic custody. Physical-field attachment of a source-selected rule
   is work in progress. A sufficiently sensitive propagation measurement can
   then rule out that branch. See the
   [screen microphysics paper](paper/screen_microphysics_and_observer_synchronization.pdf),
   the [exact receipts](code/a5_fingerprint/runtime/), and the
   [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md).

Beyond the eight receipts, the exact layer carries a set of further
results, each stated with its boundary at the link:

- A signed-graph theorem proves the screen has no free excitation at zero
  cost: on a target-clean source capture the declared signed operator obeys
  $\lambda_{\min}\geq24^{-8661}>0$. It supplies no physical clock or
  particle mass. See the
  [screen microphysics paper](paper/screen_microphysics_and_observer_synchronization.pdf)
  and the pinned
  [source-gap receipt](https://github.com/muellerberndt/oph-physics-sim/blob/d99ca548a4853e83f819a3a2c9d813f7a3429bdb/data/local_domain/source_gap_receipt.json).
- A finite capacity theorem maximizes generalized entropy at $\log M$,
  gives the exact shock shift $\log(1-f)$, and fixes the pure de Sitter
  relation $\mu^2=d-2$; the physical time-advance reading waits on its
  stated dictionaries. See the
  [focused de Sitter paper](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
  and its [Lean proof](Lean/ObserverPatchHolography/DeSitterCapacityShock.lean).
- On separately declared Maxwell, Yang--Mills, and Einstein branches, the
  quadratic kernels have zero hard mass parameters and the expected
  transverse or transverse-traceless classical modes. These are classical
  carrier statements, not quantum pole predictions. See the
  [forced-structure ledger](docs/POSTDICTION_LEDGER.md#forced-structure).
- The finite completion layer proves one public endpoint across completed
  schedules and representatives, an internal Lorentz model with a bounded
  soldering contract, and a proof-carrying regional-net interface with exact
  witnesses and obstruction theorems; the source attachments stay open and
  typed. See the
  [public-world endpoint proofs](Lean/Tower/FixedPointEndpoint.lean), the
  [geometry theorem stack](Lean/Geometry.lean), the
  [finite regional-net interface](Lean/QFT.lean), and the
  [Lean boundary notes](Lean/docs/) for the exact scopes.
- The library defines what counts as an observer through seven tests:
  bounded window, readback, stable records, action, forward prediction, and
  survival under refinement. A worked example passes all seven and three
  deliberately broken systems fail exactly one each. See the
  [operational observer proof](Lean/Tower/OperationalObserver.lean).

The supporting Lean library contains more than 4000 theorems and lemmas and no
admitted proofs. Explicit axiom reports cover the audited theorem subset.
Twenty-three finite proofs use `native_decide`; their generated native-code
evaluation axioms extend the trust base beyond kernel-only checking. See
[Lean/](Lean/).

The rest of this README is the architecture those receipts come from.

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
spacetime and Einstein paper.

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
        ├─ transportable sectors → independent Tannaka compact-group route
        └─ local 12-port carrier → exact inverse-port response theorem
                → A1/A2 theorem forcing the abstract compact Lie type
                → conditional matrix current and rank-15 matter construction
                → exact Z6 kernel on declared tensors
                source current, matter action, global form, scalar, spectrum,
                  and family attachments open
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

The Einstein branch is instrumented end to end: every clause of its
antecedent is either a proved theorem or a measured quantity with a
fail-closed instrument and adversarial controls, never an assumption. Direct
measurement supplies the strongest empirical result in this corpus. Across
runs with 16,384, 65,536, and 262,144 carriers, the held-out event forms
carry Lorentzian signature $(1,3)$ with cone margins $-5.62$, $-3.22$, and
$-1.41$ and decreasing coupling spread, while a same-size control with a
narrow support drops to signature $(2,2)$. The measurements establish
reproducible sensitivity to support structure; they do not establish a
convergence law. Primary data sit in
[evidence/einstein_convergence](evidence/einstein_convergence/), bit-for-bit
reproducible from the
[simulation repository](https://github.com/muellerberndt/oph-physics-sim).

The carrier geometry then does surprising exact work. Oriented incidence
alone derives the antipodal pairing, the proper $A_5$ action, the rank-three
icosahedral frame, and a target-blind inverse-response law. The complete
reversible response clause in A1 and endogenous transport in A2 then force
the local Lie type of the Standard Model:

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5,
\qquad
(P_{12},[\ ,\ ]_\Theta)
\cong\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

No gauge group is assumed anywhere in the premises; a centreless alternative
is excluded by the one-dimensional fixed space. The released matrices are an
exact conditional witness, and the classified bracket families carry a
forced kinetic sector: block-diagonal in the certified projectors, one
coupling ray per simple factor, and equal to the derived action of its own
source dynamics.

Inside the declared exterior-response algebra, an exhaustive 1024-subset scan
leaves exactly one chiral anomaly-free selection: fifteen states with the
charges of one Standard Model generation and a common $\mathbb Z_6$ kernel,
so the maximal faithful image is
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. Under two further named
premises the exact band-cost order $5-\sqrt5<6<5+\sqrt5$ selects the
rank-three family band.

The boundaries stay explicit throughout: the matrix current, matter action,
and global-form selection are not reconstructed from source histories,
and laboratory identification, family attachment, scalar multiplicity, and
the closure packets are open. The
[claim scoreboard](tracking/claims_scoreboard.md) carries every premise;
the exact bracket-search and Jacobi receipts sit in
[code/a5_closure/](code/a5_closure/).

## Claim Scope

The [claim scoreboard](tracking/claims_scoreboard.md) states the scope,
premises, and evidence class of each branch. This README concentrates on the
strongest exact and measured receipts.

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

about $0.63$ and $0.39$ percent below the Planck base-$\Lambda$CDM comparison
value $3.3129271\times10^{122}$. The theory does not select between the two
corrections, and both numbers were computed after the comparison value was
known, so neither is a prediction. The
[claim scoreboard](tracking/claims_scoreboard.md) states what each step
assumes and what is missing.

## Technical status

The eight receipts above are the reader-facing summary. Exact premises,
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

The bonus is quantitative: if the loop closes, $P$ and $N$ cannot be
arbitrary. They must satisfy self-referential closure conditions: the cell
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
existence. That is the twist the book saves for late in the story, where it
belongs, after the observers-first reconstruction stands on its own. None of
the results above depend on it.

## Open Proof Obligations And Falsification Boundary

The direct $N$ theorem contains a finite, source-derived simulator
public-checkpoint packet whose capacity theorem reduces to
$M_0=|X_{\rm reach}|$, with exact bounded countermodels showing which
completions the finite controls do not pin down. A physical $N$ theorem
requires a complete source antecedent, one physical zero, proof that
both sides read the same universe-level quantity, and the physical carrier
attachment; an independent finite $A_5$ control proves that raw equality at
one rung is not physical closure.

The other named obligations are:

- complete the capacity source antecedent and select one positive physical
  carrier; the horizon-record identification is not evaluable without it;
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
- [`book/`](book): the book source and downloadable PDF.
- [`cosmology/`](cosmology): dark-sector and cosmology research.
- [`physics-problems/`](physics-problems): focused applications and open-problem notes.
- [`docs/`](docs): claim policy, falsification program, and technical audit material.
- [`assets/`](assets): diagrams and public figures.

The simulation source is maintained in the companion
[oph-physics-sim](https://github.com/muellerberndt/oph-physics-sim)
repository, which produces the simulation receipts and evidence artifacts
cited here.

## Explore OPH

- [The book, web edition](https://oph-book.floatingpragma.io)
- [The book, print PDF](https://cfxrbtseaimxxqsxlrku.supabase.co/storage/v1/object/public/books/reverse-engineering-reality.pdf)
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
