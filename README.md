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
reconstructs an exact finite structural core: quantum event rules, Lorentz
kinematics on the stated global-support branch, the Standard Model gauge
Lie type, and a conditional one-generation exterior matter pair with an exact
common central kernel.

Three axioms govern the simulator architecture and how observers reach
consensus. Beside them sit two closures. One gives the pixel constant $P$,
tied to the fine-structure constant. The other gives the computational
capacity $N$, tied to the cosmological constant. Since the universe being
simulated and the universe doing the simulating are the same universe, the
simulated fine-structure and cosmological constants must equal the simulating
ones. This self-reference locks in the possible values.

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

## Seven Reproducible Physics Receipts

These public results link directly to their papers, proofs, data, and
certificates:

1. **Three-dimensional space emerges from the algebra of repair records.**
   Observers keep repair records and add them up. Under the declared
   twelve-port repair mean, those sums complete to an ordinary continuous
   three-dimensional Euclidean space, with the sixty proper carrier rotations
   acting on it as isometries. No coordinate
   grid goes in. Physically this would mean distance and direction are
   bookkeeping over comparison records, and space has three dimensions because
   the carrier has twelve ports. A separate finite event instrument measures
   held-out inertia $(1,3)$ at 16k, 65k, and 262k carriers. Calling the
   quotient physical position, gluing overlaps, and fixing physical scale are
   work in progress. See the
   [spacetime and Einstein paper](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf),
   the Lean proofs of the
   [intrinsic completion](Lean/Screen/PrimitivePortFrameQuotient.lean) and
   [repair-response limit](Lean/Screen/PortGramRepairCovariance.lean), the
   [seam-current quotient](Lean/Screen/SeamCurrentCarrierQuotient.lean), and
   [homogeneous internal action](Lean/Screen/SeamCurrentHomogeneousAction.lean),
   the independently verified receipts for the
   [metric quotient](https://github.com/muellerberndt/oph-physics-sim/blob/main/data/repair_closure/port_load_metric_quotient_receipt.json) and
   [proper carrier action](https://github.com/muellerberndt/oph-physics-sim/blob/main/data/repair_closure/port_gram_equivariant_action_receipt.json),
   and the [signature data and regeneration scripts](evidence/einstein_convergence/).
2. **Quantum rules on public records.** Consensus picks out the algebra of
   records that survive comparison, and on a separately declared finite
   algebra-state representation its projectors obey Born probabilities, Lüders
   conditioning, and the Tsirelson bound. Physically this would mean quantum
   probability is the arithmetic of what observers can jointly write down, with
   no extra postulate needed. The declared spinor branch has an exact finite
   candidate at $|S_{\mathrm{CHSH}}|=1+3/\sqrt5>2$, past the classical limit,
   for a setting family the source does not select, so it is not a physical
   Bell prediction. See
   [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf), the
   [consensus paper](paper/reality_as_consensus_protocol.pdf), and the
   [Lean Tsirelson proof](Lean/EventAlgebra/Tsirelson.lean), together with the
   [exact finite candidate receipt](https://github.com/muellerberndt/oph-physics-sim/blob/main/data/quantum/icosahedral_chsh_candidate_receipt.json).
3. **The Standard Model gauge group from twelve ports.** Complete reversible
   port response and endogenous overlap transport force
   $\mathfrak{su}(3)\oplus\mathfrak{su}(2)\oplus\mathfrak u(1)$, with maximal
   faithful image $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$ for the declared
   matter table. Physically this would mean the symmetry group behind the
   strong, weak, and electromagnetic forces is whatever a twelve-port carrier
   can do reversibly, with nothing left to choose. The conditional current
   algebra has no $X/Y$ generators, so the proton-decay channel of minimal
   grand unification is absent, which is narrower than proton stability.
   Matrix-current and
   physical-quotient selection are work in progress. See the
   [conditional current receipt](code/a5_closure/receipts/port_current_inner_reference.receipt.json),
   [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf),
   the [forced-structure scorecard](docs/POSTDICTION_LEDGER.md#forced-structure),
   and the Lean proofs of the
   [A2 holonomy bridge](Lean/Screen/A2HolonomyBridge.lean),
   [gauge trichotomy](Lean/Screen/A5OPH.lean), and
   [finite Z₆ descent](Lean/Screen/Z6Descent.lean).
4. **One generation of matter out of a finite search.** An exhaustive scan of
   the declared exterior-response algebra leaves one charge-conjugate pair of
   chiral, anomaly-free rank-15 projectors, carrying the fifteen Standard Model
   hypercharges of one generation with exact anomaly cancellation. Physically
   this would mean the quarks and leptons of one generation, with their exact
   charges, come out of a finite search instead of a table read off from
   experiment. Under the complete-band and cost premises a separate theorem
   selects rank three, which is where three families would come from. Matter
   attachment, continuum Spin/locality, and laboratory attachment require
   separate source constructions. See the
   [particle paper](paper/deriving_the_particle_zoo_from_observer_consistency.pdf),
   the [finite matter-attachment receipt](code/a5_closure/manifests/matter_attachment_receipt.json),
   the [Lean exterior-selection proof](Lean/Screen/ExteriorSelection.lean),
   and the [Lean family-band proof](Lean/Screen/A5FamilyBand.lean).
5. **A gap that is provably positive.** A target-clean source capture fixes
   the causal order, seam topology, typed sections, 38 frustrated triangles,
   and a zero twisted kernel, and the declared signed operator obeys
   $\lambda_{\min}\geq24^{-8661}>0$, with a numerical refinement of
   $0.1175367$. Physically this would mean the screen has no free excitation at
   zero cost, so disagreement always costs something. This is a signed-graph result, distinct from the
   compact-gauge repair spectrum and the continuum Yang--Mills mass gap, and it
   supplies no physical clock or particle mass. See the
   [screen microphysics paper](paper/screen_microphysics_and_observer_synchronization.pdf)
   and the pinned
   [source-gap receipt](https://github.com/muellerberndt/oph-physics-sim/blob/d99ca548a4853e83f819a3a2c9d813f7a3429bdb/data/local_domain/source_gap_receipt.json).
6. **The Koide lepton relation, with a frozen tau test.** A Hermitian $C_3$
   response obeys $Q=1/3+(2/3)(|b|/a)^2$, so $Q=2/3$ exactly when
   $|b|/a=1/\sqrt2$ in the nonnegative-eigenvalue chamber. Under the declared
   balance and ordering premises, the electron and muon masses fix the tau mass
   within 72 eV, 0.43 standard deviations from the comparison value. Physically
   this would mean the old numerical coincidence among the charged lepton
   masses is one condition on one response matrix, with two masses fixing the
   third. The rejection rule is frozen, so a shifted tau measurement can kill
   it. See the
   [Koide paper](extra/koide_identity_from_positive_c3_face_circulants.pdf),
   [Lean proof](Lean/ObserverPatchHolography/KoideCirculant.lean), and
   [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md).
7. **A frozen fingerprint in how waves travel.** Carrier symmetry fixes the
   wave action on the same three-dimensional carrier, which pins the dispersion
   relation to exact numbers, with one length $a$ and nothing to tune. For
   $q=ak$ on $0\leq q\leq1$, exact arithmetic proves
   $(19/20)q^2\leq\widehat\Lambda\leq q^2$ and
   $\lvert\widehat\Lambda-P_6\rvert\leq(7/388800)q^8$ for
   $P_6=q^2-q^4/20+(1/840-I_6/12600)q^6$, with $I_6$ the icosahedral rank-six
   harmonic. Physically this would mean the vacuum carries a grain, so waves
   run slightly slow at short wavelength and the first direction-dependent
   effect shows up at sixth order in momentum, along icosahedral axes. The same
   action carries a two-polarization massless oscillator, and its upper bound
   forbids that photon from decaying into an electron-positron pair. The
   coefficients and the decision rule are frozen under cryptographic custody
   ahead of any eligible comparison, so a propagation
   measurement with registered exclusion power can refute this branch. Maxwell
   theory, a gauge quotient, and the attachment to a physical photon are work
   in progress. See the
   [screen microphysics paper](paper/screen_microphysics_and_observer_synchronization.pdf),
   Lean proofs of the
   [Dirichlet action](Lean/Screen/SeamCurrentDirichletGenerator.lean) and
   [conditional transverse oscillator](Lean/Screen/SeamCurrentFreePhotonLift.lean),
   [conditional photon-lepton kinematics](Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean),
   the exact [coefficient-ray receipt](code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json),
   [unit-domain remainder receipt](code/a5_fingerprint/runtime/fz12_full_symbol_remainder_receipt.json),
   [joint-threshold receipt](code/a5_fingerprint/runtime/fz12_joint_threshold_receipt.json),
   [synthetic recovery receipt](code/a5_fingerprint/runtime/fz12_synthetic_recovery_coverage_receipt.json),
   [frozen custody packet](https://github.com/FloatingPragma/oph-meta/tree/25da61a800226e0232336ccc86de8dec7d6b51c6/falsification/frozen_targets/fz12_2026-08-02),
   [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md),
   and the [physical propagation](https://github.com/FloatingPragma/observer-patch-holography/issues/666)
   and [comparison](https://github.com/FloatingPragma/observer-patch-holography/issues/667) contracts.

A separate finite theorem maximizes generalized entropy at $\log M$, gives
the exact shock shift $\log(1-f)$, and fixes the pure de Sitter relation
$\mu^2=d-2$. Its physical time-advance reading requires the stated horizon,
observer-mass, gravitational, gauge-mode, and kinetic dictionaries. See the
[focused de Sitter paper](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
and its [Lean proof](Lean/ObserverPatchHolography/DeSitterCapacityShock.lean).

The structural ledger records three further action-level consequences. On
separately declared unbroken Maxwell, perturbative Yang--Mills, and pure
Einstein branches, the quadratic kernels have zero hard mass parameters and
the expected transverse or transverse-traceless classical modes. These are
classical carrier statements, not quantum photon, gluon, or graviton pole
predictions. See the
[forced-structure ledger](docs/POSTDICTION_LEDGER.md#forced-structure).

The supporting Lean library contains more than 1700 theorems and lemmas, with
per-theorem axiom reports and no admitted proofs. See [Lean/](Lean/).

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

The Einstein branch is instrumented end to end. Every clause of its
antecedent (geometric modular normalization, GNS cyclicity and modular
intersections, the Lorentzian event cone, same-source stress and coupling) has
a machine-certified fail-closed instrument with adversarial negative controls
and semantic countermodels, so each clause is either a proved theorem or a
measured quantity, never an assumption. Two clauses are theorems: coupling
universality holds with zero spread for every icosahedrally
symmetric source law, and generator positivity holds by construction for the
declared law family. Direct measurement supplies the strongest empirical result
in this corpus: the Einstein-cone scale path. The selected configurations use
$(16{,}384,128,96)$, $(65{,}536,256,96)$, and
$(262{,}144,512,384)$ for carrier count, observer count, and support width.
Their held-out event forms have Lorentzian signature $(1,3)$, with cone
margins $-5.62$, $-3.22$, and $-1.41$ and decreasing coupling spread. A
same-size control at 262,144 carriers uses support width 96. Its
cross-observer edge count is 312 instead of 1,062, and its signature is
$(2,2)$ instead of $(1,3)$. These measurements establish reproducible
sensitivity to support and cross-read structure under the archived
configurations. They do not establish a fixed-density convergence law or an
infinite-scale limit. The primary data are stored in
[evidence/einstein_convergence](evidence/einstein_convergence/) and every
number reproducible bit for bit from the
[simulation repository](https://github.com/muellerberndt/oph-physics-sim).
Two measured clauses are open: cap-state modular temperature and a
preregistered larger-rung test of the event form. Both carry frozen verdicts.

The evidence stack combines exact finite derivations, machine-checked proofs,
and deterministic measurements with primary data. Mathematical statements,
conditional physical readings, and measured properties carry separate claim
classes in the [claim scoreboard](tracking/claims_scoreboard.md).

The carrier geometry then does surprising exact work. On the certified
echosahedral lineage, the declared integer atom-counting grammar and normalized
Hilbert--Schmidt readback cost give the exact twelve-unit split and gap two.
Deriving that counting grammar and physical cost from the full three-axiom
schema is open. Oriented incidence independently derives the antipodal
pairing, proper $A_5$ action, rank-three icosahedral frame, and the decomposition
$\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5$. Incidence also fixes
the unique nonidentity central graph involution $J$. A target-blind protocol
injects an impulse at every port, reads the adjacency history through graph
diameter, and solves the common farthest-shell filter. It derives
$10J=A^3-4A^2-5A+10I$. The response $R=-J$ has exact relative sector signs;
its common sign is charge conjugation.

The complete reversible response clause in A1 and endogenous proper-carrier
transport in A2 turn the twelve-dimensional port tangent into a compact
current algebra with inner $A_5$ action. Its one-dimensional fixed space and
compact-simple classification force the abstract type
$\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. The released
charged-double-triplet matrices realize that type exactly, while ordered
source tomography and same-current holonomy are work in progress. The
registered recurrence generates a four-dimensional commutative word algebra;
it does not supply twelve current generators, their bracket, or nonidentity
proper rechartings. Twelve target-free diagonal port phases have rank twelve
and commute. Adding the connected adjacency tangent generates
$\mathfrak u(12)$ with derived rank 143, so the source current requires a
non-diagonal response law with derived rank eleven. Conditional on the
canonical oriented carrier, the complete target-free
$A_5$-equivariant alternating-bracket space over $\mathbb Q$ has dimension
fourteen and an exact rational Reynolds basis. The complete Jacobi condition
has 38 independent quadratic coefficient rows, with an exact 11+27 rowspace
decomposition. The residual solution variety, compactness, source
reconstruction, and same-current holonomy are undischarged.
See the
[exact alternating-bracket search packet](code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json)
and the
[Jacobi reduction receipt](code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json).

Inside the declared exterior-response algebra, the exhaustive 1024-subset
scan leaves one unordered conjugate rank-15 pair as the unique nonempty chiral
anomaly-free selection. Anomaly freedom gives determinant balance and
primitive charges up to charge conjugation. Exhaustive central-action
calculation gives a common $\mathbb Z_6$ kernel on those declared tensors, so
their maximal faithful image is
$(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. The cover and its
$\mathbb Z_2$ and $\mathbb Z_3$ quotients carry the same local tensors. The
six-axis calculation has order six only after its diagonal and zero-sum
coefficient relations are declared. A complete source character category and
a same-source loop-to-kernel identification are required to select the
physical quotient. Continuum fermion typing and laboratory-current
identification are open. The transportable-sector/Tannaka construction is a
separate compact-group route, and the source identification of the two routes
is open.

The carrier theorem and its declared matrix witness can be written as

$$
P_{12}\cong_{A_5}\mathbf1\oplus\mathbf3\oplus\mathbf3'\oplus\mathbf5,
\qquad
(P_{12},[\ ,\ ]_\Theta)
\cong\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2).
$$

This abstract Lie type follows from the A1 response and A2 endogenous
transport clauses together with the finite carrier theorem. It does not
follow from the module decomposition or target-blind inverse-port readback
alone. The released matrices are a conditional witness. The matter and
descent certificates prove the corresponding conditional representation and
kernel arithmetic. Source realization, physical global-form selection, and
laboratory identification are separate tests.

The exact carrier results retain explicit physical boundaries. The matrix
current, matter action, and global-form selection are not reconstructed from
source histories. Laboratory current identification,
physical family attachment, exclusion of extra light sectors, scalar
multiplicity,
the Einstein source tower, and the physical closure packets are not derived.
The CP and weak-sector clauses give the exact
conditional window $3\le N_g\le5$ without selecting within it. Under two
additional named premises, the exact band-cost order
$5-\sqrt5<6<5+\sqrt5$ uniquely selects the rank-three screen band. A finite
unitary response receipt recovers the same band at the lowest positive
generator frequency. The physical matter-family attachment, Spin/locality
data, and exclusion of extra light sectors are open. Local
icosahedral incidence constrains the carrier, while the federation nerve
requires its own construction.

## Claim Scope

The [claim scoreboard](tracking/claims_scoreboard.md) states the scope,
premises, and evidence class of each branch. This README concentrates on the
strongest exact and measured receipts.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## The Two Constants: P and N

**$P$ is the local pixel ratio**: the size of the elementary observation cell
in natural geometric units, informally the universe's **resolution**. OPH
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
declared self-read map. The detuning and inside/outside identification laws
are architectural closure premises rather than theorems of the three axioms.
The fixed-point theorem states that a self-map of the physical interval with
contraction constant below one has exactly one fixed point. Outward-rounded
interval certificates verify those hypotheses for each declared $P$ map and
exclude a second root across its full analytic domain. The
[claim scoreboard](tracking/claims_scoreboard.md)
states the root, external comparison, residual, and claim class. The
comparison uses $P_C$, which is defined from the measured endpoint.
Source-derived same-scheme hadronic transport is absent. The registered
comparison has diagnostic status, with a physical
fine-structure constant claim outside its scope.

**$N$ is the public-record capacity** of the whole observer system: how much
correctable memory the substrate carries. It sits opposite $P$, tied to the
cosmological constant rather than to the fine-structure constant.

The direct route reads $N$ off the universe itself. The self-read condition
$N=\log M_0(\mathfrak U_N)$ asks the capacity handed to a trial universe to
match the record capacity reconstructed inside it, and if both sides are
readings of one quantity, self-reference forces them to agree. The proof that
they are one quantity does not exist yet, so this route currently returns no
number. Nothing else in the reconstruction waits on it.

A second route goes through $P$. At the source-forward pixel value the
uncorrected capacity is
$N_0=\pi\exp[6\pi/(P\alpha_U(P))]=3.5321315\times10^{122}$. Two ways of
applying the finite survival correction to it give

$$
N_{\rm pres}=N_0\left(1-\frac{P}{24}\right)=3.2920979\times10^{122},
\qquad
N_{\rm Pois}=N_0e^{-P/24}=3.3000722\times10^{122},
$$

about $0.63$ and $0.39$ percent below the Planck base-$\Lambda$CDM comparison
value $3.3129271\times10^{122}$. The theory does not yet pick between the two
corrections, and both numbers were computed after the comparison value was
known, so neither is a prediction. The
[claim scoreboard](tracking/claims_scoreboard.md) states what each step
assumes and what is still missing.

## Results At A Glance

| Result | What OPH contributes | Main source |
| --- | --- | --- |
| Finite observer consensus | Terminating repair, protected readout, schedule-independent quotient normal forms, and central records | [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf) |
| Conditional quantum event surface | Consensus selects the finite commuting public-record algebra. Given a declared finite algebra-state and two-wing representation, its projectors obey Born probabilities, Lüders conditioning, and the Tsirelson bound. A declared binary-icosahedral spinor branch has an exact finite candidate with $\lvert S_{\mathrm{CHSH}}\rvert=1+3/\sqrt5>2$. Its setting family and completed two-wing instrument are not source-selected, so this is not a physical Bell prediction | [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf) and the [exact candidate receipt](https://github.com/muellerberndt/oph-physics-sim/blob/main/data/quantum/icosahedral_chsh_candidate_receipt.json) |
| Finite local action domain | One target-clean source capture carries an exact causal order on 2,304 events, six closed observer neighborhoods, a sign-frustrated seam complex, typed scalar, chiral, and gauge sections, deterministic integer operator checks, and an exact zero-kernel theorem. An isolated rerun reproduces canonical receipt content. Its declared unit-counting signed seam operator has a rigorously positive finite-domain gap; the numerical refinement is 0.1175367. This operator is distinct from the compact-gauge repair generator used in the conditional Yang–Mills branch. One neighborhood has Euclidean fitted inertia and every cone margin is negative, so the receipt does not establish a continuum spacetime, physical clock, or mass scale | [Screen microphysics](paper/screen_microphysics_and_observer_synchronization.pdf) |
| Relativity | On the certified global support branch with an independently complete algebra-state comparison on the same tower, $\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ and $H^3\cong\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$ | [Spacetime and Einstein paper](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Einstein dynamics | Typed composition from modular flow, null stress, entropy stationarity, and small-ball geometry; construction of one source-derived common-domain tower is work in progress | [Spacetime and Einstein paper](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) |
| Twelve-port Standard Model Lie-type theorem | Oriented incidence gives the proper $A_5$ action and the port module $1+3+3'+5$. Complete reversible port response and endogenous overlap transport make this a compact twelve-dimensional current with inner $A_5$ action. Its one fixed line and compact classification force $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. Target-blind readback separately derives $R=-J$. The registered recurrence has a four-dimensional commutative word algebra. Twelve diagonal port phases commute, while adjoining the connected adjacency tangent generates $\mathfrak u(12)$ with derived rank 143. Neither source lift supplies the required derived-rank-11 current. The released matrix current is an exact conditional realization; non-diagonal source tomography and same-current holonomy are work in progress | [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Conditional Standard Model faithful matter image | On the scan-selected conjugate pair of fifteen-state exterior modules, anomaly balance fixes the primitive charge pair up to conjugation. The exact common kernel on the declared tensors is $\mathbb Z_6$, so their maximal faithful image is $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. The cover and its $\mathbb Z_2$ and $\mathbb Z_3$ quotients carry the same local tensors. The six-axis menu matches $\mathbb Z_6$ only after its coefficient relations are declared. The source does not select the physical global form | [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Matter structure | Exact conditional one-generation exterior modules, hypercharge/anomaly balance, three-color carrier, and the compatible scalar-charge pair and three interaction channels. The CP and weak-sector clauses give $3\le N_g\le5$. Under separate single-band and cost-order premises, an exact finite theorem selects the rank-three screen band and a declared unitary simulator recovers its residue at the lowest positive generator frequency. Tensoring that band with the declared fifteen-state table gives a conditional complex rank-45 candidate. The table carries the nondegenerate chirality grading and exact diagonal $\mathbb Z_6$ action. A separate 8,662-node local-domain receipt checks the declared extension $D_\sigma\otimes I_{45}$ and its conditional inheritance of the positive finite-domain gap. This action is not source-selected. The twelve-port Spin packet and local operator domain have no certified source, domain, or transport bridge. Physical matter-pole identification, continuum Spin/locality, physical seam selection, scalar multiplicity, and exclusion of extra light sectors are open | [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Quantum field-theory landing | Finite-action invariance; exact finite determinant-line and Hamiltonian criteria; formal perturbative restoration and strict finite-order W/Z algebra; separate nonperturbative reconstruction and resonance implications. The exact finite and perturbative routes are parallel descendants of the local action, with source-native constructions as explicit physical gates | [Standard Model gauge paper](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) |
| Finite de Sitter screen | Exact pure-de-Sitter shock normalization, finite entropy maximum, uniform capacity-transfer law for the logarithmic sector coordinate, and analytic curvature; the physical time-advance reading is conditional on the horizon and shock dictionaries stated in the focused paper | [Finite de Sitter capacity paper](extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf) |
| Strict W/Z analytic checks | For a complete renormalized packet, the strict scalar consumer and its order and neutral-mixing rules are exact. The exact quotient of its one-loop-truncated pole coordinates has a passive common-unit scale cancellation at fixed normalized self-energy factors; its strict one-loop re-expansion is stated separately. The remaining coupling ratio and normalized self-energy factors are not source-selected, so no number follows. Interval receipts exclude scalar zeros in the declared principal-sheet boxes and isolate, for each of W and Z, one simple scalar zero with derivative and scalar-residue balls in its declared lower-half pole box on a channel-specific algebraic chart. They identify neither declared chart with the physical resonance sheet and prove no unique continuation, sign bridge, full-matrix Laurent residue, physical-current amplitude, or independent numerical replay. The external fixture is not composed with the OPH electroweak chart, so no physical W/Z pole or mass comparison follows | [Particle paper](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Local $P$ closure | $P=\varphi+\sqrt\pi/A_T(P)$; the fixed-point uniqueness schema and interval certificates give one root for each declared map; physical Thomson transport is work in progress | [Fine-structure constant paper](extra/fine_structure_constant_derivation.pdf) |
| Direct global $N$ readback | $N=\log M_0(\mathfrak U_N)$, with $M_0(q)=\alpha(G_q)$ and $M_0=\lvert X_{\rm reach}\rvert$ on the reversible branch. The fixed $D=24$ packet is exact, but an exact counterfamily shows that base agreement, positivity, and the carrier bound admit completions with different solution sets, so the condition does not yet single out a value. The direct route returns no number until the capacity source is completed | [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf) |
| Conditional common-load $N$ candidates | The screen and electroweak readings must agree once a physical bridge proves they denote one quantity. The exact conditional formulas $N_{\rm pres}=N_0(1-P/24)$ and $N_{\rm Pois}=N_0e^{-P/24}$ evaluate to approximately $3.2920979\times10^{122}$ and $3.3000722\times10^{122}$, against the weighted Planck base-$\Lambda$CDM value $3.3129271\times10^{122}$, residuals $-0.63$ and $-0.39$ percent. Exact countermodels show the finite survival datum selects neither correction, and both comparisons are retrospective, so neither row is an OPH prediction | [Deriving the Particle Zoo](paper/deriving_the_particle_zoo_from_observer_consistency.pdf) |
| Exact verification | Interval certificates, finite receipts, and reproducible simulations | [`code/`](code) |

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
closure is machine-checked in Lean. The declared $P$ map has a certified fixed
point, while its comparison with the physical fine-structure constant has
diagnostic status. The evaluation boundaries of the closure conditions and
their missing physical inputs are stated in the
[OPH Falsification Program](docs/OPH_FALSIFICATION_PROGRAM.md).

A physical closure of both constants would give a zero-continuous-parameter
branch with both values returned by the architecture. That physical
attachment is open. The fixed-point theorems certify roots of declared maps;
they do not turn an observed basin or target-defined coordinate into a
physical derivation. On the $N$ side the finite counting is exact, but the
capacity source it would close over is incomplete, so the direct condition is
not yet evaluable and the common-load route stays conditional on its physical
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
public-checkpoint packet. At fixed $D=24$, the packet has the reachable public
records, the
publicness rule, joint checkpoint kernels, carrier projections, and extension
and refinement maps. Injective checkpoint generators reduce its capacity
theorem to $M_0=|X_{\rm reach}|$, computable by exact CSP or model counting.
The target-clean all-rung counterfamily has an exact bounded verdict: base
agreement, positivity, and the carrier bound admit completions with different
slack-zero sets. Executable certificates check additional finite controls at
declared rungs. Universal all-rung membership of those countermodels in the
complete A1--A3 terminal fibers, atom maps, joint kernels, meaning maps,
feasible sets, and regulator controls has not been proved, and no
executable-to-Lean bridge supplies it. A physical \(N\) theorem requires a
complete source antecedent, one physical zero, proof that both sides read the
same universe-level quantity, and the physical carrier attachment. A
separately named stronger source law is one possible route. The independent finite $A_5$ control has $M_0=60$ and
$D_{\rm raw}=60k$; its publicly inert multiplicity proves that raw equality at
$k=1$ is not physical $N$-closure.

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
| The full observer-first synthesis | [From Observer Consensus to Standard Physics](flagship/from_observer_consensus_to_standard_physics.pdf) |
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
