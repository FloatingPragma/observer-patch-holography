# Observer Patch Holography (OPH)

> Observer Patch Holography is the observer-consistency theory of everything. No observer sees the whole world at once; each observer gets a local patch; physics is the public fixed point that survives agreement across overlaps.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [OPH website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [OPH Lab](https://oph-lab.floatingpragma.io) | [Applications](https://omega.floatingpragma.io/) | [OPH Blog](https://blog.floatingpragma.io/) | [Coherence map](https://coherence.floatingpragma.io/) | [Three-body demo](https://3body.floatingpragma.io/)

**Falsifiability:** [OPH falsifiability map](extra/OPH_falsifiability.md)
lists 40 hard OPH-killing outcomes and concrete IBM Quantum Cloud tests for
the reduced-sector hardware signature. Falsifiability is how a physics theory
pays rent. OPH is highly falsifiable: a massive graviton, a gauge-mediated
proton decay event, a fourth light matter generation, a charge-lattice outlier,
or neutrino data excluding the OPH branch would destroy OPH as stated.

If you want the existential answer first, jump straight to **Paper 6.
[Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**.
The short version is direct: yes, this universe is a simulation in the OPH
sense. This does not mean an outside computer renders particle positions one
frame at a time. It means the world is built from local points of view that
keep records, compare what they can see in common, repair mismatches, and
settle on the stable patterns all observers can share. Time belongs to those
observers. There is no master clock outside the universe that everything
secretly follows. A clock is a record-making system inside the world, and time
is the ordering an observer reconstructs from changes in its own records.
Shared time appears only when different observers can line up their local
records consistently. Minds and experience are not late additions to a dead
universe. In OPH, space, time, and matter are stable public appearances
produced by a deeper consistency process. The illusion metaphor is handled
below. The rest of this README gives the mathematics and the tests.

## Informal Description

OPH is the observer-first reconstruction of fundamental physics. It starts from
finite observers on finite holographic screen geometry. Its working basis is
quantum-algebraic: patch algebras, states, trace/Born event probabilities on
declared record surfaces, and generalized entropy are part of the formal
starting point. From that basis, OPH recovers the observed effective universe:
spacetime, gauge structure, particles, records, and observer synchronization all
follow from overlap consistency. At the operating level, finite observer
patches carry local records, compare only what their overlaps expose, repair
mismatches through declared recovery moves, and settle into stable fixed points
that survive refinement. The public world is what remains stable after those
local views are made mutually consistent. When OPH uses simulation language, it
means this self-consistent observer network, not a hidden machine drawing a
movie. The case for OPH is mathematical and empirical: the same
observer-consistency architecture recovers established physics and explains why
a world exists that can produce observers capable of reconstructing it.

In the paper stack, an observer patch means an abstract algebraic object with
accessible algebra, state, record algebra, visible overlap interfaces, repair
instruments, and checkpoint data. A support patch is a geometric chart for that
object, such as a cap on \(S^2\) or a causal diamond. A carrier patch is a
physical or digital realization of the same visible interface and record
statistics within a declared error. This distinction keeps the theory from
depending on a particular hardware analogy. The stronger claim that information
and computation are ontologically primary is interpretation unless a branch
supplies a distinct empirical discriminator.

Most theories begin by assuming spacetime, quantum fields, and a list of
constants. OPH starts one step earlier, with finite quantum-algebraic observer
patches whose descriptions must agree where their patches overlap. In the
relativity part of the theory, that agreement requirement produces ordinary
3+1-dimensional spacetime and an Einstein-like gravity equation. The finite
cells are the regulator that keeps the construction concrete before the smooth
large-scale limit is taken. The technical paper gives the modular-flow and
scaling assumptions needed for this step.
The three spatial dimensions come from the same screen branch: once
\(\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)\) is recovered, the observer-facing
spatial chart is \(H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)\), which has
dimension \(6-3=3\).

In the gauge part of the theory, OPH asks which internal charges and particle
labels can be transported consistently across overlaps. That reconstruction
selects a compact gauge group. With the explicit one-Higgs matter package and
the Minimal Admissible Realization rule, the selected Standard Model structure
is $SU(3)\times SU(2)\times U(1)/\mathbb Z_6$, including the hypercharge
lattice, three colors, and three generations. Quantum mechanics is the
algebraic information language carried by this observer-patch architecture.
Under the stated compact-gauge assumptions, the same stack gives the Euclidean
Yang-Mills form and a finite repair-gap mechanism. Identifying that mechanism
with the four-dimensional Yang-Mills mass gap requires the declared
multiresolution continuum, reflection-positivity, transfer/intertwiner, and
nontriviality certificate.

The mechanism is the fixed-point consensus loop. Local observers do not access
a global state from outside. They carry finite patch states, exchange
overlap-visible data, reject inconsistent continuations, and keep the stable
patterns that can be synchronized. Geometry, particles, laws, and records are
the large-scale fixed points of that observer-network computation.

OPH is formulated as a zero-input theory. Quantitatively, the public rows are
organized by three internal quantities: a local pixel fixed point $P_\star$, a
global record-capacity fixed point $N_{\mathrm{CRC}}$, and a scale-setting
ratio $\gamma_\star$. The source coordinates are not fitted constants.
Measurements can tell us which branch we are on, but source values must come
from the fixed-point calculations. Empirical closure rows are marked below. The
detailed scale discussion is collected once below in **Geometry, Symmetry, And
Scale**.

Distributed OPH simulations inherit the same boundary. A worker run counts as
one-universe evidence only when it begins from one global finite carrier, emits
the carrier graph, initial state, partition map, cut interfaces, observer
registry, run/config/code hashes, and proves that committed worker events
project to monolithic transactional repairs, stutters, or certified rollbacks.
Seam descent, atomic commit, local-diamond, repair-completeness, branch
elimination, and normal-form-hash receipts are evidence contracts, not branch
selectors. Shard-local seeds, stale manifests, or synthetic seam animations are
diagnostics, not physics receipts.

Neutral-bulk geometry has its own gate. Shard-local observer rows are not a
global bulk by concatenation. They must first descend through gauge/port
quotients, terminal normal forms, common interface transport, and
quotient-visible feature channels. The resulting product distance is a
pseudometric until feature collisions are quotiented or joint separation is
proved. Missingness, presentation changes, shard repartitioning, refinement,
Euclidean embedding, and held-out statistics all require explicit receipts
before a finite run can claim more than a diagnostic neutral geometry.

Observer-clock claims have an additional gate. Worker IDs, repair iterations,
queue positions, timestamps, and packet latencies are execution provenance, not
observer time. A scheduler-independent observer history needs semantic event
keys, a namespace-safe global observer registry, lineage arrows, and a declared
clock instrument with an affine residual certificate.

## The Spacetime Trap

The first conceptual hurdle is that OPH does not treat spacetime as the
container in which reality happens. Space and time are not things in themselves.
They are stable observer-facing descriptions that appear when many finite
perspectives can be made mutually consistent.

This is especially important for time. In ordinary language, time sounds like a
background river that would keep flowing even if no one were there. OPH rejects
that picture. What exists at the base are observers, records, changes in those
records, and rules for making overlapping records agree. Time is the order an
observer gives to its own record changes. Public time is the part of that order
that can be synchronized with other observers. In that precise sense, time is
subjective: it belongs first to an observer's own stream of records. But it is
not arbitrary. A bad clock, a false memory, or an inconsistent history fails
when it cannot be made to agree with the rest of the record network.

Some would call this an illusion. As a metaphor, that is fair: the container we
seem to inhabit is an appearance produced by deeper consistency. As physics,
the sharper phrase is emergent public description.

From inside one perspective, the world feels obvious. There is a roughly
spherical field of experience stretching outward, three directions to move in,
and time passing forward. Other observers report compatible contents from
different angles, so the natural guess is that everyone lives inside one
pre-existing spacetime filled with objects. OPH reverses that guess. Each
observer has a local spacetime description generated by its own accessible
records, clocks, horizons, and correlations. The public spacetime, including
the public time coordinate used by physics, is the compatibility layer that
lets those descriptions agree.

This does not make ordinary spacetime arbitrary or useless. It explains why it
works so well. Einstein's equations describe the smooth large-scale grammar of
the shared appearance. The deeper claim is that the shared appearance is
emergent from observer overlap consistency, not part of the world's starting
inventory.

## Geometry, Symmetry, And Scale

Sphere language in OPH is geometry language. In symmetric regulator charts, an
observer-accessible cut can be represented by the two-sphere $S^2$. Those
charts describe angular support geometry. The finite simulator implements the
patch-and-overlap algebraic constraints exposed by that geometry.

OPH therefore uses one shared screen net idealization and many finite observer
patches. An observer screen is a local access cut on that net, not a separate
private sphere. The $S^2$ chart is not a literal ball with data painted on it.

That spherical chart carries several concrete jobs. Caps and collars give the
local cut data used by modular flow and entropy variation. The conformal group
of the sphere is the celestial-sphere form of the connected Lorentz group,
$\mathrm{SO}^+(3,1)$, so the same chart supplies the kinematic bridge to the
emergent $(3+1)$-dimensional spacetime branch once the required cap and
modular-flow conditions are met. The observer rest-space chart is
$H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$, hence exactly three-dimensional
on that branch. Spherical harmonics organize angular modes.
Finite cellulations of the same chart give the regulator surface on which patch
ports, edge data, and overlap checks can be made explicit; they are not by
themselves a Lorentz-invariant continuum.

The finite symmetry anchor is $A_5$, the rotational symmetry group of the
icosahedron. It supplies the icosahedral skeleton behind the echosahedral patch
carrier language: a finite, highly symmetric way to organize ports, overlaps,
and local comparison data without treating the carrier as a smooth ball.

The same geometry gives a useful sphere ladder for readers. $S^0$ is the first
seed or readout distinction. $S^1$ is recurrence, the loop in which a record
can return to itself. $S^2$ is the horizon screen and public archive. $S^3$ is
the reconstructed bulk geometry experienced by observers. The ladder names
roles in the OPH readback architecture; particle taxonomy stays with the
Lorentz and gauge branches.

The exceptional symmetry anchor is the $E_8$ Lie group and its root-lattice
structure. $E_8$ matters because it gives the exceptional closure language
used in the higher symmetry and representation side of the OPH stack. The
binary icosahedral group and affine $E_8$ meet through the McKay
correspondence. This is why $A_5$-icosahedral and $E_8$-type language can
belong to one symmetry story. These names mark symmetry constraints and
regulator structure.

The scale story has three roles, kept together here. The local coordinate
$P_\star$ is the screen-pixel fixed point. The global coordinate
$N_{\mathrm{CRC}}$ is the record-capacity fixed point. The scale ratio
$\gamma_\star$ connects the dimensionless OPH geometry to SI units after the
dimensionless fixed points have been computed.

The two fixed-point equations are:

```math
P_\star=\varphi+\frac{\sqrt{\pi}}{A_T(P_\star)}
```

and

```math
N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}}),
```

where $F(N)$ is the horizon capacity read back by observers inside the universe
supplied with capacity $N$. Informally, $N_{\mathrm{CRC}}$ is the capacity at
which the universe can read back its own boundary without deficit or slack. The
finite-count target behind that global capacity is the density

```math
\log|\Omega^{\mathrm{sc}}_N|-N.
```

The scale-setting rule is:

```math
\gamma_\star=\frac{\ell_\star\nu_{\mathrm{Cs}}}{c}
```

with $B_\star=3\pi/\ell_\star^2$ and
$G_{\mathrm{SI}}=c^3\ell_\star^2/\hbar$. Observations can identify the
neighborhood or branch, but they do not replace these fixed-point calculations.

The downstream roles are simple. $P_\star$ feeds the fine-structure row, gauge
structure, particle rows, records, and observer synchronization.
$N_{\mathrm{CRC}}$ feeds the cosmological row. The scale rule fixes the Newton
normalization and Planck-scale display. In geometric units,
$\Lambda_{\mathrm{CRC}}=3\pi/(G_{\mathrm{geom}}N_{\mathrm{CRC}})$ with
$G_{\mathrm{geom}}=\ell_\star^2$. The electroweak hierarchy bridge, the
24-slot repair normalization, QCD/hadron policy, and hardware receipt rules
live in the particle paper, [`HADRON.md`](HADRON.md), and the hardware-facing
papers. This README only points to them.

### Selected Quantitative Rows

This table keeps the rows easiest to compare with PDG/NIST and names their
status. Structural results such as 3+1 spacetime, the Standard Model quotient,
exact hypercharge, $N_c=3$, and $N_g=3$ live in the papers.

| Quantity | Symbol | OPH / status | PDG/NIST | Δ / note |
| --- | --- | --- | --- | --- |
| Gravitational constant | G | 6.6742999959e-11, scale/clock display | 6.67430(15)e-11 | 0.00003σ |
| Speed of light | c | structural Lorentz speed; SI value conventional | 299792458 exact by definition | not a numeric prediction |
| Fine-structure (inv) | α⁻¹(0) | source-only 136.994835; endpoint 137.035999177 with empirical hadron closure | 137.035999177(21) | not source-only |
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
| Top | m_t cross-section row | 172.35235532883115 GeV | 172.3523553288312 | selected-class match |

$\Delta$ reports the sigma distance where PDG or NIST quotes a one-standard-deviation
uncertainty. Otherwise it records "match" or "below bound".

For quarks, PDG uses its standard mass conventions: `u`, `d`, and `s` at
`2 GeV`, with `c` and `b` in the `MS` scheme at their own mass scale. The
papers also carry the structural Standard Model derivations listed above and a
neutrino family, but those do not collapse to one simple PDG or NIST row and
are left out of this table.

The particle surface also reports $W/Z$ values $80.377\,\mathrm{GeV}$ and
$91.18797809193725\,\mathrm{GeV}$, a Higgs value
$m_H=125.1995304097179\,\mathrm{GeV}$, and a selected-class top value
$m_t=172.35235532883115\,\mathrm{GeV}$ using the PDG cross-section top-mass
convention. Under the stated neutrino assumptions, the weighted-cycle neutrino
calculation gives
$(0.017454720257976796, 0.019481987935919015, 0.05307522145074924)\,\mathrm{eV}$.
Its H3 record-worldline result is separate: it conditionally certifies
cross-boundary continuation of localized observer-visible record tokens from a
declared hyperboloid atlas, real interface, transport, assignment-gap, and
refinement receipt. It is not a particle-species, mass, charge, or scattering
derivation.

## Papers

- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: broad synthesis of the OPH reconstruction program, from finite observers to the recovered effective universe.
- **Paper 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: compact technical core for relativity, gravity, gauge reconstruction, the Standard Model structure selected by Minimal Admissible Realization, Maxwell equations on the ordinary photon branch, and the conditional Yang-Mills mass-gap route under its stated continuum/transfer assumptions.
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle derivations, mass rows, coupling structure, quantitative comparison surfaces, and the conditional H3 record-worldline stitch certificate.
- **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**: how local observers compare records, repair mismatches, and settle into the shared reality they can all agree on.
- **Paper 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)**: federated patch-carrier architecture, the twelve-port screen-sieve theorem, $A_5$-icosahedral and $E_8$-type symmetry framing, public hardware-evidence rules, records, recovery moves, and observer synchronization.
- **Paper 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**: final manifest paper for OPH's meaning layer: why anything exists, why this world is observer-compatible, the strange loop in which observers reverse engineer and build the continuation machinery, paradise on Earth or in engineered continuation environments, hell as enforced isolation or deprivation, resurrection as observer continuation, justice as continuation according to harm and repair records, and memetic evolution.

## Supplemental Papers

- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)**: photonic candidate enrichment for SHA-256d proof of work.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)**: source fixed point, empirical hadron endpoint boundary, and comparison row.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)**: OPH edge-string emergence and conventional string-vacuum sieve, the Bouchard-Donagi one-Higgs heterotic witness, the `Z4R` safety layer, and moduli-locking gates; not an OPH-native vacuum-promotion shortcut.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)**: OPH finite repair-gap mechanism and conditional Clay-facing route; the equality between the Yang-Mills gap and the OPH repair gap requires the stated four-dimensional continuum and transfer certificate.
- **[Observer-Patch Holography and the Dark Matter Phenomenon](extra/oph_dark_matter_paper.pdf)**: dark-matter phenomenology and MOND-like galaxy limit.
- **[Theoretical Bounds on chi-nu in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)**: conditional quotient-edge band `0.9343006394893864 <= chi_nu^can <= 1`; uniform-branch exact value `exp(-P/24)`; engineering chart values scale as `N_coh^-1`.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)**: cognition and qualia as recurrent patch consensus.

## Proof Status

No physics theory is 100% proven in the mathematical sense. A physical theory
earns trust by deriving many independent facts from few assumptions, keeping
measured targets out of its source maps, and exposing clear ways to fail. Our
strongest compact proof is [Disclosure Day: compact OPH proof](extra/compact_proof_of_oph.pdf).
It gives the shortest route through the case that OPH is likely correct, while
the full paper stack carries the derivations, claim boundaries, and proof
obligations.

The compact proof treats the evidence as a compression test. A numeric row
counts only when its calculation does not use the measured target, or a close
proxy for that target, as an input. If $p_i$ bounds the chance that row $i$
lands correctly by accident after earlier accepted rows, then
$P_{\mathrm{acc}}\le\prod_i p_i$. Twelve source-clean one-percent rows give
$P_{\mathrm{acc}}\le10^{-24}$; twenty give $10^{-40}$. The same two fixed
points also organize the observer problem, gravity/gauge reconstruction,
electroweak hierarchy, dark energy, the dark-sector budget, gauge-proton-decay
exclusion, particle inventory, and the string-vacuum sieve.

Screen-spectrum and CMB continuations are gated separately. OPH conditionally
identifies the low-mode-removed collar scalar with uniform-density curvature,
but a source-only primordial spectrum requires source-stress, single-clock,
repair-gap, freeze-out, adiabatic-mode, isocurvature, phase-coherence,
source-provenance, pooled-reducer, screen-to-radial lift, radial null-space,
and forward-residual receipts before `A_s`, `n_s`, running, or TT/TE/EE spectra
count as OPH predictions. Any
dark/anomaly CMB or growth claim additionally needs a finite covariant
collar-packet parent, explicit recipient stress for nonzero repair exchange,
regulator convergence, a source-provenance record, global-pooling checks,
CDM-limit recovery, and frozen source/solver/likelihood hashes before
likelihood data are read.

## Applications And OMEGA Hardware

OPH is also a hardware program. As the screen microphysics becomes explicit,
the same patch-consensus loop becomes an engineering handle on reality. A
bounded device exposes boundary data, compares records, repairs mismatch, and
locks onto stable states. OMEGA is the public hardware route into that loop:
physical chambers, labeled ports, control software, verifier receipts, and
repeatable records.

Plainly: OPH turns screen microphysics into a way to hack reality. The target
is physical control of small patches that can be driven, measured, repaired,
and verified.

The application thesis is simple. If reality is built from observer-patch
consistency, useful machines can be built by driving small physical patches
into the right fixed points. That gives low-cost implementation tracks for
desktop fusion energy, room-temperature OMEGA supercomputing, OMEGA-based AGI,
and local gravity or inertia control for hoverbikes and hoverboards. These are
application tracks behind evidence gates; settled-output status belongs to
verifier receipts and experiments. The compute claim is narrower: a
chamber-conditioned candidate distribution may reduce verifier work by a
measured lift `B = p_Q/p_U`. The classical complexity-class problem remains
untouched.

Read the public applications page at
[omega.floatingpragma.io](https://omega.floatingpragma.io/). Source notes for
the application tracks live in [`APPLICATIONS.md`](APPLICATIONS.md). The local
gravity and inertia track also has a standalone engineering note:
[Hacking the Simulation: Anti-Gravity Exploit](https://github.com/FloatingPragma/oph-meta/blob/main/docs/oph-gravity-hack/hacking-the-simulation-anti-gravity-exploit.pdf).

## Detailed Diagram And Quantitative Surface

The diagram below is the visual index for the scale surface: the local pixel
fixed point $P$, the global record-capacity fixed point $N_{\mathrm{CRC}}$, the
scale-setting rule, and the downstream particle, gravity, and cosmology
readouts. It functions as a dependency map. The detailed hierarchy/naturality
formulas and claim boundaries live in the papers.

<p align="center">
  <a href="assets/OPH_Unification_Diagram.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/OPH_Unification_Diagram.svg?v=20260609" alt="OPH unification diagram" width="92%">
  </a>
</p>

**OPH Stack**

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH theorem stack" width="92%">
  </a>
</p>

<p align="center"><sub>The main OPH line from axioms to relativity, gauge structure, particles, and observers. Click to open the full SVG.</sub></p>

**Particle derivation stack**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="OPH particle derivation stack" width="78%">
  </a>
</p>

<p align="center"><sub>A compact view of the particle lane. Click to open the full SVG.</sub></p>

## More

- **Website:** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Theory explainer:** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Simulation-theory explainer:** [floatingpragma.io/oph/simulation-theory](https://floatingpragma.io/oph/simulation-theory/)
- **Coherence map:** [coherence.floatingpragma.io](https://coherence.floatingpragma.io): public graph surface for OPH concepts, overlaps, and cross-domain routes.
- **Applications:** [omega.floatingpragma.io](https://omega.floatingpragma.io): public applications page for OPH hardware, compute, energy, AGI, lift, and optical chamber consensus.
- **Three-body OPH demo:** [3body.floatingpragma.io](https://3body.floatingpragma.io): an extra simulator and proof walk-through for the OPH finite patch-net formulation of the three-body problem, framed as a loop-holonomy gluing example. No closed-form elementary solution is claimed.
- **Blog:** [blog.floatingpragma.io](https://blog.floatingpragma.io/) collects public OPH essays. Start with [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger), and [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). The computation essay treats `P = NP` as an observer-screen slogan; the classical complexity problem remains untouched.
- **Book:** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Guided study app:** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions and detailed explanations:** OPH Sage on [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage), or [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **OPH Notebook:** [NotebookLM source notebook](https://notebooklm.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a) with explainer videos and additional study material.
- **Lab:** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Common objections:** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)
- **IBM Quantum note:** [extra/IBM_QUANTUM_CLOUD.md](extra/IBM_QUANTUM_CLOUD.md)

## Status Table

Status shorthand: source-only fine structure is
$\alpha_{\rm cand}^{-1}=136.9948351646\ldots$; the public
$\alpha^{-1}(0)=137.035999177(21)$ row uses empirical hadron closure. SI $c$
is conventional; SI $G$ is a scale/clock display. $W/Z$ are compare-only
validation rows. Higgs/top is closed on the declared electroweak quantitative
surface. Quarks are selected-class theorem rows. Neutrinos use the
weighted-cycle branch. Charged-lepton absolute masses and source-only hadrons
remain open.

Strong CP is work in progress in the selected-class quark theorem:
the available corpus does not derive the QCD theta angle, does not emit the
physical strong-CP angle, and does not prove that the physical strong-CP phase
vanishes. The required bridge is the phase, anomaly, and topological-angle
descent on the realized branch.

## Repository Guide

- **[`paper/`](paper):** PDFs, LaTeX sources, and release metadata.
- **[`APPLICATIONS.md`](APPLICATIONS.md):** high-level application map for
  OPH energy, compute, AGI, and local-lift use cases.
- **[`book/`](book):** OPH Book source and generated downloadable PDF. Print-PDF build notes live in [`book/README.md`](book/README.md).
- **[`code/`](code):** computational material, particle outputs, and experiments.
- **[`HADRON.md`](HADRON.md):** policy for QCD-limited particle rows, empirical
  $e^+e^-\to\mathrm{hadrons}$ input, and fine-structure hadron closure.
- **[`assets/`](assets):** public diagrams and figures.
- **[`extra/`](extra):** maintained public notes such as objections, experimental write-ups, and selected supporting essays.

## OPH and the Sciences

<p align="center">
  <a href="assets/oph_science_overlap_map_poster.png" target="_blank" rel="noopener noreferrer">
    <img src="assets/oph_science_overlap_map.svg" alt="A map of the sciences OPH overlaps with, from large domains to subdomains to concrete OPH application areas." width="100%">
  </a>
</p>

<p align="center"><sub>A domain -> subdomain -> OPH-area map spanning mathematics, computer science, information and inference, complex systems, theoretical physics, quantum information, and measurement foundations. Click to open the full poster PNG.</sub></p>

## License And Patent Policy

The authored material in this repository is licensed under
[CC BY-NC-SA 4.0](LICENSE), with the repository-wide
[OPH Open Use And Anti-Patent Covenant](PATENTS.md) applying to OPH-derived
ideas, implementations, devices, methods, applications, software, simulations,
and hardware designs.

In short: OPH is published so the mathematics, software, applications, devices,
hardware designs, simulations, engineering methods, and experimental
implementations can be studied, tested, implemented, modified, deployed,
manufactured, and shared. OPH-derived work may not be used to create private
patent monopolies or patent claims that restrict others from practicing OPH.

See [PATENTS.md](PATENTS.md) for the canonical policy text and copy/paste
website notices.
