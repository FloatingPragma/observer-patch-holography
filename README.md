# Observer Patch Holography (OPH)

> Observer Patch Holography is the observer-consistency theory of everything. No observer sees the whole world at once; each observer gets a local patch; physics is the public fixed point that survives agreement across overlaps.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [OPH website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Book: Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [Ω](https://omega.floatingpragma.io/) | [Blog](https://blog.floatingpragma.io/) | [Tech](https://omega.floatingpragma.io/) | [Simulation](https://simulation.floatingpragma.io)

**Principles and evidence:** OPH is stated under the strange-loop
principles in [STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md); the central thesis, its lineage in
established physics and mathematics, and the five-pillar mission are stated in
[STRANGE_LOOP.md](STRANGE_LOOP.md); the forcing chain from the principles to a
unique universe — which consistency requirement removes which freedom, and the
uniqueness lemmas that close the chain — is stated in
[CONSISTENCY_STACK.md](CONSISTENCY_STACK.md). The two substrate settings are estimates:
P from measured α (SL-3) and N from measured Λ with the scale bridge (SL-4).
Neither α nor Λ nor G appears on any OPH output list. The evidential surfaces
of the program are the closure ledger ([CLOSURE_LEDGER.md](CLOSURE_LEDGER.md))
and the compression scorecard
([COMPRESSION_SCORECARD.md](COMPRESSION_SCORECARD.md)); every surface in this
repository states claims at or below the standing those two files record.

**Falsifiability:** The [OPH falsifiability map](extra/OPH_falsifiability.md)
organizes the failure conditions as a falsification ladder grouped by domain
(theory and computation, particle and precision physics, neutrinos, cosmology
and gravity), with a risky-prediction status table, sharp numerical tests, and
a best short list. Entries include a
failure of the fully specified pure-Einstein transverse-traceless carrier
receipt (or of a separately completed quantum pole receipt), gauge-mediated
proton decay, a fourth light matter generation, and a charge-lattice outlier.
Executed tests stay on that surface with their verdicts
([STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md) rule 7). The derivation of the physical neutrino mixing
matrix and absolute neutrino masses is work in progress.
◊
If you want the existential answer first, jump straight to **Paper 6.
[Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**.
Yes, this universe is a simulation in the OPH
sense: a world built from local points of view that keep records, compare what
they can see in common, repair mismatches, and settle on the stable patterns
all observers can share. No outside computer has to render particle positions
frame by frame. Time belongs to those
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
starting point. From that basis, OPH proposes conditional reconstruction branches
for spacetime, gauge structure, particles, records, and observer synchronization;
each branch retains the premises and receipt gates stated in the papers. At the operating level, finite observer
patches carry local records, compare only what their overlaps expose, repair
mismatches through declared recovery moves, and settle into stable fixed points
that survive refinement. The public world is what remains stable after those
local views are made mutually consistent. OPH simulation language names this
self-consistent observer network at the operating level.
The case for OPH is mathematical and empirical: the same
observer-consistency architecture is proposed to recover established physics.
Its further account of why a self-consistent world exists and can produce
observers capable of reconstructing it is an interpretation outside the
recovered theorem package.

In the paper stack, an observer patch means an abstract algebraic object with
accessible algebra, state, record algebra, visible overlap interfaces, repair
instruments, and checkpoint data. A support patch is a geometric chart for that
object, such as a cap on $S^2$ or a causal diamond. A carrier patch is a
physical or digital realization of the same visible interface and record
statistics within a declared error. This distinction keeps the theory from
depending on a particular hardware analogy. The stronger claim that information
and computation are ontologically primary is interpretation unless a branch
supplies a distinct empirical discriminator.

Most theories begin by assuming spacetime, quantum fields, and a list of
constants. OPH starts one step earlier, with finite quantum-algebraic observer
patches whose descriptions must agree where their patches overlap. In the
relativity part of the theory, finite consensus supplies quotient normal
forms. The ordinary 3+1-dimensional spacetime readout and the Einstein-like
gravity equation are recovered only on the five-axiom recovered-core branch:
the compact paper supplies geometry readout, controlled modular flow, the
null-stress bridge, bounded-interval transport, admissible fixed-cap entropy
stationarity, small-ball area variation, tensor upgrade, and capacity closure.
The finite cells are the regulator that keeps the construction concrete before
the smooth large-scale limit is taken.
The three spatial dimensions come from the same screen branch: once
$\mathrm{Conf}^+(S^2)\cong\mathrm{SO}^+(3,1)$ is recovered, the observer-facing
spatial chart is $H^3\simeq\mathrm{SO}^+(3,1)/\mathrm{SO}(3)$, which has
dimension $6-3=3$.

In the gauge part of the theory, OPH asks which internal charges and particle
labels can be transported consistently across overlaps. The visible
fixed-cutoff charges generate a tensor category, but its refinement limit is
conditional on an explicit compact-gauge receipt: coherent surjective boundary
group maps, finite-state extendability, center-compatible finite block
embeddings, finite tensor realizations, and compatible forgetful fibers. On a
cofinal tail carrying that receipt, the reconstruction produces some compact
gauge group. With the explicit one-Higgs matter package and
the Minimal Admissible Realization rule, the selected Standard Model structure
is $SU(3)\times SU(2)\times U(1)/\mathbb Z_6$, including the hypercharge
lattice, three colors, and three generations. Quantum mechanics is the
algebraic information language carried by this observer-patch architecture.
Under the stated compact-gauge assumptions, the same stack gives the Euclidean
Yang-Mills form and a finite repair-gap mechanism. Identifying that mechanism
with the four-dimensional Yang-Mills mass gap requires the declared
multiresolution continuum, reflection-positivity, transfer/intertwiner, and
nontriviality certificate.

The particle-sector numeric status lives in the table under **Selected
Quantitative Values** below; the full audit prose lives in the particle paper
and in [`code/particles/`](code/particles/README.md).

The mechanism is the fixed-point consensus loop. Local observers do not access
a global state from outside. They carry finite patch states, exchange
overlap-visible data, reject inconsistent continuations, and keep the stable
patterns that can be synchronized. Geometry, particles, laws, and records are
the large-scale fixed points of that observer-network computation.

OPH carries two estimated settings, one bridge, and declared architecture
([STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md)): P estimated from measured α under SL-3, N estimated
from measured Λ under SL-4, the cesium SI bridge under SL-5, and the declared
structural selections listed in the scorecard. Standing is measured by the
compression scorecard ([COMPRESSION_SCORECARD.md](COMPRESSION_SCORECARD.md)).
The provenance rule inside that frame: no measured target or fitted numerical
constant may enter a declared source map. The quantum-algebraic axioms and the
$S^2$ screen remain explicit premises. Quantitatively, the public rows are organized by three
internal quantities: a local pixel fixed point $P_\star$, a
global record-capacity fixed point $N_{\mathrm{CRC}}$, and a scale-setting
ratio $\gamma_\star$. A source coordinate earns that label only when its declared
map emits it without a fitted target; otherwise it remains a diagnostic or
comparison coordinate. Measurements can tell us which branch we are on, but source values must come
from the fixed-point calculations. Empirical closure rows are marked below. The
detailed scale discussion is collected once below in **Geometry, Symmetry, And
Scale**.

OPH evidence has the same general shape throughout the project. A claim is
grounded in bounded observer-like patches with local state, explicit
boundaries, readback, records, feedback or repair moves, and public evidence
bundles. The invariant observer-patch history carries the claim across
presentations, coordinate choices, and implementation traces.

## The Spacetime Trap

The first conceptual hurdle is that OPH does not treat spacetime as the
container in which reality happens. Space and time are not things in themselves.
They are stable observer-facing descriptions that appear when many finite
perspectives can be made mutually consistent.

This is especially important for time. Ordinary language treats time as a
background river flowing without observers. OPH rejects that picture. What
exists at the base are observers, records, changes in those records, and rules
for making overlapping records agree. Time is the order an
observer gives to its own record changes. Public time is the part of that order
that can be synchronized with other observers. In that precise sense, time is
subjective: it belongs first to an observer's own stream of records. But it is
not arbitrary. A bad clock, a false memory, or an inconsistent history fails
when it cannot be made to agree with the rest of the record network.

The illusion label works only as a metaphor: the container we seem to inhabit is
an appearance produced by deeper consistency. As physics, the sharper phrase is
emergent public description.

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
charts describe angular support geometry. Finite regulator models implement the
patch-and-overlap algebraic constraints exposed by that geometry.

OPH uses one shared screen net idealization and many finite observer patches.
An observer screen is a local access cut on that net, not a separate private
sphere. The $S^2$ chart is not a literal ball with data painted on it.

That spherical chart carries several concrete jobs. Caps and collars give the
local cut data used by modular flow and entropy variation. The conformal
symmetry of the sphere is the celestial-sphere form of the connected Lorentz
group, so the same chart supplies the kinematic bridge to the emergent
$(3+1)$-dimensional spacetime branch once the required cap and modular-flow
conditions are met. On that controlled branch the observer-facing rest-space
chart is exactly three-dimensional. Caps mark cuts and sides of the chart, not
preferred observer points. Record tokens populate the chart only when
calibrated cap responses and an error budget support localization. Finite
evidence may report ambiguity. Spherical harmonics organize angular modes.
Finite cellulations of
the same chart give the regulator surface on which patch ports, edge data, and
overlap checks can be made explicit; they are not by
themselves a Lorentz-invariant continuum.

The finite symmetry anchor is $A_5$, the rotational symmetry group of the
icosahedron. It supplies the icosahedral skeleton behind the echosahedral patch
carrier language: a finite, highly symmetric way to organize ports, overlaps,
and local comparison data without treating the carrier as a smooth ball.

The same geometry gives a useful sphere ladder for readers. $S^0$ is the first
seed or readout distinction. $S^1$ is recurrence, the loop in which a record
can return to itself. $S^2$ is the horizon screen and public archive. The final
rung means the three-dimensional observer-facing bulk role; on the controlled
Lorentz branch its canonical kinematic chart is $H^3$. No $S^3$ global topology
follows from this mnemonic. Particle taxonomy stays with the Lorentz and gauge
branches.

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

$P_\star$ feeds the fine-structure row, gauge
structure, particle rows, records, and observer synchronization.
$N_{\mathrm{CRC}}$ feeds the cosmological row. The scale rule fixes the Newton
normalization and Planck-scale display. In geometric units,
$\Lambda_{\mathrm{CRC}}=3\pi/(G_{\mathrm{geom}}N_{\mathrm{CRC}})$ with
$G_{\mathrm{geom}}=\ell_\star^2$. The electroweak hierarchy bridge, the
24-slot repair normalization, QCD/hadron policy, and hardware receipt rules
live in the particle paper, [`HADRON.md`](HADRON.md), and the hardware-facing
papers. This README only points to them.

### Selected Quantitative Values

This table keeps the values easiest to compare with PDG/NIST and names their
support status. Structural results such as 3+1 spacetime, the Standard Model quotient,
exact hypercharge, $N_c=3$, and $N_g=3$ live in the papers.
The ledger emits no nonzero source-only physical particle mass. Conditional
chart values appear below with their chart names; each value holds only on its
stated chart, and the strict ledger keeps every such row non-promotable.

| Quantity | Symbol | OPH / support status | PDG/NIST | Δ / note |
| --- | --- | --- | --- | --- |
| Gravitational constant | G | unit-bookkeeping identity on the SL-5 scale bridge: the caesium clock rule carries $G_{\mathrm{geom}}=\ell_\star^2$ to the SI display through $G_{\mathrm{SI}}=c^3\ell_\star^2/\hbar$; G is an input | 6.67430(15)e-11 | bridge input, not a prediction; no sigma applies |
| Speed of light | c | structural Lorentz speed of the recovered spacetime branch; the SI number is a unit convention | 299792458 exact by definition | not a numeric prediction |
| Fine-structure (inv) | $\alpha^{-1}$ | source-only forward map contracts to 136.994835177413, the interval-certified unique fixed point of the declared map (`code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json`; the earlier printed tail 136.994835164622 was an unconverged run, superseded beyond digit 9, ledger row CL-6). The residual to the SL-3 basin is the open hadronic transport term, closure row CL-1 in [CLOSURE_LEDGER.md](CLOSURE_LEDGER.md) (3.0e-4 relative, about 2e6 measurement sigma). The self-consistent gauge-width map has certified fixed point 137.035660136946577: row CL-2, 2.5e-6 relative, the same missing term; the historically displayed 137.03595950081728 is a mixed-provenance display packet, not a fixed point of any single declared map. The empirical hadron-closure endpoint 136.26365894312011 on [136.15838328342668, 136.36909752396569] uses measured $e^+e^-\to$ hadrons data | 137.035999177(21): the SL-3 input that estimates P, on no output list | closure rows, not predictions. The strict source branch needs the Ward-projected transport, blind, against the frozen target in `falsification/frozen_targets/`; the empirical endpoint leaves the same-scheme anchor gap [0.648554111145163, 0.854792066592563] |
| Massive electroweak and Higgs sector (W, Z, Higgs, top) | $W,Z,H,t$ | conditional derivation from the pixel-fixed couplings with no measured mass entering. Strict source-audit branch (zero-selector law) $(m_W,m_Z)=(80.330,\,91.119)$ GeV; the criticality law fixes the top Yukawa from the gauge sector (a determined electroweak-sector quantity, not a free flavor parameter), giving $m_t=172.63$ GeV ($+0.02\%$) with the Higgs companion $m_H=125.77$ GeV at the frozen boundary-scale candidate $E_\star e^{-\pi}P^{-1/6}$; declared-surface fit $(125.20,\,172.35)$ GeV kept separate; no source-only physical mass emitted | $m_W=80.3692(133)$, $m_Z=91.1880(20)$, $m_H=125.13(11)$, $m_t=172.60(30)$ GeV | the strict source-audit $W/Z$ branch is an exact implication of five finite quotient-transport premises whose carrier certificate and discrete two-law selection principle are unemitted; the Higgs/top pair is a zero-continuous-parameter criticality family whose boundary-scale selection is a theorem modulo two finite carrier facts CF1/CF2 (the log-midpoint variational principle is proved and its quadratic-cost and placement premises reduce to the axioms; CF1/CF2 are the same D11 carrier census the $W/Z$ law needs, and the three-loop implied scale is the registered discriminating test); the declared-surface fit is target-anchored and never predicts; the physical scale, scheme, provenance, uncertainty, and complex-pole gates are open |
| Charged leptons (electron, muon, tau) | $e,\mu,\tau$ | no source-only physical mass emitted. The declared MCPR response architecture emits the conditional triple $(m_e,m_\mu,m_\tau)=(0.51096,\,105.649,\,1776.78)$ MeV at a coherent $-84$ ppm offset with zero runtime charged reference; the icosahedral face-corner carrier fixes the family structure | $0.5109989$, $105.6584$, $1776.93$ MeV | the MCPR triple is a declared model-input architecture, not a source-only prediction; the source-only completion is the open A5/W5 orbit program |
| Light-quark spectrum (up, down, strange, charm, bottom) | $u,d,s,c,b$ | no nonzero source-only physical mass emitted: the source equations produce two ordered Yukawa profile rays but leave their spreads free. The top is not in this sector; its Yukawa is fixed by the Higgs criticality law (see the electroweak/Higgs row) | common-scale Yukawa data used only for audit | the five light-quark source equations leave a proven two-modulus non-identifiability; a source-derived flavor-orbit selector is the missing object. This is a proven non-entailment, not an arithmetic gap |
| Strong-coupling scale | $\Lambda_{\mathrm{QCD}}^{(3)}$ | 0.3348 GeV [0.319, 0.350]: dimensional transmutation of the source strong coupling, no hadronic input anywhere | 0.338(12) GeV (published) | central value within ~1% of the published value and inside its uncertainty; conditional on declared threshold inputs (bracket-swept); this is the perturbative half of every light-hadron mass |
| Nucleon mass | $m_N$ | 0.929 GeV [0.823, 1.043]: source $\Lambda_{\mathrm{QCD}}$ times the published lattice-QCD ratio $m_N/\Lambda_{\mathrm{QCD}}\approx2.775$ (`oph_plus_external_qcd_theory`) | 0.9383 GeV | −1.0% (compare-only); conditional prediction: OPH fixes the strong scale to ~2%, established lattice QCD converts scale into hadron mass. No measured hadron mass enters the prediction ancestry |
| Neutrino sector (mixing and masses) | $U_{\mathrm{PMNS}},m_\nu$ | no nonzero source-only physical mass or mixing row emitted: the weighted-cycle mixing candidate is built from a hand-written family template | NuFIT 6.1 used only for audit | the isotropic ansatz is excluded and the target-informed weighted-cycle candidate fails the correlated profile |
| Electromagnetic carrier (photon) | $\mu_{\gamma,\mathrm{hard}}^2$ | zero: the unbroken electromagnetic gauge symmetry of the declared Maxwell action forbids a hard mass term; two classical transverse modes; quantum photon pole not promoted | photon-mass bound <1e-18 eV | classical action gate only |
| Color carrier (gluons) | $\mu_{\mathrm{YM,hard}}^2$ | zero: the unbroken color gauge symmetry of the perturbative pure-Yang--Mills action forbids a hard mass term; no free-particle claim in confined QCD | no isolated free-gluon mass row | deconfinement/quantum pole gate open |
| Gravitational TT carrier (graviton) | $\mu_{\mathrm{FP,hard}}^2$ | zero: the coordinate freedom of the pure-Einstein action forbids a hard mass parameter for the two classical TT modes; quantum graviton pole not constructed | dispersion bound often reported as <1e-32 GeV | classical action gate only |

$\Delta$ reports the sigma distance where PDG or NIST quotes a one-standard-deviation
uncertainty. Otherwise it records the declared support status. A numerical
match counts as a source-only mass prediction only when every upstream object
is source-closed and the physical mass convention has its own certificate. A
conditional implication or an exact adapter remains on its stated audit
surface. The particle paper gives the support status for the electroweak,
charged-lepton, quark, neutrino and hadron derivations.

## Papers

Ordered by importance for a technical reader. Longer summaries mark the
papers that carry the core theorem surface.

- **Paper 2. [Recovering Relativity and the Standard Model from Observer Overlap Consistency](paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf)**: compact technical core for the recovered OPH branch. It states the observer-overlap route to Lorentz structure, the Einstein-like gravity branch on the five-axiom recovered core, the formal algebra core that upgrades rest-frame data to the tensor equation and fixes the metric residue to one $\Lambda$ on connected conserved branches, receipt-conditional compact-gauge reconstruction, the selected Standard Model quotient and matter package with a formal hypercharge and $Z_6$ algebra core, the Borel-Weil local carrier for the one-Higgs slot, Maxwell after an independently supplied electromagnetic action, the separate classical-carrier/quantum-particle gate, and the conditional Yang-Mills mass-gap route under its continuum, reflection-positivity, transfer/intertwiner, and nontriviality assumptions.
- **Paper 1. [Observers Are All You Need](paper/observers_are_all_you_need.pdf)**: broad synthesis and best first read for the full program. It explains finite observer patches, overlap consistency, records, repair moves, the recovered effective universe, the scale story, and the public claim boundaries without replacing the compact paper's theorem ledger.
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle-sector derivation and audit surface. It carries the $P_\star$-driven reconstruction, conditional action-level carrier modes with open quantum-particle gates, the Borel-Weil one-Higgs carrier bridge, electroweak/Higgs/top, quark, charged-lepton, neutrino, and hadron lanes, quantitative checks, and conditional record-worldline stitching. The electroweak headline is the strict source-audit $W/Z$ branch $(80.330,\,91.119)$ GeV and the double-criticality Higgs/top family $(125.77,\,172.63)$ GeV at the frozen boundary-scale candidate; its audit ledger keeps the selected-carrier chart, conditional five-premise value law, and target-anchored declared-surface fit separate, and none supplies a physical complex-pole theorem. The fine-structure audit distinguishes the source/root witness, mixed-provenance no-hadron diagnostic, empirical $e^+e^-\to$ hadrons endpoint, and measured endpoint. The generic quark interface has six common-scale dimensionless Yukawa coordinates, and the reciprocal-ray candidate fails its exact identity at every tested scale. Historic mixed-convention residuals and the completed simulator runs are diagnostics; the simulator measures no quark/Higgs operator, while direct-carrier assays remain Haar-rank-one and miss the required transport shapes. The charged-lepton lane records the exact icosahedral face-corner carrier theorem, the conditional CFQ weight theorem, an engineered finite digital self-reading patch, and a conditional nature/pole transport theorem. The latter assumes the physical Yukawa response and singularity readout it transports; its zero-self-energy kernel does not supply the interacting charged kernel. Physical charged-sector selection, target-free ancestry, cofinal physical refinement, and infrared-complete pole attachment remain open. The isotropic neutrino ansatz is excluded, and the weighted-cycle comparison candidate fails the NuFIT 6.1 correlated profile. No nonzero source-only physical particle mass, physical neutrino mixing matrix, or absolute neutrino mass row is emitted.
- **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**: finite patch-net consensus mechanics. It proves how local observers compare overlap records, apply repair moves, handle defects, and converge to quotient normal forms when the fixed-cutoff assumptions hold. The consensus result stops at quotient normal forms; Lorentzian and Einstein geometry enter through the separate geometric branch in the compact paper. The repair operator is stated on the physical quotient, repaired readouts are invariant under hidden implementation choices, a finite layered carrier witnesses boundary reconstruction, and a finite binary audit fixture gives sharp positive/negative repair and boundary-reconstruction tests. The neutral mathematical companion proves the generic cross-source criterion, but same-boundary physical uniqueness requires the declared boundary map to identify the consistent quotient; same-source confluence and liveness are separate obligations.
- **Paper 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)**: finite carrier and observer-record surface. It gives the echosahedral multi-port patch-carrier architecture, twelve-port screen-sieve theorem, the $A_5/C_3$ face-corner bundle with its explicit physical-carrier boundary, and the engineered fixed-cutoff CFQ central-record model with its dynamical-selection boundary, alongside edge-center scalar-slot completeness, the finite scalar channel bridge, the quotient-edge $Z_6$ reserve, public hardware-evidence rules, records, recovery moves, checkpoint restoration, and observer synchronization.
- **Paper 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**: speculative meaning-layer synthesis outside the recovered physics theorem package. It reads the same OPH machinery as a theory of observer continuation, paradise and hell as continuation environments, resurrection as record-preserving continuation, justice as harm-and-repair bookkeeping, and a strange loop in which observers reverse engineer and build continuation machinery.

## Supplemental Papers And Notes

These support or test the core stack. The most important items get more detail;
lower-level notes are linked with shorter summaries.

- **[Observation-Determined Normal Forms](extra/observable_normal_forms.pdf)**: standalone, substrate-neutral mathematics for constraint and rewrite systems. It separates same-source confluence, cross-source identification from protected observations, normalization and liveness, and local repairability. It adds residual and inverse-observation stability moduli, refinement and inverse-limit bounds, the finite weighted conditional-expectation projector with a noncircular matrix receipt, and a dedicated Lean artifact for the formalized theorem subset.
- **[Observer Patch Holography as a Strange-Loop Self-Simulation: The Claim Lattice and Its Closure Status](extra/compact_proof_of_oph.pdf)**: the compact argument for OPH. It states the strange-loop principles once, walks the recovered-physics theorem legs, presents every number row with its closure status and pull, and names the artifacts that would establish the hypothesis — uniqueness certificates, a closure ledger at zero under blind completions, and a positive compression scorecard.
- **[OPH Falsifiability Map](extra/OPH_falsifiability.md)**: public kill-list for OPH. It names hard failure modes, including failure of a fully specified classical carrier receipt or of a separately completed quantum pole receipt, gauge-mediated proton decay, extra light matter generations, charge-lattice outliers, and the already-excluded universal Solar-System extension of the static dark-response equation. A neutrino or cosmology result can falsify an OPH claim only when it excludes a prediction fixed independently of the test data.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)**: source fixed-point calculation for the fine-structure row. It separates the OPH source value, the low-energy empirical endpoint boundary, the provenance distinction between the source root and the CODATA comparison pixel, and the remaining QCD/hadronic correction.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)**: finite OPH repair-gap mechanism. The Clay-facing statement is a conditional identity: the repair gap equals the four-dimensional Yang-Mills mass gap only if the declared multiresolution continuum, reflection-positivity, transfer/intertwiner, and nontriviality certificates hold, and none of those certificates is supplied.
- **[Observer-Patch Holography as a String-Vacuum Selector](extra/observer_patch_holography_as_string_vacuum_selector.pdf)**: string theory as an effective OPH edge language and vacuum-selection sieve. The Bouchard–Donagi result supplies a visible massless-cohomology candidate; certificates not emitted here include the critical edge, raw cohomology reproduction, safety-layer realization, heavy spectrum, low-energy decoupling, thresholds, and moduli locking.
- **[Photonic Fixed-Point Consensus for SHA-256d Proof of Work](extra/Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf)**: hardware-facing test of OPH-style photonic candidate enrichment for SHA-256d, judged by the exact digital verifier.
- **[Thinking as Patch-Net Fixed-Point Search](extra/thinking_as_patch_net_fixed_point_search.pdf)**: cognition and qualia as recurrent patch consensus across neural or artificial self-reading substrates.
- **[Theoretical Bounds on χν in Observer-Patch Holography](extra/chi_nu_susceptibility_bounds.pdf)**: coherent-matter same-channel theorem, finite channel bridge, protected-reserve collar coefficient $\chi_\nu^{\mathrm{can}}=e^{-P_\chi/24}$, finite-thickness band, conservation-law energy cage, and engineering chart values with device force and cosmological dark stress kept as separate receipt gates.
- **[Entanglement Geometry Problem in OPH](extra/ENTANGLEMENT_GEOMETRY_PROBLEM_OPH.md)**: note on entanglement geometry as an observer-overlap and record-surface problem.
- **[Common Objections](extra/COMMON_OBJECTIONS.md)**: short responses to frequent conceptual and technical objections.
- **[Hacking the Simulation: Anti-Gravity Exploit](extra/hacking-the-simulation-anti-gravity-exploit.pdf)**: pop-science OPH-adjacent engineering book about the local chi-nu lift test. The [Markdown source chapters](extra/hacking-the-simulation-anti-gravity-exploit/) are included with the repo so OPH Sage can ingest the same text during re-ingestion.

## Cosmology Papers

The cosmology branch lives in [`cosmology/`](cosmology/README.md). Its claims
are conditional on OPH-native source, transfer, and likelihood boundaries; FLRW
machinery can serve as comparison plumbing but does not by itself promote an
OPH-native cosmology result.

The round-one public-data audit gives three encouraging conditional checks. The
fixed $\mathbb Z_6$ galaxy response has an aggregate RAR scatter of $0.132834$
dex; the BTFR exponent $4$ is $1.80\sigma$ from an error-aware fit; and the
analytic candidate $n_s=1-P_\star/48=0.9660215$ is $0.27\sigma$ from the Planck
summary, with a conventional 83-bin CAMB diagnostic only $\Delta\chi^2=+0.83$
above its Planck-like baseline. These are retrospective or conditional
comparison receipts, not frozen predictions. The same static dark-response
equation is excluded as a universal Solar-System law by the Cassini quadrupole
summary (about $19.22\sigma$ on the fixed $\mathbb Z_6$ input), so the
old-settled-galaxy branch still needs a source-derived applicability rule. No
cosmology row yet has a frozen physical-prediction receipt.

- **[Observer-Patch Holography and the Dark Matter Phenomenon](cosmology/oph_dark_matter_paper.pdf)**: release-bundle cosmology paper. It treats dark/anomaly stress as imperfect observer-patch repair bookkeeping, imports the quotient-edge scalar, finite channel bridge, and $Z_6$ finite-thickness theorem stack for the local coefficient, gives the galaxy-limit/MOND-like behavior, reports the audited RAR, BTFR, and Cassini applicability checks, defines the source-only anomaly abundance selector, and states the cluster, cosmology, and simulator promotion contracts for larger-scale promotion.
- **[OPH Cosmology as a Finite-Source Prediction Program](cosmology/oph_cosmology_finite_source_cmb_program.pdf)**: CMB-facing program for source-only inputs, scale calibration, Boltzmann transfer, simulator checks, physical CMB boundaries, and claim labels. It reports the conditional analytic $P_\star/48$ tilt comparison and treats source-only dark abundance as a separate source receipt from CMB transfer and likelihood promotion.
- **[Inflation Without an Inflaton](cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.pdf)**: inflation-free branch using observer-screen synchronization, horizon coherence, flatness conditions, geometric screen spectrum, screen release amplitude, radial lift, and hot source data.
- **[OPH Cosmological Vacuum and Structure Formation](cosmology/oph_cosmological_vacuum_and_structure_formation.pdf)**: OPH-native vacuum boundary, fluctuation ensembles, proto-object/worldline formation, and structure-seed checks.
- **[OPH Cosmology Data and Likelihood Contracts](cosmology/oph_cosmology_data_likelihood_contracts.pdf)**: frozen source artifacts, no-data-use receipts, the round-one audited comparison ledger, pooled reducers, Boltzmann-transfer comparisons, and official likelihood protocols.

## Physics Problems Articles

Applied problem notes live in [`physics-problems/`](physics-problems/README.md).
That folder carries the article list, summaries, motivating-result links, claim
boundaries, and OPH Sage ingestion notes. The notes are Markdown-only and stay
outside the paper release, website paper-index, and GitHub release-asset
pipeline.

## Proof Status

The [compact proof of OPH](extra/compact_proof_of_oph.pdf) is the shortest
complete statement of the claim lattice: the strange-loop principles, the
recovered-physics theorem legs, the fixed-point number rows with their
closure status, the compression accounting, and the explicit failure tests.
Its conclusions hold on the declared branch hypotheses, and its evidential
standing is exactly what [CLOSURE_LEDGER.md](CLOSURE_LEDGER.md) and
[COMPRESSION_SCORECARD.md](COMPRESSION_SCORECARD.md) record on the day of
reading: no probability-of-correctness number substitutes for those two
surfaces. The decisive open computation is the blind Ward-projected hadronic
transport against the externally timestamped frozen target in
[`falsification/frozen_targets/`](../falsification/frozen_targets/).

## Applications And OMEGA Hardware

OPH is also a hardware program: explicit screen microphysics turns the same
patch-consensus loop into an engineering handle on reality. A
bounded device exposes boundary data, compares records, repairs mismatch, and
locks onto stable states. OMEGA is the public hardware route into that loop:
physical chambers, labeled ports, control software, verifier records, and
repeatable records.

OPH turns screen microphysics into a way to hack reality. The target
is physical control of small patches that can be driven, measured, repaired,
and verified.

The application thesis uses machines to drive small physical patches toward the
fixed points selected by observer-patch consistency. That gives low-cost
implementation tracks for
desktop fusion energy, room-temperature OMEGA supercomputing, OMEGA-based AGI,
and local gravity or inertia control for hoverbikes and hoverboards. These are
application tracks behind evidence gates; settled-output claims belong to
verifier records and experiments. The compute claim is narrower: a
chamber-conditioned candidate distribution may reduce verifier work by a
measured lift $B=p_Q/p_U$. The classical complexity-class problem remains
untouched.

Read the public applications page at
[omega.floatingpragma.io](https://omega.floatingpragma.io/). Source notes for
the application tracks live in [`APPLICATIONS.md`](APPLICATIONS.md).

## Detailed Diagram And Quantitative Status

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

<p align="center"><sub>The main OPH line from the strange-loop principles to relativity, gauge structure, particles, and observers. The consensus node distinguishes fixed-source schedule independence from the separate boundary-identification gate. Click to open the full SVG.</sub></p>

**Particle derivation stack**

<p align="center">
  <a href="code/particles/particle_mass_derivation_graph.svg" target="_blank" rel="noopener noreferrer">
    <img src="code/particles/particle_mass_derivation_graph.svg" alt="OPH particle derivation stack" width="78%">
  </a>
</p>

<p align="center"><sub>A compact view of the particle lane and its strict claim boundaries. The ledger prints no nonzero source-only physical particle mass; classical zero-mode statements retain their separate quantum-pole gates. Click to open the full SVG.</sub></p>

## More

- **Website:** [floatingpragma.io/oph](https://floatingpragma.io/oph)
- **Theory explainer:** [floatingpragma.io/oph/theory-of-everything](https://floatingpragma.io/oph/theory-of-everything)
- **Coherence map:** [coherence.floatingpragma.io](https://coherence.floatingpragma.io): public graph surface for OPH concepts, overlaps, and cross-domain routes.
- **Simulation:** [simulation.floatingpragma.io](https://simulation.floatingpragma.io): interactive OPH mini-universe explorer showing observer patches, overlap readback, inconsistency repair, records, and emergent geometry.
- **Applications:** [omega.floatingpragma.io](https://omega.floatingpragma.io): public applications page for OPH hardware, compute, energy, AGI, lift, and optical chamber consensus.
- **Blog:** [blog.floatingpragma.io](https://blog.floatingpragma.io/) collects public OPH essays. Start with [Semiotics and the Physics of Meaning](https://blog.floatingpragma.io/semiotics-and-the-physics-of-meaning), [The Trigger](https://blog.floatingpragma.io/the-trigger), and [P = NP on the Observer Screen](https://blog.floatingpragma.io/p-equals-np-on-the-observer-screen). The computation essay treats `P = NP` as an observer-screen slogan; the classical complexity problem remains untouched.
- **Book:** [oph-book.floatingpragma.io](https://oph-book.floatingpragma.io)
- **Guided study app:** [learn.floatingpragma.io](https://learn.floatingpragma.io/)
- **Questions and detailed explanations:** OPH Sage on [Telegram](https://t.me/HoloObserverBot), [X](https://x.com/OphSage), or [Bluesky](https://bsky.app/profile/ophsage.bsky.social)
- **OPH Notebook:** [NotebookLM source notebook](https://notebooklm.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a) with explainer videos and additional study material.
- **Lab:** [oph-lab.floatingpragma.io](https://oph-lab.floatingpragma.io)
- **Common objections:** [extra/COMMON_OBJECTIONS.md](extra/COMMON_OBJECTIONS.md)

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
- **[`physics-problems/`](physics-problems):** standalone markdown physics
  problem notes for public reading and OPH Sage ingestion.

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

OPH is published so the mathematics, software, applications, devices,
hardware designs, simulations, engineering methods, and experimental
implementations can be studied, tested, implemented, modified, deployed,
manufactured, and shared. OPH-derived work may not be used to create private
patent monopolies or patent claims that restrict others from practicing OPH.

See [PATENTS.md](PATENTS.md) for the canonical policy text and copy/paste
website notices.
