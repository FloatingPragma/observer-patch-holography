# Observer Patch Holography (OPH)

> Observer Patch Holography is the observer-consistency theory of everything. No observer sees the whole world at once; each observer gets a local patch; physics is the public fixed point that survives agreement across overlaps.

**French version:** [README_FR.md](README_FR.md)

**Quick links:** [OPH website](https://floatingpragma.io/oph/) | [OPH Textbooks](https://learn.floatingpragma.io/) | [Book: Reverse Engineering Reality](https://oph-book.floatingpragma.io/) | [Ω](https://omega.floatingpragma.io/) | [Blog](https://blog.floatingpragma.io/) | [Tech](https://omega.floatingpragma.io/) | [Simulation](https://simulation.floatingpragma.io)

**Falsifiability:** The [OPH falsifiability map](extra/OPH_falsifiability.md)
lists 40 empirical outcomes that rule out specific OPH claims, including a
failure of the fully specified pure-Einstein transverse-traceless carrier
receipt (or of a separately completed quantum pole receipt), gauge-mediated
proton decay, a fourth light matter generation, and a charge-lattice outlier. The derivation of the physical neutrino mixing
matrix and absolute neutrino masses is work in progress.

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

OPH is formulated as a zero-input theory. The term concerns numerical
provenance: no measured target or fitted numerical constant may enter a
declared source map. The quantum-algebraic axioms and the $S^2$ screen remain
explicit premises. Quantitatively, the public rows are organized by three
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
| Gravitational constant | G | 6.6742999959e-11: emerges from the geometric area normalization $G_{\mathrm{geom}}=\ell_\star^2$ once the caesium clock rule sets the SI display scale | 6.67430(15)e-11 | 0.00003σ |
| Speed of light | c | structural Lorentz speed of the recovered spacetime branch; the SI number is a unit convention | 299792458 exact by definition | not a numeric prediction |
| Fine-structure (inv) | $\alpha^{-1}$ | emerges from the pixel fixed point $P_\star$ through the unified-coupling map: source/root witness 136.99483516462165; separate mixed-provenance no-hadron diagnostic 137.03595950081728; empirical hadron-closure endpoint 136.26365894312011 on [136.15838328342668, 136.36909752396569] | measured endpoint 137.035999177(21) | these are four distinct coordinates. The strict source branch needs the Ward-projected self-map interval proof. The empirical endpoint uses measured $e^+e^-\to$ hadrons data and leaves the same-scheme anchor gap [0.648554111145163, 0.854792066592563] |
| Massive electroweak and Higgs sector (W, Z, Higgs, top) | $W,Z,H,t$ | conditional derivation: the pixel-fixed couplings and the declared electroweak surface fix the branch, and the mass chart then computes the pair with no measured mass entering. Target-free value law $(m_W,m_Z)=(80.3770000154,\,91.1879780779)$ GeV; declared-surface $(m_H,m_t)=(125.1995304097,\,172.3523553288)$ GeV; no source-only physical mass emitted | $m_W=80.3692(133)$, $m_Z=91.1880(20)$, $m_H=125.13(11)$, $m_t=172.1(6)$ GeV | the value-law pair is an exact implication of five finite quotient-transport premises whose carrier certificate is unemitted; the Higgs/top pair holds on the declared running, matching, and threshold surface whose core couplings and Jacobian slopes are declared calibration inputs; the physical scale, scheme, branch, provenance, uncertainty, and complex-pole gates are open |
| Charged leptons (electron, muon, tau) | $e,\mu,\tau$ | no nonzero source-only physical mass emitted: the icosahedral face-corner carrier fixes the family structure, but the map from that structure to absolute masses stays target-anchored | comparison coordinates kept outside this ledger | the exact family witness is target-anchored and lives in the audit surfaces; physical source selection is open |
| Quark spectrum (up, down, strange, charm, bottom, top) | $u,d,s,c,b,t$ | no nonzero source-only physical mass emitted: the source equations produce two ordered Yukawa profile rays but leave their spreads free | common-scale Yukawa data used only for audit | the source equations leave a proven two-modulus non-identifiability; a source-derived flavor-orbit selector is the missing object |
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
- **Paper 3. [Deriving the Particle Zoo from Observer Consistency](paper/deriving_the_particle_zoo_from_observer_consistency.pdf)**: particle-sector derivation and audit surface. It carries the $P_\star$-driven reconstruction, conditional action-level carrier modes with open quantum-particle gates, the Borel-Weil one-Higgs carrier bridge, electroweak/Higgs/top, quark, charged-lepton, neutrino, and hadron lanes, quantitative checks, and conditional record-worldline stitching. Its audit ledger keeps the selected-carrier chart, conditional five-premise value law, and target-containing freeze-once inverse adapter separate; none supplies a physical complex-pole theorem. The fine-structure audit distinguishes the source/root witness, mixed-provenance no-hadron diagnostic, empirical $e^+e^-\to$ hadrons endpoint, and measured endpoint. The generic quark interface has six common-scale dimensionless Yukawa coordinates, and the reciprocal-ray candidate fails its exact identity at every tested scale. Historic mixed-convention residuals and the completed simulator runs are diagnostics; the simulator measures no quark/Higgs operator, while direct-carrier assays remain Haar-rank-one and miss the required transport shapes. The charged-lepton lane records the exact icosahedral face-corner carrier theorem, the conditional CFQ weight theorem, an engineered finite digital self-reading patch, and a conditional nature/pole transport theorem. The latter assumes the physical Yukawa response and singularity readout it transports; its zero-self-energy kernel does not supply the interacting charged kernel. Physical charged-sector selection, target-free ancestry, cofinal physical refinement, and infrared-complete pole attachment remain open. The isotropic neutrino ansatz is excluded, and the weighted-cycle comparison candidate fails the NuFIT 6.1 correlated profile. No nonzero source-only physical particle mass, physical neutrino mixing matrix, or absolute neutrino mass row is emitted.
- **Paper 4. [Reality as a Consensus Protocol](paper/reality_as_consensus_protocol.pdf)**: finite patch-net consensus mechanics. It proves how local observers compare overlap records, apply repair moves, handle defects, and converge to quotient normal forms when the fixed-cutoff assumptions hold. The consensus result stops at quotient normal forms; Lorentzian and Einstein geometry enter through the separate geometric branch in the compact paper. The repair operator is stated on the physical quotient, repaired readouts are invariant under hidden implementation choices, a finite layered carrier witnesses boundary reconstruction, and a finite binary audit fixture gives sharp positive/negative repair and boundary-reconstruction tests. The neutral mathematical companion proves the generic cross-source criterion, but same-boundary physical uniqueness requires the declared boundary map to identify the consistent quotient; same-source confluence and liveness are separate obligations.
- **Paper 5. [Federated Echosahedral Screen Microphysics](paper/screen_microphysics_and_observer_synchronization.pdf)**: finite carrier and observer-record surface. It gives the echosahedral multi-port patch-carrier architecture, twelve-port screen-sieve theorem, the $A_5/C_3$ face-corner bundle with its explicit physical-carrier boundary, and the engineered fixed-cutoff CFQ central-record model with its dynamical-selection boundary, alongside edge-center scalar-slot completeness, the finite scalar channel bridge, the quotient-edge $Z_6$ reserve, public hardware-evidence rules, records, recovery moves, checkpoint restoration, and observer synchronization.
- **Paper 6. [Paradise as Fixed-Point Consensus](paper/paradise_as_fixed_point_consensus.pdf)**: speculative meaning-layer synthesis outside the recovered physics theorem package. It reads the same OPH machinery as a theory of observer continuation, paradise and hell as continuation environments, resurrection as record-preserving continuation, justice as harm-and-repair bookkeeping, and a strange loop in which observers reverse engineer and build continuation machinery.

## Supplemental Papers And Notes

These support or test the core stack. The most important items get more detail;
lower-level notes are linked with shorter summaries.

- **[Observation-Determined Normal Forms](extra/observable_normal_forms.pdf)**: standalone, substrate-neutral mathematics for constraint and rewrite systems. It separates same-source confluence, cross-source identification from protected observations, normalization and liveness, and local repairability. It adds residual and inverse-observation stability moduli, refinement and inverse-limit bounds, the finite weighted conditional-expectation projector with a noncircular matrix receipt, and a dedicated Lean artifact for the formalized theorem subset.
- **[Compact Proof That We Most Likely Inhabit an OPH Simulation](extra/compact_proof_of_oph.pdf)**: shortest compression-style argument for OPH. It collects the five-axiom route, the fixed-branch outputs, the failure points, and the reason numerical agreement only matters when target leakage is excluded.
- **[OPH Falsifiability Map](extra/OPH_falsifiability.md)**: public kill-list for OPH. It names 40 hard failure modes, including failure of a fully specified classical carrier receipt or of a separately completed quantum pole receipt, gauge-mediated proton decay, extra light matter generations, and charge-lattice outliers. A neutrino result can falsify an OPH claim only when it excludes a prediction fixed independently of the test data.
- **[The Fine-Structure Constant as an OPH Pixel Fixed Point](extra/fine_structure_constant_derivation.pdf)**: source fixed-point calculation for the fine-structure row. It separates the OPH source value, the low-energy empirical endpoint boundary, the provenance distinction between the source root and the CODATA comparison pixel, and the remaining QCD/hadronic correction.
- **[Explaining the Yang-Mills Mass Gap with Observer-Patch Repair Dynamics](extra/yang_mills_gap_clay_problem.pdf)**: finite OPH repair-gap mechanism and conditional Clay-facing route. The four-dimensional Yang-Mills identification requires the stated continuum and transfer certificate.
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

- **[Observer-Patch Holography and the Dark Matter Phenomenon](cosmology/oph_dark_matter_paper.pdf)**: release-bundle cosmology paper. It treats dark/anomaly stress as imperfect observer-patch repair bookkeeping, imports the quotient-edge scalar, finite channel bridge, and $Z_6$ finite-thickness theorem stack for the local coefficient, gives the galaxy-limit/MOND-like behavior, defines the source-only anomaly abundance selector, and states the cluster, cosmology, and simulator promotion contracts for larger-scale promotion.
- **[OPH Cosmology as a Finite-Source Prediction Program](cosmology/oph_cosmology_finite_source_cmb_program.pdf)**: CMB-facing program for source-only inputs, scale calibration, Boltzmann transfer, simulator checks, physical CMB boundaries, and claim labels. It treats source-only dark abundance as a separate source receipt from CMB transfer and likelihood promotion.
- **[Inflation Without an Inflaton](cosmology/oph_inflation_without_inflaton_observer_screen_synchronization.pdf)**: inflation-free branch using observer-screen synchronization, horizon coherence, flatness conditions, geometric screen spectrum, screen release amplitude, radial lift, and hot source data.
- **[OPH Cosmological Vacuum and Structure Formation](cosmology/oph_cosmological_vacuum_and_structure_formation.pdf)**: OPH-native vacuum boundary, fluctuation ensembles, proto-object/worldline formation, and structure-seed checks.
- **[OPH Cosmology Data and Likelihood Contracts](cosmology/oph_cosmology_data_likelihood_contracts.pdf)**: frozen source artifacts, no-data-use receipts, pooled reducers, Boltzmann-transfer comparisons, and official likelihood protocols.

## Physics Problems Articles

Applied problem notes live in [`physics-problems/`](physics-problems/README.md).
That folder carries the article list, summaries, motivating-result links, claim
boundaries, and OPH Sage ingestion notes. The notes are Markdown-only and stay
outside the paper release, website paper-index, and GitHub release-asset
pipeline.

## Proof Status

The [compact OPH proof](extra/compact_proof_of_oph.pdf) states a conditional
closure theorem and a protocol for an empirical compression test. Its five
axioms, branch rules, source maps, and readout maps must be fixed independently
of the test data. If the two source maps have unique fixed points, no rival
source map survives the declared rules, and the finite packet passes the
five-part simulation test, the declared observables are uniquely generated up
to observer-invisible equivalence. Physical identification requires comparison
with data. The theorem establishes uniqueness on the declared branch; its
premises carry their stated proof obligations.

The compact count includes mathematical rows for the finite observer-patch
normal form, the quantum event surface, the Standard Model structural branch,
the unified-coupling interval, and geometric gravity normalization. The
balanced-circulant charged-family result is a frozen retrodictive continuation
benchmark; its MaxEnt balance, phase, and absolute normalization are not
source-derived. The QCD-free weak hierarchy is a separate retrodictive
source-side match. The machine-checked artifact covers a finite consensus and
repair subset; the OPH-specific consensus theorem is an open formalization
obligation. Low-energy fine structure and the SI value of Newton's constant are
comparison rows. Particle, neutrino, dark-sector, cosmology, collider, and
hardware results sit outside the compact count with their own evidence classes.

The numerical coincidence budget is an illustration. Its deliberately generous
windows multiply to roughly $10^{-9}$, about 30 bits, and roughly $10^{-3}$
after the windows are widened one hundredfold. These figures are not a posterior
probability. The window sizes are judgment calls, their independence is
approximate, and discarded branch searches must be charged. The formal
empirical claim requires a frozen theory, source maps, branch selectors,
tolerances, uncertainty model, comparison theories, a held-out dataset, and a
declared compression threshold. The compact count admits a result only when its
own calculation and checks satisfy that protocol.

The conditional geometry statements require finite cap data with evidence for
order, orientation, modular behavior, and thermal normalization. Under those
conditions the data converge to the expected geometric cap flow. The compact
paper settles the production question in both directions. A producer theorem
computes an incidence complex directly from the repaired normal forms and,
when that complex passes four decidable receipts (a closed orientable surface
with the right Euler count, a nondegenerate cap mesh, Cauchy modular
cross-ratios, and an independently normalized thermal comparison), it builds
the round sphere, its conformal structure, and the whole certificate.
Explicit countermodels (a torus, a three-sphere skeleton, a wedge of two
spheres, a wrong-temperature clock) prove that bare consensus never selects
the topology, the dimension, or the normalization on its own. The receipts
are the physical input, and the geometry is derived exactly on the receipt
branch.

The same style of bookkeeping runs through the rest of the gravity chain.
The null-net packet proves the standardness and positive-translation
structure the null bridge consumes, with a worked counterexample showing
that finite-range Gibbs locality does not imply modular locality. The
collar-recovery packet likewise separates exact and conditional routes. The
declared central-interface branch has zero collar CMI. Off that branch, a
finite-range Gibbs theorem gives

```math
I(A_\delta:D_\delta\mid B_\delta)
\le c|\partial C|_{\mathrm{UV}}e^{-\delta/\xi}
```

only under uniform strong conditional matrix mixing; ordinary two-point
clustering is insufficient. The continuum schedule must satisfy
$\delta/\xi-\log|\partial C|_{\mathrm{UV}}\to+\infty$.
The weaker ratio $\delta/\ell_{\mathrm{UV}}\to\infty$ can lose to the
growing boundary count. The finite evaluator checks the complete log
envelope and recovery error, but remains a branch-instantiation proxy until
the Gibbs and mixing premises are certified on a realized tower. Scalar CMI
is never treated as a stress tensor or dark-sector source.
The event-manifold packet defines events as coincidence classes of localized
records. On its named receipts (dense population, certified separation,
rank-four response frames, consistent chart overlaps) the event set is a
genuine four-dimensional Lorentzian manifold of signature (-+++), with the
hyperbolic frame space kept strictly as the space of rest frames over each
event, never as space itself. The stress tensor is built from modular
charges by null tomography, the entropy first law carries its edge term
exactly, one uniform scaling family controls the small-diamond limit, a
vacuum-reference receipt fixes the absolute baseline, and the whole chain
composes into a single branch-entry theorem with the Einstein equation as
its output.

The [simulator side](https://github.com/muellerberndt/oph-physics-sim) carries a receipt evaluator for all of this, and the
machine state reads as follows. The collar-CMI artifact is an analytic finite
proxy with no empirical-evidence or theorem-promotion flag; it does not
certify strong conditional mixing on the realized repair tower. A genuine cyclic repair run (real
conflicts, transactional repair, confluence and schedule independence
verified at runtime on a three-stage tower) passes the sphere, mesh, and
naturality receipts on its own repaired output; the overlap net of that run
is chosen spherical, which is explicit branch selection, never a derived
result. Free-fermion collar states pass the modular cross-ratio and thermal
normalization receipts, with the interior modular profile matching the
conformal prediction to about one part in ten thousand at the finest stage.
The null-net receipts pass at one-particle level, including the half-sided
compression condition and a percent-level match of the modular commutator
against the conformal structure constants. Records taken from the actual
repair dynamics reproduce a Lorentzian causal cone: signature (1,2) on the
screen sheet, and signature (1,3) at every tested seed once records carry
the produced depth coordinate on the corpus-consistent multi-scale
dynamics. A strong-coupling variant of that dynamics measures (2,2) and is
kept as a countermodel, which is what a sharp instrument looks like. The
cone margin is small, its convergence with commit density is work in
progress, and the limit clauses, cap-interior data, and physical
identification receipts (universal coupling, vacuum reference, absolute
scale) are open. The realized geometric branch is not certified nonempty,
and the record localization and event receipts must pass on physical data
before any of this counts as more than structure.

Screen-spectrum and CMB continuations are provisional without the screen
branch's geometric scale, source dynamics, clock, refinement behavior, and
observational readout from OPH-native records. The galactic dark-sector row is
a pre-likelihood scorecard outside the compact count. Anomaly, vacuum, and
quantum-foam views are diagnostic. Physical likelihood claims require a
quotient-derived ensemble, regulator-stable reconstruction, source-only
abundance where applicable, and a frozen validation target before likelihood
data are read.

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

<p align="center"><sub>The main OPH line from axioms to relativity, gauge structure, particles, and observers. The consensus node distinguishes fixed-source schedule independence from the separate boundary-identification gate. Click to open the full SVG.</sub></p>

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
