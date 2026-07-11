# Fractional Excitons As OPH Quotient-Sector Readouts

## Motivating Result

This note entered the queue after a 2026 photoluminescence experiment in twisted
MoTe$_2$ interpreted emergent peaks as trions bound to fractional charge inside
a fractional Chern-insulator state
([Nature 651, 48–53, 2026](https://www.nature.com/articles/s41586-026-10101-w)).
A second optical report presented evidence for a fragile $\nu=-1/3$ fractional
quantum anomalous Hall state in the same material family
([Phys. Rev. Lett. 136, 056601, 2026](https://doi.org/10.1103/f4dj-7sts)).
The OPH question is how a photon-visible bound state can be attached to a
fractional topological sector without turning every optical line into a phase
identification.

**Status:** conditional fractional-material representation and diagnostic
receipt contract. The note does not establish that a reported optical line is a
fractional exciton, nor that OPH supplies a unique material selector.
Material-specific promotion remains blocked at the frozen Hamiltonian/source
law, global topological-charge accounting, and phase/optical certificates.

Date: 2026-07-08

## Problem, Standard Baseline, and Unresolved Target

Conventional many-body theory supplies projected Chern-band
Hamiltonians, $K$ matrices and braided fusion categories, exciton and polaron
spectral functions, optical selection rules, and exact-diagonalization or
tensor-network diagnostics. The unresolved target is not naming these objects.
It is showing that one frozen material Hamiltonian realizes a particular
fractional phase and that a measured optical feature couples injectively to a
specific allowed fusion channel after binding drift, disorder, and ordinary
excitonic alternatives are controlled.

## Actual OPH Contribution

OPH packages the material as a bounded self-reading patch with local state,
electromagnetic and optical ports, durable transport and spectral records,
repair or continuation maps, readback, and public receipts. A material
presentation $\mathcal F_{x,r}$ defines the finite quotient. A frozen source law
supplies candidate weights that a normal form cannot supply. A topological
ledger records conventional charge, statistics, fusion, and edge data. An
optical module records the allowed sector action of optical operators, and a
line fan is the public readout. This is an evidence architecture, not a claim
that OPH replaces standard topological or spectroscopy calculations.

## Claim Boundary

The result is a conditional specification for a sandbox branch. It says how
OPH records a fractionalized material once the material presentation, source
law, and certificates are supplied. It does not claim that an arbitrary moire
sample, Chern band, fractional Chern insulator, fractional quantum anomalous
Hall device, or optical line belongs to a specified topological phase.

A material-specific proof requires public evidence for the source Hamiltonian,
active band projector, Chern number, band geometry, many-body gap,
ground-sector degeneracy, flux-insertion pump, Hall conductance, edge spectrum,
sector ledger, optical operator, global topological-charge balance, refinement
stability, conventional controls, and no-target-leak audit. Without those
receipts, the branch remains a diagnostic analogy and not a material theorem.

The epistemic ledger is:

- **Assumed:** the physical quotient, repair/canonicalization rules, candidate
  topological ledger, optical module action, finite cutoffs, and frozen source.
- **Derived within the model:** representative-invariant readouts, spectral
  decomposition for a declared optical operator, conditional identifiability,
  and non-selection by canonicalization alone.
- **Empirical or numerical receipts:** phase realization, gap and transport,
  global charge compensation, optical matrix elements and linewidths, lever
  arms, binding drift, ordinary-exciton controls, and refinement stability.

## Fractional Material Presentation

**Definition 1** (Fractional material presentation). *At regulator $r$, a fractional material presentation for a sample or simulation target $x$ is
$$
\mathcal F_{x,r}
=
\left(
\Sigma^{\mathrm{mat}}_{x,r},
\Gamma_{x,r},
Q_{x,r},
\mathcal R_{x,r},
\mathcal U_{x,r},
\mathrm{Chk}_{x,r}
\right).
$$
$\Sigma^{\mathrm{mat}}_{x,r}$ is the finite presentation space: lattice,
geometry, active bands, filling, interactions, disorder or strain fields,
electromagnetic ports, optical ports, edge/collar charts, and finite
Hilbert-space cutoffs. $\Gamma_{x,r}$ quotients gauge representatives, basis
choices, mesh labels, orbital relabelings, edge charts, repair schedules,
hidden carrier coordinates, and inert ancillas only when they preserve every
declared public correlation and response. Physically distinct orbitals,
boundaries, or schedules are not quotient artifacts. The physical quotient is
$$
Q_{x,r}=\Sigma^{\mathrm{mat}}_{x,r}/\Gamma_{x,r}.
$$
$\mathcal R_{x,r}$ is the public record algebra, $\mathcal U_{x,r}$ is the allowed repair/update family, and $\mathrm{Chk}_{x,r}$ is the refinement and receipt data.*

**Theorem 2** (Material quotient normal form). *Suppose the accepted repair
relation on $Q_{x,r}$ is terminating, quotient-descended, locally confluent,
and repair-complete on the declared branch. Then the material presentation has
a schedule-independent normal form
$$
n_{x,r}:Q_{x,r}\to N_{x,r}.
$$
Supply an inclusion $i_{x,r}:N_{x,r}\hookrightarrow Q_{x,r}$ with
$n_{x,r}\circ i_{x,r}=\mathrm{id}_{N_{x,r}}$ before treating $n_{x,r}$ as an
idempotent canonicalizer. To represent $N_{x,r}$ by a bulk topological ledger,
a separately declared edge/collar realization, optical-module data, and public
records, the branch must additionally supply a readout map $D_{x,r}$ that is
constant along accepted repairs and complete in the sense that
$$
D_{x,r}(q)=D_{x,r}(q')
\quad\Longleftrightarrow\quad
n_{x,r}(q)=n_{x,r}(q').
$$
Under this additional premise, $N_{x,r}$ is identified with
$\operatorname{im}D_{x,r}$.*

*Proof.* The finite consensus theorem applies to the physical quotient, not to
hidden representatives. The four stated repair hypotheses give a unique chosen
representative in each repair-connected quotient class. Hidden labels collapse
under $\Gamma_{x,r}$. Constancy of $D_{x,r}$ makes it factor through the normal
form, and the reverse implication in the displayed condition makes the induced
map on $N_{x,r}$ injective. Termination and confluence establish uniqueness of
the chosen representative; completeness of $D_{x,r}$ and existence of the
physical topological phase remain input receipts. $\square$

## Source Law And Non-Selection

The canonicalizer does not choose a material phase. A material branch must declare a source law
$$
\mu_{x,r,T}(q)
=
Z^{-1}_{x,r,T}m_{x,r}(q)\exp[-S_{x,r,T}(q)]
$$
where $S$ is dimensionless. On a finite classical quotient,
$Z=\sum_qm(q)e^{-S(q)}$ must be finite and nonzero; on a continuous
quotient the base measure and integrability must be declared. Alternatively,
use a positive trace-one quantum Gibbs state or a ground-projector source built
from a frozen Hamiltonian with its temperature and boundary conditions stated.
The required tags are
$$
\mathrm{SOURCE\_LAW\_REQUIRED},
\qquad
\mathrm{NORMAL\_FORM\_IS\_NOT\_SELECTOR}.
$$

**Theorem 3** (Normal-form non-selection). *Let $n_{x,r}:Q_{x,r}\to N_{x,r}$ be a material normal-form map, and let $\mathcal X$ be a finite set of candidate topological sectors inside $N_{x,r}$. The map $n_{x,r}$ does not select a unique material sector unless the branch supplies a quotient-intrinsic source law, Hamiltonian, transfer operator, or vacuum certificate.*

*Proof.* With the inclusion $i:N\hookrightarrow Q$ declared above, the typed
self-map on laws is $\mathcal C=(i)_\#(n)_\#$. It satisfies
$\mathcal C^2=\mathcal C$. Every law concentrated on an embedded candidate
sector is fixed by the canonicalizer. Therefore canonicalization can classify
candidate sectors, but it cannot choose their material weights. Selection
requires the source object. $\square$

## Minimal Hamiltonian And Promotion Certificate

A practical fractional-Chern sandbox may use a projected interacting Chern-band Hamiltonian
$$
H^{\mathrm{FCI}}_{x,r}
=
P_C\left[
\sum_{\bm k} \varepsilon_{\bm k} c^\dagger_{\bm k}c_{\bm k}
+
\frac12\sum_{\bm q}V(\bm q)\rho_{-\bm q}\rho_{\bm q}
+
H_{\mathrm{dis/strain/gate}}
\right]P_C ,
$$
with optical extension
$$
H^{\mathrm{opt}}_{x,r}
=H^{\mathrm{FCI}}_{x,r}+H_{\mathrm{hole}}+H_{\mathrm{light}}+H_{\mathrm{bind}}.
$$
The interaction term must declare normal ordering, neutralizing background,
finite geometry, and the treatment of the $\bm q=0$ mode. This is a schema,
not a material proof. The branch must freeze the source before comparison and
emit a promotion certificate
$$
\mathrm{PhaseCert}_{x,r}
=
\left(
C,\nu,\Delta,G,\sigma_{xy},\mathrm{Pump},
\mathrm{Edge},\mathrm{Ent},\mathrm{Anyon},\mathrm{Modular},
\mathrm{Refine}
\right).
$$

Here $C$ is the active-band Chern number, $\nu$ the filling, $\Delta$
the many-body gap, and $G$ the ground-sector degeneracy/manifold receipt;
the remaining fields respectively record electrical Hall response, flux pump,
edge spectrum, entanglement diagnostics, anyon ledger, modular data when
defined, and regulator refinement.

The receipt names are:
$$
\begin{gathered}
\mathrm{SOURCE\_HAMILTONIAN\_FROZEN},\quad
\mathrm{ACTIVE\_BAND\_PROJECTOR},\quad
\mathrm{CHERN\_NUMBER},\\
\mathrm{BAND\_GEOMETRY},\quad
\mathrm{MANYBODY\_GAP},\quad
\mathrm{GROUND\_SECTOR\_DEGENERACY},\\
\mathrm{FLUX\_INSERTION\_PUMP},\quad
\mathrm{HALL\_CONDUCTANCE},\quad
\mathrm{EDGE\_SPECTRUM},\\
\mathrm{TOPOLOGICAL\_SECTOR\_LEDGER},\quad
\mathrm{REFINEMENT\_STABILITY},\quad
\mathrm{NO\_TARGET\_LEAK}.
\end{gathered}
$$

**Proposition 4** (Selection versus identification). *Let $H_{x,r}$ be frozen
before comparison, and let $\mathcal X_{x,r}$ be the candidate sector set. A
source-side calculation selects a unique candidate only if a declared
energy/free-energy or projector criterion has a unique, refinement-stable
winner before target optical data are loaded. Separately, let
$$
\pi_{x,r}:\mathcal X_{x,r}\to \mathcal L_{x,r}
$$
map each candidate to its predicted public topological ledger. If
$\pi_{x,r}$ is injective and the measured certificate agrees with exactly one
frozen prediction within declared uncertainties, the evidence identifies that
candidate. Injectivity does not make the Hamiltonian select it and must not be
manufactured after inspecting the target. If the source-side winner is absent
or the evidence map is non-injective, the correct output is
$\mathrm{SECTOR\_AMBIGUOUS}$.*

## Topological Ledger And Fractional Readout

On an Abelian electronic $U(1)^N$ Chern--Simons branch, the ledger contains a
symmetric, nonsingular integral matrix $K$ and a primitive integral charge
vector $t$. In an electronic basis it must also satisfy the spin/charge parity
condition $K_{II}=t_I\pmod 2$. The physical ledger is the stable equivalence
class $[K,t]$ under integral changes of basis and addition or removal of
explicitly trivial local sectors, not one chosen matrix presentation. Anyons
are classes
$$
[\ell]\in \mathbb Z^N/K\mathbb Z^N .
$$
For a representative $\ell$,
$$
\nu=t^T K^{-1}t,\qquad
Q_\ell/e=t^TK^{-1}\ell,
$$
$$
\theta_\ell=\pi\ell^TK^{-1}\ell,\qquad
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell' .
$$
Charge is defined modulo the charge of local particles, topological spin
modulo $2\pi$, and mutual braiding modulo $2\pi$. These formulas do not apply
unchanged to a non-Abelian phase.

On the non-Abelian branch, the ledger carries
$$
\mathcal C_r=(\mathrm{Irr},N,F,R,d,\theta;S,T),
$$
together with the local electron object, electromagnetic charge grading, and
the consistency relations appropriate to a unitary braided fusion category or
modular extension. The modular matrices $S,T$ are included only when the
category or its declared extension supplies them; a fermionic super-modular
category is not silently treated as modular.

**Proposition 5** (Fractional charge and statistics are quotient-sector
readouts). *For a certified fractional Hall or fractional Chern branch,
fractional charge, exchange phase, mutual braiding, fusion, and the universal
bulk anomaly data are functions of the quotient-sector ledger. A detailed edge
spectrum additionally depends on the boundary realization. None of these data
is a property of a bare local representative in
$\Sigma^{\mathrm{mat}}_{x,r}$.*

*Proof.* Gauge, basis, mesh, orbital, and hidden carrier relabelings can change
local representatives without changing the public sector. The fractional
observables above are invariant under those relabelings and are read through
holonomy, interferometry, transport, a flux pump, modular data, or a separately
certified edge realization. Therefore the universal data factor through
$Q_{x,r}$ and its normal form. $\square$

## Optical Module And Line Fan

**Definition 6** (Optical module). *An optical fractional-exciton ledger is a
module category
$$
\mathcal M_{\mathrm{opt}}^{x,r}
$$
over the topological category $\mathcal C_r$. Each candidate sector $m$ has a
predicted topological shadow $\tau_{\mathrm{pred}}(m)\in\mathrm{Irr}(\mathcal C_r)$,
modeled total electromagnetic charge
$Q_{\mathrm{tot}}(m)$, oscillator data, polarization data, and a binding-energy
term. The module action and optical operator must be derived from the frozen
source or supplied as assumptions; the observed line must not be used to
define them.*

A neutral composite may carry a nontrivial *relative* constituent or fusion
channel,
$$
Q_{\mathrm{tot}}(m)=0,
\qquad \tau_{\mathrm{pred}}^{\mathrm{rel}}(m)\ne 1 .
$$
This is not permission for a strictly local photon operator to create isolated
net topological charge. In a closed sample a local operator preserves the
global superselection sector: relative to a vacuum initial sector, the complete
final state must fuse to a local particle or the vacuum. A nontrivial relative
channel therefore requires a compensating anyon, a pre-existing sector, or an
explicitly modeled boundary or nonlocal process. The receipt must state which
charge-balancing mechanism is present.

**Definition 7** (Optical line fan). *Let $\mathcal Y_{\mathrm{opt}}$ be the
uncertainty-aware space of measured optical records. A minimal record is
$$
y=(\Delta E,I,\partial_g\Delta E,\mathrm{pol},\eta,\ldots)
\in\mathcal Y_{\mathrm{opt}} .
$$
The ellipsis may contain independently calibrated time, momentum, temperature,
or magnetic-field dependence. The topological shadow and total charge are
candidate-ledger annotations, not automatically measured coordinates. The
predicted line fan is therefore
$$
\widehat{\mathcal L}_{\mathrm{opt}}^{x,r}
=\{(R_m,\tau_{\mathrm{pred}}(m),Q_{\mathrm{tot}}(m)):m\in\mathcal X_{x,r}\},
$$
where $R_m\subseteq\mathcal Y_{\mathrm{opt}}$ is a frozen acceptance region
constructed from the candidate's predictive law and declared calibration
uncertainties.*

**Proposition 8** (Line-fan decomposition). *Let $O$ be a declared absorption
operator certified to descend to the quotient. At zero temperature, with $E$
denoting photon energy, its spectral function is
$$
A_O(E,g)
=-\frac{1}{\pi}\operatorname{Im}
\left\langle 0\left|O^\dagger
\frac{1}{E+E_0(g)-H(g)+i\eta}O\right|0\right\rangle
=\sum_m |\langle m|O|0\rangle|^2
\delta_\eta\!\left(E-[E_m(g)-E_0(g)]\right).
$$
Here $\delta_\eta$ is normalized in energy, so it has units of inverse energy
and $A_O$ has squared-matrix-element units per energy. If angular frequency is used instead, every
energy argument is replaced consistently by $\hbar\omega$. Photoluminescence
is not obtained merely by relabeling this absorption function: initial-state
occupations, relaxation, and radiative rates must also be supplied. A public
line fan is a quotient readout only after those operator and population
assumptions are certified.*

**Proposition 9** (Fractional optical slope boundary). *Write the transition
energy for a calibrated electrostatic potential $\phi(g)$ as
$$
\Delta E_m(g)
=-Q_{\mathrm{tot}}(m)\phi(g)
+E_{\mathrm{bind},m}(g)+E_{\mathrm{other},m}(g).
$$
Let $s_m=\partial_g\Delta E_m$, $\alpha_\phi=\partial_g\phi\ne0$, and suppose
the combined nuisance slope is independently bounded by
$$
\left|\partial_gE_{\mathrm{bind},m}
+\partial_gE_{\mathrm{other},m}\right|\le b_m .
$$
Then
$$
\left|Q_{\mathrm{tot}}(m)+\frac{s_m}{\alpha_\phi}\right|
\le \frac{b_m}{|\alpha_\phi|}.
$$
The units of $g$, $\phi$, $s_m$, and the lever arm $\alpha_\phi$ must be
reported. A neutral composite can have a vanishing leading charge slope while
its internal relative fusion channel is nontrivial.*

*Proof.* Differentiating the displayed energy gives
$s_m=-Q_{\mathrm{tot}}\alpha_\phi+\delta_m$, with
$|\delta_m|\le b_m$. Rearrangement gives the interval above. Thus an
uncalibrated lever arm or an unbounded nuisance derivative prevents a charge
inference. $\square$

**Proposition 10** (Optical identifiability). *Let
$$
\mathrm{ID}_{x,r}:\mathcal X_{x,r}\longrightarrow
\mathcal P(\mathcal Y_{\mathrm{opt}})
$$
send each frozen candidate to its declared acceptance region $R_m$ in
measured-observable space. Distinct point predictions are necessary in the
noiseless idealization, but set-map injectivity is not sufficient for a
finite-resolution experiment because different regions can overlap. A
preregistered globally identifying design requires the candidate regions to be
pairwise disjoint at the declared uncertainty level; for an individual record,
it must at minimum belong to exactly one frozen region. Ordinary
exciton/polaron alternatives must be included, and the map must be frozen
before the target line is inspected. Including
$\tau_{\mathrm{pred}}$ itself as an input coordinate would make the argument
circular. Otherwise the correct output is
$\mathrm{OPTICAL\_SECTOR\_AMBIGUOUS}$.*

## Simulator Quotient Correctness

The simulator is evidence only if it implements the declared quotient. It must not define missing objects after the target comparison has been inspected.

**Theorem 11** (Simulator quotient correctness). *Let a simulator store representatives $s\in\Sigma^{\mathrm{mat}}_{x,r}$, a quotient map $q:\Sigma^{\mathrm{mat}}_{x,r}\to Q_{x,r}$, a canonicalizer $c$, observables $O_a$, transition kernels $K$, source law $\mu$, and refinement maps $c_{sr}$. The simulation is quotient-correct for the fractional branch only if it emits:
$$
\begin{gathered}
\mathrm{CANONICALIZER\_IDEMPOTENCE},\quad
\mathrm{REPRESENTATIVE\_INVARIANCE},\quad
\mathrm{QUOTIENT\_LUMPABILITY},\\
\mathrm{DETAILED\_BALANCE\_OR\_DECLARED\_NONEQUILIBRIUM},\quad
\mathrm{REFINEMENT\_COMPATIBILITY},\\
\mathrm{NO\_ORBIT\_SIZE\_BIAS},\quad
\mathrm{NO\_TARGET\_LEAK}.
\end{gathered}
$$
If any required receipt fails, the material conclusion is not promoted.*

*Proof.* Given the declared quotient map and normal-form image,
idempotence verifies that repeated canonicalization leaves the chosen form
fixed; idempotence alone would not construct that quotient. Representative
invariance makes observables functions on $Q_{x,r}$. Lumpability makes a
stochastic transition law descend to the quotient; a deterministic update
instead requires the corresponding quotient-compatibility condition. Source
freezing and the no-target-leak audit prevent posterior fitting. Refinement
compatibility prevents a finite-regulator artifact from being reported as a
stable material sector. Orbit-size control prevents the simulator from
selecting a sector because it has more hidden representatives. $\square$

Implementation receipts exist on both simulator surfaces. The active simulator helpers live in `oph-physics-sim/oph_fractional/`, with the generated sandbox bundle at `oph-physics-sim/runs/fractional/quotient_sector_sandbox/`. The paper-stack mirror lives in `reverse-engineering-reality/code/particles/fractional/`, with the generated bundle at `reverse-engineering-reality/code/particles/runs/fractional/quotient_sector_sandbox/`. Both bundles stop at `FRACTIONAL_QUOTIENT_SANDBOX_DIAGNOSTIC` and block promotion at `MATERIAL_SPECIFIC_HAMILTONIAN_PROOF_RECEIPT`. Any hard-coded passing booleans or synthetic line records in those bundles are scaffold assumptions, not material evidence.

## Failure States

| State | Meaning |
| --- | --- |
| `SOURCE_NOT_FROZEN` | comparison occurred before the material source was fixed |
| `NOT_QUOTIENT_INVARIANT` | an observable changed under hidden representative relabeling |
| `CANONICALIZER_NOT_IDEMPOTENT` | repeated canonicalization changed the normal form |
| `KERNEL_NOT_LUMPABLE` | transition probabilities failed to descend to quotient sectors |
| `ORBIT_SIZE_BIAS_DETECTED` | hidden representative multiplicity affected sector weights |
| `NO_GAP_CERTIFICATE` | the many-body gap was not certified |
| `CHERN_NUMBER_UNSTABLE` | the active-band Chern number failed stability checks |
| `PHASE_CERTIFICATE_NONINJECTIVE` | the material certificate did not select a unique ledger |
| `SECTOR_AMBIGUOUS` | transport or topological data admit more than one sector |
| `OPTICAL_OPERATOR_UNCERTIFIED` | optical operators were not shown to descend to the quotient |
| `GLOBAL_SECTOR_UNBALANCED` | a claimed relative topological channel lacks a compensating global sector |
| `OPTICAL_POPULATION_UNMODELED` | a photoluminescence claim omitted occupations, relaxation, or radiative rates |
| `BINDING_DRIFT_UNBOUNDED` | binding or other nuisance slopes prevent a charge inference |
| `OPTICAL_SECTOR_AMBIGUOUS` | frozen prediction regions overlap or the record matches multiple candidates |
| `TARGET_LEAK_DETECTED` | measured target data entered the source or selector |
| `REFINEMENT_DEFECT_TOO_LARGE` | finite-regulator results failed refinement stability |
| `DIAGNOSTIC_ONLY` | the result is useful for visualization or debugging but not a material claim |

## First Experimental Reading

Twisted transition-metal dichalcogenide bilayers are a natural first sandbox because the same finite material platform can expose moire bands, Chern-band transport, edge or pump signatures, and optical excitonic spectra. The conditional OPH organizing hypothesis is not one fitted peak. It is a preregistered, sector-indexed line fan: the optical peaks should organize by the same quotient-sector ledger that explains charge, edge, and transport. Charged anyon-trion branches can have nonzero slopes after binding drift and other nuisance slopes are bounded. A neutral composite may expose a nontrivial relative fusion channel only with the global compensating sector stated explicitly.

The direct falsifier is equally simple. If preregistered predictions from one
frozen source and one ledger $(\mathcal C,\chi,\Theta)$ cannot jointly account
for the transport, edge, and optical records within declared errors, after
conventional alternatives are included, this proposed OPH material
identification fails.

## Outcome

This note supplies a conditional representation and audit package; it does not
close the material-identification problem. Fractional excitons and related
optical peaks may be represented as quotient-sector readouts when the material
source, topological ledger, global charge balance, optical operator, population
model, line fan, and identifiability map are independently certified on the
same finite quotient. The simulation's job is receipt production: gaps, Chern
numbers, sector data, optical spectra, refinement defects, no-target-leak
audits, quotient lumpability, and identifiability. It is not a substitute for
the material proof.
