# Fractional Quantum Hall States

## Motivating Result

This note entered the queue after the zero-field fractional quantum anomalous
Hall signal in twisted MoTe2, where a moire Chern band showed fractional Hall
structure without an external magnetic field
([Cai et al., Nature 622, 63-68, 2023](https://pubmed.ncbi.nlm.nih.gov/37315640/)).
The result was surprising because fractional Hall phenomenology had left the
Landau-level setting and entered a material-engineered Chern band. The OPH
question is which boundary-visible records identify the fractional sector, and
which material details merely choose a representative.

**Status:** conditional OPH reformulation and receipt specification. The note
recovers standard $K$-matrix and fusion-category readouts, proves limited
non-selection statements, and separates high-field FQH from zero-field
FCI/FQAH branches. It does not select a material $5/2$ order or derive a
fractional phase from an arbitrary Hamiltonian.

Date: 2026-07-08

## Introduction

The problem has three distinct layers: classify a gapped charged topological
phase, infer its public transport and defect data, and determine which phase a
specific Hamiltonian or sample realizes. A plateau alone need not distinguish
candidate orders. At $\nu=5/2$, filling and minimal charge do not by themselves
imply non-Abelian order; a noncentral sector is required only if a
non-invertible fusion package is independently established.

## Standard-Theory Baseline and Unresolved Target

Conventional fractional Hall physics has powerful local and
substrate-neutral ingredients: Landau-level
projection, trial wavefunctions, composite-fermion constructions, Chern-Simons
effective theories, $K$-matrix stable equivalence, spin topological field
theories, braided fusion categories, edge conformal field theories, and
numerical spectra. These methods classify broad Abelian and non-Abelian
branches. The unresolved target is material selection and identifiability:
several physically different phases can agree on a coarse transport plateau,
and finite width, Landau-level mixing, disorder, and edge reconstruction affect
which candidate is realized and which measurements distinguish it.

The $\nu=5/2$ sector is the sharp failure mode. Pfaffian, anti-Pfaffian,
PH-Pfaffian, Abelian, and reconstructed candidates can share filling fraction
and partial charge data while disagreeing about fusion, chirality, thermal
transport, interferometry, and edge reconstruction. A plateau alone therefore
does not select the observer-facing sector. Legacy physics can name candidate
phases and compute their invariants. No OPH argument here supersedes that
machinery or makes the many-body selection calculation easier.

In that framing, the classification target is underdetermined: the same public
plateau can sit over more than one observer-distinct sector.

## Actual OPH Contribution

OPH repackages the system as a bounded self-reading patch: it has local
many-body state, electromagnetic and edge ports, durable transport records,
repair or continuation maps, readback, and public receipts. The hidden bulk presentation enters
through the declared source law, then passes through quotienting. The physical
state is the quotient of droplet, flux, active Landau level, source Hamiltonian,
edge algebra, records, repair instruments, and checkpoint data. Gauge
representatives, basis choices, mesh labels, electron labels, edge charts, port
labels, repair schedules, hidden carrier coordinates, and ancillas are removed
only when they demonstrably preserve visible interface statistics, record
readout, repair dynamics, and checkpoint continuation.

That quotient supplies an evidence ledger around the conventional normal-form
classification. Abelian
plateaus are identified by an integral edge/holonomy ledger. Non-Abelian
branches require defect and higher-defect data. The $5/2$ branch becomes
identifiable only when the public record includes enough data to distinguish the
noncentral repair sector: charge, heat, tunneling, interferometry, fusion,
braiding, chirality, and refinement-stable gap information. The OPH step that
matters is the observer-facing equivalence relation. Hidden-label differences
collapse. Boundary-visible differences in fusion or transport receipts remain
physical.

The epistemic ledger is:

- **Assumed:** a gapped branch, a declared physical equivalence relation, an
  accepted repair/canonicalization map, and a conventional topological-order
  description on the chosen branch.
- **Derived within the model:** coordinate-invariant readouts from a valid
  $(K,t)$ presentation, the conditional non-Abelian implication of
  non-invertible fusion, and normal-form non-selection.
- **Empirical or numerical receipts:** the source Hamiltonian, gap, Hall and
  thermal Hall response, charge, fusion/braiding evidence, edge realization,
  refinement stability, and no-target-leak audit.

## Abstract

Fractional Hall systems become Hall-collar quotients in OPH. Gap-preserving
local continuation stays inside a topological phase; the public residue lives
in bulk topological data, charge holonomy, and a separately declared edge
realization. Abelian plateaus recover the standard integral ledger. At
$\nu=5/2$, an observed non-invertible fusion package would require a noncentral
repair sector; filling and $e/4$ charge alone remain compatible with multiple
Abelian and non-Abelian candidates.

## Source Branch

The source branch is

$$
\mathsf{OPH}_{\mathrm{FQH}}
=
\mathsf{ObserverPatch}
+
\mathsf{ConsensusHolonomy}
+
\mathsf{ScreenMicrophysics}
+
\mathsf{EdgeCasimir}
+
\mathsf{FiniteQuotientEnsemble}.
$$

The fixed-cutoff observer patch is

$$
\mathsf O_i=
\left(
\mathcal A_i,\rho_i,\mathcal R_i,
(\mathcal I_e,\pi_{i,e})_{e\ni i},
\mathcal U_i,\mathrm{Chk}_i
\right),
$$

where $\mathcal A_i$ is the accessible algebra, $\rho_i$ the state,
$\mathcal R_i$ the record algebra, $\mathcal I_e$ the overlap-visible
interface algebra, $\pi_{i,e}$ the visible restriction, $\mathcal U_i$ the
allowed repair instruments, and $\mathrm{Chk}_i$ checkpoint data.

The Hall droplet is a support/collar chart. Hidden carrier coordinates, edge
charts, port labels, and mesh labels are silent when they preserve visible
interface statistics, record readout, repair instruments, and checkpoint
continuation.

## Branch Assumptions: Continuum FQH Versus FCI/FQAH

Both branches assume a finite regulated two-dimensional charged system with
conserved electromagnetic $U(1)$, a declared many-body or mobility gap,
observer-visible boundary/collar records, and quotient-visible charge, heat,
tunneling, and interferometry readouts. Their source presentations differ.

The continuum Landau-level FQH branch is

$$
\mathrm{FQH}_r =
\left(
2+1\mathrm{D},
U(1)_{\rm em},
B_{\rm ext}\ne0,
N_{\phi,r},
P_{\mathrm{LL},r},
\Delta_{\rm bulk}>0,
\mathcal A_{\partial,r},
\mathcal R_{\partial,r}
\right).
$$

The zero-field fractional-Chern or fractional-quantum-anomalous-Hall branch
motivated by twisted MoTe$_2$ is instead

$$
\mathrm{FCI/FQAH}_r =
\left(
2+1\mathrm{D},
U(1)_{\rm em},
B_{\rm ext}=0,
P_{C,r},
C_r\ne0,
\nu_{\rm cell,r},
\Delta_{\rm bulk}>0,
\mathcal A_{\partial,r},
\mathcal R_{\partial,r}
\right),
$$

with a Chern-band projector, moire-cell filling, and a declared explicit or
spontaneous time-reversal-breaking source and record. $N_\phi$ and a
Landau-level projector must not be inserted into this zero-field branch as if
they were material inputs. The results below classify a supplied gapped
topological ledger; they do not prove that an arbitrary Hamiltonian enters
either phase.

## Hall-Collar Quotient

For the continuum FQH branch, let the finite presentation space
$\Sigma_r^{\mathrm{FQH}}$ consist of tuples

$$
s=
\left(
\mathcal D_r,
\partial\mathcal D_r,
B_r,
N_{e,r},
N_{\phi,r},
\mathcal H_r,
P_{\mathrm{LL},r},
H_r^{\mathrm{eff}},
\mathcal A_{\partial,r},
\mathcal R_{\partial,r},
\mathcal U_{\partial,r},
\mathrm{Chk}_r
\right).
$$

Here $\mathcal D_r$ is the finite droplet, annulus, cylinder, sphere, or torus
geometry; $\partial\mathcal D_r$ is the Hall collar; $B_r$ is flux data;
$N_{e,r}$ and $N_{\phi,r}$ are electron and flux counts; $\mathcal H_r$
is the finite many-body Hilbert space; $P_{\mathrm{LL},r}$ is the active
Landau-level projection; $H_r^{\mathrm{eff}}$ is the effective source
Hamiltonian or action; $\mathcal A_{\partial,r}$ is the collar algebra;
$\mathcal R_{\partial,r}$ is the record algebra; $\mathcal U_{\partial,r}$
is the allowed repair/update family; and $\mathrm{Chk}_r$ is
refinement/checkpoint data.

Let $\Gamma_r^{\mathrm{FQH}}$ quotient away gauge representatives, basis
choices, mesh labels, electron labels, edge charts, port labels, repair
schedules, hidden carrier coordinates, and inert ancillas. The physical
quotient is

$$
Q_r^{\mathrm{FQH}}
=
\Sigma_r^{\mathrm{FQH}}/\Gamma_r^{\mathrm{FQH}}.
$$

For the FCI/FQAH branch, replace $(B_r,N_{\phi,r},P_{\mathrm{LL},r})$ by the
lattice or moire geometry, active Chern-band projector $P_{C,r}$, unit-cell
count and filling, and declared time-reversal-breaking source/record. Denote
the resulting quotient by $Q_r^{\mathrm{FCI}}$. Later uses of
$Q_r^{\mathrm{Hall}}$ mean the disjoint declared union of these branch-specific
quotients, not a single presentation containing incompatible inputs.

For the $\nu=5/2$ branch,

$$
\Sigma_r^{5/2}\subset \Sigma_r^{\mathrm{FQH}},
\qquad
\nu_r=\frac{N_{e,r}}{N_{\phi,r}}\to \frac52.
$$

The arrow is a thermodynamic statement; finite spherical systems also carry a
topological shift, so $N_\phi=\nu^{-1}N_e-\mathcal S$ rather than an exact
finite-size ratio.

After subtracting two filled lower Landau levels,

$$
\nu_{\rm active,r}\to \frac12,
\qquad
Q_r^{5/2}=\Sigma_r^{5/2}/\Gamma_r^{5/2}.
$$

## Hall Normal Form

The accepted repair map is

$$
n_r^{\mathrm{Hall}}:
Q_r^{\mathrm{Hall}}\to N_r^{\mathrm{Hall}}.
$$

Here $N_r^{\mathrm{Hall}}$ is either embedded in $Q_r^{\mathrm{Hall}}$ as a
chosen set of representatives or supplied with an inclusion
$i_r:N_r^{\mathrm{Hall}}\hookrightarrow Q_r^{\mathrm{Hall}}$ satisfying
$n_r\circ i_r=\mathrm{id}_{N_r}$. This typing is required before calling $n_r$
idempotent or applying its pushforward twice.

A Hall normal form is

$$
n_r^{\mathrm{Hall}}(q)
=
\left(
\mathcal E_{\partial,r},
\mathcal H_{\mathrm{hol},r},
\mathcal C_r,
\mathfrak e_r,
\chi_r,
\Theta_r
\right).
$$

$\mathcal E_{\partial,r}$ is a declared edge/collar realization, not by itself
a bulk invariant: distinct reconstructed edges can bound the same bulk order.
$\mathcal H_{\mathrm{hol},r}$ is the Abelian holonomy lattice or
non-Abelian higher-holonomy datum. $\mathcal C_r$ is the defect/fusion
category. $\mathfrak e_r$ is the local electron object. The charge grading is

$$
\chi_r:\mathrm{Obj}(\mathcal C_r)\to \mathbb R/\mathbb Z.
$$

$\Theta_r$ is the spin, braid, and modular-data package.

The readout map is

$$
\mathrm{HallRead}_r:
N_r^{\mathrm{Hall}}\to\mathcal Y_r^{\mathrm{Hall}},
$$

$$
\mathrm{HallRead}_r(n)
=
\left(
\nu,
\sigma_{xy},
e^*,
c_-,
\mathcal A_{\mathrm{anyon}},
\mathcal F,
\mathcal R,
\mathcal B,
\mathcal S,
\mathcal I_{\mathrm{int}}
\right).
$$

The components are filling fraction, electrical Hall conductance, minimal
quasiparticle charge, chiral central charge, anyon set, fusion rules, exchange
data, braiding data, modular/geometric response, and interferometry response.

## Proposition 1: Conditional Hall Normal Form

**Statement.** Let $Q_r^{\mathrm{Hall}}$ be a branch-specific finite physical
quotient of a gapped $2+1$ dimensional charged Hall or Chern system with
conserved $U(1)_{\rm em}$ and an observer-visible boundary/collar algebra.
Assume the accepted OPH repair law on
this branch is terminating, quotient-descended, locally confluent, and
repair-complete. Then every stable Hall phase has a schedule-independent normal
form

$$
n_r^{\mathrm{Hall}}(q)
=
\left(
\text{edge algebra},
\text{charge holonomy},
\text{defect/higher-defect data},
\text{record readout}
\right).
$$

If the residual defect datum is Abelian, it is recorded by an integral
holonomy lattice together with its charge and quadratic/braiding pairing. If it
is genuinely noncentral, the ledger must retain full braided fusion-category
data; a generic reduction to a crossed module is not assumed.

**Proof.** OPH physical data factor through accessible algebras, records,
visible interface restrictions, repair instruments, and checkpoint
continuation. Hidden presentation data are quotient artifacts.

The consensus theorem gives schedule-independent quotient normal forms when the
repair relation is terminating, quotient-descended, locally confluent, and
repair-complete. Pairwise overlap agreement need not guarantee global
agreement; on the Abelian branch the residual can be holonomic. A general
non-Abelian Hall branch is recorded by its braided fusion data. No equivalence
between arbitrary such data and a crossed-module Cech obstruction is assumed.

The additional physical premise is that gap-preserving local continuation stays
inside the same bulk topological phase. That premise comes from the conventional
stability theory of topological order; it is not proved by confluence alone.
The bulk order, a compatible edge realization, charge holonomy,
defect/higher-defect data, and observer-facing readout then form the declared
normal form. $\square$

## Theorem 2: Abelian $U(1)^N$ Chern-Simons Readout

On the Abelian branch, let $\Lambda\simeq\mathbb Z^N$ be the stable defect
lattice, let

$$
K:\Lambda\times\Lambda\to\mathbb Z
$$

be a nondegenerate symmetric integral pairing, and let $t\in\Lambda^\ast$ be
the electromagnetic charge vector. On an electronic branch require a primitive
unit-charge vector and the fermionic parity condition

$$
n^TKn\equiv t^Tn\pmod 2
\qquad(n\in\mathbb Z^N).
$$

**Statement.** On the Abelian branch,

$$
n_r^{\mathrm{Hall}}(q)=[K,t]
$$

where the brackets identify

$$
(K,t)\sim(W^TKW,W^Tt),
\qquad W\in GL(N,\mathbb Z),
$$

and, when declared, stable addition or removal of explicitly specified
counter-propagating null-response blocks. An arbitrary unimodular block is not
automatically trivial: an invertible Hall or neutral chiral phase can change
electrical or thermal response. Without the stated quotient, $(K,t)$ is a
coordinate presentation rather than a normal form. Its readouts are

$$
\nu=t^TK^{-1}t,
\qquad
\mathcal A_{\rm qp}=\Lambda^\ast/K\Lambda
\simeq\mathbb Z^N/K\mathbb Z^N,
$$

$$
Q_\ell/e=t^TK^{-1}\ell,
\qquad
\theta_\ell=\pi\ell^TK^{-1}\ell,
$$

$$
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell'.
$$

Charges are defined modulo local integer charge and the angles modulo $2\pi$
(with the transparent-local-fermion qualification on a spin theory). The
transport units are

$$
\sigma_{xy}=\nu\frac{e^2}{h},
\qquad
\frac{\kappa_{xy}}{T}
=
c_-\frac{\pi^2k_B^2}{3h}
$$

when the edge is thermally equilibrated.
For a minimal Abelian Chern--Simons presentation $c_-=\operatorname{sig}(K)$;
any separately stacked invertible sector must be recorded rather than silently
removed.

**Proof.** In coordinates, stable quasiparticle lines are vectors
$\ell\in\Lambda^\ast\simeq\mathbb Z^N$. Local excitations form the sublattice
$K\Lambda\subset\Lambda^\ast$. Hence quotient-visible quasiparticles are

$$
\mathcal A_{\rm qp}
=
\Lambda^\ast/K\Lambda
\simeq
\mathbb Z^N/K\mathbb Z^N.
$$

The electromagnetic charge readout pairs $t$ with the inverse holonomy
pairing:

$$
Q_\ell/e=t^TK^{-1}\ell.
$$

The Abelian linking form gives mutual braiding

$$
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell',
$$

and self-statistics

$$
\theta_\ell=\pi\ell^TK^{-1}\ell.
$$

The response of the edge/collar current to the background electromagnetic
connection is

$$
\nu=t^TK^{-1}t.
$$

Thus a declared Abelian $U(1)^N$ Chern-Simons branch has the standard
$K$-matrix readout. This is a recovery of conventional mathematics in the OPH
receipt ledger, not a new classification theorem for every Abelian spin
topological order. $\square$

## Corollary 2.1: Laughlin Branch

For one generator,

$$
K=(m),\qquad t=(1),
$$

so

$$
\nu=\frac1m,\qquad
\mathcal A_{\rm qp}=\mathbb Z_m,
\qquad
Q_\ell/e=\frac{\ell}{m},
\qquad
\theta_\ell=\frac{\pi\ell^2}{m}.
$$

Odd $m$ gives the fermionic Laughlin branch.

## Corollary 2.2: Jain Branch

For the Jain branch with $n$ filled composite-fermion levels and $2p$ flux
attachments, take

$$
K=I_n+2p\,\mathbf 1_n\mathbf 1_n^T,
\qquad
t=(1,\ldots,1)^T.
$$

Using Sherman-Morrison,

$$
K^{-1}
=
I_n-\frac{2p}{1+2pn}\mathbf 1_n\mathbf 1_n^T,
$$

which gives

$$
\nu=t^TK^{-1}t
=
n-\frac{2pn^2}{1+2pn}
=
\frac{n}{2pn+1}.
$$

For the reverse-flux convention use

$$
K=-I_n+2p\,\mathbf 1_n\mathbf 1_n^T,
$$

which gives $\nu=n/(2pn-1)$ when $2pn>1$.

## Proposition 3: Abelian Hierarchy Presentation

**Statement.** In the conventional Haldane-Halperin daughter-fluid
construction, adjoining one statistical gauge field gives the presentation

$$
K'=
\begin{pmatrix}
K & \ell\\
\ell^T & m
\end{pmatrix},
\qquad
t'=
\begin{pmatrix}
t\\
0
\end{pmatrix},
$$

provided $\ell$ and $m$ obey the required integrality, fermionic parity,
nondegeneracy, and charge-conservation conditions. Iterating an admissible step
gives an Abelian hierarchy presentation.

**Proof.** The added statistical gauge field enlarges the integral lattice,
pairs with the old fields through $\ell$, and has self-pairing $m$. It carries
no independent electromagnetic charge-vector entry in this convention. The
standard Abelian readout formulas then apply to $(K',t')$. $\square$

This block extension must not be conflated with arbitrary anyon condensation.
General anyon condensation requires a condensable bosonic algebra and can
identify or confine sectors rather than merely append a generator.

## Non-Abelian Hall Normal Forms

An Abelian normal form is pointed: every simple sector is invertible and fusion
is group-like. A non-Abelian normal form contains at least one simple object
$a$ with quantum dimension $d_a>1$, equivalently a fusion product

$$
a\times b=\sum_c N_{ab}^c\,c
$$

with more than one allowed output or a non-invertible simple object.

In OPH terms, Abelian holonomy is a central cycle obstruction. Non-Abelian Hall
order requires noncentral repair-sector data, because the residue is not
captured by a single commuting charge holonomy group.

## Theorem 4: Non-Abelian Repair-Sector Criterion

**Statement.** If a Hall normal form has a simple defect $a$ with
$d_a>1$, then the OPH normal form cannot be represented by Abelian holonomy
alone. It requires a genuinely noncentral repair sector represented here by
braided fusion-category data. Crossed-module data may model restricted
higher-gauge examples but are not asserted to classify arbitrary non-Abelian
anyon theories.

**Proof.** Abelian holonomy sectors form a pointed category. Each simple sector
is invertible, fusion has a unique output, and all quantum dimensions are one.
If $d_a>1$, the corresponding sector is not invertible. Its fusion rule cannot
be encoded by a group-valued Abelian holonomy class. The finite-cutoff ledger
must therefore retain its non-pointed fusion and braid data. $\square$

## Theorem 5: Conditional $5/2$ Noncentral Implication

**Statement.** Let $q\in Q_r^{5/2}$ be a Hall quotient state whose observed
readout includes

$$
\nu=\frac52,\qquad e^\ast=\frac e4,
$$

and an Ising/Majorana fusion package

$$
\sigma\times\sigma=1+\psi,\qquad
\psi\times\psi=1,\qquad
\sigma\times\psi=\sigma.
$$

Then $n_r^{5/2}(q)$ cannot be purely Abelian. It must include a noncentral
neutral repair sector.

The Ising/Majorana fusion package is a premise to be established by an
independent receipt. It does not follow from $\nu=5/2$ and $e^\ast=e/4$.

**Proof.** In a purely Abelian $K$ matrix normal form, quasiparticle sectors
form the finite Abelian group $\mathbb Z^N/K\mathbb Z^N$, and fusion is group
addition. Fusion of two simple sectors gives a single simple sector. The Ising
rule

$$
\sigma\times\sigma=1+\psi
$$

has two distinct outputs. Thus $\sigma$ is non-invertible and has quantum
dimension $\sqrt2$. By Theorem 4, this is impossible in a pointed Abelian
holonomy group. The $5/2$ normal form must include a noncentral neutral repair
sector. $\square$

## Theorem 6: Charge-Readout Non-Injectivity

**Statement.** The readout consisting only of filling fraction and minimal
charge,

$$
(\nu,e^\ast)=\left(\frac52,\frac e4\right),
$$

does not injectively determine the $5/2$ topological order.

**Proof.** Distinct candidate orders share the same charged $U(1)$ sector. The
Pfaffian, anti-Pfaffian, PH-Pfaffian, and related reconstructed sectors all
carry a charge $e/4$ quasiparticle and the same filling fraction. They differ
in neutral sector, chirality, edge content, particle-hole realization, and
thermal Hall response. Since the map

$$
\mathrm{HallRead}_{\mathrm{charge}}:
N_r^{5/2}\to
\mathbb R^2,
\qquad
n\longmapsto(\nu(n),e^\ast(n)/e)
$$

collapses those neutral and edge distinctions, it is not injective. Charge
holonomy alone cannot select the material $5/2$ topological order. $\square$

## Theorem 7: Normal-Form Non-Selection

**Statement.** Let

$$
N_r^{5/2}
=
N_{\rm Pf}\sqcup N_{\rm APf}\sqcup N_{\rm PH}\sqcup N_{\rm Ab}\sqcup\cdots
$$

be a finite or regulated set of candidate $5/2$ normal-form sectors. The
normal-form projector

$$
n_r:Q_r^{5/2}\to N_r^{5/2}
$$

does not select a unique material sector unless the branch supplies an
additional quotient-intrinsic ensemble or action

$$
\mu_r(q)=Z_r^{-1}m_r(q)e^{-S_r(q)}.
$$

**Proof.** Let $i_r:N_r^{5/2}\hookrightarrow Q_r^{5/2}$ be the declared
inclusion with $n_r\circ i_r=\mathrm{id}$. The induced self-map on probability
laws over $Q_r^{5/2}$ is

$$
\mathcal C(\mu)=(i_r)_\#(n_r)_\#\mu.
$$

It is idempotent:

$$
\mathcal C^2=\mathcal C.
$$

Every law supported on the embedded normal forms is fixed. In particular, a law supported
entirely on $N_{\rm Pf}$, a law supported entirely on $N_{\rm APf}$, and a
law supported entirely on $N_{\rm PH}$ are all fixed points of
$\mathcal C$. Hence $n_r$ is a canonicalizer, not a selector. A unique
material sector requires extra source-side weights or action data. $\square$

## Theorem 8: Strong Independence

**Statement.** For any probability vector

$$
(p_X)_{X\in\mathcal X},
\qquad
\sum_X p_X=1,
$$

over a finite set of candidate $5/2$ sectors $\mathcal X$, there exists a
quotient-intrinsic action $S_r$ whose induced sector probabilities are
$(p_X)$. Therefore normal-form data alone cannot imply a unique sector.

**Proof.** Let $N_{r,X}\subset N_r^{5/2}$ be the normal forms in sector $X$ and
$Q_{r,X}=n_r^{-1}(N_{r,X})$. At finite regulator choose a base measure $m_r$
with $0<M_X:=m_r(Q_{r,X})<\infty$. Define an extended action on $Q_r$ by

$$
e^{-S_r(q)}=\frac{p_X}{M_X}
\quad
\text{for } q\in Q_{r,X},
$$

where $p_X=0$ means $S_r=+\infty$ on that sector.

Then the partition function is

$$
Z_r=\sum_X\int_{Q_{r,X}}m_r(dq)\frac{p_X}{M_X}
=\sum_X p_X=1.
$$

The induced probability of sector $X$ is $p_X$. Since every probability
vector over candidate sectors can be realized by some quotient action, the
normal-form data alone cannot force a unique sector. This constructed action is
a countermodel, not an admissible prediction: choosing $p_X$ from the desired
material outcome would violate the no-target-leak rule. $\square$

## Missing Selector Object

A positive $5/2$ selector requires the object

$$
\mathfrak S_r^{5/2}
=
\left(
Q_r^{5/2},
m_r,
S_r,
c_{sr},
\{P_{r,X}\}_{X\in\mathcal X},
\mathrm{HallRead}_r
\right).
$$

$Q_r^{5/2}$ is the finite physical quotient. $m_r$ is the base measure.
$S_r$ is a dimensionless source-side action (or an explicitly declared
$\beta H$, Euclidean action, or ground-state limit). $c_{sr}:Q_s\to Q_r$ are
refinement maps. $P_{r,X}$ are sector projectors. The readout map records
charge, heat, edge, tunneling, and interferometry data.

The action $S_r$ must include the material terms that can split the $5/2$
candidates: Coulomb scale, finite width, Landau-level mixing, particle-hole
breaking, disorder, spin polarization, edge confinement, edge equilibration,
and any substrate-specific term used by the claimed selector.

## Theorem 9: Conditional Selector

**Statement.** Suppose $\mathfrak S_r^{5/2}$ is supplied and the sector
weights

$$
W_r(X)
=
\sum_{q\in P_{r,X}Q_r^{5/2}}
m_r(q)e^{-S_r(q)}
$$

are positive and have a unique maximizer $X_\star$. Define the log-weight
margin

$$
g_r^W
=
\min_{Y\ne X_\star}
\left[
-\log W_r(Y)+\log W_r(X_\star)
\right]
$$

If $g_r^W\to+\infty$ for the same eventual sector $X_\star$, then its
sector probability tends to one. A finite positive margin identifies only the
maximum-weight sector at that regulator; it does not imply concentration or
refinement stability by itself.

**Proof.** The sector probability is

$$
\Pr_r(X)=\frac{W_r(X)}{\sum_Y W_r(Y)}.
$$

For $Y\ne X_\star$,

$$
\frac{W_r(Y)}{W_r(X_\star)}
\le e^{-g_r^W}.
$$

Thus

$$
\Pr_r(X_\star)
=
\frac{1}{1+\sum_{Y\ne X_\star}W_r(Y)/W_r(X_\star)}
\ge
\frac{1}{1+(|\mathcal X|-1)e^{-g_r^W}}.
$$

If $g_r^W\to+\infty$, the bound tends to one. $\square$

## Theorem 10: Refinement-Compatible Selector

**Statement.** Let $c_{sr}:Q_s\to Q_r$ be refinement maps. If

$$
(c_{sr})_\#\mu_s=\mu_r
$$

and the sector labels intertwine with refinement up to null sets,

$$
P_{s,X}Q_s
=
c_{sr}^{-1}(P_{r,X}Q_r)
\quad(\mathrm{mod}\ \mu_s),
$$

for each sector $X$, then the sector probabilities are refinement-compatible:

$$
\Pr_s(X)=\Pr_r(X).
$$

More generally, let $\epsilon_{sr}$ bound the total-variation defect of the
pushforward law and define the sector-label mismatch

$$
\delta_{sr,X}
=
\mu_s\!\left(
P_{s,X}Q_s\,\triangle\,
c_{sr}^{-1}(P_{r,X}Q_r)
\right).
$$

Then

$$
|\Pr_s(X)-\Pr_r(X)|
\le
\epsilon_{sr}+\delta_{sr,X}.
$$

**Proof.** Exact compatibility gives

$$
\Pr_s(X)
=
\mu_s(P_{s,X}Q_s)
=
(c_{sr})_\#\mu_s(P_{r,X}Q_r)
=
\mu_r(P_{r,X}Q_r)
=
\Pr_r(X).
$$

The approximate statement follows by adding and subtracting
$\mu_s(c_{sr}^{-1}(P_{r,X}Q_r))$, then applying total variation and the
symmetric-difference bound. $\square$

## Theorem 11: No-Target-Leak Rule

**Statement.** A proposed $5/2$ selector is invalid as an OPH source selector
if $S_r^{5/2}$, $m_r$, or the sector projectors are chosen using the target
material order or the measurements that the selector is supposed to predict.

**Proof.** OPH finite quotient ensembles require the base weight, action,
observables, refinement maps, and source expectations to be declared from the
source side. If $S_r^{5/2}$ is chosen because the target is PH-Pfaffian, or
because a thermal Hall measurement prefers $c_-=5/2$, then the action is
posterior fitting. Such a selector violates the finite-ensemble firewall.
$\square$

## Fractional Exciton Extension

The fractional-exciton sandbox uses the same quotient logic with one extra
public readout layer. A fractional Chern or fractional Hall material is
presented as

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
\right),
\qquad
Q_{x,r}=\Sigma^{\mathrm{mat}}_{x,r}/\Gamma_{x,r}.
$$

The source law or Hamiltonian is mandatory. A normal form classifies the
surviving sector, but it does not select which material phase a sample realizes.
The required source-side tags are

$$
\mathrm{SOURCE\_LAW\_REQUIRED},
\qquad
\mathrm{NORMAL\_FORM\_IS\_NOT\_SELECTOR}.
$$

For a fractional Chern sandbox, a projected interacting Chern-band Hamiltonian
and its optical extension are only candidate source objects until they emit a
material promotion certificate:

$$
\mathrm{PhaseCert}_{x,r}
=
\left(
C,\nu,\Delta,G,\sigma_{xy},\mathrm{Pump},
\mathrm{Edge},\mathrm{Ent},\mathrm{Anyon},\mathrm{Modular},
\mathrm{Refine}
\right).
$$

The certificate must freeze the source Hamiltonian before comparison and carry
public receipts for the active band projector, Chern number, band geometry,
many-body gap, ground-sector degeneracy, flux pump, Hall conductance, edge
spectrum, topological-sector ledger, refinement stability, and no-target-leak
audit. If the certificate map is not injective on candidate sectors, the
correct result is `SECTOR_AMBIGUOUS`.

Optical fractional excitons add a module category

$$
\mathcal M_{\mathrm{opt}}
$$

over the topological sector category. Each optically active object has a
predicted topological shadow $\tau_{\mathrm{pred}}$, modeled total
electromagnetic charge $Q_{\mathrm{tot}}$, polarization, oscillator data, and
binding-energy data. Those topological labels are candidate annotations. A
minimal measured record instead lies in an observable space
$\mathcal Y_{\mathrm{opt}}$:

$$
y=(\Delta E,I,\partial_g\Delta E,\mathrm{pol},\eta,\ldots)
\in\mathcal Y_{\mathrm{opt}}.
$$

Each frozen candidate $m$ supplies an acceptance region
$R_m\subseteq\mathcal Y_{\mathrm{opt}}$. The annotated prediction ledger is

$$
\widehat{\mathcal L}_{\mathrm{opt}}
=
\{(R_m,\tau_{\mathrm{pred}}(m),Q_{\mathrm{tot}}(m)):m\in\mathcal X\}.
$$

A gate or field slope identifies charge only when the lever arm and all
binding or nuisance derivatives are independently bounded. A neutral
composite may have

$$
Q_{\mathrm{tot}}=0,
\qquad
\tau_{\mathrm{pred}}^{\mathrm{rel}}\ne 1,
$$

so its leading charge slope may vanish while an internal relative fusion
channel remains nontrivial. A strictly local photon operator preserves the
global superselection sector. The complete state therefore needs a
compensating anyon, a pre-existing sector, or an explicitly modeled boundary
or nonlocal process. An optical identification requires frozen candidate
regions that separate at the declared resolution and a record compatible with
exactly one region; inserting $\tau_{\mathrm{pred}}$ into the measured tuple
would be circular. Otherwise the correct result is
`OPTICAL_SECTOR_AMBIGUOUS`.

The companion problem note
[Fractional Excitons as OPH Quotient-Sector Readouts](fractional_excitons_as_oph_quotient_sector_readouts.md)
records the material presentation, source-law boundary, Hamiltonian promotion
certificate, optical module, line-fan theorem, identifiability theorem, and
simulator quotient-correctness certificate. The matching implementation
surfaces are `oph-physics-sim/oph_fractional/` and
`reverse-engineering-reality/code/particles/fractional/`; their generated
bundles intentionally remain diagnostic until a material-specific Hamiltonian
proof is supplied.

## Summary

Let $Q_r^{\mathrm{Hall}}$ denote the declared continuum-FQH or FCI/FQAH
quotient, with conserved electromagnetic $U(1)$ and an accepted typed
canonicalizer $n_r$. Given a conventional gapped topological-order ledger, OPH
packages its bulk invariants, compatible edge realization, readback, and
evidence receipts as a self-reading patch. On the Abelian $U(1)^N$
Chern-Simons branch the normal form is an equivalence class $[K,t]$, with

$$
\nu=t^TK^{-1}t,\qquad
\mathcal A_{\rm qp}=\Lambda^\ast/K\Lambda,
$$

$$
Q_\ell/e=t^TK^{-1}\ell,\qquad
\theta_\ell=\pi\ell^TK^{-1}\ell.
$$

The conventional hierarchy admits block-matrix presentations subject to its
integrality, parity, and charge constraints. At $\nu=5/2$, an independently
established Ising/Majorana package with
$\sigma\times\sigma=1+\psi$ is impossible as Abelian holonomy and requires a
noncentral repair sector; $\nu$ and $e/4$ alone do not. The normal-form
canonicalizer does not select a material phase. Selection requires a frozen
source Hamiltonian or ensemble and refinement-stable evidence.

For fractional-exciton and fractional-Chern optical branches, the same selector
rule applies. The optical line fan is evidence only when it descends from a
frozen source, a certified topological ledger, a quotient-descended optical
module, and a simulator certificate with representative invariance, quotient
lumpability, no orbit-size bias, refinement compatibility, and no target leak.

## Promotion Tests and Falsifiers

A material claim must freeze the continuum-FQH or FCI/FQAH source before
comparison and report, as applicable, the many-body gap, spectral flow and
charge pump, Hall conductance, torus ground-sector structure, entanglement or
modular data, edge anomaly, thermal Hall response, and refinement uncertainty.
The OPH contribution fails as an identifying ledger if representative changes
alter those readouts, if the same frozen receipts remain non-injective on the
candidate sectors, or if the claimed selector uses the target measurement it is
supposed to predict.

## References

- R. B. Laughlin, "Anomalous Quantum Hall Effect: An Incompressible Quantum
  Fluid with Fractionally Charged Excitations", Physical Review Letters 50,
  1395, 1983. https://doi.org/10.1103/PhysRevLett.50.1395
- F. D. M. Haldane, "Fractional Quantization of the Hall Effect: A Hierarchy of
  Incompressible Quantum Fluid States", Physical Review Letters 51, 605, 1983.
  https://doi.org/10.1103/PhysRevLett.51.605
- J. K. Jain, "Composite-Fermion Approach for the Fractional Quantum Hall
  Effect", Physical Review Letters 63, 199, 1989.
  https://doi.org/10.1103/PhysRevLett.63.199
- X.-G. Wen and A. Zee, "Classification of Abelian Quantum Hall States and
  Matrix Formulation of Topological Fluids", Physical Review B 46, 2290, 1992.
  https://doi.org/10.1103/PhysRevB.46.2290
- G. Moore and N. Read, "Nonabelions in the Fractional Quantum Hall Effect",
  Nuclear Physics B 360, 362, 1991.
  https://doi.org/10.1016/0550-3213(91)90407-O
- J. Cai et al., "Signatures of Fractional Quantum Anomalous Hall States in
  Twisted MoTe2", Nature 622, 63, 2023.
  https://doi.org/10.1038/s41586-023-06289-w
