# Fractional Quantum Hall States

**Status:** solved as a conditional OPH theorem package. Standalone markdown
supplemental writeup for public reading and OPH Sage ingestion.

## Introduction

Fractional Hall physics is hard in legacy language because the measured plateau
is a macroscopic transport record, while the microscopic many-electron state,
edge channels, quasiparticle statistics, disorder, and sample geometry all
enter the answer. The problem statement is to classify which gapped charged
phases survive as public transport plateaus and to identify what extra data
selects the $5/2$ sector. In OPH the boundary is the object: a Hall droplet is
a charged collar with edge records and holonomy defects, so the surviving phases
are its repair normal forms. Abelian plateaus become the ordinary integral
edge/holonomy ledger; $5/2$ becomes the branch that requires a noncentral
repair sector.

## Why Legacy Physics Gets Stuck

Legacy fractional Hall physics has powerful local ingredients: Landau-level
projection, trial wavefunctions, composite-fermion constructions, Chern-Simons
effective theories, edge conformal field theories, and numerical spectra. The
hard part is that the laboratory plateau is a boundary-visible transport record,
while the proposed explanation usually lives in a microscopic state, an
effective topological field theory, or an edge model. Several gauge-distinct
presentations can produce the same coarse transport plateau, and several
physically different edge or defect sectors can agree on part of the readout.

The $\nu=5/2$ sector is the sharp failure mode. Pfaffian, anti-Pfaffian,
PH-Pfaffian, Abelian, and reconstructed candidates can share filling fraction
and partial charge data while disagreeing about fusion, chirality, thermal
transport, interferometry, and edge reconstruction. A plateau alone therefore
does not select the observer-facing sector. Legacy physics can name candidate
phases and compute many invariants, but it lacks a single substrate-neutral
normal-form object that says which public boundary records identify the phase
and which hidden details are gauge.

In that framing, the classification target is underdetermined: the same public
plateau can sit over more than one observer-distinct sector.

## Why OPH Makes It Solvable

OPH makes the finite Hall collar, not the hidden bulk presentation, the object
being classified. The physical state is the quotient of droplet, flux, active
Landau level, source Hamiltonian, edge algebra, records, repair instruments, and
checkpoint data by gauge representatives, basis choices, mesh labels, electron
labels, edge charts, port labels, repair schedules, hidden carrier coordinates,
and inert ancillas.

That quotient turns the problem into a normal-form classification. Abelian
plateaus are identified by an integral edge/holonomy ledger. Non-Abelian
branches require defect and higher-defect data. The $5/2$ branch becomes
solvable only when the boundary record includes enough data to distinguish the
noncentral repair sector: charge, heat, tunneling, interferometry, fusion,
braiding, chirality, and refinement-stable gap information. OPH is unique here
because it makes the observer-facing equivalence relation load-bearing. Two
presentations that differ only by hidden labels collapse to one state, while two
sectors that differ in boundary-visible fusion or transport receipts remain
distinct.

## Abstract

Fractional Hall systems fit the OPH pattern as regulated Hall-collar quotients.
On the gapped branch, local bulk mismatch repairs inside a fixed topological
sector, while the quotient-visible residue lives in the edge algebra, charge
holonomy, and defect or higher-defect data. The Abelian branch reduces to the
usual integral $K$ matrix normal form. Laughlin states, Jain sequences, and
the Abelian hierarchy appear as one-generator or iterated holonomy-refinement
normal forms. At $\nu=5/2$, Ising/Majorana fusion forces a noncentral repair
sector. Normal-form projection classifies Pfaffian, anti-Pfaffian,
PH-Pfaffian, Abelian, and reconstructed sectors, while material selection
requires a quotient-intrinsic action and refinement-stable sector gap.

## Source Branch

The source branch is

```math
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
```

The fixed-cutoff observer patch is

```math
\mathsf O_i=
\left(
\mathcal A_i,\rho_i,\mathcal R_i,
(\mathcal I_e,\pi_{i,e})_{e\ni i},
\mathcal U_i,\mathrm{Chk}_i
\right),
```

where $\mathcal A_i$ is the accessible algebra, $\rho_i$ the state,
$\mathcal R_i$ the record algebra, $\mathcal I_e$ the overlap-visible
interface algebra, $\pi_{i,e}$ the visible restriction, $\mathcal U_i$ the
allowed repair instruments, and $\mathrm{Chk}_i$ checkpoint data.

The Hall droplet is a support/collar chart. Hidden carrier coordinates, edge
charts, port labels, and mesh labels are silent when they preserve visible
interface statistics, record readout, repair instruments, and checkpoint
continuation.

## Branch Assumptions

The theorem assumes a finite regulated two-dimensional charged electron system
with conserved electromagnetic $U(1)$, nonzero magnetic flux, a bulk mobility
gap on the plateau, observer-visible boundary/collar algebra, and
quotient-visible charge, heat, tunneling, and interferometry records:

```math
\mathrm{Hall}_r =
\left(
2+1\mathrm{D},
U(1)_{\rm em},
B\neq0,
\Delta_{\rm bulk}>0,
\mathcal A_{\partial,r},
\mathcal R_{\partial,r}
\right).
```

The result fixes the normal-form shape of fractional Hall response. It carries
no claim that an arbitrary Hamiltonian enters a fractional Hall phase.

## Hall-Collar Quotient

Let the finite Hall presentation space $\Sigma_r^{\mathrm{Hall}}$ consist of
tuples

```math
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
```

Here $\mathcal D_r$ is the finite droplet, annulus, cylinder, sphere, or torus
geometry; $\partial\mathcal D_r$ is the Hall collar; $B_r$ is flux data;
$N_{e,r}$ and $N_{\phi,r}$ are electron and flux counts; $\mathcal H_r$
is the finite many-body Hilbert space; $P_{\mathrm{LL},r}$ is the active
Landau-level projection; $H_r^{\mathrm{eff}}$ is the effective source
Hamiltonian or action; $\mathcal A_{\partial,r}$ is the collar algebra;
$\mathcal R_{\partial,r}$ is the record algebra; $\mathcal U_{\partial,r}$
is the allowed repair/update family; and $\mathrm{Chk}_r$ is
refinement/checkpoint data.

Let $\Gamma_r^{\mathrm{Hall}}$ quotient away gauge representatives, basis
choices, mesh labels, electron labels, edge charts, port labels, repair
schedules, hidden carrier coordinates, and inert ancillas. The physical
quotient is

```math
Q_r^{\mathrm{Hall}}
=
\Sigma_r^{\mathrm{Hall}}/\Gamma_r^{\mathrm{Hall}}.
```

For the $\nu=5/2$ branch,

```math
\Sigma_r^{5/2}\subset \Sigma_r^{\mathrm{Hall}},
\qquad
\nu_r=\frac{N_{e,r}}{N_{\phi,r}}\to \frac52.
```

After subtracting two filled lower Landau levels,

```math
\nu_{\rm active,r}\to \frac12,
\qquad
Q_r^{5/2}=\Sigma_r^{5/2}/\Gamma_r^{5/2}.
```

## Hall Normal Form

The accepted repair map is

```math
n_r^{\mathrm{Hall}}:
Q_r^{\mathrm{Hall}}\to N_r^{\mathrm{Hall}}.
```

A Hall normal form is

```math
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
```

$\mathcal E_{\partial,r}$ is the repaired edge/collar algebra.
$\mathcal H_{\mathrm{hol},r}$ is the Abelian holonomy lattice or
non-Abelian higher-holonomy datum. $\mathcal C_r$ is the defect/fusion
category. $\mathfrak e_r$ is the local electron object. The charge grading is

```math
\chi_r:\operatorname{Obj}(\mathcal C_r)\to \mathbb R/\mathbb Z.
```

$\Theta_r$ is the spin, braid, and modular-data package.

The readout map is

```math
\operatorname{HallRead}_r:
N_r^{\mathrm{Hall}}\to\mathcal Y_r^{\mathrm{Hall}},
```

```math
\operatorname{HallRead}_r(n)
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
```

The components are filling fraction, electrical Hall conductance, minimal
quasiparticle charge, chiral central charge, anyon set, fusion rules, exchange
data, braiding data, modular/geometric response, and interferometry response.

## Theorem 1: Hall Normal Form

**Statement.** Let $Q_r^{\mathrm{Hall}}$ be the finite physical quotient of a
gapped $2+1$ dimensional charged Hall system with conserved
$U(1)_{\rm em}$, nonzero magnetic flux, a bulk mobility gap, and an
observer-visible boundary/collar algebra. Assume the accepted OPH repair law on
this branch is terminating, quotient-descended, locally confluent, and
repair-complete. Then every stable Hall phase has a schedule-independent normal
form

```math
n_r^{\mathrm{Hall}}(q)
=
\left(
\text{edge algebra},
\text{charge holonomy},
\text{defect/higher-defect data},
\text{record readout}
\right).
```

If the residual defect datum is Abelian, it is classified by Abelian holonomy.
If it is genuinely noncentral, it is classified by higher-gauge/fusion data.

**Proof.** OPH physical data factor through accessible algebras, records,
visible interface restrictions, repair instruments, and checkpoint
continuation. Hidden presentation data are quotient artifacts.

The consensus theorem gives schedule-independent quotient normal forms when the
repair relation is terminating, quotient-descended, locally confluent, and
repair-complete. The same theorem surface states that pairwise overlap agreement
does not guarantee global agreement; the residual obstruction is holonomic,
Abelian in the central branch and crossed-module Cech in the genuinely
noncentral branch.

In a gapped Hall phase, all local bulk mismatches that preserve the mobility gap
repair inside the same topological sector. The quotient-visible residue after
local repair is boundary/collar data or global obstruction data. A stable Hall
phase is represented by edge algebra, charge holonomy, defect/higher-defect
data, and observer-facing readout. $\square$

## Theorem 2: Abelian Normal Forms Are $K$ Matrix Normal Forms

On the Abelian branch, let $\Lambda\simeq\mathbb Z^N$ be the stable defect
lattice, let

```math
K:\Lambda\times\Lambda\to\mathbb Z
```

be a nondegenerate symmetric integral pairing, and let $t\in\Lambda^\ast$ be
the electromagnetic charge vector.

**Statement.** On the Abelian branch,

```math
n_r^{\mathrm{Hall}}(q)=(K,t)
```

has readouts

```math
\nu=t^TK^{-1}t,
\qquad
\mathcal A_{\rm qp}=\mathbb Z^N/K\mathbb Z^N,
```

```math
Q_\ell/e=t^TK^{-1}\ell,
\qquad
\theta_\ell=\pi\ell^TK^{-1}\ell,
```

```math
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell'.
```

**Proof.** In the Abelian OPH holonomy branch, stable defects are vectors
$\ell\in\Lambda\simeq\mathbb Z^N$. Local electron excitations form the
sublattice $K\Lambda$. Hence quotient-visible quasiparticles are

```math
\mathcal A_{\rm qp}
=
\Lambda/K\Lambda
\simeq
\mathbb Z^N/K\mathbb Z^N.
```

The electromagnetic charge readout pairs $t$ with the inverse holonomy
pairing:

```math
Q_\ell/e=t^TK^{-1}\ell.
```

The Abelian linking form gives mutual braiding

```math
\theta_{\ell,\ell'}=2\pi\ell^TK^{-1}\ell',
```

and self-statistics

```math
\theta_\ell=\pi\ell^TK^{-1}\ell.
```

The response of the edge/collar current to the background electromagnetic
connection is

```math
\nu=t^TK^{-1}t.
```

Thus the Abelian OPH holonomy normal form is exactly the $K$ matrix Hall
normal form. $\square$

## Corollary 2.1: Laughlin Branch

For one generator,

```math
K=(m),\qquad t=(1),
```

so

```math
\nu=\frac1m,\qquad
\mathcal A_{\rm qp}=\mathbb Z_m,
\qquad
Q_\ell/e=\frac{\ell}{m},
\qquad
\theta_\ell=\frac{\pi\ell^2}{m}.
```

Odd $m$ gives the fermionic Laughlin branch.

## Corollary 2.2: Jain Branch

For the Jain branch with $n$ filled composite-fermion levels and $2p$ flux
attachments, take

```math
K=I_n+2p\,\mathbf 1_n\mathbf 1_n^T,
\qquad
t=(1,\ldots,1)^T.
```

Using Sherman-Morrison,

```math
K^{-1}
=
I_n-\frac{2p}{1+2pn}\mathbf 1_n\mathbf 1_n^T,
```

which gives

```math
\nu=t^TK^{-1}t
=
n-\frac{2pn^2}{1+2pn}
=
\frac{n}{2pn+1}.
```

The opposite flux orientation gives $\nu=n/(2pn-1)$.

## Theorem 3: Hierarchy as Normal-Form Refinement

**Statement.** Suppose an Abelian Hall normal form $(K,t)$ has a stable
quasiparticle class $\ell$ that condenses as an additional repair-stable
defect generator. Then the refined OPH normal form is an enlarged $K$ matrix

```math
K'=
\begin{pmatrix}
K & \ell\\
\ell^T & m
\end{pmatrix},
\qquad
t'=
\begin{pmatrix}
t\\
q
\end{pmatrix},
```

where $m$ and $q$ encode the self-statistics and charge of the condensed
generator. Iterating this step gives the Abelian hierarchy.

**Proof.** A hierarchy step adds one stable quotient-visible defect generator to
the Abelian holonomy lattice. Adding a generator enlarges the lattice from
$\Lambda$ to $\Lambda'=\Lambda\oplus\mathbb Z e_{\rm add}$. The integral
pairing on the enlarged lattice restricts to $K$ on $\Lambda$, pairs the
additional generator with old defects by $\ell$, and assigns it self-pairing
$m$. The charge vector extends by $q$. The Abelian readout theorem gives the
usual hierarchy formulas from $(K',t')$. $\square$

## Non-Abelian Hall Normal Forms

An Abelian normal form is pointed: every simple sector is invertible and fusion
is group-like. A non-Abelian normal form contains at least one simple object
$a$ with quantum dimension $d_a>1$, equivalently a fusion product

```math
a\times b=\sum_c N_{ab}^c\,c
```

with more than one allowed output or a non-invertible simple object.

In OPH terms, Abelian holonomy is a central cycle obstruction. Non-Abelian Hall
order requires noncentral repair-sector data, because the residue is not
captured by a single commuting charge holonomy group.

## Theorem 4: Non-Abelian Repair-Sector Criterion

**Statement.** If a Hall normal form has a simple defect $a$ with
$d_a>1$, then the OPH normal form cannot be represented by Abelian holonomy
alone. It requires a genuinely noncentral repair sector, modeled at finite
cutoff by crossed-module or fusion-category data.

**Proof.** Abelian holonomy sectors form a pointed category. Each simple sector
is invertible, fusion has a unique output, and all quantum dimensions are one.
If $d_a>1$, the corresponding sector is not invertible. Its fusion rule cannot
be encoded by a group-valued Abelian holonomy class. OPH's noncentral
obstruction branch is the higher-gauge/crossed-module branch, and the
finite-cutoff edge-sector language records it as fusion-category data.
$\square$

## Theorem 5: $5/2$ Noncentral Admissibility

**Statement.** Let $q\in Q_r^{5/2}$ be a Hall quotient state whose observed
readout includes

```math
\nu=\frac52,\qquad e^\ast=\frac e4,
```

and an Ising/Majorana fusion package

```math
\sigma\times\sigma=1+\psi,\qquad
\psi\times\psi=1,\qquad
\sigma\times\psi=\sigma.
```

Then $n_r^{5/2}(q)$ cannot be purely Abelian. It must include a noncentral
neutral repair sector.

**Proof.** In a purely Abelian $K$ matrix normal form, quasiparticle sectors
form the finite Abelian group $\mathbb Z^N/K\mathbb Z^N$, and fusion is group
addition. Fusion of two simple sectors gives a single simple sector. The Ising
rule

```math
\sigma\times\sigma=1+\psi
```

has two distinct outputs. Thus $\sigma$ is non-invertible and has quantum
dimension $\sqrt2$. By Theorem 4, this is impossible in a pointed Abelian
holonomy group. The $5/2$ normal form must include a noncentral neutral repair
sector. $\square$

## Theorem 6: Charge-Readout Non-Injectivity

**Statement.** The readout consisting only of filling fraction and minimal
charge,

```math
(\nu,e^\ast)=\left(\frac52,\frac e4\right),
```

does not injectively determine the $5/2$ topological order.

**Proof.** Distinct candidate orders share the same charged $U(1)$ sector. The
Pfaffian, anti-Pfaffian, PH-Pfaffian, and related reconstructed sectors all
carry a charge $e/4$ quasiparticle and the same filling fraction. They differ
in neutral sector, chirality, edge content, particle-hole realization, and
thermal Hall response. Since the map

```math
\operatorname{HallRead}_{\rm charge}:
N_r^{5/2}\to
\{(\nu,e^\ast)\}
```

collapses those neutral and edge distinctions, it is not injective. Charge
holonomy alone cannot select the material $5/2$ topological order. $\square$

## Theorem 7: Normal-Form Non-Selection

**Statement.** Let

```math
N_r^{5/2}
=
N_{\rm Pf}\sqcup N_{\rm APf}\sqcup N_{\rm PH}\sqcup N_{\rm Ab}\sqcup\cdots
```

be a finite or regulated set of candidate $5/2$ normal-form sectors. The
normal-form projector

```math
n_r:Q_r^{5/2}\to N_r^{5/2}
```

does not select a unique material sector unless the branch supplies an
additional quotient-intrinsic ensemble or action

```math
\mu_r(q)=Z_r^{-1}m_r(q)e^{-S_r(q)}.
```

**Proof.** The induced map on probability laws is

```math
\mathcal C(\mu)=(n_r)_\#\mu.
```

It is idempotent:

```math
\mathcal C^2=\mathcal C.
```

Every law supported on $N_r^{5/2}$ is fixed. In particular, a law supported
entirely on $N_{\rm Pf}$, a law supported entirely on $N_{\rm APf}$, and a
law supported entirely on $N_{\rm PH}$ are all fixed points of
$\mathcal C$. Hence $n_r$ is a canonicalizer, not a selector. A unique
material sector requires extra source-side weights or action data. $\square$

## Theorem 8: Strong Independence

**Statement.** For any probability vector

```math
(p_X)_{X\in\mathcal X},
\qquad
\sum_X p_X=1,
```

over a finite set of candidate $5/2$ sectors $\mathcal X$, there exists a
quotient-intrinsic action $S_r$ whose induced sector probabilities are
$(p_X)$. Therefore normal-form data alone cannot imply a unique sector.

**Proof.** Let $N_{r,X}\subset N_r^{5/2}$ be the normal forms in sector $X$,
and let $M_X=|N_{r,X}|$, assuming finite regulator $r$. Choose uniform base
weight $m_r=1$. Define

```math
e^{-S_r(n)}=\frac{p_X}{M_X}
\quad
\text{for } n\in N_{r,X}.
```

Then the partition function is

```math
Z_r=\sum_X\sum_{n\in N_{r,X}}\frac{p_X}{M_X}
=\sum_X p_X=1.
```

The induced probability of sector $X$ is $p_X$. Since every probability
vector over candidate sectors can be realized by some quotient action, the
normal-form data alone cannot force a unique sector. $\square$

## Missing Selector Object

A positive $5/2$ selector requires the object

```math
\mathfrak S_r^{5/2}
=
\left(
Q_r^{5/2},
m_r,
S_r,
c_{sr},
\{P_{r,X}\}_{X\in\mathcal X},
\operatorname{HallRead}_r
\right).
```

$Q_r^{5/2}$ is the finite physical quotient. $m_r$ is the base measure.
$S_r$ is the source-side material action. $c_{sr}:Q_s\to Q_r$ are
refinement maps. $P_{r,X}$ are sector projectors. The readout map records
charge, heat, edge, tunneling, and interferometry data.

The action $S_r$ must include the material terms that can split the $5/2$
candidates: Coulomb scale, finite width, Landau-level mixing, particle-hole
breaking, disorder, spin polarization, edge confinement, edge equilibration,
and any substrate-specific term used by the claimed selector.

## Theorem 9: Conditional Selector

**Statement.** Suppose $\mathfrak S_r^{5/2}$ is supplied and the sector
weights

```math
W_r(X)
=
\sum_{q\in P_{r,X}Q_r^{5/2}}
m_r(q)e^{-S_r(q)}
```

have a unique maximizer $X_\star$. If the gap

```math
\Delta_r
=
\min_{Y\ne X_\star}
\left[
-\log W_r(Y)+\log W_r(X_\star)
\right]
```

satisfies $\Delta_r\to+\infty$ or stays above the declared refinement margin,
then the material selector chooses $X_\star$ in the refinement limit.

**Proof.** The sector probability is

```math
\Pr_r(X)=\frac{W_r(X)}{\sum_Y W_r(Y)}.
```

For $Y\ne X_\star$,

```math
\frac{W_r(Y)}{W_r(X_\star)}
\le e^{-\Delta_r}.
```

Thus

```math
\Pr_r(X_\star)
=
\frac{1}{1+\sum_{Y\ne X_\star}W_r(Y)/W_r(X_\star)}
\ge
\frac{1}{1+(|\mathcal X|-1)e^{-\Delta_r}}.
```

If $\Delta_r\to+\infty$, the bound tends to one. If $\Delta_r$ exceeds the
declared finite refinement margin, the selector is stable at that regulator.
$\square$

## Theorem 10: Refinement-Compatible Selector

**Statement.** Let $c_{sr}:Q_s\to Q_r$ be refinement maps. If

```math
(c_{sr})_\#\mu_s=\mu_r
```

and

```math
c_{sr}(P_{s,X}Q_s)\subseteq P_{r,X}Q_r
```

for each sector $X$, then the sector probabilities are refinement-compatible:

```math
\Pr_s(X)=\Pr_r(X).
```

If exact compatibility has total-variation defect $\epsilon_{sr}$, then

```math
|\Pr_s(X)-\Pr_r(X)|\le \epsilon_{sr}.
```

**Proof.** Exact compatibility gives

```math
\Pr_s(X)
=
\mu_s(P_{s,X}Q_s)
=
(c_{sr})_\#\mu_s(P_{r,X}Q_r)
=
\mu_r(P_{r,X}Q_r)
=
\Pr_r(X).
```

The approximate statement is the definition of total variation applied to the
indicator of $P_{r,X}Q_r$. $\square$

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

## Summary

Let $Q_r^{\mathrm{Hall}}$ be the finite OPH quotient of a gapped charged Hall
collar with conserved electromagnetic $U(1)$ and accepted repair normal form
$n_r$. Every stable fractional Hall phase has a quotient-visible
edge/holonomy normal form. On the Abelian branch this normal form is an
integral $K$ matrix pair $(K,t)$, with

```math
\nu=t^TK^{-1}t,\qquad
\mathcal A_{\rm qp}=\mathbb Z^N/K\mathbb Z^N,
```

```math
Q_\ell/e=t^TK^{-1}\ell,\qquad
\theta_\ell=\pi\ell^TK^{-1}\ell.
```

Hierarchy formation is normal-form refinement by adjoining surviving defect
generators. At $\nu=5/2$, any Ising/Majorana package with
$\sigma\times\sigma=1+\psi$ is impossible as Abelian holonomy and requires a
noncentral repair sector. The normal-form projector remains an idempotent
canonicalizer. Unique material selection requires a quotient-intrinsic ensemble
and refinement-stable sector gap.
