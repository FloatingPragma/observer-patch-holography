# High-Temperature Superconductivity

## Motivating Result

This note entered the queue after the nickelate surprise: high-pressure
measurements of La3Ni2O7 reported superconductivity signatures near 80 K
([Nature 621, 493-498, 2023](https://pubmed.ncbi.nlm.nih.gov/37437603/)).
Nickelates had long been plausible cousins of the cuprates, but this result
made the material-search problem concrete again. The OPH question is why a
strong local pairing clue promotes only when phase stiffness, defect clearance,
retention, and bulk evidence close on the same material record.

**Status:** conditional OPH materials-audit and inverse-design specification.
The note does not derive a microscopic pairing interaction, predict a new
material, or report a material-specific $T_c$. Its established ingredients are
standard superconductivity theory and measurement practice; its OPH-specific
contribution is to place those ingredients in an observer-patch, source-law,
and public-receipt contract. Standalone markdown supplemental writeup for
public reading and OPH Sage ingestion.

Date: 2026-07-11

The current package is a conditional material branch: candidate-local pairing,
phase ordering, defect-limited transport, material retention, and bulk magnetic
evidence are tracked separately. These are not asserted to be independent
necessary-and-sufficient laws of every superconductor. A simulator may estimate
declared observables and emit receipts, but neither the simulator nor the OPH
vocabulary defines the material phase.

## Introduction

The physical problem has two levels. First, one must explain and quantitatively
predict the pairing interaction, ordered state, and transition temperature of a
declared material. Second, one must establish that the resulting state is bulk,
stable under the stated pressure and processing history, and not a filamentary
or contact artifact. Strong correlations, lattice structure, disorder,
competing orders, phase fluctuations, and processing history can all matter.

Standard theory already supplies BCS/Eliashberg treatments for conventional
pairing, material-specific Hubbard and $`t\text{-}J`$ models, spin-fluctuation and
phonon models, density-functional and many-body pipelines, Ginzburg--Landau and
XY descriptions of phase stiffness and vortices, and established transport,
magnetic, and thermodynamic diagnostics. The unresolved part is a reliably
predictive microscopic account across unconventional material families and a
prospective search procedure that outperforms those declared baselines.

OPH does not solve the microscopic interaction in this note. It contributes an
auditable decomposition: a local-pair observable, a phase-ordering observable,
a defect/transport observable, a process-retention record, and a public evidence
bundle must be declared on the same material branch. When operational gate
temperatures are monotone, their minimum is a useful ranking heuristic, not a
universal theorem for physical $T_c$.

## Standard Physics and the Remaining Gap

Superconductivity has excellent material-class explanations. BCS and
Eliashberg theory work for conventional phonon-mediated materials, including
high-pressure hydrides in regimes where their approximations are controlled.
Hubbard
and $`t\text{-}J`$ models, density-functional pipelines, spin-fluctuation models,
phonon models, and strong-correlation numerics each explain part of the high
$T_c$ record. The obstruction is simple: high-temperature superconductivity has
more than one failure mode. Strong local pairing can coexist with weak phase
stiffness. A promising gap symmetry can be destroyed by disorder or competing
orders. A high-pressure hydride can have a large pairing scale while failing
ambient retention. A cuprate can have a strong pair branch while grain
boundaries, vortices, oxygen history, or processing scars block a reproducible
macroscopic state.

That makes the room-temperature design problem underdetermined by any single
diagnostic. A
candidate can look good in one diagnostic and fail the observer-facing claim.
Transport alone can be fooled by geometry, contacts, filamentary paths, or
inhomogeneous phases. A computed pairing scale alone does not prove a material
that carries bulk current without loss. An OPH audit can make explicit which
pairing, phase, defect, process-retention, and public evidence records are being
claimed together.

The scientific problem remains unsolved when a model cannot prospectively
predict the material, its transition, and the complete bulk evidence record.

## What OPH Adds

OPH represents the search as a finite material-quotient and repair-gate
contract.
The source branch keeps the physical facts that matter: composition, structure,
doping, strain, pressure, disorder, competing channels, and processing history.
It quotients away hidden coordinates, basis labels, gauge choices, mesh labels,
and presentation data that do not change the observer-facing material record.

The OPH-specific object is the observer-like patch and receipt structure: local
state, material boundaries and ports, gauge-invariant readback, persistent
records, physical or externally controlled update moves, checkpoints, and a
public evidence bundle. The proposed ranking score combines a declared local
pairing scale, phase-ordering scale, and defect-limited transport scale with
penalties for instability, disorder, toxicity, and synthesis failure. This
organizes an inverse-design protocol; it does not replace a microscopic
mechanism or establish uniqueness relative to conventional multi-objective
materials design.

## Abstract

High $T_c$ superconductivity is treated as a material quotient with three
conditional operational diagnostics: local pairing, phase ordering, and
defect-limited transport.
Cuprates, pnictides, nickelates, and hydrides are different branches of that
same audit pattern. OPH can rank material and processing candidates by a
declared bottleneck score, with penalties for instability, disorder, toxicity,
and synthesis failure. Such a ranking becomes a physical prediction only after
its source law, units, thresholds, uncertainty, and held-out baseline comparison
are frozen.

## Source Branch

The source branch is the finite observer-patch quotient branch. A candidate
material is represented by

```math
\mathcal B_{\rm mat}
=
\left(
\mathbf x,
G_{\mathbf x},
\{\mathcal A_i^{\mathbf x}\},
H_{\mathbf x},
\Gamma_{\mathbf x},
\mathcal C_{2e},
\mathcal C_{\rm comp},
\mathcal P_{\rm proc}
\right).
```

$\mathbf x$ encodes composition, structure, doping, strain, pressure,
disorder, and processing history. $G_{\mathbf x}$ is the finite patch graph.
$\mathcal A_i^{\mathbf x}$ are local accessible algebras.
$H_{\mathbf x}$ is the source Hamiltonian or action.
$\Gamma_{\mathbf x}$ removes gauge, coordinate, mesh, and presentation labels
only after every retained observable and transition law is shown to factor
through the equivalence. An unobserved variable is not discarded merely
because it is hidden. $\mathcal C_{2e}$ is the charge $2e$
pair channel set. $\mathcal C_{\rm comp}$ is the set of competing spin,
charge, orbital, nematic, lattice, CDW, SDW, PDW, and disorder channels.
$\mathcal P_{\rm proc}$ is the processing map.

The physical quotient is

```math
Q_{\mathbf x,r}=\Sigma_{\mathbf x,r}/\Gamma_{\mathbf x,r}.
```

A material patch is an observer-like OPH object:

```math
\mathsf O_i^{\rm mat}
=
\left(
\mathcal A_i,\rho_i,\mathcal R_i,
\{(\mathcal I_e,\pi_{i,e})\}_{e\ni i},
\mathcal U_i,\mathrm{Chk}_i
\right).
```

Here $\mathcal A_i$ is the local accessible material algebra, $\rho_i$ is
the local state, $\mathcal R_i$ is the material record algebra, $\mathcal I_e$
is the port or overlap algebra, $\pi_{i,e}$ is the readback map, $\mathcal U_i$
is the update/repair/control family, and $\mathrm{Chk}_i$ is checkpoint data.
The material record algebra contains density, current, energy, pair-correlation,
phase, fluxoid/gauge-field, defect, quasiparticle, aggregate-transport, and
process-history records:

```math
\mathcal R_i^{\rm mat}
=
\mathcal R_i^n\vee\mathcal R_i^J\vee\mathcal R_i^E
\vee\mathcal R_i^\Delta\vee\mathcal R_i^\theta
\vee\mathcal R_i^{\rm flux}\vee\mathcal R_i^{\rm def}
\vee\mathcal R_i^{\rm qp}\vee\mathcal R_i^{\rm agg}
\vee\mathcal R_i^{\rm proc}.
```

The observer-like claim is operational rather than metaphorical only after a
material implementation is supplied: local order and defect state are the
state; electromagnetic, Josephson, phonon, and material exchange are the
ports; gauge-invariant current, phase-difference, strain, and spectroscopy are
readback; metastable structural or defect configurations are records; TDGL
relaxation, vortex motion, annealing, gating, or strain control are update
moves; and synchronized structural, transport, magnetic, and spectroscopic
measurements are checkpoints. Intrinsic material relaxation and an external
laboratory controller must be labeled separately.

Normal-form projection by itself leaves the material law unselected. A
simulator or material proof must declare a source law, either as a finite
quotient law

```math
\mu_{\mathbf x,r,T}(q)
=
Z_{\mathbf x,r,T}^{-1}
m_{\mathbf x,r,T}(q)e^{-S_{\mathbf x,r,T}(q)}
```

Here $S$ is dimensionless; if an energy functional is used instead, the
Boltzmann factor is $e^{-\beta E}$ and the base measure and normalization must
be declared. In the quantum expression below, $\beta=(k_BT)^{-1}$.

or as a quantum density operator on the declared quotient physical algebra:

```math
\rho_{\mathbf x,r,T}
=
Z^{-1}e^{-\beta(H_{\mathbf x,r}-\mu N)}.
```

The source type must be explicit:
`EXTERNAL_FIXED_REFERENCE`, `OPH_NATIVE_QUOTIENT_ENSEMBLE`,
`OPH_VACUUM_BRANCH`, `EXPERIMENTALLY_FITTED_MODEL`, or
`DIAGNOSTIC_BASELINE`.

## Theorem 0: Normal-Form Repair Does Not Select the Material Ensemble

**Statement.** Let $n:Q\to Q$ be an idempotent normal-form map with image
$Q_{\rm nf}\subseteq Q$, and let
$`\mathcal C_Q(\mu)=n_\#\mu`$ be the induced map on probability laws. Then
$\mathcal C_Q$ is idempotent and has many fixed laws. Therefore canonical
repair does not select the material source law.

**Proof.** Since $n$ is a normal-form map, $n\circ n=n$. Hence

```math
\mathcal C_Q^2(\mu)=(n\circ n)_\#\mu=n_\#\mu=\mathcal C_Q(\mu).
```

Every law supported on $Q_{\rm nf}$ is fixed by $`n_\#`$, so repair leaves many
possible material ensembles. A physical branch therefore requires a declared
base measure/action, density operator, MaxEnt ledger, or other source law.
$\square$

The simulator consequence is strict: quotient normal forms classify material
records, but they do not authorize uniform representative sampling unless
representative counting is the declared physical base measure.

## Local-Pair Observable and Continuous-Instability Kernel

**Established theory.** A local pairing tendency, a broken-symmetry order
parameter, and global phase order are different objects. The anomalous
expectation below is an ordered-state quantity; it is not by itself a
gauge-invariant measure of preformed local pairs.

Let $B_{i,a}$ be a local charge $2e$ pair operator in channel $a$. Couple
it to a source field $h_{i,a}$:

```math
H_{\mathbf x}[h]
=
H_{\mathbf x}
-
\sum_{i,a}
\left(h_{i,a}^\ast B_{i,a}+h_{i,a}B_{i,a}^\dagger\right).
```

The finite-temperature partition function is

```math
Z_{\mathbf x}[h]
=
\mathrm{Tr}
\exp[-\beta(H_{\mathbf x}[h]-\mu N)].
```

The induced pair field is

```math
\Delta_{i,a}[h]
=
\frac1\beta
\frac{\partial \log Z_{\mathbf x}[h]}{\partial h_{i,a}^{\ast}}.
```

In a fixed gauge, or in the corresponding neutral effective model, the ordered
thermodynamic/source-limit order parameter is

```math
\Delta_{i,a}
=
\lim_{h\to0^+}\lim_{r\to\infty}\Delta_{i,a}[h].
```

For the local-pair gate, use instead a declared gauge-invariant pair
correlator or pair susceptibility, for example

```math
\chi^{\rm pair}_{ia}(T)
=
\int_0^\beta d\tau\,
\left\langle
B_{i,a}(\tau)B_{i,a}^\dagger(0)
\right\rangle_{\!c}.
```

A nonlocal pair-correlation matrix must include the appropriate gauge
connector between its endpoints.

The operational local-pair scale $T_{\rm pair}$ must state which eigenvalue,
spectral gap, correlation length, or other gauge-invariant observable crosses
which threshold. In a phase-fluctuating or pseudogap regime this observable can
be large while the source-free anomalous expectation remains zero.

After integrating out or jointly minimizing the competing channels $\Xi$, the
effective action near $\Delta=0$ is

```math
\Gamma_{\mathbf x,T}^{\rm eff}[\Delta]
=
\Gamma_{\mathbf x,T}^{\rm eff}[0]
+
\Delta^\dagger K_{\mathbf x,T}^{(2e)}\Delta
+
\Gamma_4[\Delta]
+
\cdots.
```

For a continuous instability of the normal state, define the inverse-pair
susceptibility or quadratic kernel

```math
K_{\mathbf x,T}^{(2e)}
=
\left.
\frac{\partial^2\Gamma_{\mathbf x,T}^{\rm eff}}
{\partial\Delta^\ast\partial\Delta}
\right|_{\Delta=0}.
```

## Conditional Lemma 1: Quadratic Instability at a Continuous Transition

**Model assumptions.** The effective action is analytic near $\Delta=0$, the
transition is continuous, the kernel is Hermitian on the declared pair space,
and higher-order terms stabilize the ordered state. Under these assumptions,
the normal state is locally unstable when

```math
\lambda_{\min}(K_{\mathbf x,T}^{(2e)})<0.
```

At the continuous onset the lowest eigenvalue crosses zero. A nondegenerate
leading eigenvector gives the symmetry of the incipient order within the
declared model space.

**Proof.** Diagonalize $K$:

```math
K v_n=\lambda_n v_n,
\qquad
\Delta=\sum_n a_n v_n.
```

Then

```math
\Gamma^{(2)}[\Delta]=\sum_n \lambda_n |a_n|^2.
```

If every $\lambda_n>0$, $\Delta=0$ is locally stable against infinitesimal
perturbations in this model space. If $\lambda_0<0$, the quadratic action
decreases along $v_0$, and the assumed stabilizing higher-order terms set a
nonzero ordered-state amplitude. This is a local-stability result, not an
if-and-only-if criterion for all superconducting transitions: a first-order
transition can develop a nonzero global minimum while the Hessian at zero is
still positive, and a degenerate leading space requires the nonlinear terms to
select the final symmetry. $\square$

## Reduced-Model Lemma 2: A Repulsive Edge Favors a Sign Change

**Model status.** This is an algebraic consequence of an assumed reduced pair
functional, not an OPH derivation of the microscopic interaction.

**Statement.** Let the dominant pair scattering between patches or pockets be
repulsive:

```math
J_{pq}>0.
```

If the pair functional contains

```math
E_{\rm rep}
=
\sum_{(p,q)}J_{pq}\mathrm{Re}(\Delta_p^\ast\Delta_q),
```

then the minimizing phase relation across a dominant edge is

```math
\theta_q-\theta_p=\pi
\quad (\mathrm{mod}\ 2\pi).
```

**Proof.** Write $\Delta_p=|\Delta_p|e^{i\theta_p}$. For one dominant edge,

```math
J_{pq}\mathrm{Re}(\Delta_p^\ast\Delta_q)
=
J_{pq}|\Delta_p||\Delta_q|\cos(\theta_q-\theta_p).
```

Since $J_{pq}>0$, the minimum occurs at
$\cos(\theta_q-\theta_p)=-1$. Repulsive mismatch repairs by phase reversal.
$\square$

## Reduced-Model Lemma 3: Cuprate-Like Bond Subspace

**Model assumptions.** On a square-lattice cuprate-like branch with strong on-site
repulsion, nearest-neighbor singlet repair, and $C_{4v}$ symmetry, the
dominant bond-pair kernel diagonalizes into

```math
A_{1g}:\quad \Delta_x+\Delta_y,
\qquad
B_{1g}:\quad \Delta_x-\Delta_y.
```

If the repulsive repair coefficient between the $x$ and $y$ bond channels is
positive, the selected leading channel is $B_{1g}$, the
$d_{x^2-y^2}$ type branch.

**Proof.** The on-site Hubbard term suppresses local $s$ wave repair. The
nearest-neighbor exchange supplies bond singlet repair. In the two-dimensional
bond subspace $(\Delta_x,\Delta_y)$, $C_{4v}$ symmetry gives the kernel

```math
K_{\rm bond}
=
\begin{pmatrix}
a & b\\
b & a
\end{pmatrix}.
```

Its eigenvectors are

```math
\Delta_s=(1,1),\qquad \lambda_s=a+b,
```

and

```math
\Delta_d=(1,-1),\qquad \lambda_d=a-b.
```

For repulsive inter-bond repair $b>0$, $\lambda_d<\lambda_s$. If the
pairing instability occurs first in this subspace, the selected normal form is
$B_{1g}$. $\square$

This reproduces the conventional symmetry analysis of an assumed two-channel
kernel. It does not show that a material realizes that kernel or exclude a
first-order transition or a larger competing pair space.

## Reduced-Model Lemma 4: Bipartite Multiband Sign Structure

**Model assumptions.** Let $P$ be the set of Fermi pockets. Suppose the full
real reduced kernel is symmetrized in the declared pocket measure and there is
a bipartite sign assignment

```math
P=P_+\sqcup P_-.
```

such that, after the sign transformation by
$D_\sigma=\operatorname{diag}(\sigma_p)$, every off-diagonal entry of
$D_\sigma K D_\sigma$ is nonpositive and the resulting matrix is
irreducible. Then a nondegenerate lowest eigenmode can be chosen with opposite
signs on the two parts. This is an $s_\pm$-type branch.

**Proof.** Let $\sigma_p=+1$ on $P_+$ and $\sigma_p=-1$ on $P_-$. Write

```math
\Delta_p=\sigma_p\eta_p.
```

By assumption the sign-transformed matrix is an irreducible symmetric
Z-matrix. After subtracting it from a sufficiently large scalar multiple of
the identity, Perron--Frobenius gives a strictly positive eigenvector for the
largest shifted eigenvalue, equivalently the lowest eigenvalue of the original
matrix. Returning to $\Delta=D_\sigma\eta$, the eigenvector changes sign
between $P_+$ and $P_-$. $\square$

Merely having *dominant* bipartite repulsive edges is insufficient if retained
subleading edges frustrate the Z-matrix sign pattern. Non-bipartite
interactions, orbital weights, density-of-states factors, and subleading
channels can change or frustrate this result. A material claim must
derive or fit the full linearized gap kernel and compare its symmetry with
spectroscopy.

## Gauge-Invariant Phase, Fluxoids, and Defect Records

**Established theory.** Let $q=2e$ denote the magnitude of the pair charge.
A gauge-invariant edge record can be stored without choosing a branch of the
angle:

```math
u_{ij}
=
\exp\!\left[
i\left(
\theta_i-\theta_j
-\frac{q}{\hbar}\int_i^j A\cdot dl
\right)
\right].
```

For a closed loop,

```math
U_\gamma=\prod_{(i,j)\in\gamma}u_{ij}
=\exp\!\left(-i\frac{q}{\hbar}\Phi_\gamma\right),
\qquad
\Phi_\gamma=\int_{\Sigma_\gamma}B\cdot dS,
```

provided the phase is single-valued along the sampled contour. $U_\gamma$ need
not equal one in a coherent, current-carrying superconductor. Reducing every
edge modulo $2\pi$ and then subtracting $2\pi n_\gamma$ would erase the integer
winding information, so the winding ledger is stored separately.

In the London regime, away from vortex cores, the fluxoid receipt is

```math
\Phi_\gamma
+\mu_0\oint_\gamma \lambda_L^2\,\mathbf j_s\cdot d\mathbf l
=n_\gamma\Phi_0,
\qquad
\Phi_0=\frac{h}{2e},
\qquad n_\gamma\in\mathbb Z.
```

For a declared contour and measurement uncertainty, a dimensionless fluxoid
residual is therefore

```math
R_\gamma^{\rm flux}
=
\frac{1}{\Phi_0}
\left|
\Phi_\gamma
+\mu_0\oint_\gamma\lambda_L^2\mathbf j_s\cdot d\mathbf l
-n_\gamma\Phi_0
\right|.
```

This is a consistency receipt, not an independent necessary-and-sufficient
gate for superconductivity. Applied fields, persistent currents, grain
boundaries, and resolved vortices are physical records rather than automatic
failures.

**OPH operational proposal.** The third design diagnostic is instead a declared
defect-limited transport predicate. Depending on dimensionality and experiment,
it may use free-vortex density, phase-slip rate, percolating weak-link
probability, and linear resistivity:

```math
\mathcal D(T;J,B,L,\tau_{\rm obs}):
\begin{cases}
\Gamma_{\rm ps}\tau_{\rm obs}<\varepsilon_{\rm ps},\\
\rho_{\rm lin}<\rho_{\max},\\
p_{\rm perc}<p_{\max},
\end{cases}
```

Here $p_{\rm perc}$ is the probability of a percolating resistive obstruction,
not the probability of a superconducting path. Every threshold, unit, field,
current, geometry, loop basis, and estimator must be
frozen before evaluation. In a two-dimensional phase-only model, vortex
unbinding and the BKT stiffness criterion are the relevant standard
description; in three dimensions no universal stiffness threshold is implied.

## Static and Public Superconducting Predicates

**Model assumptions.** For a chosen material family, dimensionality, field,
current, sample geometry, and observation window, define three operational
diagnostics. For example, let $\chi^{\rm pair}(T)$ be the declared finite
pair-susceptibility matrix and define
$A_{\rm pair}(T)=\lambda_{\max}(\chi^{\rm pair}(T))$, with its units and
normalization retained. The local-pair diagnostic is then

```math
\mathcal A(T):A_{\rm pair}(T)>A_{\min}.
```

The phase stiffness is the helicity-modulus response in a declared twist
convention,

```math
\rho_{s,\mu\nu}
=\left.\frac{1}{V}
\frac{\partial^2F(\boldsymbol\varphi)}
{\partial\varphi_\mu\partial\varphi_\nu}
\right|_{\boldsymbol\varphi=0},
\qquad
\Theta(T):\rho_s(T)>\rho_{\min}(T).
```

Here $\boldsymbol\varphi$ denotes the uniform twist density used in this
formula. The units of $\rho_s$ and $\rho_{\min}$ must be stated for that convention. In
a strictly two-dimensional phase-only model, the BKT transition instead uses
the standard universal-jump relation

```math
k_BT_{\rm BKT}=\frac{\pi}{2}\rho_s(T_{\rm BKT}^{-}),
```

when $\rho_s$ is expressed as an energy. A three-dimensional material has no
universal threshold of this form. The defect-limited transport diagnostic is
$\mathcal D(T;J,B,L,\tau_{\rm obs})$ from the preceding section.

The conjunction

```math
\mathsf{SC}_{\rm model}(T)
=
\mathcal A(T)\wedge\Theta(T)\wedge\mathcal D(T)
```

is a conditional ranking predicate for the declared model. It is not claimed
to be necessary and sufficient for every superconducting phase.

The public evidence predicate separately requires controlled transport and a
Boolean bulk magnetic/thermodynamic pass. Define the magnetic pass only after
freezing field, geometry, demagnetization correction, zero-field-cooled and
field-cooled protocol, and thresholds on quantities such as penetration depth,
susceptibility, internal field, magnetization, and shielding fraction:

```math
\mathcal M_{\rm bulk}(T)
=
\operatorname{Pass}
\left(
\lambda_L,\chi_{\rm dia},B_{\rm in}/B_{\rm ext},M,
f_{\rm shield};\ \text{declared controls}
\right).
```

Then

```math
\mathsf{SC}_{\rm verified}(T,\mathbf x,L,\tau_{\rm obs})
=
\mathsf{SC}_{\rm model}(T)
\wedge\mathcal R_{\rm transport}(T,L,\tau_{\rm obs})
\wedge\mathcal M_{\rm bulk}(T).
```

The highest verified temperature is an evidence boundary, not a redefinition
of the material's physical transition:

```math
T_{\rm verified}(\mathbf x,L,\tau_{\rm obs})
=
\sup\{T:\mathsf{SC}_{\rm verified}(T,\mathbf x,L,\tau_{\rm obs})\}.
```

## Conditional Bottleneck-Score Identity

Define the individual operational thresholds by

```math
T_{\rm pair}=\sup\{T:\mathcal A(T)\},\quad
T_{\rm phase}=\sup\{T:\Theta(T)\},\quad
T_{\rm def}=\sup\{T:\mathcal D(T)\}.
```

If all three pass sets are downward-closed over the declared temperature
interval, then the model score

```math
T_{\rm bottleneck}(\mathbf x)
=
\sup\{T:\mathsf{SC}_{\rm model}(T)\}
=\min(T_{\rm pair},T_{\rm phase},T_{\rm def}).
```

**Proof.** Under downward closure, each pass set is an interval ending at its
declared threshold. Their intersection ends at the smallest threshold.
$\square$

This is a set identity for the chosen diagnostics. It neither derives those
thresholds nor proves that the resulting score equals physical $T_c$.

## Defect Repair, Escape, and Transport Records

The defect sector is the quotient-visible sector of vortex records,
phase-slip records, branch cuts, grain-boundary obstructions, and fluxoid
consistency residuals:

```math
D_{\mathbf x,r}
=
\{
v_\gamma,n_\gamma,\text{vortex records},\text{phase-slip records},
\text{branch cuts},\text{grain-boundary obstructions},
\text{fluxoid residuals}
\}/\Gamma_{\mathbf x,r}.
```

**Model status.** The following clearance operator is an optional finite-state
kinetic model, not a material-independent law. A lowest smooth
phase eigenvalue can vanish with system size because of harmless Goldstone or
twist modes; that alone is not a superconducting failure.

Let $\mu_D$ be the induced defect-sector law and let $E_C$ be conditional
orthogonal repair projections on $L^2(D_{\mathbf x,r},\mu_D)$. Let $c_C\ge0$
have units of inverse time and be derived from a declared kinetic model or
measurement. Define

```math
L_{\rm clear}^{\rm def}
=
\sum_C c_C(I-E_C),
```

and

```math
\Delta_{\rm clear}^{\rm def}
=
\inf_{f\perp\ker L_{\rm clear}^{\rm def}}
\frac{
\langle f,L_{\rm clear}^{\rm def}f\rangle_{\mu_D}
}{
\langle f,f\rangle_{\mu_D}
}.
```

Let $S_{\rm pass}=\{q:\mathcal D(q)\text{ passes}\}$. For a time-homogeneous
continuous-time defect generator $K_D$, an upper bound on the defect escape or
nucleation rate is

```math
\lambda_{\rm esc}^{\rm def}
=
\sup_{q\in S_{\rm pass}}
\sum_{q'\notin S_{\rm pass}}K_D(q,q').
```

For a discrete simulator kernel $P_D$, use

```math
p_{\rm esc}^{\rm step}
=
\sup_{q\in S_{\rm pass}}P_D(q,S_{\rm pass}^c).
```

## Mathematical Defect-Model Decay Bound

**Statement.** If $D_{\mathbf x,r}$ is finite and
$L_{\rm clear}^{\rm def}$ is self-adjoint and nonnegative, then for
$f\perp\ker L_{\rm clear}^{\rm def}$,

```math
\langle f,e^{-tL_{\rm clear}^{\rm def}}f\rangle_{\mu_D}
\le
e^{-t\Delta_{\rm clear}^{\rm def}}\langle f,f\rangle_{\mu_D}.
```

**Proof.** Diagonalize the finite self-adjoint generator. All components
orthogonal to the kernel have eigenvalue at least
$\Delta_{\rm clear}^{\rm def}$. Expanding $f$ in that eigenbasis gives the
bound term by term. $\square$

## Conditional Finite-Window Stability Check

**Statement.** Within the finite jump model above, a large clearance gap alone
does not guarantee a quiet observation window. A sufficient declared check is

```math
\lambda_{\rm esc}^{\rm def}\tau_{\rm obs}\ll1
```

and

```math
\tau_{\rm clear}:=(\Delta_{\rm clear}^{\rm def})^{-1}
\ll\tau_{\rm allow},
```

where $\tau_{\rm allow}$ is the maximum accepted clearance time. In a
time-homogeneous finite jump process, the probability of at least one escape
from $S_{\rm pass}$ during $\tau_{\rm obs}$ is bounded by
$1-e^{-\lambda_{\rm esc}^{\rm def}\tau_{\rm obs}}$. Once an unresolved defect
is present, the clearance-model lemma controls its decay. If escape is large,
defects are repeatedly created regardless of how quickly individual defects
clear. These are sufficient model checks, not necessary conditions for every
laboratory superconductor. $\square$

Define a transport mismatch functional

```math
\Phi_{\rm tr}(q)
=
\Phi_J(q)+\Phi_V(q)+\Phi_{\rm slip}(q)+\Phi_{\rm heat}(q),
```

where $\Phi_V$ includes voltage-lead mismatch and $\Phi_{\rm slip}$ counts
unresolved phase-slip events in the observation window.

## Evidence Rule: Zero Resistance Alone Is Not a Bulk Claim

**Statement.** A branch with apparent zero resistance but no
Meissner/London/bulk-shielding receipt is classified as
`TRANSPORT_ARTIFACT`, `FILAMENTARY_PAIR_PATH`, or `INCOMPLETE_SC_CLAIM`, not
as a verified bulk OPH superconductivity claim.

**Proof.** Transport records and magnetic/bulk gauge-field repair records are
distinct public readouts. A perfect conductor can preserve magnetic flux, and
a filament can short voltage contacts without producing a bulk phase normal
form. The public evidence predicate therefore requires a controlled bulk
magnetic or thermodynamic readout in addition to transport. Shielding alone
must not be called Meissner expulsion without the declared cooling, field, and
demagnetization controls. $\square$

## Inverse-Design Functional

Let $\mathfrak X$ be the candidate family of compositions, structures, and
processing histories. Let the penalties be dimensionless probabilities or
normalized costs, and let every $\lambda$ carry units of kelvin in this
convention. Define the conditional utility

```math
U_{\rm OPH}(\mathbf x)
=
\mathbb E[T_{\rm bottleneck}(\mathbf x)]
-
\lambda_{\rm inst}\mathcal P_{\rm inst}(\mathbf x)
-
\lambda_{\rm dis}\mathcal P_{\rm dis}(\mathbf x)
-
\lambda_{\rm tox}\mathcal P_{\rm tox}(\mathbf x)
-
\lambda_{\rm synth}\mathcal P_{\rm synth}(\mathbf x).
```

The penalties record instability, disorder, toxicity or unusability, and
synthesis failure. Their weights, uncertainty model, and training data must be
frozen before a held-out ranking test.

## Mathematical Lemma: Existence of an Optimizer

**Statement.** If $\mathfrak X$ is finite or compact after regulator,
stability, and synthesis filters, and if $U_{\rm OPH}$ is bounded above and
upper semicontinuous, then a maximizer exists:

```math
\mathbf x_\star\in\mathrm{arg\,max}_{\mathbf x\in\mathfrak X}U_{\rm OPH}(\mathbf x).
```

$\mathbf x_\star$ is the highest-scoring target under the declared utility, not
an assurance that it is synthesizable or superconducting.

**Proof.** A finite filtered candidate set has a maximizer by enumeration. In
the compact case, upper semicontinuity and boundedness give a maximizer by the
Weierstrass theorem. $\square$

This lemma establishes existence only. It supplies neither an algorithm nor a
material prediction; predictive value requires prospective comparison with
declared conventional materials-design baselines.

## Conditional Lemma: Metastable Processing and Pressure Quench

**Model assumptions.** Let $K^P_T$, $\rho_s^P(T)$, and $D_{\rm def}^P(T)$ be
the identified pair kernel, stiffness, and declared defect score of a
high-pressure or
processed phase. Suppose a retained ambient phase $R$ satisfies

```math
\|K^R_T-K^P_T\|\le \epsilon_K,
\qquad
|\rho_s^R(T)-\rho_s^P(T)|\le\epsilon_\rho,
\qquad
|D_{\rm def}^R(T)-D_{\rm def}^P(T)|\le\epsilon_D,
```

after a declared identification of the two kernel spaces. Suppose these errors
are smaller than the three operational margins below
$T_{\rm bottleneck}^P-\delta T$, and suppose a measured or modeled lifetime
obeys

```math
\tau_{\rm life}\simeq\nu_0^{-1}
\exp(\Delta G^\ddagger/k_BT)\ge\tau_{\rm use}.
```

Then the retained phase passes the same model diagnostics up to that
temperature:

```math
T_{\rm bottleneck}^R\ge T_{\rm bottleneck}^P-\delta T.
```

**Proof.** The amplitude condition is stable under small Hermitian perturbations
by Weyl's eigenvalue inequality. The phase-stiffness condition is stable under
$\epsilon_\rho$. The declared defect diagnostic is stable under $\epsilon_D$
when the perturbation is smaller than its margin. The lifetime premise makes
the retained state long-lived on the use timescale. $\square$

This is a robustness statement conditional on measured closeness and lifetime;
it does not predict that pressure quenching will create such a phase.

## Design Hypothesis: Heterostructure Coupling of Pairing and Stiffness

Suppose layer $A$ has high $T_{\rm pair}$ and low stiffness,
while layer $B$ has high stiffness and compatible pair symmetry. If the
Josephson coupling

```math
E_J=-\sum_{\langle A,B\rangle}J_{AB}\cos(\theta_A-\theta_B-\phi_{AB})
```

adds stiffness without frustrating the phase relation or suppressing the pair
amplitude through inverse proximity, then a coupled calculation may raise the
declared bottleneck score

```math
\min(T_{\rm pair},T_{\rm phase},T_{\rm def})
```

above either weak isolated bottleneck.

The helicity modulus of the coupled phase functional can include positive
contributions from intralayer stiffness and Josephson locking. Whether the
benefit survives inverse proximity, competing order, interface disorder, and
phase frustration must be computed in a coupled GL or microscopic model and
tested experimentally. This is a design hypothesis, not a theorem.

## Build Surface

The build surface separates established devices, record-class cuprate
processing, and hydride metastability.

### Established HTS Devices

REBCO and YBCO-family high-temperature superconductors are established
engineering materials. REBCO coated conductors operate in some
liquid-nitrogen-range applications, with usable current depending on field,
orientation, geometry, and temperature, and appear in long-length tape and
magnet work. In OPH terms, this
tier has a closed material receipt for known families: source material, cooling
boundary, transport readout, magnetic response, and engineering load.

The build instruction is ordinary materials practice: use a known REBCO or YBCO
conductor, cool through its transition, and verify zero resistance plus magnetic
response under the declared sample geometry and contacts.

### Pressure-Quenched Multilayer Cuprates

One recent empirical case study for the retention audit is pressure-quenched
Hg-1223,

```math
\mathrm{HgBa_2Ca_2Cu_3O_{8+\delta}}.
```

This is a reported materials result, not a prediction derived by OPH. The 2026
study reports an ambient-pressure transition up to 151 K after pressure
quenching, a bulk fraction of about 78% in the reported magnetic measurement,
retention for at least three days at 77 K, and degradation after heating above
200 K. Independent replication and scalable processing remain separate
receipts. The OPH interpretation is a retained-observable problem:

```math
\mathcal O_{\rm high\ pressure}
\longrightarrow
\mathcal O_{\rm metastable\ ambient}.
```

The pressure protocol tries to retain the high-pressure structural and
superconducting observables while preserving phase stiffness and controlled
defect-limited transport after pressure release. The professional-lab
constraints are part of the receipt: mercury
chemistry, oxygen stoichiometry, pressure path, quench path, metastability
lifetime, transport contacts, magnetic response, and structural readout.

This lane belongs in a serious materials lab.

### High-Pressure Hydrides

Hydrides occupy the high-transition-temperature, high-pressure lane. LaH10 has
a reported superconducting transition around $250\,\mathrm K$ near
$170\,\mathrm{GPa}$. That transition does not independently measure the
local-pair scale used in this note. In the present audit language,

```math
T_c^{\rm reported}\ \text{is high under pressure},
\qquad
\text{ambient stability and retention remain open design constraints}.
```

The actionable hydride program searches for metastable descendants whose
retained structure keeps the high-pressure pair kernel after decompression.

### Claim Boundary for Build Claims

A build claim requires these receipts:

1. structural receipt: phase, stoichiometry, strain, pressure or quench history;
2. transport receipt: zero resistance with contact and geometry controls;
3. magnetic receipt: controlled Meissner expulsion or a shielding response
   explicitly labeled as such, with geometry and demagnetization correction;
4. thermodynamic or equivalent bulk receipt where the claim needs bulk
   promotion;
5. stiffness and gap-symmetry receipt for the OPH branch claim;
6. replication receipt by an independent lab for record-class or extraordinary
   claims.

Dias-associated room-temperature claims and the LK-99 episode are the cautionary
boundary: social momentum and partial transport anomalies do not promote a
material branch. The OPH receipt requires the whole packet.

### Illustrative Research Targets

No quantitative OPH ranking is reported here. Illustrative research lanes are:

1. pressure-quenched, oxygen/strain-optimized Hg-1223 or related multilayer
   cuprates;
2. cuprate heterostructures that combine high-pairing layers with high-stiffness
   layers;
3. hydride metastability programs that try to retain a high-pressure pair kernel
   after decompression.

The compact conditional ranking rule is

```math
T_{\rm bottleneck}
=
\min(T_{\rm pair},T_{\rm phase},T_{\rm def}).
```

Pairing strength alone is not a sufficient objective. The score targets the
weakest declared diagnostic and must not be relabeled physical $T_c$ until
validated.

## Build-Discovery Protocol

1. Choose a lane: cuprate-like, pnictide/chalcogenide, nickelate, or hydride.
2. Generate candidate structures and processing histories $\mathbf x$.
3. Reject unstable, unsynthesizable, unsafe, or inaccessible candidates.
4. Compute a declared local-pair observable and, for a continuous-transition
   model, $K_{\mathbf x,T}^{(2e)}$ and its leading eigenmode.
5. Compute the operational local-pair threshold $T_{\rm pair}$.
6. Compute phase stiffness and $T_{\rm phase}$.
7. Compute fluxoid consistency receipts, defect dynamics, controlled transport,
   and $T_{\rm def}$.
8. Rank by
   $T_{\rm bottleneck}=\min(T_{\rm pair},T_{\rm phase},T_{\rm def})$.
9. Prefer candidates where processing raises the bottleneck diagnostic while
   preserving the other diagnostics.
10. Synthesize or process top candidates under professional materials controls.
11. Verify with structure, zero resistance, Meissner, stiffness, and gap-symmetry
    receipts.
12. Feed failures back into the pair model, $\rho_s$, and $\mathcal D$.

## Verifier Receipts

The public material receipt must include:

- structural receipt: composition, phase, lattice, strain, pressure or quench
  history, disorder, and reproducible processing record;
- normal-state receipt: carrier density, resistivity, competing orders, and
  relevant magnetic or orbital response;
- pairing receipt: gap symmetry, pair channel, tunneling or spectroscopy record,
  and competing-order separation;
- superconductivity receipt: zero resistance, Meissner response, critical
  current, and field dependence;
- stiffness receipt: penetration depth, phase stiffness, or helicity-modulus
  proxy with units and controls;
- fluxoid/defect receipt: winding, flux, supercurrent, grain-boundary, vortex,
  phase-slip, disorder, and controlled-transport data;
- process-retention receipt: lifetime and relaxation margin for metastable
  phases.

## Falsifiers and Audit Failures

Implementation failures are not by themselves falsifiers of OPH. Missing bulk
evidence, hidden processing history, dimensionally inconsistent scores, or a
ranking that ignores declared diagnostics fails the public claim or protocol.

A distinguishing predictive test must instead freeze the candidate set, source
law, pair observable, stiffness and defect estimators, thresholds, weights,
uncertainties, and conventional baselines before inspecting held-out materials.
The material-family model is falsified on that domain if, under those controlled
assumptions:

- the predicted continuous-onset symmetry disagrees with spectroscopy beyond
  the preregistered uncertainty;
- predicted transition or ranking errors exceed the preregistered tolerance;
- the OPH score fails to improve the held-out ranking or calibrated likelihood
  relative to the declared Eliashberg, DFT/many-body, GL/XY, or materials-model
  controls; or
- a claimed retention law fails its preregistered lifetime and replication
  tests.

Until such a frozen comparison exists, the note is falsifiable as a specific
model contract but does not make a unique OPH material prediction.

## Scope

This note closes a conceptual audit specification, not the high-temperature $`T_c`$ mechanism
or materials-discovery problem. The remaining research protocol is to derive or
fit a material source law, compute local-pair, stiffness, defect/transport, and
retention observables with units and uncertainty, freeze a prospective ranking,
synthesize candidates, verify bulk receipts, and compare held-out performance
with standard baselines.

## References

- V. J. Emery and S. A. Kivelson, "Importance of phase fluctuations in
  superconductors with small superfluid density", Nature, 1995.
  https://doi.org/10.1038/374434a0
- D. R. Nelson and J. M. Kosterlitz, "Universal Jump in the Superfluid Density
  of Two-Dimensional Superfluids", Physical Review Letters, 1977.
  https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.39.1201
- B. S. Deaver Jr. and W. M. Fairbank, "Experimental Evidence for Quantized
  Flux in Superconducting Cylinders", Physical Review Letters, 1961.
  https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.7.43
- Hualei Sun et al., "Signatures of superconductivity near 80 K in a nickelate
  under high pressure", Nature, 2023.
  https://www.nature.com/articles/s41586-023-06408-7
- Kai Wang et al., "Advances in second-generation high-temperature
  superconducting tapes and their applications in high-field magnets", Soft
  Science, 2022. https://www.oaepublish.com/articles/ss.2022.10
- Liangzi Deng et al., "Ambient-pressure 151-K superconductivity in
  HgBa2Ca2Cu3O8+delta via pressure quench", PNAS, 2026.
  https://www.pnas.org/doi/10.1073/pnas.2536178123
- A. P. Drozdov et al., "Superconductivity at 250 K in lanthanum hydride under
  high pressures", Nature, 2019. https://www.nature.com/articles/s41586-019-1201-8
- Dan Garisto, "Superconductivity scandal: the inside story of deception in a
  rising star's physics lab", Nature, 2024.
  https://www.nature.com/articles/d41586-024-00716-2
- "Perception and argumentation in the LK-99 superconductivity controversy",
  Scientific Reports, 2025. https://www.nature.com/articles/s41598-025-98554-3
