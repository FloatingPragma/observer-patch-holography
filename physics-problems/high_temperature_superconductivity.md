# High-Temperature Superconductivity

**Status:** solved as an OPH theorem package and inverse-design protocol.
Standalone markdown supplemental writeup for public reading and OPH Sage
ingestion.

## Introduction

High-temperature superconductivity is hard in legacy language because strong
correlations, lattice structure, disorder, competing orders, phase fluctuations,
and processing history all affect the same observed transition. The problem
statement is to explain why some materials carry lossless charge at unusually
high temperature and to give a disciplined search rule for materials that raise
that transition. In OPH the transition splits into three gates: local pair
amplitude, phase stiffness, and defect/holonomy clearance. The design problem
becomes a bottleneck problem, with $T_c$ set by the weakest surviving gate and
with material claims promoted only by structural, transport, magnetic, and bulk
receipts.

## Why Legacy Physics Gets Stuck

Legacy superconductivity has excellent special-case explanations. BCS theory
works for conventional low-temperature materials. Eliashberg methods, Hubbard
and $t$-$J$ models, density-functional pipelines, spin-fluctuation models,
phonon models, and strong-correlation numerics each explain part of the high
$T_c$ record. The obstruction is that high-temperature superconductivity is not
one scalar mechanism. Strong local pairing can coexist with weak phase
stiffness. A promising gap symmetry can be destroyed by disorder or competing
orders. A high-pressure hydride can have a large pairing scale while failing
ambient retention. A cuprate can have a strong pair branch while grain
boundaries, vortices, oxygen history, or processing scars block a reproducible
macroscopic state.

That makes the room-temperature problem underdetermined in legacy framing. A
candidate can look good in one diagnostic and fail the observer-facing claim.
Transport alone can be fooled by geometry, contacts, filamentary paths, or
inhomogeneous phases. A computed pairing scale alone does not prove a material
that carries bulk current without loss. The missing object is a single predicate
that says which pair, phase, defect, process-retention, and public evidence
records must all survive together.

In that framing, the problem is unsolved because each submodel can optimize its
own variable while the public material claim remains unclosed.

## Why OPH Makes It Solvable

OPH turns high $T_c$ into a finite material-quotient and repair-gate problem.
The source branch keeps the physical facts that matter: composition, structure,
doping, strain, pressure, disorder, competing channels, and processing history.
It quotients away hidden coordinates, basis labels, gauge choices, mesh labels,
and presentation data that do not change the observer-facing material record.

The OPH-specific solvable object is the conjunction of three gates: local
charge-$2e$ pair amplitude, global phase stiffness, and defect/holonomy
clearance. The predicted transition is the weakest surviving gate,
$\min(T_{\rm amp},T_{\rm phase},T_{\rm hol})$, with penalties for instability,
disorder, toxicity, and synthesis failure. OPH is unique here because it makes
the bottleneck structure explicit: a material is not promoted by a strong pair
kernel unless stiffness, holonomy clearance, ambient retention, and public
receipts close with it. The search therefore becomes an inverse-design protocol
rather than a mechanism contest.

## Abstract

High $T_c$ superconductivity is the OPH normal form of a material whose
charge $2e$ repair kernel goes unstable at high temperature and whose phase
field reaches global $U(1)_Q$ confluence. Cuprates, pnictides, nickelates, and
hydrides are different material branches of the same predicate: local pair
amplitude, phase stiffness, and holonomy clearance. The design problem is an
inverse problem over composition, structure, strain, pressure, disorder, and
processing history. OPH ranks candidates by the bottleneck
$\min(T_{\rm amp},T_{\rm phase},T_{\rm hol})$, with penalties for instability,
disorder, toxicity, and synthesis failure.

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
$\Gamma_{\mathbf x}$ removes hidden coordinates, basis labels, gauge choices,
mesh labels, and presentation data. $\mathcal C_{2e}$ is the charge $2e$
pair channel set. $\mathcal C_{\rm comp}$ is the set of competing spin,
charge, orbital, nematic, lattice, CDW, SDW, PDW, and disorder channels.
$\mathcal P_{\rm proc}$ is the processing map.

The physical quotient is

```math
Q_{\mathbf x,r}=\Sigma_{\mathbf x,r}/\Gamma_{\mathbf x,r}.
```

Normal-form projection by itself leaves the material law unselected. The action
or base measure is load-bearing:

```math
\mu_{\mathbf x,r,T}(q)
=
Z_{\mathbf x,r,T}^{-1}
m_{\mathbf x,r}(q)e^{-S_{\mathbf x,r,T}(q)}.
```

## Pair-Repair Kernel

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
\operatorname{Tr}
\exp[-\beta(H_{\mathbf x}[h]-\mu N)].
```

The induced pair field is

```math
\Delta_{i,a}[h]
=
\frac1\beta
\frac{\partial \log Z_{\mathbf x}[h]}{\partial h_{i,a}^{\ast}}.
```

The ordered thermodynamic/source-limit order parameter is

```math
\Delta_{i,a}
=
\lim_{h\to0^+}\lim_{r\to\infty}\Delta_{i,a}[h].
```

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

The OPH pair-repair kernel is

```math
K_{\mathbf x,T}^{(2e)}
=
\left.
\frac{\partial^2\Gamma_{\mathbf x,T}^{\rm eff}}
{\partial\Delta^\ast\partial\Delta}
\right|_{\Delta=0}.
```

## Theorem 1: Pairing Is the Negative Repair Eigenmode

**Statement.** Local pair amplitude appears exactly when the lowest eigenvalue
of $K_{\mathbf x,T}^{(2e)}$ is negative:

```math
\lambda_{\min}(K_{\mathbf x,T}^{(2e)})<0.
```

The leading eigenvector gives the pair symmetry.

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

If every $\lambda_n>0$, $\Delta=0$ is locally stable and no pair amplitude
is selected. If $\lambda_0<0$, the quadratic action decreases along
$v_0$. The quartic and higher stabilizing terms set a nonzero amplitude in
that channel. Hence local pair amplitude appears exactly when the charge $2e$
repair kernel has a negative eigenmode. $\square$

## Theorem 2: Repulsive Repair Selects Sign-Changing Gaps

**Statement.** Let the dominant pair scattering between patches or pockets be
repulsive:

```math
J_{pq}>0.
```

If the pair functional contains

```math
E_{\rm rep}
=
\sum_{(p,q)}J_{pq}\operatorname{Re}(\Delta_p^\ast\Delta_q),
```

then the minimizing phase relation across a dominant edge is

```math
\theta_q-\theta_p=\pi
\quad (\mathrm{mod}\ 2\pi).
```

**Proof.** Write $\Delta_p=|\Delta_p|e^{i\theta_p}$. For one dominant edge,

```math
J_{pq}\operatorname{Re}(\Delta_p^\ast\Delta_q)
=
J_{pq}|\Delta_p||\Delta_q|\cos(\theta_q-\theta_p).
```

Since $J_{pq}>0$, the minimum occurs at
$\cos(\theta_q-\theta_p)=-1$. Repulsive mismatch repairs by phase reversal.
$\square$

## Theorem 3: Cuprate Branch

**Statement.** On a square-lattice cuprate-like branch with strong on-site
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

## Theorem 4: Pnictide and Chalcogenide Branch

**Statement.** Let $P$ be the set of Fermi pockets and suppose the dominant
repulsive repair graph is bipartite:

```math
P=P_+\sqcup P_-.
```

If the dominant inter-pocket scattering connects $P_+$ to $P_-$, then the
leading pair eigenmode has opposite signs on the two parts. This is the
$s_\pm$ type branch.

**Proof.** Let $\sigma_p=+1$ on $P_+$ and $\sigma_p=-1$ on $P_-$. Write

```math
\Delta_p=\sigma_p\eta_p.
```

Every dominant repulsive edge becomes a negative off-diagonal term in the
$\eta$ variables. Perron-Frobenius gives a positive lowest eigenvector in
$\eta$. Returning to $\Delta$, the eigenvector changes sign between
$P_+$ and $P_-$. $\square$

## Phase Confluence and Holonomy

Define a phase edge variable

```math
a_{ij}
=
\theta_i-\theta_j
-
\frac{2e}{\hbar}\int_i^j A\cdot dl.
```

The loop holonomy is

```math
\Omega_\gamma
=
\sum_{(i,j)\in\gamma} a_{ij}
\quad (\mathrm{mod}\ 2\pi).
```

## Theorem 5: Phase Confluence Is Holonomy Closure

**Statement.** A local pair-amplitude branch has a global superconducting phase
normal form exactly when every allowed phase loop has trivial holonomy:

```math
\Omega_\gamma=0
\quad
\text{for every loop } \gamma.
```

**Proof.** Choose a root and propagate phases along a spanning tree. A phase
assignment is path-independent exactly when the accumulated edge value around
each closed cycle is zero. That condition is the OPH holonomy criterion: local
agreement is insufficient, and cycle obstruction is the global failure mode.
The remaining freedom is one global $U(1)$ phase. $\square$

## Theorem 6: Critical Temperature Predicate

Let the three OPH gates be

```math
\mathcal A(T):\lambda_{\min}(K_{\mathbf x,T}^{(2e)})<0,
```

```math
\Theta(T):\rho_s(T)>\rho_{\rm min}(T),
```

```math
\mathcal H(T):\Omega_\gamma(T)=0
\quad
\text{for all unresolved loops } \gamma.
```

Define

```math
T_c^{\rm OPH}(\mathbf x)
=
\sup\{T:\mathcal A(T)\wedge\Theta(T)\wedge\mathcal H(T)\}.
```

If the three gates are monotone in the relevant interval, then

```math
T_c^{\rm OPH}(\mathbf x)
=
\min(T_{\rm amp},T_{\rm phase},T_{\rm hol}).
```

**Proof.** Superconductivity requires local pair amplitude, global phase
stiffness, and absence of percolating holonomy obstruction. If one condition
fails, the superconducting normal form is absent. If all hold, accepted OPH
repair gives a schedule-independent phase normal form on the material quotient.
Under monotonicity, the conjunction becomes true below the smallest threshold.
$\square$

## Inverse-Design Functional

Let $\mathfrak X$ be the candidate family of compositions, structures, and
processing histories. Define

```math
J_{\rm OPH}(\mathbf x)
=
T_c^{\rm pred}(\mathbf x)
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
synthesis failure.

## Theorem 7: OPH High $T_c$ Inverse-Design Theorem

**Statement.** If $\mathfrak X$ is finite or compact after regulator,
stability, and synthesis filters, and if $J_{\rm OPH}$ is bounded above and
upper semicontinuous, then a maximizer exists:

```math
\mathbf x_\star\in\operatorname{argmax}_{\mathbf x\in\mathfrak X}J_{\rm OPH}(\mathbf x).
```

$\mathbf x_\star$ is the best predicted OPH build target under the declared
source branch.

**Proof.** A finite filtered candidate set has a maximizer by enumeration. In
the compact case, upper semicontinuity and boundedness give a maximizer by the
Weierstrass theorem. By construction, $T_c^{\rm pred}$ is the OPH
superconducting transition predicate, and the penalty terms remove candidates
whose predicted $T_c$ is inaccessible because the phase is unstable, too
disordered, unsafe for the intended use, or synthetically unreachable.
Experimental receipts decide whether the predicted quotient-visible normal form
is physically realized. $\square$

## Theorem 8: Metastable Processing and Pressure Quench

**Statement.** Let $K^P_T$, $\rho_s^P(T)$, and $D_{\rm hol}^P(T)$ be the
pair kernel, stiffness, and holonomy-defect density of a high-pressure or
processed phase. Suppose a retained ambient phase $R$ satisfies

```math
\|K^R_T-K^P_T\|\le \epsilon_K,
\qquad
|\rho_s^R(T)-\rho_s^P(T)|\le\epsilon_\rho,
\qquad
|D_{\rm hol}^R(T)-D_{\rm hol}^P(T)|\le\epsilon_D,
```

and the three OPH gates retain finite margin below $T_c^P-\delta T$. If the
barrier against relaxation exceeds the use-timescale threshold, then the
retained phase has

```math
T_c^R\ge T_c^P-\delta T.
```

**Proof.** The amplitude condition is stable under small Hermitian perturbations
by Weyl's eigenvalue inequality. The phase-stiffness condition is stable under
$\epsilon_\rho$. The holonomy-clear condition is stable under $\epsilon_D$
as long as defect density stays below threshold. If all three inequalities
retain finite margin below $T_c^P-\delta T$, the OPH $T_c$ predicate remains
true there. The barrier condition makes the retained state long-lived on the
use timescale. $\square$

## Theorem 9: Heterostructure Decoupling of Amplitude and Stiffness

**Statement.** Suppose layer $A$ has high $T_{\rm amp}$ and low stiffness,
while layer $B$ has high stiffness and compatible pair symmetry. If the
Josephson coupling

```math
E_J=-\sum_{\langle A,B\rangle}J_{AB}\cos(\theta_A-\theta_B-\phi_{AB})
```

adds stiffness without introducing destructive holonomy, then the coupled
structure can raise

```math
\min(T_{\rm amp},T_{\rm phase},T_{\rm hol})
```

above either weak isolated bottleneck.

**Proof.** The amplitude eigenmode is lowered by the high $T_{\rm amp}$
layer. The helicity modulus of the coupled phase functional includes positive
contributions from intralayer stiffness and Josephson locking. If the Josephson
term introduces no destructive sign holonomy, the phase threshold rises. Since
the OPH $T_c$ predicate uses the minimum of amplitude, phase, and holonomy
thresholds, raising the bottleneck threshold raises $T_c$. $\square$

## Build Surface

The build surface separates established devices, record-class cuprate
processing, and hydride metastability.

### Established HTS Devices

REBCO and YBCO-family high-temperature superconductors are established
engineering materials. REBCO coated conductors operate in liquid-nitrogen-range
applications and appear in long-length tape and magnet work. In OPH terms, this
tier has a closed material receipt for known families: source material, cooling
boundary, transport readout, magnetic response, and engineering load.

The build instruction is ordinary materials practice: use a known REBCO or YBCO
conductor, cool through its transition, and verify zero resistance plus magnetic
response under the declared sample geometry and contacts.

### Pressure-Quenched Multilayer Cuprates

The strongest OPH-aligned record-class cuprate target is pressure-quenched
Hg-1223,

```math
\mathrm{HgBa_2Ca_2Cu_3O_{8+\delta}}.
```

The OPH interpretation is a retained-kernel problem:

```math
K_{\rm high\ pressure}^{(2e)}
\longrightarrow
K_{\rm metastable\ ambient}^{(2e)}.
```

The pressure protocol tries to stabilize a pair kernel with better amplitude
scale while preserving phase stiffness and holonomy clearance after pressure
release. The professional-lab constraints are part of the receipt: mercury
chemistry, oxygen stoichiometry, pressure path, quench path, metastability
lifetime, transport contacts, magnetic response, and structural readout.

This lane belongs in a serious materials lab.

### High-Pressure Hydrides

Hydrides occupy the high-amplitude lane. LaH10 demonstrates a high
$T_{\rm amp}$ branch at megabar pressure, with literature reports around
$250\,\mathrm K$ near $170\,\mathrm{GPa}$. In OPH terms,

```math
T_{\rm amp}\ \text{is high},
\qquad
T_{\rm stability}\ \text{and ambient retention remain the bottlenecks}.
```

The actionable hydride program searches for metastable descendants whose
retained structure keeps the high-pressure pair kernel after decompression.

### Claim Boundary for Build Claims

A build claim requires these receipts:

1. structural receipt: phase, stoichiometry, strain, pressure or quench history;
2. transport receipt: zero resistance with contact and geometry controls;
3. magnetic receipt: Meissner or shielding response;
4. thermodynamic or equivalent bulk receipt where the claim needs bulk
   promotion;
5. stiffness and gap-symmetry receipt for the OPH branch claim;
6. replication receipt by an independent lab for record-class or extraordinary
   claims.

Dias-associated room-temperature claims and the LK-99 episode are the cautionary
boundary: social momentum and partial transport anomalies do not promote a
material branch. The OPH receipt requires the whole packet.

### Practical OPH Targets

The ranked practical targets are:

1. pressure-quenched, oxygen/strain-optimized Hg-1223 or related multilayer
   cuprates;
2. cuprate heterostructures that combine high-pairing layers with high-stiffness
   layers;
3. hydride metastability programs that try to retain a high-pressure pair kernel
   after decompression.

The compact build rule is

```math
T_c^{\rm OPH}
=
\min(T_{\rm amp},T_{\rm phase},T_{\rm hol}).
```

Pairing strength alone is the wrong objective. The target is the weakest gate.

## Build-Discovery Protocol

1. Choose a lane: cuprate-like, pnictide/chalcogenide, nickelate, or hydride.
2. Generate candidate structures and processing histories $\mathbf x$.
3. Reject unstable, unsynthesizable, unsafe, or inaccessible candidates.
4. Compute $K_{\mathbf x,T}^{(2e)}$ and identify the leading pair eigenmode.
5. Compute the amplitude threshold $T_{\rm amp}$.
6. Compute phase stiffness and $T_{\rm phase}$.
7. Compute holonomy, defect density, and $T_{\rm hol}$.
8. Rank by $T_c^{\rm OPH}=\min(T_{\rm amp},T_{\rm phase},T_{\rm hol})$.
9. Prefer candidates where processing raises the bottleneck gate while
   preserving the other gates.
10. Synthesize or process top candidates under professional materials controls.
11. Verify with structure, zero resistance, Meissner, stiffness, and gap-symmetry
    receipts.
12. Feed failures back into $K_{\rm mat}$, $\rho_s$, and
    $\mathcal D_{\rm hol}$.

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
- holonomy/defect receipt: grain-boundary, vortex, disorder, and loop-obstruction
  data;
- process-retention receipt: lifetime and relaxation margin for metastable
  phases.

## Falsifiers

The OPH high $T_c$ branch fails on a material family if:

- the computed leading negative pair eigenmode has the wrong symmetry under
  controlled source assumptions;
- a candidate has strong pair amplitude and stiffness but persistent holonomy
  obstruction, and the protocol classifies it as superconducting anyway;
- a claimed material $T_c$ lacks Meissner and zero-resistance receipts under
  declared controls;
- the design functional ranks candidates by pairing strength while ignoring
  phase stiffness or holonomy bottlenecks;
- a processing claim depends on hidden or irreproducible history outside the
  public material receipt.

## Scope

The OPH explanation of high $T_c$ closes at theorem-package level. The recipe
surface remains a research protocol: compute the pair kernel, stiffness, and
holonomy obstruction for declared material candidates; synthesize top-ranked
candidates; verify the receipts; update the source model from failures.

## References

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
