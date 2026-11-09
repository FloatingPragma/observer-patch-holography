# Quantum Higgs feedback on the finite gauge carrier

The finite quantum scalar of
[`higgs_join`](../interacting_extension/higgs_join/PROOF.md) replaces a fixed
neutral Higgs number in a complete memory-controlled interaction. The
construction retains the same fifteen internal currents, physical CAR,
quantum links, native record law and declared coherent pulse clock. It
proves a finite interacting gauge-Higgs action limit. It does not select a
Higgs vacuum, remove the naive kinetic mirrors, or prove a relativistic
interacting limit.

## Quantum scalar and the exact vertex norm

The specified scalar carrier is
\(\mathcal H_H=\mathbb C|0\rangle\oplus\overline{\mathbb C^2}\), with
\(B_h^\dagger=|h\rangle\langle0|\), \(h=0,1\). This is the finite source
model; it is not asserted to be an invariant truncation of unrestricted
bosonic creation. Write \(\Gamma^h_{i,e}\) for its actual charged-Higgs
tensor, \(L_h\) for the left doublet and \(e^c\) for the charged singlet.
For a family matrix \(Y\), define the source vertex

\[
 W(Y)=\sum_{h,i,f,g,s,t}\overline{\Gamma^h_{i,e}}Y_{fg}
       \epsilon_{st}B_h^\dagger c_{L_i,f,s}^\dagger c_{e,g,t}^\dagger,
 \qquad V(Y)=W(Y)+W(Y)^\dagger.                         \tag{1}
\]

This is exactly the directed `W_fl_fe` convention of
[`interacting_carrier.py`](../interacting_extension/higgs_join/interacting_carrier.py),
without adding the transposed internal tensor a second time. The source
fixtures authenticate their own bytes; the code replays all72 Higgs Ward
substitutions and the equality of their twelve-dimensional current algebra
with the literal exterior lift. No execution of the large retained
438/1596-dimensional interacting packet is claimed here.

For a unit rank-one family projector \(P=u u^\dagger\), independent unitary
family changes on the left and singlet fields reduce (1) to one family.
The other fermions are spectators. With its six active CAR modes, choose
the literal source signs so that

\[
 A_h=c_{L_h,\uparrow}^\dagger c_{e,\downarrow}^\dagger
       -c_{L_h,\downarrow}^\dagger c_{e,\uparrow}^\dagger,
 \qquad W=\sum_h B_h^\dagger A_h.
\]

Since \(W^2=0\), \(\|V\|^2=\|\sum_hA_h^\dagger A_h\|\). Exact occupation
arithmetic on all64 fermion states gives eigenvalues
\(0,1,2,3,4\), with multiplicities \(29,16,10,8,1\). The vacuum attains4.
Consequently

\[
                      \|V(P)\|=2.                    \tag{2}
\]

The complete192-dimensional scalar/CAR matrix, including every occupied
scalar and fermion input, is checked. For positive \(Y\), spectral
decomposition and real linearity imply \(\|V(Y)\|\le2\operatorname{Tr}Y\).
Equation (2) needs no selected matter or scalar state.

## Local Gauss law and a common interacting Hamiltonian

For each anti-Hermitian source generator \(G\), the local Gauss operator is
the sum of \(d\Gamma(G)\), the conjugate-doublet scalar generator and the
incident link generators of
[`GAUGE_ATTACHMENT.md`](../standard_model_source_map/GAUGE_ATTACHMENT.md).
The tensor Ward identity and the full CAR creation commutator give
\([\mathcal G_x(G),V_x(Y)]=0\) for every family matrix \(Y\). The executable
checks all twelve commutators on the full192-state local carrier. It also
checks \([N_F-2N_H,V]=0\). This operator is CAR-even; neither its scalar
factor nor the auxiliary memory supplies fermion parity.

Take \(H_0\) to be the existing finite quantum-link electric terms, invariant
loop terms, spatial Weyl hopping, scalar/fermion number energies, and any
other explicitly proved invariant source terms. Setting the Yukawa coupling
to zero in `higgs_join.build_join` isolates precisely this common part in
its retained two-vertex model. Combining its scalar vertex with the
Weyl-hopping attachment configures the corresponding finite spatial model.
All coefficients and the controller service are supplied action data.
The proof does not replace quantum links by classical backgrounds.

The scalar is charged, so it must appear in the Gauss generator. Dropping
that generator while retaining (1) fails the weak/hypercharge Ward identity.
Every complete record-conditioned propagator below preserves the full
local Gauss kernel, including a reference carrying any required external
charge. An initial state outside that kernel is not made physical by the
unitary.

## Complete native feedback with quantum Higgs and gauge dynamics

Use the same twelve-memory source and four independent flags from
[`FEEDBACK_LIMIT.md`](../standard_model_source_map/FEEDBACK_LIMIT.md).
For duty \(D=1,C,E,CE\), apply

\[
 K_{m,r}=\sqrt{w_r}|F(m,r)\rangle\langle m|\otimes
 e^{-ih[H_0+g\chi_D(r)V(P_{F(m,r)})]},\qquad h=t/N.     \tag{3}
\]

All original events, disabled pulses and final memories are retained. The
auxiliary current legs which generate the event probabilities are separate
from the arbitrary quantum field and reference. Reset-record pulses depend
on the reset record/destination, never the old memory. Their stationary
mean is exactly \(V(gZ_D)\), with the original traces \(1,p,q,pq\); no duty
probability is divided out. The executable replays all48 coarse rows and
the source mean for each duty. Fine native labels are restored by the same
conditional completely positive refinement as in the owning theorem.

For independent memories at \(V\) finite cells, put
\(p_* =1-(1-p)^V\) and \(L=2V|g|\). The existing interaction-picture proof
applies to the entire gauge/scalar/fermion/reference carrier and yields

\[
 \delta_N\le\frac{4p_*Lt+4(1+p_*)L^2t^2}{N(1-p_*)}.   \tag{4}
\]

This is unhalved reduced trace distance, uniformly over references and
initial CQ memory correlations. The complete flagged output is within
\(\sqrt{2\delta_N}\) of the instrument which keeps the same initial
correlations and entire path law while applying
\(e^{-it[H_0+\sum_xV_x(gZ_D)]}\) on every path. The ideal archive is not a
unitary or a stationary product archive. No \(\|H_0\|\) enters (4), even
when \([H_0,V]\ne0\). Per-step complete implementation errors add by
telescoping. Finite records do not have their mean unitary exactly.

Thus full quantum Higgs and gauge dynamics can coexist with this controlled
finite mean-action limit. Equation (4) is not uniform in source volume.
The auxiliary gate overhead, coherent pulse service, attempt/archive
capacity and common physical implementation clock remain explicit costs.
No infinite tape is executed or vacuum state selected by this certificate.

## What a quantum link can and cannot replace

Each nonzero source link component is \(U=|\ell\rangle\langle0|\) on its
60-dimensional register. Its nonzero two-state block has \(U^2=0\), norm1,
numerical radius1/2, and

\[
 (I-U)^\dagger(I-U)\text{ has eigenvalues }
 (3\pm\sqrt5)/2.
\]

Every normalized link/reference state therefore obeys
\(\|(U-I)\psi\|\ge(\sqrt5-1)/2\), and every density operator obeys
\(|1-\operatorname{Tr}(\rho U)|\ge1/2\). Spectator levels have singular
value1 and do not weaken the bound. These statements rule out an
arbitrarily accurate identity-component substitution in this fixed link
carrier.

Rescaling a small mean does not make the bidirectional hopping classical.
Let \(\lambda=\operatorname{Tr}(\rho U)\),
\(q=\operatorname{Tr}\rho(P_0+P_\ell)\). Positivity gives
\(2|\lambda|\le q\le1\), and

\[
 \operatorname{Tr}\rho\{(U-\lambda)^\dagger(U-\lambda)
 +(U^\dagger-\bar\lambda)^\dagger(U^\dagger-\bar\lambda)\}
 =q-2|\lambda|^2\ge2|\lambda|^2.                       \tag{5}
\]

For nonzero mean, the combined relative RMS error, divided by
\(\sqrt2|\lambda|\), is at least1. The mean condition in (5) is essential;
the displayed estimate is not asserted for arbitrary scalar targets.
This excludes a background shortcut, not collective dressed transport.

For distinct oriented edges, the contracted operator
\(W_{xz}=U_{x x_1}U_{x_1 x_2}\cdots U_{x_n z}\) has cancelling interior
Gauss variations. The full bilinear \(c_x^\dagger W_{xz}c_z\) cancels both
endpoint variations too. A two-edge weak-block control checks all36
vertex/current identities; linear independence of CAR bilinears promotes
its coefficient identity to every occupation sector. The path without
endpoint fields transforms at its endpoints. Neither object is assumed
unitary, and a gauge-invariant bilocal does not supply an isolated local
charged-particle state or its dispersion.

## Chiral regulator compatibility

The unchanged classical-background Weyl symbol
\(-c\sum_j\sigma_j\sin k_j\) has eight corners, four of each orientation.
An invertible onsite family basis change cannot remove its zeros. Gauge
Ward identities and adding the quantum scalar vertex do not prove that an
interacting mirror sector decouples; the free eight-pole spectrum is not
asserted for the interacting Hamiltonian.

The existing [overlap construction](../common_source/overlap_regulator.py)
instead uses unitary classical links, a Euclidean Wilson operator, a
uniform sign-kernel gap and a modified chiral projector. The
[Hernandez--Jansen--Lüscher locality theorem](https://arxiv.org/abs/hep-lat/9808010)
requires smooth admissible gauge fields. The link controls above do not
meet those assumptions. [Lüscher's exact lattice symmetry](https://arxiv.org/abs/hep-lat/9802011)
changes the chiral transformation; it is not an identity for the unchanged
naive regulator.

An operator-valued Wilson kernel can admit a gap and sign functional
calculus without approaching a classical identity link. There is an exact
control on a one-cell periodic **four-direction Euclidean** Wilson stencil,
using the source \(q=6Y=6\) singlet link components. Each of its four full
60-level registers is either its active vacuum/singlet pair or one of58
retained dark levels. On a sector with \(k\le4\) active links,
\(U_\mu=(X_\mu-iY_\mu)/2\) and the literal Wilson assembly gives

\[
 A=D_W-I=B-iC,\quad B=3I-\tfrac12\sum_{\mu\le k}X_\mu,
 \quad C=\tfrac12\sum_{\mu\le k}\gamma_\mu Y_\mu,\quad H=\gamma_5 A,
\]
\[
 H^2=B^2+\tfrac{k}{4}I-\tfrac12\sum_{\mu\le k}\gamma_\mu Z_\mu
 \succeq[(3-k/2)^2+k/4-\sqrt{k}/2]I\succeq I.           \tag{6}
\]

The last spin/link sum squares to \(kI\), and the commuting \(X_\mu\)
give \(B\succeq(3-k/2)I\). All five sectors, up to64 rows, are replayed
exactly. Thus the sign kernel exists, including the dark link levels.
However, \(\operatorname{Re}A\succeq I\) and \(\|A\|\le6\).
For its unitary polar factor \(T=A|A|^{-1}=\gamma_5\operatorname{sign}H\),
an eigenvector of \(T\) with eigenvalue \(e^{i\theta}\) satisfies
\(\cos\theta\langle|A|\rangle\ge1\), hence \(\cos\theta\ge1/6\).
The overlap operator therefore obeys

\[
 D=I+T,\qquad D^\dagger D\succeq(7/3)I.               \tag{7}
\]

It has no zero, whereas the classical identity-link stencil on one site
has \(D_W=0\) and \(D=0\). This is a finite quantum-link sign construction
and a failure of that identity-link comparison, not an all-volume gap for
all M15 blocks or a Lorentzian chiral continuum theorem.

A gap alone also does not finish the full-CAR join. If \(P_{ij}\) act on a shared
quantum-link carrier and \(b_i=\sum_jP_{ij}c_j\), then

\[
 \{b_i,b_k^\dagger\}=P_{ik}
       +\sum_{j,l}[P_{lk},P_{ij}]c_l^\dagger c_j.       \tag{8}
\]

Unlike scalar background entries, the entries need not commute. The code
uses the actual charge-singlet source edge, whose one-fermion link matrix
is \(A=\begin{psmallmatrix}0&iU\\-iU^\dagger&0\end{psmallmatrix}\).
Let \(\Pi=A^2\), \(K=I-\Pi+A\), and \(P=(I-K)/2\). This retains the dark
complement with sign+1 and has \(K^2=I\). It commutes with both endpoint
Gauss generators. With \(\gamma_5=\operatorname{diag}(I,-I)\),
\(D=I+\gamma_5K\) satisfies the GW relation exactly. Nevertheless (8) has
a nonzero defect, including eigenvalues \(\pm1/4\). This is a gapped
dressed-edge control using the actual source algebra, not the cubic Wilson
kernel or a no-go theorem for other many-body chiral constructions.

The existing [chiral electric completion](../interacting_extension/chiral_electric/PROOF.md)
supplies a connection and a positive geometric term on a fixed-rank
admissible abelian coordinate patch. It does not supply the nonabelian
quantum-link many-body join needed here. A recent
[overlap-measure preprint, section10](https://arxiv.org/html/2610.09093v1#S10),
uses the same \(10\oplus\bar5\) representation. It treats admissible
Euclidean gauge fields and leaves reflection positivity and an interacting
continuum definition open. No independent replay of that measure theorem
is claimed by this finite certificate.

The concrete next regulator test is a uniform gap/locality estimate for an
explicit quantum-link Wilson kernel, followed by a compatible many-body
fermion construction and the same complete charged readouts. A pointwise
background projector or deletion of mirror records does not supply those
steps. No physical species, mass, chiral continuum or vacuum selection is
assigned by the finite interacting construction above.
