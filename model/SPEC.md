# CCF-LT-2 executable specification

Official model: **A-ccf2-local-capacity-v2**. Sections 1–8 match the supplied
CCF-2 contract, with local implementation links and state-only wording.
[Design](DESIGN.md), [runtime](ccf_lt2/README.md), [proofs](ccf_lt2/bridges/PROOF.md)
and [authority checks](AUTHORITY.json) accompany this specification.
The current native repair law is a calibrated intermediate test of this model; the search for the correct native repair law is ongoing.

## 1. Finite carrier and conventions

Choose a finite oriented complex with vertex positions, positive edge lengths, unit edge directions, and closed oriented faces. Its vertex-edge incidence $B$ has tail +1 and head −1; face boundaries $C$ satisfy $BC^T=0$.

At every vertex $v$, retain all 96 native fermion modes, with CAR ordering

$$
i(v,s,c,f)=96v+48s+3c+f,
\quad s\in\{0,1\},\ c\in\{0,\ldots,15\},\ f\in\{0,1,2\}.
$$

Internal coordinates 0–14 are lexicographic 2-subsets followed by 4-subsets of $\{0,1,2,3,4\}$; coordinate 15 is sterile. Their integer charges $n_c=3Q_c$ are

$$
(-2,-2,2,-1,-2,2,-1,2,-1,3,0,-3,1,1,1,0).
$$

Each edge carries $E_e\in\{-L,\ldots,L\}$, default $L=6$. The shift $U_e|q\rangle=|q+1\rangle$ vanishes at the upper wall; its adjoint vanishes at the lower wall. Negative powers mean adjoint powers, never an inverse or cyclic wrap. Thus $[E,U]=U$, while products retain hard-wall projectors.

The local native charge is $n_v=\sum_{s,c,f}n_c c^\dagger_{vscf}c_{vscf}$. With a declared reference charge $r_v$,

$$
G_v=(BE)_v-n_v-r_v.
$$

Let $P_G$ project onto all $G_v=0$. Keep the complete carrier and incompatible preparation outcomes: at fixed cutoff some native occupations have no compatible flux. A root reference has declared charge range and preparation domain.

Each vertex also has the three-state $\Phi$ and 16-state capacity/memory registers below. Intrinsic dimension is

$$
2^{96|V|}(2L+1)^{|E|}48^{|V|}.
$$

External references, purifications, apparatus, and tapes add declared dimensions. Matrix-free keys are
`(capacity_indices, (phi_index, (native_bitmask, flux_tuple)))`.
CAR relabeling uses graded permutations; disjoint physical algebras are even.

## 2. The inherited native source

The inherited two-cell source has 192 native modes and **114 source factors**, whose full normal-ordered export has 22,273 monomials. Its canonical export SHA-256 is

```
591949feb5120a050c49a8b4edf52db9a920bfaf0db920d8c7937089b4613013
```

[snapshot.py](ccf_lt2/snapshot.py) and [gauging_compiler.py](ccf_lt2/gauging_compiler.py) reconstruct this exact coefficient object from the [114-factor snapshot](ccf_lt2/source_snapshot.json) and verify the digest, without an external checkout. The primary generator is RER `model/openwave_common_mass_source/source.py::full_source_terms` at commit `b37ba4dae0e069dc030ceb0f130e52905fbb5465`. **Executed coefficients are the binary64 snapshot promoted to exact rationals.** The formulas and calibration inputs below describe its construction; recomputing rounded entries at another precision changes the source.

### 2.1 Internal operators and the 114 factors

Lift fundamental 5-by-5 generators by the exterior derivative representation on $\Lambda^2\mathbb C^5\oplus\Lambda^4\mathbb C^5$, and extend by zero on the sterile coordinate. Color generators $C_A$ are $(E_{ij}+E_{ji})/2$ and $-i(E_{ij}-E_{ji})/2$ for the three pairs $0\le i<j\le2$, followed by

$$
C_7=\tfrac12\operatorname{diag}(1,-1,0,0,0),\qquad
C_8=\tfrac12\operatorname{diag}(1,1,-2,0,0).
$$

Weights are $\kappa_{1\ldots7}=1,\kappa_8=1/3$; $W_a$ are lifted Pauli/2 generators on fundamental indices 3,4.

$P_{\rm col}$ selects all coordinates except 9,10,11,15; $P_Q$ selects (2,5,7,3,6,8); $P_L$ selects (10,11). The baryon matrix $b$ is +1/3 on $P_Q$, −1/3 on other colorful coordinates, zero otherwise. Write $J_x(M)=c_x^\dagger Mc_x$, with identity on unmentioned factors, $\sigma^\mu=(I,\boldsymbol\sigma)$, and $\eta_\mu=(1,-1,-1,-1)$.

The complete factor list is:

| Count | Factor |
|---:|---|
|1|Colorful quadratic source: local up/down masses and $t\sum_i(n_{0i}+n_{1i}-c^\dagger_{0i}c_{1i}-c^\dagger_{1i}c_{0i})$, over colorful modes|
|16|$U\kappa_A J_x(C_A)^2$, two cells and eight color generators|
|6|$J_s J_x(P_{\rm col}\otimes\sigma_a/2)^2$, two cells and three spin generators|
|6|$J_w J_x(P_QW_aP_Q)^2$, two cells and three weak generators|
|2|$A J_x(b)^2$|
|8|$\kappa_{\rm colour}\kappa_A[J_0(C_A)+J_1(C_A)]^2$|
|3|$\eta[J_0(P_QW_aP_Q)+J_1(P_QW_aP_Q)]^2$|
|24|$2G\eta_\mu J_x(\sigma^\mu\otimes P_QW_aP_Q)J_x(\sigma^\mu\otimes P_LW_aP_L)$|
|48|$E_jL_j^\dagger L_j$, the full charged and neutral lepton quadratic source|

The colorful up/down masses are $\bar m-\delta/2$ and $\bar m+\delta/2$, respectively. Down means a conjugate down singlet or a $Q$ coordinate with $T_3<0$. All products are actual CAR products, with all normal-order scalar terms retained.

### 2.2 Parameters and calibration ancestry

All energies in the executable source use the electron anchor $0.00051099895069$ GeV as unit. The strong-source parameter ancestry is

| Parameter | Supplied value |
|---|---:|
|$\bar m$|9343959228105 / 68719476736|
|$\delta$|129333251000 / 51099895069|
|$A$|49097929115755 / 34359738368|
|$U$|16345196922453 / 68719476736|
|$t$|1/8|
|$J_s,J_w$|1/4, 1/4|
|$\eta,\kappa_{\rm colour}$|1, 1|
|$G$|1/64|

Calibration inputs, in GeV: proton 0.93827208943, neutron 0.93956542194, charged pion 0.13957039, plus $|g_A/g_V|=1.2753\pm0.0013$ and the electron anchor. These determine strong-source parameters; remaining coefficients are supplied choices.

For charged leptons take the raw family tuple $q=(a,1,b)$, with

$$
a=\frac{15075139752763506821122969774584631}{3117091410011963900862644383689281103},\quad
b=\frac{52423389005790828147871788202091689820}{3117091410011963900862644383689281103},
$$

and energy ratios $\lambda_f=q_f/q_0$. These are source mass inputs. For neutral leptons take normal ordering, $m_0=0.001$ eV, and masses

$$
(m_0,\sqrt{m_0^2+7.49\,10^{-5}},\sqrt{m_0^2+2.513\,10^{-3}})\ \mathrm{eV}.
$$

Convert them to source units by multiplying by $10^{-9}/0.00051099895069$. The neutral mixing matrix is the standard $R_{23}R_{13}(\delta)R_{12}$ convention with

$$
\sin^2\theta_{12}=0.308,\quad \sin^2\theta_{23}=0.470,\quad
\sin^2\theta_{13}=0.02215,\quad\delta=212^\circ.
$$

Oscillation values are calibrated; normal ordering and lightest mass are declared. There is no Majorana term. The rounded $R_\nu$ is

$$
\begin{pmatrix}
.8226008752730573&.5487966836634492&-.12621394713219705+.07886722734641406i\\
-.33205011909985804+.0449778352415709i&.6536281816274402+.03000688129676683i&.6779303061524835\\
.4569094560606586+.04776255547084231i&-.5193039841954019+.03186470235275267i&.7199031184819246
\end{pmatrix}.
$$

The 48 lepton factors use $L_j=\sum_i\overline{A_{ij}}c_i+B_{ij}c_i^\dagger$. At cell $x$, spin $s$, family $f$, put $\epsilon_s=1-2s$. For each pair $(l,r)=(11,9)$ or $(10,15)$, let $a=(x,s,l,f)$, $b=(x,1-s,r,f)$. In both pairs $A_{bb}=1/\sqrt2$, $B_{ba}=\epsilon_s/\sqrt2$. For the charged pair also $A_{aa}=1/\sqrt2$, $B_{ab}=-\epsilon_s/\sqrt2$. For the neutral pair replace these latter entries by

$$
A_{(x,s,10,\alpha),(x,s,10,f)}=R_{\nu,\alpha f}/\sqrt2,\qquad
B_{(x,s,10,\alpha),(x,1-s,15,f)}=-\epsilon_sR_{\nu,\alpha f}/\sqrt2.
$$

Unspecified entries vanish; both pair columns receive the corresponding family energy.

### 2.3 Replication and gauging

For parent $H_{02}$ with scalar $C_0$, vacuum-restrict $H_{02}-C_0$ to each cell and add $C_0/2$ to each local block $L_0,L_1$. They agree under translation. Insert $L_v$ once per vertex and the reindexed residual bond $H_{02}-L_0-L_1$ once per edge, recovering the original one-edge source exactly.

Dress each charge-homogeneous monomial $O$ by $U^\ell$ with $B\ell=\Delta n(O)$. Bonds use their actual edge; local neutral words need no shift. Retain every normal-order scalar and native mode.

## 3. Added matter terms

### 3.1 Signed Dirac kinetic term

Use the following full signed kinetic term. For charged $(l,r)=(11,9)$ and neutral $(10,15)$ pairs, define

$$
a_{l,s}=(c_{l,s}+\epsilon_sc^\dagger_{r,1-s})/\sqrt2,\qquad
a_{r,1-s}=(c_{r,1-s}-\epsilon_sc^\dagger_{l,s})/\sqrt2.
$$

Let $\Psi=(a_{l,0},a_{l,1},a^\dagger_{r,1},-a^\dagger_{r,0})^T$, $\beta=\operatorname{diag}(1,1,-1,-1)$, and $\alpha_i=\sigma_1\otimes\sigma_i$. Its common integer charge is $q_l=-3$ for the charged pair and zero for the neutral pair. Keep the parent's rounded mass matrix.

For every edge $u\to v$, length $\ell$, direction $d$, pair, and family, add

$$
\frac{ic}{2\ell}\Psi_v^\dagger(\alpha\cdot d)U_e^{-q_l}\Psi_u+{\rm h.c.}
+\frac{r}{2\ell}\left[\Psi_u^\dagger\beta\Psi_u+\Psi_v^\dagger\beta\Psi_v-
\Psi_v^\dagger\beta U_e^{-q_l}\Psi_u-{\rm h.c.}\right].
$$

Defaults: $c=r=1$, with $c>0,r\ge0$. Expand into the full native algebra, retaining both bands and pair terms.

### 3.2 Leptonic and nonleptonic currents; family splitting

Put $\bar\sigma^\mu=(I,-\sigma_1,-\sigma_2,-\sigma_3)$, $\tau_a=\sigma_a/2$. The native left-lepton doublet is $L_{s,w,f}=c_{s,(10,11)_w,f}$. The quark left fields include their conjugate partners:

$$
Q_{s,w,j,f}=\frac{c_{s,q_{w,j},f}+\epsilon_s\chi_jc^\dagger_{1-s,\bar q_{w,j},f}}{\sqrt2},
$$

where $q_\uparrow=(2,5,7)$, $q_\downarrow=(3,6,8)$, $\bar q_\uparrow=(4,1,0)$, $\bar q_\downarrow=(14,13,12)$, and $\chi=(1,-1,1)$. Rotate the down quark family by

$$
R_q=\begin{pmatrix}4/5&3/5&0\\-3/5&4/5&0\\0&0&1\end{pmatrix},
$$

leaving up unchanged. For $X=L,Q$, form all twelve Hermitian currents $J_X^{\mu a}=\sum X^\dagger(\bar\sigma^\mu\otimes\tau_a)X$, summing the full family and, for quarks, color channels. Add

$$
V_{LL}=2g_L\sum_{v,\mu,a}\eta_\mu(J_{L,v}^{\mu a})^2,\qquad
V_{QQ}=2g_Q\sum_{v,\mu,a}\eta_\mu(J_{Q,v}^{\mu a})^2,
\quad g_L=g_Q=1/4096.
$$

Add $D_F$: native occupation shifts (0,2,8) on up/anti-up and (0,1,4) on down/anti-down coordinates, summed over spins, colors, vertices. These shifts and $R_q$ are supplied, not empirical quark/CKM predictions. Keep the original semileptonic $G=1/64$ term exactly once. All additions conserve electric charge, baryon number, and Dirac lepton number.

## 4. Compact field energy and complete source shards

The complete matter-and-field source is

$$
H_S=H_{\rm parent,spatial}(U)+K_D(U)+V_{LL}+D_F+V_{QQ}
+\frac{g_E}{2}\sum_eE_e^2+g_B\sum_f\left[I-\frac{W_f+W_f^\dagger}{2}\right],
$$

Defaults: $g_E=g_B=1$, both positive. $W_f$ is the actual oriented edge-shift product.

Partition $H_S=\sum_vh_v$: split each matter coefficient equally over its native-word vertices and endpoints of nonzero dressing shifts; split a pure scalar over all vertices, electric energy over edge endpoints, and magnetic energy over face vertices. Each $h_v$ is Hermitian and Gauss invariant. This energy allocation is part of the law.

## 5. Local geometry and capacity law

### 5.1 Local response register

At each vertex let $\Phi_v$ be diagonal with values (−1,0,1), and $N_v=I+\Phi_v/4$. Let $\mathcal L_3$ be the three-point path graph Laplacian: off-diagonal −1 between successive labels and diagonal equal to degree. With edge weights $w_{uv}=1$, hopping $j=1/2$, stiffness $\kappa=1/2$, regulator $\rho=1$, and boundary weights $b_v=0$, define

$$
H_\Phi=j\sum_v\mathcal L_{3,v}
+\frac\kappa2\left[\sum_{uv}w_{uv}(\Phi_u-\Phi_v)^2+\sum_v(\rho+b_v)\Phi_v^2\right].
$$

For configurable members require nonnegative $j,\rho,b_v,w_{uv}$, positive stiffness, ordered finite response labels, and strictly positive $N_v$. Split this energy as

$$
g_v=j\mathcal L_{3,v}+\frac\kappa2(\rho+b_v)\Phi_v^2+
\frac\kappa4\sum_{u:uv}w_{uv}(\Phi_u-\Phi_v)^2,
\quad\sum_vg_v=H_\Phi.
$$

For specified local references $e_v$, default zero, put

$$
k_v=g_v+N_v(h_v-e_vI)+e_vI.
$$

### 5.2 Per-vertex volume, capacity, and memory

Each vertex has independent labels $v_{\rm vol}\in\{1,3/2\}$, $k\in\{1,2\}$, and a two-bit hardware memory $m\in\{00,01,10,11\}$. The order is volume, capacity, memory. Let $x(v_{\rm vol})=4(v_{\rm vol}-1)-1$ and $M_v=I+x/4$; its eigenvalues are 3/4 and 5/4. On this 16-state register define

$$
A_v=\frac12\mathcal L_{2,\rm vol}+\frac18\mathcal L_{2,\rm cap}
+\operatorname{diag}\left[\frac1{32}\frac{v_{\rm vol}}k+
\frac1{16}(v_{\rm vol}-1)^2+\frac1{16}(k-1)^2\right],
$$

Unmentioned factors are identities. Memory is stationary under $A_v$. Active capacity $k$ exposes $2^k$ labels while retaining all four hardware states; readout records volume, capacity, and memory modulo $2^k$.

The **current whole Hamiltonian** is

$$
\boxed{H_{\rm v2}=\sum_v\left[A_v+M_v(k_v-FI)+FI\right].}
$$

$F=0$ is a common scalar reference at all vertices; independent $F_v$ are not implemented. Every auxiliary is local.

Configurable capacity members use positive ordered volume labels, $k=1,\ldots,k_{\max}$, $2^{k_{\max}}$ hardware states, volume coordinate $x\in[-1,1]$, $M=I+\gamma x$, $0\le\gamma<1$, and nonnegative hopping/penalties. Defaults above are supplied constants.

### 5.3 Energy references and currents

Jointly shifting $h_v\mapsto h_v+cI$, $e_v\mapsto e_v+c$, and $F\mapsto F+c$ changes $H_{\rm v2}$ only by $|V|cI$. Nonuniform shifts cannot generally be absorbed by common $F$; holding references fixed makes scalar energy affect the auxiliaries.

The charge current of the current law is the derivative of that law with respect to the matter link phase:

$$
J_e=\sum_vM_vN_v\,\partial_{\theta_e}h_{v,\rm matter}(\theta)\big|_{0}.
$$

Retain second phase derivatives when required. The full law preserves Gauss charge and its associated continuity identity; field commutators retain hard-wall effects.

## 6. States, preparation, and complete observations

Before a run, declare graph, parameters, initial joint state, reference charge, preparation domain, durations, complete instruments, outputs, and parser. Autonomous evolution is $e^{-itH_{\rm v2}}$, without intermediate repair.

A spanning-tree preparation may transport nonroot Gauss syndrome by $\Delta E=-B_{ve}G_v$, preserving native charge and chord flux. Retain displaced syndromes, work, and finite-cutoff/capacity failures; reversible uncomputation requires its declared blank-register domain.

Instruments obey $\sum_aK_a^\dagger K_a=I$. Retain $p_a=\operatorname{tr}K_a\rho K_a^\dagger$, conditional quantum states, and all classical outcomes, including no-arrival, mismatch, overflow, and complements. Coarse graining sums outcomes without changing the original denominator. Tapes report exhaustion or explicitly account for export/erasure.

On $P_G$, a complete joint boundary-flux/charge meter satisfies $1_S^TBE=\sum_{v\in S}(n_v+r_v)$; keep its zero-probability mismatch outcome. Charge/current are distinct contexts. Meter work is $\operatorname{tr}H_{\rm v2}(\sum_aK_a\rho K_a^\dagger-\rho)$.

Clocks compare actual energy differences through coherent alternatives and a physical phase reference. Stationary global phase is unobservable; a vacuum-to-pair gap is labeled $2m$, not $m$.

### 6.1 Explicit coherent phase instrument

Use the oriented triangle with edges (0,1), (0,2), (1,2), default parameters, zero flux/capacity indices, flattened $\Phi$ index one, and empty native occupation $|i\rangle$. Lexicographically order matter word/shift tuples. Select the first two-creator word with nonzero link shift, internal labels in $\{9,11\}$, and positive real vacuum-to-pair coefficient. Its action defines $|f\rangle$, keeping auxiliaries unchanged.

On the complete carrier set $|\pm\rangle=(|i\rangle\pm i|f\rangle)/\sqrt2$. Kraus operators are $P_\pm=|\pm\rangle\langle\pm|$ and $I-P_+-P_-$. Evolve $|i\rangle\langle i|$ with full $H_{\rm v2}$; retain all three outcomes and the complete initial Hamiltonian column. The receipt tests the exact initial probability derivative, with terminal-dephasing and Dirac-deleted controls using this same instrument. This supplied phase convention is not a local hardware or self-starting mass-clock construction.

### 6.2 Declared mathematical reference family

For the finite reference certificates choose $L=6$, $r=0$, $\beta=1/1024$, and

$$
\rho_\beta=P_Ge^{-\beta H_{\rm v2}}/Z.
$$

This is a supplied mathematical reference. For $d=\operatorname{rank}P_G$ and $\|H_{\rm v2}\|\le M$, $\rho_\beta\ge e^{-2\beta M}P_G/d$. Compute $d$ by summing all flux assignments with local charge multiplicities from $\prod_{i=1}^{96}(1+z^{n_i})$, times $48^{|V|}$.

For a declared region of ordered tensor factors let $\rho_R=\operatorname{tr}_{\bar R}\rho_\beta$ and $P_R=\operatorname{supp}\operatorname{tr}_{\bar R}P_G$. Keep all supported boundary-charge sectors. On this support the logarithm $K_R=-\log\rho_R$ exists, and marginalization composes. The finite entropy identity is

$$
S(\tau)-S(\rho_R)=\operatorname{tr}[(\tau-\rho_R)K_R]-D(\tau\|\rho_R).
$$

Physical modular-kernel and stationarity identifications remain separate obligations.

## 7. Spatial, coframe, and scaling families

### 7.1 Angular and radial refinement

The genuine finite refining shell family has level $m$, $n=2^m$, and tetrahedral directions obtained from vertices (1,1,1), (1,−1,−1), (−1,1,−1), (−1,−1,1). Use outward faces (0,1,2), (0,3,1), (0,2,3), (1,3,2). Subdivide each triangular face into $n^2$ barycentric triangles and normalize its points to the sphere. Add radii $r_j=1+j/n$, $j=0,\ldots,n$; join angular triangles and radial edge strips into prism cells. Every vertex carries the same complete 96-mode native block.

It has $(2n^2+2)(n+1)$ vertices, $4n^2$ triangles per shell, diameter bound $(8\sqrt6+1)/n$, exact boundary identities, and vertex/edge-path refinement maps. Full-source, instrument, and hard-wall cutoff compatibility under refinement are unproved.

### 7.2 Mechanical coframe and shift variation

The implemented constant-background adapter starts at the default $c=r=1$ source. Choose a real constant matrix $V$ with positive determinant and a real shift vector $b$. On an edge with unit direction $d$, put $p=V^{-T}d$. Split the exact Dirac polynomial into three directional derivative pieces $D_A$, Wilson piece $W$, and the identity-spin derivative piece $P$, obtained from the directional derivative by replacing its $\alpha$ matrix with $I_4$ while retaining the $i/(2\ell)$ coefficient and Hermitian adjoint. Replace it by

$$
\sum_Ap_AD_A+(p\cdot p)W+(b\cdot p)P.
$$

Multiply the electric coefficient on that edge by $d^TV^TVd/\det V$; multiply a face magnetic coefficient by the same expression with its unit face normal. Keep the inherited parent and added contact terms unchanged in this mechanical variation. At $V=I,b=0$ the source is exactly the original $H_S$.

For symmetric $S$ and $V=I+\epsilon S$, the directional derivative is $p'=-Sd$, the Wilson multiplier derivative is $-2d^TSd$, and the field-weight derivative is $2d^TSd-\operatorname{tr}S$. This provides six strain and three shift derivatives of the same finite source. It is not a spatially varying coframe, spin connection, or covariant Lorentz stress tensor.

### 7.3 Declared analytic scaling family

A separate analytic family, expressible through existing assembler arguments, takes spacing $a>0$, $g,c>0$, $r\ge0$, and

$$
g_E=g^2/a,\qquad g_B=c^2/(g^2a),\qquad \text{Dirac speed}=c,\quad\text{Wilson parameter}=r.
$$

On a uniform cubic chart its quadratic field dispersion is $4g_Eg_B\sum_j\sin^2(ak_j/2)$; the free Dirac symbols are $p_j=c\sin(ak_j)/a$ and $w=(r/a)\sum_j(1-\cos ak_j)$. These have bounded-momentum estimates in the proof record. The spatial refinement supplement remains unexecuted. Separate fixed-graph background and weak-coupling investigations are recorded separately in [scaling scope](ccf_lt2/bridges/SCALING.md); they do not change the default coefficients or establish a physical common speed. Irregular spherical shells require compatible Hodge/volume weights and source/reference/instrument limits.

## 8. Physical theorem interface and remaining obligations

The finite model does not identify its response with a physical metric. The retained OPH/RER Einstein route requires one common shrinking family, physical modular/stress, area, and stationarity estimates whose errors divided by the fourth power of radius tend to zero, and a positive spacetime-independent limiting area coefficient. It is detailed in [Einstein bridge](ccf_lt2/EINSTEIN_BRIDGE.md). Spectra, rates, physical continuum identification, vacuum selection, and cosmological/dark-sector attachment are tracked in the [dependency map](particle_joint_reconstruction/DEPENDENCIES.md).

## 9. Implementation map

All implementation paths below are local to this repository. The package
contains its native source snapshots and runs without the OpenWave checkout.

| Definition | Implementation |
|---|---|
| Exact native CAR operations, parent export, replication, fields, energy shards | [core.py](ccf_lt2/core.py), [source.py](ccf_lt2/source.py), [snapshot.py](ccf_lt2/snapshot.py) |
| Signed kinetic term and complete current additions | [dirac.py](ccf_lt2/dirac.py), [interactions.py](ccf_lt2/interactions.py) |
| Finite chart data | [spatial.py](ccf_lt2/spatial.py) |
| Local response and auxiliary factors | [local_geometry.py](ccf_lt2/local_geometry.py), [geometry.py](ccf_lt2/geometry.py) |
| Whole law and native column entry point | [local_capacity.py](ccf_lt2/local_capacity.py), [source.py](ccf_lt2/source.py) |
| Refining shell and constant-background derivatives | [gravity bridge](ccf_lt2/bridges/gravity/coframe.py) |
| Faithful Gauss/Gibbs and regional reference certificates | [reference bridge](ccf_lt2/bridges/reference/check.py) |
| Analytic scaling supplement | [kinematics bridge](ccf_lt2/bridges/SCALING.md) |
| Complete phase instrument | [phase instrument](ccf_lt2/bridges/kinematics/phase.py) |
| Charge-preserving preparation | [FIELD_PREPARATION.md](ccf_lt2/FIELD_PREPARATION.md) |

Run `python -m model` for the complete source calculation and strict result
comparison. Run `python tools/check_ccf_lt2_authority.py --authority-dir
../openwave/candidates/current` to compare the local contract with the supplied
specification and design. Omitting `--authority-dir` checks the local pinned
contract without requiring another checkout.
