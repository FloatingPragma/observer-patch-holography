# A finite gauge, fermion and Higgs system

## Domain

The graph has two vertices and two oppositely oriented edges. Each edge carries

\[
\mathcal L=\mathbb C|0\rangle\oplus\bigoplus_r(V_r^*\otimes V_r),
\qquad \dim\mathcal L=60.
\]

The blocks are the native nonvacuum even exterior representations:

| Block | Internal indices | Charge 6Y | Dimension |
|---|---|---:|---:|
| u conjugate | 0,1,4 | −4 | 3 |
| L | 10,11 | −3 | 2 |
| Q | 2,3,5,6,7,8 | 1 | 6 |
| d conjugate | 12,13,14 | 2 | 3 |
| e conjugate | 9 | 6 | 1 |

Only L and e-conjugate fermions are occupied. There are F families and two spin
components, with canonical anticommutation relations in the ordering
`((15*site+internal)*F+family)*2+spin`. Colored fermion occupation is zero;
up/down pair vertices are disabled. The scalar at vertex 0 has vacuum or one
conjugate-doublet quantum. This scalar cap defines the finite model. It is not
an invariant cutoff of an unrestricted bosonic creation operator.

The conserved integer is \(Q_*=N_F-2N_H\). The Q*=0 sector has no static
reference; Q*=1 has an e-conjugate dual reference at vertex 0. Summed Gauss
constraints imply

\[
N_L=N_H,\quad N_e=Q_*+N_H,\quad
q_{6Y}^{\rm bulk}=6N_e-3N_L-3N_H=6Q_*.
\]

Thus the occupied scalar accompanies one L and one or two e-conjugate
fermions. An unconjugated scalar cannot satisfy these integer charge equations
in either sector. The reference has charge −6Q*, making the combined state
neutral.

## Complete physical spaces

`native_inputs.py` recomputes the exterior lift of the source currents and
all 72 Higgs Ward identities. A rational generator basis spans exactly the
same 12-dimensional algebra. The Gauss operator acts by the fermion bilinear,
the scalar generator and the reference charge, with −T transpose on the source
of a link and T on its target.

The exact constructor partitions by occupation and irreducible link blocks.
It removes nonzero Cartan weights, then includes all images in all 24 Gauss
equations. Exact kernel vectors and an equal modular lower rank over the prime
1,000,003 exhaust each kernel. Rational orthogonalization and the entire lifted
Gram matrix are checked. Every operator is applied to every physical basis
vector and its full raw complement residual must vanish.

A separate representation count explains the dimensions. With no occupied
scalar, the Q*=0 sector has six closed flux states. For Q*=1, an e at vertex 0
has six flux states and an e at vertex 1 has two, giving 16F states. With a
scalar, the internal multiplicities are:

| Q* | e at 0 | e at 1 | L at 0 | L at 1 |
|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 2 | 1 |
| 0 | 1 | 0 | 8 | 1 |
| 1 | 0 | 2 | 0 | 1 |
| 1 | 1 | 1 | 2 | 1 |
| 1 | 2 | 0 | 8 | 1 |

Equal color representations are necessary at both vertices. A weak singlet
link block contributes one contraction when L and the scalar share a vertex;
a weak doublet contributes two because 2 tensor 2-dual contains a singlet and
triplet. When L and the scalar occupy different vertices, the two weak link
representations must differ. The central charges select the entries above.
In particular, the mixed u-conjugate/d-conjugate flux contributes to the
multiplicity two. Colorless mobile matter does not justify deleting it.

Multiplication by the fermionic occupation choices gives

\[
d_0=6+48F^2,\qquad d_1=16F+64F^3-20F^2.
\]

They are 54/60 for F=1 and 438/1,596 for F=3. The three-family producer checks
126,108 operator columns with zero complement leakage. The separate raw audit
checks all 48,816 lifted Gauss equations, 36,612 Higgs spin identities, 15,768
addition anticommutators and 12,204 reference-order signs.

## Hamiltonian and currents

The native raising operator is \(U^r_{ab}=|\bar a,b;r\rangle\langle0|\).
For the two oriented edges, form

\[
D_e=\sum_{r,a,b,f,g,s}(F_e)_{fg}
c^\dagger_{x_e,a,f,s}U^{r,e}_{ab}c_{y_e,b,g,s},
\quad B=\sum_{r:C_3(r)>0,a,b}U^{r,0}_{ab}U^{r,1}_{ba}.
\]

There is one endpoint Casimir per link:

\[
E=\sum_e(2C_3+\tfrac12C_2+q^2/12)_e.
\]

The literal charged-lepton tensor has components Gamma(0,10,9)=1 and
Gamma(1,11,9)=1. Using epsilon(0,1)=1 and epsilon(1,0)=−1, its finite scalar
raising operator gives

\[
W=\sum_{h,i,f,g,s,t}\overline{\Gamma^h_{i,9}}Y_{fg}\epsilon_{st}
b_h^\dagger c^\dagger_{0,i,f,s}c^\dagger_{0,9,g,t}.
\]

The e creator acts first. The explicit adjoint annihilates L first and e
second. Restricting the first tensor index to L prevents double-counting its
symmetric internal coefficient. Each adjoint is checked independently.

\[
H=\kappa_E E-\kappa_B(B+B^\dagger)+\epsilon_FN_F+\epsilon_HN_H+E_{\rm ref}Q_*
-t\sum_e(e^{-i\theta_e}D_e+e^{i\theta_e}D_e^\dagger)
+g(e^{-i\phi}W+e^{i\phi}W^\dagger).
\]

For currents defined by phase derivatives,

\[
i[H,N_1]=J_1-J_0,\quad i[H,N_H]=J_H,\quad i[H,N_F]=2J_H.
\]

Literal fermion spin generators commute with H. These are finite spin and
continuity identities, not a Lorentz or continuum theorem.

## Interacting response

In Q*=1 split the scalar vacuum P of dimension 48 from its full occupied
complement Q of dimension 1,548. With H blocks A,C-adjoint,C,D,

\[
P(z-H)^{-1}P=[z-A-\Sigma(z)]^{-1},\quad
\Sigma(z)=C^\dagger(z-D)^{-1}C.
\]

The off-diagonal block C is exactly the Higgs term. For Im z>0,

\[
-\operatorname{Im}\Sigma(z)=(\operatorname{Im}z)
C^\dagger(\bar z-D)^{-1}(z-D)^{-1}C\succeq0.
\]

In the declared primary control the norm at z=1+0.5i is approximately
1.793879153064. C has full column rank and the hopping/Higgs commutator has
Frobenius norm about 12.5786412. These numerical checks establish active
interaction in this configured finite system.

The reference has odd parity. With its mode ordered after all matter modes,
\(B_j^\dagger=c_j^\dagger d_{\rm ref}^\dagger\) and
\(B_j=d_{\rm ref}c_j\) are even. Their reference CAR sign is +1 on both
physical sibling maps. The six addition maps retain Pauli blocking and obey
\(B_jB_j^\dagger=I-n_j\) on Q*=0.

For a prepared stationary mixture of selected initial eigenvectors indexed by
m with weights p_m, the full response is

\[
C_{jk}(t)=\sum_{n,m}p_m e^{-i(E_{1n}-E_{0m})t}
\langle m|B_j|n\rangle\langle n|B_k^\dagger|m\rangle.
\]

Every initial/final pair is retained, including zero matrix residues.
Numerically grouped low-energy levels retain their distinct initial energies;
a near degeneracy is not replaced by an exact degeneracy. Summing residues
recovers the complete equal-time probe Gram matrix.

The primary control uses the exact prepared family Y, with
\(F_0=Y\), \(F_1=Y+(1-\cos\log(5/4))D+\sin\log(5/4)J\), and Higgs family Y.
All dynamical strengths and the scalar energy are one; the fermion and
reference energies are zero. Its Q*=0 ground is numerically unique. The first
three distinct addition frequencies are approximately

\[
2.58004293348\,10^{-8},\quad0.001006653510723,\quad0.452014700447579.
\]

Exact-input 50/80-digit residual refinement supports these values and the
positive tiny difference. It does not supply interval-certified eigenvalue
ordering or a rigorous enclosure. The matrix-rounding diagnostic includes the
union of exact and floating supports, including entries that vanish exactly.

## Physical boundary

The graph, three families, two spin components, scalar cap, placement and
strength of couplings, state preparation and measuring reference are supplied
controls. Preparing Y does not select the common interacting ground state.
No energy unit, measured particle mass, nonabelian chiral continuum or Einstein
dynamics is derived by this finite system.

The identity
\(\epsilon_FN_F+\epsilon_HN_H=\epsilon_FQ_*+(2\epsilon_F+\epsilon_H)N_H\)
shows why a bulk energy shift at fixed reference changes addition frequencies.
A consistent shift of the compensating reference cancels it exactly. This is
relative reference calibration, not a proof that intrinsic masses are arbitrary.
