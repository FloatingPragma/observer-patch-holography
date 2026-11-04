# Local currents and the interacting radial multiplet

## Source and current domain

Let the source be the complete fermion Fock space on the connected open cube
`{-R,...,R}^3`, with integer `R>=1`, two CAR modes per site, positive `v,a,m`, and
the original Hamiltonian

\[
 H=\sum_{x,y}c_x^\dagger h_{xy}c_y
   +U\sum_x(n_{x\uparrow}-1/2)(n_{x\downarrow}-1/2),\qquad
 h_{x,x+e_j}=\frac{iv}{2a}\sigma_j,
\]
\[
 h_{xx}=B_x=m\sum_j(-1)^{x_j}N_j(x)\sigma_j,\quad B_0=0.
\]

The radial profile is `N(x)=x/|x|`. The same current theorem applies to the
supplied rational-map textures of `winding_sector_search`, of any nonzero degree:
their values at `e1,e3` are nonparallel after the edge transport below. No
texture selection principle is used. The interaction is centered exactly as
written. A chemical potential or a scalar onsite mass changes the symmetry.

A local current means a **sum of arbitrary even Hermitian onsite operators**
`Q=sum_x q_x`. Each `q_x` acts on the full four-state local Fock algebra, commutes
with local fermion parity, and has no restriction to a chosen spectral sector.
This includes pair creation and annihilation, and the quartic parity operator.
Constants are retained and then quotiented as physically trivial derivations.
The classification does not cover multisite currents, emergent nonlocal
algebras, or broken gauge currents coupled to dynamical gauge links.

## Exhaustive local classification

The even Hermitian algebra at one site has real dimension eight. A basis is

\[
 I,\quad n-1,\quad S_j=c^\dagger\sigma_jc\ (j=1,2,3),\quad
 P^++P^-,\quad-i(P^+-P^-),\quad Z=(-1)^n,
\]

where `P+=c_up† c_down†`. The last basis element is quartic. The centered
contact equals `UZ/4` and commutes with every even onsite operator.

In `[Q,H]`, an onsite term and a hopping-edge term have distinct tensor support.
Every edge term has zero partial trace over either endpoint: it is odd at both
endpoints. Hence terms on distinct edges cannot cancel, and `[Q,H]=0` is
equivalent to each edge commutator and each onsite mass commutator vanishing.

The complete sixteen-state CAR calculation on an axis-j edge has sixteen real
unknown coefficients and a rank-eight linear equation. Its full kernel is:
two independent constants; equal number coefficients; spin matrices related
by `A_y=sigma_j A_x sigma_j`; opposite pair coefficients `f_y=-f_x`; and **both
quartic coefficients zero**. `algebra.edge_classification` checks all 256
matrix equations and the dimension of their entire kernel for each axis.
This is a necessary and sufficient classification, with no discarded columns.

Conjugations by the three Pauli matrices commute. Therefore edge transport
on the connected cube determines every spin matrix from `A_0`, with

\[
 A_x=aI+\sum_j b_j(-1)^{\sum_kx_k-x_j}\sigma_j,
 \qquad f_x=(-1)^{\sum_kx_k}f_0.
\]

Transporting `B_x` back to the origin gives
`m(-1)^sum(x) N(x).sigma`. Thus the remaining condition is
`[b.sigma,N(x).sigma]=0` at every nonzero site. The radial values at `e1,e3`
force `b=0`. For the rational map `w -> w^k`, the equator point `w=1` maps to
itself and the north pole maps to itself; either target-orientation convention
preserves their nonparallelity. The same rank-three constraint follows for
every nonzero signed degree. Parallel local directions without this transport
would not suffice.

Consequently every conserved local current, modulo constants, is a real
linear combination of

\[
 \eta_z=(N_{\rm total}-V)/2,\quad
 \eta_+=\sum_x(-1)^{\sum_jx_j}c_{x\uparrow}^\dagger c_{x\downarrow}^\dagger,
 \quad\eta_-=\eta_+^\dagger.
\]

Here `eta_x=(eta_++eta_-)/2`, `eta_y=(eta_+-eta_-)/(2i)` are Hermitian.
The full CAR relations give `[eta_z,eta_+]=eta_+` and
`[eta_+,eta_-]=2eta_z`. The Lie algebra is precisely `su(2)_eta`; its
number-preserving subalgebra is `u(1)`. The coefficients do not depend on
volume or interaction strength. This is an open-cube theorem; a nonbipartite
periodic identification need not admit the alternating pairing sign.

For direct source verification, put `Delta=diag((-1)^sum(x)) tensor i sigma_2`.
Every hopping edge obeys `h_xy Delta_y+Delta_x h_yx^T=0`, since its two signs
are opposite. Every traceless `B_x` obeys `B_x J+J B_x^T=0`. Thus
`h Delta+Delta h^T=0` on all source rows and `[H,eta_+]=0`. The complete
54-row mass64 matrix is checked independently. The contact commutator also
vanishes: empty and doubly occupied sites have equal centered-contact energy.

## Exact interacting spectral intertwiner

Consume the complete-Fock isolation theorem in
[`radial_interacting_core/PROOF.md`](../radial_interacting_core/PROOF.md).
For the retained source `R=a=v=1`, `m=64`, `U=1/65536`, it gives exactly four
low states: a two-dimensional ground eigenspace at particle number `V=27`,
and one state each at `V-1,V+1`, with equal energy `E_G+g`, `g>0`.
All other states have a strictly higher certified threshold. This is an
all-number theorem; no many-body complement is projected out of the source.

Since eta commutes with `H`, it preserves its eigenspaces. Applying `eta_+`
or `eta_-` to the ground space would give a same-energy state at `V+2` or
`V-2`; neither exists. Both operators annihilate that whole space, and
`eta_z=0` there. The Kramers ground pair therefore consists of **two eta
singlets**, regardless of the choice of orthonormal ground basis.

Let `C_-` be a normalized charged state at `V-1`. Its eta-z eigenvalue is
`-1/2`. Since there is no same-energy state at `V-3`, `eta_- C_-=0`. Hence

\[
 \|\eta_+C_-\|^2
 =\langle C_-|\eta_-\eta_+|C_-\rangle
 =-2\langle C_-|\eta_z|C_-\rangle=1.
\]

Define `C_+=eta_+ C_-`. This fixes a phase and is the unique normalized
charged state at `V+1`; `eta_- C_+=C_-` and `eta_+ C_+=0`. The isometry

\[
 J:(e_-,e_1,e_2,e_+)\longmapsto(C_-,G_1,G_2,C_+)
\]

intertwines the actual Hamiltonian and all eta generators exactly, including
their full-space actions: `HJ=J diag(E_G+g,E_G,E_G,E_G+g)` and
`eta_j J=J t_j`. Its eta Casimir is `diag(3/4,0,0,3/4)`. Leakage is zero.
This constructs an actual source representation. Eta is not identified with
electroweak isospin, spatial spin, Lorentz spin or a family label.

This same current acts on the complete source field algebra. At every site,
the operator pairs `(c_up^dagger,D_x c_down)` and
`(c_down^dagger,-D_x c_up)`, with `D_x=(-1)^sum(x)`, are eta doublets:
the eta-z commutators have weights `+1/2,-1/2`, and the raising/lowering
commutators have coefficient one. These are exact full-CAR identities.
They also determine the actual interacting pole amplitudes. In the phase
convention `C_+=eta_+ C_-`, define

\[
 (A_+(x))_{sa}=\langle C_+|c_{xs}^\dagger|G_a\rangle,\qquad
 (A_-(x))_{sa}=\langle C_-|c_{xs}|G_a\rangle.
\]

Using `eta_- G_a=0` gives `A_+(x)=D_x i sigma_2 A_-(x)` at every site,
for both ground columns. No free-state replacement enters this relation.
It expresses the source eta charge of the addition/removal poles, without
assigning Standard Model quantum numbers to them.

## A strictly nonzero local-spin leakage

There is a quantitative obstruction to identifying an arbitrary local spin
probe with that exact doublet. Define the center Pauli current
`S_j=c_0† sigma_j c_0` and `D0=(n_0-1)^2`. On the entire local Fock algebra,
`S_j^2=I-D0` and `S_j` preserves total particle number.

The exact antiunitary `Theta=R D (i sigma_2) K` of the radial source preserves
each unique charged eigenspace and sends the center `S_j` to `-S_j`. Therefore
`<C_±|S_j|C_±>=0`. No other member of the complete low quartet has the same
particle number. If `P_low` denotes that full quartet projector, then
`P_low S_j C_±=0` for all three axes.

Write `alpha` for the exact zero-frame weight at the center, and `t` for the
full interacting spectral-projector angle from the core theorem. The graph
frame is `(I,-D_Q^{-1}V)G^{-1/2}`. Since `V†V=(3/2)I` and
`||D_Q||<=m+2=66`, it gives

\[
 \alpha\le\frac{1}{1+(3/2)/66^2}=\frac{2904}{2905},\qquad
 t=\frac{27}{16252874}.
\]

The filled-sea Wick identity gives charged reference parity
`(1+alpha^2)/2`. The actual charged effect differs by at most `t`.
Consequently

\[
 \|(I-P_{\rm low})S_jC_\pm\|^2
 =1-\langle C_\pm|D0|C_\pm\rangle
 \ge\frac{1-(2904/2905)^2}{2}-t>\frac1{4096}.
\]

Thus the leakage norm is strictly greater than `1/64`. It concerns the actual
interacting states and entire Fock complement. A separate one-body identity
illustrates why compression alone can conceal this problem: for the normalized
zero frame `W`, `W†S_j^(one)W=alpha sigma_j`, while
`||(I-WW†)S_j^(one)W||^2=alpha(1-alpha)`. The compressed brackets differ from
compressed source brackets by `2i alpha(1-alpha) epsilon_jkl sigma_l`.
The interacting bound does not replace the true spectral states by these
one-body columns.

## Comparison with the declared native matter current

`native_matter_attachment` constructs a declared five-mode spin lift and its
nonvacuum even fifteen-dimensional exterior module from the actual twelve
native ports. This module has nontrivial `su(3)` and `su(2)` factors and a
commuting hypercharge action. Its auxiliary register is not included in the
radial two-mode source by this calculation. The complete native color and weak
current ranks and the charge/Casimir multiplets are independently retained.

Every Lie homomorphism from `su(3)` into the classified local `su(2)_eta` is
zero: its kernel is an ideal, simplicity would make a nonzero map injective,
and an injection of the eight-dimensional algebra into dimension three is
impossible. A nonzero homomorphism from the weak `su(2)` onto eta is an
isomorphism; its centralizer inside eta is zero. It leaves no nonconstant
commuting onsite hypercharge current. Constants generate no CAR derivation.
Hence no faithful restriction of the complete native matter-current action
exists inside this source's **conserved onsite-even current algebra**.

This conclusion uses the exhaustive CAR classification, source texture and
actual current brackets, rather than the number of low modes alone. It does
not exclude gauge interactions with additional dynamical fields, a different
source action, broken currents, nonlocal symmetries or an emergent continuum
current. Native repair preservation, source-state preparation, physical
particle identification and an interacting continuum need separate theorems.

## Executable checks

The receipt retains the exact local eight-dimensional basis, all three full
edge kernels, the literal source pairing equation, actual low-current matrices,
strict rational leakage bounds and the complete target current multiplicities.
Independent three-site controls retain all 64 Fock states and every number
sector, verify the interacting eta ladders and leakage, and reject unstaggered
pairing and an uncentered contact. Floating eigensystems are diagnostics; the
theorem admission uses exact algebra and inherited rational isolation bounds.
