# Induced regional modular interaction on the finite scalar source

The declared interacting source has computable regional modular operators.
Restricting its faithful state to a fermion family produces a pair term
inside the two-particle sector. This term is visible to a gauge-invariant
local probe and changes entropy without changing any local number-sector
weight. The calculation uses the source state itself, rather than assigning
a desired geometric modular kernel.

## One source, preparation and clock

The [owning source](../particle_observer_outputs/HIGGS_RECORDS.md) has
`Y=diag(1,2,3)`, pair cost `d=1`, and complete 307-dimensional conserved
sector `N_F-2N_H=0`. Its basis is

\[
|\Omega\rangle=|h=0,\mathrm{mask}=0\rangle,\qquad
|h,\mathrm{mask}\rangle,\quad h=1,2,\quad
\operatorname{popcount}(\mathrm{mask})=2
\tag{1}
\]

on all eighteen fermion modes. The normalized bright vector is
`b=W(Y)Omega/g`, with `g^2=56`; the Hamiltonian is
`H=V+N_H`, with `V=g(|Omega><b|+|b><Omega|)`. Its energies are
`(-7,8,1)` with ranks `(1,1,305)`.

The [faithful preparation](../particle_observer_outputs/THERMAL_CLOCK.md)
chooses `beta=log(2)` and derives the full density
`rho=exp(-beta H)/Z`, `Z=71809/256`. Its per-state spectral weights are
`(32768,1,128)/71809`. This is a declared reproducible preparation. Its
canonical purification retains a conjugate 307-dimensional reference.
The global field density is gauge invariant and includes nonsinglet states;
the joint field/reference purification is a singlet. The compensating
reference is part of that construction and is not discarded from the
program when a local marginal is computed.

Put `D0=1077135`. In the vacuum/bright/dark decomposition,

\[
\begin{aligned}
\rho={}&\alpha P_\Omega+wP_{\mathrm{occ}}+(\beta_b-w)P_b
       +c(|\Omega\rangle\langle b|+|b\rangle\langle\Omega|),\\
\alpha={}&262151/D0,\qquad \beta_b=229384/D0,\qquad
w=128/71809,\qquad c=-32767g/D0.
\end{aligned}\tag{2}
\]

The symbol `beta_b` denotes a population, not inverse temperature.
The same clock is `U(t)=exp(-itH)`. The exact global modular relation
`-log rho=log(2)H+log Z I` is inherited; no separate Gibbs proof is needed.

## The actual regional map and its algebra

Select the six modes of family `f`, with coefficient `y_f=f+1` for
`f=0,1,2`. Move these modes to the front of the fermionic tensor ordering.
For every basis vector in (1), the sign is the parity of the occupied-mode
inversions in that permutation. The map

\[
J|h,\mathrm{mask}\rangle
 =(-1)^{\mathrm{inversions}}|\mathrm{local\ mask}\rangle
                    \otimes|h,\mathrm{outside\ mask}\rangle
\tag{3}
\]

is an isometry on all 307 input columns. The reduced density is
`R_f(rho)=Tr_out J rho J^dagger`. The executable uses these actual
labels and signs. The same expression tensor identity on any reference
defines a completely positive, trace-preserving map for joint local and
reference outputs. No outcome probability is renormalized.

The local six-mode Fock space has dimension 64. Occupations zero, one and
two form a 22-dimensional positive support; the other 42 basis outcomes
have probability zero in every input from this retained source sector.
Those zero outcomes are retained in the complete ambient detector.
The natural regional algebra preserving this source sector is

\[
\mathcal A_f=\mathbb C\oplus M_6\oplus M_{15},\qquad
\mathcal H_{307}\simeq
(\mathbb C\otimes\mathbb C^{133})\oplus
(\mathbb C^6\otimes\mathbb C^{24})\oplus
(\mathbb C^{15}\otimes\mathbb C^2).
\tag{4}
\]

The multiplicities count all exterior scalar and fermion configurations:
`133=1+2 binom(12,2)`, `24=2*12`, and `2=2*binom(12,0)`.
The sum is `133+6*24+15*2=307`.

The dual injection `j_f(A)=J^dagger(A tensor I)J` is a unital star
homomorphism on (4). It is not a homomorphism on the full 64-dimensional
CAR algebra: a local creation operator takes this fixed sector to a
different total `N_F-2N_H`, so its compression vanishes while the
compression of its positive product does not. A larger source carrier
is required for that different local algebra. This domain distinction
prevents a partial trace from being misidentified with a multiplicative
restriction.

There is also a trace-preserving conditional expectation onto `j_f(A_f)`.
Take the regional operator partial trace, divide its number blocks by
`(133,24,2)`, and inject back. The resulting map is completely positive,
unital, trace preserving, idempotent and fixes the regional algebra.
It realizes the single-region expectation contract of
[`CPRestrictionNet.lean`](../../Lean/QFT/CPRestrictionNet.lean), including
`expect_idem` and `expect_retracts`. No full generating spatial net,
coverage or spacetime interpretation is inferred. Preservation of the
trace does not mean preservation of the interacting state (2).

## Exact regional spectrum and induced pair term

Let `P_n` be the local number-n projector. The two scalar components of
`W(Y)Omega` restrict to orthogonal local pair vectors of squared norm
`2 y_f^2`. After normalization call them `b_1,b_2`, and set
`B_f=sum_h |b_h><b_h|`. This rank-two projector lies inside `P_2`;
`D_f=P_2-B_f` has rank thirteen. Direct partial trace of (2) gives

\[
\rho_f=p_0P_0+p_1P_1+p_bB_f+p_2D_f,
\tag{5}
\]

where

\[
\begin{aligned}
p_0&=\alpha+132w+(\beta_b-w)(1-4y_f^2/56),\\
p_1&=24w,\qquad p_2=2w,\qquad
p_b=2w+(\beta_b-w)2y_f^2/56.
\end{aligned}\tag{6}
\]

All four values are strictly positive. Their ranks are `(1,6,2,13)`;
`p_0+6p_1+2p_b+13p_2=1`. For the first family they are

\[
\left(\frac{5087653}{7539945},\quad
\frac{3072}{71809},\quad
\frac{83746}{7539945},\quad
\frac{256}{71809}\right).
\tag{7}
\]

The exact modular Hamiltonian on positive support is therefore

\[
K_f=-\log p_0\,P_0-\log p_1\,P_1-\log p_2\,P_2
       +\kappa_f B_f,
\qquad \kappa_f=\log(p_2/p_b)<0.
\tag{8}
\]

The last term is an induced regional pair interaction. It is calculated
from the actual coupled state; it is not an added action parameter.
For the first family `kappa_0=log(13440/41873)`. The canonical modular
flow on the faithful support is

\[
\sigma_s(A)=\sum_{a,b}(p_a/p_b)^{is}P_aAP_b,
\quad (P_a)=(P_0,P_1,B_f,D_f).
\tag{9}
\]

Restricting (9) to the accessible algebra (4) preserves its number blocks.
The formula also defines a mathematical flow on the full positive-support
matrix algebra; implementing cross-number operators requires the larger
carrier mentioned after (4).

The pair term has an admitted gauge-invariant probe. On each bright pair,
apply left-spin `Z`: it changes the spin-singlet combination to its
orthogonal spin-triplet partner `d_h` in the same weak component. Then
`A=sum_h(|b_h><d_h|+|d_h><b_h|)` preserves number and commutes with the
source's twelve gauge generators. The executable checks these commutators
against the owning one-family current matrices. It also verifies

\[
[K_f,A]=\kappa_f\sum_h
(|b_h\rangle\langle d_h|-|d_h\rangle\langle b_h|)\ne0.
\tag{10}
\]

The apparatus has a chosen spin axis; it is not declared invariant under
all spin rotations. It is a finite gauge-invariant matrix measurement
with the complete spectral apparatus supplied by the
[readout contract](../particle_observer_outputs/READOUT_CONTRACT.md).

The bare local fermion onsite energy for this fixture is
`epsilon_F N_f=0`. Even allowing arbitrary constant, `N_f` and `N_f^2`
terms cannot reproduce (10), since those commute with A. Equation (8)
provides the missing induced operator explicitly. Calling it a geometric
diamond energy would require an independent spatial and stress dictionary.

## A source-admitted entropy response

The tangent

\[
X_f=\tfrac12B_f-\tfrac1{13}D_f
\tag{11}
\]

has zero trace in every number block. It is gauge invariant, has norm
`1/2` and Hilbert-Schmidt squared norm `15/26`. Its entropy first variation
is the source-specific value

\[
DS_{\rho_f}(X_f)=\operatorname{Tr}(X_fK_f)
                 =\log(p_2/p_b)=\kappa_f.
\tag{12}
\]

This changes entropy while preserving all local number probabilities.
It measures precisely the induced pair contribution that a number-only
description loses.

The tangent has a faithful global realization on the same source and
Hamiltonian: `X_global=j_f(X_f)/2`. The factor two is the full exterior
multiplicity of the local two-particle sector. Its reduction is `X_f`,
its trace is zero and its operator norm is `1/4`. Choose
`epsilon=1/71809`, the global density's exact minimum eigenvalue. Then

\[
\rho+\epsilon X_{\rm global}
 \succeq\frac3{4\cdot71809}I>0.
\tag{13}
\]

This is an admitted different preparation on the same carrier, action and
clock. Its canonical conjugate-reference purification is a computable
finite preparation using the owning primitive or compiler contract.
The gauge-invariant tangent preserves the joint-reference Gauss route.
No passive repair preparation of this varied state is assumed.

The existing [regional entropy remainder theorem](../native_gravity_identification/PROOF.md#regional-first-law-and-its-controlled-remainder)
applies to these concrete matrices. Their support floor is
`m=p_2=256/71809`. It gives

\[
S(\rho_f+\epsilon X_f)-S(\rho_f)
=\epsilon\kappa_f-D(\rho_f+\epsilon X_f\Vert\rho_f),\qquad
0\le D\le\frac{15\epsilon^2}{52(m-\epsilon/2)}
          =\frac{15}{26\cdot71809\cdot511}.
\tag{14}
\]

This instantiates the finite theorem instead of repeating its proof.
The complete number-center weights are
`(p_0,6p_1,2p_b+13p_2)`. Their Shannon entropy and the weighted normalized
block entropies must both be retained when evaluating `S(rho_f)`. The
42 unsupported local configurations remain zero outcomes. The center
labels here are local occupation records; they are not a geometric area.

## Why the global modular clock does not restrict directly

The local number `N_f` is central in (4), so its regional modular flow is
constant. In the actual global source,

\[
\|[H,j_f(N_f)]\Omega\|^2=16y_f^2>0.
\tag{15}
\]

Thus no nonzero rescaling of regional modular time can identify it with
the restriction of this global Hamiltonian evolution on this algebra.
The global evolution does not preserve that regional algebra. This is
compatible with the global relation `-log rho=beta H+log Z I`.

There is also a direct state-preservation obstruction for regional
conditional expectations. Normalize the part of `W(Y)Omega` in family f
as `b_f`, and put `X=|Omega><b_f|+|b_f><Omega|`. It joins the central local
number-zero and number-two sectors and has exact mean

\[
\operatorname{Tr}(\rho X)=-131068y_f/1077135\ne0.
\tag{16}
\]

Every bimodule conditional expectation onto `j_f(A_f)` kills this
off-diagonal central block, because `P_0 E(X)P_2=0`. If it preserved rho,
its expectation of X would have to equal (16), a contradiction. The
tracial expectation after (4) exists, but fails to preserve rho. The
regional partial trace deliberately removes these global scalar/family
coherences from the local marginal; the global program and reference
retain them. This is not a failure of complete positivity or of the
finite regional-net theorem's trace-based expectation fields.

The scalar-only region supplies an exact comparison. Its density is

\[
\rho_H=\operatorname{diag}(262151,407492,407492)/1077135.
\tag{17}
\]

On its retained algebra `C direct_sum M2` this density is central, so
modular flow is the identity. A formal representation
`K_H=constant+log(262151/407492)N_H` has a negative number coefficient,
despite the positive global beta. Actual Hamiltonian evolution changes
`N_H`, with `||[H,N_H]Omega||^2=56`. The scalar region is therefore not
the positive-temperature global clock restricted to one factor. Cross
empty/occupied scalar operators also leave this retained sector and
require the larger ambient carrier.

## Geometric interpretation requires additional data

Equations (5)-(14) supply actual finite regional modular operators and
source-admitted entropy responses, with reference, support and center
specified. The family and scalar cuts are internal subalgebras of one
cell. They carry no spatial ball radius, proper volume, boundary area or
boosted observer foliation. They cannot directly instantiate the
[matched regional Einstein theorem](../matched_einstein_limit/PROOF.md).

That consumer needs the modular response of spatial regions in the same
gauge-coordinate field and state family, attached to its action-derived
stress with a diamond kernel. It also needs fixed-proper-volume area
response, generalized-entropy stationarity and matched shrinking-ball
errors. The induced term (8) shows why tracing an interacting source must
be performed before assuming a local modular-energy formula. Other
regions, state families, corrected stress operators and continuum limits
are not excluded by these finite controls. Stationary repair entropy is
not used as geometric entropy, and no desired Einstein contraction is
introduced as a definition.
