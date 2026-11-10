# Complete finite Peter-Weyl links

## Source and retained representations

The declared group is the connected compact group
\(G=S(U(3)\times U(2))\). Its twelve Hermitian current matrices are the actual
native matrices returned by `physical_attachment.chiral.representation`.
They act on the native exterior representation, with the five invariant
blocks in the following order:

| Primitive charge \(6Y\) | Native indices | Dimension |
|---|---|---|
| \(-4\) | \(0,1,4\) | 3 |
| \(1\) | \(2,3,5,6,7,8\) | 6 |
| \(6\) | \(9\) | 1 |
| \(-3\) | \(10,11\) | 2 |
| \(2\) | \(12,13,14\) | 3 |

The retained representation set \(\mathcal R\) consists of the trivial
representation, these five blocks, and their duals. The dual current is
\(-T^T\). These eleven irreducible representations have distinct hypercharges.
Each scalar commutant is verified by exact linear algebra on all twelve
matrices. Hypercharge labels alone are not the representation test.

In normalized Haar measure, let
\(|r,a,b\rangle=\sqrt{d_r}D^r_{ab}\in L^2(G)\). The link carrier is

\[
\mathcal L=P L^2(G)
=\bigoplus_{r\in\mathcal R}V_r\otimes V_r^*,\qquad
\dim\mathcal L=1+2(3^2+6^2+1^2+2^2+3^2)=119.
\]

The index conventions use the identical finite basis and endpoint currents
as `oriented_links`, including all 119 states. Dual link sectors describe flux
orientation; they add no matter modes. The operator law is literal compressed
group multiplication and is specified independently of the minimal link law.

## Complete Clebsch derivation

For every native field block \(r\), source \(s\in\mathcal R\), and retained
target \(t\), `clebsch.py` solves

\[
(T_r^A\otimes I+I\otimes T_s^A)C_{rs}^{t,\alpha}
=C_{rs}^{t,\alpha}T_t^A,\qquad A=1,\ldots,12.
\]

The hypercharge equation forces \(q_t=q_r+q_s\), leaving at most one retained
target. Exact diagonal Cartan equations eliminate zero coefficients; the
remaining equations include every matrix entry of every current. The full
nullspace is normalized using
\(C_\alpha^\dagger C_\beta=\delta_{\alpha\beta}I_{d_t}\).
All 55 field/source pairs are exhausted, giving 37 nonzero channels of
multiplicity one. Ten channels meet the vacuum and 27 have both source and
target excited. The receipt gives every normalized embedding and every
excluded pair. The five native/native exterior-product channels agree, up to
an irrelevant embedding phase, with the independently constructed wedge
projectors using `particle_gauge_band.quantum_link_completion.wedge_sign`.

The full multiplication formula is

\[
U^r_{ab}=P M_{D^r_{ab}}P
=\sum_{s,t,\alpha}\sqrt{\frac{d_s}{d_t}}
\sum_{m,n,i,j}
(C_{rs}^{t,\alpha})_{am,i}
\overline{(C_{rs}^{t,\alpha})_{bn,j}}
|t,i,j\rangle\langle s,m,n|.
\]

The factor follows by multiplying
\(\sqrt{d_s}D^s_{mn}\) by \(D^r_{ab}\), decomposing the tensor product with
the isometric embeddings, and expressing the result in normalized target
matrix coefficients. The representation-basis construction is described in
[Zohar and Burrello, New Journal of Physics 18 (2016) 043008](https://doi.org/10.1088/1367-2630/18/4/043008).
The finite matrices here are derived from the stated native currents.

The vacuum creation and dual absorption terms have coefficient
\(1/\sqrt{d_r}\). They agree with the corresponding minimal-link terms.
Their agreement does not remove the other 27 channels. In particular,
\(e^c\otimes L\simeq\overline L\), and the actual coordinate obeys

\[
U^{e^c}_{00}|L,0,0\rangle=|\overline L,1,1\rangle.
\]

This unit-amplitude excited-to-excited transition is absent from the minimal
link law. A rescaling by \(\sqrt{d_r}\), used for a different coordinate
normalization in `quantum_link_completion`, is excluded here.

## Gauge currents, orientation and bounded paths

The endpoint currents on a forward block are \((-T^T,T)\), and on its dual
block are \((T,-T^T)\). Exact full-matrix checks prove

\[
[L_0^A,U_{ab}]=-\sum_c T^A_{ac}U_{cb},\qquad
[L_1^A,U_{ab}]=\sum_c U_{ac}T^A_{cb}.
\]

There are 1416 endpoint identities, covering every native coordinate and all
119 input states. Inversion is the unitary involution
\(J|r,a,b\rangle=|\bar r,b,a\rangle\), with vacuum fixed. Its exact identities
are \(J^2=I\), \(JL_0J=L_1\), and \(JU_{ab}J=U_{ba}^\dagger\).

On \(V_r\otimes L^2(G)\), multiplication by the internal matrix \(D_r(g)\)
is unitary. Its compression
\(V_r=\sum_{ab}E_{ab}\otimes U^r_{ab}\) therefore has norm at most one.
The finite matrices additionally satisfy
\((V_r^\dagger V_r)^2=V_r^\dagger V_r\) exactly. The respective projection
ranks are \(45,72,24,80,90\), on dimensions \(357,714,119,238,357\).
The blank-link internal subspace has unit norm, so every transporter has norm
exactly one. Products of these complete operators are contractions. Every
distinct-edge blank path also attains norm one. The path executor checks all
238 single-link columns in the forward and adjoint branches.

The excitation parity \(K=\operatorname{diag}(1,-1,\ldots,-1)\) does not
anticommute with these coordinates. Its excited-to-excited transitions remove
the minimal-link local star symmetry
\(C_x=(-1)^{N_x}\prod_{e\ni x}K_e\). The complete Hamiltonian has an explicit
nonzero \([H,C_x]\) basis witness; the physical propagation calculation below
establishes this distinction inside the gauge-singlet sector as well.

## Complete periodic quantum Hamiltonian

`model.column` acts on all
\(2^{30N^3}119^{3N^3}5^{N^3}\) basis states. It keeps every native CAR mode,
five scalar states per site, and every link state. It uses all axial kinetic
paths of `periodic_family.kinetic_terms`, the three actual quantum Yukawa
tensors of `periodic_family.local_yukawa`, and scalar/link excitation
penalties. The coordinates above replace the path link operators. No
composite hopping term is added and no scalar or link expectation is inserted.
The action is explicitly Hermitian, with every kinetic and Yukawa word
accompanied by its adjoint, and each word has even CAR degree.

The Ward identities and actual scalar/Yukawa identities prove local Gauss
invariance on the full domain. The geometry and spin lifts in
`oriented_links.rotations` apply because the carrier, endpoint currents and
inversion are identical. All 24 proper cubic rotations preserve the complete
path words for \(N=2,4\). The antiperiodic seam phases and spin lifts are kept;
matter internal charges are not conjugated. Actual occupied-column tests
include 33 Hamiltonian/rotation comparisons with unequal rational couplings,
three translation comparisons and three Gauss commutators. The full-word
argument applies to every basis state; these executions test its implementation.

The coefficients and regulators are source design choices. This change of
link law transfers no ground state, spectrum, chirality limit or continuum
theorem from another link architecture. No implementation by the native idle
link source is asserted.

## Physical neutral transfer and complete records

Take \(N=2\), spacing \(a=\pi\), and all three Yukawa coefficients and both
penalties equal to one. The state \(\Omega\) has empty CAR occupations and
blank scalar/link memories. It is gauge invariant and nonstationary under
this Hamiltonian. Define

\[
D_{x,s}=\frac{\phi_{x,0}c_{x,s,11}-\phi_{x,1}c_{x,s,10}}{\sqrt2},\qquad
|x,s\rangle=D_{x,s}^\dagger\Omega.
\]

Each packet is a normalized two-component exact gauge singlet. The composite
has operator norm one and is noncanonical. The physical even-parity
preparation retains a two-mode, gauge-neutral fermionic reference of spin
one-half and raw lepton label minus one:

\[
B^\dagger=\frac{D_{0,\uparrow}^\dagger r_\downarrow^\dagger
-D_{0,\downarrow}^\dagger r_\uparrow^\dagger}{\sqrt2},\qquad
K=B^\dagger+B.
\]

On the empty source/reference input, \(K\Omega=\Phi\) and \(K\Phi=\Omega\),
so \(e^{-i\pi K/2}\) prepares the normalized spin singlet. The complete
320-state active preparation verifies all twelve currents, three spin
generators and total parity, and is extended with the full spectator CAR
parity convention. The reference is retained with probability one. Evolution
acts on the source with the same \(H\) specified above. The two source spin
branches have identical arrival probabilities by a proper rotation.

For \(|s\rangle=|0,\uparrow\rangle\) and
\(|t_\sigma\rangle=|(1,0,0),\sigma\rangle\), the full first and second images
have respectively 433 and 62,799 basis terms, with

\[
\|Hs\|^2=\frac{1651}{8},\qquad
\|H^2s\|^2=\|H^2t_\sigma\|^2=\frac{15494801}{128}.
\]

The code verifies all 96 cell/current Gauss constraints on both \(s\) and
\(Hs\). The first two powers retain quark, lepton, scalar and nonblank link
states. Exact full-source moments give

\[
\langle t_\sigma|H^k|s\rangle=0\quad(0\leq k\leq4),\qquad
\langle t_\downarrow|H^5|s\rangle=-\frac{151i}{32},\qquad
\langle t_\uparrow|H^5|s\rangle=0.
\]

The fifth moment is evaluated as \(\langle H^2t|H|H^2s\rangle\), retaining
both complete second images. For each actual Hamiltonian summand, matching
its unchanged spectator coordinates excludes only identically zero matrix
elements. Every onsite and diagonal middle term has zero intersection. Of
the 24 stored directed edges, only the two joining the origin and destination
can contribute. Every native field block on those edges is included; the
receipt gives their separate contributions. The method is also checked
against an unfiltered full-Hamiltonian bilinear on a physical input.

The exact \(\pi\) rotation about \(z\), including the antiperiodic seam, gives
eigenvalue \(-i\) to \(s\), \(+i\) to \(t_\uparrow\), and \(-i\) to
\(t_\downarrow\). Since it commutes with \(H\), the same-spin transfer
amplitude is zero at every time. Thus

\[
A(t)=\langle t_\downarrow|e^{-itH}|s\rangle
=-\frac{151}{3840}t^5+R_6(t).
\]

Each scalar matrix, CAR monomial and coordinate entry has norm at most one.
The complete kinetic/onsite word triangle bound gives \(\|H\|\leq1896\).
Applying three powers on each side of the sixth derivative yields

\[
|A^{(6)}(t)|\leq\|H^3t_\downarrow\|\,\|H^3s\|
\leq1896^2\frac{15494801}{128}=:M,
\qquad |R_6(t)|\leq\frac{Mt^6}{720}.
\]

Let \(m=151/32\), \(t_*=3m/M\), and \(L=mt_*^5/120\). The remainder is
at most \(L/2\). The exact rational values are printed in the receipt.
The complete two-outcome readout is

\[
P=\sum_\sigma|t_\sigma\rangle\langle t_\sigma|,\qquad I-P.
\]

These are gauge-invariant even operators. This exclusive arrival detector
also checks that every other source memory is blank; its support is the
complete finite source. Both records are retained. Its probabilities obey

\[
\Pr(P)\in[L^2/4,9L^2/4],\qquad
\Pr(I-P)\in[1-9L^2/4,1-L^2/4].
\]

The local detector of `neutral_mobility` is the rank-two onsite packet
projector times \((I-C_{(1,0,0)})/2\), with all other fields unrestricted.
Its support is the destination and incident links. Both exclusive target
states lie in its image, as the executor verifies. Thus its complete
arrival/no-arrival probability intervals are respectively \([L^2/4,1]\) and
\([0,1-L^2/4]\). The narrower upper bound above applies to the exclusive
detector. No source component is postselected away during evolution.
This finite nonstationary transfer result supplies no
stationary particle spectrum, mass, relativistic dispersion, macroscopic
speed or continuum limit.

## Reproduction and custody

`--write` and `--check` recompute every CG matrix, current identity, complete
source power and probability bound. The receipt contains rational/algebraic
strings, with no floating-point evidence. Declared inputs and the transitive
import closure are hashed; missing inputs fail explicitly. Deliberately
incorrect orientation, normalization and missing-input controls are tested.
No numerical receipt from another architecture is consumed.
