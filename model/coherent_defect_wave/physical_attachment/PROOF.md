# Spin, currents and chiral masses on declared field extensions

The [joined native clock](../PROOF.md), the [spatial CAR occupation law](../species/PROOF.md)
and the [Clifford routing law](../spatial/ROUTE.md) are distinct specified
programs. This construction attaches spin and currents to the actual routing
matrices, and computes the mass restrictions of the fifteen-state internal
representation. It does not identify the native wave's internal labels with
spatial spin or infer a physical mass from their multiplicity.

## Spin carried by the routing amplitudes

Use exactly the Hermitian matrices `alpha_i,beta` consumed by the routing
program. They obey the Clifford relations. For cyclic `(i,j,k)` set

\[
 J_i=-\frac{i}{2}\alpha_j\alpha_k,\qquad
 \chi=-i\alpha_1\alpha_2\alpha_3.
\]

Direct substitution gives

\[
 [J_i,J_j]=i\epsilon_{ijk}J_k,\quad
 [J_i,\alpha_j]=i\epsilon_{ijk}\alpha_k,\quad
 [J_i,\beta]=0,\quad \sum_iJ_i^2=\frac34 I.
\]

Thus `S_i(phi)=exp(-i phi J_i)` rotates the actual principal matrices,
`S alpha_j S* = sum_i R_ij alpha_i`, fixes beta and gives `S_i(2pi)=-I`.
The action is spin one half. The chirality operator is Hermitian, has two
eigenvalues of each sign, commutes with every alpha and anticommutes with
beta. The massless ordered route preserves each chirality exactly; a beta
mass term exchanges the two. The source's triplet/doublet label j is on a
different internal space and does not acquire spin through this identity.

For `H_D(p)=alpha.p+m beta`, the preceding identities imply
`S_R H_D(p) S_R*=H_D(Rp)`. The routing proof gives

\[
 \|D_a(p)^n-e^{-itH_D(p)}\|\le ta C(p),\qquad
 C(p)=\sum_{r<s}|h_rh_s|,\quad h=(p_1,p_2,p_3,m),\quad t=na.
\]

Consequently

\[
 \|S_R D_a(p)^nS_R^*-D_a(Rp)^n\|
 \le ta\{C(p)+C(Rp)\}.
\]

This supplies rotation covariance of the controlled free window with a
quantitative error. It does not assert exact rotations of an arbitrary
finite momentum grid or the ordered finite tick. In fact, at
`p=(pi/2,pi/3,0)`, `m=0`, a quarter turn about axis 3 gives squared
Hilbert-Schmidt residual **12** for the one-tick covariance equation.
The discrepancy comes from the fixed order of noncommuting directional
shifts. It vanishes in the stated small-spacing, bounded-momentum limit.

## A current measured by the same directional shifts

For an axis set `P+=(I+alpha_i)/2`, `P-=(I-alpha_i)/2`. Its actual finite
shift sends

\[
 \psi'_x=P_+\psi_{x-e_i}+P_-\psi_{x+e_i}.
\]

Orthogonality of the projectors gives an exact lattice continuity equation
with no long-wavelength approximation:

\[
 \rho'_x-\rho_x=J_{x-e_i/2}-J_{x+e_i/2},\qquad
 J_{x+e_i/2}=\psi_x^\dagger P_+\psi_x
             -\psi_{x+e_i}^\dagger P_-\psi_{x+e_i}.
\]

Every position outcome is retained. The executable checks the complete
four-site position distributions and all directed edge currents on each
axis for a coherent complex input. Concatenating the three directional
layers gives the full tick balance; its onsite mass coin preserves density.

For an internal Hermitian generator T commuting with the kinetic coin,
insert T in all three bilinears. The same proof gives its charge continuity.
Exterior second quantization turns the one-body identities into exact CAR
operator identities on every occupation sector. The current/creation
commutators and occupied-state signs are those of the existing
[fermion attachment](../../standard_model_source_map/FERMION_ATTACHMENT.md).
The declared neutral native-controlled gate commutes with local routing
number and these internal charges, so it adds no charge source at its cell.
This proves conservation through that finite coupling, without replacing
native link dynamics by a fixed background.

For the free continuum fields, the corresponding current is
`j_T^mu=psi^dagger ((I,alpha) tensor T) psi`. The spinor boost
`B_i=cosh(u)I+sinh(u)alpha_i` gives the usual Lorentz-vector transformation
of these four matrices and preserves `B_i* beta B_i=beta`. Exact controls
use `cosh(u)=5/3`, `sinh(u)=4/3`; the induced vector boost has diagonal
entries `41/9` and mixed time/axis entries `40/9`. All twelve current
identities are checked, including preservation of the Lorentz quadratic form.
This is the actual routing version of the existing two-component
`standard_model_source_map.kinetic_attachment.lorentz_currents` calculation.
A pointwise boost is not a unitary operation on the finite equal-time CAR
register; the full field and hypersurface transformation are separate.

## Which native coupling could be a Lorentz scalar?

The finite controlled-occupation interaction is a valid unitary. Its number
density limit is of the form `psi^dagger psi K`, where K is a native operator.
The displayed boost gives `B_i* B_i != I`, so number density alone is not
a Lorentz scalar. Rotational scalar behavior is insufficient to establish
Lorentz covariance of this interaction. It may instead belong to a vector
coupling if its coefficient is the time component of a covariant field.

A Dirac scalar candidate has the form `psi^dagger beta psi K`. Its spinor
factor is invariant under the displayed boosts. For example,
`K=(U-U*)/(2i)` is a source-defined Hermitian local operator commuting with
the native interior constraint. A local Hamiltonian term using this K can
therefore be specified without a unitary-logarithm choice. The complete
term is a Lorentz scalar only if the continuing field associated with K
has a compatible scalar transformation and dynamics. Identical independent
native cells do not supply that conclusion by themselves. The free routing
bound is not an interacting Lorentz-limit theorem.

## Internal representation and anomaly traces

The primitive joined source has exact pure-internal commutant `U(3) x U(2)`.
Choose its trace-free part and the declared normalization
`Y5=diag(-1/3,-1/3,-1/3,1/2,1/2)=-Q5/6`. The separate auxiliary
representation is

\[
 R=\Lambda^2\mathbb C^5\oplus\Lambda^4\mathbb C^5.
\]

The executable consumes the retained exterior basis and invariant tensors
from `common_source.chiral_candidate.action._representation`, checks that
its twelve one-body currents span this same block algebra, and uses its
exterior derivative to construct an exact rational basis. No trial source
coefficients, scalar vacuum or mass kernel from that package are consumed.

| Representation | Hypercharge | Dimension |
|---|---:|---:|
| `(3bar,1)` | `-2/3` | 3 |
| `(3,2)` | `1/6` | 6 |
| `(1,1)` | `1` | 1 |
| `(1,2)` | `-1/2` | 2 |
| `(3bar,1)` | `1/3` | 3 |

These are the weights of one all-left-handed Standard Model family without
a right-handed neutrino. The right singlets are written as left-handed
conjugates; the [PDG electroweak review, section 10.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf)
states the conventional doublet/singlet description. The representation
comparison supplies neither three generations nor measured mass values.

For the twelve Hermitian generators `T_A`, exact evaluation gives

\[
 \operatorname{Tr}T_A=0,\qquad
 \operatorname{Tr}\{T_A(T_BT_C+T_CT_B)\}=0.
\]

The cubic tensor is symmetric by trace cyclicity. All 364 independent
unordered triples are checked. The weak-doublet count, including color,
is four. These are the algebraic anomaly cancellation conditions of this
declared representation. They do not construct a regulated nonperturbative
chiral gauge measure or remove additional lattice modes.

To put R on propagating fermions, declare a one-particle Weyl factor
`C2 tensor R` and its full exterior Fock algebra. Its Cauchy dimension is
30. The even auxiliary masks labeling R do not themselves supply odd
fermionic creators. The complete exterior-current transfer from the cited
fermion attachment applies to this declared factor.

Putting the same R on all four Dirac routing components gives Cauchy
dimension 60. Its left and right chiral blocks carry the same internal
representation. This is a vectorlike doubling, rather than the fifteen
left-handed family. At zero mass a chiral projection commutes with the
free route; admitting that projection to a local gauge theory requires its
own locality, regulator and mode-counting proof.

## Exact mass selection rules and conditional spectrum

A Lorentz-scalar bilinear of left Weyl fields is
`psi_i^T epsilon M_ij psi_j`, with `M=M^T`. Hypercharge invariance requires

\[
 Y^T M+MY=0,\qquad (Y_i+Y_j)M_{ij}=0.
\]

No pair of the displayed hypercharges sums to zero. All 120 independent
symmetric entries vanish. Therefore the unbroken fifteen-field family has
**no constant gauge-invariant Weyl mass matrix**. This statement does not
forbid a broken phase, a larger matter representation, or an emergent
composite mass. It shows why a generic beta coin copied over R cannot serve
as the chiral Standard Model mass mechanism.

For a declared Higgs doublet with `Y_H=1/2`, the only charge origins satisfying
`Y_i+Y_j+/-Y_H=0` are the up, down and charged-lepton pairs. The color singlet
contraction and weak doublet invariant fix one tensor line for each pair.
The retained tensors `C_channel,k` are tested against all twelve generators:

\[
 G_A^T C_k+C_k G_A+\sum_j(H_A)_{jk}C_j=0.
\]

The up channel uses H and the down/lepton channels use its conjugate. All
72 generator/component identities hold. The internal symmetric pairing
combines with antisymmetric Weyl spin contraction to give the required
fermionic sign. This reuses the existing invariant-line construction and
does not infer coefficients from symmetry.

For a background written as `H=(0,v)` in a selected gauge, take `v>=0` to
denote that component amplitude, with no hidden square-root-two convention.
This does not assert a gauge-variant expectation in a gauge-invariant state.
The conditional one-family kernel is

\[
 K=v(y_uC_{u,1}+y_dC_{d,1}+y_eC_{e,1}).
\]

Its squared singular values are

\[
 v^2|y_u|^2\;(6\text{ copies}),\quad
 v^2|y_d|^2\;(6),\quad
 v^2|y_e|^2\;(2),\quad 0\;(1).
\]

The doubled nonzero multiplicities pair conjugate Weyl components into
three color copies each of an up and down Dirac field and one charged
Dirac field. The remaining neutral Weyl component has no mass in this
minimal renormalizable construction. These are coefficients of the
declared tree-level kernel, not interacting pole or threshold masses.

Gauge covariance leaves three complex coefficients for one family and
three arbitrary complex three-by-three matrices for three declared
families. To exhibit the freedom, `v=1`, `(y_u,y_d,y_e)=(1,1,1)` gives
`Tr K*K=14`, while `(1,2,3)` gives 48; both obey all current identities.
Thus those identities do not determine a charged mass ratio. The source
interactions, scalar state and mass matching must supply those values if
they are to be predictions.

## The extra phase symmetry and a finite compensating memory

The native internal commutant is the thirteen-dimensional `U(3) x U(2)`.
Its trace-free twelve-dimensional subgroup supplies the gauge-current
identities above, but the remaining determinant character matters for the
source dynamics. Let `N_j` count the occurrence of auxiliary index j in an
exterior basis mask. The Higgs fundamental has the corresponding unit
characters, and its conjugate has their negatives. Direct evaluation gives,
for every j and each of the six retained Higgs-channel components,

\[
 N_j^T C_k+C_kN_j+\sum_l(H_j)_{lk}C_l=C_k.
\]

All thirty identities hold. Thus each complete Yukawa vertex transforms by
`det(g_3) det(g_2)`: its five diagonal phase charges are `(1,1,1,1,1)`.
Under a common phase on the five auxiliary modes its charge is five. The
up vertex has auxiliary degrees `2+2+1=5`; the down and lepton vertices
have `2+4-1=5` or `4+2-1=5`. Restricting to the determinant-one subgroup
annuls this character, which explains the twelve-current invariance.

These auxiliary degrees are representation labels. They are not the
physical fermion number or grading of the independently declared Weyl CAR
fields. In particular, a degree-two auxiliary mask can label a physical
odd field. The scalar field and the compensating memory below are even
under that physical grading.

The native occupation dynamics and its neutral repair registers preserve
every `N_j`. If an exact source-to-field intertwiner preserves these five
charges, a constant nonzero charged Yukawa vertex cannot be generated by
that dynamics: commutation with each charge is preserved by the intertwiner,
whereas the displayed vertex has a nonzero charge commutator. This excludes
that specific symmetry-preserving assignment with neutral references.
It does not exclude other composite assignments, a charged reference, or
a different declared law.

Retaining this extra symmetry is a choice of the source-to-field map, not
a universal physical Higgs requirement. On the primitive fixed-one-load
`C5` space, a common internal phase is an overall phase and acts trivially
on its observable algebra. The mixed `Lambda2 C5 + Lambda4 C5` lift gives
that phase different charges on its two blocks. The determinant condition
therefore applies when the map retains this lifted action. A declared
extension can instead retain only `S(U(3) x U(2))` and break the extra
accidental phase symmetry. A projective or quotient description must
specify its observable algebra and show that its dynamics descends to it;
the mixed-degree relative action is not removed merely by calling it an
overall phase. Neither the trace-free gauge group nor its mass tensors
require the following memory as a unique construction.

There is a concrete finite controller that removes this algebraic obstacle
while retaining the full lifted symmetry.
Adjoin the even two-state memory representation

\[
 \mathcal M=\mathbf 1\oplus\det^{-1},\qquad
 X=|\det^{-1}\rangle\langle\mathbf 1|.
\]

The transition X has charges `(-1,-1,-1,-1,-1)`. Multiplying each complete
Yukawa tensor vertex by X and adding the adjoint gives a Hermitian interaction
invariant under the full `U(3) x U(2)`, including its determinant character.
The Weyl spin contraction and full Higgs/internal tensor are retained; a
single selected Higgs/color component alone is not that invariant vertex.

`compensator.py` executes three finite component controls to test the charge
exchange and physical parity explicitly. Two fermionic modes, a two-state
Higgs component, and the two-state memory have dimension sixteen. For the
three allowed degree pairs, write

\[
 V=a_0^\dagger a_1^\dagger\otimes |1_H\rangle\langle0_H|
       \otimes X,\qquad K=V+V^\dagger.
\]

The auxiliary charges of the matter pair and Higgs transition sum to five;
the memory contributes minus five. Therefore `[K,Q_total]=0`, and K
preserves physical fermion parity. Exact matrix evaluation gives `K^3=K`.
For the declared control angle `pi/4`,

\[
 e^{-i\pi K/4}=I+(\cos(\pi/4)-1)K^2-i\sin(\pi/4)K.
\]

Starting from the empty component state, the complete sixteen-outcome
detector distribution has probability one half at the empty state and
one half at the pair/Higgs/memory state; every other outcome has probability
zero. Both populated states have total charge zero. This is a finite
charge-conserving interaction, not a Standard Model scattering calculation.

A charge eigenstate of the memory has `expectation(X)=0`; a declared coherent
superposition has `expectation(X)=1/2`. Treating that expectation as a fixed
Yukawa coefficient requires a controlled approximation to the joint memory
dynamics. In this finite controller `X^2=0`, so the transition saturates;
it does not supply an inexhaustible constant background. The coupling
strength, memory preparation, Higgs state and their long-time response are
inputs of this extension. No charged-lepton ratio or GeV mass follows from
charge compensation alone.

In particular, the coherent memory superposition is not prepared from a
fixed-total-charge neutral state by a charge-conserving unitary acting on
that memory alone. A preparation relative to another charged reference must
include that reference and its charge records. The complete charge-zero
pair/Higgs/memory evolution above has no such omitted charge reservoir;
it is not replaced by the coherent-memory expectation in the certificate.

## Finite theorem interfaces

`Geometry/CanonicalLorentzModule` supplies the Hermitian two-by-two Lorentz
quadratic module. The routing Clifford matrices and current transformations
give a concrete free-field representation of that structure, with a stated
momentum window and error. `Geometry/MassShellKinematics` takes a mass as
input; it does not choose v or a Yukawa coefficient. `Screen/FermionSectorAssembly`
supplies algebraic charge-compatible pairings on its declared bundle; the
source-specific spin/current and mass tensors here identify which algebraic
premises a field extension can instantiate. None of those results promotes
an auxiliary occupation mask to a native propagating particle without a
state, clock and dynamical intertwiner.
