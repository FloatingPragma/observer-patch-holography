# One finite chiral gauge and Higgs Hamiltonian

The declared program `coherent-spectral-quantum-link-Higgs-v1` uses the
[sitewise native source](../chiral_source/PROOF.md), the coefficients of the
[spectral derivative](../chiral_regulator/PROOF.md) and the
[quantum Higgs tensors](../chiral_higgs/PROOF.md). It has one Hamiltonian,
one physical Gauss subspace and a complete filter preparing a stationary
state of that same Hamiltonian. Its state and parameters are declared
architecture. A free background spectrum is not substituted for this
quantum gauge system.

## 1. Complete finite carrier and exact Hamiltonian

For an even side N take the N-cubed spatial cells of the spectral regulator,
with its antiperiodic matter spin structure and physical length L. Each cell
contains the full `F(C^2 tensor R15)` source code, dimension `2^30`, and the
five-state scalar memory `|0>,|H_1>,|H_2>,|antiH_1>,|antiH_2>`.
Every coordinate nearest-neighbor edge `(x,x+e_i)` with `x_i<N-1` carries
the sixty-state quantum link from the same source proof. The gauge edges
form a finite coordinate grid. The spectral fermion boundary convention is
in its coefficients; no extra periodic gauge edge is assumed. At side two,
the full logical dimension is `2^240 * 60^12 * 5^8`. All these factors and
occupation sectors are retained in the definition.

Let P_i(x,y) be the actual dense spectral momentum matrix along axis i,
tensored with Kronecker deltas in the other coordinates. For each unordered
pair with `x_i<y_i`, use the monotone coordinate path from x to y. If the
path has l edges, define in each actual R block

`W_xy,ab=sum_(c_1...c_(l-1)) U^1_(a,c_1) tensor ... tensor U^l_(c_(l-1),b)`.

The U factors are the literal link raising matrices `|a,b><0|`, not unitary
background values. The reverse transporter is `W_yx=W_xy^dagger`. Define

`H_kin=sum_(axis,unordered x,y,spin u,v,a,b)
 c_(x,u,a)^dagger [-sigma_i[u,v] P_i(x,y)] W_xy,ab c_(y,v,b) + h.c.`.

The matrix P is Hermitian. The formula includes each unordered pair once
and retains the adjoint term, so H_kin is Hermitian without replacing or
normalizing any link. It is spatially nonlocal: its support includes every
cell and edge along each path. Substituting classical identity transporters
as a separate parameter calculation reproduces the spectral derivative,
but the quantum link has no state with identity raising operators. Its
blank link state is not an identity-background vacuum.

At each cell take the full quantum Yukawa operator

`H_Y=sum_(channel,k) y_channel A_channel,k Phi_channel,k + h.c.`,

`A_channel,k=(1/2)c^T(epsilon tensor C_channel,k)c`.

Here Phi is the quantum field phi for the up channel and phi-dagger for the
down and charged-lepton channels, using the actual retained invariant tensors.
The source code transfers this even local operator with its spatial support.
With `N_H=diag(0,1,1,1,1)` and `P_link=I-|0><0|`, the complete selected clock is

`H=H_kin+sum_x H_Y,x+omega sum_x N_H,x+kappa sum_edges P_link`.

The receipt specifies the unit parameter example
`(y_u,y_d,y_e,omega,kappa,beta)=(1,1,1,1,1,1)`. Other exact parameter choices
or cell-dependent switches are declared controls, not source-selected mass
predictions. No plaquette magnetic term is implicit in this definition.
Adding such a term would require its own specified coefficients and proof.

## 2. Actual local Gauss invariance

The total Gauss generator at a cell is the sum of its native matter current,
its quantum scalar current, and the oriented endpoint currents of every
incident quantum link. A link transforms with the dual representation at
its left endpoint and the direct representation at its right endpoint.
At an internal path vertex, the two generator actions cancel in the contracted
index. At an endpoint the remaining link action cancels the commutator with
the matter creation or annihilation operator. Hence every path term commutes
with all local total Gauss generators.

The executable expands the actual one- and two-edge raising tensors into
59 and 279 matrix-unit terms and checks all twelve currents at every
endpoint and intermediate vertex: sixty exact Ward identities. The same
index cancellation proves any finite path length. Independence of the
normally ordered matter bilinears extends the identity to full Fock space.
The local Higgs Ward identity and scalar-number invariance are those executed
on the same Cauchy order and current matrices in the Higgs proof. The excited
link projector commutes with both endpoint currents. Thus the entire H is
Hermitian and gauge invariant on one finite carrier. Sitewise encoding
transfers it to actual native source currents. The grading and all of its
couplings retain their declared architectural ancestry.

## 3. A physical reference prepared by its own clock

Let P be the orthogonal projector onto the joint kernel of all local total
Gauss generators. It exists on this finite-dimensional carrier. The empty
matter, blank links and scalar vacuum form a nonzero vector in this kernel,
so `rank(P)>=1`. This vector need not be stationary: the quantum Yukawa term
creates gauge-compatible matter/scalar excitations. Gauge invariance gives
`[P,H]=0`.

Choose beta positive. The physical Gibbs state is

`rho_beta=P exp(-beta H) P / Tr(P exp(-beta H))`.

It is faithful on P, stationary for exactly `exp(-it H)`, and generally
not faithful on the unconstrained ambient carrier. This definition asserts
neither an efficient state-preparation algorithm nor a unique physical
vacuum. The zero-temperature state is the normalized full ground projector
of H restricted to P; no uniform gap or unique ground is inferred.

A complete finite preparation is available entirely within the physical
Gauss sector. Set `r=rank(P)`. Since P and H preserve matter parity, choose
a parity-homogeneous orthonormal basis p_j of P. Use a neutral graded
r-dimensional ancilla whose |j> has the same parity as p_j. Its initial
|0> is even, matched to the empty source vector. The known empty gauge singlet Omega can be
sent by a specified finite unitary to
`(1/sqrt(r)) sum_j p_j tensor |j>`, with the ancilla retained. Both input
and output lie in the trivial gauge representation and have even total
parity. Extend this map unitarily within the total-even part of
`P tensor H_ancilla`, and as the identity on its orthogonal complement.
This is a gauge- and parity-preserving purification of the input `P/r`,
with no assumption of access to unphysical source states. Computing the
basis and controller is a finite but potentially expensive specified service.

Let B be a proved bound with `||H||<=B`. Use two Kraus operators

`L=P exp[-beta(H+B I)/2]`,

`F=sqrt(I-L^dagger L)`.

Since `[P,H]=0` and `-B I<=H<=B I`, the effects obey
`0<=L^dagger L<=P<=I` and `L^dagger L+F^dagger F=I`.
The success branch is exactly rho_beta. Its original probability is

`p_success=Tr(P exp[-beta(H+B I)])/r >= exp(-2 beta B)>0`.

If an ambient `I/D` input is admitted in a comparison controller, the same
success probability instead has denominator D and lower bound
`rank(P) exp(-2 beta B)/D`. That ambient control is not required by the
physical-sector preparation.

Failure is retained as an ordinary record. Since L and F are commuting
positive operators, the pointer block `[[L,-F],[F,L]]` is an explicit unitary
dilation of `|0> L+|1> F`. Every block commutes with the Gauss generators
and matter parity, and the pointer is neutral and even.
It acts through the declared source code and has no implicit postselection
renormalization of complete records. This is a mathematical finite controller
using spectral functions of the specified H and Gauss projector, without a
claimed efficient circuit or a source-pulse count. The global exponential
and its diagonalization are not numerically materialized in the receipt.

For unit parameters a safe explicit B follows from the triangle inequality.
An l-edge kinetic path in an irreducible block of dimension d has
`d^(l+1)` elementary internal-index terms. Each raising matrix and each
fermion operator has norm at most one. Two Pauli entries and the adjoint
give bound `4 |P_i(x,y)| sum_r d_r^(l+1)` per unordered path. The local
Yukawa bound is `2 sum_(channel,k,i,j)|C_channel,k[i,j]|`, and each scalar
number or link electric projector contributes one. Their sum at side two,
L=2pi, is B=1884. It is conservative and valid on the entire logical carrier,
including all occupations and link/scalar alternatives. In particular it
also bounds H restricted to P.

## 4. An executed encounter of quantum creation and propagation

A small control uses one actual side-two spectral edge. Its right cell has
the charged-lepton quantum Yukawa coupling one, its left coupling is zero,
and both quark Yukawa couplings are zero. The full internal representation
is retained; empty quark occupations and unexcited colored links form a
reducing spectator sector because their hopping conserves quark number.
The left scalar vacuum is also reducing. The active source has twelve
lepton CAR modes, six link states and five right-cell scalar states, hence
122880 dimensions. Its link states are the blank, the e^c singlet channel
and four weak-doublet raising channels. Both spin components and every
active fermion occupation are retained.

The code applies that literal Hamiltonian as sparse exact columns. Let
Omega be its empty matter/link/scalar vector and R the complete projector
onto at least one lepton at the left cell. Direct execution gives

`R Omega=R H Omega=0`, `||H Omega||^2=4`, `||R H^2 Omega||^2=3`.

The first term creates a matter pair and scalar at the right cell; the
second lets that same Hamiltonian transport a created lepton through the
quantum gauge link. All reverse source coefficients in the populated
neighborhood are checked for Hermiticity. This calculation does not reset
a link or insert a supplied classical Higgs background between the steps.

On the entire active carrier, its ten raising hopping terms and their
adjoints have total norm bound ten, its four pair terms and their adjoints
have bound eight, and the electric/scalar projectors add two. The bound
24 therefore covers the full H. At t=1/1000000, the second-order unitary
Taylor remainder is at most `delta=t^3*24^3/6`. The complete arrival
probability obeys

`p_yes >= (t^2 sqrt(3)/2-delta)^2 >0`, `p_no=1-p_yes`.

The receipt retains the exact algebraic lower bound and the complementary
outcome. These records are a controlled finite quantum encounter, without
an asymptotic scattering or particle mass interpretation. The empty input
here differs from the physical Gibbs reference; the claim identifies that
preparation explicitly and does not relabel it a ground state.

## 5. Scope of the common construction

The selected H, physical subspace, Gibbs filter and clock define a single
finite interacting gauge/Higgs source family. Its cutoffs, arbitrary-range
paths, coefficients and temperature are specified inputs. The same source
supports complete observables `Tr(rho_beta U^* E_j U)` for any declared
complete instrument, with Gauss-compatible source operations and apparatus.
Calculating these huge full-family distributions is separate from the exact
local Ward identities and the finite encounter executed here.

No proof transfers the free spectral filled sea to this quantum-link
Hamiltonian. A controlled interacting Lorentz limit, positive zero-temperature
field reference, nonzero connected correlations under refinement, physical
Higgs phase and source-selected couplings require their own estimates on
this H. Its finite Gibbs reference is not a vacuum mass theorem. Likewise
its antiperiodic finite-volume energy levels are not rest masses. The
conditional constant-background particle dictionary uses a different scalar
specification and retains that boundary. This construction supplies no
measured Standard Model or GeV prediction.
