# One finite defect, scalar and charged-observer construction

This construction chooses one quantum law. It joins the actual compact
defect hopping to the native thirty-mode matter/Higgs algebra and complete
charged observer instruments. It supplies no claim that those architecture
choices are forced, and uses no individual particle masses or family matrix.

## Complete source and clock

The primary geometry is the open cube `[-3,3]^3`, with 343 sites, the compact
texture and hopping of [compact_defect_pair](../compact_defect_pair/PROOF.md),
and texture amplitude `M=64`. Each site has the actual spin-major native
CAR module `C2 tensor R15`. Mode order is `30*x+15*spin+internal`. The full
fermion factor is `F(C^(10290))`; every occupation remains in the source.

The scalar factor retains the vacuum and both particle/antiparticle doublets
at every site, with total scalar occupation at most one across the finite
member. Its dimension is `1+4N=1373`. This global occupation cutoff is a
declared architecture, and differs from a product of one-particle cutoffs at
each site. Its fields are the literal matrix units

`phi_xk=|0><particle_xk|+|antiparticle_xk><0|`.

Use the same defect amplitude and graph to define the positive stiffness
`K=I+L_graph+M^2 diag(|n_x|^2)`. The scalar energy is `Omega=sqrt(K)`, the unique
positive matrix square root. Thus

`H_H=diag(0,Omega tensor I2,conjugate(Omega) tensor I2)`.

The positive square root is an explicit finite algebraic service. Large
displayed columns evaluate it numerically and retain the matrix residual;
that residual is a diagnostic, not an outward error certificate. The
[scalar package](../defect_scalar_modes/README.md) supplies every frequency,
its complete vacuum response and its separately identified reference.

At site x, let `P_channel,x,k=sum_ij C_channel,k[i,j]
c_(x,0,i)c_(x,1,j)` use the exact retained native tensors. Set
`D_x=(-1)^(x1+x2+x3)`. The full Hamiltonian is

`H=dGamma(h tensor I15)+H_H+
 sum_x,channel,k D_x(P_channel,x,k Phi_channel,x,k+h.c.)`.

Here Phi is phi for up, and phi-dagger for down and charged lepton. Every
primitive coefficient has modulus one. The fixed stagger phase makes the
spatial pair exactly the existing compact kernel `C=D i sigma2`; it is not a
family coefficient or a profile fitted to a mass. The sole clock is
`exp(-itH)` in these source units.

The executable supplies any complete occupation/scalar column, including
all Pauli-blocked alternatives and scalar transport. It uses no truncation
of fermion number. Full thirty-mode normal-ordered identities verify
Hermiticity and the twelve native current commutators. The hopping is the
identity on native internal rows, while H_H is the identity on the doublet
indices. The local tensor Ward identities therefore prove invariance under
the same twelve GLOBAL native currents on the entire spatial Fock space.
This is global charge covariance; quantum links, local Gauss constraints
and the interacting spectra of the separate chiral quantum-link source are
not imported.

If d is the maximum graph degree and L is the sum of absolute coefficients
of the actual local Yukawa matrix-unit words, a complete norm bound is

`||H|| <= 30N(M+d/2)+sqrt(1+M^2+2d)+N L`.

When the texture is switched off, its two M terms are omitted. The bound
uses `||dGamma(h)||<=30N||h||`, the graph row-sum bound for K and the norm-one
CAR/scalar factors of every pair word. It applies to all occupations and
all scalar states. Its loose size does not claim efficient simulation.

## One preparation, reference and complete readout

Let `P=1_{0}(sum_A Q_A^2)` project onto the global native singlets. The empty
fermion/scalar state proves `rank(P)>=1`. Since H commutes with every Q_A,
it preserves P. Let P0 be the full lowest-energy projector of H restricted
to P, and `rho=P0/rank(P0)`.

A declared finite preparation service first prepares `vec(P)/sqrt(rank(P))`
on the source and retained conjugate reference. The source and reference
are separately native singlets, and matching their fermion parities gives
even combined parity. The reference Hamiltonian is zero. A complete source
energy measurement has probabilities

`p(E)=rank(P P_E)/rank(P)`.

Every energy is an ordinary record. The ground branch has its original
probability `rank(P0)/rank(P)` and joint vector `vec(P0)/sqrt(rank(P0))`.
The other preparation records are retained. This definition is a finite
spectral preparation algorithm for this H; the full spatial ground and
those numerical ranks are not executed.

For an orthonormal set of m charged CAR coordinates a_mu, the observer
isometry is

`V=(sum_mu a_mu^dagger tensor |mu dual,+>
       +a_mu tensor |mu,->)/sqrt(m)`.

Its Gram is identity by the complete CAR algebra. Its dual addition and
direct removal reference representations cancel the actual native charge,
and their odd grading cancels the matter parity change. The two branch
labels and every energy outcome are retained while the internal reference
coordinates remain quantum. The whole two-point matrices are

`C+_mu,nu(E)=Tr rho a_mu P_E a_nu^dagger`,

`C-_mu,nu(E)=Tr rho a_nu^dagger P_E a_mu`.

Their summed CAR weight is exactly identity. The P_E are the spectral
projectors of the same full H, including states reached through the
fermionic exterior. A charged reference can leave the source-only singlet
sector. Energies relative to its constrained source ground are therefore
not assumed positive before the charged sectors are calculated.

The scalar instrument uses every `phi_xk/sqrt(4N)` and
`phi_xk^dagger/sqrt(4N)`, with coherent compensating reference coordinates,
and the no-operation Kraus operator

`K0=sqrt(1-1/(4N))*P_scalar_occupied`.

Indeed `sum_xk(phi_xk^dagger phi_xk+phi_xk phi_xk^dagger)
=4N P_vacuum+P_scalar_occupied`, so all effects sum to identity. A genuine
charged scalar pulse is
`X=sum_xk(phi_xk^dagger tensor |x,k dual><0|+h.c.)`.
It commutes with the total native currents including the dual apparatus
representation. Neither the scalar cost nor a coupling constant is
differentiated in this field insertion. The original preparation reference
is untouched and retained.

## Exact local member and the necessary control switch

The one-cell member at the origin has the texture explicitly disabled,
`h=0`, `K=Omega=1`, and `D=1`. Its full thirty-mode normal-ordered polynomial
is exactly the all-three-unit local source in
[local_reference](../coherent_defect_wave/chiral_joint/local_reference/PROOF.md).
The identification preserves every scalar coordinate and CAR generator;
it imports no many-cell ground or particle spectrum.

Only the distinct parameter member `y_u=y_d=0,y_e=1` has the empty quark
factor as a reducing subspace. The embedding of the six lepton modes is
`(9,10,11,24,25,26)`. Every one of its 320 columns equals the
[finite Higgs readout](../yukawa_defect_vertex/PROOF.md) source. An explicit
all-unit blank-column control produces nonzero quark pairs, rejecting the
same reduction when the other channels are enabled. The `14/5`
scalar moment and its ground reference are not transferred to another
parameter or geometry member.

## The complete quantum contact vertex

Normal ordering on the actual thirty-mode algebra verifies all 900 pair
identities. On the spatial source they read

`{[c_(x,s,a),H],c_(y,t,b)}
 =-delta_xy D_x epsilon_st sum_channel,k
   conjugate(C_channel,k[a,b]) adjoint(Phi_channel,x,k)`.

The kinetic and pure scalar terms vanish from this anomalous
anticommutator. Thus it is the actual scalar-valued equation-of-motion
vertex of the full Hamiltonian, independent of a supplied family matrix.

Let W be the complete four-dimensional compact one-body band, including
both signs. The same-source CAR coordinates are
`a_(alpha,i)=sum_xs conjugate(W_xs,alpha)c_(x,s,i)`. Their metric is identity
on all sixty coordinates. This is an injective CAR map, but the local
quantum interaction does not preserve its dynamics by itself.

With `Y_x=W^dagger P_x C conjugate(W)`, the full band contact is exactly

`K_band=-sum_x,channel,k Y_x tensor conjugate(C_channel,k)
                          adjoint(Phi_channel,x,k)`.

Every site tensor is retained. Expanding all scalar modes as
`phi_xk=sum_j u_xj phi_jk` gives `Y_j=sum_x u_xj Y_x`. The numeric packet
contains all 343 mode tensors and all site tensors; no single mode is
selected as a substitute for the quantum field.

For the separately specified free scalar vacuum, the complete response is
`D_phi(z)=2Omega(z^2-Omega^2)^(-1)`. Its canonical coordinates satisfy
`phi_j=sqrt(2omega_j) q_j`, so the corresponding contact tensor is
`sqrt(2omega_j)Y_j`. This normalization is derived on that free reference.
On the interacting ground, both `<[phi,phi^dagger]>` and
`<[phi,[H,phi^dagger]]>` must instead be obtained from the complete
interacting records before assigning a spectral normalization or
amputating a vertex. No free-vacuum moment is substituted for them.

The complete exterior contribution is retained through

`L[f]=(I-WW^dagger) diag(f) C conjugate(W)`,

`L[f]^dagger L[f]
 =S^dagger [sum_x f_x^2 B_x-(sum_x f_x B_x)^2]S`,

where `B_x=W^dagger P_x W` and `S=W^dagger C conjugate(W)`. The executable
checks the literal full-row expression for every scalar profile. The
nonzero leakage rules out an exact replacement of this interacting source
by its isolated band compression. Full two- and three-point functions
must propagate through the same complete source Hamiltonian.

These identities join the source and determine its contact tensors. They
do not evaluate the full spatial interacting ground, its charged spectral
poles or the frequency-dependent one-particle-irreducible vertex. A fixed
texture also breaks continuous translation: its trapped energies are not
rest masses without a separate propagation and observer-clock calculation.


## A complete finite-time propagation control

A second declared member has three sites `(-1,0,0),(0,0,0),(1,0,0)`,
the texture disabled, all ninety fermion modes, thirteen scalar states and
all three native Yukawa coefficients one. Its exact positive scalar root is

`Omega = [[4+3sqrt2,-2,4-3sqrt2],[-2,10,-2],[4-3sqrt2,-2,4+3sqrt2]]/6`.

Squaring gives the literal `I+L` of this source and its three eigenvalues
are `1,sqrt2,2`. The source columns execute integer arithmetic in
`Q(sqrt2,i)`, with common denominator `6^n` for `H^n`. Coincident rows
combine exactly. No fermion-number sector, Pauli-blocked alternative or
scalar output is dropped. The same full source columns are independently
evaluated numerically and compared through the third power.

Start with the fermion/scalar vacuum and a blank apparatus. The invariant
pulse

`X=sum_k(phi_left,k^dagger tensor |dual k><blank| + h.c.)`

acts on its complete reached space as
`[[0,1,1],[1,0,0],[1,0,0]]`. At pulse angle `pi/(2sqrt2)` it prepares
`(sum_k |particle_left,k> tensor |dual k>)/sqrt2` with probability one,
up to a common phase. The source and dual reference form a global native
singlet. The apparatus then has zero Hamiltonian and every reference
coordinate is retained. This pulse is an explicit preparation intervention
in source clock units; the prepared state is not a ground state.

The receiver projects onto exactly two fermions at the right site and a
blank scalar, summed over every internal pair and retained apparatus
coordinate. That projector and its complete complement are the two
classical outcomes. The projector commutes with every global native
current. Weak covariance makes its probabilities equal for the two
normalized prepared components, so one component suffices for the reported
calculation.

For that component, the exact supports of `H^n psi` for `n=0,1,2,3` have
`1,11,664,7799` coordinates. Their squared norms are

`1,10,827,(289759+2792sqrt2)/3`.

The receiver sees zero at powers zero and one, and

`||P_receiver H^2 psi||^2 = 8(2/3-sqrt2/2)^2 > 0`.

This is a scalar-transport and pair-creation path on the actual Hamiltonian.
The all-unit evolution also reaches the quark coordinates; the six-mode
lepton reduction is not used. At `t=1/10000`, evaluate
`psi_2(t)=psi-itHpsi-t^2 H^2psi/2`. The complete Hermitian spectral calculus
gives the state-specific remainder

`||exp(-itH)psi-psi_2(t)|| <= r=t^3 ||H^3psi||/6`.

With `a=t^2||P_receiver H^2psi||/2`, the exact yes probability lies in
`[(a-r)^2,(a+r)^2]` because exact algebra verifies `a>r`. The displayed
numeric interval is approximately `[2.70e-19,3.90e-19]`; the receipt
retains the exact radical endpoints. The complementary no probability
has endpoints `1-(a+r)^2` and `1-(a-r)^2`, even where machine arithmetic
rounds them to one. Neither branch is discarded or normalized separately.

This executes a finite-time encounter in the specified source family. It
does not supply a spatial ground, an asymptotic particle pole, a rest mass
or a relativistic propagation law. In particular the square-root scalar
energy couples the two outer sites directly; no finite-speed assertion
is inferred from this event.

## Bosonic operators and the complete free scalar control

The bosonic carrier supplied by this source is the scalar carrier described
above. The pure scalar suboperator has one energy-zero blank state and every
one of the 343 eigenfrequencies of Omega with multiplicity four, for the two
particle and two antiparticle components. These ranks total 1,373. Its 686
complex field coordinates and their adjoints give the complete free response
`D_phi(z)=2Omega(z^2 I-K)^(-1)`, equal-time commutator expectation zero, and
first energy moment `2Omega`. The code reconstructs the response and moment
from every retained mode and compares them with the literal joined K and
Omega. The normalized free field is `q=(2Omega)^(-1/2)phi`, with response
`D_q(z)=(z^2 I-K)^(-1)`.

These statements use the scalar blank for H_H alone. The all-unit Yukawa
Hamiltonian acts nontrivially on the matter/scalar blank, so that free
reference is not substituted for the full interacting ground. The scalar
suboperator frequencies are not identified with physical rest masses.

The twelve native currents are global charge matrices. This source has no
dynamical gauge-link factor, gauge-field oscillator or gauge-field energy.
The separate finite quantum-link program uses 60-state link factors,
covariant quantum paths, link energy and tensor-product local scalar
carriers. Those are different source operators and preparation domains.
No map preserving their complete Hamiltonian, retained reference and clock
into this joined member is supplied. The exact one-cell native-polynomial
identity does not transfer a spatial quantum-link spectrum. Consequently
this source supplies no W, Z, photon or gluon mass record; global charge
covariance is not a derivation of those dynamical sectors.
