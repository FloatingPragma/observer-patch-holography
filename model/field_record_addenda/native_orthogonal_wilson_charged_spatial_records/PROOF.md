# Charged local records with a receding compensated reference

This result extends the frozen orthogonal Wilson spatial family at fixed
spacing a=1/m and g=1/64. Its actual Psi, full native occupation code,
Gamma regulator, raw beta-sea scalar, photon completion, Gauss rule and
all references remain. In each open box the root is static during the
whole experiment. The observer stays far from every boundary face, and
the root r_n recedes from the observer. The original corner-root family
is one example. Comparing another root is a comparison of finite sources
with a different static reference location; no finite-volume equivalence
of their charged spectra is asserted.

The accepted initial source is a localized charged excitation of the
prepared perpendicular free sea. It is a nonstationary input for the
actual coupled H_n. We prove a reference-independent limit of complete
local gauge-invariant records at bounded time and a reference-subtracted
**prepared-state energy excess**. Neither result is a charged ground-state
threshold, isolated pole, asymptotic particle mass or scattering theorem.
The full-native global code guard remains, so locality is a statement on
the admitted code and its complete preparation flags.

## A uniformly admitted local charged packet

Use the parent's complete code-bank assembly, all original references,
physical routing and algebraic perpendicular-sea rotation. On the copied
assembly-success flag the active state is exactly Omega_n, the free
perpendicular Slater sea, times transverse photon vacuum, with the
normalized Gauss labels. The original assembly failure probability is at
most f_n=V(1-p0^4)^R, on the prescribed schedule tending to zero. No
interacting-ground projection is supplied.

Choose one interior vertex x0 and a fixed unit spinor v_+ with beta v_+=v_+.
Let f_+ be v_+ at x0 and zero elsewhere. Choose similarly beta v_-=-v_-.
The compensated probes on the physical code are

`D_+=Y c_(f_+)^dagger Y^dagger`,
`D_-=Y c_(f_-) Y^dagger`,

and zero on its retained complement. Y here means the same normalized
integer-Gauss construction Y_(n,r); changing r changes which longitudinal
label is omitted, not the canonical field or density frame.

Put P_n=1_(h_perp,n<0). In dimensionless units h=beta+A, {beta,A}=0,
h^2=I+A^2 and I<=|h|<=10I. Therefore |h|^-1 commutes with beta and its
spectrum lies in[1/10,1]. The A|h|^-1 term is off-diagonal in beta grading.
Consequently the two literal local pump probabilities obey

`w_+=<f_+,(I-P_n)f_+>>=11/20`,
`w_-=<f_-,P_n f_->>=11/20`.                          (1)

These are exact source bounds, with no orbital fitted to a desired mass.
Conditional on pump success the normalized matter state is the Slater
state with projector

`P_n,q=P_n+q u_n,q u_n,q^dagger`,
`u_n,+= (I-P_n)f_+/sqrt(w_+)`,
`u_n,-= P_n f_-/sqrt(w_-)`, q=±1.                    (2)

It has total source charge q and reference R=-q, with actual cost4|R|=4.
It is localized in the stated canonical observer coordinates, not asserted
to be a microscopic-native local operation through the global FW map.

The parent's sign-polynomial bound eps_J=100(99/100)^(J+1), of range2J+1,
controls localization: outside that radius from x0,
`||u_n,q||<=eps_J/sqrt(11/20)<2eps_J`. Local projectors and w_n,q converge,
and this uniform tail implies norm convergence of the packets after the
ordinary one-particle coordinate embeddings. It does not identify distinct
sharp physical Gauss fibers as nested Hilbert subspaces.

Implement the probes using the complete two-mode fermionic pointer and
source-first graded even pump of the charged-cluster package. The fair
recorded setting q uses `exp[(pi/2)(D_q f_q^dagger-f_q D_q^dagger)]`.
The pointer starts empty and has zero Hamiltonian. Its four occupation
outcomes are retained. On the even sea (N=2V), the success Kraus map is D_q
and the empty-pointer map is I-D_q^dagger D_q. The physical Gauss generator
Q+R and the distinct dual global flavor Q-n_f+ +n_f- are conserved, as is
joint fermion parity. The pointer is not the bosonic compensator, and its
label does not remove the source's reference cost4.

All assembly failures and empty-pointer outcomes execute recorded no-op
subsequent tests. Thus the original charged success union is at least
`(1-f_n)11/20`, and each fair signed success at least `(1-f_n)11/40`.
Pump failures are part of the limiting law; their probability need not
vanish. They are not discarded or included in f_n.

## The three-dimensional monopole estimate needed for the reference

Let L_n=D_n^T D_n and G_n=L_n^+ for the actual open Neumann box, with
N vertices per axis and V=N^3. Write p_t=exp(-tL_n). The exact product
cosine basis has eigenvalues at least4|k|^2/N^2 and squared coefficients
at most8/N^3. For u=4t/N^2 and
`S(u)=sum_(k>=1)exp(-u k^2)`,

`0<=p_t(x,x)-1/V<=24 S(u)(1+S(u))^2/N^3`.

Indeed the nonzero three-index sum is at most3S(1+S)^2. For k>=1,
`exp(-u k^2)<=exp(-u/2)exp(-u k^2/2)`, and monotone integral comparison
implies `S<=c u^-1/2 exp(-u/2)`, c=sqrt(pi/2). Hence the heat bound is
at most

`3c t^-3/2 [u exp(-u/2)+2c sqrt(u)exp(-u)+c^2 exp(-3u/2)]`.

Using u exp(-u/2)<=2/e<1, sqrt(u)exp(-u)<=1/sqrt(2e)<1/2,
c<4/3 and c^2<11/7 gives a coefficient below328/21<16. Thus, uniformly
in the box, vertex and positive time,

`0<=p_t(x,x)-1/V<=16 t^-3/2`,
`|p_t(x,y)-1/V|<=16 t^-3/2`.                         (3)

The second inequality follows because exp(-tL)-P0 is positive. The
integral tail beyond T is at most32/sqrt(T). On[0,1] the diagonal is at
most1, so

`0<=G_n(x,x)<=33`, `|G_n(x,y)|<=33`.                 (4)

No bounded infinite-volume monopole inverse on all l2 vectors is claimed.
Its point kernel exists in three dimensions; its vector need not be l2.

For distance d=dist(r,x)>0, all powers L^j_(r,x) vanish for j<d. Since
||L||<=12, a direct exponential-series remainder proves

`|G_n(r,x)|<= T exp(12T)(12T)^d/d! + T/V +32/sqrt(T)` (5)

for every positive T. Thus G_n(r_n,x) tends to zero for every fixed
observer vertex as r_n recedes and V grows. This remains true when the
root is at an open face, edge or corner; no translation-invariant
boundary assumption is made.

For two fixed local vertices, finite-range powers agree with the full
cubic lattice before reaching the boundary. If their boundary distance
is at least b, a sufficient point-kernel error is

`2T exp(12T)(12T)^b/b! + T/V +64/sqrt(T)`.            (6)

First increase T, then b and V. This proves the local Green limit and
its root-independence. It also supplies a constructive, if conservative,
finite accuracy budget. The inherited dipole l2 theorem is stronger for
local electric tests: `G_n D_n^T s -> G_infinity D_infinity^T s` in l2.
Evaluating that convergence at any receding root makes its value tend to
zero. Equations(5)-(6) give a separate point-kernel control when the
energy expectation involves a monopole rather than a dipole.

## Charged equal-time and positive-time records

In the fixed Q=q success sector, the actual longitudinal electric test is

`g s^T D G rho_ref = X(a_n)-q a_n(r_n)`,
`a_n=g G_n D_n^T s`, `rho_ref=rho-q e_r`.             (7)

The additional term is a scalar. It tends to zero by the dipole theorem,
and its exponential differs from identity by at most |q a_n(r_n)|. It is
never omitted in a finite member. Full transverse E/B/Q covariances and
Weyl phases use the same finite cochains and the same three-dimensional
uniform heat estimates as the parent.

Charged-state expectations of a word A equal
`<Omega_n,D_q^dagger A D_q Omega_n>/w_n,q`. All matter insertions are
finite local CAR words. For a k-factor local CAR word A, the parent
centered density estimate extends directly to the extra endpoint:

`||X(delta a) A c_f^(dagger) Omega_n||`
` <=(3+k)||delta a||_2`.                             (8)

Here f is a fixed unit local spinor; the density commutator with it has
norm at most max_support|delta a|. Onsite spin rotations preserve this
identity. Normalizing divides by at most sqrt(11/20). All finite-word
expectations then converge by the local P_n limit, bounded finite Slater
determinants and the controlled l2 tails. No infinite determinant or
classical replacement of rho is made.

For positive times use exactly the parent's spatial comparison: remove
only hopping terms outside a buffer, retaining the entire quadratic
field, Coulomb energy, onsite matrices and every raw reference constant.
The uniform physical commutator majorant and its bounded-hop truncation
are independent of source state and of the fixed value q. They therefore
apply on this charged point fiber. Its gauge-parent use is a direct-integral
comparison, not a normalizable sharp-longitudinal state in L2(R^E).

In each retained free-baseline hopping word the Coulomb phase is exactly

`g^2 t[(D G rho)_e-q(D G e_r)_e+(D G D^T)_ee/2]`.     (9)

The new root-dependent factor is a scalar of modulus1 whose exponent
tends to zero. Finitely many retained edges and bounded times give a
uniform sum of these phase errors. The scalar root self-energy and cost4
commute with every gauge-invariant history and cancel from its record
probabilities; their raw phases remain in the source clock. The matter
word has the same additional endpoint pair D_q^dagger,D_q, so (8) controls
all density tails. Transverse Q covariance, free-time Weyl phases and
bounded current words converge by the same parent theorem.

At fixed spatial buffer R and finite Dyson order J, every retained word
therefore converges. The exact norm majorants bound all omitted hopping
and Dyson terms uniformly in volume. The literal finite-cap comparison
starts with an arbitrary normalized matter vector and photon vacuum, so
it applies to (2), not only the neutral sea. The pump is photon-number
preserving. Its inclusion in the complete instrument adds no displacement
radius and requires no replacement of exp(ig Q_K) by a compressed full
exponential. Use the parent's cap schedule with the observer's actual
total wait and field/current test strengths.

If d_n is that complete cap error, a sufficient two-volume bound has the
same form

`TV(law_n,law_l)<=f_n+f_l+d_n+d_l+2delta_R+2delta_J`
`                  +epsilon_words(n,l;R,J)`.        (10)

The last term now includes root phases, local projector/pump weights and
the two additional endpoint insertions. The law is on the fixed finite observer alphabet and includes the fair
setting, assembly-success copy, every pump outcome and every detector
outcome. The remote VR fine code flags remain retained in each global
program; their different-dimensional registers are not identified across
volumes or silently claimed to have this fixed-alphabet limit.
Unnormalized word expectations prove this directly, so nonvanishing pump
failures are included. Conditional signed-success laws also converge,
because their original probabilities stay at least(1-f_n)11/40.

Thus the actual complete charged local records have a root-independent
bulk limit at fixed spacing and bounded source time. The local physical
net consists of gauge-invariant observables. Charged endpoints require
the stated reference dressing and are not relabeled as local physical
observables. The result neither repairs full-native global-guard locality
nor substitutes this input for an interacting vacuum.

## A strictly nonzero complete local charge record

The zero-wait limiting law is locally distinguishable from the neutral
prepared sea. This density test is performed immediately after the pump,
with tau=0; its positive lower bound is not asserted after an arbitrary
positive autonomous wait. The separate bounded-time convergence theorem
does allow such waits.
For either packet, its density at the creation vertex is q p_x0. Equation(2)
gives |<f_q,u_q>|^2=w_q, so p_x0>=w_q>=11/20. Use the full finite density
spectrum[-4,4] to define complete effects

`E_±=(I±q rho_x0/4)/2`,

with their positive square-root Kraus maps and both outcomes retained.
They preserve Gauss, total charge and parity. On the successful code packet,
the signed outcome imbalance is p_x0/4>=11/80; in the neutral sea it is0.
This is a physical gauge-invariant local detector, not a bare charged CAR
observable. Finite functions of rho are finite local CAR polynomials, so
these instruments also satisfy the stated finite-word limit argument.

Before conditioning on any flag, each fair signed preparation contributes
at least `(1-f_n)121/3200` to the difference between its pump-success/plus
and pump-success/minus record probabilities. Summing both q settings gives

`Pr(assembly,pump success,plus)-Pr(assembly,pump success,minus)`
` >=(1-f_n)121/1600`.                               (11a)

All failed preparations and all other pointer outcomes are still present
in the law. The positive limiting imbalance is not an interaction witness:
it detects the source charge already present at g=0. The inherited
nonzero-coupling current witness is a different experiment.

## A reference-subtracted prepared energy, not a spectral mass

Let E_n,q^r be the expectation of the actual H_n,r in the successful
charged input(2) with photon vacuum, and E_n,0 its expectation in the
neutral free sea with photon vacuum. The source preparation uses exactly
these states, but they are not stationary coupled ground states. Native
baselines, the old sea scalar, photon zero point and code completion all
cancel in their difference. The retained reference contributes exactly4.

For p_x=||u_x||^2, Slater Wick contraction gives the exact identity

`<rho_x rho_y>_(P+q uu^dagger)-<rho_x rho_y>_P`
` =q[delta_xy p_x-2 Re(u_x^dagger P_xy u_y)]`.        (11)

The rank-one covariance term -p_x p_y cancels the new mean-density product
+p_x p_y. This cancellation is important; treating the packet as a classical
charge density would give a different answer.

Let G odot P denote the block matrix G_xy P_xy, and let hbar_n,K be the
actual photon-vacuum compression of the one-body kinetic operator. Finite
photon parity makes each link expectation r_e,K real, so hbar has the
original onsite matrices and links r_e,K B_e. Direct contraction yields

`E_n,q^r-E_n,0-4-(m g^2/2)G_n(r,r)`
` =q<u,hbar_n,K u>`
`  +(m g^2 q/2)[sum_x G_n(x,x)p_x-2<u,(G_n odot P_n)u>]`
`  -m g^2 sum_x G_n(r,x)p_x`.                        (12)

The subtractions are displayed reference contributions to an expectation
comparison, not changes in the Hamiltonian. Without the self-energy
subtraction, different corner/face/interior reference environments need
not yield the same energy excess. No claim of reference-independent
finite raw charged eigenvalues is made.

The last term tends to zero: truncate the packet to radiusL, apply(5) on
that finite support, then bound its probability tail by33||u_out||^2 using
(4). The diagonal and exchange terms converge by(6), the local projector
limit and the uniform packet tail. In particular Schur multiplication by
positive G is positive and

`0<=G odot P<=diag(G_xx) tensor I_4<=33I`.

If D=sum G_xx p_x and A=<u,(G odot P)u>, then0<=A<=D<=33, so the absolute
value of the intrinsic Coulomb term in(12) is at most(33/2)mg^2. This
uniform operator estimate avoids an unjustified l1 tail for the packet.

The same field heat estimate gives every edge's vacuum coordinate variance
v_e<13 uniformly in volume. For every literal cap K>=1, its second moment
is exactly v_e and parity gives
`0<=1-r_e,K=<1-cos(gQ_e,K)><=g^2 v_e/2`.
Thus `||hbar_n,K-h_perp,n||<=39mg^2`, by degree at most6. Also
`||hbar_n,K||<=10m`, since |r_e,K|<=1. Equation(2) lies in the positive
or negative free band, respectively, with |h|>=m. The limiting intrinsic
prepared energy therefore obeys

`m[1-(111/2)g^2] <= lim(reference-subtracted excess)`
`                  <=m[10+(33/2)g^2]`.              (13)

The lower endpoint is0.9864501953125m at g=1/64. It is a strictly positive
source-prepared energy cost, not a lower bound on every charged spectrum
or an identified rest mass.

For the cap limit, with r=g sqrt(13), the inherited literal-displacement
bound gives

`|r_e,K-r_e,infinity|<=b_K(r)`
` =r^(K+1)/sqrt((K+1)!)*(1+sqrt(K+1)r)`.

The kinetic energy difference in(12) is at most6m b_K(r), uniformly in
volume. Coulomb terms and the matter packet do not depend on K. Therefore
this energy limit is compatible with the actual increasing-cap schedule.
The K=3 allowance is explicitly evaluated as an upper bound; no finite
CCR or exact finite Gaussian characteristic function is assumed. The
quantity hbar is used only to compute this prepared expectation. It is
not promoted to an interacting self-energy or charged dispersion.

For each finite member these expectations also have complete bounded
source instruments: if B_n bounds the full raw H_n, use the two effects
`(I±H_n/B_n)/2`, with their positive square-root Kraus maps. They preserve
Gauss and all source sectors, retain both outcomes and have expectation
B_n(p_+-p_-). Run them on the neutral or successfully pumped preparation,
retaining every earlier flag. This is an exact finite matrix-service
energy readout. No uniform-resolution local spectral detector or charged
pole is thereby implemented; the global B_n and statistical accuracy
needed to resolve a finite difference are explicit resource boundaries.

## Complete finite clocks and capacity

There are three inherited preparation slots, one added fair-setting/even-
pump slot, and r declared detector slots, plus the actual autonomous wait
T=tau/m on every bank and every failure branch. Total source time is
`4+r+tau/m`. Every driven slot cancels the actual H_n,r background with a
specified finite Hermitian logarithm. Neither the reference energy nor a
partner bank's clock is paused or reset for free. The optional single
energy readout has r=1 and tau=0, hence time5.

For R banks the parent's full source capacity is retained. Add three bits
for the complete two-mode pointer and fair sign. Including local records,

`C=R[4V(256ln2+8ln5+24ln119+ln(2^31+1))`
` +(V-1)ln9+ln(8V+1)+ln binom(K+d,d)]`
` +VR ln2+(r+4)ln2`.

The finite observer alphabet has2^(r+4) labels; the remote VR fine code
flags remain in the full program. With the inherited revised raw bound
B_raw,new, each control has norm at most R B_raw,new+pi and gross switching
work at most2(4+r)(R B_raw,new+pi). Exact matrix services and additional
coefficient workspace remain declared architecture. No memory-budget-
selected volume, efficient passive repair circuit or interacting ground
volume limit is claimed. A separate c(V) completion changes complement
energies but leaves these admitted code laws unchanged under its stated
new background-cancelling controls.
