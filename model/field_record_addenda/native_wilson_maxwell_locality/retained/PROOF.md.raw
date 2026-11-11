# A gauge-invariant propagation bound on the same cube

The source is the positive macrocell-completed native Wilson–Maxwell cube
of `model/native_wilson_maxwell_ground_refinement`, with the same
`g=1/64`, source-selected mass `m`, spacing `a=1/m`, full native
complements, reference charges, code and autonomous clock. Its photon
cap changes only when explicitly stated. The comparison uses the same
five unbounded transverse oscillators admitted by that package. No
additional physical gauge coordinate or autonomous Hamiltonian is added.

The result concerns bounded gauge-invariant electric, magnetic, density
and current observables. Their labels refer to the declared macrocell
cube. Their Foldy-Wouthuysen/native implementations are spatially
extended and their global code guards are retained. Locality in this
observer algebra is not locality of every microscopic native tensor
factor. A charged endpoint is not an observable of this local physical
algebra without its gauge dressing or boundary reference.

## The physical algebra and its coordinate parent

Let `D` have edge rows with tail coefficient -1 and head coefficient +1,
let `C` be the actual oriented face boundary, and let G be the cube
Laplacian pseudoinverse. The source has `C D=0`,
`P_T=I-D G D^T`, and five transverse coordinates. On the code its
gauge-invariant fields are

`E=P_T p+g D G rho_ref`, `B=C q`,
`rho_ref=(I-e_root ones^T)rho`, `rho_x=Psi_x^dagger Psi_x-2`,
`F_e=Psi_head^dagger B_e^dagger Psi_tail exp(i g q_e)`.

The full 32-component Psi satisfies the canonical anticommutation
relations on the complete `16^8` occupation code. This follows from
the source's complete four-constituent code and its specified global
Foldy-Wouthuysen CAR action. No five-state hard-core carrier is identified
with a complete CAR algebra.

The exact identity `P_T+D G rho_ref-map D^T=I` gives
`[E_e,F_f]=g delta_ef F_f`. Also `[B_f,E_e]=i C_fe` and
`[rho_x,F_e]=D_ex F_e`. Thus the apparent Coulomb nonlocality of the
reduced coordinates cancels in the physical electric/current algebra.
Keeping only the transverse coordinate commutator would give a
different, nonlocal coefficient. The receipt checks every edge and face
with the actual source orientation.

For the proof, introduce a mathematical coordinate parent with twelve
canonical link pairs `(Q_e,P_e)` on `L2(R^12)`, the same finite CAR
factor and the same total-charge reference cost. In dimensionless
time `tau=m t`, its Hamiltonian, up to retained scalar constants, is

`H_par=(P^2+||C Q||^2)/2`
`      +sum_x Psi_x^dagger Omega_x Psi_x`
`      +sum_e(F_e+F_e^dagger)+4|R|/m`,                 (1)

where `Omega_x=(1+degree(x)/2)beta`,
`F_e=Psi_head^dagger B_e^dagger Psi_tail exp(i g Q_e)` and
`B_e=-(beta+i alpha_axis)/2`. Each B is a rank-two partial isometry,
so `||B||=1` and the complete many-body `||F_e||<=2`.
These identities are checked on the actual three spin matrices.

The reference cost is a function of total charge. It commutes with
all gauge-invariant fields, currents and densities and with the other
terms of (1). Their dynamics is therefore unchanged when this central
term is omitted during the following auxiliary commutator calculation.
It is retained in the physical source, energy preparation and clock.
This observation is needed because charged Psi itself does not commute
with the reference cost.

Choose orthonormal frames V and L for the transverse and longitudinal
edge spaces, respectively. Put

`Q=V q_T+L q_L`, `chi=G D^T L q_L`,
`T=exp[-i g sum_x chi_x rho_ref,x]`.

This is a well-defined unitary multiplication operator on the parent
coordinate space, with a finite matter matrix at each coordinate. It
satisfies

`T Psi_x T^dagger=exp[i g(chi_x-chi_root)] Psi_x`,
`T P_L T^dagger=P_L+g L^T D G rho_ref`.

The longitudinal phase cancels the fermion phase in every F. Under T,
the Gauss generator `D^T P-g rho_ref` becomes `D^T P`. Fourier
transforming `q_L` therefore gives a direct integral over `k in R^7`,
with fiber

`H(k)=H_T+H_matter(q_T)+|k+g L^T D G rho_ref|^2/2`
`     +4|R|/m`.

At k=0 the last quadratic term is the exact source Coulomb energy.
The electric observable becomes
`E(k)=V p_T+L k+g D G rho_ref`, and magnetic/current/density observables
become their source operators. The physical comparison is this point
fiber. Sharp longitudinal Gauss values are not normalizable states or
subspaces of `L2(R^12)`, and no isometry from a parent Gauss eigenvector
is asserted.

A uniform parent operator bound for decomposable gauge-invariant
observables holds almost everywhere in k. It also holds at k=0:
after removing the scalar `|k|^2/2`, `H(k)-H(0)` is a bounded finite
matter operator linear in k. Unitary Duhamel makes the evolution of
each bounded observable norm-continuous in k. Electric Weyl operators
have only the additional continuous scalar phase `exp(i s L k)`.
Currents, magnetic fields and densities have no extra k dependence.
Consequently each commutator norm is continuous at zero. An
almost-everywhere bound then extends to the physical point fiber.
This is a direct-integral comparison argument, not an extra physical
longitudinal register.

## A closed interacting commutator majorant

Fix a bounded unitary intervention U in the coordinate parent with
bounded initial commutators with Q, P and Psi. For positive time define

`q_e(tau)=||[Q_e(tau),U]||`,
`p_e(tau)=||[P_e(tau),U]||`,
`z_x(tau)=||column_alpha [Psi_x,alpha(tau),U]||`.

The four-component CAR column has norm at most 2. The corresponding
adjoint row has the same commutator norm: its row commutator equals
`U [column(Psi),U]^dagger U`. This uses U unitary, without replacing
row and column norms for arbitrary operator matrices.

The phase derivative is `J_e=i(F_e-F_e^dagger)`. The oriented charge
flow is `I_e=-J_e`, consistent with the parent's incidence convention.
Leibniz, the CAR column norm and `||F_e||<=2` give

`||[J_e(tau),U]||`
` <=4[z_tail(tau)+z_head(tau)]+4g q_e(tau)`.           (2)

Here `||[exp(i g Q),U]||<=g||[Q,U]||` follows by the exponential
Duhamel identity. The exact equations are

`dot Q=P`, `dot P=-C^T C Q-g J`,
`dot Psi_x=-i Omega_x Psi_x`
`           -i sum_(y adjacent x) B_xy W_xy Psi_y`.

Removing the onsite unitary spin rotation preserves the column norm.
Every hopping spin matrix has norm one. Thus, with S the unsigned
edge/vertex incidence and A the vertex adjacency,

`q(tau)<=q(0)+integral p`,
`p(tau)<=p(0)+integral[(|C^T C|+4g^2 I)q+4g S z]`,
`z(tau)<=z(0)+integral[A z+2g S^T q]`.                (3)

All inequalities are componentwise. Let M be the resulting nonnegative
`(2*12+8)` square matrix in the order `(q,p,z)`. Iterating (3), or
the elementary positive Volterra comparison, proves

`(q,p,z)(tau)<=exp(tau M)(q,p,z)(0)`.                 (4)

The unbounded-coordinate step in (3) has a domain justification. Split
(1) into its quadratic field Hamiltonian plus bounded matter hopping
and onsite matrices. Free Q/P evolve linearly. The interaction-picture
hopping has bounded commutators with every free linear Q/P, uniformly
on compact time intervals, because its P derivative is the bounded
current and its Q commutator vanishes. The norm-convergent bounded
Dyson series therefore also converges in the Q/P graph norms: commute
Q or P through each finite Dyson word, sum the inserted bounded
commutators, and use the exponential scalar majorant. This proves
domain invariance and the integrated Heisenberg equations. In
particular Q/P at time tau equal their free linear expressions plus
bounded current integrals, so the displayed commutators are bounded.
The bounded CAR equations hold as strong integrals with the displayed
operator-norm bounds. Exponential Duhamel can then
be applied on these invariant domains and extended by continuity.
No finite-cap canonical commutator is used in this proof.

For a field intervention
`U=exp[i(s_U.P+t_U.Q)]`, the initial data are
`q(0)=|s_U|`, `p(0)=|t_U|`, `z(0)=0`.
A later field Weyl V with coefficients `(s_V,t_V)` obeys

`||[V(tau),U]||<=|s_V|.p(tau)+|t_V|.q(tau)`.         (5)

There are bounded gauge-invariant matter extensions. For a density
test `exp(i theta rho_x)`, initial z is at most `2|theta|` at x and
zero elsewhere, with initial q=p=0. For a current test
`exp(i theta J_e)`, initial z is at most `2|theta|` at each endpoint,
initial p is at most `4g|theta|` at e, and q=0. At the receiving end,
`||[rho_x(tau),U]||<=4z_x(tau)` and (2) controls J. Exponentiating
these bounded generators and telescoping finite products gives
propagation estimates for their finite unitary words. Only the
gauge-invariant words are transferred to the Gauss fiber. A bare
charged Psi endpoint is an auxiliary calculation, not a local physical
observable.

The field/current/density algebra therefore has a same-source
geometric commutator majorant on this fixed cube. Matrix powers in
(4) are weighted paths in its edge/vertex incidence graph. The estimate
has continuous-time tails, not exact finite-speed support. It supplies
no Lorentz velocity or spatial continuum limit.

## An executed nonincident electric/magnetic setting

Use the actual edge `e=0=(000,100)` and actual face `f=3`, the xy face
at z=1. The source orientation gives `C[f,e]=0`. Let

`U=exp[i B_f/64]`, `V=exp[i E_e/64]`.

The initial majorant has only `p(0)=|C[f,:]|/64`. The receiving bound
is `p_e(tau)/64`. Every entry of M is rational at `g=1/64`.
`locality.py` computes its nonnegative Taylor terms through order 30
and bounds the omitted tail by the scalar infinity-norm exponential
series. The first two coefficients vanish exactly; the degree-two
coefficient before division by `2!` is `1/4096`. The coupled matrix,
including all matter-current terms, gives

`||[V(1/64),U]||<2.981081774e-8`.                     (6)

The table also gives bounds at dimensionless times 0, 1/16, 1/4 and 1.
No full `2^32` fermion matrix or infinite oscillator matrix is
enumerated. The certificate is a finite rational majorant for the
complete interacting operator, rather than a free-field numerical
trajectory.

For a complete Hadamard detector of V, its two effects are
`(2I+V+V^dagger)/4` and their complement. Inserting U before the
autonomous wait changes the detector law by total variation at most
`||[V(tau),U]||/2`, for every code state. This proves an influence
upper bound, not a positive signal or a measurement of a propagation
speed.
The corresponding Kraus maps are `(I+V)/2` and `(I-V)/2`; their
adjoint products sum exactly to I. A final common source factor
multiplies both maps on the left and preserves completeness. The outer
Gauss-orthogonal unitary extension is identity, so its plus outcome is
deterministic and its full source vector is retained. No undefined
readout action or omitted constraint sector is needed.

## Transfer to the actual full source and all preparation records

Every physical E/B field test, current and density is block diagonal in
the global native code P. This follows from the source's verified
`[P,rho]=0`, photon factorization and globally guarded currents. The
autonomous completed Hamiltonian reduces P. Thus every code history
intertwines exactly with the corresponding restriction of the full
source. Its complement is retained with its actual dynamics and
actual field instruments; it is not replaced by a fictitious local
microscopic net.

The ground-refinement package prepares the actual physical ground with
four independent native trials, three complete energy-agreement filters,
and all `M^3` pointer records, where `M=315179904752653412`.
Its all-pass main state has ground fidelity F and hence code weight
at least F. Because all subsequent effects are block diagonal in P,
code/complement coherences do not contribute to any recorded
probability. For any such history the law is the convex combination
of its normalized code and complement laws. Consequently a code-state
influence bound b transfers to the actual successful preparation as
`b+(1-F)`, without requiring a square-root state-error subtraction.
Here `1-F<5.583e-14`. This support argument concerns the chosen
physical observable algebra; it does not prove full-native operator
locality.

Choose the finite photon cap `K=2^64` within the admitted same-source
family and use its certified finite energy-clock preparation. On filter
failure, both the intervention and detector
slots execute declared no-op records, retaining the failure flag and
all source/partner states. On all-pass, compare two settings: magnetic
kick U versus identity, an autonomous wait `tau=1/64`, and the complete
electric V detector. All original-attempt records consist of every
preparation flag and the electric detector bit. No outcome is removed.
The intervention choice is a classical program input. The two laws are
conditional on that specified input and use identical labels for every
random preparation/detector outcome. A retained input bit may be
appended after canonically relabelling its fixed setting; without that
relabeling, two opposite deterministic input bits trivially have total
variation one. No selection of measured successful data is involved
in the original-attempt law comparison.

The exact finite kick is a duration-one driven matrix service with
Hermitian generator `i Log(U_controlled)` and control
`i Log(U_controlled)-H_copies`; the identity setting uses the same
duration and explicit cancellation control. All four physical sources
and references are retained. This specified control cancels the
background evolution during that slot; it is not an implicit pause
or an autonomous local implementation. The final detector service
has its own common free factor after the bit is recorded. Every
branch has source elapsed time `6+1/(64m)`: one trial slot, three filter
slots, one kick/noop slot, the autonomous wait, and one readout slot.
All finite controls have the parent's raw cap-dependent norm/work
bounds. There are six duration-one pulses, each with control norm at
most `4B_raw(K)+pi`; their gross switching work is at most
`12[4B_raw(K)+pi]`. A retained intervention-setting bit and Bob's bit
add `2ln2` nats to the parent's full four-source/filter capacity.
Coefficient workspace and efficient implementation are separate from
these finite register counts. The code, ground, mass and autonomous
Hamiltonian are unchanged.

For the cutoff comparison, embed the actual finite prepared density
matrix itself, with every flag and partner, into the unbounded
mathematical carrier. The initial states then agree exactly. The
comparison does not repeat the preparation at infinity or substitute
an infinite-source ground. The code commutator estimate holds for
every code state, so the finite source's proved code weight F suffices.
The finite source's energy/moment estimate also applies to its embedded
initial state. All subsequent moment estimates concern this infinite
comparison branch. Put `n=K+1`, `L=floor(K/2)`,
`r=g sqrt(v)` and let `M_N` bound `||sqrt(N)psi||`. The inherited
complete filter tree has an all-branch bound
`M_N^2<0.284617466`. The low-sector displacement lemma and a photon
tail split give

`d_W(M_N)<=n^-2+2M_N/sqrt(L+1)`.

Under an autonomous comparison wait,
`||sqrt(N)psi(tau)||<=M_N+48r tau`. Unitary Duhamel and the same split
therefore give

`d_U(M_N,tau)<=48tau/n^2`
`             +[96tau M_N+2304r tau^2]/sqrt(L+1)`.     (7)

The magnetic kick increases M_N by at most `sqrt(2)/64`. Both actual
field radii are below `1/32`, within the inherited low-sector bound.
Four source copies continue evolving during the autonomous wait, so
all four contributions to (7) are included, including the retained
partners. No preparation discrepancy is needed: `d_p=0` for this
same-state embedding. The two complete protocol discrepancies obey

`d_0<=d_p+4d_U(M_N,tau)+d_W(M_N+48r tau)`,
`d_1<=d_p+d_W(M_N)+3d_U(M_N,tau)`
`       +d_U(M_N+sqrt(2)/64,tau)`
`       +d_W(M_N+sqrt(2)/64+48r tau)`.

The intervention is coherently controlled by retained success flags;
the unconditioned moment bounds apply before any branch conditioning.
The final common source evolution cannot change the stored bit law.
Combining both comparisons, the code influence (6), and the complete
complement allowance proves the actual finite-source, original-attempt
bound

`TV(law_kick,law_identity)`
` <=d_0+d_1+(1-F)+||[V(1/64),U]||/2`
` <1.810490e-8`.                                      (8)

The sum of the two post-preparation cutoff allowances is less than
`3.199429e-9`. This bound uses the actually prepared finite state and
its retained failures, with no approximation to a mathematical ground
as an initial-state premise.

This compares the same complete preparation records and final bit in
the two intervention settings. It is not a bound on the distance
between terminal partner quantum states after an additional uncounted
evolution. Conditional detector laws require their explicit success
denominators. The finite photon cap is global and its full-carrier
linear commutator has the parent's boundary counterexample; (8) does
not replace that counterexample with an exact finite canonical algebra.

The construction supplies a fixed-cube local physical field/current
algebra, an interacting commutator majorant, and a complete prepared
record consequence. It does not identify the global Foldy-Wouthuysen
map with a microscopic local map, assert locality of every native
complement observable, or prove a spatial/ultraviolet continuum or
relativistic cone.
