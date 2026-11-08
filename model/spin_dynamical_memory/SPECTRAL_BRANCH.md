# Source spectral windows from bounded resolvent preparations

## 1. Actual response hierarchy and covariance

Use a chosen source ground Omega and its positive generator K, on one
specified state and clock. Let A_a be admitted bounded creation operators,
||A_a||<=1, and z>0. Put S=z(K+z)^(-1) and psi_a=A_a Omega. For r>=1 set

    H_r=span{S^j psi_a: a=1,...,m, 0<=j<r},
    Pi_r=projection onto H_r,  T_r=Pi_r S|_(H_r).            (1)

Zero or dependent columns can be omitted after a rank certificate. These
are mathematical response spaces, not presumed particles or executable
projection detectors. H_infinity=closure(union_r H_r) reduces S and K:
S preserves the dense generating union, and S is bounded selfadjoint.
For every fixed vector in H_infinity, Pi_r converges strongly to identity.
The finite reducing defect of `RESOLVENT_CLOSURE.md` need not converge to
zero in operator norm as the response dimension grows.

If spatial translations commute with the actual clock, include their
orbits of the seed creators. The closed orbit-generated response then
also reduces those translations. For a finite translation group, each
H_r is translation invariant; for the continuum, an increasing dense
translation orbit gives the covariance of the closed response, with no
asserted covariance of an arbitrarily truncated orbit. On the joint
positive-field spin family use the even microscopic translations common
to all four fields. No equality to another state or clock is used.

For any two seed vectors and every integer 0<=j<=2r-1,

    <psi_a,S^j psi_b>=<psi_a,T_r^j psi_b>.                   (2)

For powers through r-1, T_r^j psi=S^j psi. Split a power at most2r-2 into
two such powers and use selfadjointness. At2r-1 use
<S^(r-1)psi_a,S S^(r-1)psi_b> and insert Pi_r between the two retained
vectors. This proves (2) for a block hierarchy, including rank deficiency.

For a seed psi of original norm squared w<=1 and a real L-Lipschitz
function f on [0,1], its degree-N Bernstein polynomial B_N f has uniform
error at most L/(2sqrt(N)): it is E[f(X/N)] for a binomial variable X of
mean Nu and variance Nu(1-u). With N=2r-1, (2) implies

    |<psi,f(S)psi>-<psi,f(T_r)psi>| <= w L/sqrt(2r-1).      (3)

This is a uniform bounded-score estimate independent of source volume and
of an assumed dispersion. The compressed spectral measure has positive
weights summing to the original w. Its finite Ritz values are quadrature
nodes, not certificates of physical spectral atoms. On the cyclic union,
continuous bounded scores converge; discontinuous windows need the
boundary estimate below.

## 2. Gram stability and complete preparation probabilities

For the d retained source columns phi_j, let G_ij=<phi_i,phi_j>.
A strictly positive floor G>=mu I must be certified before inversion.
An exact positive rational matrix has the conservative floor

    mu=det(G)/(trace G)^(d-1)>0.                            (4)

Every eigenvalue is at most trace G, which proves (4). For measured entries
with individual errors at most eta, the operator error is at most d eta.
Thus a certified floor mu_hat of the measured Hermitian matrix yields
mu=mu_hat-d eta whenever this is positive. A source-specific uniform
floor as r or volume grows is not inferred.

An admitted coherent combination sum_j c_j phi_j is prepared by a linear
combination of their complete source creation instruments, normalized by
alpha=sum_j|c_j|. Its original success probability is

    c*G c/alpha^2 >= mu/d.                                (5)

The normalization and all failed branches are retained. Equation(5)
constructs finite response vectors; it does not implement the abstract
Gram projection Pi_r. Interference of admitted coherent creations supplies
cross Gram entries, with their actual normalization errors. A small Gram
floor is an experimental cost and cannot be removed by conditioning.

## 3. Bounded moments from local source histories

The required moment m_j=<psi,S^j psi> is the original success probability
of S^(j/2) A Omega. These vectors have an explicit same-source preparation.
For alpha>0 and E>=0, use

    (z/(E+z))^alpha = integral_0^infinity gamma_(alpha,z)(u) exp(-u E) du,
    gamma_(alpha,z)(u)=z^alpha u^(alpha-1)exp(-zu)/Gamma(alpha).

Closing the elementary Fourier contour in the appropriate half-plane
gives the Cauchy transform exp(-u|E|) from density u/[pi(u^2+t^2)]. Consequently the even nonnegative normalized weight

    w_(alpha,z)(t)=integral_0^infinity gamma_(alpha,z)(u)
                                   u/[pi(u^2+t^2)] du

satisfies

    S^alpha A Omega=integral_R w_(alpha,z)(t) tau_t(A)Omega dt,
    integral_(|t|>T) w_(alpha,z)(t)dt <= 2alpha/(pi z T).    (6)

Here tau_t(A)=exp(itH)A exp(-itH). Groundness gives
exp(itK)A Omega=tau_t(A)Omega exactly, so the scalar ground phase is not
silently changed or omitted. The tail bound follows from the Cauchy tail
bound2u/(pi T) and the Gamma mean alpha/z. Alpha=0 is the identity.

For a local source A, C_A=||[H,A]|| is finite and
||tau_t(A)-tau_s(A)||<=C_A|t-s|. Truncate(6) at T and use cell midpoints
of width at most h with their integrated positive weights. A positive
rational coefficient approximation of l1 error epsilon_c gives a selected
creation-column error at most

    epsilon_A=||A||[2alpha/(pi z T)+epsilon_c]+C_A h/2.     (7)

The missing coefficient weight goes to a rejected branch; it is not
renormalized. Coherent control of the time-translated complete creation
instruments gives a finite dilation. Each selected creation is a
contraction, so its original success score error is at most2epsilon_A.
This selected-column estimate is not an arbitrary-complement channel
bound. Actual native complete-instrument errors, preparation errors and
sampling errors must be added from their own complete-history estimates.

The positive weights are computable definite integrals with the explicit
Gamma and Cauchy tails. Finite rational approximations suffice. For an
actual finite source, the native controlled times and creator dilations
are the admitted services. For growing spin supports, truncate each
source time history using its proved local-history bound. The largest
time is T; the corresponding spatial buffer and positive-field marginal
preparation modulus are separate computable terms. At the joint family
scale, the physical duration is T/(Jp), and C_A is divided by Jp.
These are finite resource statements. They prove no effective convergence
rate to the continuum history state.

A useful positive-field refinement is exact. If the source ground is
unique and the full excitation gap is at least g>0, then

    ||S^N-P_ground|| <= (z/(z+g))^N.                       (8)

The filter I-S^N rejects the ground and retains an excitation of energy
at least g with squared amplitude at least[1-(z/(z+g))^N]^2.
This is a contraction-operator bound. Direct coherent subtraction of
I and S^N has coefficient normalization2, so that LCU implementation
retains the displayed weight divided by4. Unit-normalized acceptance
requires its own compiled dilation; it is not attributed to that LCU.
For the actual spin source g=B/8. On K_B/(Jp), the joint prescription gives
g>=1/(8n^3), which tends to zero. Equation(8) is useful for finite
preparation and readout; it is not a positive continuum mass gap.

## 4. Original interval probabilities and finite resolution

For a closed interval I in [0,1] and delta>0, take ramp functions
f_-<=1_I<=f_+, valued in [0,1], with Lipschitz constants1/delta.
f_- rises inside the interval and vanishes at its endpoints; f_+ equals
one on I and falls to zero outside its delta extension. Then(3) gives

    q_r(f_-)-w/[delta sqrt(2r-1)] <= mu(I)
              <=q_r(f_+)+w/[delta sqrt(2r-1)],             (9)

clipped to [0,w]. Here q_r is the original compressed spectral measure.
The difference between the two exact ramp scores is bounded by the
actual weight of the delta strips around the endpoints. An atom at an
endpoint prevents that boundary error from vanishing. No uniform
convergence of discontinuous CDFs is asserted.

One can also evaluate Bernstein scores directly from measured moments.
For b_k=f(k/N), the power coefficients are

    a_j=binomial(N,j) sum_(k=0)^j (-1)^(j-k)binomial(j,k)b_k.

Entry errors epsilon_j contribute at most sum_j |a_j|epsilon_j, in
addition to the analytic Bernstein error. Thus no unverified Gram
inversion is required for this alternative, and cancellation/statistical
costs are explicit. Neither method conditions away low original weight.
The monotone map E=z(1/u-1) converts finite energy windows to S windows.
An energy tail above M has original weight at most C/M if the chosen
source preparation has original first energy moment at most C. Such a
moment bound must belong to the same preparation and clock.

## 5. Exact source example with joint finite symmetry labels

Use the unchanged two-site source J=1,B=3/8 from the closure control,
K=H_B+9I/8, z=1/8, and

    psi=(I-S)(X_0+Z_0)Omega/2.

The exact site permutation T exchanging the two sites commutes with
this source Hamiltonian and fixes Omega. It is a finite Z2 graph symmetry,
not the continuum spatial momentum of a physical particle. The joint
original success weights are

    (E=9/8,T=+1):81/4000,
    (E=9/8,T=-1):729/4000,
    (E=5/4,T=+1):16/121.                                  (10)

Preparation failure1/2 and filter failure7999/48400 are retained. The
energy Krylov vectors psi,S psi have a positive exact Gram matrix and
span the two energy components; their generalized resolvent eigenvalues
are1/10 and1/11. Adding the translated seed T psi yields the exact
three-dimensional joint response. Its Gram floor, complete polynomial
preparation probabilities and closure are checked from the full source
matrix. The response contains both sharp finite levels and different
finite symmetry labels. It is not a continuum pole or a particle census.

## 6. What a limiting spectral branch would require

Let actual original joint energy-momentum measures mu_n converge weakly
on one compact observation window K, with complete errors tending to zero.
For a continuous proposed dispersion omega(p), suppose shrinking closed
tubes |E-omega(p)|<=delta_n, delta_n->0, retain original weight at least
q>0 after their full errors. For any fixed delta>0 the eventual tubes lie
in the fixed closed delta tube. The portmanteau bound and continuity from
above imply

    mu({(E,p) in K:E=omega(p)}) >= q.                     (11)

This is a sufficient positive spectral-branch certificate, not an outcome
of the finite hierarchy. Weight may instead remain in a continuous
threshold; the finite window tools apply to that case without an atom
inference. Irreducibility, particle scattering and physical labels are
additional questions even if(11) holds.

Given actual energy and spatial-momentum intervals [E_-,E_+] and
[p_-,p_+] and an independently fixed positive speed c, the invariant
E^2-c^2|p|^2 lies in

    [E_-^2-c^2 p_+^2, E_+^2-c^2 p_-^2].                  (12)

Only if the lower endpoint is nonnegative does this give a timelike
squared-mass interval by division by c^4. Calling it a physical mass
also requires a common Lorentz action and the appropriate particle or
threshold interpretation; no speed or mass is selected by fitting(12).
The finite dimer's Z2 label is not substituted for |p| in this test.
