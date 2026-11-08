# Exact finite color action under spectral coarse-graining

The complete gauge–Higgs source determines a positive Euclidean transfer
family and an exact action on its six scalar-vacuum flux states. Eliminating
the other 432 states generates additional gauge-invariant operators. Exact
intervals resolve those operators and the scale dependence of the projected
electric and magnetic coefficients. The source Hamiltonian, Gauss space,
scalar cutoff and clock are those of the [complete color-loop
readout](../color_loop_readout/PROOF.md). No particle data enter.

## 1. The literal link and the Wilson-coordinate test

The link register is the source space

    L = C|0> plus direct sum_r (V_r* tensor V_r), dim L = 60.

Its entries are U_ab=|r,a,b><0|. For a block of dimension d, regard U as
an operator on C^d tensor L. Direct multiplication gives

    U² = 0,       U*U = d I_d tensor |0><0|.                    (1)

Thus V=U/sqrt(d) is a partial isometry with orthogonal initial and final
spaces. In particular

    ||V*V-I|| = 1.                                            (2)

The same identity holds after tensoring any spectator or increasing the
number of untouched links. On every unit vector psi,
|<psi,V psi>|<=1/2: write its orthogonal components in the initial space,
final space and their complement and use ab<=(a²+b²)/2. Consequently
||(V-I)psi||>=1/2. No statewise limit of this normalized local operator is
the identity. This is a direct algebraic test of the specified raising
coordinate; it is not a prohibition on collective gauge coordinates or on
continuum limits of finite-dimensional lattice systems.

An ordinary Wilson link matrix is unitary on its full group-coordinate
space. Equations (1)–(2) exclude substituting the source raising matrices
for that coordinate without a separate reconstruction. A larger tensor
product alone does not change the tested local identities. The native
Casimir law and Gauss covariance therefore do not determine a Wilson bare
coupling by taking a ratio of the existing coefficients. The source's
electric color coefficient 2e remains a defined finite action coefficient.

## 2. Positive transfer on the same source

Let H be the exact 438-dimensional Hermitian Hamiltonian at the identity
family controls. The source provides an outward rational bound B>=||H||.
Set H_+=H+B I and T(t)=exp(-tH_+) for t>=0. Functional calculus gives

    0<T(t)<=I,       T(t+s)=T(t)T(s).                         (3)

This uses the existing generator, with a scalar reference shift, rather
than another evolution law. Every T(t) commutes with the Gauss action.
For vectors v_i and nonnegative times t_i, its reflected kernel satisfies

    sum_ij <v_i,T(t_i+t_j)v_j>
      = ||sum_i T(t_i)v_i||² >= 0.                           (4)

More generally any finite ordered source-operator histories F_i have the
reflected Gram <F_i chi_0,F_j chi_0> on the original empty source input.
This proves finite reflection positivity with that boundary state. It does
not turn chi_0 into a stationary ground state or prove the spacetime
hypotheses of a relativistic reconstruction theorem.

The contraction has a complete finite measurement realization at the
operator level: K_ok=T(t), K_fail=sqrt(I-T(t)²). Both outcomes, their
output states and untouched references are retained. This is a declared
source-derived instrument, without an asserted native repair compiler.
Its successful branch is never divided by its probability in (4).
Imaginary time is the analytic argument of this transfer family; it is
not a second physical clock calibration.

## 3. An exact action retaining the whole eliminated sector

Let P be the scalar-vacuum projector and Q=I-P in Q*=0. Their dimensions
are six and 432. Define A=PHP, C=QHP and D=QHQ in the actual orthonormal
source basis. For a real spectral resolution mu>B, the positive operator
mu+H has the exact Schur action

    K(mu) = mu I_P + A - Sigma(mu),
    Sigma(mu) = C*(mu I_Q+D)^(-1)C.                           (5)

For every coarse vector phi,

    <phi,K(mu)phi>
      = min_xi <phi+xi,(mu+H)(phi+xi)>,   xi in QH,
    xi_min = -(mu+D)^(-1)C phi,
    P(mu+H)^(-1)P = K(mu)^(-1).                              (6)

Completing the square proves all three statements. This integrates the
full complement; it does not delete its probability or assert that P is
invariant under H. In particular

    (mu-B)I <= K(mu) <= (mu+B)I,
    K'(mu)=I+C*(mu+D)^(-2)C >= I,
    K''(mu)=-2C*(mu+D)^(-3)C <= 0.                          (7)

All higher derivatives of Sigma have alternating signs. These positivity
and monotonicity properties supply a controlled finite family of effective
actions and a resolvent observable dictionary. The transfer representation
is (mu+H)^(-1)=integral_0^infinity exp(-mu t)exp(-tH)dt; convergence follows
from mu>B. Neither (5) nor this transform identifies mu with a laboratory
energy or a chosen continuum renormalization scale.

Under a common clock change H->lambda H and mu->lambda mu, with lambda>0,
K->lambda K. Thus K(mu)/mu and dimensionless ratios among its coefficients
have no dependence on the overall clock unit. No physical reference scale
is generated by this covariance.

There is also a complete finite energy-record dictionary for (6). Write
H=sum_j epsilon_j Pi_j using all its spectral projections, including
degeneracies. For a normalized coarse input phi, score every energy outcome
epsilon_j by 1/(mu+epsilon_j). Its mean is
<phi,P(mu+H)^(-1)P phi>. In the real source basis,
the six basis inputs and the thirty normalized sums and differences
(e_i plus or minus e_j)/sqrt(2) reconstruct every entry of this resolvent:
the difference of the two scored means is twice its off-diagonal entry.
Every energy outcome is retained. General complex sources also require
the corresponding imaginary superpositions.

The scaled score (mu-B)/(mu+epsilon_j) lies in [0,1], so complete preparation,
instrument and sampling errors can be propagated as bounded-score errors.
If each reconstructed resolvent entry has error at most delta, the operator
error is at most eta=6delta. Its exact eigenvalue floor is g_*=1/(mu+B).
For eta<g_*, inversion contributes at most

    ||G_est^(-1)-K(mu)|| <= eta/[g_*(g_*-eta)].              (7a)

This gives an observable definition of the whole effective action on the
same finite Hamiltonian. The energy projectors and coherent coarse inputs
are specified mathematical instruments; no native Higgs-source pulse
compiler is claimed. Observers with access to the source matrices can
also calculate this dictionary directly, as the exact receipt does.

## 4. The source generates more than an electric/magnetic ratio

Use the normalized scalar-vacuum order (empty, u-conjugate, L, Q,
d-conjugate, e-conjugate). Their raw source squared norms are
(1,9,4,36,9,1), so normalization here is exact and rational. Write

    A=e E-m M,
    E=diag(0,8,9/4,25/4,6,6),
    M_01=M_10=3, M_03=M_30=6, M_04=M_40=3,                 (8)

with all other entries of M zero. The full source, including all scalar
and matter states, gives

    C*C = 12g² I,
    C*DC = 12g²(A+sI),
    C*D²C = 12g²(A+sI)² + g²h² R.                           (9)

The only nonzero rows and columns of R are (empty,L,e-conjugate):

                   [36 12 12]
    R_restricted = [12  6  0].                              (10)
                   [12  0 12]

It is positive, with nonzero eigenvalues 27 plus or minus 3sqrt(33).
In source operators C=g QSP; QSQ=0 at the specified scalar cap. The
electric and loop pieces intertwine S with the same scalar-vacuum action,
whereas C*T C=0. Expanding D C and its Gram gives (9), with
R=P S T² S P. The exact replay evaluates these products in all 438 source
coordinates at several rational controls. It never evolves a truncated
six-state Hamiltonian in place of the full source.

For mu>||D||, (5) therefore has the expansion

    K(mu)=mu I+A-12g²I/mu+12g²(A+sI)/mu²
          -[12g²(A+sI)²+g²h²R]/mu³+O(mu^-4).              (11)

The empty-to-L and empty-to-e-conjugate entries of (A+sI)² vanish.
Both corresponding entries of R equal 12. For gh nonzero, the exact
effective action therefore contains transitions absent from every
operator in span{I,E,M}. A running value of a single gauge coefficient
cannot represent this source coarse-graining. This conclusion is about
the complete finite source; it does not exclude removal of extra operators
in a separately controlled continuum limit.

The program retains all six-by-six entries rather than forcing such a
closure. It also reports the Hilbert–Schmidt projection onto I,E,-M as a
defined finite scheme. Its Gram constants are

    tr E=57/2, tr E²=1441/8, 6 tr E²-(tr E)²=537/2,
    tr M²=108.                                              (12)

The residual is retained. Projected electric and magnetic coefficients
are not called alpha_s, since this finite scheme has no proved QCD
observable matching map.

The positive projected electric/magnetic ratio has disjoint certified
intervals at mu=64,128,256, near 1.000091738891, 1.000013363087 and
1.000001810753 respectively. Thus this specified finite scheme has a
resolved scale dependence even though the source defaults e=m=1 are
fixed. It is accompanied by the nonzero residual in (11).

## Exact numerical errors and source boundary

Let B_D be an outward rational row bound for the normalized full Q block.
For M_k=C*D^kC, the Neumann approximation through degree n has error

    ||Sigma(mu)-sum_(k=0)^n (-1)^k M_k/mu^(k+1)||
      <= 12g² B_D^(n+1)/[mu^(n+1)(mu-B_D)].                 (13)

This follows by summing the operator geometric tail, using the exact
identity ||C||²=12g². It bounds every matrix entry and the entire operator,
including inputs entangled with a reference. Finite weighted inner
products compute every M_k exactly. The retained values mu=64,128,256,
n=24 give strictly negative intervals for both induced entries in (11).
No numerical diagonalization or decimal residual supplies this proof.
Independent full-matrix solves are tests only.

For an actual family of this source's admitted coefficients, with
mu_n>B_D,n and a positive comparison energy Z_n, the same proof gives

    ||Sigma_n/Z_n|| <= 12g_n²/[Z_n(mu_n-B_D,n)],
    ||[Sigma_n-M_0,n/mu_n+M_1,n/mu_n²]/Z_n||
       <= 12g_n² B_D,n²/[Z_n mu_n²(mu_n-B_D,n)].             (14)

These are sufficient uniform error criteria for dropping the full
self-energy, or only its terms beyond the first two moments. They are
explicit hypotheses on the given family, rather than a statement that
the extra operators are irrelevant. Under a common clock rescaling all
coefficients, mu, B_D and Z multiply by the same positive factor. Both
bounds in (14), and the normalized exact action, are unchanged. Increasing
a capacity label or changing only the clock unit therefore does not
suppress the induced dimensionless action. Taking mu to infinity at fixed
coefficients does suppress it, but that is a high-resolution resolvent
limit on the fixed finite graph, not a spatial continuum construction.

Equations (3)–(13) provide a finite positive transfer, its exact resolvent
action and controlled scale-dependent source coefficients. A physical
strong coupling requires a spatial refining family, a common state and
clock, a surviving color action and an observable matching map to a
specified physical scheme. The raising-link test and the induced operators
specify which direct substitutions fail. The construction supplies neither
the scaling of these extra operators in a continuum limit nor a physical
alpha_s value. Its physical coupling and particle-mass fields are null.
