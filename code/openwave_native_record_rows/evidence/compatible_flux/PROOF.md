# Compatible native record sums: a conditional Gaussian theorem

This is a written construction proposal, with no new source execution or
row award. Its native ternary instruments, record-compatibility filter,
accepted and rejected histories, and normalization must all be declared
as a new N+ apparatus. The theorem does not say that arbitrary conditioning
is selected uniquely by the endpoint law. It proves what the specified
record apparatus would measure, including the cost of its conditioning.

The finite Hodge/Maxwell target identities are already present in
`reverse-engineering-reality/lean/Screen/WhitneyMaxwellDynamics.lean`:
`curl_stiffness_gradient_null`, `scalar_potential_action_expansion` and
`action_derivative_iff_recurrence`. That module supplies its mass and
stiffness pairings. The construction below proposes a native statistical
source for particular pairings; it does not relabel those supplied forms
as existing native evidence. `SeamMaxwellContinuum.lean` similarly provides
finite symbol bounds with its stated continuum and clock premises.

## 1. The finite native record family

Fix m record coordinates and an integer matrix A of rank r. Write
L=ker(A:R^m -> R^k), d=m-r, and Lambda=L intersect Z^m. It is a full
rank-d lattice in L: integer Gaussian elimination gives a rational basis,
and multiplying denominators gives d independent integer vectors. The
lattice of all integer solutions, obtainable by Smith normal form, is
used; it is not replaced by an arbitrary finite-index sublattice.

For coordinate i, an independent native gate accepts with probability
p_i in [p_min,1/2], where p_min>0 is fixed on this finite member. A further
independent fair native record sets the sign. The resulting record is

    P(X_i=+1)=P(X_i=-1)=p_i/2,  P(X_i=0)=1-p_i.

Different coordinates and repetitions use independent prepared probes and
fair records in the ideal comparison. Every bit, probe, failure, hold and
rejected word is retained. Finite persistent-probe reset implementations
are covered by section5, not assumed exactly independent without proof.

A concrete fair record needs no hidden coin-only oracle: prepare one unit
on a disclosed two-port pair and repeatedly call the pinned native endpoint
rule on that pair. The public bit |current| is0 or1 with conditional
probability1/2 for either possible current location. Chain-rule multiplication
therefore proves that these movement bits are IID fair, although the unit
persists and its location is reconstructed from all previous currents.
Three movement-bit gates and a fourth movement sign bit give p=1/8.
Changing the endpoint event to an identity hold makes every movement bit0;
if the gate requires three1 bits, the entire flux variance vanishes. This
control changes the native event law on the observable's own object. A
coin-only decoder on empty pairs would not detect that particular mutation.

Let Z_i=sum_{j=1}^n X_ij. An observer accepts a complete word exactly when
A Z=0, and reads F_n=sqrt(8/n) Z. This is an integer test on public records.
No Gaussian density, energy, field action or Hamiltonian is executed to
produce the word. A is supplied compatibility apparatus and must be named.

## 2. Conditional full-moment limit and acceptance cost

Set Sigma=diag(p_i) and let covol(Lambda) mean Euclidean d-volume of a
fundamental lattice cell in L. For d=0 use the counting convention that
the integral at the origin and the zero-dimensional covolume are one.
Define

    g_Sigma(x)=(2pi)^(-m/2) det(Sigma)^(-1/2)
               exp[-x^T Sigma^(-1) x/2].

Then

    P(AZ=0) ~ n^(-r/2) / covol(Lambda) * integral_L g_Sigma(x) dx.  (1)

Conditioned on AZ=0, F_n converges in every joint polynomial moment to
the centered Gaussian measure on L whose density, with respect to
d-dimensional Lebesgue measure, is proportional to

    exp[-F^T Sigma^(-1) F/16].                              (2)

All mixed linear exponential moments also converge. In particular this
is a full Gaussian hierarchy, not an inference from the covariance alone.
The limiting covariance is

    C=8[ Sigma-Sigma A^T (A Sigma A^T)^+ A Sigma ].           (3)

Dependent rows of A do not change (1)-(3); the exponent in the cost is its
rank r. Each accepted word requires asymptotically a positive constant
times n^(r/2) complete attempts on this fixed finite member. Multiplying
by the actual n*m coordinate-instrument cost is mandatory. A resource
claim cannot omit rejected words or replace their cost by the accepted
word's length.

### Proof of the local estimate used in (1)

For one coordinate, write q_n(z)=P(Z_i=z). Its characteristic function is
f_p(theta)=1-p+p cos(theta), which lies in[0,1] for p<=1/2. Fourier inversion
and 1-cos(theta)>=2theta²/pi² on[-pi,pi] give

    q_n(z)=(1/2pi) integral_{-pi}^pi f_p(theta)^n exp(-iztheta)dtheta,
    f_p(theta)^n <= exp[-2p_min n theta²/pi²].

For |theta|<=1, Taylor's bound and the power-series remainder for log give

    |log f_p(theta)+p theta²/2| <= theta^4/16.               (4)

Indeed y=p(1-cos theta)<=1/4; the logarithm remainder is at most
y²/[2(1-y)], while |1-cos theta-theta²/2|<=theta^4/24.
Their sum is at most (p/24+p²/6)theta^4<=theta^4/16.

On |theta|<=theta0=min(1,sqrt(4p_min)), (4) leaves at least
p_min*theta²/4 in the Gaussian decay exponent. After u=sqrt(n)theta,
the L1 Fourier error against exp(-p u²/2) is bounded by a constant times
n^-1 integral_R u^4 exp(-p_min u²/4)du, plus an exponentially small tail.
The complementary theta region is controlled by the preceding global
bound. Thus, uniformly for all integers z and p in[p_min,1/2],

    sqrt(n) q_n(z) = (2pi p)^(-1/2) exp[-z²/(2np)] + O(n^-1),  (5)

with computable constants from the displayed integrals. The complex phase
has modulus one, so this estimate is uniform in z without a central-window
assumption. Taking the finite product gives the local multivariate estimate
on every compact set of z/sqrt(n).

### A uniform local tail bound, needed after conditioning

A Chernoff bound alone lacks the n^-m/2 prefactor needed here. Use exponential
tilting with t=z/n for |z|<=n. In the tilted three-point law all three atoms
are at least epsilon=p_min/(2e²), since |t|<=1 and its moment-generating
function M(t)<=e. For any such law, the squared characteristic-function
modulus is at most

    1-4epsilon²(1-cos theta),

by keeping just the zero-to-plus and zero-to-minus cross terms. Fourier
inversion therefore bounds every n-fold tilted atom by C_epsilon/sqrt(n).
Also M(t)<=cosh(t)<=exp(t²/2), because the original law is centered and
log cosh(t)<=t²/2. Undoing the tilt yields

    q_n(z) <= C_epsilon n^-1/2 exp[-z²/(2n)]                (6)

uniformly on its support; outside[-n,n] it is zero. This proof keeps both
the Gaussian tail and the local prefactor.

Choose a bounded fundamental parallelepiped of Lambda. The scaled lattice
Lambda/sqrt(n) has cell volume covol(Lambda)n^-d/2 and a cell diameter
tending to zero. Applying (5) on compact sets, the lattice Riemann sums
converge to the integral in (1). Product estimate (6) bounds their tails
by an integrable Gaussian on L, uniformly for all sufficiently large n.
The same argument with any polynomial weight or exp(ell dot x) gives
uniformly vanishing tails for that weight. This proves (1), convergence of
all moments and exponential moments, and the normalized density (2).
Completing the square on L gives (3). This also supplies a finite,
computable convergence procedure: choose a Gaussian tail cutoff, bound
the compact Riemann-sum error using a lattice basis and bounded derivatives,
then use the explicit Fourier error above. No noncomputable convergence
rate needs to be used by a diagonal construction.

## 3. Affine compatible preparations

One may instead condition on A Z=A z_n, where z_n is an integer vector
and z_n/sqrt(n) tends to a fixed x0. The accepted coset is z_n+Lambda.
The same local estimates and translated Riemann sums give the Gaussian
on sqrt(8)(x0+L), with the same quadratic density and covariance on its
tangent space. Its mean is the unique minimum of F^T Sigma^-1 F on that
affine space. The normalization still scales as n^-r/2 with a strictly
positive coefficient. Taking z_n by rounding sqrt(n)x0 is one concrete
choice, but x0 and the physical boundary interpretation must be disclosed.
This theorem does not turn a supplied boundary amplitude into a derived
constant, nor identify the conditioning apparatus with a charged particle.

## 4. Maxwell-shaped finite Gaussian target

Let D be the integer edge-to-plaquette coboundary of a specified finite
cubical complex, with specified boundary conditions. Choose A so that
ker_R A=im_R D. The accepted integer lattice is exactly im_R D intersect
Z^m. On a contractible box, Bianchi compatibility supplies this subspace;
for other topology, harmonic sectors and any extra period constraints
must be analyzed explicitly. No reflection-positivity claim follows merely
from imposing arbitrary global zero-mode constraints.

If 8p_f=1/w_f, then (2) is the Gaussian flux measure

    dmu(F) proportional exp[-(1/2) sum_f w_f F_f²] dF,
    F in im D.                                             (7)

Writing W=diag(w_f), its covariance is exactly

    C_F=D (D^T W D)^+ D^T.                                  (8)

This identity is obtained by using D on the orthogonal complement of its
kernel as coordinates on im D. The Jacobian is constant in F and cancels
in the normalized law. Gauge zero modes are quotient coordinates, not
unobserved extra positive photon polarizations.

For the static metric ds²=-u^-2 dt²+u² dx², the continuum target weights are
w_E=u² and w_B=u^-2. The plus-hit controller supplies b=u/2; the separately
weighted native first-hit response supplies q=1/u. Two independent b flags
and one fair flag give p_B=u²/8. Two independent q flags and three fair
flags give p_E=u^-2/8. The sign uses another independent fair native bit.
Since1<=u<=2, all probabilities lie in[1/32,1/2], exactly the bounds used
above. No inverse-probability estimate or potential table is needed in
this ideal nested instrument. Its actual recursive native costs, caps and
reset implementation are additional obligations before an execution claim.

A discrete electric edge pairing can more precisely use M_xy=u_x u_y and
p_E,xy=q_x q_y/8. Then the same source-derived reciprocal v_x=1/u_x obeys
sum_y M_xy(v_x-v_y)=sum_y(u_y-u_x)=0 in the original harmonic exterior.
This is the finite Dirichlet mean target for an affine boundary preparation
of (7); it is not yet an executed boundary preparation or stress measurement.

## 5. Errors must be tightened before acceptance

Let P_n be the ideal complete-word law and Q_n the fully capped, reset,
native implementation on the same retained record alphabet. Suppose their
total-variation distance is at most delta_n, proved by sequential coupling
of the complete word. Put alpha_n=P_n(AZ=0). If delta_n<alpha_n, the actual
acceptance probability is positive and

    TV(P_n(.|AZ=0), Q_n(.|AZ=0)) <= 2delta_n/alpha_n.         (9)

This follows by adding the numerator-law error and the acceptance-
normalization error. Since alpha_n is of order n^-r/2, delta_n->0 alone
does not suffice. For a degree-s polynomial, |F_n|<=sqrt(8mn) bounds its
supremum by a constant times n^(s/2). A sufficient full-moment schedule is

    delta_n=o(n^(-(r+s)/2)).                                (10)

For instance delta_n<=n^(-(r+s)/2-2) works on this fixed graph. Per-event
reset and cap allowances must be divided by the full number of attempted
coordinate instruments to reach that complete-word bound. Rejected words
multiply the expected resource cost but do not change the conditional-law
bound for one successful word. To certify an entire experiment with K
attempted words, the full archive error is bounded separately by Kdelta_n.

The conditional CLT is first a fixed-complex theorem. A spatial/time
continuum diagonal must recompute r, the acceptance lower bound and all
moment/cap budgets for each finite complex; it cannot treat r as fixed
while the number of plaquettes grows. Source first-hit approximation,
Gaussian-word limit and geometric continuum error are distinct entries.

## 6. Reflection positivity and remaining construction

Equation (7) alone is not yet a complete curved Maxwell reconstruction.
For a time-independent spatial complex and conducting boundary conditions,
the expected proof uses a positive electric pairing M on spatial edges
and K=C^T beta C. On the gauge quotient, integration over the temporal
potential produces a positive Gaussian convolution in the difference of
neighboring spatial potentials. Its Fourier multiplier is a positive
Gaussian. Sandwiching that convolution by the square root of the positive
magnetic multiplication weight gives a positive self-adjoint transfer.
The same oriented plaquette readout then gives the reflection action on F.

That proof must specify the exact finite complex, boundary state, quotient,
all harmonic sectors, normalization and physical field algebra. In particular
an unexamined global torus projection is not justified by this paragraph.
The zero-mode caution is substantive: the supplement of Borsanyi et al.,
[arXiv:1406.4088](https://arxiv.org/pdf/1406.4088), distinguishes a global
zero-mode subtraction that breaks reflection positivity from a per-slice
treatment. Their physical QED calculation is not imported as OPH evidence.

Once those finite conditions are proved, smooth-background continuum and
stress convergence, the shared scalar clock, canonical normalization and
native boundary-state preparation remain to be joined. Target inertial
energy, universal gravitational response and a cosmological term are not
consequences of the conditional Gaussian theorem alone. This file supplies
the complete conditional full-moment step and its error/cost obligations;
it is presently an unexecuted proposal awaiting independent written review.
