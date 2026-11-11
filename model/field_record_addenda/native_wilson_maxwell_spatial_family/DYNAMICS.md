# Positive-time local records at fixed lattice spacing

This theorem uses the exact source family, complete code-bank
preparation and photon-cap schedule in [PROOF.md](PROOF.md). Fix a
finite observer window, finitely many electric/magnetic Hadamard
instruments in that window, and a bounded total autonomous time in
the same dimensionless clock `tau=m t`. All detector outcomes and the
local success/failure flag are retained. The intervening physical
Hamiltonian is the full coupled source at `g=1/64`, not a free-field
replacement. The initial accepted state is the source's nonstationary
free Wilson/code vacuum. No larger-volume interacting ground is used.
The observer remains a diverging distance from all six open-box
boundary faces, as specified in the equal-time theorem.

The resulting local record laws have a volume limit, with the
vanishing failure and finite photon-cutoff allowances included. This
is a fixed-spacing local-record limit. It neither constructs a
Lorentz/ultraviolet limit nor identifies an observed particle field.

## Truncate bounded hopping, keeping the full field baseline

On the admitted code of an unbounded-oscillator comparison box write,
in units of m,

`H_n=H_0,n+sum_e V_e,n`, `V_e=F_e+F_e^dagger`,

where `||F_e||<=2`. The baseline contains the full completed photon
quadratic operator, full Coulomb energy, all onsite Wilson matrices,
all scalar zero points, and every reference cost. Nothing in the
quadratic field or its Gauss boundary representation is removed.

Let S contain the coordinates of the fixed observer's initial
commutator support. Keep hopping only on edges whose tail is within
Manhattan distance R of S:

`H_n,R=H_0,n+sum_(e in Omega_R) V_e,n`.

This is a mathematical bounded-perturbation comparison on the same
finite-box physical Hilbert space. It is not an alternative physical
source, preparation, or clock. Its gauge-parent proof has the same
full quadratic operator and the same longitudinal point fiber as the
source. Setting some hopping coefficients to zero only decreases
the nonnegative commutator majorant.

The q/p node coordinate is its edge's tail and the z coordinate is
its vertex. Every majorant entry spans Manhattan distance at most
two. At weight `mu=ln2`, the row bound is `4 kappa`, where
`kappa=16+8g+4g^2`. If C_A bounds the initial commutator components
of a local unitary A, the current Leibniz estimate gives

`||[V_e,alpha_s^(n,R)(A)]||`
` <=(12+4g) C_A exp(4 kappa |s|) 2^(-dist(tail(e),S))`.

The factor 12 accounts for the two endpoints, whose coordinate
distances can differ by one. There are at most
`18 |S| (r+1)^2` positive-direction edges at distance r from S.
The exact dyadic tail identity is

`sum_(r>R)(r+2)^2 2^-r=2^-R(R^2+8R+18)`.

Bounded-perturbation Duhamel therefore proves, uniformly in n,

`||alpha_tau^n(A)-alpha_tau^(n,R)(A)||<=delta_A(R,tau)`,
`delta_A=18|S|(12+4g) C_A |tau| exp(4 kappa |tau|)`
`        *2^-R(R^2+8R+18)`.                           (1)

No commutator with a removed E-squared or B-squared boundary term
occurs. The removed operator is a finite sum of bounded hoppings.
This is the reason the comparison avoids an unproved unbounded
boundary-Hamiltonian estimate.

A complete Hadamard instrument has a controlled unitary followed by
fixed pointer gates. Replacing its Heisenberg unitary A by A_R changes
the entire instrument isometry in norm by at most `||A-A_R||`.
Telescoping the finitely many instruments gives a complete-history
bound `delta_R=sum_j delta_Aj(R,tau_j)`. It applies to every input
state, including its retained references, and tends to zero with R
independently of box volume. The initial state's volume limit is not
being assumed in this step.

## Transverse-coordinate covariance, including the infrared tail

Peierls words require the covariance of Q itself. Equal-time E/B
square-root covariances alone are insufficient. The actual one-form
Hodge decomposition gives the required additional estimate.

For a fixed oriented edge vector e, its squared product-basis
coefficients are at most `8/N^3`. Its component's own-axis sine
index is nonzero. With `u=8t/N^2`, the same product sums as in the
vertex proof imply

`||exp(-t L_1,n)e||^2<=8 t^(-3/2)`, t>0,             (2)

uniformly in n. Indeed
`sum_(k>=1)exp(-u k^2)<=exp(-u/2)sqrt(pi/(2u))` and the
other two sums are at most `1+sqrt(pi)/(2sqrt(u))`.
After substitution the coefficient is bounded by
`2sqrt(pi)exp(-4z^2)(z+sqrt(pi)/(4sqrt(2)))^2<8`.
For an analytic uniform check, `(z+c)^2<=2z^2+2c^2`,
`z^2 exp(-4z^2)<=1/(4e)`, `c^2=pi/32`, `e>2`, and
`pi<22/7` bound this coefficient by
`2sqrt(22/7)[1/4+(22/7)/16]<8`.

P_T commutes with L_1 and its restriction is the transverse K=C^T C.
Consequently

`Sigma_Q,n=K_T^(-1/2)P_T/2`
` =(2sqrt(pi))^-1 integral_0^infinity t^(-1/2)`
`                                  exp(-tK)P_T dt`.

Using (2) at t/2, the diagonal tail beyond heat time T is less than
`12/T` for every local edge; off-diagonal tails obey the corresponding
positive-form Cauchy-Schwarz bound. There is no volume-dependent
infrared constant in this estimate.
Indeed the integral coefficient is `8sqrt(2/pi)<8sqrt(2)<12`.

The dipole Green theorem in PROOF.md gives
`P_T,n=I-D_n G_n D_n^T -> P_T,infinity` strongly on local edge
vectors. Bounded finite-range powers give strong finite-time
convergence of K_n and its continuous functional calculus. Splitting
the displayed heat integral into a finite interval and the uniform
tail proves convergence of every local Q covariance entry.

Free-time coordinate covariances use
`K_T^(-1/2) cos((t-s)sqrt(K))P_T/2`.
The cosine commutes with the positive tail form and has norm at most
one, so the same tail bound applies. The exact Weyl commutator phases
use `sin((t-s)sqrt(K))/sqrt(K) P_T`, a bounded continuous function
times P_T. Electric and magnetic covariances use the corresponding
regular powers or curls. All these local entries converge uniformly
when the finitely many times range over a compact interval. For
example the derivative of the Q covariance cosine factor is bounded
by one half, and the derivative of the sine quotient by one. Thus
the complete free Weyl phases, not only their real covariances, have
controlled limits.

The construction concerns finite local test forms and their finite
free-time evolutions. No bounded inverse square root on every
infinite-lattice l2 vector is assumed, and no extra photon carrier is
put into a finite program.

## Exact free hopping words and their source-state limits

The density is constant under H_0,n: the onsite spin rotations preserve
local charge, and the Coulomb operator is a function of commuting
densities. The actual freely evolved hopping has the ordered form

`F_e^0(t)=Psi_head^dagger exp(i Omega_head t) B_e^dagger`
`          exp(-i Omega_tail t) Psi_tail`
`          *exp[i g Q_T,e(t)]`
`          *exp{i g^2 t[(D G rho_ref)_e+(D G D^T)_ee/2]}`.    (3)

The final half-phase is essential. It follows directly from
`H_C(rho+D_e^T)-H_C(rho)` and equals the Weyl ordering phase in the
mathematical gauge parent. All raw scalar and total-charge reference
terms are retained in the physical clock; they commute with these
gauge-invariant words. The initial sea and every hopping/E/B word
have total charge zero, so the root subtraction may be removed only
on this stated sector.

Each fixed word in finitely many F_e^0(t), their adjoints and freely
evolved E/B detector Weyls reduces to a finite sum of terms consisting
of a local CAR word, a dipole-density exponential and a photon Weyl
word. Commuting a density exponential past a local fermion only
produces its exact endpoint phase. Its density weights are finite
linear combinations of `g G D^T s` and `g^2 t G D_e^T`. They
converge in l2 uniformly on compact time intervals by the proved
three-dimensional Green estimate.

This convergence also holds with the fermion insertions. If A is a
k-factor local CAR monomial, its norm is at most one and

`||X(delta a) A Omega||`
` <=||A X(delta a)Omega||+||[X(delta a),A]Omega||`
` <=(2+k)||delta a||_2`.                              (4)

The first term uses the centered source Slater variance; the second
uses the exact density/CAR commutator at the finitely many endpoints.
Onsite spin rotations preserve this charge identity. A finite linear
combination has the sum of the corresponding finite coefficients.
Thus density tails can be truncated before evaluating a mixed word.
On fixed finite support, Wick's finite-dimensional formula makes all
remaining expectations polynomials or finite determinants in the
local projector entries and bounded density phases. The source
Wilson projector convergence proves their limit. No infinite
determinant, substituted classical scalar, or unproved tail bound on
an excited insertion state is used.

The photon word expectation is its exact Gaussian characteristic
function times the exact Weyl phase. Its entries converge by the
preceding section. Therefore every fixed free-baseline Dyson word
has a volume limit, uniformly over its finitely many times in a
compact simplex. All operator norms are bounded by the finite
product of hopping norms, independently of volume.

## A complete coupled finite-window record limit

Fix R. There are finitely many retained edges in Omega_R, independent
of n once the box contains that buffer. Their perturbation has norm
at most `4|Omega_R|`. The strong interaction-picture Dyson expansion
around H_0,n has its usual absolute norm bound, uniform in n. For
total autonomous time tau, put `x=4|Omega_R| tau`; a sufficient
two-sided omitted-tail allowance is

`delta_J=2 x^(J+1)/(J+1)!/[1-x/(J+2)]`, `J+2>x`.      (5)

All retained terms are the finite free words just controlled. The
bounded hopping strong integrals and their scalar norm majorants
justify dominated convergence; norm continuity of a Weyl path in
time is not assumed. Complete detector isometries may be inserted
between the intervals. Their joint amplitudes have the same total
Dyson weight, and their probabilities are finite combinations of
these retained word expectations. Consequently the complete local
record laws for H_n,R converge as n grows.

Let f_n be the actual code-assembly failure bound, d_n the literal
photon-cap history error, and epsilon_free(n,m;R,J) the finite-word
convergence bound obtained above. The original finite-program local
laws in two boxes satisfy the sufficient estimate

`TV(law_n,law_m)<=f_n+f_m+d_n+d_m`
`                 +2delta_R+2delta_J`
`                 +epsilon_free(n,m;R,J)`.           (6)

All errors concern the same fixed observer program, with every local
detector outcome and its success/failure flag. The remote fine code
flags and all unused banks remain retained. For any positive target
error choose R by (1), then J by (5), then box sizes by the explicit
Green, covariance and Wilson-projector limits. The prepared failure
and cap schedules make f_n,d_n tend to zero. Equation (6) proves a
Cauchy sequence on the fixed finite local record alphabet. Its limit
has nonnegative probabilities summing to one.

For r readout times the finite physical clock is three preparation
slots, r declared cancelling readout slots, and the stated total
autonomous wait `tau/m`, on every success and failure branch. A
same-time group of tests can share one readout slot. Every bank
continues under its actual source clock during the autonomous waits.
The unused-bank unitaries do not change the active local marginal;
their states are not erased. Exact matrix controls, workspace and
raw work costs retain the source-family declarations.

The limit describes complete interacting local records of this
nonstationary preparation at fixed lattice spacing. It does not
transfer the cube's interacting ground, require an infinite physical
universe, or correct the original Wilson squared-dispersion mismatch.
There is no assertion that different finite-box sharp Gauss fibers
are literal nested Hilbert subspaces. The proof compares probabilities
through bounded-hopping truncations and explicitly convergent source
word expectations.
