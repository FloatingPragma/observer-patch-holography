# Hydrogenic standing levels from the native first-hit feedback source

Written theorem on the declared N+ `MFR-native-probe-binding` architecture.
This extends the controller proof in `../2026-11-14_codex_probe_binding`;
it does not replace native events by a Coulomb solver. No new numerical run
is asserted here. The native controller uses target/escape letters, never
the function or strength appearing in its derived description below.

## 1. Fixed source and statement

Let a=1/m, delta=a², m a positive integer. Outer external seams have weight
m² and internal seams weight1. The inner probe uses unweighted, port-
preserving cubic spatial routing, target cell0, escape at the boundary of
[-B,B]^3, with B=m³. The outer spatial graph is the periodic cube of side
2B+1, and its Euclidean time period is H=2B+1 as well. Its feedback at coordinate z accepts an internal seam precisely when
the native inner first-hit probe escapes. Relocation uses native attempts
until successful; the finite probe stops almost surely. These are the
exact-preparation/uncapped production terms. Finite cutoffs are charged in
section7, rather than silently replacing this marginal by an exact Markov
law after a history-dependent failed preparation.

The source-derived scalar internal bands are

    lambda = 5-sqrt(5), 6, 5+sqrt(5),

with multiplicities3,5,3. Fix one component of one positive band. Write
L for the ordinary spatial cubic Laplacian, h_B for the target-before-
escape probability, and

    S_a = a^-2 L + lambda(1-h_B),
    H_delta = (S_a-lambda)/a²
            = delta^-2 L - lambda h_B/delta.                 (1)

The Bohr coordinate is y=delta z. Its periodic side tends to infinity:
delta(2B+1) ~2m. For every fixed principal shell n, the corresponding n²
ordered eigenvalues of H_delta converge, counting multiplicity, to

    e_n = -(lambda c_*)²/(4n²),
    c_* = 1/(4pi G(0)),                                    (2)
    G(0) = integral_[-pi,pi]^3 [d³p/(2pi)³] /
                       [4 sum_j sin²(p_j/2)].

No further eigenvalues persist below a fixed negative threshold between
these levels. The total source degeneracy includes the above internal
multiplicity; n² is the spatial orbital shell multiplicity per component.
At finite a the cubic source may split the members of a shell.

The reconstructed standing-wave energy is

    E_a(s)=2 asinh(a sqrt(s)/2)/a.

Writing s_(n,k)(a) for a squared-energy level, k=1,...,n², the conclusions
are

    s_(n,k)(a)=lambda-a²(lambda c_*)²/(4n²)+o(a²),
    E_a(lambda)-E_a(s_(n,k)(a))
       = a² lambda^(3/2)c_*²/(8n²)+o(a²).                   (3)

The subtraction in(3) is essential: E_a(lambda) itself differs from
sqrt(lambda) by -a² lambda^(3/2)/24+O(a^4). The unproved target paragraph
in the earlier controller proof omitted this same-order vacuum lattice
shift when describing energy relative to sqrt(lambda). Equation(3)
corrects that paragraph and measures binding against the same source's
free threshold.

## 2. Native first-hit probability and its asymptotic coefficient

Conditional on a spatial move the six cell directions are equally likely;
native holds and internal moves only change the stopping cost. Therefore
the cell hit probability is the simple cubic one. Strong Markov at the
first visit to0 gives h_infinity(z)=G(z)/G(0). Here G is the inverse of
the cubic Laplacian; its normalization is one sixth of the discrete-time
simple-random-walk Green function.

The simple-cubic Green asymptotic is G(z)=1/(4pi|z|)+O(|z|^-3), and
G(z) <= C/(1+|z|). These follow from the finite-range symmetric-walk Green
theorem, with covariance I/3 and the factor1/6 just specified.
Reference: Lawler–Limic, [Random Walk: A Modern Introduction,
Theorem4.3.1](https://www.math.uchicago.edu/~lawler/srwbook.pdf), and its
first-hit identity in section4.6. The rest of this proof supplies the
singular-potential spectral limit needed here; it is not imported from a
bounded-potential continuum theorem.

Let tau_B be first boundary hit. Splitting paths according to first target
or boundary gives

    0 <= h_infinity(z)-h_B(z)
       = E_z[h_infinity(X_tau_B); tau_B<tau_0] <= C/B.        (4)

Consequently the potential error from the finite native probe in(1) is
uniformly at most C lambda/(delta B)=O(1/m). It vanishes. On every fixed
annulus rho<=|y|<=R, the sampled potential

    V_delta(y)=lambda h_infinity(z)/delta, y=delta z,

converges uniformly to V(y)=beta/|y|, beta=lambda c_*. Indeed the far-tail
error is O(delta²/rho³); moving a cell representative contributes
O(delta/rho²). Globally V_delta<=C_lambda/(|y|+delta).

## 3. The singular target cannot create an extra deep level

Use norm ||u||²_delta=delta³ sum_z |u_z|² and kinetic form
T_delta[u]=delta sum_{z,j}|u_(z+e_j)-u_z|². On a periodic cube of physical
side D>=1 there is a constant independent of delta,D such that

    ||u||_(6,delta)² <= C(T_delta[u]+D^-2||u||²_delta).       (5)

Here and below only finitely many bounded-size m can be omitted from an
asymptotic statement. For completeness, interpolate nodal values
trilinearly on cubes. Norm equivalence for polynomials of degree at most1
in each variable bounds the discrete sixth norm by the interpolant's
sixth norm. Integrating its derivatives bounds its gradient norm by the
edge form. The periodic Sobolev inequality on a unit cube, scaled to side
D, proves(5). On the infinite lattice the D^-2 term is absent. All norm
equivalence constants are finite-dimensional constants, independent of
the mesh.

Counting lattice shells and using V_delta<=C/(|y|+delta), for delta<=rho

    ||1_(|y|<rho) V_delta||_(3/2,delta) <= C rho.            (6)

The origin cell contributes O(delta), not a nonvanishing point defect, to
this norm. Holder and(5) therefore bound the singular-region form by
C rho(T_delta+D^-2||u||²). On the complement it is bounded by
C/rho ||u||². Thus the potentials are uniformly infinitesimally form
bounded. The same bounds apply to the continuum Coulomb potential.
In particular H_delta has a uniform lower bound, and the lattice target
value h(0)=1 cannot produce a hidden level diverging to minus infinity.

## 4. Compact spectral convergence below zero

Here are details sufficient to keep both eigenvalue multiplicities and
the absence of additional negative levels. Embed lattice vectors as
piecewise-constant functions on mesh cubes, using the isometry J_delta.
The periodic cubes expand to R³; extend their functions by zero outside.
Fix kappa>0 and consider the positive Birman–Schwinger operator

    K_delta(kappa)=V_delta^(1/2)
                      (T_delta+kappa)^-1 V_delta^(1/2).     (7)

First restrict V_delta to the annulus rho<|y|<R using its indicator.
Inner, annular and outer supports are disjoint, so the square roots sum
exactly. The annular indicator can be approximated by smooth bounded
cutoffs; the thin-shell L^(3/2) norm tends to zero and(5) controls the
resulting resolvent error uniformly. This also justifies using smooth
cutoffs in the local compactness argument below. The free lattice
resolvents converge strongly under J_delta to(-Delta+kappa)^-1. This can
be seen on compact smooth data by consistency of finite differences and
the resolvent bound1/kappa; density extends it to L². Periodic boundaries
are eventually outside every compact set. A cutoff energy estimate, or
the same resolvent comparison on growing boxes, removes boundary terms.

On a fixed bounded set, the images of the resolvents of bounded inputs
have uniformly bounded discrete gradient and norm. The interpolants used
in(5) have bounded H¹ norm, and their L² difference from the step embedding
is at most C delta sqrt(T_delta). Rellich compactness makes their local
restrictions collectively compact. Strong convergence, self-adjointness
and this compactness imply norm convergence of the doubly localized
resolvents: if norm convergence failed, a sequence of unit inputs and its
compact output subsequence would contradict the strong limit tested
against fixed vectors. Uniform convergence of the annular potential then
gives norm convergence of the annular versions of(7).

The discarded inner part has

    ||V_inner^(1/2)(T_delta+kappa)^-1/2||²
          <= C rho(1+1/(kappa D²)),

by(5)–(6). The discarded outer part has the same squared norm at most
C/(kappa R). The norm of the full potential-resolvent square root is
uniformly bounded by the form estimate. Writing V as the sum of these
three nonnegative pieces, Cauchy–Schwarz bounds cross terms by the product
of the corresponding square-root bounds. Hence removing the two tails
changes(7), after embedding, by a quantity tending to zero as rho->0,
R->infinity, uniformly for small delta. The continuum bounds are identical.
Together with the annular result this proves norm convergence of(7) to

    K(kappa)=(beta/r)^(1/2)(-Delta+kappa)^-1(beta/r)^(1/2),

a compact operator. Replacing h_infinity by h_B changes H_delta by the
operator norm O(1/m) in(4), so it has no effect on this conclusion.

For an elementary version of the last spectral step, an eigenvalue of
T_delta-V_delta below -kappa corresponds to a variational direction of
the dual operator (T_delta+kappa)^(-1/2) V_delta
(T_delta+kappa)^(-1/2) above1: substitute
v=(T_delta+kappa)^(1/2)u in its quadratic form. This dual is A*A, while
the displayed K_delta is AA*, with A=V_delta^(1/2)
(T_delta+kappa)^(-1/2). Their nonzero spectra and multiplicities agree.
Thus counts below -kappa agree with counts above1 in(7), with multiplicity.
Norm convergence of compact self-adjoint operators preserves this count
whenever1 is not a limiting eigenvalue. Applying this at thresholds on
both sides of each isolated continuum negative level proves convergence
of the ordered levels and their multiplicities. It also excludes extra
levels in any fixed compact negative energy interval. Levels accumulating
at zero are treated by fixing a finite shell menu before taking the limit.

This argument also proves the assertion for the expanding periodic boxes:
their global constant mode contributes no nonvanishing extra localized
Birman–Schwinger mode. Local compactness and the uniform inner/tail bounds
already control it. No unqualified fixed-box or arbitrary aspect-ratio
limit is used.

## 5. The limiting levels and standing modes

For -Delta-beta/r separate spherical harmonics Y_lm and write the radial
wavefunction as u(r)/r. At energy -k², the radial equation is

    -u''+[l(l+1)/r²-beta/r]u=-k²u.

Put t=2kr and u=t^(l+1)e^(-t/2)v(t). Direct substitution gives

    t v''+(2l+2-t)v'+[beta/(2k)-l-1]v=0.

The origin-regular solution is square-integrable at infinity precisely
when its series terminates: beta/(2k)-l-1=q is a nonnegative integer.
Otherwise the regular confluent-hypergeometric solution has an e^t term,
while its exponentially decaying independent solution violates the
origin form-domain condition. Thus n=q+l+1, k=beta/(2n), and l=0,...,n-1.
Summing 2l+1 gives n². This proves(2); section4 transfers these ordered
spatial shells to the native-source description. The conclusion proved
here is ordered eigenvalue and shell-multiplicity convergence. A chosen
basis of degenerate eigenfunctions need not converge; no numerical
profile-convergence claim is needed or made.

The source occupation theorem applies to the derived symmetric completed-
macro-attempt transfer. On each positive internal band its Euclidean
operator is -Delta_time,a+S_a, with S_a strictly positive at finite volume.
Its scalar slice covariance has positive transfer multiplier
exp[-2asinh(a sqrt(s)/2)] for a spatial eigenvalue s. Reflection-positive
Gaussian reconstruction therefore yields a standing mode with that spatial
profile and the energy E_a(s). This is a description derived from native
full occupation moments, not a primitive Hamiltonian phase installed in
the event. The same source-clock interpretation as the free-field bundle
is used, with the following explicit finite-time error. Section3 gives
H_delta>=-C uniformly, hence S_a>=lambda/2 for sufficiently small a.
For any fixed shell, omega=2asinh(a sqrt(s)/2)>=c_lambda a. On the
declared cycle H=2m³+1, the relative image correction to its slice
covariance is O(exp(-c_lambda m²)). The exact first-hit multiplier is
cosh(omega(H/2-1))/cosh(omega H/2); its reconstructed energy differs
from omega/a by at most2omega exp(-omega(H-2))/a, also exponentially
small for these bands. Both errors are o(a²), smaller than binding.
The weaker uniform massless clock bound alone would not justify that
scale; the positive-band gap supplies the needed stronger estimate.
Formula(3) follows from(2) and
E_a'(lambda)=1/[2sqrt(lambda)sqrt(1+a²lambda/4)].

The lambda=0 internal band remains unchanged by feedback, because the
internal Laplacian annihilates it. Its free Maxwell sector can coexist
with these bound scalar bands on the same outer-source architecture.
This alone does not identify the binding as electric charge coupling.

## 6. A computable source diagonal; no fitted Coulomb strength

G(0) in(2) is a computable native constant. Its Fourier integrand follows
from the same six spatial neighbor counts. Since
4 sum sin²(p_j/2)>=4|p|²/pi², a ball of radius rho contributes at most
C rho; outside it the integrand has a computable derivative bound.
Rational Riemann sums and this bound give certified intervals for G(0),
c_* and every target e_n. No measured Rydberg or fitted gamma is used.

For each requested finite shell menu and tolerance, enumerate m=1,2,...
with B=m³. The finite h_B is the unique solution of a rational Dirichlet
linear system obtained from the native transition counts. The controller
does not read that solution. The finite S_a has algebraic entries; rational
interval isolation or exact inertia with algebraic coefficients computes
its ordered eigenvalues, including multiplicities. Use interval precision
tied to the outer enumeration and a strict interior tolerance margin;
do not wait forever for equality at one precision threshold.

Sections2–5 prove termination of the search for each finite shell menu.
The search chooses a declared source resource, never a coupling parameter.
An increasing menu n<=r and tolerance2^-r gives one explicit source family.
This is not a claim of practical polynomial complexity, nor a claim that
these increasingly large eigensystems have already been executed.

## 7. Bounded native execution and charged costs

Exact preparation/probes terminate almost surely with finite expected cost.
An inner relocation path has at most6B edges, each taking mean2 attempts.
For an incident probe, each of the six spatial increments has probability
1/22 per event. Consequently |X_t|²-(3/11)t is a martingale up to absorption.
Stopping on the target or cubic boundary and using |X_tau|²<=3B² gives
E[tau]<=11B² uniformly over starting points. This also holds after any
unabsorbed history. Therefore

    Pr(tau>22B² j)<=2^-j.

For a finite observation of K outer macro-attempts and allowed history
error epsilon, choose probe cap22B² ceil(log2(2K/epsilon)), and per-edge
relocation cap ceil(log2(12BK/epsilon)). A union bound couples the bounded
protocol to the exact one with discrepancy at most epsilon over its full
outer history. Every failed preparation/probe is retained as a hold.
Charge all inner calls and route attempts; the displayed caps also supply
a deterministic finite resource bound.

For finite moments of degree d, a bound on the finite record observable
times this history discrepancy controls the moment error. Choose epsilon
after the block length and observable bound, so the cutoff is small enough
for that menu. Then use the existing fixed-graph all-moment occupation
limit with the actual contact term L_eff^+-Pi/(4W), followed by volume,
spacing and the spectral limit above. This is a constructive diagonal,
not an unjustified exchange of infinite observation and cutoff limits.

## 8. Exact scope

The theorem supplies hydrogen-like discrete spatial shells and standing
energies for a specified native first-hit feedback source, with an unfitted
source-derived coefficient and controlled full limit. It does not identify
an electron, proton, physical electric charge, fine structure or laboratory
Rydberg constant. The target and feedback architecture are declared. Owner
pin, another-group review and an OpenWave registered task remain external
acceptance steps. The finite calibration verifies the local interaction;
it is not silently promoted to an executed large-shell production run.
