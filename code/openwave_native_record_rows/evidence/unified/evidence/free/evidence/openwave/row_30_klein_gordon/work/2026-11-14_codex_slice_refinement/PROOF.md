# A uniform massive slice limit from the native record resolvent

This is a conditional theorem extending the E0/E1 space-slice calculation. It is
not a new native execution or an operational physical-time identification. Its
source-derived input is the one-unit native transfer on a specified routed graph.

## 1. Prior theorem and new content

`openwave-push/b/E1_matter/SOLUTION.md` sections 3 and 4 and `MEMBER.md` define the
port-diagonal four-axis routing and establish the slice equation. Its companion
`VERDICT.md` already observes that massive bands approach a common speed when the
external seam weight grows. E0 gives the finite-weight derivative
omega/sinh(omega). Those results are reused here.

The extension is a uniform, explicit O(a²) bound for every native internal band,
a common-speed bound, an exact group-speed bound, and the normalized covariance
limit along one declared refinement family. It shows exactly why the unit-weight
failure is not a failure of all refined members. It does not remove the routing,
refinement, ruler or time-interface premises.

## 2. Source family and the observable

Each cell has the twelve internal ports and the thirty internal seams of the
pinned native graph. Copies are indexed by z in Z^4. For every port p, join
(z,p) to (z+e_mu,p), mu=1,2,3,4. The internal seam weight is 1; each external
seam weight is w=a^(-2). Both directions occur through undirected seams. On a
finite periodic quotient, choose a seam independently with probability
proportional to its weight and apply the pinned fair native endpoint event.

For a=1/m, m a positive integer, the weights can be realized by m² parallel
addresses for each external seam and one address per internal seam, all chosen
uniformly. This changes the scheduler family, not the endpoint function. It is
declared architecture, not a fitted value of a physical coupling.

Let W be the total seam weight on a finite quotient. For one unit, the exact
branch-count identity is

    P_a = I - L_a/(2W),
    L_a = L_int + a^(-2) sum_(mu=1)^4 (2 - S_mu - S_mu^(-1)).

The first equality follows directly from the native half-exchange probability on
each seam. An observer holding the known initial reference and signed current
tape reconstructs the occupation at every event. The accumulated occupation has
expectation sum_n P_a^n. Project out its constant mode, or use a nonzero internal
band/spatial momentum. Then its infinite-event limit is 2W L_a^(-1). A finite
record horizon leaves the ordinary geometric remainder. No oscillator step or
mass coin has been inserted.

The fourth coordinate is a declared Euclidean slice index. It is not the source
event order: a native event can move either way in that coordinate. Converting
the slice index to a physical time requires a separate same-history interface.
This theorem takes three momenta p conjugate to x=a(z1,z2,z3) and a slice
separation t=a z4. Those coordinate assignments are stated inputs.

## 3. Exact native internal spectrum

For the pinned icosahedral incidence,

    spec L_int = {0 (multiplicity1),
                  5-sqrt(5) (multiplicity3),
                  6 (multiplicity5),
                  5+sqrt(5) (multiplicity3)}.

This is the already-certified twelve-port spectrum. One source is the exact
characteristic polynomial of 60P in the record-foundations review:

    det(xI-60P) = (x-60)(x-54)^5((x-55)^2-5)^3,

together with P=I-L_int/60. No laboratory mass appears. The three nonzero
square-root gap ratios in the limit are

    1 : sqrt(6/(5-sqrt(5))) : (1+sqrt(5))/2.

They are internal-mode ratios, not an identification of electron, muon or tau.
The prior lepton comparison fails and is not repaired by this theorem.

## 4. Exact finite-spacing dispersion read from slices

For an internal eigenvalue lambda, define

    q_a(p) = (4/a²) sum_(j=1)^3 sin²(a p_j/2),
    A_a(lambda,p) = lambda + q_a(p),
    E_a(lambda,p) = (2/a) asinh(a sqrt(A_a)/2).

Away from A_a=0, the decaying solution of the slice resolvent recurrence is
exp(-a E_a |n|). Indeed,

    (A_a+2/a²)g_n - (g_(n+1)+g_(n-1))/a² = 0, n!=0,
    4 sinh²(a E_a/2)/a² = A_a.

The native record supplies the resolvent. E_a is its measured decay parameter,
not an assertion that a unitary runs in event time. On the infinite slice line,
normalizing the resolvent by the slice cell measure gives

    C_a(t,p;lambda) = a exp(-E_a |t|)/(2 sinh(a E_a)),  t in aZ.

For raw accumulated occupations, this normalization is division by 2W a after
the internal-band and spatial-Fourier projections. These are explicit record
statistics with supplied coordinate normalization. At A_a=0 the infinite-line
covariance diverges; that zero mode is excluded from the covariance assertion,
not silently discarded from the source. The dispersion formula still has E_a=0.

## 5. Uniform quadratic refinement theorem

For every a>0, every p in R^3 and every native internal eigenvalue lambda>=0,

    0 <= (lambda+|p|²) - E_a(lambda,p)²
      <= (a²/12) [sum_j p_j^4 + (lambda+|p|²)²].                 (1)

Proof. The elementary cosine remainder gives

    0 <= |p|²-q_a(p) <= (a²/12) sum_j p_j^4.

For x>=0, 0<=x-asinh(x)<=x³/6. To prove the upper bound, integrate
1-(1+x²)^(-1/2)<=x²/2, which follows from the tangent-line bound of the convex
function u -> (1+u)^(-1/2). Since asinh(x)<=x,

    0 <= x²-asinh(x)² <= x^4/3.

Set x=a sqrt(A_a)/2 and multiply by4/a². This gives
0<=A_a-E_a²<=a² A_a²/12. Add the cosine bound and use
A_a<=lambda+|p|². QED.

Let Lambda=5+sqrt(5). On |p|<=K, all twelve native bands therefore satisfy

    sup |E_a²-(lambda+|p|²)|
      <= (a²/12) [K^4+(Lambda+K²)²].                            (2)

Thus the same declared family has the continuum scalar dispersion
E²=|p|²+lambda for every band, uniformly on each fixed momentum ball. The rest
gap is sqrt(lambda), determined by the internal native response in these units.
A lattice-spacing-dependent mass coin was not used; the relative scheduler
scaling w=a^(-2) is the explicit remaining construction choice.

## 6. A common infrared speed, and a separate group-speed bound

For lambda>0 let omega0=a E_a(lambda,0). The coefficient of |p|² in the rest
dispersion is

    c_a(lambda)² = omega0/sinh(omega0)
                 = asinh(x)/(x sqrt(1+x²)), x=a sqrt(lambda)/2.

For lambda=0 its continuous value is1. For all lambda>=0,

    0 <= 1-c_a(lambda)² <= a² lambda/6.                        (3)

Proof. Differentiate x sqrt(1+x²)-asinh(x): its derivative is
2x²/sqrt(1+x²). Its value is at most2x³/3. Divide by
x sqrt(1+x²)>=x. QED.

All native bands share the speed1 limit in the supplied rulers and slice clock.
For a=1/16, equation(3) bounds every deficit below1/200, since
5+sqrt(5)<7.68. This is a theorem-level example, not a new empirical run.
More generally w>=Lambda/(6 epsilon) suffices for a speed-squared deficit at
most epsilon. The old numerical threshold119.2 at one percent is compatible
with the sufficient upper threshold Lambda/0.06 from this bound.

There is also an exact finite-spacing bound wherever E_a is differentiable:

    |grad_p E_a|²
      = sum_j sin²(a p_j) /
        [a² A_a (1+a² A_a/4)]
      <= q_a/[A_a(1+a² A_a/4)] <= 1.                          (4)

Here sin²(2u)<=4sin²(u) was used. A maximal graph-step speed is a different
quantity and need not equal either this group-speed bound or the infrared slope.
None of these unitless statements yet identifies laboratory c.

## 7. Normalized two-point covariance limit

Let S=lambda+|p|² and delta_a=(a²/12)[K^4+(Lambda+K²)²]. Fix S>=sigma²>0,
|p|<=K, |t|<=T, a<=a0 and delta_a<=sigma²/2. Put
Emax=sqrt(Lambda+K²) and Emin=sigma/sqrt(2). Then

    |C_a(t,p;lambda) - exp(-sqrt(S)|t|)/(2sqrt(S))|
      <= delta_a/sigma [T/(2Emin)+1/(2Emin²)]
         + a² Emax cosh(a0 Emax)/12.                          (5)

Proof. Equation(1) gives E_a>=Emin and
|E_a-sqrt(S)|<=delta_a/sigma. The derivative with respect to E of
exp(-E|t|)/(2E) has magnitude at mostT/(2Emin)+1/(2Emin²).
For z>=0, sinh(z)-z<=z³ cosh(z)/6 by Taylor's integral remainder.
Therefore the extra prefactor z/sinh(z) changes exp(-E|t|)/(2E) by at most
a² E cosh(aE)/12. Substitute the uniform bounds. QED.

The right side is O(a²). The limiting covariance has Euclidean Fourier symbol
(p4²+|p|²+lambda)^(-1). Its usual analytic continuation is the free scalar
Klein-Gordon description with mass squared lambda. This is a mathematical
description of the limiting one-unit record covariance. It is not yet an
operational Lorentzian process, an implemented complex-effect instrument, or a
proof about many-unit conditional records. The E1 two-unit counterexamples to
the tested reflection-positive product reading remain on file.

## 8. What advances, and what remains

The unit-weight band-dependent speed is repaired mathematically along this
specific family with the explicit bound(3). The source determines the internal
gaps; an independent scalar mass term is unnecessary for this conditional
one-unit slice result. These are useful input statements for an observer-first
construction and do not need an installed Maxwell or Dirac Hamiltonian.

Still supplied: four-dimensional routing, external/internal scheduling ratio
as a function of refinement, spatial ruler a, the choice of a Euclidean slice,
the accumulation/projection decoder and infinite-volume/horizon limit.
Still owed: selecting or pinning that one candidate; a source observer's
same-history physical time and intervention correspondence; common matter/light
comparison; and the full joint conditional law beyond the one-unit covariance.
The theorem supplies neither two transverse polarizations for row29 nor a
source interaction and standing bound-state phases for row31.
