# A causal bound uniform in physical lattice refinement

The same gauge-coordinate field has a finite propagation bound uniform
as its assigned lattice spacing tends to zero, under an explicit scaling
of its action. The estimate concerns complete local observer records,
all occupations and quantum references. It does not construct a limiting
state, prove convergence of observables across spacings, or identify a
Lorentz light speed. Spatial lengths here are the declared lattice length
convention; no emergent metric is inferred from a commutator estimate.

## 1. Reused proof and fixed architecture

Use the [full-Fock gauge-coordinate locality proof](../particle_field_causality/PROOF.md)
and its bounded clipping of unbounded onsite generators. It applies to
the native four Dirac components and fifteen internal rows at every cell,
with coordinate links owned by their tail cells. All sixty CAR modes per
cell are retained. The [smooth interaction expansion](../particle_smooth_sign_field/PROOF.md)
supplies explicit exponentially summable interactions. Only its assigned
spacing and magnetic coefficient vary in the family below.

The [finite causal observer theorem](../../Lean/QFT/FiniteCausalObserverNet.lean)
describes exact remote repair under its locality premises. The present
Hamiltonian estimate supplies approximate remote-influence control with
explicit tails. The [refining lattice cone theorem](../../Lean/Geometry/RefiningLatticeCausalCone.lean)
concerns a different layered move law. Its exact reachability statement
is not transferred to continuous Hamiltonian evolution. No clock is
redefined as a scheduler layer or an apparatus gate count.

On an open cubic box, or a periodic box with graph distance including
wraparound, let a>0 be the assigned cell spacing. Fix delta>0 and kappa0>=0.
Declare

```
H_a = H_on(a) + dGamma(h_delta,a)
                  - kappa(a) sum_p Re Tr_M15(U_p)/15,
h_delta,a = [beta+sqrt(1+delta²) H_W(H_W²+delta² I)^(-1/2)]/a,
0 <= kappa(a) <= kappa0/a.                              (1)
```

H_on(a) contains the source electric Casimirs and any additional declared
coordinate-independent onsite even terms. The onsite part of dGamma(h)
is transferred into the interaction-picture generator, counted only once.
Their coefficients may vary without a common bound;
each finite-spacing onsite operator must have its stated self-adjoint
domain. The source electric coefficients are nonnegative. The locality
argument does not use their magnitude. A source member with kappa(a)
growing faster than kappa0/a is outside this certificate; failure of this
majorant is not a no-go result for that member's actual propagation.

This is a declared refinement family. It does not assert that a fixed
source preparation or magnetic coefficient from another calculation
automatically belongs to the family. The smoothing parameter remains
fixed; no discontinuous-sign limit is involved.

## 2. The spacing cancels from the physical propagation rate

Write the spatial weight of the smooth proof as lambda>1, and its
weighted interaction bound at spacing1, coefficient kappa0 as W_star.
Its kinetic terms scale exactly as 1/a and its magnetic bound is linear
in kappa. Consequently

```
sup_x sum_(Z contains x) |Z| ||Phi_Z(a)|| lambda^diam(Z)
 <= W_star/a.                                         (2)
```

The onsite interaction picture preserves every support, norm and parity.
Bounded spectral clippings of each onsite generator yield the same (2),
uniformly in the clipping ceiling. The finite-volume strong-unitary limit
therefore preserves the norm inequality, as in the owning proof. No
electric-energy expectation of the input state is needed.

For disjoint homogeneous operators A_X and B_Y, set
`M_xy=sum_(Z contains x,y)||Phi_Z||`. Equation(2) gives
`sup_x sum_y M_xy lambda^dist(x,y)<=W_star/a`. The triangle inequality
for graph distance implies the weighted row sum of M^n is at most
`(W_star/a)^n`. Hence
`sum_(x in X,y in Y)(M^n)_xy<=|X|lambda^(-D)(W_star/a)^n`.
The interaction chain expansion therefore yields

```
||[alpha_a,t(A_X),B_Y]_graded||
 <= 2 ||A_X|| ||B_Y|| |X| lambda^(-D)
                         [exp(2 W_star |t|/a)-1],
D=dist_graph(X,Y).                                    (3)
```

The weighted row sum bounds the sum over every destination cell in Y;
there is no factor proportional to the cardinality of the remote support.
The zero-length chain vanishes because the supports are disjoint.
When either operator is even the bracket is the ordinary commutator;
for two odd operators it is the anticommutator. Local physical detector
effects and admitted interventions are even, including any required
parity-compensating apparatus.

Put d=aD. Dropping the minus-one term in (3) gives

```
||[alpha_a,t(A_X),B_Y]_graded|| /(||A_X||||B_Y||)
 <= 2|X| exp[-(log(lambda)d-2 W_star|t|)/a].             (4)
```

The cone speed `v=2 W_star/log(lambda)` is independent of a. It is a
conservative propagation bound in the declared coordinates and clock,
not a measured light speed, an isotropic dispersion, or a boost action.
Euclidean separation can be used as a weaker lower bound because
Euclidean distance does not exceed the assigned Manhattan distance.
At finite spacing the tails need not vanish exactly.

## 3. A sharper finite-range architecture

The directly declared ordinary Wilson action is also allowed:

```
h_Wilson,a=(beta+H_W)/a.
```

It is not identified with the smooth action or a chiral overlap theory.
Its onsite contribution goes into H_on. A nearest-neighbor bond has
one-body off-diagonal block `-beta P_mu^- R(U)/a`, with
`P_mu^-=(I-gamma_mu)/2`. The native projector has spin rank2, so the full
block has thirty singular values1/a and thirty zeros. On the complete
two-site one-body space of dimension120 the bond has eigenvalues
`+1/a` and `-1/a`, each of rank30, and zero of rank60. Its full-Fock
second quantization consequently has norm30/a: the largest occupied sum
is the sum of all thirty positive one-body eigenvalues. This counts all
occupations, not a selected one-particle sector. Coordinate unitarity
preserves these singular values on the actual quantum link space.

There are six incident bonds, each of support size2 and diameter1, and
twelve incident plaquettes, each of support size4 and diameter2. Thus
for any fixed lambda>1,

```
W_star = 360 lambda + 48 kappa0 lambda².               (5)
```

The exact source spin matrices and all fifteen internal columns are
checked in `native_wilson_bonds()`. The smooth family's larger bound
remains valid on its own architecture; (5) does not replace it silently.

## 4. Complete apparatus records and growing physical supports

Compare two experiments on the same full state and Hamiltonian member:
one applies a remote even unitary U_Y at time0, and the other applies
identity. The remote operation can include a private Stinespring
apparatus, which is traced in the observer-output comparison. No remote
record is communicated directly to the observer. Arbitrary untouched
references correlated with the input are retained in the output.

The observer uses k complete local Stinespring interventions, a shared
feedback policy, and a final local output over duration T. Every local
record is retained coherently until readout. All field supports lie in X.
Move U_Y past the Heisenberg intervention product. Each of the k
interventions contributes once on each side of a dual output effect;
the final effect contributes once. Applying (3) at times no greater than
T and effect duality gives the conservative bound

```
half_trace(output_U,output_I)
 <= min(1,2(2k+1)|X| exp[-(log(lambda)d-2 W_star T)/a]). (6)
```

The proof is uniform in input state, fermion number, spectator dimension
and remote-support cardinality. Shared input correlations do not violate
this bound. A global energy projector, global Fourier momentum menu or
random-site effect n_x/N is not a local intervention merely because its
reported value is assigned a spatial label.

Assume 0<a<=1 in the chosen length unit and X is contained in a physical
coordinate cube of half-width R. An alignment-safe bound is
`|X| <= (floor(2R/a)+3)^3 <= (2R+3)^3/a^3`.
Assume also `k(a)<=K/a^p`, for integer K,p>=0. Set

```
mu=(lambda-1)/lambda <= log(lambda),
g=mu d-2 W_star T,
C=2(2K+1)(2R+3)^3,    q=p+3.
```

Equation(6) is bounded by `min(1,C a^(-q) exp(-g/a))`.
At every fixed positive margin g, the exponential dominates the support
and apparatus growth as a tends to zero. This is uniform in exterior
volume and holds for every real spacing below the certified cutoff.

The geometric separation premise concerns the actual assigned cell
supports: `a dist_graph(X_a,Y_a)>=d`. If the supports are instead chosen
by intersecting cells with continuum regions, enlarge those regions by
their cell radius before using d. For example, a continuum separation
d0 and at most3a total Manhattan enlargement permit d=d0/2 whenever
a<=d0/6. The code does not replace this support condition by an unchecked
continuum distance.

## 5. Exact finite resolution certificates

The derivative of `a^(-q)exp(-g/a)` has the sign of `g-q a`.
It therefore increases on `0<a<=g/q`. For rational data, use
`exp(-g/a)<=2^(-floor(g/a))` and
`C a^(-q)<=2^ceil(log2(C/a^q))`. A sufficient binary exponent is

```
ceil(log2(C/a^q))-floor(g/a).                           (7)
```

The dyadic planner decreases a until both a<=g/q and (7)<=-bits.
Its termination follows from exponential domination. Its output is valid
for every smaller positive real spacing by monotonicity, not only for the
tested dyadic members. The code rejects g<=0 as insufficient for this
certificate; it does not call an inside-cone experiment acausal.

For the Wilson choice lambda=2, kappa0=1/100, R=1/10, T=1/1000,
d=4 and at most two intermediate interventions, W_star=18048/25 and
g=1738/3125. The planner gives a<=1/128 for a requested20-bit bound.
The complete record error is at most2^-41 at that cutoff and at every
smaller spacing. An example on the positive-smoothing family keeps T=1
and chooses a large d giving g=1; its separation is stated explicitly.
Neither example allocates the large field or executes its full apparatus.

## 6. Representation and physical scope

The estimate is proved on the coordinate/Fock parent before a global
coordinate-word cutoff. Such a cutoff can be nonlocal. If finite-box,
polynomial, representation, preparation or numerical implementations
have complete-experiment half-trace errors epsilon_left and epsilon_right,
their implemented output difference is bounded by (6) plus both errors,
capped at one. A common small error in only one experiment is insufficient.
Each record comparison retains the same source clock; implementation
overhead is a separate resource unless explicitly included in that clock.

At each fixed a, the [smooth approximation and box planners](../particle_smooth_sign_field/PROOF.md)
and [finite-word history planner](../particle_local_history_family/PROOF.md)
can supply these errors when their explicit initial-state, apparatus and
filtration premises hold. For example, allocating a total error at most a
to each complete experiment gives an implemented influence bound of
`min(1,C a^(-q)exp(-g/a)+2a)`, which also tends to zero for g>0.
This permits polynomial order, simulation box and word cutoff to grow with
refinement. It does not presume a single fixed cutoff has local dynamics,
or that an arbitrary unspecified input has a known preparation-tail bound.

If independently identified observables and states have convergent joint
record laws under refinement, (6) forces the limiting operational remote
influence to vanish outside this conservative cone. The bound alone
supplies neither those limits nor positivity of a selected energy vacuum,
Lorentz transformations, physical chirality, a mass shell or Einstein
dynamics. Its output is one uniform causal premise for such a construction.
