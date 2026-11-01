# A finite-density interacting ground representation at fixed lattice spacing

The original zero-chemical-potential Weyl/contact action on the infinite
three-dimensional lattice has a **translation-invariant, gauge-invariant
ground state and a positive-energy GNS representation**. For fixed
`c>0,U>0`, one such state has occupation density between `1/4` and `1` per
site. Its represented local dynamics retains the nonzero contact
interaction. No half-filled reference or chemical-potential shift is used.

This is an existence construction through compactness, with source-bound
finite premises and quantitative boundary errors. It does not compute or
prepare the limiting state, select a unique subsequence, or give a rate of
convergence of finite-volume states. The result is at **fixed spacing**;
it does not establish a relativistic continuum, a mass gap, scattering,
physical species or Einstein dynamics. Ground-state functional faithfulness
and stationarity under the repair instruments are not asserted.

## The original action and its available dynamics

At each `x in Z^3` there are two CAR modes. On a finite open rectangular box,

```
H_Lambda = sum_<x,y> (a_x^* h_xy a_y + h.c.) + U sum_x D_x,
h_(x,x+e_j) = i c sigma_j/2,       D_x=n_xup n_xdown.
```

Keep the entire local matrix algebra of dimension `4^|Lambda|`, including
every particle-number sector and odd operators. Coefficients are fixed,
finite `c>=0,U>=0`; the nonzero-density statement below requires both to
be positive. Each interaction term is even. A spatial bond has norm at
most `c`, and a contact term has norm `U`, giving
`-c E_Lambda I <= H_Lambda <= (c E_Lambda+U|Lambda|)I`.

The existing [full-Fock construction](../full_fock_field/README.md) already
proves a strongly continuous automorphism group `tau_t` on the quasi-local
CAR algebra, with convention `tau_t(A)=exp(itH)A exp(-itH)` in finite
volume. For `||A||<=1` supported on `s` sites in `[-r,r]^3`, its outward
boundary bound for a containing cube of radius `R>r` is

```
B_R(t) = min(2, 6 c s (2R+1)^2 |t| 2^(-(R-r)) [exp(48c|t|)-1]).
```

It is uniform in larger volume, number sector, onsite `U`, and bounded
time intervals. It applies to odd `A`: each exterior Hamiltonian term is
even, so the relevant disjoint commutator is zero. `local.py` independently
checks the full sixteen-element one-site Clifford basis against the
actual source on nested `3^3` and `4^3` boxes, without allocating their
full Fock matrices. Removing an incident bond is a nonzero negative
control. The exact local derivation is

```
Q(A)=[H,A]=sum_(Z intersects supp A) [Phi(Z),A],
||A^*Q(A)|| <= 2 ||A||^2 s (U+6c).
```

Its support lies in the one-neighbor enlargement of `supp A`. All other
even terms commute exactly, not approximately. The bounded calculators
are finite witnesses of these formulas; their evaluation caps do not
restrict the mathematical lattice or exhaustion.

## Full-ground references, even extension and spatial averaging

Let `Lambda_L={-L,...,L-1}^3`, with `V_L=(2L)^3`. Define `rho_L` to be
the normalized spectral projector of the **absolute** lowest eigenvalue
`E_L` of `H_Lambda_L`, taken over all `4^V_L` states. Finite-dimensional
spectral theory supplies this state for every `L`. No number sector,
Gaussian approximation, gap, unique eigenvector, numerical eigensolver
accuracy or thermal oracle is assumed.

Since `[H_Lambda_L,N]=0`, this projector is gauge invariant and even.
Extend it to the entire CAR algebra by the even product with the exterior
vacuum. One explicit finite-region definition orders inside modes before
outside modes and uses `rho_L tensor |vac><vac|`; evenness makes these
restrictions consistent with the graded CAR inclusions. The resulting
state `omega_L` restricts exactly to `rho_L` inside and is invariant under
the finite-box dynamics. Every exterior zero-energy mode remains present;
this extension is not a conditioning or a claim that the exterior vacuum
is the interacting ground.

With spatial translations `theta_x`, form

```
nu_L(A) = (1/V_L) sum_(x in Lambda_L) omega_L(theta_x(A)).
```

The local CAR algebra has a countable dense collection of rational
matrix entries. Their state values are bounded. Successive subsequence
extraction and a diagonal argument produce a subsequence `nu_(L_j)`
converging on every local matrix entry. Positivity, normalization and
compatibility pass to the limit; the bound `|nu_L(A)|<=||A||` extends
the limit uniquely to a state `omega` on the norm-completed quasi-local
algebra. Thus this is a sequence-level weak-star construction, not a
presupposed infinite-volume state. Compactness gives no computable
subsequence or state-convergence rate.

Gauge invariance passes immediately. For a fixed integer translation
`z`, cube averaging gives

```
|nu_L(theta_z(A))-nu_L(A)|
 <= ||A|| |(Lambda_L+z) symmetric_difference Lambda_L|/V_L ->0.
```

Consequently `omega` is translation invariant. `folner_translation_error`
retains the exact analogous ratio for odd-sided cubes. Both exhaustions
have the same vanishing-boundary property.

## Ground inequality and stationarity, including odd probes

For a local `A` whose one-neighbor enlargement lies inside the finite box,
the finite-volume spectral minimum gives exactly

```
omega_L(A^*[H,A])
 = Tr rho_L A^*(H_Lambda_L-E_L) A >=0.
```

This argument uses the **full** ground minimum, so it holds for charged
odd operators and arbitrary mixtures of parity and particle-number
changes. It would not follow for an arbitrarily selected canonical
ground when `A` changes number.

For support in `[-r,r]^3`, the fraction of translates whose one-neighbor
enlargement fails containment is at most
`1-(max(0,L-r-1)/L)^3`. Each bad contribution can be complex: bound its
real part below by `-2||A||^2 s(U+6c)` and its imaginary magnitude by
the same constant. Good contributions are real and nonnegative, so the
averaged lower defect and imaginary part both vanish. The
local operator `A^*[H,A]` is fixed as `L` grows, so weak-star convergence
proves

```
omega(A^*[H,A]) >=0       for every local A, including odd A.
```

Stationarity can be proved directly without assuming that this local
inequality is already a theorem about an unbounded generator. Fix
`R>r`. A fraction at least `(1-R/L)^3` of translates have their entire
radius-`R` cube inside `Lambda_L`. For each good translate, compare
finite-box and infinite dynamics to the same intermediate cube: the
operator error is at most `2B_R(t)||A||`. For the remaining fraction
use `2||A||`. Since `omega_L` is exactly invariant under finite-box
dynamics,

```
|nu_L(tau_t(A))-nu_L(A)|
 <= ||A|| min(2, 2B_R(t)+2[1-(1-R/L)^3]).
```

First let `L_j` tend to infinity, then `R` tend to infinity. This proves
`omega(tau_t(A))=omega(A)`, uniformly on each compact time interval for
fixed local `A`. Density extends invariance to every quasi-local
observable. `averaged_dynamics_error` is an exact finite bound on this
stationarity comparison, **not** on the unknown difference `nu_L-omega`.

## Positive-energy GNS representation

Let `(pi,H_omega,Omega)` be the GNS representation of `omega`. Invariance
defines the strongly continuous unitary group

```
U_t pi(A)Omega = pi(tau_t(A))Omega,     U_t=exp(it L_omega),
U_t Omega=Omega.
```

Here Stone's generator `L_omega` is an energy above the ground
representation. It is not formed by adding a chemical potential or
changing the source action.

To prove `L_omega>=0`, for local `A` consider the finite correlations

```
F_L(t) = (1/V_L) sum_x omega_L(theta_x(A)^*
                                  tau_t^Lambda_L(theta_x(A))).
```

Each summand has a positive spectral measure supported on `[0,infinity)`:
on a finite containing algebra the exterior Hamiltonian is zero and the
inside density is supported at its absolute minimum `E_L`, so the
correlation is a sum of nonnegative weights times
`exp(it(E_n-E_L))`. The same argument includes exterior components of
`theta_x(A)`, all charged states and odd `A`.

The preceding good/bad-translate estimate, multiplied by `||A||`, proves
that `F_(L_j)(t)` tends to
`omega(A^*tau_t(A))=<pi(A)Omega,U_t pi(A)Omega>` uniformly on compact
time intervals. Weak-star convergence is uniform on the norm-compact
set `{A^*tau_t(A):|t|<=T}` by a finite epsilon-net and the uniform state
norm bound. All correlations are bounded by `||A||^2`.
For a smooth nonnegative function `chi` compactly supported on negative
energies, let `g(t)=(1/2pi) integral chi(E) exp(-itE)dE`. It is integrable.
Every finite correlation has `integral g(t)F_L(t)dt=0`. Compact-time
convergence and the uniform bound pass this equality to the limit.
The spectral theorem then gives
`<pi(A)Omega,chi(L_omega)pi(A)Omega>=0`.
Every such negative-energy test vanishes. Since local GNS vectors are
dense, the negative spectral projection is zero and `L_omega>=0`, with
`L_omega Omega=0`. This supplies the positive-energy condition for the
actual three-dimensional interacting lattice action.

## A nonzero density, without half-filling or a chemical potential

The existing exact eight-site cube vector `G` has four hidden-up
particles, no contact, and energy `-2sqrt3 c`. Tile every even-sided
box by disjoint `2^3` cubes and use their graded product. Each tile has
fixed number, hence zero odd expectation; every cross-tile hopping has
expectation zero. This is a trial vector, not a proposed large-box
eigenstate. It proves

```
E_L/V_L <= -sqrt3 c/4.
```

For rectangular hopping, write `h=c sum_j sigma_j D_j` with Hermitian
Cartesian difference matrices `D_j`. They commute, each has norm at
most one, and the Pauli anticommutators cancel the mixed terms:
`h^2=c^2 sum_j D_j^2<=3c^2 I`. Therefore
`H>=-sqrt3 c N` for `U>=0`. The full-ground comparison at `c>0` yields
`Tr(rho_L N)/V_L>=1/4`.

For the upper bound, the already proved all-box onsite SU(2) frame
turns the original action into real, spin-independent bipartite hopping
while preserving `N,D`. The simultaneous staggered particle-hole map
on both hidden components gives

```
P H P^*=H+U(V_L-N),
E_L(2V_L-N)=E_L(N)+U(V_L-N).
```

At `U>0`, any `N>V_L` state at the absolute minimum would have a
strictly lower-energy partner. Thus no such ground component exists,
and `Tr(rho_L N)/V_L<=1`. `particle_hole_certificate` checks the full
two-site CAR identity, with its all-box extension following bond by
bond from real bipartite hopping. No half-filling conclusion follows.

Spatial averaging is essential here:
`nu_L(n_0up+n_0down)=Tr(rho_L N)/V_L` exactly. Hence the same cluster
state has

```
1/4 <= omega(n_0up+n_0down) <= 1.
```

This fixed-spacing state has finite, nonzero density. It is outside the
fixed-particle-count domain of the earlier free Weyl continuum theorem.
The density bounds neither select a density nor control its behavior
under `a->0`.

## Ground energy density and finite-box surface brackets

Let `E_m` denote the absolute full-Fock minimum on an even-sided cube
of side `m`, at the same fixed `c,U`. If even `ell` divides `m`, tile
the large cube by side-`ell` cubes. There are exactly
`3(m/ell-1)m^2` spatial bonds between tiles. Deleting them yields a
sum of independent tile Hamiltonians, each bounded below by `E_ell`;
their total operator norm is at most `c` times this count. Conversely,
the graded product of gauge-invariant tile ground states has zero
cross-tile hopping expectation. Hence

```
E_ell/ell^3 -3c/ell+3c/m <= E_m/m^3 <= E_ell/ell^3.
```

For general even `m`, tile the side-`ell floor(m/ell)` inner cube.
The remaining `m^3-(ell floor(m/ell))^3=O(ell m^2)` sites meet at most
six bonds each. All terms touching them have norm at most
`(6c+U)` times that remainder count. This controls both lower and upper
comparisons; the normalized remainder vanishes for fixed `ell`.
Extensive stability bounds the sequence. For every even `ell`, its
limsup and liminf differ by at most `3c/ell`, so the energy density
exists and satisfies

```
e_infinity = lim_(even m->infinity) E_m/m^3,
E_ell/ell^3-3c/ell <= e_infinity <= E_ell/ell^3.
```

Define `h_x=U D_x` plus one half of every incident spatial hopping
bond. The exterior-vacuum extension has zero cross-boundary hopping
expectation, so the spatially averaged reference obeys
`nu_L(h_0)=E_(2L)/(2L)^3` exactly. The same ground subsequence therefore
satisfies `omega(h_0)=e_infinity`. The earlier trial and number bounds
also give `-sqrt3 c<=e_infinity<=-sqrt3 c/4` when `c,U>0`.

`energy.py` evaluates the surface bracket given a separately certified
finite full-ground energy interval. The retained side-two input is
the proved original cube value `E_2=-sqrt3/4` at `c=1/8,U=1`;
it is not extrapolated as the exact cube minimum at arbitrary couplings.
Larger finite-ground energies have not been computed here. These are
conditional computable energy-density brackets, not a rate of convergence
of the states, a unique reference, or particle rest masses.

## The interaction survives representation; the free sea does not suffice

The exact even local witness is

```
[H_U-H_0,a_0^*a_2] = U a_0^*a_2(n_1-n_3) !=0.
```

Its Clifford expansion has degree-four terms. A nonzero unital GNS
representation restricts to a unital representation of each full local
matrix algebra; matrix algebras are simple, so this restriction is
injective. Thus this local operator remains nonzero in the ground
representation. Equivalently, the quasi-local CAR algebra itself is
simple. This does **not** require the state functional to be faithful.
On the dense invariant domain of local GNS vectors,
`[L_omega,pi(A)]pi(B)Omega=pi([H,A])pi(B)Omega`.
Consequently the represented derivation contains the original quartic
interaction. This is not a proof of nontrivial scattering, a connected
vacuum four-point function, or an interacting relativistic continuum.

There is also an exact countercontrol to replacing this construction by
the filled free negative band. For a pure finite Slater covariance `P`
with onsite two-mode blocks `P_xx=I2/2`, Wick's determinant identity gives

```
Var(D_x)=3/16,
Cov(D_x,D_y)=-Tr(P_xy P_yx)/4+|det P_xy|^2  (x!=y).
```

Since `P^2=P`, `sum_(y!=x)Tr(P_xy P_yx)=1/2`, and therefore

```
Var(D)=V/16+sum_(x!=y)|det P_xy|^2 >=V/16.
```

If this Slater vector is an `H_0` eigenvector, then
`Var(H_0+UD)=U^2 Var(D)>0` for `U>0,V>0`. It cannot be a stationary
pure interacting reference. This does not exclude correlated grounds.

`filled_sea.py` constructs all 3,844 nonzero occupation amplitudes of
the source cube's filled negative band over `Q(sqrt3)/1296`. Exact
complete-Fock arithmetic gives `N=8`, `<D>=2`, `Var(D)=2/3`,
`<H>=2-sqrt3/2` and `Var(H)=2/3`. An independent covariance calculation
gives the same variance. Wick is used only for this known Slater input;
no interacting Gaussian evolution is assumed.

## An analytic local retarded response on the same ground subsequence

For bounded local operators `A,B`, including odd CAR probes, define

```
G_AB(z)=-i integral_0^infinity exp(i z t)
                         omega({tau_t(A),B^*}) dt,       Im z>0.
```

Its finite comparator averages the corresponding anticommutator over
the same finite full-ground states, with finite-box dynamics. The
norm bound on the anticommutator is `2||A||||B||`. Consequently, at
`Im z>=eta>0`, each time tail after `T` is bounded by
`2||A||||B|| exp(-eta T)/eta`. Differentiating under the integral gives
`|d_z^n G_AB(z)|<=2||A||||B|| n!/eta^(n+1)` for every nonnegative
integer `n`; the response is analytic in the upper half-plane.

The same good/bad-translate estimate controls each finite-to-infinite
dynamics replacement. The remaining state discrepancy on `[0,T]` is

```
epsilon_L(T)=sup_0<=t<=T
 |nu_L({tau_t(A),B^*})-omega({tau_t(A),B^*})|.
```

It tends to zero along the **same** ground subsequence: the displayed
local-in-time family is norm compact, so weak-star convergence and a
finite epsilon-net suffice. It has no supplied numerical convergence
rate. The complete response comparison is bounded by

```
2||A||||B||T min(2,2B_R(T)+2bad_fraction)
  + T epsilon_L(T) +4||A||||B|| exp(-eta T)/eta.
```

Take `T` large, then `R` large, then the ground subsequence index large.
This proves uniform convergence on compact subsets of `Im z>0`.
`retarded.py` retains the exact tail/derivative bounds and a comparison
calculator that requires an explicitly supplied `epsilon_L(T)` bound;
it never silently sets that unknown discrepancy to zero. Its rational
exponential upper bound is the reciprocal of a positive Taylor partial
sum, hence rigorous without a floating exponential premise.

This supplies an analytic interacting local response in the proved
ground representation. It computes no infinite-volume response values
and establishes no real-axis continuation, isolated pole, mass gap,
spatially local inverse, metric or relativistic dispersion.

## Replay and remaining scope

```
.venv/bin/python -W error -m pytest -q model/thermodynamic_ground
.venv/bin/python -W error -m model.thermodynamic_ground --write
.venv/bin/python -W error -m model.thermodynamic_ground --check
```

The receipt binds the action, previous thermodynamic dynamics, exact
local/source checks, density inequalities, countercontrols and this proof.
It records an existence theorem, not retained amplitudes or local density
matrices of an infinite-volume state. Native preparation/readout, a
specific ground functional, repair-preserving continuation, uniqueness,
gap or clustering, finite-density continuum estimates, gauge/species
identification and geometric Einstein premises still require separate
constructions.
