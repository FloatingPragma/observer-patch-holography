# Gauge response of the same interacting pair

This package couples the existing Weyl/contact CAR action to a declared
external U(1) link field. It joins the complete number current to the
energy-dependent bound-pair response, including the separated-pair norm.
The leading pair hopping has twice the link phase; an exact Schur complement
and a remainder bound specify what that approximation leaves out.

The links, action, coordinate speed and Hamiltonian model time are supplied.
The integer charge here is CAR particle number. Neither its normalization
nor its name identifies electromagnetic charge, a dynamical gauge field,
Lorentz spin, a native defect, a physical mass or a chiral continuum.
All tastes and the entire two-particle scattering space remain present.

## Local coupling, current and source implementation

For an oriented neighboring site pair `x<y` in direction `a`, write

```
h_xy(A) = exp(i A_xy) i c sigma_a/2,
B_xy = sum_(s,t) c_(x,s)^* (h_xy)_(s,t) c_(y,t),
H(A) = sum_(x<y) (B_xy+B_xy^*) + U sum_x n_(x,0)n_(x,1).
```

This is the same action as `interacting_field` at zero link field. The
`GaugedField` constructor uses a finite induced cubic site set with open
boundary; absent bonds are not periodic wrapping links. It acts on every
state of the exact `N=2` sector.

Let `G(alpha)=exp(i sum_x alpha_x N_x)`. CAR identities give

```
A'_xy = A_xy + alpha_x - alpha_y,
H(A') = G(alpha) H(A) G(alpha)^*,
J_(x->y) = partial H/partial A_xy = i(B_xy-B_xy^*),
i[H,N_x] + sum_y J_(x->y) = 0.
```

Reversing orientation negates the current. The contact term commutes with
every `N_x`; internal currents cancel when a region is summed. These are
operator identities, including nonflat links with nontrivial square
holonomy. Site phases do not remove such flux. Their CAR derivation applies
to every finite Fock sector, although this package's global evaluator is
the declared complete `N=2` sector.

The actual two-mode source word in `gaussian_source` accepts the complex
hopping coefficient `exp(i A_xy) z`. All four local Fock input columns,
including vacuum and double occupation, satisfy the same gauge covariance.
The contact word and occupation pointer commute with local number phases.
Consequently every finite split evolution and every complete occupation
history is gauge covariant, with arbitrary untouched reference factors.
This uses the existing fermionic routes and source compiler; it introduces
no additional source primitive and discards no zero outcomes. Link phases
change the supplied word parameters. The emitted primitive counts, loads,
exports and fermionic routes are unchanged for the same bond schedule.
Their floating residual and hardware errors retain the existing separate
admission conditions.

The exact certificate constructs fourth-root link values directly as
Gaussian integers at `c=1/8`, `U=1`, rather than rounding trigonometric
values. On the complete 1,431-state radius-one cube it checks all 27 local
Ward identities, 54 bond currents, full-sector gauge covariance and the
pair-transfer identity below. Independent tests use arbitrary complex links,
full local source columns and complete source/reference pointer outputs.

## Exact contact reduction and its controlled leading term

On `N=2`, let `P` project onto all onsite pairs `|d_x>` and `Q=I-P` retain
every separated pair. With `T=dGamma(h)`, one-body hopping gives `PTP=0`.
Put `B=PTQ` and `D=QTQ`. For real `E>||D||`, solving the second block of
the eigenvalue equation gives

```
v = (E-D)^(-1) B^* u,
F(E) = U I_P + B(E-D)^(-1)B^*,
H (u,v) = E (u,v)  iff  F(E)u = E u.
```

Where both inverses exist, the exact projected resolvent is
`P(E-H)^(-1)P = (E-F(E))^(-1)`. No scattering mode is removed from `D`.
The bounded numerical routine solves this full resolvent; the block identity
does not require constructing a global exponential. Dense-column resolvent
execution is limited to 8,192 total two-particle states; larger admitted
fields retain the algebraic formulas and leading operator but are rejected
by this host evaluator before factorization.

For any link configuration a Schur bound is
`b=2 max_p sum_q |h_pq| <= 6c`, giving `||T||,||D||<=b`.
If `U>b`, min--max gives exactly `rank P` eigenvalues at least `U`, while
the remaining eigenvalues are at most `b`: the Rayleigh quotient on `P`
is exactly `U`, and `U P` has rank `rank P`. The pair sector is thus
isolated by at least `U-b`. This is a sufficient strong-coupling bound;
it is not asserted to be sharp. The sharper anticommuting-symbol estimate
for a flat field is not transferred to arbitrary flux.

The leading pair operator is

```
H_pair^(2) = U I_P + B B^*/U.
```

Two CAR hopping orders give, for a neighboring pair of sites,

```
(B B^*)_(x,y) = 2 det(h_xy),
(B B^*)_(x,x) = sum_(y adjacent x) Tr(h_xy^* h_xy).
```

For the off-diagonal identity, the two-particle exponential is the exterior
square of the one-particle exponential. Its pair-to-pair matrix element is
`det((exp(t h))_xy)=t^2 det(h_xy)+O(t^3)`; comparing the coefficient of
`t^2` with `exp(t T)` gives the factor two and its sign. On the diagonal,
the first hop separates the pair; summing the squared one-hop amplitudes
gives the displayed Frobenius norm. No intermediate separated state is
omitted in either calculation.

For the Weyl bond this is
`(H_pair^(2))_(x,y) = (c^2/(2U)) exp(2i A_xy)` and a diagonal shift
`degree(x)c^2/(2U)`. A direct full-16-state calculation with an unrelated
spin-mixing bond checks the determinant sign and the two virtual paths.
Under site phases the contact basis transforms by `exp(2i alpha_x)`.
Both the exact `F(E)` and its leading term therefore carry the charge-two
site representation; the leading nearest-neighbor link phase is doubled.

Let `C` be a Schur upper bound for `||B B^*||`. The resolvent identity gives
the explicit error

```
||F(E)-H_pair^(2)||
 <= C [ |1/E-1/U| + b/(E(E-b)) ],       E>b, U>0.
```

Indeed `||(E-D)^(-1)-I/E||<=b/(E(E-b))`, and the replacement of `E`
by `U` is separately charged. The bound applies uniformly to every contact
vector and untouched reference. At any actual pair eigenvalue it bounds
the distance to the spectrum of the leading operator. It does not make
the leading pair evolution the exact full history.

In particular, doubling a link phase does not imply that every exact
full-field response has half the single-fermion flux period. A square
control has identical leading pair matrices at fluxes differing by `pi`,
but different full-field spectra. The omitted separated-pair paths matter.

## The current needs the full-state norm

Set `X(E)=(E-D)^(-1)B^*`. A contact vector `u` reconstructs the full vector
`(u,Xu)`, whose norm is

```
u^* M(E) u,       M(E)=I+X(E)^*X(E)=I-partial_E F(E).
```

For a link parameter `A`, differentiate `F` at fixed `E`. If `J=partial_A H`,
then block differentiation gives

```
partial_A F = J_PP + J_PQ X + X^* J_QP + X^* J_QQ X.
```

It follows directly that the current expectation of the normalized full
eigenstate is

```
<J> = (u^* (partial_A F) u)/(u^* M u).
```

For a simple eigenvalue this is also its Hellmann--Feynman derivative
`dE/dA`. The formula remains an expectation identity for each reconstructed
eigenvector at degeneracy; following an individual degenerate eigenvalue
branch requires the usual diagonalization inside that eigenspace.
If `u` has contact norm one, its full eigenstate's contact probability is
`1/(u^* M u)`. Thus the same separated-pair contribution fixes both residue
and current normalization. Dividing only by the contact norm gives a wrong
current, as an independent full-eigensystem and finite-link-derivative
control demonstrates.

This attaches a concrete interaction-dependent gauge response to the same
finite pair. Selection of the link field and action, preparation and clock
calibration, a native localized protection mechanism and a controlled
physical continuum are additional obligations.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m model.gauged_field.certificate
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/gauged_field
```

The exact integer identities carry their stated finite conclusions.
Resolvent, eigensystem and source-contraction residuals are independent
floating diagnostics with explicit tolerances, not interval or hardware
certificates. No claim depends on postselection or numerical state
renormalization during the source evolution.
