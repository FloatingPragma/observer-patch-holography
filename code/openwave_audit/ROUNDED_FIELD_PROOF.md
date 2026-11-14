# Native rounding and a controlled classical field

Let A be the actual integer edge-by-face boundary matrix. Fix nonnegative
integers a,b, alpha=2^-a, beta=2^-b and fresh independent fair native coins
at every rounding operation and coordinate. The declared synchronous source is

    E' = E + A h^a(B),     B' = B - A^T h^b(E').

Its mean is the linear shear M(E,B)=(E+alpha A B,
B-beta A^T(E+alpha A B)). This construction is Fable's reversible field
amendment. Its schedule, field memory and stencil are declared; this proof
does not derive them from the unchanged one-load process.

## Exact rounding noise

Repeated native halving of an integer x is an unbiased random rounding of
x/2^a to its adjacent integers. Induction through ceiling/floor halving
proves both support and mean. Its variance is r(1-r)<=1/4, where r is the
fractional part. For a=0 it has zero variance.

Write xi=h^a(B)-alpha B and eta=h^b(E')-beta E'. Conditional on the old
state, xi has mean zero and independent coordinates. Conditional on the
old state and the first kick, eta has mean zero and independent coordinates.
In particular every cross moment E[xi_f eta_e] vanishes even though E'
depends on xi. The complete update is

    z' = M z + L xi + H eta,
    L = [A; -beta A^T A],       H = [0; -A^T].

The two Gauss constraints are unchanged because each increment lies in an
actual boundary or coboundary image. Recorded coins give the exact inverse.

## Shadow energy and its error

For P with blocks beta I, alpha I and off-diagonal alpha beta A/2, let

    Q(E,B) = beta |E|^2 + alpha |B|^2 + alpha beta E^T A B.

Block multiplication gives M^T P M=P, L^T P L=beta A^T A and
H^T P H=alpha A A^T. Q is positive definite when
alpha beta ||A||^2<4, including its harmonic and constraint sectors.
Conditional means remove the old-state and xi/eta cross terms. Thus the
one-step conditional increase is exactly

    beta sum_f ||A_f||^2 Var(xi_f | old)
      + alpha E[sum_e ||A_e||^2 Var(eta_e | first kick) | old],

and is bounded by

    C = ||A||_F^2 [beta 1_(a>0) + alpha 1_(b>0)]/4.

For deterministic initial integer z0, put delta_t=z_t-M^t z0. Induction
gives E Q(delta_t)<=t C. No state-dependent variance is replaced by an
independent unconditional variance in this argument.

The process Y_t=M^-t delta_t is a martingale in the positive P inner
product. M is an isometry for this norm. The nonnegative submartingale
Q(Y_t) therefore obeys

    Pr(max_(t<=T) Q(delta_t/R) >= epsilon^2)
        <= min(1, T C/(R^2 epsilon^2)).

For the average of K independent replicas with the same deterministic
initial state, replace R^2 by K R^2. Independence is an explicit premise;
shared coins need a covariance bound and do not receive the factor K.
The replicas are distinct experiments, not an interacting shared field.

## Spatial refinement

For the actual LOCAL12 complex, the frozen integer positive-weight
certificate gives ||A||^2<=1504/17 for every periodic size. At a=2,b=3,
alpha beta=1/32 and alpha beta ||A||^2/4<=47/68<25/36. Hence

    Q(E,B) >= (beta |E|^2 + alpha |B|^2)/6.

This is a uniform conditioning bound. It does not use a numerical largest
eigenvalue. The elementary per-face bound ||A_f||_1<=16 and 26 faces per
cell gives ||A||_F^2<=6656 L^3, even allowing repeated edge traversals.
The receipt computes the sharper exact trace on sizes 2,3,4.

For h=1/L on a fixed physical periodic box, let the observed field be the
integer field divided by R. Multiply the diagonal norm by h^3. The maximal
noise probability is bounded by a constant times T/(K R^2 epsilon^2),
uniformly in L. With physical elapsed time t and T<=ceil(t/h), noise
vanishes in probability if h K R^2 tends to infinity. Deterministic
initial discretization error and deterministic evolution error remain
separate terms in the triangle inequality.

Consequently any controlled deterministic Maxwell convergence theorem on
the same retained sector and observation map transfers to these native
records under that scaling. The two acoustic branches have leading speed
sqrt(2 alpha beta/3)=1/sqrt(48) cells per macro, equal to the scalar
LOCAL12 walk's leading speed. This statement neither removes optical
branches nor proves a deterministic spatial approximation error by itself.
Native service duration and the chosen macro clock remain distinct.

Amplitude scaling changes the admitted integer preparation; a charge
scaled by R must be declared if a fixed rescaled charge is intended.
Replica averaging measures an ensemble field. Neither operation makes
the unit-amplitude single-run field noiseless or derives matter feedback,
the electromagnetic strength, or quantum particle physics.
