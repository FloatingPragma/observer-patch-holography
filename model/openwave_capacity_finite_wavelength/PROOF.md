# Finite-wavelength cosmological records with an interval remainder

The source is the supplied flat Einstein/dust/radiation/canonical-capacity-
scalar model of `openwave_capacity_transfer`. Its Hubble normalization,
density fractions, scalar attachment and horizon choice retain their input
ancestry. This calculation concerns a specified preparation from log scale
factor N=-1 to N=0. It does not identify that preparation with the primordial
mode evolved from the radiation era.

## Comoving reduction

Use the parent's variables and write

    U=U_m, v=U_r-U_m, R=Phi+U,
    w=y-sU, j=z-s_N U-s[(h-1)U+R],
    Delta_m=delta_m+3U, Delta_r=delta_r+4U_r.

The complete seven-equation system is exactly equivalent to

    Delta_m'=-Q U+6 Omega_r v+(3/2)s w,
    Delta_r'=Delta_r-(6 Omega_m+2s²+4Q/3)v+2s w-4QU/3,
    v'=h v+Delta_r/4, w'=j,
    j'=-(3+h)j-(M-3s²/2+Q)w+6 Omega_r s v-QsU,
    U'=(h-1)U+R, R'=2 Omega_r v+s w/2.

Primes mean N derivatives. Symbolic substitution retains the scalar
background derivative, its potential gradient and the exact gravitational
constraint. For Q=0 the regular solution has the first five coordinates
zero, R=3/2 and U=T, with T given by the proved single positive quadrature.
The apparent k=0 perturbation is residual adiabatic coordinate freedom.

At N=-1 supply the finite-wavelength initial correction
Delta_m=-2Q Phi/(3 Omega_m), with the other comoving corrections zero.
This cancels the full Hamiltonian constraint exactly. The parent's
constraint propagation theorem carries it through the complete evolution.
This late preparation has an explicit density adjustment; it is not a
derived early-Universe state-selection law.

Set epsilon=(k/H0)² and use coordinates
Y=(Delta_m,Delta_r,4v,2w,j,U,R). The system has the form
Y'=(A0+epsilon B)Y. Its solution has

    Y=Y0+epsilon Y1+Z,
    Y1'=A0 Y1+B Y0,
    Z'=(A0+epsilon B)Z+epsilon² B Y1, Z(-1)=0.

Only the known scalar T enters B Y0. Appending 10T therefore gives an
eight-coordinate affine equation for Y1 and T; its matrix and its exact
background derivative are written directly in `matrices`. The coefficient
of the metric potential is R1-U1. The coefficient of comoving dust is the
first coordinate. These are distinct observables.

## Uniform norm and complete remainder

The infinity-norm logarithmic norm is the maximum of each diagonal entry
plus the sum of the absolute off-diagonal entries. On the entire past
branch, 0<=s²<=8/27, -2<=h<=0 and 2<=M<=6. The scalar restoring coefficient
M-3s²/2 is positive. In the displayed weighted coordinates these facts give
mu_infinity(A0)<=4. For example the radiation row is bounded by
1+3/2+4/27+11/20<4; the scalar-derivative row by -1+3+33/40<4.
Every other row is smaller. No physical sign is replaced before taking
the logarithmic norm.

Write f=Q/epsilon=a²(1-Omega_dark)/(Omega_m0 a+Omega_r0).
For a<=1, f<=1/Omega_m0. The supplied Omega_m0>5/18 gives
||B||_infinity<=5/(3 Omega_m0)<6. The exact regular solution has
||Y0||<=3/2, and the interval initial correction gives ||Y1(-1)||<8.
Thus a coarse independent bound is

    ||Z(0)|| <=102 epsilon² exp(4+6epsilon).

The reported remainder uses a sharper cellwise version. On each cell the
enclosed Y1 supplies an upper bound F_i on the complete forcing B Y1,
including its U, relative-velocity and comoving-scalar terms. A positive
recurrence

    z_(i+1)=exp[(4+6epsilon)h] (z_i+epsilon² h F_i)

bounds the exact variation-of-constants integral. The code also compares
it to the independent coarse bound. Metric and dust errors are at most
2z and z, respectively. There is no omitted k-dependent complement.

## Interval transport

The inherited monotone logit Taylor map gives the background at every
endpoint. Monotonicity encloses its whole cell. The coefficients and their
N derivatives use outward rational arithmetic, including square roots and
exponentials. For an exact centre c, enclose its local trajectory in
c+[-d,d]^8, where

    d=h [||A|| (||c||+h||b||) exp(h||A||)+||b||].

Taylor's integral remainder encloses the next point by

    c+h(A_start c+b)+(h²/2)[A' Y+A(A Y+b)].

The second derivative is evaluated on the whole trajectory tube and whole
background cell. Its interval midpoint is the next exact rational centre;
the interval radius is added to the existing error propagated by
exp(h mu_infinity(A_cell)). Rounding is outward after every operation.
The all-cell tube used in the k^4 forcing bound also retains this propagated
initial error. The background uncertainty, initial quadrature tail and
quadrature remainder are included. The 1024- and2048-step calculations are
independent enclosures with overlapping outputs and a smaller refined
coefficient radius.

The finite nonzero wave numbers are k/H0=1/100 and1/50. Numerical values
and outward displays are in `results.json`; no observational target enters
this selection. The supplied linear normalization is R=3/2, not a physical
primordial amplitude. A simultaneous H0->lambda H0, rho->lambda²rho,
L->L/lambda, k->lambda k leaves every dimensionless record unchanged.
There is no native capacity-to-horizon/clock map, absolute cosmological
normalization, early thermal population, Boltzmann hierarchy or CMB fit.
