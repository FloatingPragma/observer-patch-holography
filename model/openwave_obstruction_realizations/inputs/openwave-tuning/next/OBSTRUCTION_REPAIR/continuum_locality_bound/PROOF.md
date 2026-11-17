# An explicit full three-dimensional observer-locality bound

This is a root derivation for the massless flattened source in
`../gravity/spectral_flattening/`, using the complete continuum multiplier
and the observer construction in `../gravity/boost_orbit_observer/`.
It provides a conservative numerical constant for the existing qualitative
bound. It changes neither the source nor its observer map.

Set R=1 initially. For real p, put r=|p|, q=p dot p and
omega=sqrt(1+q^3). The source's exact unitary multiplier can be written

    V(p) = sqrt((omega+1)/(2 omega)) I
           + beta (alpha dot p) q / sqrt(2 omega (omega+1)).

This is the same cos(theta) I + beta(alpha dot p/r) sin(theta) with
theta=arctan(r^3)/2, including its regular extension at p=0.

## A complex directional disk with a uniform norm bound

Fix any real unit vector n. Extend p to p+w n, where w is complex and

    |w| <= epsilon(r) = max(1,r)/16.

All square roots below use their analytic principal branches. We show
that V is holomorphic on a neighbourhood of this closed disk and its
operator norm is less than 3 there. The Clifford matrices are Hermitian
with unit directional norm, so

    ||alpha dot (p+w n)|| <= r+|w|.

For r<=1, the change in q is at most 33/256. Since the real central q is
r^2<=1, the change in q^3 is at most

    (289/256)^3-1 < 1/2.

Thus Re(1+q^3)>1/2 throughout the disk. It follows that
|omega|>=1/sqrt(2), Re omega>=1/sqrt(2), and
|omega+1|>=1+1/sqrt(2). The scalar coefficient has magnitude less than
3/2. The numerator norm of the second term is at most

    (17/16)(289/256) = 4913/4096 < 3/2,

and its denominator magnitude is at least sqrt(1+sqrt(2))>1. Consequently
the total norm is less than 3. The real-part bounds also exclude every
square-root branch singularity in this disk.

For r>=1, |q-r^2|<=33 r^2/256. Thus

    (223/256) r^2 <= |q| <= (289/256) r^2,
    |arg q| <= arcsin(33/256) < pi/12.

In particular Re(q^3)>0, so |omega|>=1. Moreover

    |omega|^2 = |1+q^3|
              >= Re(q^3)
              >= |q|^3/sqrt(2)
              >= r^6/4.

The last inequality follows from (223/256)^6>1/8. The scalar coefficient
has magnitude at most 1. Since Re omega>0, |omega+1|>=|omega|, and the
second term has norm at most

    sqrt(2) (4913/4096) < 2.

Again the total is less than 3. The strictly positive real parts ensure
holomorphy of the required square roots. The argument is independent of
the direction n and includes all real momenta.

## Cauchy estimate and the full continuum tail

Cauchy's integral formula, applied to the matrix-valued holomorphic
function, gives

    ||(n dot grad)^6 V(p)||
        <= 6! * 3 / epsilon(r)^6
         = 2160 * 16^6 * min(1,r^-6).

Its three-dimensional integral therefore obeys

    M6 <= 4 pi * 2160 * 16^6
          * [integral_0^1 r^2 dr + integral_1^infinity r^-4 dr]
        = 5760 pi * 16^6.

Using the integration-by-parts and Schur argument in the boost-orbit
proof, separated spatial sets at distance d obey

    ||chi_A V_R chi_B|| <= C3 (R/d)^3,
    C3 <= 960 * 16^6 / pi < 320 * 16^6 = 5,368,709,120.

The same estimate holds for V_R dagger. This is a bound on the entire
continuum operator, including all ultraviolet modes. Distributional
terms supported at zero separation do not contribute between these sets.

For the boost-orbit apparatus with rapidity z and positive geometric
spacelike gap g as defined in the source proof,

    overlap <= 8 C3 (1+exp(3|z|)) (R/g)^3.

For example, if |z|<=log(2) and g>=100,000 R, the operator norm is less
than 386,547,056,640 / 10^15 = 0.00038654705664. This deliberately loose
bound establishes a quantitative macroscopic regime without dropping
bands or subtracting a finite-cutoff control. It does not assert that
this large separation is necessary, or that the bound is sharp.

## What this establishes

The observer-orbit construction has an explicit full-space approximate
locality tolerance at sufficiently large geometric separation. Exact
microscopic locality still fails. There is no physical identification of
R here, no measured G, no native observer execution, no curved transport,
and no Einstein stationarity law. The bound is a component repair of the
observer construction, not a solution of gravity.
