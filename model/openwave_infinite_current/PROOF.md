# Infinite-photon current and Larmor certificate

## Coordinate-matrix derivative bound

Keep the complete twenty-dimensional spin/spatial block, the same mass16,
Gaussian transverse profile, field coefficient1/8 and trap as the declared
infinite-photon source. In photon coordinates x in R³,

    B(x)=B0+(sqrt2/8)sum_i M_i x_i,
    K(b,x)=sqrt(16²+(B(x)+bE_z)²).

The global B is unbounded. Every coordinate fibre is a20-by20 Hermitian
matrix. Exact LDL factorizations give ||E_z||,||M_i||<4/5, and
tr(E_z²)=5. The scalar f(y)=sqrt(16²+y²) has |f'|<=1. In a spectral basis
of B, the Frechet derivative has entries f[1](lambda_i,lambda_j)E_ij;
all divided differences have magnitude at most one. Consequently

    ||K'(b,x)|| <= ||K'(b,x)||HS <= ||E_z||HS=sqrt5

uniformly in x and b. Thus K(b)-K(0) is bounded on the full field space by
sqrt5 |b|. Differentiating K²=16²+(B+bE_z)² twice and using the Sylvester
inverse for K>=16 gives

    ||K''|| <= (||E_z||²+||K'||²)/16 <=141/400.

The domination is uniform, so these fibre derivatives define the derivatives
of the full multiplication operator with its common unbounded domain. No
finite-cap global norm is used for the infinite field.

## Global differentiated polynomial remainder

Let p6(t) be the sixth Taylor polynomial of sqrt(16²+t) at zero. The scalar
positive remainder and its derivative obey, for t>=0,

    0<=sqrt(16²+t)-p6(t)<=33 t^7/(2048*16^13),
    0<=d/dt[sqrt(16²+t)-p6(t)]<=7*33 t^6/(2048*16^13).

Both follow from the integral remainder and the decreasing positive next
derivative. For g(y)=sqrt(16²+y²)-p6(y²), this implies
|g'(y)|<=14*33 |y|^13/(2048*16^13). Its divided differences on each fibre
are bounded by the same coefficient times ||B(x)||^13. The Hilbert-Schmidt
argument therefore bounds the derivative error in a normalized trial phi by

    |<phi,Dg(B)[E_z]phi>|
       <=14*33 sqrt5 <||B(x)||^13>/(2048*16^13).

The exact bounds on B0 and M_i give ||B(x)||<=8/5+(sqrt6/10)|x|. For a trial
supported in total photon number at most C,

    <|x|^26> <= product_(j=0)^12 [2(C+j)+3].

Indeed for a vector of total number at most n,
sum_i||x_i v||²<=2n+3, and every coordinate multiplication increases the
number support by at most one. Apply this bound to every ordered monomial
and sum. Minkowski's inequality for the scalar L26 norm then gives the
explicit remainder used by the certificate. The finite support is a property
of the trial only; the actual dressed ground has the full photon carrier.

## Exact line selection and enclosed current

A cap-four numerical ground vector supplies integer dyadic trial coordinates.
The exact spectral polynomial in total Jz projects them onto+1/2, using all
half-integers between−13/2 and+13/2. Interval normalization defines an exact
unit vector; the eigenvalue calculation itself is not assumed exact.
Fourteen additional photon degrees evaluate B^k phi exactly for k<=14.
The inherited global square-root remainder, positive excited floor17 and
Temple inequality enclose the actual infinite ground energy and the trace
distance between the trial and its exact+1/2 line.

The polynomial current is evaluated with every ordered Frechet term:

    sum_(k=1)^6 c_k sum_(j=0)^(2k-1)
        <B^j phi,E_z B^(2k-1-j)phi>.

The same sparse outward arithmetic contract as `openwave_infinite_photon`
propagates exact algebraic entry and trial intervals. The global derivative
remainder above encloses the true trial current. The exact current operator
has norm at most sqrt5, so its ground-line expectation differs by at most
sqrt5 times the certified trace distance. Combining these bounds gives a
literal infinite-photon current interval and g=4E0|mu| in the same energy,
current and clock conventions. The spatial regulator is fixed throughout.

## Full-reference finite-time precession

Rotation covariance makes P K'(0) P=mu sigma_z on the exact dressed spin-half
band. For gamma=17-E0_upper, set

    v=sqrt5 |b|, r=(141/800)b²,
    D=r+v²/(gamma-2v),  gamma>2v.

The exact Jz sectors, perturbation theory and the same Temple floor give
|E_s(b)-E0-sbmu|<=D. Phase-aligned perturbed vectors differ from their
unperturbed lines by at most sqrt2 v/(gamma-2v). Therefore

    ||U_b(t)P-exp[-it(E0+bmu sigma_z)]P||
       <=2sqrt2 v/(gamma-2v)+|t|D.

Using the displayed current centre adds |bt| times its certified error.
Twice this norm bounds output trace distance with arbitrary entangled
spectator references. Every photon outcome, retained filter failure and
complete subsequent instrument is covered by channel contraction. The rows
are executable bounds on the infinite-source dynamics, not a numerical
execution of infinitely many photon outcomes. No physical Tesla conversion,
anomalous electron moment or trap-free magnetic coefficient is claimed.
