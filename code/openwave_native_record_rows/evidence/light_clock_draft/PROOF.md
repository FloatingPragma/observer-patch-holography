# A common native clock and quantitative light-propagation route

This is a new theorem/construction draft, with no new source execution.
It uses the GR7 scalar controller and compatible native Maxwell flux
source on the SAME static chart and coefficient u=1+h. It targets the
relative clock, slow-body, delay and deflection predictions with one
capacity coefficient. Identifying that capacity with G times an
independently measured source mass is a separate strength join; it is
not silently required to test every light/clock relation below.

## 1. Existing source attachments and their exact scope

GR7 compiles temporal/internal/spatial scalar weights u^4/16,u^2/16,1/16
using four plus-hit slots. The fixed field normalization Y/4 yields the
scalar density -u^4 partial_t^2-Delta+lambda u^2. Its Lorentzian metric is

    ds^2=-u^-2 dt^2+u^2 dx^2.                             (1)

The compatible-flux theorem supplies a full-moment Gaussian Maxwell
limit with electric pairing M_xy=u_xu_y and magnetic pairing tending to
u^-2. Its positive quotient transfer and finite-spacing CCR are proved
in native_compatible_flux_research/MAXWELL_ATTACHMENT.md. These are new
source-defined weights, not the supplied forms of the earlier Lean
WhitneyMaxwellDynamics module. That module's exact gauge-null and energy
algebra is reused; its integrated-prism time step is not substituted for
the source's nearest-time OS transfer.

SeamMaxwellContinuum proves useful uniform flat-symbol bounds but supplies
its chart and reversible evolution. Its translation-invariant1/20 bound
is not simply asserted for variable u. CommonWorldMaxwellClockJoin is
explicitly a same-index product without a physical calibration relation.
The present join instead uses the same native u, chart coordinate step a,
scalar density and Maxwell electric/magnetic precision. The temporal
coordinate is the reconstructed field/clock coordinate; the much longer
native sampling/reset resource clock stays separately counted.

The ideal static continuum Maxwell density is

    (1/2)[u^2 |E|^2-u^-2 |B|^2].                          (2)

It is Maxwell on(1), with canonical field amplitude fixed by sqrt(8/n)
of the native record sums. No later amplitude or coordinate-time fit is
permitted to make its rays agree with the scalar clocks.

## 2. Source clock and slow-body calibration

Let m=sqrt(lambda)>0 be any native internal scalar band. Freezing u in a
small detector region gives local source mass shell

    omega^2=|k|^2/u^4+m^2/u^2.

Thus every band's rest phase has coordinate rate m/u. Proper time inferred
from the same rest-clock family is d tau=dt/u. For the source OS time step a,
the finite rest rate is

    omega_a(u)=(2/a) asinh(a m/(2u)),                      (3)

with 0<=m/u-omega_a(u)<=a^2 m^3/(24u^3). This follows from
0<=x-asinh x<=x^3/6. Clock ratios retain the denominator correction in(3);
no exact finite-a equality to1/u is asserted. Spatial variation across a
clock packet and the packet preparation/readout error are additional
controlled terms. Native scalar WKB phases obey the massive Hamiltonian

    H(x,p)=sqrt(|p|^2/u^4+m^2/u^2).

At coordinate rest p=0 their ray acceleration is exactly

    d^2x/dt^2=grad u/u^5.                                (4)

Indeed p_dot=m grad u/u^2 and partial_p^2 H=I/(m u^3) at rest. The band
mass cancels. This is a reconstructed massive-packet statement with its
WKB detector error, not a claim that an individual native random-walk
path is a timelike trajectory.

For stationary emitter and receiver, conservation of coordinate frequency
and the same clocks give

    omega_receiver/omega_emitter = u_receiver/u_emitter. (5)

Here the frequencies are measured per local proper time. Finite spacing,
localization and detector errors must accompany(5) in an executed test.

## 3. Exact continuum Maxwell energy estimate

Set n=u^2 and H=B/n. In a source-free region, the continuum equations are

 n partial_t E=curl H,       n partial_t H=-curl E,
 div(nE)=div(nH)=0.                                       (6)

Their energy is (1/2) integral n(|E|^2+|H|^2). Use a finite conducting box
whose boundary lies outside the signal's domain of dependence, or whole
space with smooth positive bounded n and decaying fields. On the latter,
Psi=E+iH obeys

    i n partial_t Psi=curl Psi.                           (7)

The curl realization is self-adjoint in the n-weighted Hilbert space;
this also follows directly by integration by parts and the zero boundary
energy flux. If an approximate field obeys(7) with residual R, its exact
solution with initial mismatch e0 satisfies

 ||Psi(t)-Psi_app(t)||_n
    <=||e0||_n+integral_0^t ||n^-1/2 R(s)||_2 ds.          (8)

This is the unitary variation-of-constants estimate. It remains valid
without a geometric-optics approximation and gives a directly computable
error from a chosen finite wavepacket ansatz. Electric/magnetic boundary
conditions must be imposed in the real two-field formulation; a conducting
wall is not silently treated as a helicity-invariant boundary for Psi.

## 4. A constructive WKB packet, rather than a named ray assumption

Choose a real smooth eikonal phi on a time tube without caustics, with
k=grad phi nonzero and positive omega=-phi_t. Put

 Psi_app=exp(i phi/epsilon)(a0+epsilon a1).

The principal matrix of(7) is B=n omega I-i k cross. It has a one-dimensional
chosen-helicity kernel when n omega=|k|. Its other eigenvalues are |k| and
2|k|, so its inverse on the orthogonal complement has norm<=1/|k|.
Choose a0 in that kernel and solve the projected transport equation

    P (i n partial_t-curl) a0=0.

It is a linear first-order equation along the smooth eikonal rays, with
known coefficients from n and phi. With P a1=0 choose

    a1=-B_perp^-1 (I-P)(i n partial_t-curl)a0.             (9)

The order1 term then vanishes, and the exact residual is

    R=epsilon exp(i phi/epsilon)
                     (i n partial_t-curl)a1.            (10)

All constants in(8) are thus derivative norms of explicitly constructed
phi,a0,a1 on a compact tube. No unproved uniform WKB remainder is needed.
For fixed smooth envelope and tube, the error is O(epsilon), with its
actual constant displayed by(8)-(10). A shrinking beam width changes
those derivative norms and must be inserted before taking limits.

An exact Gauss initial field is available without assuming that a0 alone
has zero divergence. Choose g=i k cross(n a0)/|k|^2 at the initial time.
Then i k cross g=n a0, and

    n Psi_initial=epsilon curl[exp(i phi/epsilon)g]

has zero divergence exactly. Its difference from the WKB initial datum
is epsilon times the explicit smooth quantity curl g-n a1. Real and
imaginary parts give the two real initial fields. This supplies the e0
in(8) and a source-admissible transverse coherent displacement target.
The weighted Gaussian source's finite Weyl/OS preparation must be compiled
or reconstructed with its own error; it is not replaced by a new unrelated
wave source. GR5's finite Gaussian exponential/phase preparation strategy
extends by using the actual weighted covariance, with all normalization
and phase costs retained.

The eikonal Hamiltonian is omega=|k|/n. Therefore

 x_dot=k/(n|k|),     k_dot=|k| grad n/n^2.                 (11)

These are null rays of(1), with optical index n=u^2. The equality of the
Maxwell and scalar characteristic speeds u^-2 is an operator consequence
of the same native u, not an independent clock calibration.

For a bounded detector window D at fixed source time, use the positive
local energy readout W_D(Psi)=||1_D Psi||_n^2. If the norm error in(8) is d,
then

 |W_D(Psi)-W_D(Psi_app)|<=d(2||Psi_app||_n+d).              (12)

If the ansatz is supported in a chosen ray tube, the exact energy outside
that tube is at most d^2. A detector arrival/angle interval follows from
separated spatial/time windows and an envelope whose intended signal has
a declared positive energy margin. The window resolution and pulse width
must be charged; a peak location cannot be inferred from L2 convergence
without such a margin.

## 5. The finite source OS correction is quantitatively separate

Let L_a=M_a^-1/2 K_a M_a^-1/2 /a^2 be the positive physical spatial operator
on the finite transverse quotient. The actual source transfer frequency is

 E_a=(2/a) asinh(a sqrt(L_a)/2),
 Omega_a=sqrt[L_a(I+a^2 L_a/4)].                           (13)

On the exact source spectral band L_a<=K^2,

 ||E_a-sqrt(L_a)||<=a^2 K^3/24,
 ||exp(-it E_a)-exp(-it sqrt(L_a))||<=|t|a^2 K^3/24,
 ||sqrt(L_a)/Omega_a-I||<=a^2 K^2/8.                      (14)

The first two bounds use functional calculus and Duhamel; the last uses
1-(1+x)^-1/2<=x/2. The last term is the plaquette-electric versus canonical
momentum correction. It cannot be dropped at finite spacing. A component
outside the chosen band contributes its explicitly measured norm tail;
the spectral projector is a reconstructed functional of the SAME native
weighted operator, not an imported flat Fourier cutoff.

For the spatial continuum join, sample a fixed smooth continuum packet
with the stated edge/face geometry and compute its residual in the actual
discrete Maxwell equations. The discrete positive energy gives the same
estimate as(8), with the discrete residual norm. Midpoint Yee differences
have O(a^2) consistency for sufficiently smooth fields and coefficients;
a corner-sampled magnetic weight has only O(a) consistency. The actual
sampling rule must be declared. Source-hitting coefficient error enters
that residual through the explicit M,beta differences. In particular,
flat SeamMaxwellContinuum constants are not imported into a variable-u
stencil. A rigorous execution can bound the residual directly by interval
Taylor derivatives of the selected packet and u, without a general
spectral-convergence theorem.

For a practical joint limit, first fix a smooth WKB tube and a desired
error, then choose epsilon, detector width/margin, band K, grid a, source
Gaussian-word length and native caps so that each displayed term is below
its allocated error. Since WKB derivatives grow with1/epsilon and inverse
beam width, a must be chosen AFTER those quantities. Compatibility
acceptance and source moment-error amplification are exactly those of the
conditional native theorem, with a new rank/error budget on every complex.
This paragraph gives an ordered construction; it is not a completed finite
residual certificate or an executed beam measurement.

## 6. One capacity coefficient predicts acceleration, delay and bending

In a radial exterior take u(R)=1+C/R, C>0. The coefficient C is the same
harmonic hitting capacity that appears in the scalar controller and the
Maxwell precision. A finite escape boundary, lattice target shape and
nonradial multipoles require their own bounds before this ideal radial
profile applies. Neither C nor its relation to source mass is fitted from
the light output.

Equation(4) gives the exact coordinate-rest inward acceleration magnitude

    a_rest(R)=C/[R^2(1+C/R)^5].                           (15)

For light launched parallel to the z axis at y=b>0 from z=-L, let
alpha=C/b, Lambda=L/b>=1, Y=y/b, s=z/b and q=dY/ds. Fermat's equation for
the same optical index n=(1+alpha/sqrt(Y^2+s^2))^2 is

 Y'=q,
 q'=(1+q^2)[partial_Y log n-q partial_s log n],
 Y(-Lambda)=1, q(-Lambda)=0.                              (16)

This equation is obtained by varying integral n sqrt(1+q^2) ds. It is
independent of a separately selected light speed.

Here are deliberately conservative explicit weak-field bounds. Assume
alpha<=1 and alpha Lambda^2<=1/240. On |Y-1|<=1/2, |q|<=1/2, R>=1/2,
|grad log n|<=8alpha and ||Hess log n||<=32alpha. The right-hand side F of
q' satisfies |F|<=15alpha, |F_Y|<=60alpha, |F_q|<=22alpha. A bootstrap gives

 |q|<=30alpha Lambda,        |Y-1|<=60alpha Lambda^2,

which stay inside the stated tube. At Y=1,q=0 the difference from the
linear forcing -2alpha/(1+s^2)^(3/2) is at most2alpha^2. Thus

 |q(Lambda)+4alpha Lambda/sqrt(1+Lambda^2)|
   <=alpha^2 [7200Lambda^3+1320Lambda^2+4Lambda].          (17)

The outgoing angle is atan q, so adding22alpha^2 bounds its additional
small-angle error. Replacing the finite-distance leading angle by the
asymptotic4C/b costs at most2alpha/Lambda^2. These explicit constants are
large but finite; they are not a numerical fit. A sequence with
Lambda=alpha^-1/8 and alpha->0 satisfies the bootstrap and makes the error
relative to4alpha tend to zero.

The coordinate travel time to the output plane is

    t=b integral_-Lambda^Lambda n(Y,s) sqrt(1+q^2) ds.

Using the same tube bounds, the leading excess over2L is

    Delta t=4C asinh(L/b)+error,
    |error|<=10000 b alpha^2 (Lambda^3+Lambda).            (18)

For example, the path-length term is bounded using sqrt(1+q^2)-1<=q^2/2
and n<=9, while |1/R(Y,s)-1/R(1,s)|<=4|Y-1|. Their integrated bounds are
below the displayed constant. Unequal input/output distances give the
leading2C[asinh(L_minus/b)+asinh(L_plus/b)], with the symmetric bound
applied using the larger distance. A stationary detector's measured proper
time is dt/u_detector. Convert using its independently calibrated source
clock(3), retaining that clock's finite-band/localization error.

The finite-radius acceleration(15), clock ratio(5), delay(18) and
bending(17) consequently share one independently read capacity C. Their
leading relative coefficients are1,2 and4 in the usual acceleration,
Shapiro and deflection conventions. This tests the observer model's
relative gravitational predictions even before the separate C=G M_source
strength normalization is established.

## 7. Requirements for an acceptable source execution

The next packet must freeze a finite common target/escape chart, its
native hit/coefficient error certificate, the electric/magnetic sampling
rule, packet and detector windows, internal clock bands, scalar and
Maxwell source readouts, all word/cap schedules and all numerical error
allocations. An independently checked residual must join the finite
native operator to the specified WKB packet. A ray ODE simulation alone
would test only(16), not the native Maxwell propagation claim.

The positive weighted Gaussian field and its finite record characteristic
source provide the relevant starting source. An executed comparison must
retain its acceptance/preparation and negative-phase costs, not install
a Maxwell solver as if its fields were already native data. It must also
avoid confusing the source sampling resource clock with the reconstructed
proper clock. Those are concrete remaining work items for the full row27
packet. Source inertial energy/G and general Einstein dynamics stay
separate, rather than being hidden inside the relative light/clock test.

## 8. No-caustic and detector construction on the weak ray bundle

The smooth tube in section4 can be constructed for the radial weak-field
experiment, rather than supplied without a test. Use a parallel incoming
bundle whose transverse initial points lie within b/8 of the central
impact b, and strengthen the weak condition to alpha Lambda^2<=1/10000.
The vector version of(16) has the same derivative bounds on this bundle:
||F_Y||<=60alpha and ||F_q||<=22alpha while its distance from the target
is at least b/2. The previous displacement bootstrap keeps it in that
region.

Let J be the derivative of transverse ray position with respect to its
incoming transverse position, and K=J'. Then

 J(-Lambda)=I, K(-Lambda)=0,
 K'=A(s)J+B(s)K,  ||A||<=60alpha, ||B||<=22alpha.

On a bootstrap interval ||J||<=2, variation of constants gives
||K||<=240alpha Lambda exp(44alpha Lambda) and therefore

 ||J-I||<=480alpha Lambda^2 exp(44alpha Lambda)<1/2.       (19)

The strict improvement closes the bootstrap, so J stays invertible.
The longitudinal derivative also stays positive because |q|<1/2.
The bundle has no caustic before the output plane. A smooth optical phase
is obtained by solving |grad S|=n with constant incoming-plane S and
positive outgoing longitudinal derivative; the characteristic solution
is single-valued on this invertible tube. Its phase phi=S-t satisfies
the eikonal relation used above. The helicity projector is smooth there,
since k never vanishes and stays in a fixed direction chart.

Choose smooth compact transverse and pulse envelopes on that incoming
plane, a smooth transverse helicity, and solve the linear transport in
section4 along this explicit flow. The invertible flow and bounded smooth
coefficient derivatives give finite computable norms for a0,a1 and their
residual. Those norms can be bounded by interval ODE/Taylor calculations
before source execution. There is no assertion that their constants are
independent of a subsequently narrowed detector pulse or beam.

Arrival is read by an array of positive energy windows near the output
plane at declared source-clock times. Use a smooth pulse whose transported
energy has a fixed positive margin in its predicted arrival window and
zero ansatz energy in disjoint earlier/later windows. Equation(12) then
turns an explicit sufficiently small residual into a separated signal
interval. Compare the same windows and clock convention in the u=1
changed-controller reference. The pulse width, window spacing, detector
clock error and norm-error margin are included in the delay interval.

For angle, use two downstream position-window arrays separated by a fixed
coordinate distance D. On independently prepared copies of the same beam,
their occupied windows bound the displacement direction. The envelope
width contributes its geometric angular resolution, of order width/D;
ongoing curvature between the arrays contributes the integral of the
known ray force, rather than being called zero. The source norm error is
converted to an energy-window confidence margin by(12). This avoids
claiming a peak or direction from unqualified norm convergence.

To resolve a weak angle of order C/b and a delay of order C, first fix
window/pulse widths below those target separations; next choose the WKB
wavelength so that(8)-(12) preserve the signal margin; then choose a and
the native coefficient/word/cap budgets so that(C6),(C8),(C10) and(14)
fit the remaining margin. Since the macroscopic target size is resolved,
there is no constraint forcing the signal itself to shrink at the same
rate as a. This gives a constructive order of finite choices for a
relative light/clock test; its actual chosen values and costs still belong
in the next frozen executable packet.
