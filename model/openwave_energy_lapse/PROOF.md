# Energy-lapse coupling and finite attractive response

## Complete native matter and the declared geometry law

Let M0 be the complete eight-dimensional native three-oscillator cutoff-two
Hamiltonian at geometry one, in the energy unit sqrt(14). Rebuilding the
original coefficient expressions gives

\[
 M_0=D/2-(XXI+IXX)/7,\qquad D_{bb}=2\operatorname{popcount}(b)+3.
\]

Its eight simple eigenvalues are enclosed by exact characteristic-polynomial
isolation and independent rational matrix inertia. Write E0 for the actual
lowest eigenvalue and K=M0-E0 I. The complete operator K is nonnegative;
no level is deleted. On a finite cubic graph retain one such matter factor
per site and a geometry register q_x on the dyadic grid of spacing 2^-40
within [-1/1000,1/1000]. Exterior geometry is zero. The specified Hamiltonian is

\[
 H(q)=\frac\kappa2q^TLq+\sum_x e^{q_x}K_x,\qquad \kappa=10^4,
\]

where L is the positive three-dimensional Dirichlet graph Laplacian. The
geometry register is diagonal and has no kinetic update. Every q configuration
and every matter eigenline belongs to the full carrier. This changes the
geometry coupling from the reciprocal scalar action. In particular the
q-dependent subtraction e^q E0, the spatial gluing, stiffness, massless elastic
coordinate and minimum-energy preparation service are declared architecture.
No equivalence to the unchanged reciprocal source is asserted.

For every admitted density matrix, including states entangled with a reference,

\[
 \partial_{q_x}H_x=H_x=e^{q_x}K_x\geq0.
\]

Thus the geometry response of a ground-referenced excitation is its actual
local energy. This is an operator identity on all eight levels, including
coherences. The identity does not infer this coupling from the native
reciprocal derivative. That derivative has different eigenvectors and the
retained positive-energy, opposite-charge counterexamples apply to it.

## Continuous minimum and its finite-carrier enclosure

Prepare native energy eigenlines at specified sites and the ground elsewhere.
Let E_x be their exact ground-relative source energies. Subtracting the fixed
baseline sum E_x gives the geometry energy

\[
 F(q;E)=\frac\kappa2q^TLq+\sum_x E_x(e^{q_x}-1).
\]

All E_x are nonnegative. The Hessian

\[
 A(q;E)=\kappa L+\operatorname{diag}(E_xe^{q_x})
 \succeq \mu I,\qquad \mu=12\kappa/(n+1)^2,
\]

is positive definite. The continuous extension has a unique global minimum
q*. The discrete maximum principle gives q*<=0; its stationarity equation
then gives -L^-1 E/kappa<=q*<=0. On the negative window, the map
q -> -L^-1(E exp(q))/kappa is a contraction whenever
||L^-1 E/kappa||_infinity<1. The executed sources satisfy the stronger bound
that this norm is below 1/1000.

This admissibility statement is rationally checked. The exact one-axis
barrier b_i=i(n+1-i)/2 obeys Lb>=1 on the cube, including the other-axis
boundary edges. Hence ||L^-1||_infinity<=(n+1)^2/8. Rounding the numerical
linear solution to a dyadic vector and bounding its residual with the exact
native spectral intervals encloses ||L^-1 E/kappa||_infinity. The floating
fixed-point iteration proposes a witness, while the bound certifies the
window and contraction independently of its convergence trace.

For the proposed dyadic q, every non-source residual kappa Lq is an exact
integer multiple of kappa2^-40. At the occupied sites, rational exponential
series and native spectral intervals enclose

\[
 r=\nabla F=\kappa Lq+E e^q.
\]

Strong convexity gives

\[
 F(q)-\frac{\|r\|^2}{2\mu}\leq \min_{\mathbb R^{n^3}}F
 \leq \min_{\rm dyadic\ carrier}F\leq F(q).
\]

The code uses n^3||r||_infinity^2 as a certified norm upper bound and
retains both energy endpoints. It therefore encloses the actual minimum of
the fully finite geometry carrier without enumerating its Cartesian product.
The native spectral and exponential uncertainties are included. A preparation
controller implementing this minimization is an explicit service; the displayed
diagonal source has no autonomous relaxation law.

## Attraction and its observation domain

The Hessian A is an irreducible positive M-matrix, so every entry of A^-1 is
strictly positive. For the continuous minimized energy W(E), the envelope
and implicit-function identities give

\[
 \frac{\partial^2 W}{\partial E_x\partial E_y}
 =-e^{q_x^*}(A^{-1})_{xy}e^{q_y^*}<0.
\]

Integrating this mixed derivative over two positive source strengths proves
negative mutual energy after subtracting both one-source self energies.
This derivative theorem concerns the continuous extension. Finite-grid
quantization can suppress changes below its resolution; the exact finite
pair differences are separately enclosed by the preceding certificate.

On the 31-cube, two first-excited native lines at separations 4,5,6,7,8 have
strictly negative finite-carrier mutual energies. The outward central
energy difference at separation 6 is strictly positive, giving attraction.
Its magnitude is approximately 1.53313e-7 in the stated source units; the
full rational endpoints appear in the receipt. Its ratio to
E1^2/(4 pi kappa6^2) lies between 1.05989 and1.05990. The lattice and boundary
shift is retained. This is the interaction after subtraction of the two
self energies, not a claim about total boundary forces.

At leading small source strength, stationarity gives q=-L^-1 E/kappa and
W=-E^TL^-1E/(2kappa). The long-distance continuum Green comparison therefore
has model coefficient G_model=1/(4 pi kappa). Stiffness and a permitted
positive geometry mass change this response; neither a numerical physical
Newton constant nor protection of infinite range is selected here.

The source positions are external preparation choices. The energy difference
is an operational static-force observation using complete energy readouts and
separation interventions, without a freely moving body or geodesic trajectory.
Matter spectra, geometry records and reference outputs retain all complementary
outcomes. The full source ground has every native factor in its ground and q=0;
prepared source energies, their work and the conditional geometry minimum are
distinct from that global ground.

The argument uses finite spectral isolation, the exact positive quadratic
residual bound and the discrete maximum principle. It supplies no covariant
stress dictionary, Einstein entropy-stationarity premise, physical length or
laboratory clock calibration.
