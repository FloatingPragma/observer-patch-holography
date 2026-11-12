# Separated clocks and finite scalar-wave observations

## Geometry and source ancestry

The geometry is a specified static dyadic preparation from the
[native-charge elastic extension](../openwave_native_force/PROOF.md).
It has 31 cells in each direction, zero exterior geometry displacement,
stiffness one and a central charge \(9/140000\). The sine-transform
solution of the declared quadratic energy is rounded to the same
\(2^{-48}\) geometry grid. Its exact residual and energy certificate
are regenerated. This is a retained variational preparation; it is not
asserted to be the exact minimizer of the full nonlinear energy or to
be prepared by an autonomous native controller. The geometry stays fixed
throughout the probes, with no probe backreaction claim.

Set \(r_x=1+q_x\). The reciprocal scalar coefficient density and its
metric dictionary are unchanged:

\[
C=\operatorname{diag}(-1/r,1/r,r^2,r^2),\qquad
g=\sqrt{-\det C}\,C^{-1}
=\operatorname{diag}(-r^2,r^2,1/r,1/r).
\]

The normalized local frame axes are a declared fixed spatial frame of
this extension; they are not a source-derived spherical coordinate
attachment. The lapse is \(r\), radial null coordinate speed is one,
and transverse null speed is \(r^{3/2}\). On the 31-cell transverse
central slice, a piecewise-constant metric ray has travel time

\[
T_{\rm ray}=\sum_j r_j^{-3/2},\qquad
T_{\rm ray}-31\in(0.0000641973322511292,
                         0.0000641973322511293).
\]

The receipt encloses this value with rational square-root intervals.
A radial metric ray has no corresponding delay. These are metric
readouts. Identifying them with the arrival time of the finite dispersive
field below would require a separate geometric-optics estimate.

## Complete native clocks do not share one lapse ratio

At the central cell and a cell eight edges away, retain the complete
eight-level native Hamiltonian \(H(r)=\sqrt{14}M(r)\), with

\[
M(r)=\left(\frac r4+\frac1{7r}+\frac{3r^2}{28}\right)D
-\frac{XXI+IXX}{7r},\qquad D_{mm}=2\operatorname{popcount}(m)+3.
\]

It is rebuilt from the original compressed oscillator coefficients.
All eight eigenvalues at each rational geometry are enclosed by exact
characteristic-polynomial isolation and independent rational inertia
counts. For each of the seven ground-to-excited transitions, compare
\(\omega_k(r_a)/\omega_k(r_b)\). The common factor \(\sqrt{14}\)
cancels. Every comparison interval is disjoint from \(r_a/r_b\), and
different transitions have disjoint frequency-ratio intervals. For
example, the first three ratios are approximately

\[
0.999984282048,\quad0.999990992272,\quad0.999994626553.
\]

The comparison is an operational finite-clock experiment: each selected
spectral preparation is \((|0\rangle+|k\rangle)/\sqrt2\), followed by
the unchanged source evolution for \(t=1/\sqrt{14}\). A Y-basis Ramsey
measurement gives plus-Y probability
\((1-\sin(E_k-E_0))/2\) and its complement, where the displayed \(E\)'s are eigenvalues
of \(M\). The six other local level outcomes have probability zero and
are retained. The two-clock carrier has all 64 levels. Exact interval
arithmetic and a bounded Taylor remainder enclose the probabilities.
Spectral preparation and the common phase apparatus are declared services.

The identity \(e^{-itH_r}=e^{-i(rt)(H_r/r)}\) therefore cannot by itself
establish universal gravitational redshift. In this source, proper-frame
transition frequencies depend on geometry. No clock is rescaled separately
to remove that dependence, and no energy exchange between clocks is
claimed. A physical clock standard and its universal metric coupling
are distinct requirements.

These oscillator devices are held in the fixed coordinate frame. Their
proper radial and transverse dimensions change with \(r\). The result
therefore rejects an uncorrected universal lapse reading from these
particular transitions; it does not refute the equivalence principle or
establish a comparison of identical devices held at fixed proper size.

## A declared spatial action with complete finite records

Use the same coefficient density for a scalar on the transverse slice:

\[
L=\frac12\sum_j\frac{\dot\phi_j^2}{r_j}
 -\frac12\sum_e a_e(D\phi)_e^2,\qquad
a_e=\frac{r_{e,-}^2+r_{e,+}^2}{2}.
\]

The scalar has zero Dirichlet boundary values and exterior geometry
\(r=1\). This spatial gluing and edge averaging are declared architecture
choices. It is a spatial extension of the same scalar density, with no
asserted implementation by the original uncoupled eight-level clock.
The canonical change \(\phi=R^{1/2}Q\), \(p=R^{-1/2}P\) gives

\[
H=\tfrac12 P^TP+\tfrac12Q^TKQ,\qquad
K=R^{1/2}D^TAD R^{1/2}.
\]

All 31 positive normal frequencies are retained. For uniform geometry,
\(K(r)=r^3K(1)\), with dispersion
\(\omega^2=4r^3\sin^2(k/2)\). Its low-momentum speed is the metric
transverse speed, while its finite group speed has the extra factor
\(\cos(k/2)\). The executed wave is not identified with a null ray.

Quantize the normal modes on the complete vacuum-plus-one-quantum
carrier, dimension 32, with diagonal Hamiltonian \(0,\omega_1,\ldots,
\omega_{31}\). The common oscillator zero-point constant is removed
explicitly; it changes no record. The finite cutoff theorem used in
[field observables](../openwave_field_observables/PROOF.md) gives
\(i[H,Q_k]=P_k\), \(i[H,P_k]=-\omega_k^2Q_k\) on the entire cutoff
carrier. It does not replace the compressed Hamiltonian by the square
of truncated coordinates.

Prepare an equal coherent superposition of vacuum and a normalized
right-moving Gaussian in the orthogonal one-quantum site frame. Evolve
with that Hamiltonian and measure vacuum plus all 31 site-frame
projectors. These are specified observer effects, without a physical
photon-position interpretation. All 32 probabilities, including the
vacuum weight one-half, sum to one. Against the flat source, the complete
detector-law total variations at source times 8, 12 and 16 are
approximately \(4.01534,7.40718,7.63711\) parts per million. The
conditional mean-position shifts are negative in all three experiments.

Evolution uses sixteen equal steps of the degree-36 exponential
polynomial. Since \(r_j\leq1\) and \(a_e\leq1\), the incidence norm
bound gives \(\|\sqrt K\|\leq2\). Using the larger bound three, let
\(x=3|t|/16\). The one-step operator error is at most

\[
\epsilon=\frac{x^{37}}{37!}\frac1{1-x/38},\qquad
\epsilon_{\rm total}\leq16\epsilon(1+\epsilon)^{15}.
\]

The source frequency decomposition independently evaluates each full
evolution. Its discrepancy, norm defect and square-root residual are
reported as numerical diagnostics. The analytic Taylor certificate
does not include floating-point rounding or eigensolver error; the
recorded numerical wave probabilities are not interval-certified.
The native clock and metric-ray intervals above use exact bounds.

## Testing the Einstein attachment

For the same metric at first order in \(q\), direct linearized Ricci
evaluation gives

\[
G_{00}=\partial_x^2q-\tfrac12\partial_y^2q
                    -\tfrac12\partial_z^2q.
\]

A harmonic exterior scalar profile does not make this vanish. For
\(q=-a/\sqrt{x^2+y^2+z^2}\), the value at \((3,0,0)\) is \(-a/9\).
On the actual dyadic geometry, the discrete Laplacian at \((21,15,15)\)
is \(2^{-47}\), while the discrete linearized \(G_{00}\) is
\(-2708207/35184372088832\). The pinned matter source vanishes there.
This challenges a pressureless, first-order Einstein interpretation of
this scalar geometry. A complete covariant stress dictionary is not
constructed; stress is not defined from the Einstein tensor to force
agreement. The smooth counterexample and finite-grid diagnostic retain
their distinct domains.

The [regional Einstein decoder](../matched_einstein_limit/PROOF.md)
requires matched source, clock, area, stress and curvature premises.
Those premises are not discharged by this calculation. The original
failure of universal positive-energy scalar coupling also applies.
The metric and finite-wave observations therefore support a controlled
geometry-dependent source response, not an Einstein or physical
gravitational-redshift theorem.
