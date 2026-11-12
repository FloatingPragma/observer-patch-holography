# Native clocks and wave records on one lapse profile

## One source and complete clock observations

Use the same complete native K=M0-E0 I, finite geometry register, stiffness,
boundary and prepared central source as the
[energy-lapse construction](../openwave_energy_lapse/PROOF.md).
The dyadic 31-cube geometry q is bound to its exact variational certificate.
It is prepared before the probe intervention. The full source contains the
geometry energy, all controlled native factors exp(q_x)K_x, and the finite
wave field specified below. Geometry has no kinetic update and stays constant
during each probe. Probe work is recorded; a self-consistent dynamical
backreaction claim is not made.

At two sites a and b, prepare (|0>+|k>)/sqrt(2), for each of the seven native
excited levels k. The exact native gap Delta_k is the same at both sites.
Evolution under the actual source Hamiltonians gives phases
omega_a=exp(q_a)Delta_k and omega_b=exp(q_b)Delta_k, in one common source
time. Therefore

\[
 \frac{\omega_a}{\omega_b}=e^{q_a-q_b}
\]

for every transition. This is measured with the two Y quadratures and all
six orthogonal native levels. At duration one the plus-Y probability is
(1-sin(omega))/2. The two full eight-output instruments form 64 complete
joint records, retained for each k. Characteristic-polynomial/inertia
spectral intervals, rational exponentials and a bounded sine series enclose
every probability. Preparation work is omega/2 for each clock; its mean
energy is conserved by the free source evolution. The mean energy change
of this particular Y readout is zero, although the general readout channel
does not preserve energy on arbitrary inputs.

The actual source zero E0, its q-dependent counterterm and the coherent
spectral preparation are explicit. A common controlled phase reference is
part of the apparatus. The calculation establishes universal simulation
clock ratios on this architecture, without a laboratory clock, a physical
proper-rod standard or an uncontrolled oscillator onset.

## Isotropic field gluing and its metric

Declare a scalar field with coefficient density

\[
 C=\operatorname{diag}(-e^{-4q},1,1,1).
\]

The algebraic metric dictionary g=sqrt(-det C) C^-1 gives

\[
 g=\operatorname{diag}(-e^{2q},e^{-2q},e^{-2q},e^{-2q}).
\]

Its lapse exp(q) agrees exactly with every native transition ratio. Its
proper spatial scale per coordinate length is exp(-q), and null coordinate
speed is exp(2q). The isotropic gluing is a declared choice. On the central
31-site slice, the piecewise metric ray time is sum_j exp(-2q_j), with a
strictly positive delay relative to 31. Rational exponential bounds enclose
the complete value. By cubic symmetry the same central ray result applies
to each coordinate axis. No claim equates this metric ray with a finite
lattice wave-packet arrival without a geometric-optics estimate.

The finite field action is

\[
 \tfrac12\sum_xe^{-4q_x}\dot\phi_x^2
 -\tfrac12\sum_{\rm edges}(D\phi)^2.
\]

Canonical mass weighting gives the positive frequency-square matrix
K_wave=exp(2q)D^TD exp(2q). Vacuum and every one-quantum normal mode form a
32-dimensional carrier with its exact normal-mode Hamiltonian. No mode of
that finite construction is omitted. At constant q, the exact dispersion is
omega^2=4exp(4q)sin^2(k/2). The field is a declared scalar extension rather
than an identification of its 32 states with the native eight-state matter.

Prepare an equal coherent superposition of vacuum and a normalized moving
one-quantum Gaussian. The detector retains vacuum and all 31 orthogonal
site-frame outcomes. Actual evolution at times 8,12,16 uses sixteen degree 36
Taylor steps and is checked against independent spectral evolution. Since
q<=0, the incidence norm gives ||sqrt(K_wave)||<=2. The conservative norm 3
Taylor remainder bounds the exact operator approximation. Floating-roundoff
and numerical diagonalization errors are reported as diagnostics, not
included in that analytic bound. The receipt retains source energy before
and after evolution, and the nonzero work of the site-frame measurement.
Deleting q gives a distinct complete probability law.

## First-order curvature and its boundary

Write the weak metric as h=diag(-2q,-2q,-2q,-2q). Independent differentiation
of the linearized Ricci tensor gives

\[
 G_{00}^{(1)}=2\Delta q,\qquad G_{0i}^{(1)}=0,
 \qquad G_{ij}^{(1)}=0.
\]

The equation kappa Lq+rho_q=0 follows from the energy minimum, where
rho_q is independently defined as the coordinate energy exp(q_x)E_x of
the prepared native matter. It is not defined from curvature. On the dyadic
witness, exact finite differences therefore give

\[
 G_{00}^{(1)}-2\rho_q/\kappa=-2r/\kappa,
 \qquad r=\kappa Lq+\rho_q.
\]

The inherited rational residual certificate bounds this equality at the
source and exterior points. It supplies the static first-order comparison
8 pi G_model=2/kappa in the stated coordinate-energy convention. Matter
coordinate energy versus a full covariant stress tensor, nonlinear metric
terms and geometry stress require separate treatment. No finite regional
entropy or area-stationarity premise is discharged by this comparison.

A control keeps the same lapse and clock/force source while choosing spatial
metric exp(-2 gamma q) with gamma=0. Direct curvature at q=-1/r gives
G11=2/27 at (3,0,0), while gamma=1 has no such first-order vacuum residual.
Clock and static-force consistency alone therefore does not select gamma=1.
The displayed scalar action makes that choice explicit, and no Einstein
field equation is installed as the microscopic source law.

The construction supplies complete source-bound clock and wave records,
positive ray delays and an action-derived first-order density comparison.
It supplies no physical Newton constant, laboratory unit system, freely
moving source body, full nonlinear Einstein dynamics or measured
physical gravitational-redshift prediction.
