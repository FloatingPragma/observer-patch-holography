# Finite static Maxwell reconstruction from compatible native flux records

This is a theorem-only attachment to PROOF.md. It introduces no source
execution or row award. It supplies the finite reflection-positive gauge
quotient and a boundary mean join for the proposed native compatible-flux
apparatus. A growing-complex continuum limit, source stress measurement,
target inertial energy, G and general gravitational backreaction remain
separate obligations. All gates, compatibility constraints, boundaries
and amplitudes below are declared apparatus, not consequences uniquely
selected by the endpoint law.

## 1. A finite spatial complex with no hidden zero-mode subtraction

Take a finite cubical three-ball X with its entire outer boundary as the
relative subcomplex. Its real relative cochains C0,C1,C2 have coboundaries
G:C0->C1 and C:C1->C2, with CG=0. Relative C0 potentials vanish on the
boundary and relative C1 potentials have zero tangential boundary value.
The integer incidence matrices are fixed geometric apparatus. Since
H1(X,boundary X;R)=0, ker C=im G. G is injective: a relative function with
zero gradient is constant and its boundary value is zero. These facts
can also be checked by exact ranks on each finite complex. No torus mode
or global-time zero-frequency component is subtracted.

Choose positive diagonal, time-independent electric and magnetic pairings
M on C1 and beta on C2, and set K=C^T beta C. Let

    T={q in C1: G^T M q=0}.

Then C1=im G direct-sum T in the M inner product. K is positive definite
on T, since its only null vectors belong to im G. These are finite
identities on the supplied complex, and reuse the exact gauge-null and
positive-Hodge structure of WhitneyMaxwellDynamics. The time discretization
below is the nearest-time Gaussian record action; it is not asserted to
be the different integrated-prism recurrence of that Lean module.

For one time link, use q_t in C1 and scalar potential s_t in C0. Its
oriented electric cochain and slice magnetic cochain are

    e_t=q_(t+1)-q_t-G s_t,       b_t=C q_t.                 (A1)

On the spacetime product slab these are exactly F=D a. Conversely every
compatible flux in im D has such a representation. Gauge changes are
q_t ->q_t+G chi_t, s_t->s_t+chi_(t+1)-chi_t. Thus the flux Gaussian
in PROOF equation(7) is the quotient of the Gaussian in these variables,
with only constant linear Jacobians. All magnetic slice weights, including the two ends, can stay equal to
beta. A product of the transfer kernels below contributes half a magnetic
slice at each end; multiplying each boundary vector by exp(-q^T K q/4)
restores the full endpoint weight. Thus no doubled endpoint ternary
probability, outside the stated p<=1/2 range, is required.

## 2. Positive transfer, including the integrated scalar potential

Fix the M-orthogonal gauge q_t in T. Integrating s_t yields

 integral_C0 exp[-(Delta q-Gs)^T M(Delta q-Gs)/2] ds
    = constant * exp[-Delta q^T M Delta q/2],              (A2)

because G^T M Delta q=0. The constant is finite, independent of q, and
contains det(G^T M G)^(-1/2). It cancels in normalized flux probabilities.
Do not mistake this integration for deleting electric-gradient flux:
those variables have their own Gaussian contact correlations on each
link. They are retained in the original compatible-flux measure.

On L2(T,dq), define

 T(q',q)=c exp[-q'^T K q'/4]
             exp[-(q'-q)^T M(q'-q)/2]
             exp[-q^T K q/4].                             (A3)

Here M and K denote their restricted quadratic forms, using any fixed
linear coordinates on T; c is a positive normalization constant. The
middle convolution has Fourier multiplier a positive constant times
exp[-k^T M_T^(-1)k/2]. The two exterior factors are the same bounded
positive multiplication operator. Therefore T is positive self-adjoint.
K_T>0 also makes T compact with a strictly positive Gaussian ground state.
This is a constructive positive transfer, not an inference from the
positivity of a Euclidean quadratic form alone.

Take a symmetric finite slab with equal real boundary vectors in L2(T),
or the infinite-time ground-state limit. For observables supported strictly
on one side of a reflection plane, integrate that half's scalar potentials
and internal q variables first. Reflection maps spatial q evenly and an
oriented electric time-link cochain oddly. The reflected expectation is
the norm square of the resulting boundary vector (a vertex plane), or
its nonnegative quadratic form with T (a plane crossing a link). Thus
all polynomial flux observables on their stated half-slab supports obey
Osterwalder-Schrader reflection positivity. Electric insertions exactly
on the reflection plane need a declared boundary convention and are not
included in this argument. Contact correlations do not change this
half-slab proof. Finite-n ternary records conditioned on compatibility are
not claimed reflection-positive; this is a theorem for their full-moment
Gaussian limit.

The quantum physical configuration space is the gauge quotient T. In a
large regular interior the principal symbol has two transverse modes;
no positive scalar or longitudinal photon is inferred from the temporary
potential variables. On a finite bounded complex, the exact number of
oscillator coordinates is dim C1-rank G, not an assumed count of two per
site. Gauss is G^T M e=0 for separated-time physical states; contact terms
in Euclidean time-ordered electric correlations are kept distinct.

## 3. Exact finite-spacing canonical reconstruction

Choose an M-orthonormal basis of T. In coordinates x=M_T^(1/2)q the slice
operator is Z=M_T^(-1/2) K_T M_T^(-1/2)>0. The infinite-time covariance is
that of the discrete action

    (1/2) sum_t [|x_(t+1)-x_t|^2+x_t^T Z x_t].            (A4)

Functional calculus gives

 E=2 asinh(sqrt(Z)/2),
 Omega=sinh E=sqrt[Z(I+Z/4)],
 C_x=(2 Omega)^(-1),       C_pi=Omega/2.                  (A5)

The canonical CCR is [x_j,pi_k]=i delta_jk. The transfer Hamiltonian,
up to its vacuum constant, is

 H=(1/2)[pi^T(E/Omega)pi+x^T(E Omega)x].                   (A6)

All functions of Z commute. The frequencies are eigenvalues of E, not
sqrt Z at finite spacing. If a is the common physical coordinate step,
the measured transfer rates are E/a. Spatial coboundaries are the integer
cochain differences; their low-frequency eigenvalues scale as a^2.
Temporal electric plaquettes use sqrt Z, while canonical momentum uses
Omega. Their ratio sqrt Z/Omega must be retained, just as in the common
flat native field. No finite-spacing CCR is silently replaced by its
continuum limit. Undoing the fixed M-orthonormal coordinates recovers the
original canonically normalized electric and magnetic flux cochains.

## 4. Native weights, canonical amplitude and the scalar clock

Let u=1+h be the harmonic native hit response used in GR7, and b=u/2
its plus-hit flag. Let q=1-H_c/2=1/u be the independently weighted-walk
response. A q flag is native-compilable: accept on the first fair bit0;
on bit1 accept iff an independently reset weighted walk escapes. Its
ideal acceptance is1-H_c/2. The entire weighted walk must still execute
and retain its inner records in either branch. Caps, failed resets,
conditioned rarity and all accumulated costs are charged by PROOF section5.

For an oriented electric spatial edge xy use

    M_xy=u_x u_y,       p_E,xy=q_x q_y/8.                 (A7)

For magnetic faces choose an explicitly declared time-independent local
sample u_f and use beta_f=u_f^(-2), p_B,f=b_f^2/2. If u_f is a face-center
native probe value on a refined chart, that extra chart and interpolation
error must be supplied. One may instead choose a fixed corner sampling
rule for this finite theorem; that rule gives consistent smooth-background
weights with O(a) error, rather than claiming a symmetric O(a^2) scheme.
Either prescription must be frozen before execution.

The common fluctuation readout F=sqrt(8/n)Z fixes all amplitudes. The
conditional Gaussian precision is then M for electric flux and beta
for magnetic flux. No field-dependent or fitted amplitude rescaling is
performed after identifying the geometric target. For slowly varying u,
the continuum principal density is

    (1/2)[u^2 |E|^2+u^(-2)|B|^2].                        (A8)

Its Lorentzian continuation is Maxwell on the same GR7 metric
-g00=u^-2, gij=u^2 deltaij. Freezing coefficients gives both transverse
characteristic speeds sqrt(beta/M)=u^-2, exactly the scalar characteristic
speed. This principal-symbol match is a positive clock/cone join, not
already a quantitative continuum theorem for arbitrary gradients or
boundaries. The temporal spacing, cochain normalization, smooth-weight
consistency, stable continuum limit and the same detector protocol must
all be retained in such a theorem. Uniform ellipticity1<=u<=2 gives
finite energy stability, but does not by itself supply the remaining
continuum or stress error estimates.

## 5. Affine boundary preparation gives the reciprocal static mean

The exact reciprocal identity can also be the mean of the very same
Gaussian flux state, rather than an unrelated electrostatic proxy.
Declare a stationary boundary potential v=1/2 on the target boundary
and v=1 on the outer escape boundary, with zero tangential spatial
potential. For this subsection the spatial exterior is a finite shell
with two boundary components. Its relative cohomology is not that of
the preceding ball: a boundary-voltage harmonic mode can occur. Fix its
voltage by the affine boundary data, and keep any remaining harmonic
coordinates explicit. The proof here is a static finite Gaussian mean
identity, not an unproved application of the ball's K_T>0 transfer proof
to a shell with zero modes.

On a finite symmetric time slab, fix the spatial potential lift to zero
at both temporal ends, and add the linear compatibility and boundary
constraints specifying this potential difference in each time slice.
This temporal boundary condition is explicit; it makes the static solution
stationary also against all allowed endpoint variations. It is used for
this finite mean theorem, not asserted to be a vacuum boundary state.
Use integer cochains and rational boundary-voltage approximants so the
accepted set is a nonempty affine integer coset. The target normalized
boundary voltage is fixed before source sampling, with rounding error
O(n^-1/2), as in PROOF section3. This amplitude is a declared preparation,
not a derived source charge or universal constant.

The Gaussian conditional mean minimizes the positive quadratic flux
action on that affine set. The unique static electric minimizer has
b=0, e_xy=v_y-v_x, and satisfies

    sum_y M_xy(v_x-v_y)=0                                (A9)

at unconstrained exterior vertices. Stationarity in temporal q also
holds because the electric mean is time independent; stationarity in
spatial potential holds because b=0 and the adjacent-time electric
contributions cancel. Strict positivity of the flux action gives
uniqueness within its affine flux class. Therefore the full Gaussian
mean is this static solution, with the same fluctuation covariance on
the tangent compatibility space.

For M_xy=u_xu_y, original harmonic u gives exactly

    sum_y M_xy(1/u_x-1/u_y)=sum_y(u_y-u_x)=0.

The boundary values select v=1/u. Consequently the independently measured
weighted native hitting response, the mean of the proposed Euclidean
conditional flux measure, and the scalar metric lapse have the same finite
graph value. This is not yet a positive Lorentzian charged-state assertion;
section6 addresses that distinction and supplies a separate state target.
In the unit lattice convention the mean electric flux is the exact edge
difference of this v; its conductance-weighted flux is the original
u-gradient with reversed sign. On a common physical refinement of step a,
set the temporal edge boundary potential to a v, so the electric plaquette
mean is a(v_y-v_x). Dividing the plaquette cochain by a^2 then gives the
physical mean derivative of v. Holding the raw edge voltage fixed instead
would produce an A0 amplitude divergent as1/a and is not the claimed
continuum preparation. This explicit a scaling is a declared boundary
preparation convention, not a fitted stress normalization.
This is an amplitude-preserving finite state/metric join conditional on
the declared affine preparation. It is not yet an executed affine filter.

To prove a full charged-shell quantum reconstruction one must fix or
quantize its harmonic boundary sector and prove reflection positivity
with that time-local choice. It cannot be achieved by subtracting an
arbitrary global-time zero mode. The target's independently measured
inertial energy and the renormalized stress of this prepared state remain
necessary before interpreting the geometric identity G=2T as a physical
Einstein equation or assigning Newton's constant.

## 6. A positive charged harmonic sector, with the Wick factor kept

There is a crucial distinction in section5. A real stationary Euclidean
F_0i mean changes sign under time reflection. The affine real Gaussian
with that mean is not automatically a reflection-invariant positive
quantum state. Its real mean cannot simply be named a real Lorentzian
electric field and inserted into an Einstein stress tensor. Section5 is
an exact conditional-record mean identity. The following finite quantum
construction supplies a different, positive state target on the same
weighted operator, and identifies its remaining preparation obligation.

For a cubical exterior shell with two boundary components, relative H1
has dimension one. Set epsilon=dv for the M-harmonic boundary potential
v=1/u. Then C epsilon=0 and G^T M epsilon=0, where G includes only gauge
functions vanishing on both boundaries. Fix the normalization by the
actual boundary values1/2 and1, not by rescaling epsilon after computing
stress. Write kappa=epsilon^T M epsilon>0. The M-orthogonal gauge quotient
splits into the span of epsilon and a complement on which K is positive.
All statements follow from the finite relative cohomology ranks and can
be checked by exact linear algebra on a specified shell.

Write the harmonic spatial potential as x epsilon. Its Euclidean kinetic
action is kappa*(Delta x)^2/2, with no magnetic potential. The transfer is
Gaussian convolution on L2(R,dx), with Fourier multiplier

    exp[-p^2/(2 kappa)].                                  (A10)

It is positive self-adjoint even though it has no normalizable zero-energy
ground state. Retain this degree of freedom; do not subtract its global
zero mode. On the transverse complement retain the oscillator transfer
and its Gaussian vacuum from sections2–3. The product transfer is positive.
Finite symmetric boundary states yield the same reflection-positivity
proof. A vacuum limit for the free harmonic coordinate is not asserted.

For L>0 the normalized harmonic boundary state

 psi_L(x)=(pi L^2)^(-1/4) exp[-x^2/(2L^2)] exp[i kappa x]   (A11)

has canonical mean <p>=kappa and variance Var(p)=1/(2L^2).
The Lorentzian harmonic electric operator is

    E_harm=(p/kappa) epsilon,       B_harm=0.              (A12)

It follows exactly that

 <E_harm>=epsilon,
 Cov(E_harm)=epsilon epsilon^T/(2 kappa^2 L^2),
 <H_harm>=kappa/2+1/(4 kappa L^2).                         (A13)

Its mean and momentum variance are constant under real-time evolution.
Thus the same native M and reciprocal epsilon select a positive quantum
state with the desired real static electric mean and a controlled energy
error tending to zero. The limiting fixed-momentum sector need not be a
normalizable vector in L2(R); the finite-L states already give explicit
positive approximants. This avoids assuming a normalizable charged vacuum.

In its Euclidean representation the canonical momentum insertion carries
the usual imaginary derivative factor. Equation(A11) is a phase boundary
preparation, not the positive affine electric conditioning of section5.
A native operational attachment must compile or reconstruct this boundary
vector on the same compatible-record source, retain its normalization and
phase readout, and account for any finite potential-lift/reference record.
No such preparation has been executed by this theorem. The flux harmonic
increment is readable as kappa^-1 epsilon^T M e_t; it determines Delta x
but not the initial x origin. A declared boundary reference is therefore
needed for a position-space implementation of(A11). Merely drawing an
extra unconnected Gaussian x and naming it a physical state would not
supply that missing interface.

The physical state target and the classical geometric exterior identity
can now be compared without an arbitrary field-amplitude fit. However,
(A13) supplies a finite integrated harmonic energy and its error, not a
complete local renormalized stress theorem. Curved-vacuum stress, the
continuum local stress limit, target inertial energy and common G remain
to be established. A coherent or charged state's excess classical stress
must not be substituted silently for the full quantum stress.
