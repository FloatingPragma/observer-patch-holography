# What a calibrated source clock does and does not fix about light bending

The existing source-addressed scalar construction fixes a positive spatial
action form, and its response identifies the normalized action measure.
The filled-source construction supplies a readable clock and a positive
stroboscopic excitation generator. Neither construction identifies a physical
gravitational stress response. We resolve a specific consequence of that
boundary: even the complete calibrated clock-field response and its Gaussian
fluctuations do not select the light-deflection coefficient in the positive
local two-field source-action class below.

This is a constructive nonidentifiability theorem on an explicit mathematical
source class. Its local field registers, source coupling, constitutive action
and probe dictionary are declared extensions of the existing scalar support.
They are not asserted to be selected by A1–A3, or to satisfy the additional
physical Einstein-branch premises. The same dictionary is used for every
member; no measured gravitational result chooses a member.

## 1. Fixed support and the full constitutive family

Let M be the positive dual-cell mass matrix and A the actual 64-site golden
scalar operator of `code/source_scalar_execution`. Remove its explicit unit
mass term and put K=M(A-I). Thus

    q^T K q = sum_edges c_ij (q_i-q_j)^2 + sum_boundary c_ib q_i^2.

Every coefficient is positive. The interior support is connected and touches
the Dirichlet boundary, so K is positive definite. This is a new massless
source response on the existing spatial support, not a claim that the parent's
massive field was massless. Its coordinates and conserved twelve-port addresses
are inherited; its constitutive law and clock remain specified inputs.

Place two real records u,v at each site. For a real gamma and kappa>0 define

    E_gamma(u,v;j) = [u^T K u + kappa (v-gamma u)^T K (v-gamma u)]/2 - j^T u. (1)

The source j=M rho, boundary values u=v=0, units, and external preparation are
identical across the family. Every edge uses the same two-by-two coefficient

    S = [[1+kappa gamma^2, -kappa gamma],[-kappa gamma,kappa]].               (2)

Its determinant is kappa>0. Conversely, EVERY symmetric positive definite
two-field coefficient S with Schur complement S11-S12^2/S22=1 has exactly
this form, with kappa=S22 and gamma=-S12/S22. Normalizing this Schur complement
is the one common clock/source susceptibility calibration. It leaves gamma
arbitrary. Spatial rotations and support automorphisms act on sites, commute
with S, and cannot select its internal cross coefficient. The construction
works with the inherited Gram form and does not replace it by a chosen sphere.

Writing w=v-gamma u separates the two positive Dirichlet forms; the source
term remains -j^T u. For u_star=K^-1 j the common minimum is
`E_min=-j^T K^-1 j/2`, and
`E_gamma-E_min=[(u-u_star)^T K (u-u_star)+kappa w^T K w]/2>=0`.
Thus positivity refers to stiffness and energy above this common minimum,
not to the unshifted sourced energy. Stationarity gives

    Ku=j,   Kw=0,   v=gamma u.                                             (3)

This is a unique stationary state for every source, not an imposed potential
profile. Its clock response du/dj=K^-1 is independent of gamma and kappa.
At any declared inverse temperature beta>0 the complete Gibbs marginal is

    u ~ Normal(K^-1 j, (beta K)^-1),                                       (4)

independent of both constitutive coefficients. Integrating w proves (4);
the determinant of the change (u,v)->(u,w) is one. Unconstrained
relative-entropy minimization against each declared Gibbs reference returns
that reference and hence the same u marginal. With identical constraints on
u alone, the relative-entropy chain rule again gives the same minimizing u
marginal whenever a minimizer exists. Constraints involving v need not do so.
This does not derive or equate those different full references from A3.

A dynamical completion with kinetic form
`(dot u^T M dot u + kappa dot w^T M dot w)/2` has nonnegative energy above
the same minimum for a static source, and equations
`M ddot u+Ku=j`, `M ddot w+Kw=0`. Consequently all u-response functions and
its entire Gaussian process, with matching initial data, agree too. Static
calibration, additional u samples, and a better u clock do not remove this
ambiguity. A read that resolves v or its coupling to propagation can.

This is not restricted to a linear source or a Gaussian clock distribution.
Replace the u energy by a C^2 coercive H(u;j) independent of gamma and kappa,
with static j and a finite partition integral at the chosen beta>0. The same
change of variables gives the normalized Gibbs marginal proportional to
exp(-beta H), and leaves its complete u dynamics unchanged with the same
kinetic form. Coercivity alone is insufficient for the statistical statement:
`H(x)=log(1+x^2)/2` is coercive but its beta=1 partition integral diverges.
In particular, use the resting-matter coupling of the clock action below:

    H_rest(u;j)=u^T K u/2 + sum_i j_i[exp(-u_i)-1],   j_i>=0.
    Ku=j exp(-u),   Hessian H_rest=K+diag(j exp(-u))>0.                    (5)

Since `H_rest>=lambda_min(K)||u||^2/2-sum_i j_i`, its Gibbs law is integrable
at every beta>0. Strict convexity gives a unique minimum. Thus the same
proper-time coupling for source matter and test clocks preserves the ambiguity even with
source backreaction. The finite experiment executes both (1) and (5).
The resting source is held by external support: this is not an autonomous,
generally covariant matter theory, nor a derivation of that coupling from OPH.

## 2. Local operations, retained records and finite error

The exact block coordinate minimizer at site i reads its incident neighbours
and prepared j_i and writes

    u_i <- (j_i - sum_(j!=i) K_ij u_j)/K_ii,
    w_i <- -sum_(j!=i) K_ij w_j/K_ii,    v_i <- gamma u_i+w_i.               (6)

The energy decrease is K_ii[(Delta u_i)^2+kappa(Delta w_i)^2]/2.
Repeated complete cyclic sweeps converge to (3): this is coordinate descent
of a strictly positive quadratic form. Site records have bounded neighbour
degree, local readback, versioned writers and an authenticated public chain.
The finite execution uses explicit dyadic rounding; it retains every read,
write, source preparation and sweep. Its recomputation recipe replaces a
large event tape, not omitted events. It includes both a nonzero w seed and
zero w controls, both source strengths, all stated gamma values and schedules.
The implemented residual, rather than an assumed convergence rate, certifies
the final state. Arithmetic, storage, scheduling and calibration are supplied;
this is not a native integer-repair or autonomous apparatus claim.

For the nonuniform golden grid, the quadratic barrier
`b(x,y,z)=[x(1-x)+y(1-y)+z(1-z)]/6` satisfies `(A-I)b>=1` and `0<b<=1/8`.
The missing boundary values of this positive barrier only increase the
Dirichlet operator's result. The discrete maximum principle therefore gives

    ||u-K^-1 j||_infinity <= ||rho-(A-I)u||_infinity/8,                     (7)
    ||w||_infinity <= ||(A-I)w||_infinity/8.

The certificate recomputes these residuals in exact Q(sqrt(5)) arithmetic.
It also checks the source Gram addresses, every reciprocal edge and the
source-independent action measure. A changed record-production rate changes
the event bill, not K, (3), or the calibrated observables. No event count is
silently identified with physical volume or the field-update clock.

For (5), each local update is the unique root of
`K_ii u_i+sum_(j!=i) K_ij u_j=j_i exp(-u_i)`.
The producer brackets the rounded root with rational exponential enclosures;
the independent verifier uses integer bisection and a different alternating
series enclosure. Nonlinear residuals are enclosed, not represented as exact
algebraic numbers. Subtracting the stationary equation adds a nonnegative
diagonal divided difference to K. The same barrier therefore proves (7)
with residual `rho exp(-u)-(A-I)u`. No fitted stopping rate is used.

## 3. Source solution and one probe dictionary

The corresponding declared continuum action uses the same field coefficient
S and the source `4 pi G rho_m/c^2` in the three-dimensional Gram metric.
Its stationary equations are

    -Delta u = 4 pi G rho_m/c^2,     -Delta(v-gamma u)=0.                   (8)

For a uniform ball of radius R and independently measured mass M, regularity
at the centre and decay at infinity give, with mu=GM/c^2,

    u(r)=mu(3R^2-r^2)/(2R^3) (r<=R),   u(r)=mu/r (r>=R),  v=gamma u.       (9)

Integrating the source equation fixes the radial flux; continuity fixes both
constants. The inverse-square exponent agrees with the existing source
carrier shell-law theorem. The coupling G, the mass readout and physical
position identification are retained inputs. They are not fixed by the
calculated native quasienergy. A finite outer Dirichlet radius L subtracts
mu/L from u; differences of u and gradients are unchanged. Absolute units
still need the same reference convention.

For a refinement control, use the nonnegative compact C^1 density
`rho_m=105 M(1-r^2/R^2)^2/(32 pi R^3)` inside the ball and zero outside.
The exterior remains mu/r; the interior is
`mu[35-35(r/R)^2+21(r/R)^4-5(r/R)^6]/(16R)`.
The potential is C^3: it matches the exterior through three derivatives.
Its normalization and Poisson equation are independently checked, rather than using the desired
potential as a source input.

The continuum resting-matter variant is equally definite:
`-Delta u=(4 pi G/c^2) rho_bare exp(-u)`.
For nonnegative smooth compact radial rho_bare, its energy on D^(1,2)(R^3)
is coercive: exp(-u)-1>=-u and the Sobolev inequality bounds the linear term
by a constant times ||grad u||. The direct method gives a minimizer, and the
strictly convex gradient energy gives uniqueness. Truncation and the maximum
principle give `0<=u<=u_linear`; regularity follows from the bounded smooth
source. Rotational invariance and uniqueness make u radial. Outside the
source it is exactly `mu_eff/r`, where

    mu_eff=(G/c^2) integral rho_bare exp(-u) dx,
    exp(-max u_linear) mu_bare <= mu_eff <= mu_bare.                       (10)

Newtonian source calibration measures mu_eff. The same bare source and
calibration give the same mu_eff for every gamma. The linear interior
polynomial above is not substituted into this nonlinear equation. Both
source laws produce the same exterior redshift/deflection relation at fixed
calibrated mu; this removes the unrelated-source-coupling loophole.
Here rho_bare counts prepared resting masses per coordinate volume on the
fixed background support. It is not held equal to proper metric-volume
density as v varies. External supports fix the source positions; their
physical lengths and stresses are not silently derived from this action.

This family has a controlled finite-volume realization. On a source-addressed
nonuniform grid in a cube of side D, let H be the largest adjacent gap and
M3 a bound on each pure third derivative of this C^3 potential. Taylor's
theorem bounds the three-axis consistency residual by `2 H M3`. The scaled
barrier in (7) bounds interior potential error by `D^2 H M3/4`, plus
`D^2 r_solver/8` for the executed solver residual. A zero boundary on a cube
of half-side L adds at most mu/L relative to the decaying exterior solution,
by the maximum principle. The source is resolved inside that cube. Choose
H and solver tolerance so those terms vanish as L grows; the source-address
family permits such refinement. This is an analytic sequence, not a claim
that the small 64-site run has astronomical resolution.

A potential sup-norm error epsilon alone is NOT a ray-direction certificate.
For off-grid rays, multilinear interpolation adds at most `3 H^2 M2/8`
to that potential error, where M2 bounds the pure coordinate second
derivatives. This follows by three consecutive one-dimensional interpolation
estimates; include it in epsilon before taking directional differences.
For a mesoscopic directional difference of spacing ell, its scalar derivative
error is at most `2 epsilon/ell + ell M2/2`; here M2 must also bound the second
derivative in the direction used. The corresponding scalar component of the
first-order deflection integral on a segment of length T has
error at most `|1+gamma| T (2 epsilon/ell+ell M2/2)`, plus the separately
bounded quadrature error. If both orthonormal transverse components are
approximated this way, multiply this component bound by sqrt(2) to bound the
Euclidean norm of their combined error. For the exterior straight ray, each
omitted tail beyond |z|=Z is at most `|1+gamma| mu b/(2Z^2)`. Taking ell to zero more
slowly than epsilon and refining the quadrature gives the optical limit.
The logarithm of a clock ratio has error at most 2 epsilon. These bounds
separate source discretization, solver, detector and path-tail error; the
finite execution does not supply a physical solar measurement history.

Use for EVERY member the same physical interpretation and test-probe action:

    ds^2 = -exp(-2u)c^2 dt^2 + exp(2v) dx.dx,
    S_clock = -mc^2 integral d(tau),    S_photon = Maxwell[g].              (11)

The proper detector clock is the clock of (11); laboratory coordinate time is
not substituted for it. Pointlike identical two-level clocks with a local
energy gap Delta E have proper angular frequency Delta E/hbar and coordinate
phase rate exp(-u) Delta E/hbar. For a rest-mass gap Delta m, Delta E=c^2 Delta m.
These shared factors cancel in the frequency ratio. The earlier native clock
is an implementation capability for such a declared phase; its spacetime
and physical-unit identification is not a consequence of this extension.

A static observer has d(tau)=exp(-u)dt. Conservation of the Killing frequency
of a photon exchanged between emitter e and receiver r then gives

    nu_r/nu_e = exp(u_r-u_e).                                              (12)

The round-trip reversal contrast in a static vertical resonance experiment is
`2 sinh(Delta u)` in the corresponding fractional convention. In a small
laboratory region `Delta u=g h/c^2` to first order. Here g and h are measured
in the same local proper units: (11) gives `g/c^2=exp(-v)|grad u|` and
`dl=exp(v)|dx|`, so `g dl/c^2=|du|`. There is no separate fitted gravitational
redshift normalization.

Null curves obey Fermat's principle with index

    n=exp(u+v)=exp((1+gamma)u).                                           (13)

Its variation supplies the spatial contribution; using exp(u) alone changes
the model. The first variation about a straight exterior ray of impact b is

    theta_infinity = 2(1+gamma)mu/b + O((mu/b)^2).                         (14)

For a distant source seen at elongation chi by an observer at r_o, the
first-order angular deflection is

    theta_observer = (1+gamma)mu/r_o cot(chi/2).                           (15)

Both follow by integrating the transverse gradient of log n along the same
unperturbed ray; the endpoint integral is `(1+cos chi)/b`.
For gamma=-1, n=1 exactly and EVERY null path is straight at all source
strengths. Its clock redshift is nonetheless exactly (12). For gamma=1,
the first-order deflection is the usual Einstein value. Gamma=0 gives half,
gamma=2 gives three halves. The complete family has no selected coefficient;
the four displayed members are controls, not a search menu from which a
successful member is promoted.

There is also a controlled nonlinear ray check. Write p=(1+gamma)mu/r_min,
where r_min is the actual closest approach. Require r_min>=R so the whole
ray lies in the exterior mu/r field, mu>0 and 0<=p<=1/1000. Merely requiring
the asymptotic impact b>R would not ensure this at finite strength. On this
domain `d[n(r)r]/dr=n(r)[1-(1+gamma)mu/r]>0`, giving a unique exterior
turning point. The exact infinity-to-infinity angle of (13) is

    theta(p)=2 integral_0^1 2/[2-t^2+(exp(-2pt^2)-1)/t^2]^(1/2) dt - pi.    (16)

The t=0 value is its continuous limit. This follows from the conserved optical
impact `b=r_min exp(p)` and the substitution `r_min/r=1-t^2`. The denominator
squared is at least 1-2p. Let I(p) denote the displayed integral, so
theta=2I-pi. Differentiation under the integral gives I'(0)=1
and bounds |I''(p)| by
`4/(1-2p)^(3/2)+6/(1-2p)^(5/2)<11`. Consequently

    |theta(p)-2p| <= 11p^2.                                               (17)

The analytic inequality supplies the rigorous approximation error. Separate
high-precision quadratures in bounded and unbounded radial coordinates check
the recorded numerical angles; their agreement is a numerical cross-check,
not a certified interval enclosure of quadrature rounding error.

This controls the weak-field remainder for that precisely stated geometry;
an infinity-to-infinity bound is not silently reused for the finite-observer
VLBA geometry. Exact zero bending of gamma=-1 needs no expansion or such
transfer. The published VLBA coefficient comparison below remains a
compressed first-order-template comparison, not a reanalysis with our full
nonlinear metric and the original station geometry.

The conformal gamma=-1 member is directly continuous with the supplied
conformal action control already present in the source/action-density work.
The difference is now that a positive local source action actually produces
the profile. It remains a declared constitutive extension, not a derived
physical source law. The other members use exactly the same calibration,
source, boundary, reference clock and universal matter/photon dictionary.

## 4. An identifiable missing condition

At first order write `g00=-1+2u`, `gij=(1+2v)delta_ij`. Direct linearized
curvature gives

    R00=-Delta u,
    Rij=partial_i partial_j(u-v)-delta_ij Delta v,
    G00=-2 Delta v,
    Gij=(partial_i partial_j-delta_ij Delta)(u-v).                          (18)

Outside a nonzero compact source, Delta u=0 and the Hessian of mu/r is
nonzero. Thus the linearized physical vacuum tensor condition `Gij=0` is
equivalent to gamma=1 in this entire static family. This is a single precise sufficient
selector for this paired weak-field coefficient, not a full derivation of
the Einstein equation or all of OPH's physical source law. In particular,
the full exponential gamma=1 metric is not asserted to solve the nonlinear
vacuum Einstein equation. Equivalently, a source-to-stress identification
prohibiting the extra exterior scalar stress
selects gamma=1. Positivity, symmetry, scalar reciprocity, the calibrated
Newtonian force and every u-only statistical test do not imply it.

Every member already has vanishing linearized exterior scalar curvature
and a divergence-free linearized Einstein tensor. The latter is a geometric
identity, not an independent source selector. A scalar vacuum check or Ward
identity therefore cannot replace the physical exterior tensor condition.

There is no contradiction with OPH's existing conditional Einstein theorem:
its physical curvature/stress matching, vacuum and null-balance premises
exclude the other members. The native source results used here have not
established those premises. This is the exact point at which importing the
Einstein branch would supply the requested answer.

## 5. Measurement comparison and scope

The exterior weak-field formula is for fixed finite gamma, with
`|u|, |gamma u| << 1`, impact parameter outside the source and geometric-optics
test photons. The gamma=-1 zero-bending statement and static clock ratio are
exact for the declared metric. There is no uniform finite-strength weak-field
expansion over unbounded gamma. The infinitesimal bending coefficient ranges
over all real numbers while the clock coefficient stays fixed.
The common source-strength calibration is Newtonian (leading weak-field
order). Complete u-response equivalence does not mean every matter observable
agrees: finite-strength proper accelerations and physical ruler distances can
also resolve the spatial-response freedom. No exact higher-order equality of
those different observables is used in the theorem or measurement chart.

The consumed inputs and uncertainties have different, explicit statuses:

| Input or error | Status in this result |
| --- | --- |
| Twelve-port addresses, Gram metric, positive scalar support | Reconstructed from the existing primitive tables; no new A1-A3 physical selection asserted |
| Massless two-field action, boundary and static source | Declared constitutive class; both linear and proper-time source laws retained |
| G, c, source mass/strength, positions and time units | Shared calibration/empirical inputs, not native quasienergy identifications |
| gamma and kappa | Residual universal constitutive freedom; kappa>0, no per-channel fit |
| Matter clock and Maxwell coupling to (11) | One supplied universal dictionary; neither clock nor spatial factor omitted |
| Finite source/solver and detector read error | Maximum-principle and outward exponential bounds on the actual retained records |
| Continuum and spatial resolution | C^3-potential consistency and growing-boundary bounds above; no solar precision assigned to the coarse run |
| Timing | Stationary reads at declared final versions; supplied update ticks are not laboratory seconds; static propagation conserves Killing frequency |
| Weak-field optical approximation | Proved asymptotic-ray remainder; finite-observer comparison explicitly uses the published first-order template |
| Laboratory detector and propagation reduction | Published statistical error, bounded redshift systematic and retained overlapping radio reductions; raw likelihood not rerun |
| Physical source/clock/metric identification | Unproved and not assigned a fictitious small numerical error; positive physical exit is not claimed |

For the finite program the clock log-error is at most twice the u bound;
the optical-index log-error is at most `|1+gamma| epsilon_u+epsilon_w`.
Changing record multiplicity cannot change those fields or their stationary
predictions. A nonzero w initialization does change the finite transient;
its final bound is retained even for gamma=-1, rather than calling a
nonstationary index exactly conformal.
The independent verifier additionally requires the stationary log-index
intervals for gamma=-1 and gamma=1 to be disjoint in every one of the sixteen
histories. Thus the executed optical distinction is larger than the certified
solver errors, not just a difference between nominal branch labels.

The evidence package records published compressed measurements, the original
source URLs/versions and file digests, all transcribed numerical alternatives,
and the exposure history. All comparisons are retrospective. It independently
maps (12)–(15) into the measurement conventions; it does not claim a new fit of
the original photon counts, radio visibilities, satellite orbits or clocks.
Its uncertainty treatment and the distinction between a reported standard
deviation and a systematic bound are specified in `code/paired_gravity/DATA.md`.

Agreement of gamma=1 with those measurements is a check of the supplied
metric-coupling branch. The complete source class also contains an exactly
non-bending conformal member with the identical calibrated clock law, which
does not reproduce the observed solar deflection. Choosing gamma using that
comparison is an empirical selection, not an independent OPH prediction.

The negative exit is the constructed source nonidentifiability (1)–(18).
The retained physical obligation is to select and identify the source's
spatial response, or justify the physical vacuum/null-stress condition above,
on the same source/clock history. It belongs to #779/#729, with calibration
under #736 and common-history identification under #740. A larger run of
the same u-only source cannot decide this missing coefficient.

## Prior results and credit

Schur complements, Gaussian marginalization, coordinate descent, the discrete
maximum principle, gravitational clock shifts and the PPN light-bending
coefficient are standard. The result here composes them on the existing
source-addressed action, proves equivalence of the entire calibrated u
response/statistics across the constitutive family, and identifies a concrete
physical consumer and selector. It does not claim a new theory of lensing.
