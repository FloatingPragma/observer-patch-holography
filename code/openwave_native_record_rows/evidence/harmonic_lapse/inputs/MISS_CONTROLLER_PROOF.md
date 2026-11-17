# A native first-hit controller for common clocks and weak metric response

This is a new, unexecuted N+ controller proposal. It extends the same
endpoint law and common occupation source. The acceptance bits are actual
first-hit records; no potential-value table is installed. The powers2 and4
below are declared integer features of the controller, not a uniqueness
claim about OPH. No physical Newton constant or source mass is supplied
by the target's geometric capacity.

## 1. The record rule and a common proposal clock

Use the native first-hit instrument of the row20 adaptive-record packet
on a declared finite three-dimensional routing, with target H and escape
boundary B. A probe starts at the disclosed point x. Native incident
attempts, including all idle coins, continue until hit, escape or cap T.
An escape before target and before the cap is a miss bit1. A target hit
or timeout is bit0. Pad a terminated probe to T calls with native events
on a disclosed zero-load seam. Retain every event and terminal outcome.

Let h(x) be the eventual target-hit probability and r(x)=1-h(x) the escape
probability. Let r_T(x) be the capped miss probability. Then

    0<=r-r_T<=Pr_x(tau_absorb>T).

Fresh probes use independent fair coins and disclosed fresh preparations.
Their reset/preparation implementation and resource cost must be stated
in an executable packet; returning a used unit to its start is not a free
native operation. A source-family experiment may supply independent
prepared probe copies, but must count those preparations and copies.

The finite compiler proposed here instead retains one probe unit throughout.
It returns that unit from its current-known endpoint to the requested start
along a deterministic Manhattan path. In the small box below this needs
at most9 unit moves. Each move gets R fair native attempts on the desired
incident seam, padded to R slots after success. A failed move terminates
the reset, and every remaining slot is a retained zero-pair hold. Thus
the reset uses9R calls and fails with conditional probability at most
9*2^(-R), uniformly in the past. A failed reset produces a refused probe
and a fully padded body. The four-probe reset discrepancy is at most
36*2^(-R) per main proposal.

For every globally scheduled main edge proposal run four independent
T-slot probes. For an internal edge in cell x, accept the main native event
if the first two probes started at x both miss; the remaining two probes
are retained but unused. For a spatial edge x-y, require two misses at
each endpoint. For a temporal edge, ignore the four retained probe bits
and accept. On refusal run a native zero-pair hold instead of changing the
main unit. Every proposal has exactly4T+1 native calls, apart from separately
declared preparation/reset cost. Proposal marks and all rejected records
stay in the archive.

This padding matters: using random probe duration as the occupation clock
would introduce state-dependent holding times and change the record law.
The common occupation decoder samples the main unit at the declared
proposal marks. Their relation to actual calls is fixed and disclosed.
Independent probes do not read hidden main-source loads. The same
controller applies to every internal eigenmode and field channel.

## 2. The exact averaged native kernel

Under the independent preparations the edge acceptances are

    a_internal(x)=r_T(x)^2,
    a_spatial(x,y)=r_T(x)^2 r_T(y)^2,
    a_temporal=1.                                          (1)

They are independent of the main placement coin and symmetric under edge
orientation. A rejected edge contributes an identity transition. If W0
is the fixed total proposal-address weight, endpoint averaging therefore
gives exactly

    P_main=I-L_eff/(2W0),
    L_eff=sum_e w_proposal,e a_e b_e b_e^T.                 (2)

Uniform population remains stationary. Main endpoint movement, probe
words and hold outcomes are different retained parts of one compound
experiment. Equation(2) is its marginal main-record law, not a claim that
the archive factors from the endpoint. All reset records remain retained.

The native occupation theorem now uses W0, not the sum of accepted
weights. Its covariance target is

    C_infinity=L_eff^+-Pi/(4W0).                            (3)

Equations(1)-(3) are exact for the ideal independent-start protocol. The
finite reset compiler can retain a small memory dependence through failed
resets. Its conditional main transition is within36*2^(-R) of this ideal
one; it is not silently called an exactly Markov main chain. At any finite
n-proposal horizon a sequential coupling bounds the difference of complete
main histories by n[36*2^(-R)+4 sup_x Pr_x(tau_absorb>T)]. To attach the
long-field limit, choose the retry and probe caps with n so this bound
tends to0. For the small box, the latter tail is at most(5/6)^T. For
example R_n=ceil(log2(72n^(p+2))) and
T_n=ceil(log(8n^(p+2))/log(6/5)) make the history bound at most
n^(-p-1). At fixed graph a degree-p normalized occupation moment has
maximum absolute size at most a constant times n^(p/2), so this stronger
schedule also makes its coupling error vanish. Merely making history
total variation tend to0 would not suffice for unbounded moment transfer.
Each proposal then costs4(9R_n+T_n)+1 calls. This growing resource schedule
and its proposal marks are declared parts of the controller.

This avoids an unreported normalization change when acceptance changes.
If every probe timeout probability is at most epsilon, the internal and
spatial acceptance errors relative to r are at most2epsilon and4epsilon,
by telescoping products of numbers in[0,1]. The resulting finite Laplacian
and inverse-covariance errors follow from those actual edge errors and a
declared spectral floor. A first-hit holding modification that preserves
h changes probe cost; the padded proposal clock keeps that cost distinct
from the main transfer coefficients.

## 3. Common source symbol and finite-spacing clocks

For a region where r is constant and temporal/spatial proposal weights are
a^(-2), the same native internal Laplacian eigenvalue lambda has

    A=lambda r^2+r^4 q_a(k),
    E_a=(2/a) asinh(a sqrt(A)/2).                           (4)

This is derived from the main native occupation/first-hit transfer, with
the temporal proposals unchanged. It follows by the same one-plane
harmonic recurrence as the sealed common clock theorem. In the
refinement limit,

    omega^2=r^4 |k|^2+lambda r^2.                           (5)

The finite-a rest-clock ratio is

    E_a(lambda,r)/E_a(lambda,1)
      =asinh(a r sqrt(lambda)/2)/asinh(a sqrt(lambda)/2).

It is not exactly r. For each fixed positive band it tends to r, with
species-dependent spacing corrections bounded before refinement. Thus
all native internal rest clocks have the same limiting ratio, even though
their native lambda values differ. No measured rest masses are inserted.

Proper local frequency omega/r and proper local wave number r k turn(5)
into the standard massive relation with rest mass sqrt(lambda). The
effective metric selected by these clocks and phase speeds is

    ds^2=-r(x)^2 dt^2+r(x)^(-2) dx dot dx.                   (6)

Equation(6) applies in an exterior domain with r bounded away from0.
The target itself has r=0 and is outside this regular metric argument.
This is a metric inferred from common clock and propagation records; an
extra conformal factor cannot be chosen independently while retaining
the same rest-clock dictionary.

## 4. Variable r and the exact scalar amplitude correction

For smooth r(x), symmetric edge weights r(x)^2 r(y)^2 have continuum
Dirichlet form integral r^4 |grad psi|^2. After the common internal
projection, the scalar equation has principal operator

    partial_t^2 psi-div(r^4 grad psi)+lambda r^2 psi=0.      (7)

Its static Euclidean operator is positive on the declared finite domain,
so the scalar OS reconstruction is preserved. This claim uses the massive
internal scalar sectors; it does not silently change the common massless
vector reflection to a scalar reflection.

Set psi=r^(-2) phi. Direct differentiation gives

    div(r^4 grad(r^(-2)phi))
       =r^2 Delta phi-2(|grad r|^2+r Delta r)phi.

Consequently(7) becomes

    partial_t^2 phi-r^4 Delta phi+lambda r^2 phi
      +2r^2(|grad r|^2+r Delta r)phi=0.                     (8)

For metric(6), sqrt(-g)=r^(-2). Its minimally coupled massive KG equation,
multiplied by r^2, is the first three terms of(8). The actual source has
the additional effective mass potential

    V_extra=2(|grad r|^2+r Delta r).                        (9)

The same result follows by substituting psi=r^(-2)phi in the measured
quadratic form and integrating its cross derivative by parts, with
specified boundary data. It is not permissible to drop(9) for a generic
profile. In an exterior continuum first-hit region h is harmonic,
Delta r=0. Then V_extra=2|grad h|^2, quadratic in a weak h. Thus the
scalar metric equation agrees through first order in the weak target
response, with an explicit second-order residual. Lattice derivative,
finite-boundary and timeout errors are separate and must also vanish.

The earlier direct metric-density inversion of a raw conductance tensor
assumed an unchanged scalar normalization. Here a measured r-dependent
normalization and the native rest clocks select(6); the residual(9)
records the difference. This is not a contradiction of the metric-density
identity or permission to drop an unwanted potential.

## 5. Test-body rays and the inverse-square limit

In a smooth background the local positive-energy scalar Hamiltonian from
(5) is H(x,p)=sqrt[r^4|p|^2+m^2 r^2]. Its Hamilton-Jacobi characteristic
equation is the mass shell of(6). Smooth narrow-band packets have these
rays in the geometric-optics limit; envelope, lattice and background
gradient errors must be controlled on the propagation interval. The
primitive random unit path is not identified with a Lorentzian body path.
The relevant detector observable is the reconstructed positive-frequency
packet response from the same source correlations and clock.

For a slowly moving packet, H=mr+r^3|p|^2/(2m)+higher velocity terms.
At zero initial velocity its coordinate acceleration is

    d^2x/dt^2=-r^3 grad r=r^3 grad h.                       (10)

All positive internal masses cancel from this leading acceleration.
For an exterior target with h(x)=C/|x| plus controlled multipoles and
lattice/boundary corrections, equation(10) has the attractive leading
inverse-square term -C x/|x|^3. This would be a source-derived
test-body response of the compound record controller. The coefficient C
is first-hit target capacity in the selected normalization. Equating it
to G times the target's inertial energy is a separate essential statement.

## 6. Optical and redshift targets, with the vector limitation explicit

The null characteristic speed of(6) is r^2, so its optical index is
n=r^(-2)=1+2h+O(h^2). A static reconstructed frequency is conserved;
proper detectors compare omega/r. Thus a common spectral line emitted
at x and received at y has proper-frequency ratio r(x)/r(y).

For h=C/sqrt(b^2+z^2) the first-order optical travel-time excess is
2C integral dz/sqrt(b^2+z^2), and the small-angle ray deflection is
2C integral b dz/(b^2+z^2)^(3/2), tending to4C/b for distant endpoints.
These are direct consequences of(6), and give fixed normalization targets
for a native packet/detector comparison. They have not yet been measured
in this proposal.

The flat common Maxwell proof cannot simply be reapplied with variable
r: spatial differentiation no longer commutes with the scalar Green
operator. The four-component curl's global OS positivity, vector
connection terms and detector response require a joint variable-background
proof or a specified new readout/controller. In a sufficiently small
approximately constant patch the existing transverse principal symbol
has the same characteristic cone; this alone proves no full curved
Maxwell transport theorem. Scalar clock and massive-packet claims above
remain distinct from the outstanding electromagnetic clause of row27.

## 7. A small exact finite hit example and next test

As an independently derivable finite controller ancestor, take probe box
max|x_i|<=2, target0 and absorbing escape at max|x_i|=2. Symmetry and the
native first-step equation give the target hit probabilities for interior
axis, face-diagonal and corner starts:

    h_axis=5/22, h_face=1/11, h_corner=1/22.

Indeed a=(1+4b)/6, b=(a+c)/3, c=b/2. Thus the limiting internal acceptances
are(17/22)^2,(10/11)^2,(21/22)^2, and adjacent spatial acceptances are
products of these squared misses. A complete native finite absorbing
kernel can verify these rational numbers and certify the cap-T tails.
Actual independent probe records then decide main repair or hold, with
all failures retained. This finite example is not a far-field test.

A preregistered continuation should first verify this compiled compound
kernel and common clock normalization, including one-miss and ungated
controls. Then use a growing exterior domain to compare multiple internal
clock ratios, packet accelerations, metric volume/area readouts and the
source energy of the target. It must not declare a universal gravitational
constant from capacity alone. Generalized entropy stationarity, the
full tensor Einstein response, cosmological dynamics and a physical
cosmological constant remain further common-source obligations.
