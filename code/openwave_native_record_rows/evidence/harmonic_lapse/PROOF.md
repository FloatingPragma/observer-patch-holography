# A native hit controller with an exact scalar metric density

This is a later, unexecuted proposal. It does not amend the preregistered
miss controller. All endpoint laws stay native; accepted proposals are
selected by additional recorded hit experiments and fair native bits.
The following theorem supplies an exact scalar metric density and a
separate static potential readout on the same target. It does not yet
identify a physical electromagnetic stress, target inertial energy or G.

## 1. Plus-hit bits and three conductance factors

Let h(x) be a native first-hit probability with target value1 and escape
value0, and put u=1+h, so1<=u<=2. Generate a fair native bit on a known
empty seam. Accept immediately for one bit value; for the other value
accept iff an independently reset native probe hits the target. A probe
is still executed and padded in the immediate branch so the record and
proposal cost have fixed length. The resulting ideal probability is

    b(x)=(1+h(x))/2=u(x)/2.

Use four independent such compound slots per main proposal. For a temporal
edge at spatial x require all four plus-hit flags. For an internal edge
at x require the first two plus-hit flags and a specified value of each
of the last two independent fair bits. For a spatial edge require a
specified value of all four fair bits, retaining but ignoring their hit
records. Their acceptance factors are exactly

    a_time=u^4/16,
    a_internal=u^2/16,
    a_space=1/16.                                          (1)

All factors lie in[0,1]. There is no division by an estimated probability
or potential table in the controller. The powers and fair-bit choices
are declared parts of this new candidate, not a claim of unique selection
by the unchanged endpoint law. A no-target control has h=0 and all
three factors equal1/16, a common slowdown of the existing source.

Capped hits and failed resets are retained. If h_T is the actual completed
hit probability before cap T, the ideal capped b_T=(1+h_T)/2 differs
from b by at most half the timeout probability. Product errors telescope.
The finite native reset compiler and its growing-cap whole-history
coupling bound carry over from the miss-controller proof. The extra four
fair-bit events per main proposal must be counted. No unbounded repeat
loop is required by(1).

## 2. Exact scalar metric operator

The unchanged binary one-unit endpoint kernel, averaged over the gates,
has the symmetric graph Laplacian with the weights(1). Use the common
spatial/time proposal spacing a and the same native internal eigenvalue
lambda. In a smooth static refinement its Euclidean operator is

    L_eff=(1/16)[-u^4 partial_t^2-Delta+lambda u^2].          (2)

The time coefficient depends only on space, so there is no hidden time
derivative of u. Spatial edges have constant weight, so no extra spatial
gradient of u is generated. Internal projection commutes with u because
every internal seam in a cell uses the same gate.

For the Lorentzian metric

    ds^2=-u^(-2) dt^2+u^2 dx dot dx,                       (3)

one has sqrt(-g)=u^2, -sqrt(-g)g^00=u^4, and
sqrt(-g)g^ij=delta^ij. Its minimally coupled massive scalar equation is
therefore exactly the continuation of the bracket in(2). This agreement
does not require a position-dependent field rescaling, and has no
additional gradient mass potential. It is stronger than the first
miss-controller's agreement through first order in a weak h.

The factor1/16 is fixed globally by the four fair bits. If Y denotes the
native occupation Gaussian readout, then phi=Y/4 has the covariance of
the bracket operator, with its explicit remaining contact subtraction.
This fixed amplitude convention does not fit a gravitational coefficient.
It agrees with the no-target reference and is used for every internal
band. Finite lattice, covariance and cap errors must still be bounded.

The static Euclidean operator has positive time coefficient D=u^4 and
positive spatial-plus-mass operator B=-Delta+lambda u^2. The canonical
time form is obtained with D^(1/2)phi, whose frequency operator is
sqrt(D^(-1/2) B D^(-1/2)). This proves scalar OS positivity and defines a
common positive time transfer even when D and B do not commute. The
native plane-hitting recurrence is the corresponding generalized
eigenvalue recurrence; it must use this weighted inner product rather
than replacing D by a constant.

Locally, the mass shell is omega^2=u^(-4)|k|^2+lambda/u^2.
Thus the universal rest-clock ratio is1/u and the characteristic speed
is1/u^2. These observations independently select the same metric(3).
As before, finite-a asinh corrections and the main/probe resource clock
are retained. A particle trajectory claim needs a reconstructed packet
and detector comparison; it is not a claim that a raw random-walk path
is a Lorentzian geodesic.

## 3. A separately measured reciprocal potential, exactly on a graph

The new temporal/internal gates use u, but(3) has lapse v=1/u. The lapse
also admits an independent native static readout. Let the original
unweighted first-hit graph have conductances w_xy, so its source record
identity in the exterior is

    sum_y w_xy (u_x-u_y)=0.

Give a second native probe walk conductances

    c_xy=w_xy u_x u_y.

A symmetric factor u_x u_y/4 is compiled using two independent plus-hit
flags at the endpoints; the common1/4 scale does not alter its hitting
probabilities. It is another record-driven compound walk, with all inner
and outer resets, timeouts and holds to be retained in an actual compiler.
Direct substitution gives

    sum_y c_xy(v_x-v_y)
      =sum_y w_xy (u_y-u_x)=0,     v_x=1/u_x.               (4)

Since v=1/2 at the target and1 at escape, uniqueness of the finite
Dirichlet problem implies

    v(x)=1-H_c(x)/2,                                       (5)

where H_c is the second walk's actual target-hit probability. Thus the
reciprocal is obtained from a second normalized first-hit response,
independently of naming it a field potential. On an edge its static
conductance flux is

    c_xy(v_x-v_y)=w_xy(u_y-u_x).                            (6)

Native mean-current/residence identities attach this flux to records once
the boundary-source experiment and normalization are specified. Equations
(4)-(6) are exact finite-graph identities, not a continuum Coulomb fit.
They remain new unexecuted controller claims until their nested native
apparatus is compiled and tested.

For the small unweighted probe box of the preceding proposal, the three
interior u values are27/22,12/11,23/22. The independent weighted target-hit
values predicted by(5) are10/27,1/6,2/23. These are fixed rational targets
for a later native test, with no fitted amplitude.

## 4. The geometric exterior Einstein-Maxwell identity

Suppose a controlled common refinement has smooth positive harmonic u
outside the target. Metric(3) is then the conformastatic electrovac
geometry often called Majumdar-Papapetrou. For an independent check of
this geometric target and its conventional normalization, see
[Varela, equations16-30](https://arxiv.org/abs/gr-qc/0201073).
This citation supplies a target identity, not native experimental evidence.

Here is the normalization in the canonical Maxwell convention of the
common free-field reconstruction, with
T_mu nu=F_mu alpha F_nu^alpha-g_mu nu F^2/4. Put r=1/u and tentatively
A_0=r,A_i=0. Direct differentiation gives

    Delta r=2|grad r|^2/r,
    G_00=r^2 |grad r|^2,
    G_0i=0,
    G_ij=r^(-2)[-2 r_i r_j+delta_ij |grad r|^2].            (7)

For this static field the same coordinate components of T are

    T_00=(r^2/2)|grad r|^2,
    T_0i=0,
    T_ij=r^(-2)[-r_i r_j+(delta_ij/2)|grad r|^2].           (8)

Hence the geometric pair obeys G=2T. Its source-free Maxwell equation
is div(u^2 grad r)=0, which follows from u^2 grad r=-grad u. The trace
is zero. If the physical field were instead A_0=a r, the tensor
coefficient would be2/a^2. This exposes exactly why a normalization join
cannot be replaced by an arbitrary rescaling after the calculation.

Equation(5) makes the proposed scalar potential shape source-readable,
while(1)-(3) make the scalar metric response source-readable. To promote
(7)-(8) to a physical common-source gravity result, the same native
Maxwell state/readout must realize the potential mean and its canonical
stress with the fixed amplitude convention. A scalar conductivity field
does not by itself prove the variable-background Maxwell algebra,
reflection positivity, two transverse modes or charged-source Ward law.
The target's mass/charge/energy must also be measured in the same clock.
None of those joins is discharged by the geometric identity alone.

In particular no value G=1/(4pi) is claimed merely by reading G=2T and
writing2=8piG. That inference would require the missing independent stress
normalization and source inertial-energy identification. This is a large
exterior solution family of a specified new controller, not a derivation
of arbitrary Einstein dynamics or a cosmological-constant value.

## 5. Independent proper area and volume

Metric(3), once established by the scalar records, has spatial volume
element u^3 d^3x. The metric's area form on a coordinate sphere is
u^2 times its Euclidean area element; for a radial profile,
A(r_coord)=4pi r_coord^2 u(r_coord)^2. These are independent geometric
counts/quadratures of the same measured u, not definitions from entropy
balance. Geodesic-ball volume and fixed-volume variations still require
the appropriate distance construction and small-radius error estimates.

The new positive route therefore supplies concrete simultaneous candidates
for scalar clocks, geometry and an independently probed static potential.
It leaves the generalized-entropy stationarity and full dynamical matter
response to further tests, while preserving their separate meanings.
