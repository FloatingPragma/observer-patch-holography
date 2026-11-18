# Clock agreement does not select a spatial response in the coherent source

## Result and its meaning

A two-field calculation of a calibrated clock-field response leaves open
whether its constitutive alternatives belong to an A1--A3 source. Here the
alternatives are **programs of the already constructed
coherent source**, with its existing port response, instruments, reference,
overlaps and refinement. No new reversible direction or transport primitive
is added. Their complete stationary-clock channels agree at every source
strength, but a finite, physically timed propagation experiment distinguishes
them with an exact rational probability gap.

This closes that particular admissibility objection to nonselection. It does
not turn a programmable source into a law of gravity. In particular, it is
not a prediction that nature implements either program, a derivation of
Maxwell photons, or a completion of #729/#779/#740. The stronger conclusion
is that A1--A3, even with the existing proposed RG and CCG services and an
unchanged clock calibration, **do not entail a unique clock-to-propagation
relation**. Selecting one requires a restriction on source-conditioned
transport, not another improvement in computational capability.

The exact finite result below does not assume a metric, a continuum limit,
Einstein equations, a ray approximation or successful low-band filtering.
The subsequent long-wavelength calculation identifies the spatial freedom
that this leaves in the paired-gravity interface. Its optical interpretation
is kept separate from the finite proof.

## 1. Fixed source and membership, before choosing a program

Use the finite net and directed completion of
[Coherent source clocks, sections 1--5](COHERENT_SOURCE_CLOCKS.md), which
extends the explicit [scalar-seam axiom model](../../code/source_selection_model/DERIVATION.md).
At each carrier the algebra is C^12 tensor M_6 tensor its finite private
buffers. The twelve central projections remain primitive. The complete
primitive response is the reconstructed image

    D_p = diag(hat(v_p/2) + i Q_p, hat(sigma(v_p)/2)),
    Q_p = v_p v_p^T/2 - |v_p|^2 I/6 + I/12,
    im D = u(3) direct_sum so(3).

Same-response Cayley words implement all proper carrier rechartings. The
full Hermitian coordinate grammar supplies all accessible constraints,
including the cross-triple effects needed to observe the twelfth direction.
The normalized matrix/counting trace and the injective positive observer
cover give the inherited unique tracial ontic A3 minimizer. The two-basis
CCG transition problem separately selects identity on the proper code.
These are the **same two typed optimization problems** for every program.

The added work here is a restriction to particular words of that grammar,
not an enlargement of it. A two-mode one-occupation gate is packed in
coordinates |01>, |10> of the four-state proper occupation code. Its U(2)
matrix acts on coordinates 1,2 of M_6 and is identity on the other four
coordinates. This lies in the first U(3). Diagonal phases, real two-mode
rotations and their imaginary counterparts have exact twelve-port control
coefficients. The executable certificate reconstructs these coefficients
from ordered source-response data; its checker instead reconstructs D from
the icosahedral formula. All gates used below factor into these controls.

Packing, moving and unpacking are the inherited complete CP instruments:
receiver overwrite, sender reset and all complement/failure outcomes remain.
They preserve a valid code entangled with other registers. Their extensions
to the whole physical algebra are irreversible. Therefore a
record-controlled gate does not add a thirteenth infinitesimal reversible
response. The rank-deficient-Kraus and closed-unitary-limit argument in the
parent proof applies to these very instructions, including zero angles.
The verifier tests the complete proper-code transfer on every matrix unit,
not just successful code inputs.

All local source values, pulse requests, program labels and detectors are
transported with their charts. Conjugating by the inherited U_a transports
each compiled gate; ordered products and instruments then intertwine by
induction. Selecting a calculation basis is not selecting an absolute port
label. Logical positions are RG addresses, not extra vertices of the
icosahedral screen. All finite programs fit in the parent's support fibres;
its degree-one support cycles, retained old factors, collinear flight
refinements and natural restrictions are unchanged. No independent
spectator supplies the Lie response or quantum dynamics.

Thus these are admitted finite processes in one existing mathematical
source satisfying the stated A1--A3 clauses. This is **not** an assertion
that the separate integer-repair simulator implements this source, nor that
RG or CCG follows from the three axioms. A counterexample that also satisfies
these stronger proposed services is a counterexample to the weaker axioms.

## 2. One stationary source and one clock, common to both laws

To make source ancestry executable, use a finite connected nearest-neighbour
graph with a nonempty grounded boundary, and minimize

    E_s(u) = (1/2) sum_edges (u_i-u_j)^2 - s u_o,
    u|boundary = 0,                 s >= 0.                 (1)

The grounded graph Laplacian is positive definite: its quadratic form
vanishes only on an edgewise constant vector, and grounding forces that
constant to zero. Hence u=s K^(-1)e_o uniquely. The maximum principle gives
u>=0. Equation (1) is a fixed, disclosed source law, not something claimed
to be selected by A3. A classical rational controller can solve this finite
system from its source record s and publish the resulting local records
through RG before preparation of either probe. Both candidates use the
same finite preparation and publication protocol. Its earlier computation,
source loading and waiting lie outside the reported probe ledger, which
starts from the supplied stationary source; this receipt does not certify
their full preparation cost.
No source value is inferred from the desired detector answer.

The primary finite instance is the five-site interval with endpoints
grounded, central source o=2, and unit edge weights. Its actual solution is

    u = (0, s/2, s, s/2, 0).                              (2)

Both laws use the same positive lapse read coefficient

    A(u) = 1/(1+u).                                        (3)

This rational choice avoids importing the preceding exponential metric.
More generally the argument works for any supplied positive clock function
A with A(0)=1. The resting code clock has the native phase channel

    rho -> exp(-i m t A(u) beta) rho exp(i m t A(u) beta),
    beta = diag(1,1,-1,-1).                               (4)

Each clock qubit can run continuously on the native generator
i(E_11-E_22) in its active M6 processor; two co-located copies give the
displayed twofold degeneracy. Packing and preparation precede the
interrogation interval, and readout follows it. These common preparation
costs are outside the reported probe ledger. The clock may run on its own
processor throughout
all probe service and waits: t is physical elapsed interrogation time,
not a count of only successful flights. The finite certificate also
reconstructs the same endpoint matrix by its native phase factors.
The reference gap is 2m in units hbar=1;
restoring hbar multiplies the calibrated energy gap by 1/hbar. The same
code, mass parameter, preparation, measurement and time ledger are used in
both candidates. Their clock channel, not merely one fitted frequency, is
identical for **every** m,t,u, input entanglement and clock read effect.
Clock-only adaptive protocols consequently have identical distributions by
induction. A static frequency comparison uses A(u_e)/A(u_r).

Preparing a source or probe is an actual local intervention, as in the
parent clock construction; a pure clock input is not relabelled the
unconstrained ontic entropy minimum. The stationary source is held fixed
during the probe. No autonomous stress backreaction is asserted.

## 3. A transport freedom already present in the response

Let beta and alpha_j (j=1,2,3) be Hermitian 4 by 4 Clifford matrices, with
beta=Z tensor I and alpha_j=X tensor sigma_j. Define

    R_theta(x)=exp(i theta(u_x) beta),
    S_j(k)=exp(-i a k_j alpha_j),
    W_j=R_theta S_j R_theta^* S_j,
    U_m=exp(-i m T A beta) W_3 W_2 W_1.                   (5)

These matrices specify a finite local program, not a new primitive U(4)
response. In the eigenbasis of alpha_j each S_j is a permutation of four
occupation registers, two flying each way. The change of basis and R are
native two-mode/phase words just described. On a finite interval, outward
registers reflect locally into the opposite channel and wait for the same
flight duration; inward registers fly one edge. The resulting permutation
is bijective. There is no periodic wrap taking a long spatial displacement
in one short tick. Empty registers also move and are charged.

Use two fixed policies, selected before either experiment:

    flat-optical:       cos theta(u) = 1/2,
    curved-optical:     cos theta(u) = A(u)^2/2.            (6)

All strengths use these unchanged functions. Both agree at u=0; neither is
fitted to an optical outcome. The local controller computes the requested
angle from the actual published source record. The angle can equivalently
be specified by the pair (r,sqrt(1-r^2)), r=cos theta. Finite precision gives
ordinary controlled pulse errors, quantified below. The full source admits
both policies before any gravity target is imposed.

**What is being compared.** Merely inserting different controlled optical
elements in one experiment would not disprove a gravitational prediction.
Here form two marked source completions: after the same source/probe
preparation, the default propagation update is the repeated word (5),
with respectively the first or second hard-coded rule (6). The rule is
fixed data of each candidate, not a detector-dependent intervention or an
adjustable fit. A finite controller expands the word into its local
instructions; a stage counter advances regardless of the quantum outcome,
and the published local source value is its only variable angle input.
Finite prefixes use the parent's retained controller/buffer inventory.
Controlled preparations and readouts have the same meaning in both.

Marking a permitted composite as the default update adds no operation to
the A1 grammar. Both updates are natural by section 1, and their complete
matrix-coordinate histories satisfy A2 by composition. The actual A3
ontic and code-transfer objects are unchanged; neither has a constraint
identifying a unique default composite. Thus both markings extend the
same inherited model. This is the precise sense of a source-law
counterexample, rather than an assertion that two different laboratory
interventions should give the same answer. A separately supplied physical
default law could exclude one marking; it would be the additional
selection information that is absent from these antecedents.

Each three-dimensional macrostep has six sequential flights of duration
a/c_raw. Choose the same service budget T_service and pad both programs to

    T = 6a/c_raw + T_service.                             (7)

The witness using one axis has two flights and its own common declared
service budget. Every basis change, angle calculation, code transfer,
reflection and record operation uses that budget or the separately charged
preparation/report interval. On a finite cutoff a finite maximum over the
two programs and all source values bounds the required work and control
strength. A larger budget can always be used. No zero-time operation or
identification of raw flight speed with measured probe speed is needed.
A cofinal continuum family uses the parent's explicit increasing control
resources; this is not a fixed-strength-hardware convergence claim.

For the retained unit-spacing execution c_raw=1 and T_service=2. If N is
the number of serial native pulse layers and F the number of flight layers,
choose delta=2/(3N+2F). Every pulse has a positive delta for loading, a
positive delta for the response, and a positive delta for unloading.
Every flight has delta for loading, duration one for flight or boundary
waiting, and delta for unloading. Distinct sites and flight-mode owners
operate in parallel; each site's pulse workspace executes its word in
order. The finite source inventory includes these distinct processors.
The primary artifact retains every controller stage with its start,
load-completion, action-completion and end time. It then charges one unit
for receiver readout and two for the classical returned record. The
verifier reconstructs the entire order and all times. Source preparation
is a supplied stationary initial condition; its earlier computation and
publication are not claimed to fit in this probe interrogation budget.

## 4. Exact finite experiment: no metric interpretation required

On (2), prepare one excitation at site 1 in an alpha_1 positive eigenstate.
Run one W_1, then make the local occupation read at site 3. The rightmost
S_1 carries the input to site 2. There R_theta^* splits it into positive
and negative alpha_1 eigenspaces with squared amplitudes cos^2 theta(u_2)
and sin^2 theta(u_2), because beta exchanges those eigenspaces. The next
flight sends them to sites 3 and 1. The final on-site R does not change
either site's total occupation. Thus, exactly,

    P_flat(site 3)   = 1/4,
    P_curved(site 3) = 1/[4(1+s)^4],
    gap(s) = ((1+s)^4-1)/[4(1+s)^4] > 0  for s>0.          (8)

The receiver reads its own register. Its central click/no-click record is
then returned by a charged classical flight; the same deadline applies to
both laws. No pooled simulator log is an instantaneous detector.
The complement failure effect is retained with probability zero for this
explicit valid preparation; it is not discarded from the channel.
The five-site receiver is at site 3 and returns its record to site 1,
distance two. The three-dimensional corroborating executions instead
read the receiver at (2,0,0) and return to (0,0,0), also distance two.
Their complete output vectors are independent software diagnostics; only
the specified local occupation effect is the receiver's quantum read.

Equation (8), together with (4), is the nonselection witness. It is a
finite source-operation/read theorem, not a fit, assumed dispersion or
postselected asymptotic result. Strengths 1/20, 1/10 and 1/5 are retained;
the s=0 control has zero gap. The independent replay additionally checks
every input column of the complete reflecting evolution and all three
Clifford-axis implementations. An additional pulse-order control confirms
which actual source writer the split consumes.

If a valid-code preparation has trace distance at most epsilon_p from the
specified input, and the ordered physical channel differs by half diamond
distance at most epsilon_g, any detector probability changes by at most
epsilon_p+epsilon_g+epsilon_d, where epsilon_d bounds its effect error.
Telescoping bounds epsilon_g by the sum of per-instruction half diamond
errors; a unitary operator-norm error is a sufficient bound for a unitary
instruction. Consequently the two certified ranges remain disjoint when
the sum of both experiments' budgets is smaller than gap(s). This finite
robustness statement needs no ray or continuum error estimate.

## 5. What spatial geometry is and is not reconstructed

There is also an exact dispersion statement before taking any limit. On
one homogeneous axis set p=a k, r=cos theta, q=sin theta, and mu=m T A.
The Clifford identities give

    W = (1-2r^2 sin^2 p) I
        -i [2r cos p sin p alpha' + 2rq sin^2 p beta].     (9a)

Multiplication by the native mass phase shows that the complete four-mode
characteristic polynomial of U_m is

    [z^2 - 2 f z + 1]^2,
    f = cos(mu)(1-2r^2 sin^2 p) - 2rq sin(mu) sin^2 p.   (9b)

In particular U_m(0)=exp(-i mu beta) exactly, independent of r, while
its nonzero-momentum spectrum depends on r. With m=0 a signed band phase
is omega(p)=2 arcsin(r sin p), so away from branch crossings

    d(omega/T)/dk = (2ar/T) cos p / sqrt(1-r^2 sin^2 p).  (9c)

At k=0 this is 2ar/T, and its absolute value never exceeds 2ar/T for
0<r<1. Thus the clock-equivalent velocity difference is present in the
finite-duration unitary itself. It is not an artefact of truncating a
continuum series. Both eigenvalue multiplicities and the entire symbol
are retained in the certificate, including nonzero mass. No selected
Floquet logarithm is being equated to a physical stress tensor or
unwrapped laboratory energy; the earlier energy-interpolation boundary
is unaffected.

For constant u, W_j(0)=I. Differentiating the *actual* word (5) at k=0 gives

    W_j = I - 2 i a k_j cos theta alpha'_j + O((a k_j)^2),
    alpha'_j = exp(i theta beta/2) alpha_j exp(-i theta beta/2).

The primed matrices obey the same Clifford relations and anticommute with
beta. At fixed bounded k,m, as a,T tend to zero with a/T fixed, (5) yields

    H = m A beta + v(u) sum_j alpha'_j k_j,
    v(u)=2a cos theta(u)/T,
    E^2 = m^2 A^2 + v(u)^2 |k|^2.                         (9)

This is a first-order expansion of a finite word, not an uncharged Trotter
simulation of arbitrarily many spatial flights in time T. Its remainder
is bounded by (sum_l ||X_l||)^2 exp(sum_l ||X_l||)/2 for the small
generators X_l after conjugating out the fixed R factors; the certificate
tests a stricter unitary telescoping bound without the exponential. All
finite quasienergies and all input modes remain; (9) describes only the
stated long-wavelength branch.

More explicitly let L=|m T A|+2a sum_j |k_j|. The second derivative of
the product with every small exponent scaled by t has norm at most L^2:
each differentiated factor is bounded by its generator norm and all
remaining factors are unitary. Taylor's integral remainder gives
||U-I+i T H|| <= L^2/2. Since ||T H||<=L, the same bound applies to
exp(-i T H), hence ||U-exp(-i T H)||<=L^2. A telescoping product gives
||U^N-exp(-i N T H)||<=N L^2. At fixed physical elapsed time N T and
bounded k,m, this tends to zero linearly in T when a/T is fixed. This
proves the stated homogeneous dynamics limit, not just a derivative fit.

For smooth varying coefficients, the lower-order term can be calculated
explicitly. Acting on a smooth spinor, write

    S_j psi = psi - a alpha_j partial_j psi + O(a^2),
    C_j = alpha_j + R alpha_j R*,
    R alpha_j (partial_j R*) = (1/2) partial_j C_j.

The last identity uses R=exp(i theta beta) and {beta,alpha_j}=0. Hence
the actual two-flight word, including its intermediate source writer, is

    W_j psi = psi - a [C_j partial_j psi
                       + (1/2)(partial_j C_j) psi] + O(a^2).

With V_j=(a/T) C_j=v alpha'_j, the local differential generator is

    H psi = m A beta psi
            - i sum_j [V_j partial_j psi + (1/2)(partial_j V_j) psi].

The gradient term makes this expression formally Hermitian on compactly
supported spinors in the coordinate counting limit. It is generally
nonzero and must not be discarded from the dynamics. It has differential
order zero, so the principal cone is still (9). The mass multiplication
adds no derivative of A at this order. The local one-step remainder is
O(a^2) for fixed a/T and bounded second spatial derivatives of the
coefficients and spinor. The regression test executes the native words
on a varying profile in all three Clifford directions, checks quadratic
one-step error, and detects omission of this gradient term. This is not
a uniform geometrical-optics estimate for the finite detector (8).

Calibrate once at u=0, where v_0=a/T. If a static isotropic metric is used
to encode **both** the massive rest frequency and this principal cone,
its coefficients are forced by those two reads:

    ds^2 = -A(u)^2 v_0^2 dt^2 + B(u)^2 dx^2,
    B(u)=v_0 A(u)/v(u),       n(u)=v_0/v(u).              (10)

The metric is reconstructed here, not supplied to compute the finite
reads. It is an operational principal-symbol metric; it neither identifies
the microscopic record cone with this cone nor derives Einstein dynamics.
For (6), respectively,

    flat-optical:   B=A,      n=1,
    curved-optical: B=1/A,    n=(1+u)^2.                  (11)

Both have the same nonconstant lapse. The first has straight principal
rays; the second has a nonconstant refractive index. More generally any
smooth positive v(u)<=2a/T can be realized by theta=arccos(vT/(2a)); values
strictly below the upper endpoint avoid derivative singularities. Thus
clock calibration leaves a **function**, not just one fitted constant,
undetermined in this source.

For the spherical gravity interpretation, the continuum source scaling
must also be stated. The common three-dimensional finite-difference
action is (a/2) sum_edges (u_i-u_j)^2 - a^3 sum_i rho_i u_i,
so its Euler equation is K u=a^2 rho. It approximates the separately
disclosed continuum action integral (|grad u|^2/2-rho u) d^3x, with
Dirichlet grounding taken to infinity. Keeping an unscaled point-record
strength fixed while changing a would not keep the same physical mu.
For a smooth compact spherical rho, the continuum equation -Delta u=rho
and decay at infinity give exterior u=mu/r, mu=integral rho/(4 pi).
This source normalization is common to both policies; it is not a derived
laboratory stress map or an input to the exact finite counterexample.
Equation (11) then gives, at first weak-field order, respectively
zero and 4mu/b bending of the principal rays, whereas the clock shift is
the same at that order and exactly identical under (3). This is the
gamma=-1 versus gamma=1 spatial-response distinction of the preceding
paper, now reachable inside an admissible source process. The bending
statement is the ordinary principal-ray weak-field interpretation of
(9)--(11), not a second physical measurement comparison and not a proof
that these massless code modes are electromagnetic photons. The decisive
finite distinction is already (8).

## 6. Why A3 cannot silently choose the curved policy

The source's actual ontic and code-transfer A3 problems in section 1 do
not optimize the program, source law or matter species. Both policies
use their same solutions. Reinterpreting an optimizer over states as an
optimizer over Hamiltonians would change the typed axiom.

There is also a sharper, properly typed test. Put the two program labels
in two ports of an existing local central record and allow their mixtures on one
fixed algebra. Condition on the same source and probe preparations and
all the common stationary-clock data. Those clock effects have equal
expectations on both labels by (4); neither branch is excluded. The
actual positive effect reading the flat label has value one in a feasible
state. The inherited [finite A3 support theorem](../../Lean/Geometry/SourceStateSelection.lean)
therefore assigns it strictly positive probability at any attained A3
minimum with faithful references and a scored positive pullback for this
record. This statement does not require a uniform optimizer. In the
simple initial-register problem with the symmetric tracial reference,
the unique label weights are exactly 1/2,1/2. The other ten ports are
excluded by the specified preparation constraint, not deleted from the
primitive centre or from the faithful twelve-port reference. Treating
these labels as a state-selection object is a separate test of an attempted
selector; it does not reclassify the fixed-law comparison as random
external interventions.

Running the same controlled two-branch experiment after that preparation
gives P(site 3)=p/4+(1-p)/[4(1+s)^4], strictly larger than the curved
policy's probability for s>0 whenever p>0. To select the latter exactly,
the feasible constraints must eliminate the flat label (or the reference
must cease to be faithful there). Finer clock calibration alone cannot do
so. This application supplies the actual same-algebra feasible witnesses
and positive program-record effect; it does not invoke a theorem about
an unspecified space of candidate theories.

## 7. The additional physical information actually needed

Within this fixed family the necessary and sufficient restriction for the
curved operational cone is the source-to-transport response

    v(u)/v_0 = A(u)^2, equivalently cos theta(u)=A(u)^2/2. (12)

It is a functional condition over source strengths, not one normalization
that a clock experiment can determine. The exact local statistic
P(site 3)=cos^2 theta(u_2) independently measures this response and is a
falsifier of (12). Because cos theta is positive on this family, it
determines theta on the specified arccos branch 0<theta<pi/2, without
claiming to determine all possible signed coin angles outside this family.
Source-conditioned propagation measurements, or a separately justified
physical action principle that implies (12), can remove the freedom.
Simply declaring (12) to be A3, equivalence, or covariance cannot.

The finite read also gives a positive reconstruction statement within this
specified family. If p(u) is the receiver probability in (8)'s protocol,
then v(u)/v_0=2 sqrt(p(u)), B(u)=A(u)/(2 sqrt(p(u))) and
n(u)=1/(2 sqrt(p(u))). These are identified functions when the clock and
source-conditioned local transport reads are supplied, without assuming
a metric beforehand. A finite list of strengths does not determine an
arbitrary function between them. These reads would be additional physical
information, not a second optical fit hidden in a clock-only prediction.

In the two-policy selector test, write Delta=P_flat-P_curved>0. An exact
additional constraint P=P_curved forces the flat-label weight to zero
because P-P_curved=p_flat Delta. A tolerance epsilon instead gives only
p_flat<=epsilon/Delta. This is the exact selection cost; finite precision
must not be promoted into an exclusion of every competing policy. More
generally positive probability bounds p_minus<=p(u)<=p_plus and lapse
bounds A_minus<=A(u)<=A_plus give the monotone metric interval
[A_minus/(2 sqrt(p_plus)), A_plus/(2 sqrt(p_minus))]. No laboratory
probabilities or confidence levels are fabricated here.

This result does not adopt (12) as a new axiom or pretend that it has been
independently motivated. It identifies exactly where an additional law
would act and proves why the currently stated antecedents cannot replace
it. No further conditional compiler is needed to establish this negative
answer. The physical parent still needs a selected source/action and a
matter/photon identification; our admissibility theorem does not close it.

## Reproduction and attribution

See [the finite contract](../../code/source_gravity_admissibility/CONTRACT.md)
and [execution instructions](../../code/source_gravity_admissibility/README.md).
The independent checker validates exact source solutions, twelve-port
controls, full reflected-channel matrices, the analytic finite detector
gap, common clocks, a complete error ledger and the Clifford expansion.
Analytic model membership and the all-strength functional theorem are
proved above; a finite test is not labelled their formal proof.

The Dirac principal-symbol interpretation is standard; compare
[Obukhov, *On gravitational interaction of fermions*](https://arxiv.org/abs/gr-qc/0112080).
We do not import its gravitational field equation. Product expansions and
unitary telescoping estimates are standard; see
[Childs et al., *Theory of Trotter Error with Commutator Scaling*](https://arxiv.org/abs/1912.08854).
The source-specific contribution is the admitted, charged native word,
its exact clock-equivalent finite counterexample and its typed A3 test.
