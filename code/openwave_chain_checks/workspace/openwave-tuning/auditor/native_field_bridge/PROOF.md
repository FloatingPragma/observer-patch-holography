# Native current feedback with a finite local electric solver

The executable produces an electric-energy response by changing the seam
scheduler of the native process. The integer endpoint repair and fair coin
are unchanged. The energy coefficient, inverse temperature, fixed charge
reference and state-dependent Poisson rates are declared controller inputs.
The result is a finite dissipative response, with complete current records
and waiting times. It supplies neither transverse waves nor orbital levels.

## Source and controller

The development carrier is the existing periodic cubic routing at side two:
96 ports and 288 seams. Its internal twelve-port cells and external pair
lists are reconstructed independently in `field_feedback.py`. Every port
has degree six. A single native unit load is prepared at a known port. The
reference apparatus supplies a fixed compensating charge at port zero.
This reference does not undergo repair. Its pinning is a construction input.

Orient each seam from its lower to higher port. Let B have +1 at the first
endpoint and -1 at the second; L=BB*. With initial load n0, the observer
updates its counter only through the native C record j:

    n_record = n0 - B F,       F_e <- F_e + j.

This equals the actual source load after every event, including j=0. The
source state is supplied solely to the literal endpoint routine extracted
from `model/observer/native.py`. The scheduler receives the record counter.
The complete one-load domain has 96 states. All 55,296 combinations of
source state, seam and coin were compared against an independent endpoint
floor/ceiling calculation.

Put rho=n_record-delta_0 and G=L+ on the mean-zero subspace. The ideal
electric energy is U=kappa rho*G rho/2. Its graph dependence is computed
from the native routing. The coefficient kappa is supplied; the construction
does not determine an electromagnetic unit or fine-structure constant.

## A local circuit replacing a Green inverse

For a configured neutral rho, initialize phi0=0 and repeat K rounds of

    phi_(k+1)(x) = [6 phi_k(x) + sum_(y adjacent x) phi_k(y) + rho(x)] / 12.

Each round uses only local and neighbouring coordinates. The exact result is
phi_K=G_K rho with G_K=(1/12) sum_(j<K)(I-L/12)^j. The executable computes
these columns by that stencil. Its controller never calls an inverse or
eigensolver. The dense inverse occurs only in the independent error check.

At side two the eight Fourier symbols have rational entries. Positive
rational LDL decompositions of their twelve-dimensional matrices, after
lifting the sole constant mode, certify lambda2>1/2. The degree bound gives
lambda_max<=12. All exact diagonal certificates are in the receipt. Hence

    ||G-G_K|| <= 2 (23/24)^K = epsilon_K.

Since ||rho||^2<=2, |U-U_K|<=kappa epsilon_K and the error in one-hop energy
differences is at most 2 kappa epsilon_K. At K=512, epsilon_K is
6.879015378248898e-10. The measured maximum energy error is
8.092415626492766e-12. The spectral truncation bound is an exact-arithmetic
statement; this floating-point execution is checked by the stated numerical
residual and does not supply an additional rigorous rounding certificate.

The finite controller caches one energy per admitted source state, obtained
from those local circuits before the Poisson experiment. Online state
selection consumes the known preparation and current records. This avoids
a Green oracle, while retaining an explicit centralized scheduler and a
compiled lookup service. The 512 local arithmetic rounds are preparation
work, not secretly counted as native Poisson time. A finite-speed distributed
implementation with online relaxation latency would define a further clock
adapter. No relativistic propagation claim is made for this controller.

## Unchanged native moves with an electric response

For a source at x and an incident seam to y, choose its attempt rate as

    r_e(x) = nu exp[-beta (U_K(y)-U_K(x))/2].

Nonincident seams have rate nu. On every scheduled event the fair native
coin is still drawn. One coin moves the load and the other holds it.
Consequently the actual position jump rate is r_e(x)/2. No favourable coin
is selected and no zero-current event is discarded.

For every bond, the forward/reverse jump-rate ratio is exp[-beta Delta U_K].
Thus the full position generator has the exact stationary law proportional
to exp[-beta U_K], with detailed balance. Its semigroup, including diagonal
holding probability, is a normalized stochastic matrix.

An observable response is the energy carried by the actual native current
records. Starting from a uniform position preparation, write u for the
unit-kappa energy table and Delta u=u(y)-u(x). Then

    d E[u]/dt at t=0 = -(nu/N) sum_edges Delta u sinh(beta kappa Delta u/2).

This is strictly negative when beta*kappa>0 and u is nonconstant; it is zero
with the same carrier and detector at kappa=0. The formula follows directly
from the two native transition rates. The receipt verifies it against the
complete generator and propagates all 96 output positions at time 1/2.
For beta=4 and nu=1 the initial rates at kappa=0,1,2 are respectively zero,
-0.005318327866635085 and -0.011060236987688573 in the declared source units.
They are response coefficients, not physical values fitted to observation.

The relation is dissipative and describes motion in the supplied electric
energy. It does not establish Newton's second law, a particle mass, radiation,
or a hydrogen spectrum. The Green construction can use larger cubic graphs,
but this receipt executes only the stated finite development carrier.

## Complete clocks and probabilities

At state x, Lambda(x)=sum_e r_e(x). A waiting time has density
Lambda exp(-Lambda t), the seam has probability r_e/Lambda, and either coin
has probability 1/2. A complete fine path therefore has density

    product_i [exp(-Lambda_i Delta t_i) r_(e_i)(x_i)/2]
      times exp(-Lambda_final remaining_time).

Every simulated record retains its elapsed time, total rate, seam rate,
coin and current. The terminal no-event interval remains in the density.
A finite attempt cap has a distinct failure output. The 36 development
paths contain 4,732 native events, including 4,704 zero currents; all were
independently replayed from their records. Explicit cap-failure and
zero-duration cases are retained. These paths are individual executions;
the complete output-law claim comes from the full generator and its
stochastic semigroup, not a small-sample frequency fit.

For the exact finite circuit versus the ideal Green scheduler, only the six
incident seam rates differ. With beta=4, kappa=nu=1 and duration T=1/2,
a coupling of complete marked Poisson paths gives

    TV <= T*6*exp(4)*[exp(4 epsilon_K)-1] < 4.507e-7.

The bound uses 0<=U,U_K<=2, so each ideal rate is at most exp(4). It includes
the seam and coin records, their holds, and terminal waits. A common attempt
cap is a measurable output map and cannot increase this bound. As above,
this is the exact-arithmetic circuit bound, separate from floating rounding.

## Transfer boundaries and next adapter

The uniform native attempt clock has T=I-L/(2M). The cubic acoustic symbol
lambda(k)=|k|^2/8+O(|k|^4) therefore gives diffusion coefficient 1/(16M)
per attempted event. Poissonizing each seam at rate nu gives generator
-nu L/2 and coefficient nu/16. The separately installed coherent H_m has
speed 1/sqrt(8) in its own Hamiltonian time. These clocks are distinct.

The next narrow adapter is a certified finite-precision controller compiler:
enclose the stencil energies and exponential rates, then propagate their
rate intervals into the complete path bound. A distributed controller also
needs to retain its finite computation delays. A transverse electromagnetic
model requires an additional source-executed field action and its coupling;
the present response supplies no such action by changing a clock label.

The frozen row-one spatial charge and readout files are unchanged. No official
platform status is promoted by this development receipt.
