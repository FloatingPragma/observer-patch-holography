# Controlled native transfer-ground response

## Source and operation domain

The carrier is the weighted binary total-occupation-two source on twelve ports,
with weights `(25,25,8,8,8,8,8,8,8,8,8,8)` and the uniform selector on its thirty
seams. Its original mixture of C and E attempts has population matrix P and
stationary diagonal pi. The complete 721-slot instrument, including its zero
slots, comes from [the owning producer](../native_record_response/response.py).
The recorded-current reversal and positive scalar-history reconstruction in
[TRANSFER.md](../native_record_response/TRANSFER.md) give

    S = Pi^(-1/2) P Pi^(1/2),  a I <= S <= I,  a = 296/495,
    Omega_x = sqrt(pi_x),  S Omega = Omega,  H = -log S >= 0.                 (1)

The connected source has a unique population stationary law and hence a
one-dimensional transfer ground. The diagonal reference Pi and the pure state
|Omega><Omega| are different states. The logarithm of Pi is also different
from H. The fair operator-history logarithm in
[native_quantum_history](../native_quantum_history/README.md) acts on operator
space and supplies neither of these carrier identifications.

This construction identifies the source configuration basis with the scalar
transfer coordinates and implements H by controlled pulses. It consumes the
existing nine-bit address chart, proper-code coherent services, blank
processors, load/export operations and complex-linear quantum instrument
semantics of [effect synthesis](../native_effect_synthesis/PROOF.md).
These are declared operational premises of that source-service branch.
The fixed compiled protocol supplies this state-and-intervention dictionary;
it does not assert that unread C/E repair is unitary. The full native repair
contraction and its stationary mixed state are unchanged.

All time inputs are exact rational t in [-1,1], with at most 64 bits in each
numerator and denominator. This is a dimensionless controlled transfer
parameter, measured in the unit of the one-attempt transfer logarithm. It is
not elapsed native-attempt time or a physical clock. The concrete replay uses
t=1/10. Every emitted stored primitive must pass the owning interval check.
No theorem about all possible SciPy outputs replaces that acceptance test.

## An exact residual certificate without eigenvalue gaps

The fixed [rational proposal](proposal.json) has 66 real columns v_j and 66
positive rational numbers lambda_j, all on the 2^-44 grid. A numerical
symmetric eigensolver proposes these numbers; it is not evidence for their
accuracy. `spectral.certify` independently checks, using exact rational
arithmetic and directed square-root intervals,

    ||V^T V-I|| <= g = 10^-9,
    ||S V - V Lambda|| <= r = 10^-9,
    a <= lambda_j <= 1.                                                   (2)

The source entries lie in Q(sqrt(2)). Its three unnormalized square-root
weights are 25, 10 sqrt(2), and 8. Integer square root at 100 binary places
bounds sqrt(2); rational interval multiplication and division enclose every
entry of S. For the symmetric Gram error the largest absolute row sum bounds
the operator norm. For the residual, the maximum of its absolute row and
column sums bounds the operator norm. Both gates must pass for the actual
source matrix. No eigenvalue spacing or eigenvector uniqueness is assumed.

Let W=V(V^T V)^(-1/2). This is a mathematical comparison basis. No polar gate
is requested or executed. Singular-value calculus gives ||W-V||<=g, since
|sqrt(x)-1|<=|x-1| for x in [1-g,1+g]. Thus

    ||S-W Lambda W^T|| <= r+2g.                                          (3)

For positive A,B>=aI, the resolvent representation of log gives

    log A-log B = integral_0^infinity [(B+sI)^(-1)-(A+sI)^(-1)] ds,
    ||log A-log B|| <= ||A-B||/a.                                        (4)

The resolvent identity proves the inequality directly, because the integral
of (a+s)^(-2) is 1/a. Consequently (3) bounds the corresponding logarithms by
(r+2g)/a. Duhamel's formula bounds their unitary evolutions by
|t|(r+2g)/a.

The actual compiler uses the normalized rational rays u_j=v_j/||v_j||.
They need not be orthogonal. The bounds ||u_j-w_j||<=2g and
||P_(u_j)-P_(w_j)||<=4g show that replacing one rank-one phase costs at most
8g in operator norm. Unitary telescoping over all 66 phases gives

    ||product_j exp(i t log(lambda_j) P_(u_j)) - exp(-itH)||
       <= |t|(r+2g)/a + 528g.                                           (5)

The product is chronological in the proposal's column order. For the W rays,
the factors commute; no commutation is assumed for the implemented rays.
The scalar logarithm is enclosed by its 96-term rational series in 1-lambda.
Rounding t log(lambda_j) to the 41-bit angle grid costs at most one grid unit
per phase, included in the compiler allowance below. The source ground and
transfer are fixed by (1); the proposal only chooses a finite compilation.

## Preparation with every native record retained

Run the eleven-edge E tree in
[the native configuration detector](../native_effect_coverage/detector.py).
Every original 25-slot command is retained. The complete lazy tree has
25^11 formal raw words and 66*2^11 live words. A live fine branch is
sqrt(w_h)|y_h><x_h|. Both x_h and y_h are decoded from its incoming contrasts
and signed currents, without inspecting an undisclosed input state.

Conditioned on this retained record, apply the existing X routing to move the
work address of y_h to literal address zero, then the existing amplitude-tree
word with weights pi on the 66 source addresses and zero on all padding
addresses. Its target amplitudes are exactly Omega; pi is rational, so this
uses the amplitude-tree input domain directly. The source-address export is
after the entire word. Literal address zero may be a padding address during
the word. The word has at most 12,309 pulses, including at most eighteen
routing pulses, and has ideal fine map

    K_h = sqrt(w_h)|Omega><x_h|.                                         (6)

The full archive is retained. For arbitrary joint rho_SR, every successful
ideal source factors as |Omega><Omega| from archive and untouched R. After
tracing only the preparation archive, the ideal map is
|Omega><Omega| tensor Tr_S rho_SR. The archive can remain correlated with R.
No product input, normalization of an outcome, or free pure-state reset is
assumed. The two eleven-E readouts are ordinary weighted native commands;
the conditional routing and preparation are declared observer control.

## Complete implemented evolution and intervention

Each rational ray phase is the literal word W_j^(-1), addressed phase at zero,
W_j from [the coexistent-filter compiler](../native_coexistent_effects/menus.py),
with zero padding weight. Each has 58,392 pulses, below the unchanged
65,536-pulse source-word admission ceiling. There are 66 separately admitted
words, with a source export after each. The ideal phase preserves the source
chart. The rounded implementation retains its padding as a terminal observed
abort, with 446 internal Kraus rows, and retains every first processor failure
at its original word, pulse and complement label. Subsequent words run only
on the successful source output. These maps preserve arbitrary untouched
references without renormalization.

The intervention is the existing diagonal-phase word assigning i to the
source configuration a=(0,1), index 65. It has 11,270 pulses; a control uses
the same word grammar with zero phases. Its target angle pi/2 is rounded by
at most one 41-bit grid unit. No ancillary address or primitive gate is added.
After evolution, a second complete eleven-E detector supplies the initial
configuration at that readout. Score one when it is a, and zero otherwise.
All processor and padding aborts have score zero and remain in the complete
output. The final detector's disturbance is included, and its remaining source
and reference are retained. The decoded initial configuration is measured;
the terminal configuration after the tree is not substituted for it.

For N pulses and target-grid allowance c, put b=2*10^-11 and Q=2^41. The
existing primitive certificate and product estimate give

    e = Nb/(1-Nb) + c/Q,
    f = N b^2/(1-Nb)^2,
    delta = [e(2+e)+f]/2,
    G = (1-Nb)^(-2)+f.                                                   (7)

Delta bounds complete half-diamond error; G bounds the induced trace norm of
the rounded completely positive map. The exact ideal service is trace
preserving. Stored numerical maps have bounded gain and are not asserted to
be exactly trace preserving. These bounds concern exact products of the
stored primitive matrices. Floating propagation arrays are diagnostics;
hardware error is outside the theorem.

Use c=9 for ground preparation, c=1 for the intervention and c=23 for each
transfer phase. The latter includes the two amplitude-and-phase preparation
errors and the scalar logarithm grid rounding. Complete composition uses
sum_j delta_j product_(k<j) G_k, plus (5) with its conservative gain factor.
`protocol_budget` checks for every allowed t

    complete half-diamond error < 1/11000,   complete gain < 1001/1000.     (8)

The protocol uses 68 coherent words and at most 3,877,451 native pulses, plus
22 fixed-E commands. The direct nine-bit compiler and the lazy complete tree
instrument are its execution domain. The generic word evaluator has a
six-bit/256-pulse ceiling and returns RESOURCE_STOP unchanged. The expanded
tree evaluator also retains its resource stop; it is not relabeled as a
completed tree. These resource boundaries are distinct from mathematical
completeness of the lazy instrument.

## A certified observable response of the prepared ground

Let U(t)=exp(-itH), D=I+(i-1)|a><a| and p=pi_a=625/7505. Since U(t)Omega=Omega,

    p_a(t) = |<a|U(t)D Omega>|^2
           = p |1+(i-1)u(t)|^2,       u(t)=<a|U(t)|a>.                   (9)

The control, with D=I, gives exactly p. Put h=H_aa and L=103/200.
The rational logarithm bound proves ||H||<L and
h>=1-S_aa=32/495. Spectral Taylor expansion gives
u(t)=1-it h+r_t, |r_t|<=t^2 Lh/2, since H^2<=LH. Using sqrt(2)<3/2 and retaining the
nonnegative second-order term in the squared modulus yields

    p_a(t)-p >= p*(32/495)[2t
                  -(3/2)(1+(3/2)tL)t^2 L] > 1/1000,  t=1/10.         (10)

The exact lower fraction is stored in the receipt. Each complete implemented
experiment differs from its ideal instrument by less than 1/11000. A score
in [0,1] changes by at most half the trace distance for equal-trace ideal
channels, but rounded maps can have unequal trace. The witness allowance
therefore uses the full trace-norm error, 2/11000 per experiment, unless the
score is centered and a normalized completion is supplied. The stored
implementation certificate subtracts 4/11000 for the two experiments and
proves an implemented, unnormalized score increment greater than 1/2000.
The ideal output has physical probability semantics under the declared
instrument model; the stored bounded-gain output is a certified approximation.

This response is a controlled finite transfer-ground experiment. A physical
vacuum, autonomous real-time field action, a selected physical clock, a
capacity-compatible coherent embedding and particle masses require separate
proofs. The scalar transfer spectrum is not identified with a relativistic
particle spectrum.
