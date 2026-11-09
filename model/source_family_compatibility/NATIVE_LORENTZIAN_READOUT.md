# A complete native long-time spectral experiment on the CAR dimer

The source is the full sixteen-row two-site CAR action in the declared
signed Hubbard chart, at c=1/8 and U=1.
The scalar-shifted centered clock is Q=T+(U/2)sum_x(N_x-1)^2. Its unique
normalized ground has energy e0=(4-sqrt(17))/8 and direction

    |3>+(4+sqrt(17))|6>-(4+sqrt(17))|9>+|12>.

`local_spectral.dimer_source` reconstructs Q and its complete projectors
from the source. Its local mode 0 addition and removal measures both have
energies sqrt(17)/8 minus/plus 1/16 and original weights
1/4 plus/minus1/(4sqrt(17)). No charge is normalized after preparation.
The coherent clock is Q-e0, equivalently the same centered K referenced
to its ground; Q differs from K only by a scalar.

## A stopped native experiment

Supply a blank four-bit field and prepare an integer approximation to the
displayed direction using the existing full native Gray tree. It uses 57
literal source pulses. A phase-exact occupation copy and X0 implement the
complete charge instrument; read its tag and retain both classical charge
columns and all four proper-pointer zero maps. Supply a fresh blank clock
bit and apply the native Hadamard.

At stage j=0,...,m-1 supply a fresh blank selector. Its ideal amplitudes
for stop and continue are 1/sqrt(257),16/sqrt(257). Measure that selector,
retaining all export failures. If it stops, apply a Hadamard to the clock,
read it and record charge q, time j and sign b. If it continues, apply one
controlled clock step exp(-i Delta(Q-e0)) and proceed. At the final
continue outcome retain q as a dummy record, with the unmeasured clock
kept in its output environment and score zero. Blank supply is explicit;
no measurement is assumed to reset a physical register.

Take Delta=1/256 and m=2048. With r=256/257 and w_j=(1-r)r^j the ideal
proper-record probabilities are

    p(q,j,b)=w_j [M_q+b Re<a_q g,exp(-ij Delta(Q-e0))a_q g>]/2,
    p(q,dummy)=r^m M_q,                M_q=1/2.                 (1)

Thus there are 8194 proper terminal records. The other raw charge and
clock labels are explicit zero maps: each six-level pulse has exported
its complement before the successful code readout. The controller has
the same declared continuation menu for all six charge labels, including
the four zero prefixes. There are 36m+6 raw terminal pointer addresses,
of which only 4m+2 can be nonzero. Every first-export failure from
preparation, charge, clock, selector or readout is an additional absorbing
terminal record. Continue labels remain in the prefix of each terminal
history. The dummy output retains its clock environment; its quantum
output is not identified with an unevolved dummy CAR map.

The source instrument is complete on all inputs and arbitrary untouched
references. Each native source pulse includes its full six-level unitary
and every complementary export. Each selector separates stop, continue
and failure. Composition of these complete instruments, followed by the
terminal dummy branch, proves trace preservation without conditioning on
success. The finite execution contracts all sixteen source rows and both
classical charge columns. It uses neither number nor spectral truncation.

## Native clock and complete error

Expand the actual Q into its seven rational Pauli terms. Subtract the
rational midpoint e_ref of the certified e0 interval. Controlled terms
are h_a P_a tensor(I-Z)/2; the identity term and its relative phase are
retained. A symmetric product of their exponentials gives a 404-pulse
native step. Put M=sum_a|h_a|<37/32. The symmetric product and exact
exponential have equal derivatives through order two at zero; their
third derivatives each have norm at most M^3. The unitary Taylor integral
therefore bounds their difference by Delta^3 M^3/3.

The existing rational source-port and literal-angle certificates add the
native error epsilon_n, and the ground reference adds Delta delta_E.
Thus one complete controlled step differs by at most

    epsilon_step=epsilon_n+Delta^3 M^3/3+Delta delta_E.          (2)

The certificate bounds the complete export isometry, including failure
amplitudes, rather than only its successful matrix. The supplied source
preparation error epsilon_p includes its exact integer-direction error;
epsilon_q, epsilon_H and epsilon_coin denote the corresponding complete
charge, Hadamard and selector errors. A half-trace bound for the entire
terminal record/source/reference output is

    epsilon_total <= epsilon_p+epsilon_q+epsilon_H
        +(sum_(j=0)^(m-1) r^j) epsilon_coin
        +(sum_(j=0)^(m-1) r^(j+1)) epsilon_step
        +(1-r^m) epsilon_H.                                  (3)

The ground comparison in (3) concerns the specified blank preparation.
The subsequent channel bounds allow arbitrary untouched references.

To prove (3), telescope with the ideal prefix and actual suffix at each
intervention. Complete channels contract trace distance. The ideal prefix
reaches the jth selector with weight r^j and its following clock step with
weight r^(j+1); stopping readouts have weights w_j. This is a weighted sum
of unnormalized blocks. No bound is divided by an observed branch weight.
The clock sum is less than 256 and the selector sum is less than 257,
although the protocol retains all 2048 time addresses. Exact rational
evaluation of the certificates in (2)-(3) gives epsilon_total<1/100000.

Every ideal pulse has no exported amplitude. The squared norm of the
actual failure component of a complete native word is therefore bounded
by its squared native isometry error. Summing over all visits, even using
the coarse bound one for each prefix weight, gives

    p_failure <= epsilon_p^2+epsilon_q^2+epsilon_H^2
        +m(epsilon_coin^2+epsilon_n^2+epsilon_H^2) < 10^(-20). (5)

The unitary product-formula error is excluded from epsilon_n in (5).
The two charge columns have total weight at most one, so no extra factor
two is needed. This is a rigorous upper bound; the smaller floating
failure weights are diagnostics rather than certified lower values.

The implementation checks every native probability and all absorbing
failure slots. Its working source/clock register has 32 rows. Failure
weights are retained by their stage, charge, time, pulse and physical
complement slot; unused padding slots are distinguished from outcomes.
The source program and unnormalized prefix reconstruct each failure
output vector for this prepared input. The complete index domain has
17,449,226 physical first-failure addresses. Of these, 5,816,554 have one
of the two potentially nonzero charge prefixes or precede charge
readout; all others are exact zero vectors. Keeping the 2049 source/clock
prefixes uses about 2 MiB. No dense failure archive is allocated. The
longest branch has 829,504 native pulse applications; repeated local
blocks are compiled from the literal six-level matrices and reused.

## Original spectral score and finite boundary

The signed sum of (1) approximates

    Q_1=sum_q integral (1+E^2)^(-1) dmu_q(E).

The exact source local second moment is 65/256. The inherited global
geometric-time estimate therefore gives

    |Q_1-ideal_score| <= (Delta/2)(9/16)+r^m.                  (4)

A complete-output half-trace error epsilon_total adds 2epsilon_total to
the signed score error. All source failures and dummy outcomes have score
zero. The lower bound on original weight in [0,1] is
max(0,2Q_lower-1); it includes any zero-energy contribution. The exact
finite spectrum independently checks the measured score. This check is
not used to remove any detector outcome or to select a fitted energy.

The rigorous score intervals do not rely on floating contractions or on
their small normalization residual. Enclose sqrt(17) by consecutive
90-bit rational endpoints and use the inherited rational sine/cosine
Taylor bound to enclose each one-step source spectral phase. Store each
complex recurrence at 60-bit dyadic precision. If its one-step error is
e and its stored modulus is at most two, its jth phase error is at most
j(2e+2^(-60)): each integer rounding contributes at most 2^(-60) in
complex norm, and the unit-modulus exact factor does not amplify the
previous error. All stored moduli are checked by integer arithmetic.
Exact rational interval products with the source spectral weights and
geometric weights enclose every ideal record and the signed ideal sum.
Equation (3) then encloses every actual native record and its signed sum.
The floating execution is checked against these independently certified
intervals; it is not a premise in their derivation.

The actual native signed-score enclosure is

    [0.799269041434, 0.799309041435].

Allowing the native-to-ideal difference again when inferring the original
spectral score from this native enclosure, together with (4), gives

    Q_1 in [0.797809677050, 0.800768405819],
    mu_add([0,1])+mu_remove([0,1]) >= 0.5956193541.

These are original probabilities. Both charge branches, every sign,
dummy, continue prefix, raw zero label and native service failure are
retained. The numerical host calculation evaluates the specified native
source instrument; no experimental hardware shots are asserted. The
exact interval and channel arguments certify its mathematical output
law independently of the host's floating arithmetic.

The executed source is a dimer, not the selected infinite degree-six
ground. The 8194-record thermodynamic planner's collar and its pressure
search are separate constructions and are not executed by this program.
The native pulse calculation supplies no continuum shell, physical
species or mass identification, and no hardware-shot sample is claimed.
