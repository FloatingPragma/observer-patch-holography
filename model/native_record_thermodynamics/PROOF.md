# A native record thermometer and conditional repair

This construction acts on one observer's finite record register. It uses
actual native repair attempts to produce coin-prefix records, then a
specified comparison controller to refresh that register. It is not a
claim that the unchanged occupation process has the same stationary law.
Its energy is an operational prefix cost in dimensionless units, not total
dissipated work or a laboratory Hamiltonian.

## A nonconstant cost fixed before the reference

On a fixed committed seam, execute one fair native placement attempt and
read its coin. Stop for coin zero. After coin one, execute one more attempt
and read its coin. The three terminal public prefixes are

    h_0=0,       h_1=10,       h_2=11.

Every coin is part of an actual native endpoint/current record, including
no-movement attempts. Their lengths are n=(1,2,2). Charge each attempted
placement `log(2)` units and define `E(h)=n(h)log(2)` before selecting any
reference. Fair-coin path counting gives the prefix probabilities. The label
reads only coins, so they are the same on all 78 total-two inputs by construction:

    w=(1/2,1/4,1/4)=exp(-E).

The lengths therefore also equal record surprisal in these units. The
seam choice and stopping rule are declared. The probability law relies on
independent fair native coins; deterministic repair alone does not choose
that probability measure.

Generate b prefixes consecutively, where b is a supplied positive integer.
Independent coins make these labels independent even on one continuing
carrier: the conditional next-prefix law is w for every source input and
every preceding record. No independently prepared carrier replicas or
source reset are required. Retain all fine records and compare the b labels.
An agreement on h has probability

    a_b(h)=w(h)^b,       Z_b=sum_h a_b(h)=2^(-b)+2*4^(-b).

On agreement the register is replaced by h. On disagreement it is kept.
Thus comparison selects, rather than assumes, the faithful reference

    pi_b(h)=a_b(h)/Z_b=exp(-b E(h))/Z_b.

The same nonconstant E is used for every b. The relation
`pi_b(0)/pi_b(1)=2^b` identifies b within the implemented family. This is
an internal thermometer for a controller-selected parameter, not a
contact-dynamics theorem or the selection of a physical temperature.

## The actual finite register kernel and the exact heat bath

Let `P_b(x,y)=pi_b(y)`. One completed comparison proposal has register kernel

    Q_b=(1-Z_b) I+Z_b P_b.

All disagreement outcomes and their source tapes are retained. Both kernels
are stochastic, preserve pi_b and satisfy detailed balance. The kernel
P_b is precisely the conditional-resampling kernel for the constant visible
datum. Its rows are information projections onto that whole register fibre.
The same pi_b is the state-side reference optimizer with no extra moment
constraints. This shared reference comes from the counted source words.

The controller tries at most M independent proposals and stops at the first
agreement. If every proposal disagrees, the register is kept and a timeout
is reported. Its unread register kernel is

    Q_(b,M)=(1-Z_b)^M I+[1-(1-Z_b)^M] P_b=Q_b^M.

Its exact timeout probability is `(1-Z_b)^M`, independent of the incoming
register and source. Conditional on agreement, the accepted register law is
pi_b, also independent of the incoming register. The complete flagged
comparison with an ideal always-successful P_b has total variation at most
the timeout probability. Forgetting the tape is a marginal, not a physical
memory-erasure operation. The growing archive is not claimed stationary.

For b=1, Z_b=1: one prefix proposal implements the exact nontrivial heat bath
P_1 in at most two native attempts. For b>1 and finite M, the actual kernel
is generally not idempotent; its nonconstant eigenvalue is `(1-Z_b)^M`.
The exact P_b is realized by retrying until agreement, almost surely, with
finite expected cost. It is not a uniformly bounded executor. In particular,
pi_2=(2/3,1/6,1/6) has a nondyadic probability, so no bounded fair-coin tree
can implement its exact unconditional sampler. The finite implementation
keeps its timeout instead of claiming an exact projector.

## Native resource costs and failures

One prefix costs one or two attempts with expectation 3/2. Each proposal
therefore costs at most 2b attempts and has expectation 3b/2. Sequential
renewal, including the cost of every rejected proposal, gives

    expected attempts through M proposals
      = (3b/2) [1-(1-Z_b)^M]/Z_b,
    expected attempts until agreement = (3b/2)/Z_b.

The cap is 2bM attempts. Agreement need not be independent of proposal
length; the expectation uses independence of a fresh proposal from the
previous failure history, which is sufficient. All records have a finite
maximum length at fixed b and M. Exact unlimited retry has no finite worst
case record-memory capacity. Any fixed memory implementation uses the
explicit cap and timeout clause.

Perfectly copying one prefix b times instead of generating fresh records
has success one and accepted law w. For b>1 this differs from pi_b. This
control makes the independence requirement testable. Changing the stopping
rule changes E and the family; neither is selected by an unqualified
statement of deterministic repair.

## Instantiating the finite thermodynamic core

For each b, `RepairLawData` consists of the finite three-state register,
reference pi_b, and constant visible datum. Set calibration energy E,
positive multiplier b, and `logZ=log Z_b`. Directly from the source word law,

    -log pi_b(h)=b E(h)+log Z_b.

These are the two mathematical inputs to the existing finite four-law core.
The native bounded implementation supplies its exact P_1 instance. The
other P_b instances have the stated almost-sure retry implementation, while
the bounded Q_(b,M) inherits stationarity, detailed balance, relative-entropy
contraction and supported fluctuation relations, rather than projection-
specific exact dissipation. The conditional core is a mathematical theorem;
calling its calibration field a physical energy requires additional evidence.

For P_b, with q=p P_b=pi_b, the exact entropy split is

    H(q)-H(p)=b (<E>_q-<E>_p)+D(p||q).

For the bounded Q_(b,M), the same Gibbs reference and data processing give
`H(q)-H(p) >= b (<E>_q-<E>_p)`. A fixed-energy register refresh has the exact
energy bookkeeping `Delta U=Delta Q`, with work zero by definition of this
model stroke. This defines a dimensionless register-energy balance; the
actual source apparatus, complete tape and memory-reset costs are separate.

The ground label is h_0, the gap is log(2), and the excited reference mass is

    q_exc(b)=2/(2^b+2) <= 2^(1-b).

Every finite b has full support, whereas the b-to-infinity limit is the
ground point mass. The finite Gibbs entropy tends to zero by this bound.
A finite number of these positive-reference refreshes cannot turn a
strictly positive register law into a rank-deficient ground state. This is
a finite three-state statement, with no cofinal spacetime or physical
cooling conclusion. No infinite-temperature continuum, laboratory units,
physical Landauer cost, or occupation-spectrum calibration is inferred.
