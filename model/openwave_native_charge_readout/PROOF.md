# Source-compiled charge feedback and complete records

## Shared event and memory service

The carrier and original C/E record grammar are those of
`openwave_coherent_native_architecture`. The source is encoded in nine bits,
with all 512 addresses retained. A three-bit memory contains six cyclic
addresses and two padding addresses. The active charge law is
`openwave_active_native_charge`. No physical unit or particle identification
is added by this encoding.

For each seam and committed current, the post-source feedback is a full
permutation V_r of all 4096 source/memory addresses. On a valid moving pair it
adds the declared voltage to the six-state memory. It fixes source padding,
collision addresses, both memory padding addresses and source configurations
incompatible with the stated current. The holding service is the identity.
The compiler decomposes each permutation into Gray transpositions and
reversible Toffoli ladders. It executes the actual proper-code source words
for its reusable flip and Toffoli routines, including scalar phases. Exact
integer execution of every emitted gate checks all 4096 admitted columns.
Equal permutations may share computation within this invocation; their event
and pulse costs are charged at every use. No stored scientific result is an
execution input.

The source primitives are the complete six-level exponentials of the inherited
native control generators. Their exact proper-code identities imply the
logical reversible gates. The binary64 executions are diagnostics of those
identities, with a fixed routine admission ceiling 10^-10. If each routine has
operator error at most epsilon and L routines execute, the product error is
at most (1+epsilon)^L-1, bounded by L epsilon/(1-L epsilon) when L epsilon<1.
The displayed full-program bound uses this implication; it is neither an
outward certificate for the numerical eigensolver nor a physical hardware
precision claim. No polar projection is used.

All scratch registers are tested for the declared blank initialization.
The success effect is P_blank and the separate failure effect is I-P_blank.
On success, the exact emitted circuit restores the scratch and has action
V_r. On failure the instrument retains the input and its failure record.
Thus arbitrary dirty scratch is not projected away or asserted to implement
the same logical map. Complete active operator identities extend to arbitrary
quantum references. A dense 21-qubit floating matrix is not evaluated.
The smaller 16-address source-code/memory witness is executed numerically on
every clean and dirty scratch column, giving an independent full-operator
check of the coupled phase experiment.

For rectangular code load/export maps, charge compatibility is the typed
intertwiner Q_out K=K Q_in, including identity factors for pointers and scratch.
For an endomorphism this is commutation with Q. Computational character
initialization is a separate supplied preparation operation: its pulses need
not conserve the charge of an isolated memory. A charge-conserving physical
initialization would require a specified reservoir or a neutral joint
preparation. The conservation theorem concerns the declared operations after
initialization.

## Actual coherent source preparation and complete quadratures

Let a=(0,1), b=(0,6), a'=(1,5), b'=(5,6). The native C moving record on seam
(0,5) has probability 1/2 and maps both input alternatives to their primed
partners. Its memory voltages are five and zero. Prepare the source state
(|a>+|b>)/sqrt(2) and a memory character |k>, with

    X6 |k> = exp(2 pi i k/6) |k>.

The complete eight-address source-word preparation of each character is
executed, starting from a tested blank memory. All seven nonblank input
addresses have a retained preparation failure. The source superposition is
also prepared by its actual 512-address native word. No ideal preparation
matrix replaces the source execution.

On the moving record the normalized source is

    (exp(5 pi i k/3)|a'> + |b'>)/sqrt(2).

Let P project onto the two primed states, X=|a'><b'|+|b'><a'| and
Y=i(|a'><b'|-|b'><a'|). The complete five-effect readout is

    (P+X)/4, (P-X)/4, (P+Y)/4, (P-Y)/4, I-P.

A fair classical quadrature choice, a code-complement test and the actual
inverse native eigenstate-preparation words implement these effects. All
512 source output addresses are retained. On the admitted code the plus
outcome is computational zero after the inverse word; the minus outcome is
its full complement. This measured complement is distinct from the initial
source-code failure. The calculation checks both source code columns through
the actual X and Y readout words, and compares their complete probabilities
with the independently executed smaller witness.

The four joint probabilities, including the moving-record probability, are
(1±cos(5 pi k/3))/8 and (1±sin(5 pi k/3))/8. The other native holding outcome
has probability 1/2. The receipt retains all 25 native slots times five
readout outcomes, including zeros and code failure. An initial projector admits exactly the source pair, six active memory addresses and blank scratch; its complement is a separate failure record before any native attempt. The moving output has a second complete memory/scratch admission. Actual source-preparation leakage and every numerical padding amplitude enter these failure records. The holding weight is half the actual admitted norm. Preparation and subsequent numerical norm defects are reported separately, without renormalizing any state or record. These are unconditional
joint probabilities. No successful branch is normalized into an unconditional
claim. The six phase pairs distinguish the six memory characters. Source
coherence within a single native raw branch suffices; different committed
histories are not recombined or erased.

For a concrete interaction statistic, prepare k=0 or k=2 with equal
probabilities and keep label C=+1 or -1. Score the moving X+ and X- outcomes
as +2 and -2, and every other outcome as zero. Then

    E[C]=0, E[J]=1/8, Cov(C,J)=3/8.

Removing the declared order-three voltage while retaining the transported
ordering parity gives covariance zero. This is a complete source-record
backreaction test of the active memory. The electric-charge unit, localization
of a particle core and a rest-energy dictionary are not inferred.

## Exact compatibility with the native clock recovery

The six logical clock labels use the actual six matching supports of the
512-address continuation. In a raw branch let c be the admitted support size
and chi_j its memory voltages. In charge sector k, the recovered code amplitude
is proportional to

    sum_j exp(2 pi i k chi_j/6) / sqrt(6c).

The original complete matching mixer gives this formula on all admitted
columns. Consequently the branch leakage is c/6 minus the squared amplitude,
multiplied by the original seam and coin weights. The executor evaluates all
six charge sectors, logical labels, raw records and the identity timeout in
exact arithmetic. For this voltage all code-population leakage is zero.
The phases between distinct logical labels can nevertheless depend on the
record. Summing the complete records yields exact coherence multipliers in
Q(i sqrt(3)). In sector k=1 under C, the (0,2) logical coherence multiplier is
85/96. The loss 11/96 prevents transfer of the uncompensated clock theorem,
even though the code populations and memory charge are preserved.

There is a complete compensated observation policy within the same hardware.
After each original raw record, apply the actual feedback inverse V_r^dagger
before the original source recovery R_r. Then, branch by branch,

    (R_r tensor I6) V_r^dagger V_r (K_r tensor I6)
        = (R_r K_r) tensor I6.

The emitted forward and reversed reversible programs are checked on every
4096 admitted input, including source and memory padding. Blank scratch is
restored, and the dirty-scratch failure is retained. The identity extends to
all references and preserves the original raw transcript. Source-only H2
phases commute with the memory, so the complete compensated clock protocol
has the original channel tensored with the unchanged memory. Its error bound
and unit ancestry transfer only for this compensated policy, with two feedback
invocations and their complete pulse/service costs per event. Approximate
physical gates require their own error bounds.

The covariance 3/8 is observed before inverse feedback. Undoing that feedback
also undoes its phase signature. This construction does not establish a
simultaneous charge-dependent rest-mass clock or a physical particle spectrum.
