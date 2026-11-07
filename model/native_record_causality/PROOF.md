# Independent repairs and complete records

## Domain and observation map

Use the literal `native_program.interpreter` on either the full 66-dimensional
binary total-two carrier or the full 78-dimensional nonnegative total-two
carrier. Commands have fixed native seams and fixed current/contrast record
laws, or the declared mixture. All port weights, repair probabilities and
mixture weights are parameters of this same admitted program. The selector
does not depend on a load, a previous record or the command position.
Each command retains its hold outcome, both coins, all eight current-record
addresses and all sixteen contrast-record addresses. Addresses with zero
operator or zero mixture weight are retained. The total alphabet has 25 slots.

Consider two fixed commands with disjoint endpoint sets. They execute in the
orders `AB` and `BA`, starting with empty record memory. Both executions must
finish. The finite memory either has enough space or exports its oldest
record; exported records remain part of the observed output. Resource stops
and intrinsic memory timeouts are outside this domain.

The comparison concatenates the exported archive and the final memory,
interchanges the two command-instance slots, and splits the tape using the
same capacity. This is a bijection with its own inverse. It preserves every
raw address, including holds and reserved zero/abort addresses. Slot identity
comes from the declared program, so a hold need not encode its seam. The
source and any untouched reference are unchanged by this comparison. Known
program identifiers are mapped from `BA` to `AB`; their literal hashes are
different. The comparison is a relabeling of completed outputs, with no
additional repair operation or preparation.

## Exact channel identity

For a slot `a` of command `A`, write its native Kraus operator as
\(K_a=\sqrt{w_a}V_a\). The nonnegative scalar \(w_a\) is fixed by the
declared repair probability, record-law mixture and endpoint coin bias.
The matrix \(V_a\) is the actual partial permutation on the full carrier.
Hold has \(V_a=I\); reserved addresses have \(V_a=0\).

Membership in the domain of \(V_a\), its recorded current or contrast, and
its update all depend only on the endpoints of `A`. A disjoint repair `B`
leaves those coordinates unchanged. Hence the two domain conditions commute
and the two endpoint updates commute, including when either condition fails:

\[
V_bV_a=V_aV_b,\qquad K_bK_a=K_aK_b.
\]

Both updates preserve the fixed total, so restriction to either complete
carrier preserves this identity. It holds for all 315 unordered disjoint
seam pairs and all \(25^2\) ordered raw-slot choices. The certificate checks
all 196,875 operator pairs on each carrier exactly, without a floating
tolerance or a selected invariant subspace.

For every input operator \(X\) on source and any finite reference `R`, the
`AB` output block is

\[
(K_bK_a\otimes I_R)X(K_bK_a\otimes I_R)^*.
\]

The corresponding `BA` block, with its records permuted, is identical.
Thus the complete classical-quantum channels are equal after the stated
bijective comparison. The diamond distance between those compared channels
is exactly zero. The proof includes all coherences, arbitrary mixed or
entangled preparations, and arbitrary untouched references. Real and
imaginary parts may be evaluated separately because the interpreter acts
complex linearly with real rational coefficients. No conditioning on a
successful outcome or division by a branch probability occurs.

An arbitrary final source/reference effect accompanied by a record test has
the same expectation when that test is pulled back through the record
permutation. A chronological record test without this pullback can separate
the two runs. The code checks that their literal chronological outputs
differ on the executed entangled preparations.

## Finite fixed schedules

Assign an instance label to each command of a finite fixed schedule. Swapping
adjacent commands with disjoint endpoints permutes their complete output
slots and leaves the corresponding product of quantum operators unchanged.
Repeated adjacent swaps prove the same statement for any two fixed schedules
related in this way. This is a finite partial-commutation equivalence of
observation maps. The complete tape, command-instance labels and all scalar
weights travel with the commands. The induction requires no execution to
halt and no selector or later operation to consume a reordered record.

The single memory queue is written in chronological order in every literal
run. Its full implementation provenance has storage edges between attempts.
The comparison theorem does not remove those edges from the implementation
and does not apply the disjoint-footprint theorem to overlapping memory
writes. It establishes a separate equivalence of completed instruments.

## Exact boundaries

The actual adaptive interpreter supplies a control signal between commands
with disjoint endpoint sets: two source preparations give receiver
probabilities \(1/4\) and \(0\). Its record-dependent selector cannot be
discarded in a physical causal reconstruction.

For two fixed disjoint commands on `(0,5)` and `(1,7)`, capacity one and the
intrinsic `halt` policy execute only the first repair before the timeout.
Starting with occupation `(0,1)`, the final occupation probability at port
five is \(8/33\) or \(0\), depending on order. Both timeout outcomes remain
in the output. For fixed overlapping seams `(0,5)` and `(0,1)`, starting
with occupation `(0,2)`, the complete final occupation probability at port
five is \(8/33\) or \(4/33\). These controls separate failure of the
endpoint hypothesis from failure of complete execution.

## Spatial and field interpretation

The construction supplies a source-bound equivalence for a specified family
of finite observation protocols. It does not select which observers or
record effects define physical events. It assigns no spatial coordinates,
dimension, metric, light cone, physical clock, field algebra or particle
energy. In particular, commutation of independent finite operations alone
does not imply three spatial dimensions, Lorentz symmetry or a relativistic
dispersion relation.

A physical use requires a common source family with an operational choice of
local record algebras and spacetime readout, compatible refinement and clock
maps, and controlled limits of the same prepared-state correlations and
responses. A Lorentz null-ray reconstruction can then constrain an actual
transition map if that map satisfies its null-preservation premises. Finite
gauge covariance constrains a supplied action but does not fix that action,
its couplings, a chiral continuum or particle poles. No such premises are
discharged by this finite comparison alone.
