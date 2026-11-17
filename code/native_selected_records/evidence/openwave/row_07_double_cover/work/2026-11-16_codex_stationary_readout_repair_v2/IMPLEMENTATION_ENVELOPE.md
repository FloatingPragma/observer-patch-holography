# Finite precision and finite sign assay

The nominal v5 program already uses finite rational matching phases. The exact
nearest-cell certificates include selector leakage and every source-dependent
norm. Those coefficient approximations are part of the program evaluated by
the independent source certificate. They are not subsequently treated as
ideal unit-gain filters, and finite swap dark-port leakage is retained.

There is a separate implementation question: how much can a realization of
that specified program differ from it? Give every specified source/inverse
access and coherent memory service an operator-norm error at most delta on
its complete finite workspace. Give the admitted initial state trace error
at most delta. This is an explicit implementation hypothesis, not an empirical
hardware measurement. A service is the complete declared address, phase,
Hadamard, current-add, controller or source operation; it is not an uncounted
physical elementary gate. The finite source extension and all orthogonal
workspace outcomes make the norm statement meaningful.

The cap64 observed calibration, cap16 marker retries, finite counters and
chosen finite amplifier words give a uniform upper bound N=2^22 services per
trial, including one chosen Ramsey setting at the largest retained stride.
`error_budget.py` itemizes a stricter case-specific bound. Prefix failures
cannot increase it. Multiple tomography settings require separate prepared
copies; their complete batch cost is not hidden in this per-probe bound.

Telescoping unitary/Stinespring products gives output trace error

    epsilon <= delta + 2 N delta.

Adaptive readout/controller maps and retaining or discarding a specified
outcome are contractions. If the ideal event has probability p>=p0, its
implemented probability is at least p0-epsilon. The normalized conditional
state differs in trace norm by at most2 epsilon/p0. Therefore a conditional
binary probability changes by at most epsilon/p0, and an interference
visibility by at most2 epsilon/p0. This bound is applied to the **complete
source before conditioning**, including failures, rather than to a normalized
accepted map with its cost forgotten.

For row7 use p0=.03 and delta=2^-40. This yields a nonzero finite implementation
interval far smaller than the primary sign/return margins, subject to the
independent outward certificate establishing the nominal thresholds. The
control's nominally tiny selected-minus continuation is not guaranteed to
remain nonzero under that envelope. Its nearly+1 both-port visibility remains
the relevant robust control observable. Neither a rare control conditional
clock nor the old tiny stationary fringes receives a noise guarantee by fiat.

For a fully joint two-copy readout with arbitrary bounded implementation
errors, use the whole joint acceptance p0² and twice the service count. Do not
silently assume independent noise or propagate only single-copy marginals.
`envelope(bits=52,copies=2)` supplies a stricter example whose binary error is
below1e-5. This is separate from the source impurity certificate and the fixed
pair-reader envelope. It does not claim a finite-shot experimental mass
measurement; obtaining very precise empirical probabilities has its own cost.

For the narrower sign assay,64 accepted both-port records per condition and
8192 admitted raw attempts per condition suffice for the displayed Hoeffding
and binomial-shortage bound over224 primary conditions. The reader records all herald
failures and both interference outcomes, and decides sign without a law label,
seed, filename or expected phase. This is a finite implementable assay with
an explicit workload bound. It is not a claim that these laboratory samples
have already been collected. The complete source run enumerates/model-computes
the laws that the assay would sample.
