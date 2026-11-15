# A moving native code with an exact stopping boundary

CP-03 / [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139).
This construction removes the fixed-core reset from the actively
corrected qubit. Its classical core follows the complete native current,
and each successful correction encodes the unknown logical state at that
destination core. An exact calculation also shows where this construction must stop.
It is moving, actively serviced memory; autonomous particle protection and
a reference-like exterior are not achieved.

The baseline is [PR #1151](https://github.com/FloatingPragma/observer-patch-holography/pull/1151)
at `67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c`, used provisionally for this
follow-up. The fixed source weights are `(25,25,8,...,8)`, exactly those of
its record-derived clock. The recovery controller does not need to know or
infer these weights: they multiply complete flags without changing their
code geometry. They remain parameters of the native source law. The added
active operations change the complete process, so its marginal is not
claimed to retain PR #1151's stationary reference or preparation guarantee.

The [contract](campaign.json), [implementation](transport.py),
[retained results](results.json) and [independent tests](test_transport.py)
specify this finite construction. Neither CP-03 nor a dependent gate is
accepted by this packet.

## Compatibility and nearest existing theorem

For each nonopposite two-port core `c`, use the existing
[sharp-q=1 code](../particle_screening/README.md). There are 30 edge cores
with four safe exterior ports and 30 distance-two cores with two. The two
logical columns are disjoint, equally weighted crossed rectangles; every
source pair is a nonedge. The admissible state is an arbitrary logical
qubit, including an untouched entangled reference, encoded by `V_c`.
Compatibility is the independently specified condition `q(x,c)=1`, with
the stated full native flag menu and supplied correction service.

The nearest result is the existing
[all-seam flagged-code theorem](../common_source/ENCODED_MEMORY.md#optimal-native-carrier-code).
It supplies the rank and scalar-compression argument on this exact carrier;
the preceding [compiled fixed-core service](../screened_recovery/README.md)
implements a special case but resets the label. This packet checks the
actual native maps at **all 60 cores and 43,200 raw C/E slots**, and changes
the correction target to the core determined by the recorded current.

Let `h` include seam, placement coin and full raw register. If its two
native image columns contain `m_h` terms each, and `V_c` has `n_c` terms per
column, then

```text
V_c* K_h* K_h V_c = p_h I_2,
p_h = P(coin | seam) m_h/n_c.
```

The two images are disjoint; the flag determines a unique transported
core `c'`. When `c'` is nonopposite, an emitted full-register permutation
decodes these images, the inverse `m_h` mixer removes the syndrome, and
the emitted `c'` encoder gives

```text
R_(c,h) K_h V_c = sqrt(p_h) V_c'.
```

Every unknown logical/reference state is therefore preserved in this
changing encoding, with flag probabilities independent of its amplitudes.
The identity composes for any finite history-controlled seam/C/E schedule
up to the declared stop, under ideal supplied services. Holds preserve it.
There is no successful-step reanchoring and no quantum-state normalization.
The position label is classical: this does not construct a coherent
superposition of positions or a dispersion relation. The supplied recovery
can rearrange source occupations away from the just-repaired seam. Its
source compilation does not prove locality in the port graph or an
autonomous propagating excitation.

## The boundary is reached with probability one

Six opposite core pairs have no safe exterior ports. The existing
all-seam theorem forces an exact correctable qubit to avoid native-edge
source pairs and to compress every one-port occupation to a scalar. For
an opposite core, every exterior port leaves at most one allowed core
partner. Its occupation projector then has rank at most one on the sharp
sector. Its compression to a rank-two code must vanish. All exterior
occupations would vanish, contradicting the required exterior total one.
This excludes **every exact all-seam sharp-q=1 qubit at that core**, not
merely the chosen crossed rectangles.

The moving service therefore retains the native output on an opposite-core
flag, records a stop and halts. This output preserves the logical
state through a known isometry, but it is not certified for another native
attempt. No stopped output is discarded or sent back to the initial code.
There are **11,460 positive continuing flags and 120 positive stopping
flags** across the two laws; all raw abort and valid zero slots remain.

With independent uniform seam choice and the fixed biased endpoint law,
the complete core kernel has 60 live rows and 66 target columns. C and E
have exactly the same kernel, so any mixture of those two laws has the
same core marginal. At `(0,1)` the core moves with probability **16/495**
per attempt. A distance-two core can reach an opposite pair on the next
attempt. The minimum probability of stopping within two attempts, over
all live initial cores, is **49/245025**. Thus

```text
Pr(T > 2n) <= (1 - 49/245025)^n.
```

There are no closed nonstopping classes, and eventual stopping has
probability one. The exact solution of `(I-Q)t=1`, checked in the original
60 equations, gives the initial mean

```text
97909756773858989776401485931436511086245
------------------------------------------------- = 203.025378... attempts.
482253783605192403282205691964646745488
```

The fair-coin control independently lumps by edge/distance-two/opposite
class, with live rows `(29/30,1/30,0)` and `(1/30,19/20,1/60)`. Its mean
is exactly **150 attempts** from an edge and **120** from distance two.
These are native attempt counts, not modular parameters or physical lifetimes.
Different adaptive seam schedules need their own stopping analysis. The
general branch identities do not by themselves imply this uniform-scheduler
absorption bound.

## The actual reference and clock constrain the candidate

This packet keeps the same **fixed** product-pair reference as PR #1151.
It never moves the heavy weights along with the core. For every sharp-q=1
state at a specified core, measuring the full exterior particle number gives

```text
D(rho_exterior, pi_exterior) >= 1 - pi(q=1 | core).
```

The conditioned diagonal reference attains this bound. It is therefore an
exact minimum over all quantum states in the sharp sector. At the initial
core the minimum is **701/1501**, while the chosen code's worst single
exterior-port distance is **976/7505**. The latter is a different, weaker
observation. Moving the core has not produced a reference-like full exterior.

Only **19 of the 60** code subspaces are invariant under every parameter
of the fixed reference's modular group. Every logical column uses pair
weights 64 and/or 200, and both columns have the same weight distribution.
If `r` is the fraction of weight-200 terms, applying the modular flow for
dimensionless parameter `s` has exact code leakage

```text
1 - || V_c* U_s V_c psi ||^2
  = 4 r(1-r) sin^2(s log(25/8)/2),   ||psi||=1.
```

A scalar compressed clock is not a preserved code subspace when this
leakage is nonzero. Better record inference alone cannot supply the missing
clock/code compatibility. The modular service and moving recovery need a
joint protocol before their capabilities can be composed. No rest
energy, stress tensor, spin/statistics or physical particle is assigned.

## Source execution and cost

All 60 encoders and all positive continuing corrections have emitted source
schedules. An exact integer execution checks **60 forward and 1,550 distinct
inverse Gray schedules on all 512 addresses**, including padding. The
two-/three-/four-term mixers and source primitives reuse the existing
compiler. This is an exact layout check, not a claim that every long word
was replayed numerically or physically implemented.

A fresh floating source replay executes all 60 encoders and every raw flag
on a C-then-E history with a complex full logical basis (equivalently a
maximally entangled reference after normalization). The witness moves
`(0,1) -> (1,10) -> (9,10)`, where the last core causes the explicit stop.
Its exact flag probability conditioned on the two selected seams is
`16/1089`; it is not a survival-conditioned experiment. Source branch
amplitudes and total output mass must pass fixed `1e-8` numerical admission
ceilings. The full raw menu is compared before selecting a witness flag.
Fresh numerical admission is not a directed-rounding or hardware-fault
certificate, and is not a floating replay of all 43,200 cases.

The initial encoder costs **65,484 source pulses**. The initial-core mean
native-plus-recovery cost under uniform seams is `1240015193/1980` pulses
for C and `405163927/660` for E. The maximum across all live cores, laws and
positive flags is **864,887 pulses per attempt**. Stopping flags incur the
native word but no recovery. The receipt also solves the exact absorbing
reward equations for cumulative source pulses through the stopping attempt,
with the initial encoder added separately; it does not discard failed work.
A fixed archive uses 24 C or 25 E bits per
attempt, including law, seam, coin, raw flag, old/new core and stop status.
The nine-qubit code, up to three raw-record qubits, six-level processor,
coherent load/export, blank preparation, storage and classical control
retain their supplied apparatus contracts. Cached host computations never
remove per-invocation pulse or archive costs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.moving_defect.experiment --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/moving_defect
```

The source-bound receipt owns a fresh replay; it does not accept a supplied
summary of success. Independent tests reconstruct endpoint records and the
core kernel by separate arithmetic, verify the exact stopping equations,
check the exterior measurement and modular leakage, and reject damaged
Gray schedules, a fixed-core-reset substitute, nonfinite diagnostics and
altered receipt types or stop labels.

The resulting exclusion belongs to CP-12 as a **finite canonical moving
memory obstruction**. The opposite-core exact-code theorem has its wider
specified local scope. Approximate correction, other charges or cores,
different native operations, coherent position encodings and refining
carriers are outside this finite stopping calculation. No whole-program
exclusion or empirical comparison is added.
