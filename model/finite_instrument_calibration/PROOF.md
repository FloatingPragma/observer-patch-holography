# Selecting a complete instrument with finitely many probes

The finite theorem in `Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean`
selects the Lüders operation from its projective effect and invariance of
every outcome-certain state. The following calibration criterion replaces
those two conditions by finitely many joint preparation/readback probabilities.
It also gives a uniform error bound for arbitrary untouched references.
The quantum state, Kraus-map and probability semantics are hypotheses.

## Complete calibration protocol

Let `H` have dimension `n` and a trusted orthonormal basis `e_i`. Partition
that basis into nonempty blocks `I_y`, of ranks `d_y`, with orthogonal
projectors `P_y`. The candidate instrument has finite Kraus operators
`K_ya:H->H`, every outcome label `y` is retained, and

```
sum_(y,a) K_ya* K_ya = I.
```

Its induced effects are not assumed to equal `P_y`. For each basis input
`e_i`, record the joint event: the instrument reports its assigned label
`y(i)` and a trusted subsequent rank-one test returns to `e_i`. Its
probability is

```
p_i = sum_a |<e_i,K_(y(i),a)e_i>|^2.
```

For every block of rank greater than one, also prepare
`u_y=sum_(i in I_y)e_i/sqrt(d_y)` and test the joint event with the same
label and return to `u_y`. Write its probability as `p_y^+`. There are
`n + #{y:d_y>1}` distinct preparation settings. Rank-one blocks require
only their basis setting. These are unconditional joint probabilities:
wrong labels, failed readback, and every other outcome count as failures.

**Finite selection theorem.** All the listed probabilities equal one if
and only if every outcome map is

```
Phi_y(X) = P_y X P_y.
```

In particular the protocol determines both the induced effects and the
full state update. No independently supplied projective-effect equality or
invariance over an infinite set of states is needed for this implication.

For the forward direction, trace preservation and `p_i=1` exhaust the
total squared Kraus norm on each input `e_i`. Thus every wrong-label
operator annihilates that vector, and every correct-label image is
`K_ya e_i = c_(ya,i) e_i`. The coherent input has output entirely in its
own pure ray when `p_y^+=1`. Positivity of the Kraus sum therefore gives
`K_ya u_y=c_ya u_y` for each branch. Expanding in the orthonormal basis
shows `c_(ya,i)=c_ya` throughout that block. Hence `K_ya=c_ya P_y` and
`sum_a |c_ya|^2=1`. Substitution proves the displayed map identity.
The converse follows by direct projection. This is the finite Kraus
rank-one argument underlying the parent certain-state theorem, with an
explicit finite set of sufficient probes.

The coherent probe can be any trusted vector whose coefficients in the
block basis are all nonzero for exact selection. Uniform coefficients are
used for the quantitative bound below. More generally, several pure probes
suffice when their supports connect the basis indices. If they leave two
disconnected components, a relative phase between the components passes
every such test and changes the instrument. No minimum-setting claim is
made for arbitrary calibration designs.

## Stability with arbitrary references

Suppose certified probability bounds give
`p_i>=1-epsilon_i` and `p_y^+>=1-epsilon_y^+`, with all errors in `[0,1]`.
Let

```
E_y = sum_(i in I_y) epsilon_i,
eta = sqrt(sum_i epsilon_i)
      + sqrt(sum_(y:d_y>1) d_y (sqrt(epsilon_y^+)+sqrt(E_y))^2),
delta = min(1,2 eta).
```

For the complete labeled candidate channel `Phi` and complete Lüders
channel `L`,

```
(1/2)||Phi-L||_diamond <= delta.
```

No inverse outcome probability occurs. The bound includes terminal quantum
states and hence every subsequent quantum or classical readout. It is
conservative and is not asserted to be optimal.

To prove it, use a common Stinespring isometry

```
V e_i = sum_(y,a) |y> tensor K_ya e_i tensor |y,a>.
```

The environment carries the outcome label so tracing it produces a
classical direct sum, including every cross-input coherence inside each
outcome. Put `o_i=|y(i)> tensor e_i`, decompose
`V e_i=o_i tensor a_i+r_i`, and choose `r_i` orthogonal to the output
ray `o_i`. The squared norm of `r_i` is the actual basis failure
probability, at most `epsilon_i`. Each `a_i` belongs to the environment
subspace carrying `y(i)`. Define `W e_i=o_i tensor a_i`; then

```
||V-W||_op <= ||V-W||_F <= sqrt(sum_i epsilon_i).
```

For one block, project orthogonally to the output ray
`|y> tensor u_y`. On `V u_y` that residual has norm at most
`sqrt(epsilon_y^+)`. Replacing `V` by `W` costs at most `sqrt(E_y)`.
Writing `a_bar_y=(sum_(i in I_y)a_i)/d_y`, the exact orthogonal variance
identity gives

```
sum_(i in I_y)||a_i-a_bar_y||^2
 <= d_y (sqrt(epsilon_y^+)+sqrt(E_y))^2.
```

For a rank-one block its variance is identically zero. Replacing each
`a_i` by its block mean defines `W0`, with `||V-W0||_op<=eta`.
Each mean has norm at most one. Since `V e_i` has norm one,
`1-||a_bar_y||<=eta`. Normalize each nonzero mean to a unit `b_y`.
The replacement `W0->W1`, where `W1 e_i=o_i tensor b_y`, costs at most
`eta` in operator norm, because the output rays `o_i` are orthogonal.
The `b_y` for different labels are orthogonal. If a mean vanishes, then
`eta>=1`, and the universal channel bound `delta=1` already suffices.

Thus `W1` is a Stinespring isometry for the complete Lüders channel and
`||V-W1||_op<=2 eta`. For any pure input with an untouched reference,
the trace distance of the two isometric output states is bounded by their
vector distance. Partial trace contracts that distance. Purification and
convexity give the same bound for every reference/input dimension, proving
the half-diamond estimate. The universal bound of one supplies the cap.
The executable computes all square roots outward with exact rational
arithmetic; rounded decimal inputs are rejected.

## Complete histories and finite data

For at most `N` uses interleaved with the same arbitrary adaptive CPTP
controls, hybrid replacement and contractivity give

```
TV(complete record laws) <= min(1,epsilon_prep+N delta).
```

The stronger half trace-norm bound also holds for the complete terminal
classical/quantum outputs. This is a consequence of a uniform instrument
bound, so no independence between successive experimental states is needed.
The conditional uniform bound must hold for every instrument instance or
admitted history; a stationary calibration does not cover an unconstrained
history-dependent device. Conditioning on a successful outcome requires a
separate bound and is not used here.

The inputs `epsilon` are probability error bounds, not observed failure
fractions interpreted as exact probabilities. For example, under independent
identically distributed calibration trials, `N_c` trials with no failures
for each of `m` probes yield the simultaneous one-sided allowance
`1-(alpha/m)^(1/N_c)` at confidence at least `1-alpha`, by the binomial
zero-failure bound and a union bound. Finite zero-failure counts do not prove
zero error. Correlated calibration, untrusted preparation or untrusted
readback requires its own model and uncertainty contract.

## Exact controls and scope

The certificate checks maps independently on all matrix units. It includes
a rank-two block plus a rank-one block, the two-dimensional Pauli-Y context,
and the swap-twisted maps used by the parent instrument construction.
Dephasing within a rank-two block passes every basis test and fails the
coherent test. A swap inside the block retains the same effects and the
parent repeatability property, passes its uniform coherent test, and fails
the basis tests. Exchanging outcome labels leaves an unread channel
unchanged while failing the joint labeled protocol. These controls distinguish
effect agreement, repeatability, unread dynamics and complete instrument
selection.

The module `Lean/EventAlgebra/FiniteProbeInstrumentSelection.lean` proves the
full-block, single-outcome core. Coordinate pure-matrix fixation and uniform
pure-matrix fixation force every Kraus operator to be scalar identity and
the complete channel to be identity. For a normalized channel in nonzero
dimension, unit Born return probability on every coordinate state and the
normalized uniform state implies those fixation conditions and hence the
committed Lüders map of `I`. Its eleven theorem audits use only `propext`,
`Classical.choice` and `Quot.sound`. This formalization does not include
the multi-label reduction or the quantitative stability estimate.

The examples are exact algebraic models, not source-authenticated calibration
records. The protocol assumes the finite quantum instrument class and trusted
preparations and rank-one readback tests. It does not derive Born probabilities
from scalar consensus, choose a physical state space, implement its own trusted
tests or prove that the native repair instrument is Lüders. It turns an
all-state invariance premise into a finite conditional identification and
error contract. The parent certain-state theorem and the stated single-channel
finite-probe theorem are Lean checked; the complete multi-label theorem and
analytic stability bound are not presented as Lean-kernel theorems.
