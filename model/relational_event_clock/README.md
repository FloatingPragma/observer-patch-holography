# Native records on a local event clock

This packet predicts complete continuing record laws when a clock ticks on an
actual occupation change at one native port. Remote native repairs continue
between ticks. The exact stopped quantum instrument removes independently
inserted identity holds from event-indexed predictions, retains the effect of
relative seam weights, and distinguishes the two native quantum completions
with a phase-sensitive readout. These are finite model predictions; no new
experimental dataset or physical clock calibration is used.

## Source and observer access

The carrier is the existing 66-dimensional binary total-two sector on the twelve
source ports, with the original thirty seams and native fair coin. A source
attempt chooses a seam uniformly. The two candidate completions C and E retain
their original committed Kraus labels; hidden labels are probability-summed,
not coherently recombined. Every native branch is a partial injection, and the
sum of its exact effects is the identity, including with an arbitrary spectator
reference. Floating Kraus execution is compared with exact rational targets;
its residual is a numerical diagnostic, not a hardware or arithmetic certificate.

The clock port is 2, the first source port outside the existing phase pair
(0,1). Its neighbours are 3, 4, 6, 10 and 11. The observer reads the source seam
current and records `(neighbour, change at port 2)` when that change is nonzero.
The other outcome is `no tick`. No direct occupation measurement or hidden-state
read is added. This detector grouping and the seam schedule are declared
observation premises. Tick number gives record order; attempt count measures
computational resources, not elapsed physical time.

The [local-tick packet](../local_tick/SPEC.md) ticks at every incident
attempt, including no change. Here a tick requires actual local movement, and
complete marked two-event laws, an event-triggered intervention and coherent
terminal readout replace an occupation-moment-only calculation.

## Exact stopping and its boundary

Let Q be the no-tick CP map and R_m the CP map for a signed local mark m, using
output-row/input-column conventions. The completed event instrument is

\[
S_m=\sum_{n\geq0}R_mQ^n=R_m(I-Q)^{-1}.
\]

Q contains genuine remote changes; it is not a diagonal holding probability.
Fifteen valid states have zero immediate local hazard, yet all eventually tick.
Dropping remote motion or dividing by the immediate hazard therefore gives the
wrong event map. The executable computes exact rational blocks on native
signed-difference sectors, rather than forming a dense 4356-dimensional solve.

For every initial native state, including correlations with a reference, the
probability of no tick after two attempts is at most 1799/1800. Products of
partial injections make the survival effect diagonal, so the checked classical
maximum is also the quantum/reference-uniform bound. Iteration gives
`Pr(no tick after 2m attempts) <= (1799/1800)^m`. Thus the CP series converges to
a complete event instrument. The report also gives sharper exact bounds at
specified finite caps. No infinite execution is claimed: actual runs keep the
unresolved native density as a distinct timeout outcome.

Insert an independent identity hold with activity a, where `0 < a <= 1`:
`Q_a=(1-a)Id+aQ`, `R_(m,a)=aR_m`. Then `S_(m,a)=S_m` exactly. This preserves any
complete event-triggered protocol built from these maps. It does not preserve
attempt timestamps, work costs or finite-cap success, and does not cover other
work evolving during the holds. Doubling the relative weight of remote seam
(0,1) changes the two-event law, providing a scheduler control.

## Continuing histories and predictions

The first experiment starts from the uniform native state and performs one
specified native attempt. A local movement herald succeeds with probability
5/198 and fails with probability 193/198. The successful normalized state is the
local Palm reference: the distribution seen at an event sampled from the
stationary attempt process. Its support has 51 addresses and its mean intertick
cost is 198/5 attempts. It is stationary for the completed event map. It differs
from stopping at the first tick after uniform initialization by total variation
1/3; those two preparation protocols must not be substituted.

Two continued ticks give all 100 signed-mark pairs and their full native
terminal states. The clock signs alternate. Conditional on the initial herald,
the probability that both ticks use the same neighbour is
`5428819044974590805373 / 19780089058627344141475`, approximately 0.274459,
versus 1/5 for independent neighbours with the required sign alternation.
After the first tick, an optional native attempt on remote seam (0,10) retains
all three current outcomes and gives 300 joint records. This intervention changes
the pooled two-event distribution by total variation approximately 0.00617194.
It changes the continuing state; neither protocol resets between ticks.

A second experiment uses the existing
[clocked quantum preparation](../clocked_quantum_readout/README.md): a native
(0,1) herald followed by the balanced phase mixer. Its total success is 5/33,
with two herald signs of weight 5/66 each; failure 28/33 remains in the complete
law. After one or two local ticks, the existing balanced mixer and native
pointer read phase match relative to the retained herald sign.

| Native completion | Match after one tick, given herald | Match after two ticks, given herald |
| --- | --- | --- |
| C | 3804737/7050540 = 0.539638 | 3588793/7050540 = 0.509010 |
| E | 289535633/552465980 = 0.524079 | 278113337/552465980 = 0.503404 |

C and E have identical diagonal clock records, but different phase probabilities
on the same local event clock. The one-tick separation is approximately
0.0155590 conditional on successful preparation, or 0.00235743 in the pooled
unconditional match probability after charging preparation success. These are
predictions of the specified native completions and detector protocol; neither
candidate is selected by an empirical observation here.

The exact stationary inputs above are targets. The
[native equilibration packet](../native_equilibration/README.md) supplies a
separate finite preparation attachment: 374 native attempts give trace-distance
error at most 33/131072 from the tracial native target, uniformly with a retained
reference. Any subsequent complete instrument inherits this unconditional
error by contraction. Conditioning on a rare herald requires an additional
postselection bound; the exact conditional numbers are not asserted unchanged
for that finite preparation.

## What runs

`python -m model.relational_event_clock` runs the actual native Kraus instruments
on continuing, unnormalized 66-by-66 matrices. The local herald is executed with
the primitive runtime and its serialized checkpoint is restored. Two stopped
blocks with cap two retain all 422 baseline records or 1262 intervention records,
including preparation failure, both timeout positions and zero-probability
outcomes. Every terminal matrix is compared with its exact target; exact full
terminal states are hashed in the report.

Each C/E coherent run uses cap three and retains 123 complete records: failed
preparation, signed local events, both phase results and each unresolved tail.
The runtime applies the source Kraus maps; resolvent answers are not substituted
for the evolving state. All native attempts and balanced mixer calls are charged
per branch. The eventual one- and two-tick phase probabilities are exact CP-series
predictions, separate from these bounded executions.

The coherent runtime uses the existing ideal 66-dimensional balanced-service
target. [Compiled native interface](../compiled_native_interface/README.md)
supplies the separate full-512-dimensional source gate comparison. This packet
does not claim a new error certificate for the complete hardware protocol or
silently discard implementation leakage.

## A local modular reference, with a strict distinction

The local Palm state is faithful on its corner `P M66 P`, with 51-dimensional
support (matrix-algebra vector-space dimension 2601). Its weight census is 30 entries of 1/100, ten of 1/50, five of 1/25 and six
of 1/20. This supplies the finite modular automorphism
`alpha_s(X)=sigma^(is) X sigma^(-is)` on that corner. Zeros outside the support
are not an obstruction to this supported construction. A completed local event
maps the entire native carrier into the corner; a no-tick branch can leave it,
so unresolved states must remain on the full carrier.

The event map is not the modular automorphism. An input matrix unit between
occupied configurations (0,2) and (2,3) has modular weight ratio 5/4. Its image
contains the matrix unit between (0,4) and (3,4), of ratio 1/2, with coefficient
1/14 for C and 1/24 for E. A single marked source attempt contributes 1/60.
The unequal frequencies and nonzero coefficients prove noncovariance with the
corner modular flow. Here S acts on states; the corresponding inverse-sign
state flow has the same commutation obstruction. No modular gate or
identification with physical duration is executed by this witness.

## Existing theorem connections and remaining obligations

[Regional continuity](../../Lean/Screen/RegionalContinuity.lean), specifically
`regional_continuity`, supplies the exact singleton-region relation between
load change and signed seam current used to define the tick. The implementation
checks that relation against every original native transition. Its finite
update index carries no physical duration.

[TimeOrderLedger](../../Lean/Time/TimeOrderLedger.lean) distinguishes repair
order, observer records, modular parameter and calibrated time.
[ClockReadout](../../Lean/Time/ClockReadout.lean) provides clock regrading and
affine comparison algebra with explicit correspondence/calibration premises.
The local event count is an ordered record readout; none of those conditional
comparison theorems supplies missing seconds or spacetime geometry. The
[relative modular dynamics](../relative_modular_dynamics/README.md) packet's
native preparations and signed-difference sectors are reused, with a distinct
local-event Palm reference and a directly checked modular obstruction.

The advance is a complete source event instrument and intervention law with an
exact irrelevant-hold invariance, a phase discriminator and an attached faithful
local modular corner. Relating this clock's records to observer modular flow,
calibrated physical comparisons or an effective continuum spacetime requires
additional structure and evidence. No larger primitive, new field action,
preferred universal clock or continuum QFT conclusion is introduced.

Run the focused checks with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/relational_event_clock`.
The CLI accepts `--write PATH` for a fresh receipt and `--check PATH` for a fresh
canonical replay against that receipt.
