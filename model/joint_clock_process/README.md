# Joint local clocks and quantum fields with explicit memory

This packet supplies an exact finite correspondence between native repair
histories and a declared observer's effective field, two local clocks and
interventions. It preserves complete outcome probabilities and conditional
quantum outputs, including an untouched reference. It extends the
[first-moment memory identity](../federated_observer/MEMORY.md) to quantum
instruments. It also finds a simpler, exact description at this experiment's
completed clock events: a four-mark classical predictor, together with the
retained conditional reference operators. An instrument-valued change from
the union clock to the A clock derives their relative count law,
including feedback, remote transitions and full reference outputs. No assumed
deterministic rate conversion is needed for this finite comparison.

The construction is a finite operational reduction. It does not establish a
continuum field theory, a relativistic spacetime, or an identification of native
event order with modular flow. The field below is an observer-chosen accessible
algebra; its coordinates are not fundamental fields or a derived spatial
dimension.

## One source, two local clocks

The existing [federated quantum packet](../federated_quantum_sector/README.md)
supplies the source path with active ports `(0,1,72,82)` and seams
`(0,1), (1,82), (72,82)`. Ports 0 and 1 belong to one elementary twelve-port
carrier; ports 72 and 82 belong to another. The controlled experiment chooses
these three seams uniformly and uses the original fair native placement coin.
It retains one particle on the four active ports, with all other 236 source
ports empty. Every selected operation preserves that support, including all
superpositions and reference correlations. This schedule is a declared service;
it is not the uniform schedule on the entire federation.

Clocks A and B read the signed native occupation changes at ports 0 and 72.
Their joint mark is one of

```text
no event, A exit, A entry, B exit, B entry.
```

Only actual nonzero repair currents tick a clock. The seam and coin remain
hidden CP alternatives. The no-event map Q includes genuine remote motion
between ports 1 and 82. The declared joint record interface makes both clocks'
records available for comparison; message delivery times and physical
simultaneity are not inferred. Native attempt number is a resource counter.
The two tick counts are observer-local records. Neither parameter is asserted
to be a modular parameter or calibrated seconds.

Let R_m be a nonzero-mark CP map. Direct source evaluation on all sixteen
matrix units gives `Q*(I)=(5/6)I`. Thus the chance of no union event within k
attempts is exactly `(5/6)^k` for every initial state and untouched reference.
The eventual instrument is

```text
S_m = R_m (I-Q)^(-1) = sum_(n>=0) R_m Q^n.
```

The uniform survival bound proves convergence and completeness. The expected
cost of the next union event is six native attempts. Dividing R_m by the
one-step event probability is wrong: Q contains remote motion and is not
`(5/6) Id`. Independently inserted identity holds cancel from the completed
instrument, by the same resolvent identity as the
[single local-clock construction](../relational_event_clock/README.md).
They change finite-cap probabilities and resource costs.

## An exact event-boundary simplification

In the ordered active basis, the four completed maps have the exact form

```text
S_m(rho_WR) = |d_m><d_m| tensor Tr_W[(E_m tensor I_R) rho_WR],

d = (1, 0, 82, 72),
E_Aexit  = diag(1, 0,   0, 0),
E_Aentry = diag(0, 2/3, 0, 1/3),
E_Bexit  = diag(0, 0,   1, 0),
E_Bentry = diag(0, 1/3, 0, 2/3).
```

The addresses in d are native port addresses, not matrix indices. Each final
nonzero-current Kraus branch has rank one on this one-particle support. The
source calculation checks the full measure-and-prepare identity, not only
traces. It therefore holds with every untouched reference.

The next-mark transition matrix, with columns denoting the preceding mark in
the order above, is

```text
       A exit  A entry  B exit  B entry
A exit     0       1       0       0
A entry  2/3       0     1/3       0
B exit     0       0       0       1
B entry  1/3       0     2/3       0
```

The stationary mark probabilities are each 1/4. Consecutive union events use
different clocks with probability 1/6, rather than the 1/2 from independent
uniform clock identities. Starting after any event, the probability of no A
tick during the next 2k union events is at most `(2/3)^k`; the same bound holds
for B. Thus each local clock continues almost surely. These counts and costs
do not supply physical clock rates. This chain is periodic; a stationary law
does not mean convergence of every initial mark distribution.

Four marks suffice for future work-record probabilities under fixed or
record-controlled work-only instruments. Such an intervention updates the next
transition using its action on the known state `|d_m><d_m|`. A controller that
consults more history must retain that classical history. Complete reference
outputs additionally require the conditional reference operator in the formula;
four scalar probabilities cannot replace the reference assemblage. Operations
on a correlated reference or additional hidden controller memory are outside
this four-mark closure.

C and E have exactly the same completed event instrument here. Each completed
event destroys work coherence in this chosen one-particle experiment. This is
a property of this instrument and sector, not a claim about all observer
clocks. Coherence survives in unresolved waits and can be read before an event;
discarding those waits would remove the discriminator below.

## A complete change between observer clocks

The [clock conversion](clock_change.py) predicts how many B ticks occur before
the next A tick, including the conditional field/reference state. Let `Q` be
the no-union-event map, `A_a` the two native A-event maps, and `B_b` the native
B-event maps followed by any declared B-triggered intervention. Define
`B=sum_b B_b` and the completed union maps

```text
S_Aa = A_a (I-Q)^(-1),       S_B = B (I-Q)^(-1).
T_Aa(z) = A_a (I-Q-z B)^(-1)
        = S_Aa (I-z S_B)^(-1).
```

The coefficient of `z^k` is the CP map for the next signed A tick after
exactly `k` B ticks. It sums genuine remote transitions, rather than treating
them as identity holds. The two expressions are equal by factorizing
`I-Q-zB=(I-zB(I-Q)^(-1))(I-Q)`. The executable checks the symbolic rational
matrix identity on the full sixteen-dimensional operator space for both
native laws and both feedback settings. Thus direct elimination of native
attempts and staged elimination through the union clock preserve the same
complete count instrument, not just the expected rates. Retaining individual
B signs uses the ordered products `S_Aa S_Bb_k ... S_Bb_1`; their sum gives
the coefficient above. Reference systems are untouched identity extensions.
History-dependent controllers require their classical memory to be included
in this map; the demonstrated feedback depends only on the current B sign.

The exact effect inequality `(S_B^*)^2(I) <= (2/3) I` holds both without
feedback and with a field Hadamard after each B entry. Consequently the
probability that A has not ticked after `2k` union events is at most
`(2/3)^k`, uniformly over work/reference preparations. This proves completion
almost surely and gives a complete finite censoring bound. An independent
identity hold with any fixed activity `0<a<=1` cancels from every completed
union map and hence every relative-count coefficient. Identity holds do
change native-attempt caps and resource costs. Changing the relative seam
weights, allowing evolution during an alleged hold, or timing feedback by
attempt count is outside this invariance.

An actual A exit leaves the known work state at native port 1; an A entry
leaves it at port 0. No reset or independently selected hidden state is needed
to predict from these observed boundaries. Without feedback the B-count
probability-generating functions are

```text
after A exit:   g(z) = (2-z^2)/(3-2z^2),
after A entry:  g(z) = 1.

P(B count=0 | A exit) = 2/3,
P(B count=2k | A exit) = (1/9)(2/3)^(k-1),  k>=1.
```

The next A sign alternates in the passive experiment. Starting just after an
A exit, each two-A-tick cycle returns to that same known work state. Under
the declared independent Markov scheduler these cycles are independent for
work-only records, even when an untouched reference is retained. Each cycle
has mean B count **2** and variance **16**. After `N` cycles the mean ratio
of B ticks to A ticks is **1**, its variance is **`4/N`**, and the probability
of having seen no B ticks is **`(2/3)^N`**. Equal average rates do not
give an exact affine relation between the two finite clock records. The
expected native resource costs of the two passive local intervals are 18
and 6 attempts; these are not calibrated durations.

A field Hadamard after a recorded B entry changes the first generating
function to `(4+z-z^2)/(2(3-z^2))`, with mean B count **`3/4`**, variance
**`35/16`**, and probability **`1/6`** of exactly one B tick. The feedback can
change occupation between native ticks; its action is explicitly charged and
is not itself counted as a native-current tick. The two signed A maps retain
that changed sign law. C and E have the same *completed* count
instrument under this event-triggered feedback. This is a controlled
prediction about a changed experiment, not a regrading of the same clock.

Finite executions run the actual primitive source instruments until A ticks
or three native attempts have occurred, then read the field-X detector on
the continuing state. All **106 physical records**, including timeout, are
retained. The fixed runtime menu has 250 branches: after completion it repeats
a cost-free administrative label, executes no further repair and adds no
clock tick. Impossible menu branches are zero and are removed only by an
explicit record coarse-graining that sums their full states. The source
runtime, independent rational stopped-tree oracle, source-bound checkpoint
and all terminal matrices agree. Two-attempt entangled-reference runs retain
34 physical records per law and include the B-triggered intervention.

| Capped protocol | A timeout probability | Unconditional field-plus probability, C | Same probability, E |
| --- | --- | --- | --- |
| Passive, three native attempts | `85/108` | `35/54` | `14/27` |
| Hadamard after B entry, three native attempts | `679/864` | `547/864` | `437/864` |

These field-plus probabilities include completed and unresolved A intervals.
Each completed A branch destroys work coherence in this support; the
unresolved branches retain the C/E distinction. Replacing the complete
capped record by the successful-clock subset would discard that information.
The finite execution starts with the disclosed coherent probe preparation;
the boundary-based count law applies after any actual A event. Neither
requires exact stationary preparation, a supplied clock-rate ratio or an
identification with modular flow. The sector, three-seam schedule, record
collector and coherent apparatus are the stated inputs. No laboratory
data are identified with this clock protocol.

## Full quantum field histories require memory

Identify the ordered four-dimensional support with `field qubit tensor hidden
qubit`: the field bit names the carrier and the hidden bit names its within-
carrier address. The effective field algebra is `M2 tensor I2`. Coherent field
operations are declared joint apparatus, not inferred spatially local gates.
Define the state restriction R by tracing the hidden factor and choose the
trace-preserving lift `J(x)=x tensor I2/2`. Then `RJ=I`. For a retained record
word, let L_j be its actual outcome map, including its declared intervention.
With `F=I-JR`, define

```text
A_j = R L_j J,    B_j = R L_j F,
C_j = F L_j J,    D_j = F L_j F,
x_0 = R rho_0,   y_0 = F rho_0.
```

Block substitution gives, for every finite word,

```text
x_(n+1) = A_(n+1) x_n + B_(n+1) D_n ... D_1 y_0
          + sum_(j=0)^(n-1) B_(n+1) D_n ... D_(j+2) C_(j+1) x_j.
```

An empty product is the identity. This follows by inserting
`y_(j+1)=C_(j+1)x_j+D_(j+1)y_j` repeatedly into
`x_(n+1)=A_(n+1)x_n+B_(n+1)y_n`. No Markov, Gaussian, equilibrium or unitarity
assumption is introduced. Outcome maps are unnormalized; their traces are the
full word probabilities. The identity includes arbitrary CP interventions
whose selection depends on retained records. Applying it to each reference
matrix block proves the same full quantum output identity with any untouched
reference. The finite executable guards the recurrence at eight maps; the
algebraic identity itself holds at every finite length.

The memory terms and hidden forcing are generally not positive maps. Their
complete sum is the restriction of the original CP history. This predictor
uses source operator blocks and the declared initial hidden component;
it does not prove an efficient finite memory window or reconstruct an unknown
preparation from field records alone.

Two checks make that boundary concrete. Preparations at native ports 0 and 1
have the same reduced field state. Their completed first-event laws are
`(1,0,0,0)` and `(0,2/3,0,1/3)`, with total variation one. Even if the initial
hidden factor is maximally mixed, removing the memory kernel changes the
word `(no event, A exit)` from exact probability 5/72 to 1/16. The full memory
recurrence restores every reduced matrix entry, not just these probabilities.

## What executes and what it predicts

The shared primitive runtime executes three actual native attempts on continuing
four-by-four matrices, then a complete field-X measurement. The chosen coherent
preparation is `(port 0 + port 72)/sqrt(2)`. A controlled variant applies a field
Hadamard after each recorded A exit. This supplied apparatus can change
occupation; the clock is the specified native-current detector and does not
silently count the apparatus as another native tick. Its visible trigger and
resource charge remain in the experiment. Identity feedback makes no change.

Every five-symbol attempt word and both detector outcomes are retained:
`5^3 times 2 =250` terminal records per candidate and feedback choice. Records
with fewer than one or two events explicitly retain their timeout status and
full native state. No state is reset at an event. A checkpoint after the first
attempt and its feedback is serialized, imported by prefix replay and
continued. Two-attempt entangled-reference runs retain 50 complete records
each. The exact memory recurrence is checked against every terminal reduced
field/reference matrix, and every full work/reference matrix is compared with
the actual Kraus runtime.

Without feedback, the field-plus probabilities after three attempts are
**35/54 for C** and **14/27 for E**. The common probability of no clock event
is 125/216. Within that no-event branch the corresponding conditional plus
probabilities are 189/250 and 133/250. The complete law keeps that branch and
all successful-event branches; these are not postselected success scores.
Both candidates have the same completed clock law, so the field measurement
adds information absent from the clocks alone.

These are exact laws of the declared model. The declared support, preparation,
schedule, coherent apparatus and detector need an experimental correspondence
before comparing them with laboratory data. No public dataset is identified
with this protocol, and no empirical confirmation is asserted. Floating runtime
residuals are numerical diagnostics; this packet supplies no compiled pulse
words or hardware-error certificate.

## Theorem attachments and remaining obligations

The [native signed-difference reduction](../observer/QUANTUM_MIXING.md) and
[federated quantum packet](../federated_quantum_sector/README.md) supply the
source Kraus maps, invariant support and arbitrary-reference extension. The
event resolvent and irrelevant-hold identity reuse the
[local-event construction](../relational_event_clock/README.md).
[Regional continuity](../../Lean/Screen/RegionalContinuity.lean) identifies each
native signed current with the recorded port change; it does not describe the
separately supplied coherent apparatus. The
[time-order ledger](../../Lean/Time/TimeOrderLedger.lean) and
[clock readout algebra](../../Lean/Time/ClockReadout.lean) preserve the
distinction between event order, modular parameter and calibrated comparison.

The extension of the existing block-elimination identity is from first moments
to record-conditioned quantum maps and their complete output states. Its proof
is the finite substitution above, not a Lean theorem or an invocation of a
Gaussian theorem on this non-Gaussian process. The four-mark event closure is
an additional source-specific simplification checked on the complete operator
space. The clock-conversion factorization reuses that stopped-instrument
construction and gives a compatible local-count comparison with explicit
fluctuations and feedback. It does not instantiate the affine-clock theorem
by replacing the stochastic comparison with its mean. These identities
discharge the finite joint-law correspondence for this declared experiment.
Compatible refinement, controlled continuum limits, physical
calibration and relativistic spacetime recovery are separate obligations.

```sh
python -m model.joint_clock_process --check model/joint_clock_process/results.json
python -m pytest -q -W error model/joint_clock_process
```

For the [complete local-history truncation bound](../native_action/README.md),
the observed region must contain all four active ports: the declared field
Hadamard and phase readout span both carriers. The locality theorem does not
treat those services as operations confined to one elementary observer.
