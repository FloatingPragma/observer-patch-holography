# Continuing native records and their modular protocol algebra

This packet joins the existing native first-current channel, its stationary
Palm reference, the finite modular/KMS identity and two continuing observer
experiments. The existing active-corner phase target implements the
**full reference's modular state action on the
entire algebra used by these protocols**, even though it does not implement
that action on all of `M66`. No extra phase gates are needed for this equality.
The exact target equality and the nonzero compiled-source error are separate.

The positive continuing-record prediction is

```
P(degree at event 1, degree at event 2 | Palm herald succeeds)
             degree 8    degree 10
degree 8        1/5          1/5
degree 10       1/5          2/5.
```

Both degree-eight marginals are `2/5`; their connected correlation is `1/25`.
The observer continues its actual postmeasurement state between events.
Resetting to the reference would instead give `[[4/25,6/25],[6/25,9/25]]`.
These are finite model predictions under the declared scheduler and read
apparatus, not new laboratory observations or a recovery of a continuum QFT.

## Source channel, reference and parameter

The source is the existing binary, total-two native sector with 66
configurations, uniform selection of thirty seams and a fair placement coin.
Both existing record laws C and E are checked. A native attempt produces
either zero or nonzero absolute current. Its zero-current CP map `T0` is the
Schur multiplier

```
T0(rho)[i,j] = W[i,j] rho[i,j],  W[i,j] = integer_count[i,j]/60.
```

Every count is between zero and 52. The actual moving-current map is `T1`.
The first-event channel is the existing stopped construction

```
S = sum(n >= 0) T1 T0^n.
```

For each moving partial injection, its contribution to the matrix unit
`|i><j|` is exactly `1/(60-integer_count[i,j])`. We enumerate all 4,356 input
matrix units using rational arithmetic. The CP limit is trace preserving:
the unresolved trace after `N` attempts is at most `(13/15)^N`, also for an
arbitrary untouched reference. This follows from the uniform escape bound
on `T0* (I)`; it is not an inference from a selected initial state. The
infinite limit is analytic. Every executed finite experiment retains its
timeout outcome.

The source-derived degree census is `{8:30, 10:36}`, with sum 600: an
occupied pair has 8 moving branches when its ports are adjacent and 10
otherwise. The frequency `log(5/4)=log(10/8)` is therefore a combinatorial
identity of the carrier. Hence

```
sigma = diag(degree/600),
S(sigma) = sigma,
-log sigma = log(60) I + log(5/4) P_degree8.
```

This is event stationarity. One ordinary native attempt changes this
reference by total variation `1/150`. A movement herald on one specified
attempt starting from `I66/66` prepares `sigma` with probability `5/33`;
failure `28/33` remains in the record. Waiting for the first movement from
`I66/66` gives a different distribution. Thus preparation, first-event
sampling, attempt stationarity and event stationarity are kept distinct.

For observables the convention is
`alpha_s(A)=sigma^(is) A sigma^(-is)`. The dual state action implemented by
`modular_action` is `rho -> sigma^(-is) rho sigma^(is)`. Its exact quarter-turn
target has `theta=s log(5/4)=pi/2`. The finite KMS boundary identity is
`omega(A alpha_(-i)(B))=omega(B A)` with `omega(A)=Tr(sigma A)`; its potentially
nonzero matrix-unit cases are checked over all of `M66`. The opposite sign
is required when parameterizing the observable flow by `H=-log sigma`.

The source event counter and the modular parameter `s` are separately named
observer parameters. A modular parameter may be proposed as an observer's
time without introducing a universal clock. This construction checks how its
programmed interventions relate to continuing native records; it does not
identify their two dynamics.

## An observer-relative modular phase readout

The successful branch of the existing native coherent preparation has
nonzero coherence between the two eigenspaces of `sigma`; its modular
frequency support is exactly `0` and `+/-log(5/4)`. The reference `sigma`
itself remains stationary under its own modular flow. The coherent branch is
a distinct state and therefore supplies a possible relational clock hand.

Apply the ideal target modular action for parameter `s`, pass that state
through the existing first-moving-current stopped channel, then use the
existing Y readout and condition on the original `5/33` herald. The exact
conditional plus probabilities are

```text
law C:  P(plus | s, herald) = 1/2 - sin(s log(5/4))/140
law E:  P(plus | s, herald) = 1/2 - sin(s log(5/4))/220.
```

The coefficients follow from exact rational evaluation at modular phases
`0`, `+/-pi/2` and `pi`; the two-eigenspace source has no higher modular
frequencies, so linearity fixes the whole one-harmonic response. The quarter
turn points are checked in the producer, including the executed
positive quarter-turn phase-order experiment. `s` has period
`2*pi/log(5/4)`. The nonzero derivative at `s=0` makes this an observer
readout of modular phase locally near that point. One probability is periodic
and globally ambiguous, so it is not a unique absolute timestamp.

This is a concrete bridge from the native observer's coherent preparation and
readout to a modular parameter. The continuous curve is an exact prediction
for the ideal modular target, not a continuous native gate execution: the
repair process has not been shown to generate or advance `s`. In fact, the
current native rule cannot generate this clock hand from its existing
successful coherent preparation. That prepared density is exactly real, every
native Kraus branch is a real partial injection, and the Y readout detects the
imaginary coherence generated by modular phase. The exact attempt channel
therefore maps real densities to real densities: its conditional Y-plus
probability stays `1/2` after every integer number `n` of native repair
attempts. The executable boundary test checks `n=0,...,8` against the actual
primitive runtime for both laws C and E; real-linearity proves the same result
for every `n`. At a modular quarter turn, the ideal target instead gives
`69/140` for C and `109/220` for E. Thus the ideal modular intervention is an input, and replacing it with the current native repair dynamics would erase
the phase signal. This is an obstruction for the committed finite observer
implementation, not a theorem that no other OPH construction could generate
modular evolution.

The real-state obstruction also survives conditioning on any retained native
record, record-dependent finite stopping and probability mixtures of those
histories, provided every intervening operation is real in the native basis.
Each branch maps a real density to a real unnormalized density; normalization
by a positive real probability preserves that property. A convergent stopped
CP sum does too. Thus selecting a different native event boundary or a
nonzero-probability record cannot turn this same Y protocol into a moving
modular hand. An alternative detector can read native change: the
[single-carrier clock comparison](../native_clock/COMPARISON.md) gives actual
current-count agreement and uncertainty without a federation. Those clocks
do not require modular phase advancement.

The source does supply a local carrier interpretation for this frequency.
Each of the four cross-degree coherent modes in the existing phase preparation
is supported on three named ports. Restricting the Palm reference to the
binary states with both loads inside that triple gives three positive weights
`(5/13, 4/13, 4/13)` in an order fixed by the mode. This is a faithful
three-state conditional reference, and the coherence joining weights `8` and
`10` has the same source-derived modular frequency `log(5/4)`. The condition
has Palm mass `13/300`. The executable audit derives this from the carrier and
degree census for all four modes.

This restriction does not produce three autonomous port qubits or a closed
three-port evolution. The current carrier is the 66-state binary portion of a
fixed-total-two load space (78 states when double occupancy is included).
The 12 port labels are coordinates of those load states; the current model
does not use a `2^12` tensor-product Hilbert space. The existing event-Palm
state is the available faithful nontracial reference. By contrast, the
forgotten-record native attempt channel has unique fixed operator `I66/66`,
whose modular flow is trivial. Its other existing complex preparation uses an
untouched reference qubit; tracing that qubit leaves a nonfaithful diagonal
work state. Neither alternative supplies continuing modular phase motion.

## Passive repair cannot act as a modular automorphism

The obstruction extends beyond the real C/E preparations. The
[law-uniform warm-up theorem](../record_laws/README.md#preparation-by-finite-native-warm-up)
gives `tau=I66/66` and, after `n` unread attempts,
`D(Psi_n(rho),tau)<=epsilon_n`, where `D` is half trace distance and
`epsilon_n=min(1,33 r_alpha^n)`,
`r_alpha=1-alpha*(1-19079/20000)`. The domain is every admitted current table,
including its phases and the theorem's restricted law-only feedback, under
exogenous seam coverage at least `alpha/30`. Discard the warm-up record and
controller memory; do not insert other source operations.

For every faithful matrix reference `sigma`, every modular parameter `s`,
and the dual modular action `U_s rho U_s*`, unitary invariance gives
`D(U_s rho U_s*,tau)=D(rho,tau)`. The reverse triangle inequality therefore
proves

```text
D(Psi_n(rho), U_s rho U_s*) >= max(0, D(rho,tau)-epsilon_n).
```

The same statement holds with an arbitrary untouched reference `R`, replacing
`tau` by `tau tensor rho_R` and acting with `U_s` only on the source. Any
deterministic reparameterization `s=s(n)` obeys it. A record-conditioned
mixture of phase actions, postselection, and active coherent controllers are
different comparisons. At 221 uniform attempts, `epsilon_n<=1/1000`:
every pure source input is at least `32467/33000` from its modular orbit.
For the native Palm reference, `D(sigma,tau)=3/55`, so even its stationary
modular orbit has separation at least `589/11000`. The tracial reference's
modular action is identity; passive repair is not identity on general inputs.

Restricting to an observer algebra does not evade the uniform dynamical
obstruction. Let `A` be a unital star subalgebra invariant under the chosen
modular group, and `E_A` its trace-preserving conditional expectation. The
warm-up proof gives `||Psi_n-Replacement||_(2->2)<=r_alpha^n`; the adjoint
has the same Hilbert--Schmidt norm. For any traceless `X` in `A`,

```text
||E_A Psi_n*(X)-alpha_s(X)||_2 >= (1-r_alpha^n)||X||_2.
```

Indeed `E_A` is an orthogonal contraction and the modular action is an
isometry. Every nonscalar observer algebra contains such a nonzero `X`.
This excludes agreement as maps on all admitted states, while allowing
agreement on a single stationary state or scalar readout. The scalar algebra
passes but has no moving hand.

[`passive.py`](passive.py) independently enumerates actual C/E branches and
provides an exact finite checker for a rational-complex operator basis and a
faithful diagonal reference. It checks the ambient identity, adjoints,
products, modular invariance, leakage and norm loss. Modular invariance is
checked without logarithms: group matrix units by equal `sigma_i/sigma_j`
and test that each spectral component lies in the algebra. The native degree
centralizer and a coherent `M2` block both exhibit strict norm loss. Passing
this necessary isometry test alone is not a modular-flow certificate; the
generic callback must separately be shown to be an admitted channel.

There is also a direct stopped-event witness, without replacing an event by
its mean waiting time. Start from `tau`, make one native movement herald,
read degree, wait for one native movement with cap `h`, and read degree again.
Keep the failed herald, every successful waiting count, and both degree
timeouts. A matched control uses each initial degree's same waiting law and
timeout mass but preserves that degree under an arbitrary modular phase.
The degree projectors commute with `sigma`, so every such phase has the same
control records. Exact C/E enumeration gives the unconditional record distance

```text
TV = [2-(13/15)^h-(5/6)^h]/33.
```

At cap two this is `499/29700`; the limit is `2/33`, or `2/5` conditional
on the `5/33` herald. Failed heralds and timeouts match exactly. This control
is a declared comparison law, not an implementation claim. The witness uses
the existing degree-read service and a selected preparation; the general
reference-safe bound above concerns attempts, not arbitrary stopping rules.

The [thermal-time proposal of Connes and Rovelli](https://arxiv.org/abs/gr-qc/9406019)
motivates a state-relative modular automorphism as time. These finite matrix
arguments need no imported infinite-algebra theorem. They test whether the
declared passive native channel realizes that automorphism, independently of
choosing units for its parameter.

## The modular readout has no lab-time scale

The modular response fixes a frequency `log(5/4)` per unit of the
dimensionless parameter `s`; it does not fix a frequency per second. The
executable `modular_lab_time_scale_identifiability` audit substitutes
illustrative maps `s = kappa*t` with `kappa` values `1`, `2` and `5/3` per
second. They preserve the exact probability at the corresponding modular
quarter-turn, while their lab periods differ by factors `1`, `1/2` and `3/5`.
No source quantity selects a value of `kappa`.

The native attempt chain does not generate `s`: its real preparation and real
repair maps keep the Y outcome at `1/2` for every integer attempt count. Thus
calibrating a lab clock alone would not implement the modular flow. A test
against lab data first needs a physical operation that realizes this phase
action, an identified detector/preparation mapping, and an independently
measured nonzero pair `(delta-s, delta-t)`; a constant `ds/dt` must be justified
over the tested range. No matching public dataset is identified for
this protocol, so the response curve is not an empirical lab prediction.

The related [`native_clock`](../native_clock/README.md) packet independently
shows that native attempt/event counts also have no seconds scale: arbitrary
positive seconds-per-attempt values preserve the same exact record law.
Continuous-gate and detector implementation errors are not bounded here.

## The invariant algebra closes the active-corner gap

The ten fixed-seam orientation pairs define the unital protocol algebra

```
A = (direct sum of ten M2 blocks) direct sum C46,
complex dimension(A) = 86.
```

Its matrix-unit basis consists of every native diagonal and the twenty
oriented pair coherences. Exhaustive support checks show both `S(A) subset A`
and `S*(A) subset A`. Every individual native moving Kraus branch, the
fixed-seam herald branches, their adjoints and the declared preparation,
phase, degree and phase-read instruments preserve this algebra. Zero-current
Kraus maps are diagonal. The source no-current map commutes with the full
modular action on every matrix unit.

The full reference's modular action preserves `A`. The compiled
target phases eight degree-eight coordinates in the twenty-dimensional
active corner. It omits 22 other degree-eight coordinates. Those omissions
only change scalar phases on one-dimensional blocks of `A`, which cancel in
conjugation. Consequently the **ideal active phase target equals the full
modular state action on all of `A`**, throughout these continuing protocols.
This is stronger than restricting the whole state to the active corner:
native events can populate the other 46 diagonal blocks. On all of `M66` the
equality remains false; the explicit control `|0><1|` distinguishes the two
actions. Coherences outside `A` are not silently projected away.

The reused complete-source certificates compare every coordinate of
the compiled phase and mixer, including the 512-dimensional spectator and
padding extension. Restricting the exact equality to `A` does not permit
discarding source leakage outside it.

## What executes and what is measured

The two-event experiment starts with an actual primitive-runtime
`Program`: `I66/66`, one native attempt, and its current herald. On success,
an explicit finite stopping driver repeatedly calls the existing runtime
instrument executor, retaining the unnormalized zero-current density until a
movement occurs. It applies the declared degree projector, continues that
postmeasurement matrix through a second stopping block, then reads degree
again. The complete loop is an explicit driver of the same CP instruments;
it is not represented as one shared-runtime `Program` or a serialized
checkpoint. At cap two per event, all 22 public records are retained:
failed preparation, first-event timeout, first degree with second-event
timeout, and both completed degree readings with their attempt counts.
Hidden source Kraus alternatives are summed by the declared instrument.

The phase-order experiment starts from the existing fixed-seam coherent
preparation program: failure `28/33`, plus herald `5/66`, minus herald
`5/66`, with the balanced mixer only on success. It applies a quarter modular
phase either before or after the first native event. A declared Y detector
on the equal-degree pair `(64,54)` reads the result; the herald sign controls
bit decoding. The detector gives an unbiased bit on the other blocks. At
cap four, each order retains all 19 public records, including both herald
signs, every event count/read bit and both unresolved tails. Failed branches
incur no mixer, phase or detector cost; phase-after-event timeouts incur no
phase cost.

Every terminal 66-by-66 matrix is compared with an independently assembled
exact rational-complex oracle. The degree and Y detectors and phase gate in
these continuing runs are ideal native-sector target services. The native
attempt Kraus arrays come from the existing source enumeration. Floating
matrix propagation is numerical evidence; neither it nor detector hardware
errors receive an exact implementation certificate here. No terminal state
is conditionally normalized, reset, projected to a target state or dropped.

The probability measured by the phase probe distinguishes intervention
order. The exact infinite-event conditional plus probabilities are:

| Native law | Phase before event | Phase after event | Original-herald joint probability difference |
|---|---:|---:|---:|
| C | `69/140` | `1/2` | `1/924` |
| E | `109/220` | `1/2` | `1/1452` |

The last column includes the original `5/33` total herald yield. The
finite executions predict their own capped probabilities, rather than
relabeling the infinite values as finite observations. Analytically, the
full stopped channel mixes modular frequencies in 3,360 nonzero matrix
coefficients for each law. This operational obstruction exists even on `A`:
stationarity and KMS do **not** make the irreversible event channel the
modular automorphism group, or even a channel commuting with that group.

For two successive event blocks, the conditional unresolved mass is bounded
by `min(1,2(13/15)^N)`; multiply by `5/33` for its contribution before
conditioning on preparation. Exact finite probabilities and unresolved
masses are retained at caps 2, 8, 32, 64 and 128. For example, cap 64 has
conditional unresolved mass about `9.45e-5`. The short cap-two execution is
a complete execution cross-check, not a tightly converged infinite-event
experiment.

## Existing theorem attachment and remaining scope

The stopped CP construction is reused from
[`native_relational_readout`](../native_relational_readout/README.md), with
the Palm law and native event semantics from
[`event_circulation`](../event_circulation/README.md).
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean) supplies the
finite inner-flow setting;
[`StructuralInheritance.lean`](../../Lean/QFT/StructuralInheritance.lean)
supplies the finite Gibbs-KMS identity. Here its Gibbs matrix is explicitly
`H=-log sigma`, obtained from the native Palm reference. The existing Lean
Heisenberg parameter has the opposite sign to `alpha_s` above; the exact
boundary identity is the same after that substitution. This packet checks
the concrete matrix mapping and adds the invariant protocol algebra,
continuing measurements and operational covariance witness. It claims no
new abstract Lean theorem.

[`TimeOrderLedger.lean`](../../Lean/Time/TimeOrderLedger.lean) keeps causal
order, modular parameter and clock identification distinct. The reference
used here is a native density on `M66`. The transfer-Hilbert-space reference
`K/Tr(K)` in the effective-field reconstruction is a different construction;
equal carrier dimensions would not identify their operational algebras or
preparations. Neither reference is asserted to be selected uniquely by OPH.

Fresh replay recomputes the complete source phase and mixer comparisons
using [`native_modular_reference`](../native_modular_reference/README.md)
and [`compiled_native_interface`](../compiled_native_interface/README.md).
The fixed operator-error gates are `1e-10` and `1e-11`, respectively. These
are separate source-gate bounds, **not** an end-to-end certified error bar
for this packet's native-sector prediction. A complete source implementation
of the degree/Y instruments, native evolution on the full encoded
carrier, and their composite error accounting are not supplied.

The next physical-time attachment must explain the joint continuing record
and intervention statistics, and connect them to the accessible effective
fields. Simply renaming the native event counter as this modular flow would
contradict the measured-order prediction above. The present result provides
an exact closed algebra for making that comparison and continuing observer
correlations without a fitted action or an extra fundamental clock.

## Replay

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.relative_modular_dynamics --check model/relative_modular_dynamics/results.json
python -m model.relative_modular_dynamics.passive --check model/relative_modular_dynamics/passive_results.json
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/relative_modular_dynamics
```

The retained receipt is software evidence. Exact rational predictions,
record counts, error envelopes and source bindings replay strictly; only
explicit floating diagnostics use the existing numerical comparison rule.
