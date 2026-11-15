# Same-code current and sector-action ambiguity

This CP-04/#1140 packet tests the native repair records against the sharp
`q=1` code and the [screened code](../particle_screening/README.md)
of the particle lane. It
replays the unchanged C and E laws on the binary total-two shell. Its
[campaign](campaign.json), [source-bound receipt](energy_bridge_results.json),
[implementation](energy_bridge.py) and [independent checks](test_energy_bridge.py)
keep source activity, correction resources and energy interpretation separate.
No CP-04 physical-energy gate is accepted.

## A common current rate for both code supports

For a uniformly chosen native seam and coin, the complete population kernel
has probability `1/60` on each allowed directed movement edge. Its other
outcomes, including all zero-current flags, are holds. Every basis pair in
the **two-support source-executed code**

```text
|0_L> = (|{0,2}> + |{1,3}>)/sqrt(2)
|1_L> = (|{0,3}> + |{1,2}>)/sqrt(2)
```

has degree ten. The **screened code** has four supports per column,

```text
|0_S> = (|{0,2}> + |{1,3}> + |{0,4}> + |{1,6}>)/2
|1_S> = (|{0,3}> + |{1,2}> + |{0,6}> + |{1,4}>)/2,
```

and each of its eight basis pairs also has degree ten. Its exact all-branch
code check is owned by [the particle source](../particle_screening/screening.py).
The nonzero-current
effect is diagonal on these supports, so **every logical input in either
code** has rate `10/60 = 1/6` per attempted uniform-seam repair for C and E.
This supplies the same candidate's native activity for the screened
particle lane. This packet's recovery implementation below belongs to the
two-support code; in this packet screened-code recovery is an algebraic
existence result. The separate [screened recovery](../screened_recovery/README.md)
executes the screened encoder and all 208 positive source-recovery branches;
[screened transport](../screened_transport/README.md) joins that recovery to
the same uniform-seam schedule and an explicit core-label reset.
The screened code's transported core changes value on `1/15` of uniformly
selected seam/coin attempts and leaves the native seam set on `1/30`, for
either law and either logical input. These are native metadata probabilities,
not an energy cost. The screened packets supply the reset protocol and compiler
resource counts; a physical rest-energy interpretation is open.

The charge-matched stationary co-moving `q=1` reference has uniform source
marginal, so its current rate is `5/33`; the code excess is exactly `1/66`
per attempt. That stationary benchmark averages over all core labels,
whereas both codes start at fixed core `{0,1}`. The excess is therefore a
diagnostic comparison, not an established defect maintenance cost. Two
other references make the comparison choice visible:

| Comparison ensemble | Current per attempt | Code excess |
| --- | ---: | ---: |
| Stationary co-moving `q=1` ensemble | `5/33` | `1/66` |
| Fixed-core `q=1` conditional reference | `23/150` | `1/75` |
| Full fixed-core transported reference | `341/2250` | `17/1125` |

The last two fixed-core ensembles are **not stationary** under the bare
unconditioned C/E law. Thus they cannot be silently used as its equilibrium
action measure. A hypothetical cost `epsilon` per additional current would
yield `epsilon/66` per attempted repair against the stationary comparison.
Only after supplying an attempt interval `delta_tau` in a calibrated observer
clock could that be called power `epsilon/(66 delta_tau)`. Neither expression
is a rest-energy measurement, and no energy unit or proper-time action is
supplied by the current count.

## An exact ambiguity even after demanding native stationarity

On the recorded source pair `x` and transported classical core `c`, the
existing positive measure is `pi(x,c)=sigma_c(x)/66`. Every complete C/E
record block preserves `q=|x intersection c|`. Its connected sectors have
`2970`, `1320`, `66` states, with probabilities `(2/5,5/9,2/45)`. For any
positive factors `f(q)`, the normalized measure

```text
pi_f(x,c) = f(q(x,c))*pi(x,c) / sum_{x,c} f(q(x,c))*pi(x,c)
```

is stationary for **the same reduced diagonal source/core population law**
under uniform-seam C/E repair. The growing retained archive is not assigned
a stationary law. Every
movement edge and its reverse have equal probability, and both endpoints
have the same `q`. Self-loops add no constraint. The receipt checks all
`300 × 66 = 19,800` edge/core pairs for `f=(1,1,1)` and `f=(1,2,1)`,
plus the complete raw-record charge premise from the reference-balance packet. The
second factor changes sector probabilities to `(9/35,5/7,1/35)` while the
conditional distribution, and all native records, inside sharp `q=1` stay
the same. For any `U(q)`, the native gradient form is zero because
`U(q_after)-U(q_before)=0` on every record.

If each stationary reference were **separately declared** Gibbs at one
common inverse temperature `beta`, doubling the relative `q=1` weight
would shift the modular contrast `K(q=1)-K(q=0)` by `-log(2)` and the
associated declared energy contrast by `-log(2)/beta`. This is a precise
identifiability limit for the bounded family of sector-only potentials and
the native current form. It does **not** say arbitrary `U(q)` is a physical
energy, nor does it erase a separately specified modular reference: once a
particular Gibbs reference and `beta` are fixed, its relative modular gap
is fixed. Native stationarity and the active code's own records do not select
that reference or identify its energy with an action.

This is the actual-model reason the normalized symmetrizer result in
[SourceActionMeasure](../../Lean/Geometry/SourceActionMeasure.lean) must be
applied **within each connected `q` sector**: the full joint graph is
disconnected, so sector weights are free. The formal Gibbs identity in
[GibbsReferenceEnergyIdentification](../../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean)
requires the declared `H` and `beta`; it cannot select them from native
stationarity. These statements coexist with the separately specified
transported clock reference, whose sector weights are one
possible declared choice.

## Corrected-code work and the missing inertia test

The source-executed recovery and counted routing below belong **only to the
two-support logical code** on seam `{0,5}`. All positive raw flags are
retained; the code has zero abort probability. Source encoder cost is
`28,635` native pulses once per preparation. The repair word and
record-selected recovery have these exact per-attempt counts:

| Law | Repair pulses | Mean recovery pulses | Mean repair plus recovery |
| --- | ---: | ---: | ---: |
| C | `386,739` | `52,670` | `439,409` |
| E | `386,742` | `89,999/2` | `863,483/2` |

This one-seam executed protocol is a different schedule from the uniform-seam
activity diagnostic above. These numbers cannot be added into one implemented
uniform-seam maintenance action or used as a compiled all-seam recovery cost.

The classical core value changes on reset with probability `1/4` after
recovery. One 15-bit status/old/new core entry is appended per non-abort
attempt, including entries with no value change. These are source compiler
resource counts, not action, joules, dissipated work or an energy per
calibrated observer-clock interval. The two-support code's source packet also shows that omitting
the explicit core reset leaves only `7/8` in sharp `q=1` after correction.
These one-seam counts are the boundary of this diagnostic. The
separate screened-code protocol executes its own encoder, routing, recovery
and core reset, with uniform-seam mean native-plus-recovery counts
`37,499,681/60` for C and `36,870,437/60` for E per attempted repair.
Its [resource ledger](../screened_recovery/README.md#core-transport-and-service-cost)
and [full corrected cycle](../screened_transport/README.md#the-full-corrected-cycle-has-a-flat-logical-channel)
belong to that screened protocol and are separate from the one-seam values above.
The retained `source_pulse_recovery_executed: false` field describes this
packet's screened-code implementation boundary, rather than the status of
the screened construction. Neither packet assigns rest energy to pulse counts.

The separate 240-port product-field control has an ideal rest gap and an
independently differentiated radial-symbol curvature satisfying its own
mass shell. It is not this 66-state corrected code. For this code there is
no implemented momentum-parametrized dispersion, independent kinetic
curvature, proper-length action, or source stress tensor. The
[InternalEnergyInertia](../../Lean/Geometry/InternalEnergyInertia.lean)
theorem distinguishes coupled from additive action shapes;
[ProperTimeInternalAction](../../Lean/Geometry/ProperTimeInternalAction.lean)
states refinement and proper-time premises needed to choose the coupled
form. Neither premise follows from counted pulses. Recorded occupation
currents satisfy the finite continuity map of
[RegionalContinuity](../../Lean/Screen/RegionalContinuity.lean), but the
stress realization and refinement in
[WardLimitManifest](../../Lean/Dynamics/WardLimitManifest.lean) are open.
Thus `E=mc^2` is a same-candidate test after energy, kinetic response,
stress and a common calibrated observer clock have been independently
specified; it cannot be used to define the code's mass from these counts.

```sh
python -m model.active_code_energy.energy_bridge --check
python -m pytest -q -W error model/active_code_energy/test_energy_bridge.py
```
