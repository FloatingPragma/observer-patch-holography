# Native localization, finite contrast lifetime, and exact-code limits

This packet constructs a faithful localized coherence on the **same biased
66-state carrier and product reference used by the continuing record clock**.
Its complete exterior agrees with the reference, including the joint exterior
and native record history. It has a nontrivial modular phase and a rigorously
computed finite contrast lifetime. It is a decaying, immobile coherence with
zero modular-energy excess. It does not identify a protected moving particle,
a topological charge, or a mass.

The closest prior results are the
[signed-difference quantum-sector reduction](../federated_quantum_sector/README.md),
the [general native quantum mixing theorem](../observer/QUANTUM_MIXING.md), the
[biased clock's ten-spectator response](../local_repair_family/modular_clock.py),
and the [arbitrary-code capacity theorem](../common_source/record_code.py).
Signed-difference preservation and the existence of a killed spectator kernel
are predecessor results. The additions here are the particular faithful state
with an exactly invisible exterior, its complete-history localization and
resolvent lifetime, a new source preparation on the existing register, the
same-state energy check, and general exterior/clock consequences of the
arbitrary-code theorem. None assumes a sharp core occupation.

## The actual state and reference

The source law has positive endpoint weights
`theta=(25,25,8,8,8,8,8,8,8,8,8,8)`, all 30 seams, and binary total load two.
The fixed reference is

```text
pi_{ij}=theta_i theta_j/7505.
```

Choose the actual native seam `A={1,7}` as the core. For each spectator
`k` outside this core, replace the reference block on
`{|1,k>,|7,k>}` by

```text
rho_k = (theta_k/7505) [[25,8],[8,8]].
```

Every other matrix entry is unchanged. The determinant of each displayed
block is `136(theta_k/7505)^2`, so the entire state is faithful and normalized.
Its diagonal equals `pi`. Its off-diagonal difference has zero partial trace
over the core. Consequently the **full exterior state**, including all its
many-port correlations, equals the exterior of the same faithful reference.
The global half trace distance is `776/7505`, approximately `0.1034`.

This is one known mixed state with a coherence parameter. It is not a qubit
encoding accepting arbitrary logical/reference inputs. Its construction
therefore does not contradict the exact-code restrictions below.

## Complete native history and a finite lifetime

A complete native C/E record maps matrix units by a partial injection.
Equal records preserve `x-y`. Every surviving `|1,k><7,k|` is consequently
another `|1,l><7,l|`: the common spectator can move, while the two sites of
quantum disagreement stay fixed. The diagonal histories of `rho` and `pi`
are identical. The difference in each unnormalized raw-history block has
zero exterior partial trace. This proves, inductively for every finite
history length, that the joint **full native tape plus exterior** is identical
for the two initial states. It does not equate a conditioned exterior with
the unconditioned reference. It also does not preserve the preparation
archive: the declared preparation mixture has its own retained classical
record.

The executable checks all 720 raw seam/law/coin/value slots, including 240
abort slots and 180 valid zero slots, against an independent raw-menu source.
It checks 7,200 transfer matrix units; 920 survive individual slots. The
induction includes every history without selection of successful records.

Let `v_n(k)` be the surviving real coherence coefficient after `n` actual
native attempts. The full, uniform-seam scheduler gives `v_{n+1}=v_n Q`,
where the exact ten-by-ten matrix `Q` is retained in
[results.json](results.json). Its entries include the biased coins. It is
reversible with weights `theta_k`, and a positive diagonal-dominance bound
places all eigenvalues strictly between zero and one. The contrast is

```text
D_n = sum_k v_n(k),
T_contrast = sum_{n>=0} D_n/D_0
           = theta^T (I-Q)^(-1) 1 / sum_k theta_k.
```

This half trace distance is the same with the entire native tape retained
or discarded: all surviving blocks have the same coherence sign and phase.
`T_contrast` is an integrated visibility in **native attempt count**, not a
particle survival time measured by a physical clock.

| Native law | Exact or approximate integrated relative visibility |
| --- | ---: |
| C | `86947767665473931268090/12912620562057708632977` ≈ 6.73354934 attempts |
| E | `10/3` attempts; `D_n/D_0=(7/10)^n` exactly |
| Equal C/E mixture | `1106546929049854606584045/248438522107974898878437` ≈ 4.45400705 attempts |

The receipt includes rational eigenvalue enclosures, all ten floating
spectral diagnostics, exact resolvent potentials, and visibility at attempts
0, 1, 2, 4, 8, 16 and 32. The response is a positive mixture of real decays.
There is no oscillatory native-attempt dispersion or moving coherence core.
Deleting all incident seams makes this coherence survive, but violates the
positive all-seam contract; a test retains this countercontrol.

## Clock readout and the energy distinction

The derived reference flow `U_s=pi^(-is)` multiplies each oriented
coherence by `exp(-i s log(25/8))`. The parameter `s` is dimensionless modular
parameter; it has not been identified with native attempt count `n` or a
calibrated physical time. Modular covariance of the actual branches permits
this phase before, between, or after repairs.

After `n` native attempts, the supplied bare-core Ramsey analyzer and one
additional native current read give the **complete-trial** prediction

```text
E[j] = -1649/15010 + D_n cos(s log(25/8)).
```

C and E have the same current marginal for this terminal read, while their
quantum instruments differ. The entire three-outcome law is in the receipt.
The analyzer is required: native records alone are diagonal effects and do
not reveal the coherence. The same state has

```text
Tr[(rho-pi)(-log pi)] = 0,
D(rho||pi) = S(pi)-S(rho) > 0.
```

The displayed numerical entropy deficit is approximately `0.0587094` in
natural-log units. A nonzero clock phase and information deficit therefore
do **not** establish nonzero excitation energy or mass.

For comparison, the existing compiled local rotation `R pi R*` has initial
modular-energy excess `(1649/15010) log(25/8)` and zero initial entropy
change. Its initial exterior equals the reference, but the population
perturbation reaches the exterior after one native attempt. The receipt
propagates its exact exterior distance and modular-energy coefficient
through 16 attempts on the same source. These are different states; their
favorable properties cannot be combined into a single particle claim.

## What arbitrary exact codes cannot fix

The existing full-78-state, arbitrary-complex-code theorem says that exact
flagged correction for every positive native seam first excludes all double
occupations, then excludes every one of the 30 native-edge pair states.
It also forces every compressed number operator `V* n_i V` to be scalar.
Positive endpoint biases multiply effects by nonzero numbers and leave this
proof unchanged. It covers arbitrary CPTP flagged recovery, not only the
existing rectangle codes.

For a proposed two-port core, group all pair states by their exterior
configuration. If **every** extension of an exterior configuration to the
core is a native edge, an exact code gives that exterior event probability
zero, while the faithful reference gives it positive probability. Sum these
reference probabilities. Their full 66-core census gives

```text
D(rho_code|outside, pi|outside) >= b_core >= 96/395 > 0.243.
```

The event includes forbidden two-particle exterior edges, single exterior
loads with both core extensions forbidden, and an empty exterior when the
core pair itself is a native edge. The minimum occurs for cores `(0,3)` and
`(1,2)`. The bound is attained over unrestricted diagonal states on nonedge
support; no claim that the attainer is an exact code is made. The support
obstruction also applies to a correctable subsystem with a nontrivial
quantum logical factor and full gauge support: a nonzero rank-one effect
cannot compress to `I_logical tensor B_gauge`.

For the actual clock, let `N_H=n_0+n_1`. On the allowed nonedge support,
`N_H` is a projector because `(0,1)` is a forbidden native edge. Every exact
subspace code satisfies `V* N_H V=nu I`, giving, up to a global phase,

```text
V* U_s V = [1-nu + nu exp(-i s log(25/8))] I,
P(clock leaves the code) = 4 nu(1-nu) sin^2(s log(25/8)/2).
```

This identity holds for every logical input and untouched reference. Clock
invariance requires `nu=0` or `nu=1`; the logical clock is then a scalar
phase. It is not a general prohibition on approximate codes, other
references or carriers, history observables, or supplied parameter-dependent
logical controls. Those are separate construction obligations.

## Source execution and scope

The new state's exact rational mixture uses the factorization

```text
[[25,8],[8,8]] = (25,8)(25,8)^T/25 + diag(0,136/25).
```

The classical mixture and native-address initialization are supplied.
Actual compiled rotations prepare the ten coherent components on all 512
original source addresses, with all 66 preparation columns retained. The
same source compiler then executes the existing modular word and the native
C-then-E dilations on seam `(0,1)`. All 128 two-attempt raw histories are
retained, including numerical abort and padding amplitudes. This fresh
requested-seam replay complements the exact 30-seam scheduler calculation;
it is not a pulse replay of all 30 routes. The implementation reports a
floating error diagnostic with gate `1e-8`, not a hardware or outward
interval certificate. Proper-code load/export, pulse actuation, classical
selection, storage and record export are declared services.

Run:

```sh
python -m model.native_localization.experiment --write
python -m model.native_localization.experiment --check
python -m pytest model/native_localization/test_localization.py
```

Replay compares every retained field, including current probabilities, mean
readouts, spectra, source budgets and source hashes. Counts, types and exact
fractions compare exactly; floating diagnostics have a fixed `2e-10` replay
tolerance and nonfinite values are rejected. Both CLI options accept an
optional receipt path. Hostile receipt tests cover changed observables,
source counts, budgets, physical status and nonfinite values.

The package discharges no canonical physical premise and does not close
CP03 or CP04. It removes a concrete uncertainty: an exterior-matched,
clock-sensitive native coherence can be constructed and propagated, but
its full source dynamics supplies neither motion of that coherent core nor
a nonzero energy or protected boundary charge. Exact arbitrary-qubit
correction has the separately quantified exterior and clock restrictions.
