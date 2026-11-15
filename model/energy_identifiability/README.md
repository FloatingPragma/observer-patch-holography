# Energy offsets before mass fitting

This diagnostic identifies generator energy differences that a declared finite
protocol can distinguish **before** attaching a mass to a protected mode or fitting a
dispersion curve. It adds an operational use of the existing
[Gibbs energy identities](../../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean)
and [conserved-sector result](../conserved_sector_action/README.md).
It supplies no physical clock, native field identification, particle species
or mass. The source/core carrier and the interacting field are distinct.

Run the exact controls and, optionally, the bounded full-field comparison:

```sh
python -m model.energy_identifiability.certificate
python -m model.energy_identifiability.certificate --field
python -m pytest -q model/energy_identifiability/test_certificate.py
```

## Which energy offsets can observations distinguish?

Fix the carrier, accessible algebra, reference, generator and clock convention
independently of the fitting data. For diagonal generators `K` and `beta H`,
their derivations agree on a matrix unit `E_ij` precisely when

```text
K_i - K_j = beta (H_i - H_j).
```

Thus `K - beta H` is constant on each connected component of a graph whose
edges are **admitted coherent matrix-unit comparisons**. Diagonal population
observations or classical movement edges alone do not supply those comparisons.
For `v` vertices and `c` components the incidence rank is `v-c`; there are
`c` offset freedoms, of which `c-1` are relative. A spanning tree determines
the relative constants in a connected graph. Further edges test cycle
consistency. `reconstruct_offsets` checks this exactly for rational differences
and rejects inconsistent cycles. Choosing zero independently in disconnected
components does not measure a common energy origin.

The same statement applies to blocks when the full generator inside each
block is known up to a scalar: one coherent comparison between two
blocks fixes their relative scalar. This is a conditional reduction in the
number of comparisons, not a construction of the comparisons. A supplied
positive `beta` is also a premise. Determining its value from a slope requires
a nonzero admitted energy contrast, as in the existing finite Gibbs theorem.

## A separate clock does not remove the native sector ambiguity

Let `Q_q` be the transported source/core charge projectors. Suppose:

- the initial joint source/clock/reference state is block diagonal in `q`;
- each evolution commutes with these projectors;
- every instrument branch, including feedback and clock gates, intertwines
  the same projectors; all classical outcomes are retained.

Keep the **same initial preparation** and replace `H` by
`H + sum_q a_q Q_q`, for arbitrary scalar sector offsets.
On each block the additional evolution is only a scalar phase. It cancels in
every density update. Induction gives the **same complete quantum/classical
output state at every finite history length**, including arbitrary correlations
within each block with a separate clock or untouched reference. Adaptive
choices based on the same retained history preserve this conclusion. The
time parameter here is the declared generator parameter, not calibrated seconds.

The native check verifies `Q_q K_r = K_r Q_q` on every raw C/E source/core
matrix entry: 261,360 entries per law. This operator check extends the existing
stationary-sector ambiguity to a protocol-level energy-offset ambiguity.
It does not assert that every operation in the larger source gate grammar
preserves `q`. Coherent gates or ancillary couplings that mix sectors violate
the premise. The exact clock control uses sector-dependent clock rotations,
retains eight three-step histories and finds identical outputs under a
nontrivial sector phase. The countercontrol prepares and reads a cross-sector
`|+>` state: a relative phase of `pi/2` changes its `+` probability from `1`
to `1/2`. That protocol can observe an offset because it supplies the missing
comparison. Neither control supplies a physical energy calibration.

This is a statement about those generator interventions at fixed preparation.
It does not forbid inferring relative sector energies from measured Gibbs
populations **if** a global Gibbs preparation law and its temperature are
independently established. Changing `H` in that law also changes the prepared
sector weights; that is a different experiment. The native conserving law
does not select those relative stationary weights, so declaring a Gibbs law
there supplies the missing link rather than proving it.

The initial block-diagonal condition is needed for equality of **full** output
states. With initial cross-sector coherence, its phase can change even when
sector-preserving observable probabilities cannot reveal it. This diagnostic
does not identify those two output contracts.

## The complete native tape does not give the localized coherence an energy

For the actual [localized mixed-state coherence](../native_localization/README.md),
`rho` and its faithful reference `pi` have identical initial diagonals. Each
complete native raw branch is a weighted partial injection in that basis.
Products of these maps remain partial injections, and source-independent
branch weights multiply both states identically. Therefore every history
has the same branch probability and diagonal under `rho` and `pi`.

At history length `n`, use **the same instrument pushforwards**
`R_n = Phi_n(rho)` and `T_n = Phi_n(pi)`, including the entire classical tape.
`T_n` is diagonal and faithful on its image support. Consequently

```text
Tr[(R_n - T_n) (-log T_n)] = 0.
```

This holds for every finite `n`, including record-adaptive selection among
the admitted native maps. It is an extension of the existing source energy
statement to the full archive; no stationary density on a growing tape is
assumed. Zero-probability records have zero support and are retained as zero
blocks. Positive-probability conditional branches satisfy the same identity.
The certificate checks all 720 encoded raw slots and their partial-injection
property. A test deliberately inserts a coherent merger and rejects it.
Additional coherent analyzers are outside this induction.

This is an energy statement, not absence of information: the initial state
differs from its reference and has positive relative entropy. In general
`D(R_n || T_n) = S(T_n) - S(R_n)` here; retaining or discarding the archive
need not preserve relative entropy even though the preceding trace identity
is zero in both cases. It cannot be converted into a positive mass by
relabeling the information deficit or native-attempt decay rate.

## The actual 64-fiber field control

The control uses all 8,128 states of the existing fixed-`N=2` interacting field,
with the same supplied common `beta`. Write its fiber partition functions as
`Z_k`. The global reference has modular generator `beta H + log Z`; normalizing
each fiber separately instead gives `beta H_k + log Z_k`. Intrafiber generator
responses agree, while cross-fiber coherences detect `log Z_k - log Z_l`.
The complete numerical spectrum gives a nonzero partition-log spread. Given
the intrafiber response, a mathematical tree of 63 cross-fiber matrix units
suffices to determine all relative constants. Its graph rank is exact;
partition logs and residuals are floating diagnostics. The existing moving
packet retains cross-momentum coherence, so replacing its global reference
by separately normalized fibers is consequential.

**The tree matrix units are not certified source preparations/readouts.**
The diagnostic specifies what an implementation must compare; it does not
make the global Gibbs reference native or supply its preparation.

There is a further independent boundary for mass attachments. On this fixed
`N=2` carrier, adding `c N` to `H` adds `2c I`; it leaves every field channel
and the normalized Gibbs reference unchanged. Connecting all momentum fibers
therefore does not measure the two-particle rest-energy offset relative to a
vacuum sector absent from the carrier. Calibrated dispersion curvature may
constrain an effective mass under declared kinematic assumptions;
that is a different contract from measuring an absolute vacuum-sector gap.
The existing eta-pair/dark-state contrast `U` is an internal `N=2` energy
difference and survives this scalar shift. Its measurability is unaffected;
it is not a comparison to the vacuum.

## Relation to the finite proof stack

`GibbsReferenceEnergyIdentification.lean` proves the equal-mass
partition-constant cancellation, additive-energy gauge and slope uniqueness
when a nonzero energy contrast exists. `ConservedSectorAction.lean`
proves invariant sector reweighting and zero native charge-gradient cost.
The additions here are the accessible-algebra graph criterion, the full
clocked-history induction, their actual raw-map checks and the full-archive
energy application. These matrix/instrument arguments are given above and
tested here; they are not represented as Lean theorems.

`ModularEnergyAdditivity.lean` requires the declared reference factorization
or its explicit binding correction; adjoining a clock does not select it.
`CurlStoneClockBridge.lean` similarly requires its operator map and clock
calibration. Neither theorem converts a scalar-sector gauge, a killed
coherence mode or a supplied field spectrum into a native physical mass.
