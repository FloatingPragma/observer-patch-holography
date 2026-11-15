# Conserved-sector action theorem on the native kernel

This CP-13/#1137 and CP-04/#1140 packet binds the finite
[conserved-sector theorem](../../Lean/Geometry/ConservedSectorAction.lean)
to the same native source/core process as the
[active-code energy calculation](../active_code_energy/README.md).
It extends the two retained reweighting examples to every positive sector
factor, by checking the premises of one quantified theorem.

The model domain is the 66 binary source pairs and 66 classical transported
core pairs. The compatibility condition is specified independently of the
update: `q(x,c)=|x intersection c|`. Every raw C/E record preserves `q`.
The population schedule chooses each of the 30 seams and two coins uniformly
and independently on every native attempt. The complete kernel retains all
zero-current holds. The growing record archive is outside the stationary
state space.

Exact enumeration gives three connected sectors of sizes 2970, 1320 and 66.
The supplied faithful reversible reference has masses `2/5`, `5/9`, `2/45`.
The Lean proof establishes that, for **every** positive `f(0),f(1),f(2)`,

```text
pi_f(x,c) = f(q(x,c)) pi(x,c) / Z
```

is faithful, normalized and stationary under this same reduced population
law. The same statement holds for nonreversible finite kernels whenever
their reference is stationary and their support preserves the supplied
sector map. Conditioning on a fixed sector cancels its factor, so every
instrument applied to that same conditional density has the same record
law. For every potential `U(q)`, the native gradient form is zero.

This uses [SourceActionMeasure](../../Lean/Geometry/SourceActionMeasure.lean):
its reciprocity ratios propagate along walks inside each connected sector.
One overall normalization leaves two independent sector probabilities.
Its uniqueness theorem for a connected full support therefore cannot select
these probabilities. The extension is elementary finite algebra; it adds an
explicit quantified model binding and no law of physics. A separately
declared Gibbs reference and inverse temperature retain their own values;
stationarity and the native gradient do not determine them or a rest mass.

The [receipt](sector_action_results.json) distinguishes exact Python checks
of the actual 4,356-state model from Lean's proof of the general implication.
The frozen [campaign](campaign.json) retains its energy-bridge input
under [an exact predecessor capture](predecessor/energy_bridge_results.json),
byte-identical to that input at the campaign's recorded predecessor commit.
Its specifically named historical path resolves to this capture; every other
campaign source is checked against its live bytes. The current receipt binds
the capture's actual path separately from the current energy-bridge receipt,
so documentation changes do not rewrite the campaign or its evidence.
The full native kernel is not encoded as a Lean value. Mutation tests reject
charge leakage, missing holds, an invalid reference, a disconnected
purported sector and a changed predecessor capture. No clock calibration, physical energy or integrated
completion gate is accepted.

```sh
python -m model.conserved_sector_action.sector_action --check
python -m pytest -q -W error model/conserved_sector_action
cd Lean
lake build Geometry.ConservedSectorAction
```
