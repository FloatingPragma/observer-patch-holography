# Testing particle candidates

The [candidate table](RESULTS.md) groups eighteen named comparison targets
by their source representations. The calculations test conserved charge,
internal mixing and which detector can create each finite source channel.
They supply discriminators for source excitations. Physical identification
requires their propagation and interaction records in the same field
description; no GeV mass or generation assignment is inferred here.

```sh
python -m model.particle_candidate_confirmation --write
python -m model.particle_candidate_confirmation --check
python -m pytest -q model/particle_candidate_confirmation
```

The consumer keeps its producer receipts separate, with deduplicated source
hashes. It does not repeat the full spectral certificates or change an
accepted physical dependency.

The [confirmation steps](CONFIRMATION_STEPS.md) attach six physical tests
to each target within the common field-and-particle dependency map.
The [labelled catalogue](CATALOGUE.md) groups protected objects and source
excitations separately, with their original identifiers and domains.
The [per-object mass paths](MASS_PATHS.md) state the first missing attachment
for every protected entry and the finite results that can be reused.

## Charged elementary families

The [local operator calculation](../particle_candidate_charged/PROOF.md)
tests the actual sixty-mode spin and gauge coefficients. Bare family mixing
obeys exact conserved-occupation selection rules. A spin-scalar composite
with Higgs-conjugate internal charges gives an explicit dressed charged
singlet channel. A complete compensated interference experiment proves
positive original-state overlap with the bare charged-hole channel. Its
full-occupation residual also gives original weight greater than `1/40`
in the reported pair-transfer interval `[77/8,271/8]` source units at
member zero, using the same complete energy clocks. The opposite charge
and connecting flux are included. Its physical vacuum, chirality and mass are separate
requirements; a single matched internal row does not distinguish three
generations.

The [time-filtered pair construction](../particle_candidate_charged_propagation/PROOF.md)
retains original probability above `3/80` in `[311/32,1081/32]` at member
zero and admits this local seed to the canonical energy/momentum algebra.
The bulk-aperture member-zero apparatus is inactive, so the local bound
does not prove a positive continuum residue.

The [separated endpoint calculation](../particle_candidate_endpoint_threshold/PROOF.md)
uses all shortest cube paths of lengths 1, 2 and 3. Its original complete
filtered dressed probabilities exceed `4/25`, `2/25` and `7/250` in the
reported positive bands. These bands overlap; the pair data do not
determine a charged endpoint threshold or mass.

## Neutral composites

The [composite identity counter](../particle_candidate_composites/PROOF.md)
measures the change of `S=N_-4+N_2+N_6` in the selected intervention. Seven
record values distinguish the six signed weak, color and even-composite
channels and failure, preserving their original energy histories. The
counter records an occupation increment, without assuming a physical
baryon or lepton name. Multiple excitations outside the admitted alphabet
can alias modulo seven.

The [spatial/final-energy experiment](../particle_candidate_composite_spectrum/PROOF.md)
has nonzero original neighbor-versus-nonneighbor contrast in all forty
components. Its complete final-energy records have explicit, extremely
small probability bounds. They supply no rest mass or two-point dispersion.

The [original two-point moment bound](../particle_candidate_neutral_dispersion/PROOF.md)
uses occupied-pair estimates on the weak-pair/singlet trion. Within each
G1 or H multiplet, the four Walsh spectral means lie in `(16/3,38/3)`,
span less than `5/9` and have ratio less than `53/48`. They may coincide.
A complete two-clock estimator has controlled bias; its histogram and
the signed neighbor coefficient are uncomputed. These are spectral means,
not physical mass ratios or a dispersion relation.

## Neutral scalar candidate

The [radial response](../particle_candidate_scalar_gauge_response/PROOF.md)
of the Higgs-conjugate composite has total detection probability at least
`3481/160000` throughout the source family and above `383/500` at member
zero. A separate two-time experiment has nonzero original dynamical
response. Total detection includes elastic records; a surviving inelastic
branch and Higgs-phase interaction dictionary require separate proofs.

The [resolved positive transfer](../particle_candidate_scalar_spectrum/PROOF.md)
has original weight at least `c*9/(128*a^2*D^2)` in `[1/40,D]`, with
`rho >= c I` and full energy enclosure width `D`. Complete finite clocks
retain at least three quarters of that floor in an expanded band. Exact
elastic weight has the uniform floor `3481/160000`; the positive inelastic
floor can vanish and D grows with capacity. The enormous clock resources
are specified symbolically; their histogram is not executed.

## Gauge-sector candidates

The [complete loop menu](../particle_candidate_bosons/PROOF.md) resolves
`A1g, Eg, T1u, A1u, Eu, T1g`: twelve source components and one failure
outcome. Every component has a certified interacting transition with the
original rare-sector probability retained. Even scalar and vector probes
cannot excite the isolated `A1u` line from an `A1g` input. This separates
candidate detector channels, without assigning intrinsic Lorentz spin or
physical photon, W, Z or Higgs identities.
