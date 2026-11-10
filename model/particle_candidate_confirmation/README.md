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

## Neutral composites

The [composite identity counter](../particle_candidate_composites/PROOF.md)
measures the change of `S=N_-4+N_2+N_6` in the selected intervention. Seven
record values distinguish the six signed weak, color and even-composite
channels and failure, preserving their original energy histories. The
counter records an occupation increment, without assuming a physical
baryon or lepton name. Multiple excitations outside the admitted alphabet
can alias modulo seven.

## Gauge-sector candidates

The [complete loop menu](../particle_candidate_bosons/PROOF.md) resolves
`A1g, Eg, T1u, A1u, Eu, T1g`: twelve source components and one failure
outcome. Every component has a certified interacting transition with the
original rare-sector probability retained. Even scalar and vector probes
cannot excite the isolated `A1u` line from an `A1g` input. This separates
candidate detector channels, without assigning intrinsic Lorentz spin or
physical photon, W, Z or Higgs identities.
