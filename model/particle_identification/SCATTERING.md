# From finite responses to collision records

A CERN-like comparison requires prepared incoming beams, an interaction
region and complete outgoing detector records. Addition frequencies and
a nonzero self-energy do not by themselves specify that experiment.
The finite gauge–Higgs source supplies interacting response and complete
Ramsey records; its two vertices do not provide spatially separated beams
or an asymptotic outgoing region.

## Reusable calculations

The [contact scattering theorem](../contact_scattering/ONSHELL.md) constructs
complete two-fermion lattice wave operators and a unitary scattering matrix.
Its [finite source control](../contact_scattering/README.md) is an interaction
echo on two sites, including an untouched reference and every export failure.
The echo is an executed interaction witness, not an executed asymptotic collision.

The [moving contact construction](../moving_contact_scattering/PROOF.md)
has a nontrivial three-dimensional two-particle continuum scattering limit
at bounded total momentum, with the same critically tuned contact coupling
and clock as its moving bound branch. All eight edge valleys are retained.
This is nonrelativistic. Its additive extension does not make the empty Fock
state a ground vacuum: separated bound pairs lower the energy.

The [scattering packet readout](../scattering_packet_readout/PROOF.md)
connects that theorem to an explicit incoming source packet, finite collision
time, spatial cutoff, preparation and complete interference detector.
Cook tails, packet truncation, coefficient approximation, preparation, finite
volume, product formula and native port errors have separate budgets. The
result distinguishes interaction from the matched free experiment by more
than `3/10` in a retained detector event. The very large fine-spacing program
is a finite constructive schedule, not an executed full collision simulation.
Its separate five-qubit control executes all 32 rows and checks the free
control, clock orientation and complete native outcomes.

These controls can be replayed without building another scattering package:

```sh
python -W error -m pytest -q model/contact_scattering/test_collision.py
python -W error -m pytest -q model/scattering_packet_readout/test_native.py
python -W error -m model.scattering_packet_readout --check
```

The first executed control to reuse is the five-qubit echo and its matched
free run. Replacing its preparation by independently
prepared separated incoming packets and reading outgoing momentum bins
requires a corresponding preparation, propagation and detector error bound;
the existing echo must not be renamed such an experiment.

## A common interacting collider model

The physical route joins one selected source state and Hamiltonian to a
spatial family, prepares incoming wave packets in that family, evolves them
through a collision, and retains every outgoing detector result. Failure,
soft radiation and unobserved channels contribute to the complete law.
Bounded detector bins and finite resolutions are legitimate effective
observables; each approximation must bound their joint record error.

A controlled infinite-volume or finite-time window must prevent boundaries
and recurrences from replacing the intended experiment. A relativistic
effective description additionally needs the same clock and spatial metric,
the relevant symmetry and locality bounds, surviving interactions, the full
regulator-sector census, and a quantitative domain of validity. A formal
ultraviolet limit is not necessary for a finite-resolution comparison, but
a free dispersion bound is insufficient for an interacting collision law.

Particle identification then depends on the physical channel. Stable neutral
particles may admit isolated mass shells. Charged sectors require infrared
dressing or inclusive threshold observables; a universal isolated-pole test
would exclude legitimate charged physics. Unstable particles require
resonance definitions and line-shape conventions. Confined quark masses
are parameters in specified renormalization schemes and scales, inferred
from hadronic observables. These distinctions are supported by
[Buchholz's charged-sector result](https://www.sciencedirect.com/science/article/pii/037026938691110X)
and the [Particle Data Group quark-mass review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-masses.pdf).
The corresponding cross sections, branching
fractions, angular distributions and detector response must follow from the
same model before comparison with measured collider records.

Reading simulator specifications is allowed throughout. It fixes which
source calculation is intended; it does not substitute a scattering
amplitude, state preparation or detector map for a missing one.
