# How source loops combine

Two neighboring Wilson loops cancel their shared edge and form a six-edge
loop. Opposite face loops keep eight edges. On the actual interacting
source, their complete energy records occupy disjoint certified windows.
The calculation retains all gauge and fermion modes in its complementary
spectral bounds. These are static loop-creation records in a rare
empty-matter sector; moving collisions and a physical particle mass
require a different populated-sector attachment.

The [proof](PROOF.md) gives the full-source cluster and original probability
bounds. Replay with `python -m model.particle_primary_loop_fusion --check`.

## Which collision calculations are executed?

The separate [Weyl/contact experiment](../contact_scattering/README.md)
executes22,531 native pulses on a full four-fermion-mode two-site carrier,
with45,062 load/export events, an untouched reference and every failure.
Its interaction echo has an empty-pointer probability about0.046915;
the free control gives zero. Its supplied action is
`dGamma(h)+U sum n_up n_down`, with `c=1/8`, `U=1` and duration4.
It uses a different action and preparation from the primary Wilson source.

The [moving packet experiment](../scattering_packet_readout/PROOF.md)
constructs incoming packets and a complete finite collision detector with
certified probability separation greater than3/10 in a nonrelativistic
two-particle family. The fine-spacing packet experiment is unexecuted;
its separate32-row five-qubit source control is executed. Neither these
controls nor the static loop-fusion experiment identifies CERN particles
or their physical scattering cross sections.
