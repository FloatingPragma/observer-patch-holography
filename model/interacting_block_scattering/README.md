# Finite collisions between neutral excitations

The paired native blocks support a controlled two-packet interaction.
The three-dimensional theorem uses an ideal free-symmetry-forbidden
channel. The [proof](PROOF.md)
bounds the interaction, finite-capacity corrections, packet preparation
and detector errors before squaring the signal. Two three-dimensional
members have positive collision probabilities in a fixed momentum/time
window. An independent 1024-state calculation gives complete records
for a five-block cycle. Its response is compared with full free evolution
of the same admitted input and detector; both responses are nonzero.

The [receipt](results.json) contains five observation times, hit and
complement probabilities, and native block occupation marginals with
rigorous intervals. These are quantum records, with explicit preparation
errors. They describe no classical particle trajectories.

From the repository root, with `model/requirements.txt` installed:

```sh
python -m model.interacting_block_scattering --write
python -m model.interacting_block_scattering --check
python -m pytest -q model/interacting_block_scattering
```

[Source ancestry](PROVENANCE.md) identifies the immutable parents and
declared preparation. The interaction vanishes in the controlled
infinite-capacity limit. No asymptotic scattering matrix, interacting
relativistic ultraviolet limit, Standard Model identification or efficient
implementation is asserted.
