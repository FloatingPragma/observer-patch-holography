# Cooling the full spatial field

The spatial gauge-coordinate model admits deterministic stationary states
as close to its finite ground energy as requested, without assuming an
excitation gap. The construction keeps all sixty fermionic modes per
cell, its quantum links and the conjugate reference. It gives an explicit
finite preparation budget and a fixed-cutoff spatial ground-state
subsequence. Its resources can grow exponentially with volume.

The exact 208-state native control reaches certified mean ground-energy
excess below 0.157 source units with every occupation and flux retained.
These are energies of the declared Hamiltonian; a relativistic particle
interpretation requires its field and detector limit.

[Proof and domain](PROOF.md), [implementation](ground.py), and
[complete charged detector records](../particle_ground_spectral_records/README.md)
give the source-to-preparation-to-measurement connection.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m model.particle_spatial_ground_state
python -m pytest -q -W error model/particle_spatial_ground_state/test_ground.py
```
