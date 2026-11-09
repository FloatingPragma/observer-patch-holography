# Physical gauge states on the spatial field

The full spatial Wilson source admits a deterministic preparation supported
in the common kernel of its local Gauss constraints. A finite polynomial
of the same Hamiltonian reaches any requested positive mean ground-energy
tolerance, with all allowed physical states and references retained.
At fixed spacing and link cutoff, these states supply a ground-state
subsequence on the gauge-invariant observable algebra.

The native quantum-flux control has a complete 78-dimensional physical
carrier inside its 208-dimensional ambient carrier. Its direct preparation
has certified mean ground-energy excess below 0.105 source units. Charged
measurements retain explicit compensating boundary references; a bare
charged field does not preserve the source-only physical space.

[Proof and precise domain](PROOF.md), [implementation](gauss.py), and
[the common field source](../particle_spatial_ground_state/PROOF.md)
separate the finite construction from physical refinement and particle
identification.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m model.particle_gauss_ground_state
python -m pytest -q -W error model/particle_gauss_ground_state/test_gauss.py
```
