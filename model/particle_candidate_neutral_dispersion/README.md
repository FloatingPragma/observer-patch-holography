# Neutral two-point energy constraints

The exact occupied-pair derivative bounds the unresolved original signed
neighbor moment of the two G1 and four H gauge-singlet trion probes. On
member zero, each multiplet's four normalized Walsh energy centroids form
an arithmetic progression with total spread less than5/9 source units.
Each centroid lies between16/3 and38/3; their ratio within one multiplet
is less than53/48. These are spectral means, without a particle-mass or
conserved-momentum interpretation.

The complete source-clock construction measures the signed neighbor
coefficient with a certified deterministic bias below1/2048. It retains
all failures, both energy records and the original reference. Its finite
resource plan is calculated; the original full histogram is not executed.

See [the proof](PROOF.md) for the exact CAR identities, occupation bound,
normalization, clock error and limits on physical interpretation.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_candidate_neutral_dispersion --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest model.particle_candidate_neutral_dispersion.test_model
```

`results.json` binds the fixed source, selected polynomial state, reference
and clock. It supplies no physical particle assignment or GeV mass table.
