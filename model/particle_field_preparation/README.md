# Preparing a gauge field with a controlled Wilson gap

A finite polynomial of the source's own link coordinates prepares a
gauge-invariant state with an explicit upper bound on small Wilson gaps.
The [proof](PROOF.md) controls its normalization and electric energy, then
uses the same interacting Hamiltonian clock to bound how long the gap
estimate survives. The construction chooses a preparation; it does not
identify that preparation as the physical vacuum.

The executable checks a finite seed on the actual one-cell coordinate
field. Its character preparation has an exact 1,331-dimensional code and
1,728 complete electric-irrep records, including every zero complement
record. It also verifies a native weak-link configuration where the
Wilson gap closes, showing why a state condition is necessary.

The higher-degree gap recipe is a finite representation algorithm with
an explicit resource bound. Its conservative certificate is enormous and
its small-plaquette-cost control has a short certified dynamical interval.
Those large coefficients and the preparation tape are not materialized.
The small seed has no useful
small-gap guarantee by itself.
The executable word codes and evolution cutoff are finite; the Haar
coordinate Hilbert space is the analytic parent used for their error bounds.

The [conserved-energy construction](REFINEMENT.md) gives a distinct fixed
high-plaquette-cost program whose gap-probability bound holds for every
time. It uses the existing plaquette-admissibility theorem and the same
finite polynomial preparation. Its coupling and resource bounds are
explicit and enormous; the result is not a uniform relativistic limit.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_field_preparation --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_field_preparation --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/particle_field_preparation
```

All calculations use the [declared coordinate field](../particle_observer_outputs/FIELD_ROUTE.md)
and its periodic spatial boundary condition. The result supplies no
uniform interacting relativistic limit or physical particle mass.
