# Causal influence of the interacting field

A local intervention in the declared gauge-coordinate program has an
explicitly bounded effect on distant detector records. The bound keeps
all fermion occupations and unbounded electric fluxes, and is independent
of the total lattice volume. It applies to the same Hamiltonian and clock
as the [finite field construction](../particle_observer_outputs/FIELD_ROUTE.md).

The [proof](PROOF.md) derives the local constants from the actual native
spin projectors and the retained overlap polynomial. It then bounds both
remote influence and the error from simulating only a finite neighborhood
of a complete local experiment. The construction reuses the repository's
CAR locality and complete-instrument arguments; it does not identify a
causal influence bound with Lorentz invariance.

For one declared polynomial control, a remote gauge-invariant operation
32 cells away changes the complete origin detector distribution by less
than `10^-7` after duration `1/64`. This is an analytic certificate, with
an exact rational exponential enclosure, not execution of a large Fock
array. A separate exact counterexample shows why the total representation
cutoff needs its own detector error: that global compression can make
two distant link operators fail to commute.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_field_causality --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_field_causality --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/particle_field_causality
```

The fixed-spacing locality theorem does not give a uniform interacting
relativistic refinement, physical chirality, or a particle mass assignment.
The [source/output boundary](../particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md)
states the distinct spectral and physical identification requirements.
