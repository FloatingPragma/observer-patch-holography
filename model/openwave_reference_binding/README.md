# Finite charge binding to a retained reference

The complete finite source assigns a strictly positive energy cost to moving
one unscreened native charge away from its fixed compensating reference.
For the weak-doublet removal channel, the constrained minimum increases by
exactly one source energy unit per graph edge. The charged singlet channel
has strictly positive certified intervals. The calculation uses physical
Wilson-string states and keeps the screened-cloud outcome separate.

These are constrained preparation energies. The [proof](PROOF.md) retains
the actual Hamiltonian, its interacting neutral ground and the reference
cost; it supplies no free charged-particle dispersion or observed rest mass.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_reference_binding --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_reference_binding --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/openwave_reference_binding
```
