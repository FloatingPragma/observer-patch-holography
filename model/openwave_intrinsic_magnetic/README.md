# Intrinsic finite magnetostatics

A stationary spin-polarized bound packet generates a nonzero current
through the derivative of its actual gauged Hamiltonian. That current
and the finite Maxwell field are solved together. The weak-field energy
response agrees with the independently integrated current dipole.

The source is a declared positive-energy envelope construction using
all sixteen native kinetic/mass rows, a scalar mass well and a Wilson
regulator. Its empty Fock vacuum defines the energy reference. It is
separate from the PW119 interacting source. Longitudinal electrostatics,
physical screening and a free self-bound electron are not supplied.

The [proof](PROOF.md) distinguishes exact finite existence and operator
identities from the numerical stationary branches. Their bound-energy
magnetic coefficients are 1.47517, 1.68416 and 1.78875 over three grid
refinements. A zero-Wilson control has zero intrinsic current. The
[receipt](results.json) retains both controls and all declared inputs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_intrinsic_magnetic --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_intrinsic_magnetic --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_intrinsic_magnetic
```
