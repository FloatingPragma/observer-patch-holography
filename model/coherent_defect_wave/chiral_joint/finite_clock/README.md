# Finite clocks for interacting records

An exactly unitary rational function of the declared Hamiltonian implements
a finite clock. Its error relative to the same Hamiltonian's exponential
is bounded on the full state space. The clock preserves the certified
neutral arrival on the Peter–Weyl source. A separate exact execution on the
complete 72-state periodic physical block gives both loop outcomes and
their complementary record.

The [proof](PROOF.md) states the algebraic clock service and its scope.
This service supplies no relativistic dispersion, physical clock scale or
observed particle identification.

```sh
python -m model.coherent_defect_wave.chiral_joint.finite_clock --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/finite_clock
```
