# Completion and component proofs

The official executable model is [CCF-LT-2](../SPEC.md), with independent local
capacity factors at every vertex. Its [completion programme](../../docs/COMPLETION_PROGRAM.md)
tracks the current source, observable, continuum and gravitational obligations.
A finite construction result does not discharge physical identification.

The [assumption inventory](ASSUMPTION_CLOSURE.md) and
[current theorem map](../../model/PREMISES.md) distinguish the official source
from conditional results proved for other carriers, states or operations.
The [CCF-LT-2 bridges](../ccf_lt2/bridges/README.md) provide the checked finite
Gibbs, symmetry, coherent phase, local current, mechanical response and
stationarity evidence for the actual Hamiltonian.

This directory also contains independent finite component mathematics:
[recorded observers](OBSERVER.md), [stationary references](observer_reference.py),
[modular access](MODULAR_ACCESS.md), [branch algebras](defect_algebra.md),
[radial excitations](RADIAL_SOURCE.md), [defect geometry](DEFECT_GEOMETRY.md),
and [spectral readouts](PARTICLE_MASS_TABLE.md). Each result retains its stated
carrier, action, reference and clock. A C/E repair carrier, radial source or
supplied field action is not the CCF-LT-2 Hamiltonian. The component class
review in `program_class_v2.json` is not current model adoption; its signed
review is not re-used to accept CCF-LT-2.

Run the current dependency and premise checks from the repository root:

```sh
python tools/build_completion_program_current.py --check
python tools/build_ccf_lt2_claims.py --check
python -m model.particle_joint_reconstruction.run --check
python -m model.completion.assumption_closure --check
python tools/renew_finite_attachment_navigation.py --check
python tools/completion_test_scope.py --run
```

The owning packages continue to execute their independent mathematical tests.
Only a checked same-source identification can transfer a component conclusion
to the official model.
