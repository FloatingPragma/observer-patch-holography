# Fixed-source photon cutoff

This package bounds complete records of the frozen
`model/native_wilson_maxwell` source as its total photon cap grows on the
same finite cube. Native matter, full complements, preparation, gauge
reference, current and source clock are fixed. The mathematical comparison
uses five unbounded oscillators with the same bounded Peierls coupling.

At `K=3`, `g=mT=1/64`, the vector error is below `1.863e-9`. The
comparison retains an unconditional pair-plus-photon probability lower
bound above `4.75047e-8` and connected binary response above `4.71123e-8`.
All count outcomes, native complements, outer-Gauss failure and comparison
overflow are retained. The actual parent executes 1161 labels; padding
its law with zero-probability overflow labels makes the comparison
alphabet explicit without changing that executed apparatus.

The proof controls literal finite link exponentials, completed photon
Hamiltonians, arbitrary fixed bounded electric/magnetic test histories
with their complete pointer records, and exact local incidence relations
of the gauge-invariant comparison algebra. It does not assert commuting
finite-cap links, a full-native local propagation cone, a spatial
continuum, an interacting relativistic limit or SM identification.

From the repository root:

```sh
python -B -m model.native_wilson_maxwell_refinement --write
python -B -m model.native_wilson_maxwell_refinement --check
python -B -Werror -m pytest -q model/native_wilson_maxwell_refinement
```

`bounds.py` provides exact rational enclosures. `algebra.py` checks the
same cube's longitudinal/transverse cancellation. `producer.py` verifies
the immutable parent receipt and source hashes, transfers its actual
preparation and response bounds, and emits `results.json`. The tests
include independent five-mode finite-matrix comparisons and a retained
controlled-measurement history. Those floating matrix tests are
diagnostics; the authoritative error bounds are rational inequalities
proved in `PROOF.md`.

The producer also checks every scientific field against the immutable
local `research_reference.json`. Source closure and that comparison hash
are described in `PROVENANCE.md`. No research workspace or external
absolute path is required to replay this package.
