# Source color normalization and finite coupling readout

The actual port-admitted color currents have fundamental normalization
`Tr(T_a T_b)=delta_ab/2` and the selected fifteen-state matter Casimir
`C3=0 or 4/3`. On the configured native gauge link, complete endpoint
records determine its hopping and electric coefficients, with rational
finite-time error bounds. The [proof](PROOF.md) retains all five matter
blocks and all sixty link states.

The actual source defaults are `J=1`, `epsilon=1`, hence
`epsilon/abs(J)=1` in its specified model clock. Reading these values from
the simulator is permitted; the complete detector witness recovers them
independently. They do not carry an identified QCD scheme or scale.

The source's uniform nonvacuum electric energy cannot equal any combination
`a C3+b C2+c Y^2`. Its Gauss constraints admit different electric-to-hopping
ratios. These exact results distinguish a measurable coefficient of the
configured finite model from a source-selected QCD coupling. The inherited
running calculation uses a stipulated D10 input and external threshold
masses. No physical `alpha_s` value is predicted here.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.source_strong_coupling --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.source_strong_coupling --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/source_strong_coupling
```

The [receipt](results.json) binds the actual source and complete rational
readout enclosures. It records mathematical identities and inherited
architecture separately. No empirical strong-coupling measurement enters
its producer.
