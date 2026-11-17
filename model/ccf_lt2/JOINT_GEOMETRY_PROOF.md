# Complete local-capacity tensor action

Each source energy oracle returns the full column `h_v|s>` on the native
fermion and flux carrier. Each vertex carries independent Phi and capacity
factors. With the scalar references specified in [the model](../SPEC.md),

```
k_v = g_v + N_v(h_v-e_v I) + e_v I
H = sum_v [A_v + M_v(k_v-F I) + F I].
```

`local_capacity.geometry_shard_column` evaluates the complete local Phi path
kinetic, onsite regulator/boundary energy, and half of each incident edge
gradient. The sum over vertices is exactly `H_Phi`.
`local_energy_column` applies that geometry column, the entire source shard
multiplied by its positive diagonal `N_v`, and `(1-N_v)e_v` on the input.
`apply_local_capacity_column` applies `A_v`, every column of `M_v` to the
entire `k_v` column, `-F M_v` and `F I` for each vertex. It combines equal
output keys and removes only exact zeros. `apply_native_column` supplies
all actual source-shard callbacks; `source.apply_total_column` delegates
this same operation.

Equality on the full tensor basis establishes the operator identity. The
256-dimensional finite fixture checks every column against an independently
assembled Kronecker matrix, retaining all tensor factors. Native triangle
columns exercise the complete 288-mode source without allocating its state
vector. Their complete support, Gauss identities and first/second phase
response are checked. These columns do not constitute arbitrary native
state propagation or an interacting spectrum calculation.

All auxiliary factors commute with source Gauss operators. Each `h_v` is
Hermitian and preserves Gauss, and `N_v`, `M_v` commute with their source
factors. Thus `H` is Hermitian and preserves every Gauss eigenspace, including
mismatch complements. Finite free evolution is unitary and conserves the
total energy. Instrument work is a separate complete-channel energy change.

## Disconnected components

Each retained source word and field face is supported in its graph component.
The Phi-gradient terms have the same property, while the capacity factors
are individual vertex factors. Physical disjoint fermion operators are even.
Consequently terms in distinct components commute, `H` is a sum of component
Hamiltonians, and its exponential factors exactly. A trace-preserving local
operation on one component cannot alter the other component's reduced state.
The checker verifies the component sum on every basis column of a complete
256-dimensional fixture, including unequal local references and nonzero
common capacity reference. A separate complete 1024-dimensional diagnostic
retains both default 16-state capacity registers and all source outcomes.
Its receiver probabilities are independent of the disconnected sender label.

Bounded degree, face size, finite cutoffs and uniformly bounded coefficients
give bounded finite interaction range. Uniform physical velocity under
refinement requires additional scaling, reference and metric premises.

## Reference convention

Replacing `h_v` by `h_v+c_v I` and `e_v` by `e_v+c_v` changes `k_v` by `c_v I`.
For uniform `c_v=c` and common `F` shifted by `c`, each vertex summand changes
by `c I` and the total changes by `|V|c I`. At fixed `F`, the change is
`sum_v c_v M_v`. A single common `F` cannot generally absorb nonuniform
shifts. Exact column checks exercise both identities; no independent `F_v`
parameter is present in this law. Scalar energy and calibration inputs are
retained as supplied coefficients, with no identification as unfitted
physical parameters.
