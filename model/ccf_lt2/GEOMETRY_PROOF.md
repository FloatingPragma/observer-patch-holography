# Per-vertex volume, capacity and memory

`GeometryConfig` defines one finite capacity register at each source vertex.
The default basis is volume `(1,3/2)`, capacity `(1,2)`, and the four hardware
memory words. Its dimension is 16. With normalized volume `x` ranging from
minus one to plus one, the positive lapse is `M=I+x/4`, and

```
A = (1/2) L_volume + (1/8) L_capacity
    + diag[(1/32) volume/capacity
           +(1/16)(volume-1)^2 +(1/16)(capacity-1)^2].
```

Identities act on all unmentioned factors, including memory. The configurable
member retains ordered positive volumes, a positive integer maximum capacity,
all `2^capacity_max` memory states, nonnegative finite hopping/potential
coefficients and lapse coupling in `[0,1)`. `joint_geometry.auxiliary_columns`
returns every exact rational snapshot entry of `A` and `M`.

`build_geometry_extension` is a single-vertex matrix primitive for
`A + M(h-F) + F`. It validates the whole supplied finite `h`, optional component
sum, and its Gauss commutators. The official multi-vertex law applies this
primitive independently to each `k_v` as specified in the
[complete tensor proof](JOINT_GEOMETRY_PROOF.md). The standalone witness is a
complete three-CAR diagnostic of this one-register primitive.

## Complete capacity instrument

At capacity `k`, an outcome records the volume, capacity and hardware-memory
residue modulo `2^k`. Its projector sums over every inactive bit. These
projectors are orthogonal and sum to identity on the entire capacity register.
Their tensor products over vertices give the complete multi-vertex capacity
instrument. Reducing the exposed alphabet never erases a hardware bit.
A tensor product with a complete source instrument retains every source
outcome, including zero-probability and complementary sectors.

For a supplied normalized full finite state the matrix primitive evolves
with its Hermitian Hamiltonian. The instrument routine retains branch
probabilities, conditional full states, and measurement work
`Tr H(sum_a P_a rho P_a-rho)`. Energy conservation applies to autonomous
propagation; dephasing generally changes energy. The tested apparatus is an
explicit supplied instrument, without an identification as a derived
physical detector or self-starting clock.
