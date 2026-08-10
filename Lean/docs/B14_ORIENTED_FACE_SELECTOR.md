# B14 oriented-face bracket discriminator

## Exact incidence result

`Screen/OrientedFaceBracketSelector.lean` uses the twenty oriented faces in
the pinned Reynolds port frame. `face_bracket_eq_sixty_r13` kernel-checks the
complete tensor-coordinate identity

\[
B_{\mathrm{face}}=60R_{13}.
\]

The resulting alternating bracket is not a Lie bracket:
`jacobi_failure_witness` gives one exact Jacobiator coefficient equal to
`-1`. The independently replayed Python certificate checks the full 2640
independent-coordinate census: 240 nonzero coefficients, split into 120 at
`+1` and 120 at `-1`.

## Conditional three-norm discriminator

On the same 792 upper-triangular structure-constant coordinates, the exact
certificate compares the face bracket with the complete classified compact
`P`, `F`, and `G` families under three supplied coordinate-edit rules:

| edit | `G` | `F` | `P` |
|---|---:|---:|---:|
| total absolute (`L1`) | `30(sqrt(5)-1)` | `60` | `60` |
| Hilbert--Schmidt squared (`L2^2`) | `(615-123sqrt(5))/22` | `(615+123sqrt(5))/22` | `45` |
| worst coordinate (`Linfinity`) | `(5-sqrt(5))/10` | `sqrt(5)/5` | `1/2` |

The `F` value in the `L1` row is an infimum over its compact stratum, not an
attained minimum. Exact primal-dual Python replay proves the 792-coordinate
optimizations and their feasibility conditions. Lean checks the serialized
radical values and strict family order; `three_norm_unique_nearest_G` proves
that `G` wins every one of the three encoded comparisons.

## Open boundary

The face incidence is pinned input, and each norm is a basis-dependent added
premise. Neither a norm, a minimum-distance Jacobi repair rule, nor a physical
bracket producer is selected by the source. The comparison ranges only over
the compact locus classified under the three named compact-Lie inputs; it is
not a classification of the full Jacobi variety or its noncompact
components. Ordered source tomography, source bracket selection, and
same-current holonomy remain open under issues #705 and #697.
