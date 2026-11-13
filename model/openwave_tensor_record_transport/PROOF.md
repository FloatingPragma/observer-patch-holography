# Finite tensor records and the Einstein decoder

The input is a jointly available classical record table on the existing
staggered lattice. This construction supplies the missing co-location map
for the finite null-tensor theorem. It does not supply a quantum instrument
that measures a physical stress tensor or curvature field.

Coordinates are in half-spacing units. Diagonal entries lie on the
scalar lattice. For each even site x, retain C_aa(x)=B_aa(x) and define

```
C_ab(x) = sum_{i,j=+-1} B_ab(x+i e_a+j e_b)/4,  a!=b.
```

All outside halos and zero records are included. Each off-diagonal entry
uses the four actual surrounding records. Define
`d_i(ab)=sup_y |B_ab(y+2e_i)-B_ab(y)|` on its own sublattice. For the corner
`y=x+e_a+e_b`, the four distances require respectively zero, one a-step,
one b-step, or both. The triangle inequality gives

```
|B_ab(y)-C_ab(x)| <= (d_a(ab)+d_b(ab))/2.
```

The program obtains these variations from the complete finite table. No
continuum derivative bound is assumed. The coefficient is attained by
affine data locally. A periodic checkerboard has zero corner average and
nonzero original entries: averaging is not invertible, and its variation
budget cannot be dropped.
`TensorRecordTransport.lean` proves the four-corner inequality and the
addition of its error to an existing decoder remainder.

Apply the existing nine-null inverse to C. If each null contraction has
absolute size at most epsilon, including its observation error, the existing
sharp theorem gives every entry of
`C-lambda eta` at most `(2+sqrt(3))*epsilon`, where
`lambda=tr(eta*C)/4`. The diagonal entries of B equal those of C; the metric
has zero off-diagonal entries. Hence the original staggered table satisfies

```
B=lambda eta+R,
|R_ab| <= delta=(2+sqrt(3))*epsilon+tau,
tau=max_{a!=b}(d_a(ab)+d_b(ab))/2.
```

Its existing staggered Bianchi/error theorem applies directly. With
divergence bound w, the actual nearest-neighbour scalar edge obeys
`|lambda(x+2e_j)-lambda(x)| <= a*w+8*delta`. Sum these budgets on actual
graph paths. The nine-null inverse and its error inequalities are reused
from `EinsteinBranch/Tensor.lean` and `EinsteinResidualBound.lean`; no additional
matrix inversion or repair-order premise is introduced.

For bounded finite input words, rational averages and the fixed quadratic
number field Q(sqrt(3)) admit finite exact encodings. The complete classical
record map has a reversible realization
`|x,z> -> |x,z XOR code(f(x))>`, retaining all input labels. This is a
permutation on the entire finite record carrier, including unused output
labels. It copies orthogonal classical labels, not unknown quantum states.
Initialized output memory is an explicit resource. No physical acquisition
work or free erasure theorem follows. The collecting observer is a common
causal successor of its inputs; the map does not require simultaneous
access to spacelike or incompatible quantum observables.

The exact examples use nonstatic compact metric variations and all their
curvature halos. Their divergences vanish identically. The arbitrary-table
decoder conclusion follows from the existing exact inverse identity, with
six direct numerical-algebraic checks per example. Nonlinear curvature,
physical stress/metric identification, continuum control and the value of
the cosmological metric term remain separate obligations.
