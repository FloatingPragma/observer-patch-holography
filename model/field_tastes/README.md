# Interactions, all Weyl pockets, and finite thermal positivity

This CP-06 / [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142)
and CP-13 / [#1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137)
packet uses the **same supplied** Weyl-contact action as
[`interacting_field`](../interacting_field/README.md) and
[`bound_pair`](../bound_pair/README.md):

```text
H = T + U P,  T=dGamma(h),  h(k)=-c sum_j sigma_j sin(k_j),
c=1/8, U=1, N=2, total momentum K=0.
```

In a momentum fiber, `P=|d><d|` is the on-site singlet contact projector.
Spatial lattice, periodic boundary condition, CAR algebra, coupling constants,
input states and momentum measurements are supplied. Model Hamiltonian time
is used throughout. No value here is a physical mass, measured cross section
or calibrated laboratory rate. A periodic cube is not identified with the
open cube at finite volume.

## Why dropping the other pockets changes the interacting theory

The open-cube field packet retains all eight Weyl nodes and their balanced
chiralities. That free-symbol count alone does not prove whether a selected
interacting observable sector closes. This packet supplies an operator-level
answer for the declared momentum-pocket projection.

Let `Q` select any relative-momentum set closed under `k -> K-k`, so it
preserves the physical fermion exchange constraint. It acts on every spin
component at those momenta. Because `T` is diagonal in relative momentum,
`[T,Q]=0`. Put `w=<d|Q|d>`. For `0<w<1`, define

```text
a = Qd/sqrt(w),       b = (I-Q)d/sqrt(1-w),
d = sqrt(w) a + sqrt(1-w) b.
```

The commutator `[P,Q]` vanishes off `span{a,b}` and on that space is
`sqrt(w(1-w)) [[0,-1],[1,0]]`. Therefore

```text
||[H,Q]|| = U sqrt(w(1-w)),
||(I-Q) H a||² = U² w(1-w).
```

Both statements hold on the physical antisymmetric fiber: `a` and `b`
belong to it. Thus a proper pocket with nonzero contact weight is not a
reducing subspace of this interacting Hamiltonian. The statement is stronger
than a free node count and narrower than a no-go for chiral field theories.
It does not exclude other regulators, dressed invariant sectors, interactions
with a derived selection rule, or a controlled limiting decoupling map.
Dark vectors whose contact overlap vanishes are not by themselves a complete
interacting single-pocket theory.

The retained torus has side six. Classify each component of `k` by the sign
of `cos(k_j)`. There are no `cos(k_j)=0` grid points, each of eight cells has
27 momenta, and the cells are exchange symmetric. Each cell contains one
Weyl node `pi*(b_x,b_y,b_z)` of chirality `-(-1)^(b_x+b_y+b_z)`.
Every pocket has `w=1/8`, giving exact squared coupling `7/64`. A union of
two pockets has `w=1/4` and squared coupling `3/16`. The union of all pockets
and the `U=0` deletion provide the exact zero-commutator controls.

Prepare the normalized singlet contact vector restricted to the first pocket
and evolve under the complete action. All eight unnormalized measurement
outputs are retained. At model time two the total probability in the seven
other pockets is approximately `0.267060`; all seven have positive weight.
These are supplied momentum measurements of the finite field state, not a
compiled native field detector. The independent test diagonalizes the
complete **424-dimensional physical antisymmetric fiber**, embeds its answer
back into all 864 spin/momentum coordinates and compares the resulting
evolution. No positive spectral or spin component is omitted.

A strictly positive early-time leakage bound does not depend on those
floating results. The action satisfies `||H|| <= U+2sqrt(3)c < 3/2`.
The integral remainder for self-adjoint evolution gives

```text
||exp(-itH)-I+itH|| <= t² ||H||²/2.
```

Using `5/16 < sqrt(7)/8 < 1/3`, the outside probability at `t=1/32` is
between `5041/67108864` and `80089/603979776`. The replay yields about
`0.000106799`, strictly between these bounds. The rational bounds concern
ideal Hamiltonian evolution; numerical contraction and hardware errors are
not absorbed into them. Renormalizing the chosen-pocket outcome would remove
this observable failure and is explicitly rejected by the complete-output
tests.

## A positive finite thermal premise on that same action

On the complete physical `N=2,K=0` matrix algebra, choose the faithful density
`tau=exp(-beta H)/Z`, with the declared value `beta=log(25/8)`. The choice is
compatible as a numerical ratio with the separate two-level clock corner;
it is **not** preparation of this full thermal state by native repair.
The full fiber includes all nodes, contact and scattering sectors.

For arbitrary operators `A_i` on this algebra and `0<=t_i<=beta/2`, define

```text
G_ij = Tr[ exp(-(beta-t_i-t_j)H) A_i* exp(-(t_i+t_j)H) A_j ] / Z.
```

In an energy basis with `p_n=exp(-beta E_n)/Z`, this is the Gram matrix
of the vectors

```text
F_i(m,n)=sqrt(p_n) exp(t_i(E_n-E_m)) (A_i)_(m,n).
```

Hence `G` is positive semidefinite, including its complex off-diagonal
entries. The numerical packet independently checks the direct trace formula
for the contact projector, a momentum-pocket projector, and `i[H,P]`, at
three Euclidean model times. These are even, number-preserving observables;
no odd full-Fock field is represented by an unjustified N=2 compression.

The same spectral sum proves the finite KMS identity

```text
C_A(t) = Tr[tau exp(tH) A exp(-tH) A*]
       = C_(A*)(beta-t).
```

The test uses the non-Hermitian product of contact and pocket projectors,
which do not commute. A tracial density substituted at the same nonzero beta
fails the identity, exposing the reference premise. The finite positivity
and KMS facts are ordinary matrix identities, attached here to this action;
they are not general thermal theorems.

This does not establish spatial locality of the momentum observables,
reflection positivity of a continuum field measure, a physical vacuum,
Lorentz covariance, all Osterwalder--Schrader reconstruction hypotheses,
Higgs/current matching, or decoupling of the mirror tastes. It supplies a
concrete positive finite thermal premise and an exact failed taste-removal
map on the same interaction. CP-06 is open.

The mobile-memory PR #1153 supplies active native code transport
under its declared fair-coin law. The triplet-attachment PR #1128
studies a different operational dictionary and scoped decoder obstructions.
Neither provides a chiral-regulator or full-field clock map for this action;
this packet uses no code from either PR.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.field_tastes --check model/field_tastes/results.json
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/field_tastes
```
