# Two plaquettes with a shared gauge link

The source is the positive Hamiltonian

`H = alpha sum_e C_e^2 + beta (2-W1-W2)`,

on two adjacent square plaquettes, with six vertices, seven links and the
full Haar space on every link. Here `alpha=g²/(2a)`, `beta=2/(g²a)` and
`Wi=Tr(Ui)/2`. Both coefficients are declared positive inputs. Each loop
boundary has zero divergence at every vertex. The two faces share one link.
All local SU(2) Gauss transformations are imposed.

## Complete physical basis

At a bivalent vertex the two incident spins agree and its invariant tensor
is unique. The three exterior links of each face therefore carry one spin,
`j1` or `j2`. The shared link carries `J`. At either trivalent vertex the
invariant tensor exists precisely when `|j1-j2|<=J<=j1+j2` and
`j1+j2+J` is integral, and it is unique. Peter–Weyl completeness gives one
orthonormal physical state for every such triple. No upper spin is part of
the source definition.

Equivalently, choose the same base vertex and orient both loops through the
shared link in the same direction. Reversing a face orientation leaves its
real SU(2) Wilson trace unchanged. Tree gauge fixing gives
`L²(SU(2)²,Haar)^{Ad SU(2)}`. For `di=2ji+1` and `dJ=2J+1`, its normalized
basis functions are

`Psi_(j1,j2,J)(U,V) = sqrt(d1 d2/dJ) Tr[P_J (D_j1(U) tensor D_j2(V))]`.

Here `P_J` is the orthogonal Clebsch–Gordan projector. Matrix-coefficient
orthogonality gives norm squared `Tr(P_J)/(d1 d2)=dJ/(d1 d2)` before the
displayed normalization, and orthogonality of different triples. The three
outer links of each loop contribute `3 ji(ji+1)`. A shared-link derivative
rotates both based holonomies on the left; its Casimir is `J(J+1)`.
Consequently the exact electric diagonal is

`3 j1(j1+1) + 3 j2(j2+1) + J(J+1)`.

The `J` dependence is retained. Replacing it by `j1(j1+1)+j2(j2+1)` gives the
separate independent-loop control, whose ground is a product of the two
one-loop ground functions.

## Wilson multiplication

Multiplication by `W1` leaves `j2` fixed and permits
`j1'=j1+-1/2`, `J'=J+-1/2`. Its matrix element is

`(1/2) sqrt(d1 d1' dJ dJ') { j1 j2 J ; J' 1/2 j1' }²`.

This identity follows by expanding the two projectors in Clebsch–Gordan
coefficients, multiplying a matrix coefficient by the fundamental trace,
and applying Haar orthogonality. The remaining recoupling contraction on
each matrix index is the same 6j coefficient; their product is its square.
The normalization factors give the displayed dimensions. The 6j convention
is the finite 3j sum in [NIST DLMF §34.4](https://dlmf.nist.gov/34.4).
`W2` interchanges the two loop labels. An independent program expands all
Clebsch–Gordan indices and integrates them by matrix-coefficient
orthogonality, without calling a 6j symbol. It checks the emitted exact
matrix elements. Complete successor sums also check the commutation of
the two multiplication operators, including successors outside a displayed
finite projection.

Each Wilson multiplier is real with norm at most one. Thus the magnetic
potential is positive and bounded. The elliptic link Casimir sum on the
compact group has compact resolvent; its Gauss restriction does also. A
bounded potential preserves self-adjointness on the electric domain and
compact resolvent. The positivity-improving heat kernel supplies a unique
positive gauge-invariant ground. These are finite-graph statements on
infinite-dimensional Haar spaces.

## Complete complement

The diagnostic projection retains `2j1,2j2<=N` and every admissible `J`.
An omitted state has at least one loop spin at least `(N+1)/2`. Positivity
of the other electric terms and of the magnetic potential gives

`Q H Q >= d_N Q`, with `d_N=3 alpha (N+1)(N+3)/4`.

Only the magnetic term crosses the projection, and its off-diagonal block
has norm at most `2 beta`. The complement is infinite and is bounded as a
whole. A finite coefficient matrix never replaces the source Hamiltonian.

The graph is a finite interacting gauge regulator. It supplies no
three-dimensional thermodynamic limit, four-dimensional field
reconstruction, physical glueball identification or GeV scale.
