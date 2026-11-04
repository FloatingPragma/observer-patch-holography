# Controlled spin interactions on the bipartite source

The source's bipartite hopping removes its third-order virtual spin
correction. On every fixed finite cubic graph, this gives a fourth-order
bound on the actual low-energy generator and a complete-output time-evolution
error. It sharpens the [finite spin construction](../strong_coupling_refinement/README.md)
without changing the repair source, state space or clock. The bound is
extensive in volume. A separate [local construction](LOCAL_LIMIT.md) proves
a volume-independent interacting spin limit for dressed inputs at fixed
lattice spacing. Neither result supplies an interacting relativistic
continuum or a measured particle assignment.

## Source and parity

Use the centered source `K=T+UD-(U/2)N`, with signed hopping `t_xy=±c/2`
on a finite simple bipartite graph. At half filling `N=V`, write
`H=K+UV/2=T+UD`. Let `P` be exactly one fermion per site, `Q=I-P`,
`D0=QUDQ`, `B=QTP`, `C=QTQ`, and `tau=c|E|>=||T||`. All intermediate
empty and doubly occupied states are retained. The parent construction
proves `PTP=0`, `D0 B=UB`, and

```
S0=-(B*B)/U = (c^2/U) sum_edges(S_x.S_y-1/4).
```

Here the star is the adjoint and `S_x` is the source's internal spin.
For either sublattice A, the parity `R=(-1)^(N_A)` commutes with `D0`
and anticommutes with `T`. Its restriction to `P` is the scalar
`(-1)^|A|`. Therefore

```
B* C B = 0,
B* R_z C R_z B = 0,       R_z=(D0-z)^(-1), z<=0.       (1)
```

This is an operator identity, independent of a selected spin state or
ground-state degeneracy. `algebra.py` verifies the identity with integer
matrices on the complete half-filled dimer, chain, square and eight-site
cube. The latter has 12,870 states; its spin space has 256 states.

## Every low energy with an explicit error

Assume `U>tau`. The parent min-max and Feshbach argument gives exactly
`2^V` eigenvalues `z<=0`, with `|z|<=tau^2/(U-tau)`; the other energies
are at least `U-tau`. For `A=D0+C`, two exact resolvent expansions give

```
(A-z)^(-1) = R_z-R_z C R_z+R_z C R_z C(A-z)^(-1).
```

In `F(z)=-B*(A-z)^(-1)B` the middle term vanishes by (1). Since
`R_z B=B/(U-z)`,

```
||F(z)-S0||
 <= tau^2 |z|/[U(U+|z|)]
      +tau^4/[(U+|z|)^2(U-tau+|z|)]
 <= 2 tau^4/[U^2(U-tau)].                             (2)
```

Each ordered eigenvalue of `F(z)` is nonincreasing in `z`, so its
difference from `z` is strictly decreasing and has one root. The root
ordering agrees with the ordering of the low eigenvalues of `H`.
The variational eigenvalue bound applied to (2) pairs all of them with
the ordered eigenvalues of `S0` within the same error, including degeneracies.

## Canonical generator and all-output dynamics

Assume `U>=4tau`, and put `r=tau/U`. The parent spectral-graph proof
constructs `X:P->Q` with

```
x=||X|| <= r/sqrt(1-2r) <= 3r/2,       x^2<=2r^2,
W=(P+X)a,       a=(I+X*X)^(-1/2),       ||W-P||<=2r.
```

The range of `W` is the entire low spectral subspace. Its graph equation
is `D0 X+B=X B*X-CX`. Set `Y=X+B/U`. Then
`||Y||<=r(x+x^2)` and `B*D0^(-1)=B*/U`. Writing `M=B*X`, (1) yields

```
M-S0=(B*/U)X B*X-(B*/U)C Y,
||M-S0|| <= (tau^2/U)[x^2+r(x+x^2)] <=4 tau^4/U^3.   (3)
```

The self-adjoint canonical generator is `Hlow=W*HW=a^(-1)Ma`.
Use `||a||<=1`, `||a^(-1)||<=9/8`, and
`||a-I||,||a^(-1)-I||<=x^2/2`. Equation (3) gives

```
||Hlow-S0|| <= (9/8)4 tau^4/U^3+(17/16)(tau^2/U)x^2
             <= (53/8)tau^4/U^3 <=7tau^4/U^3.        (4)
```

Exact intertwining, two uses of `||W-P||<=2r`, and Duhamel imply

```
||exp(-itH)P-P exp(-itS0)||
 <= min(2,4tau/U+7|t|tau^4/U^3).                     (5)
```

This holds uniformly for every spin input, including entanglement with an
arbitrary spectator. The output includes every virtual-charge leakage
branch. To restore `K`, multiply the half-filled output by its known
common phase `exp(itUV/2)`. The original full-Fock charged-sector gap
`U/2-tau` is unchanged.

## Refinement boundary

At fixed graph and `U=c^2/J` with `J>0`, (5) is
`4|E|J/c+7|t||E|^4 J^3/c^2`. It controls the same interacting Heisenberg
generator, whose adjacent bond terms do not commute. In a fixed physical
box with L sites per side, `c` proportional to L gives `tau=O(L^4)`.
For `U` proportional to `L^p`, the two errors have powers `4-p` and
`16-3p`. This estimate tends to zero for `p>16/3`. The exchange then
has power `2-p` and vanishes. Fixed nonzero exchange requires `p=2`,
outside the certified domain. This is a limitation of the extensive
bound, not an exclusion of a local expansion or another continuum limit.

The [local theorem](LOCAL_LIMIT.md) removes this volume factor for a
different, explicitly dressed state/readout domain. Five exact local
integer-penalty transformations give a residual of order c^-4 and a
spin-generator correction of order c^-1 at U=c^2/J. A local propagation
estimate then gives error of order c^-1 for bounded times, uniformly in
volume. Its limit is the interacting fixed-spacing Heisenberg system.
It does not transfer the bare repair instrument or identify its prepared
states with finite-c source ground states.

The spectral-subspace method is a finite instance of the general
[Schrieffer-Wolff framework](https://arxiv.org/abs/1105.0675).
The parity cancellation and constants used here are derived above;
no unproved local or thermodynamic estimate is imported.

## Verification

Run `python -m pytest -q -W error model/bipartite_spin_limit` and
`python -m model.bipartite_spin_limit --check`. The receipt retains exact
integer identities and rational constants. Tests independently construct
Jordan-Wigner matrices, compare the complete low spectral subspace and
time evolution, and reject graph or scale inputs outside the theorem.
Numerical comparisons check finite controls; the operator proof establishes
the finite-graph result for all parameters in its stated domain.

The parent receipt's own finite checks pass and its source hashes match.
Its aggregate `all_pass` is false because its `field_symmetry` input has no
`all_pass` field. This package records that inherited status separately;
it does not reclassify or replay the parent's continuum experiments.
