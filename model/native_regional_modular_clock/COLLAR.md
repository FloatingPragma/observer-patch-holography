# Interior compatibility of the regional modular clocks

Fixed-cut modular covariance fails by the exact 4/3 example in PROOF.md.
A growing collar gives a different positive result. All operators below use
the same binary source graph and the same declared density p. The result is
conditional on the common rooted port graph; it does not resolve the
all-level numerical routing boundary or prepare a global transfer state.

Write J=1/[1-p(1-p)], mu=1-2p. For a finite vertex region A let
D_A=G_A+mu N_boundary,A, with the actual crossing degrees. Let O have support
S inside A, s=|S|, and suppose the graph distance from S to every inside
endpoint of an A-crossing seam is at least r. Let the maximum graph degree
be Delta and the crossing count be m_C.

## Finite path bound

For nearest-neighbour exchange with coefficient J, scalar identity terms
may be removed, so each interaction is -J Swap_e and has norm J. Define
T_r(x)=sum_(n>=r) x^n/n!. A finite commutator estimate is

```
||[alpha_(JG_A)^t(O),B]||
 <= 2s ||O||||B|| T_r(4 Delta J |t|)                 (1)
```

when B's support is at graph distance at least r from S. This estimate is
independent of the total volume and permits B on a crossing edge.

For completeness, the finite interaction-picture commutator recursion is
C_B(X,t)<=C_B(X,0)+2 sum_(e meeting X) J integral_0^|t| C_B(e,u) du,
where C_B is the supremum over norm-one operators supported in X. It follows
by Jacobi's identity, moving the norm-preserving commutator term into its
unitary interaction picture, and integrating the remaining inhomogeneous
term. Initial values are bounded by2||B|| when the supports intersect and
zero otherwise. Iteration gives an ordered chain of n overlapping edges.
There are at most s Delta choices for its first edge and2Delta choices for
each later edge, and a chain reaching B has at least r edges. Its simplex
integral is |t|^n/n!. These counts yield (1), with a harmless extra factor two
also covering r=0. The recursion is the finite-matrix argument of equations
(2.9)-(2.12) in [Nachtergaele, Ogata and Sims, *Propagation of Correlations
in Quantum Lattice Systems*](https://arxiv.org/abs/math-ph/0603064); the native
swap specialization and explicit constants are given above.

Compare both the larger-region exchange generator JG_B, B containing A,
and the regional generator JD_A with the common free-cut generator JG_A.
Duhamel's identity requires only crossing seam operators and the boundary
occupation potentials; all exterior terms commute with the free-cut
observable. Their total operator-norm coefficients are at most
J m_C(1+|mu|). Integrating (1) therefore gives

```
||alpha_(JD_A)^t(O)-alpha_(JG_B)^t(O)||
 <= 2s J m_C(1+|mu|)|t| ||O|| T_r(4 Delta J |t|).    (2)
```

Cap the right side at2||O||. No state assumption, strict support cone or
relativistic interpretation is used. This is an approximate propagation
bound in the declared exchange time.

## Join to the actual regional density

Let E_A be the maximum generator error from equations(3)-(4) of PROOF.md
across its supported local occupation sectors. For this modular join O must
preserve local occupation number; otherwise the block scalar terms in the
modular generator cannot be discarded. Combining those equations
with (2) gives, in the same time orientation,

```
regional modular vs larger exchange flow error
 <= 2|t| E_A ||O|| + right side of(2).               (3)
```

It suffices to compute E_A without enumerating even the local sectors. Put
v=|A|, b=V-v, L=2m_A+m_C and suppose L<M. Define

```
d_r=(V/b)|q/V-p|+v/b,
d_a=d_r+(m_A+3m_C)/(4M)+1/[4(b-1)].
```

Since a_k>=1/2 and a_p>=3/4, an all-sector upper bound is

```
E_A <= L^2/[M(1-L/M)] + 4m_C d_r + (8/3)L d_a.       (4)
```

This requires only source counts and the authenticated finite region.
`collar.py` evaluates (2)-(4) using exact rational tail enclosures.
An admitted local observation has its separate preparation, detector,
coherent-service and archive errors; none is silently set to zero.

For two regions containing S, compare each to their finite union or a
common larger window and add their bounds(3). Thus compatible interior
modular actions follow on a rooted graph along any sequence for which
E_A tends to zero and m_C T_r(4Delta J|t|) tends to zero, uniformly on the
chosen bounded time interval. Bounded degree gives at most exponential
boundary growth, while the factorial tail beats that growth for fixed time.
The capacity and collar must be chosen jointly so L^2/M and the density
terms vanish. For q/V=p exactly, or q=floor(pV), v^2/M ->0 is sufficient for the
finite source family with M proportional to V. For a merely convergent
density sequence, also require v|q/V-p| ->0; convergence without a rate
does not control a growing regional algebra. The fixed-cut counterexample is
unaffected because its r stays bounded.

The literal rooted-source attachment supplies finite certified windows at
levels up to14. These finite certificates give executable error bounds;
they are not evidence that floating port routing converges at every level.
An all-level compatible routing or a specifically justified rooted
subsequence remains a premise of an infinite net. The resulting local net,
where that premise holds, is a graph exchange dynamics. The celestial
Lorentz equivalence, Lorentz vacuum, physical count/volume dictionary and
responsive Einstein geometry require their own actual attachments.
