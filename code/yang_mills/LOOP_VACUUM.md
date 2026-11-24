# Two-loop measurements in the full SU(2) vacuum

Use normalized Haar measure on each link of an open cubic graph, impose gauge invariance at every vertex, and retain all representations. In electric Casimir units write the Hamiltonian as
$$
H=\sum_e C_2(e)-b\sum_{p\in B} W_p,\qquad W_p=\tfrac12\operatorname{Tr} U_p .
$$
Adding the scalar $b|B|$ restores the positive magnetic convention. The physical multiplier used below is $\alpha=g^2/(2a)$, with $b=4/(g^2)^2$. This is an independently specified pure-gauge source.

## A noncommuting signal in the electric vacuum

Take two edge-disjoint squares with the same base point and holonomies $U,V$. In the electric vacuum their holonomies are independent Haar variables. For the normalized trace $W=\operatorname{Tr}(UVU^{-1}V^{-1})/2$, Schur averaging gives $\mathbb E(W\mid U)=(\operatorname{Tr}U)^2/4$ and $\mathbb E W=1/4$. This order signal exists at zero magnetic coupling.

Write $x,y$ for the real parts of unit quaternions $U,V$, and $c$ for the cosine between their axes. Each real part has normalized semicircle measure and $c$ is uniform on $[-1,1]$. Then $W=1-2(1-x^2)(1-y^2)(1-c^2)$ and $\mathbb E W^2=1/3$. As a function of each holonomy the word contains a fundamental and its conjugate, so its only spin components are zero and one. The centered conditional means have squared norms $1/16$ each. The both-spin-one component has squared norm $7/48$. In a loop of $l$ links each spin-one factor costs $2l$. Hence the connected Euclidean vacuum correlation for equal-length loops is
$$
C(t)=\tfrac18 e^{-2l\alpha t}+\tfrac7{48}e^{-4l\alpha t}.
$$
The executable independently multiplies the quaternion words and integrates them by Haar quadrature. The angle odd terms in the rectangle moment cancel in the symmetric quadrature.

## A magnetic interaction with two independently switchable faces

Choose coplanar adjacent elementary squares $p,q$ in the $3\times2\times2$ vertex box. Let $R$ be their six-link outer rectangle. With magnetic couplings on just $p$ and $q$, either single drive leaves $\langle W_R\rangle=0$ exactly. Multiplying an undriven exterior link of $R$ by the central $-I$ preserves that Hamiltonian and reverses $W_R$. Turning on both terms changes the gauge-invariant rectangle in the interacting ground state.

At the electric point a normalized elementary trace $W_p$ has norm $1/2$ in product-Haar space and electric energy 3. The rectangle has energy $9/2$. Transport its two constituent face holonomies to a common base. Integrating a private link of each face makes the resulting holonomies independent Haar elements, and their oriented product is the rectangle. Character orthogonality gives $I=\int W_R W_p W_q=1/16$. Let the electric ground be the constant function. Its first derivative under the perturbation $-W_p$ is $W_p/3$. Its mixed derivative, orthogonal to constants, is the reduced electric inverse applied to $2W_pW_q/3$ (the constant coefficient has no effect on this expectation). Differentiating the expectation of $W_R$ gives
$$
\partial_p\partial_q\langle W_R\rangle\big|_0=
2I/9+4I/[3(9/2)]=7/216 .
$$
The same coefficient multiplies $b^2$ when all eleven faces have equal coupling. The two selected faces are the only pair whose edge support modulo two gives the perimeter of $R$. Single-link center flips make every other quadratic coefficient zero. Reversing the sign of every plaquette is a unitary change of link variables. On a link in coordinate direction $i$ the center sign is $(-1)^{\sum_{j<i}x_j}$, where $x$ is its starting vertex. The rectangle retains its sign, so its expectation is even in $b$.

For $F$ driven faces use the electric-vacuum Riesz contour of radius $3/2$ and $\rho=3/(4F)$ in complex coupling. Every nonconstant spin network contains a cycle of at least four occupied links, so its energy is at least 3. The bounded magnetic perturbation has norm at most $F|b|$. The rank-one Riesz projection has norm at most two on $|b|=\rho$, by the resolvent series. Its trace pairing with $W_R$ has modulus at most two. Cauchy's estimate then bounds the remainder after the quadratic coefficient by $2(b/\rho)^4/[1-(b/\rho)^2]$. The complete interacting rectangle mean is in $[1.6010,2.3550]\times 10^{-6}$ with the two drives at $1/128$. With all eleven faces at $1/4096$ it is in $[1.6028,2.2605]\times 10^{-9}$. These are analytic enclosures in the full finite-volume vacuum.

The lowest excitation in electric units has energy at least $3-Fb$ above the interacting vacuum by min-max and the constant Haar trial state. Resolvent subtraction bounds the ground projection distance from its electric value by $D=Fb/(3/2-Fb)$. The rectangle variance is between $1/4-2D-\langle W_R\rangle^2$ and $1/4+2D$. The first energy moment of its centered vector is $6(1-\langle W_R^2\rangle)/4$ by the electric double commutator and the fundamental Casimir normalization. Jensen's inequality and the gap bound then enclose the same vacuum correlation at every positive time. [Code and exact rational enclosures](loop_vacuum_controls.py) · [Measurements](receipts/loop_vacuum_controls.json)

## What survives removing the cutoff?

For an open spatial cubic box let $E$ be the number of links and $F$ the number of faces. A cofinal schedule with $n=2^j$, $a=1/n$, $n^2+1$ vertices per direction and $b=1/(512F)$ takes physical side length to infinity and sends the spacing to zero. Its gap is at least $(3-1/512)g^2/(2a)$ at all stages, which diverges. The connected vacuum correlation of any unit-bounded observable at fixed $t>0$ is bounded by the exponential of minus this floor times $t$. Every probe with bounded variance loses its finite-positive-time spectral weight on this trajectory. This conclusion uses the physical clock and an upper bound on the correlation, and does not require selecting a favorable fitted mass.

There is also an explicit finite-volume certificate at **every finite positive coupling**, including on a trajectory that weakens toward short distances. At electric time 4 the single-link heat kernel lies between $3/4$ and $5/4$. Its departure from one is bounded by $\sum_{m\ge1}(m+1)^2e^{-m(m+2)}<1/4$ using the character expansion. The first term is below $1/5$ and successive ratios are at most $9/400$. For the full lattice the electric kernel is their tensor product. Feynman–Kac bounds the positive interacting transfer, in the positive magnetic convention, between the electric kernel times $e^{-8bF}$ and the electric kernel. Applying these inequalities to the strictly positive normalized ground state bounds its maximum/minimum ratio by $e^{8bF}(5/3)^E$. The unique positive ground state is gauge invariant by symmetry. Its squared density therefore has minimum/maximum ratio at least
$$
d=e^{-16bF}(3/5)^{2E}.
$$
The ground-state transform gives the electric gradient form weighted by this density. Comparing its variance and Dirichlet form with the product-Haar gauge-invariant Poincaré inequality gives $\Delta_{\rm phys}\ge3\alpha d$. The variance of a fundamental plaquette is at least $d/4$ by comparison with Haar and its first energy moment is at most $\alpha$. Consequently $\Delta_{\rm phys}\le 4\alpha/d$. All estimates concern the exact same ground state and full gauge-invariant spectrum.

For any compact connected simple Lie group with a fixed Casimir normalization, let $c_G$ be the least nonzero Casimir and choose heat time $s_G$ with $\sum_{\lambda\ne 1}(\dim\lambda)^2 e^{-s_G c_\lambda}\le1/4$. Such a time exists by convergence of the compact-group heat trace and connectedness. Nonempty invariant spin networks have at least four nontrivial links on these simple bipartite boxes. The identical argument with real normalized Wilson characters gives the finite physical gap floor $4c_G\alpha e^{-4s_G bF}(3/5)^{2E}$. This group-uniform argument retains group-dependent constants and volume dependence. A smooth nonconstant local trace retains positive variance at each finite stage.

For the illustrative bare trajectory $g^2=4/(1+j)$ on the same growing boxes, the floor is positive at every finite stage and tends to zero. Its ceiling diverges. The comparison loses control through the whole-box factors $F$ and $E$. It resolves neither the limiting spectrum nor a physical running coupling. Local, refinement-uniform estimates on a dense set of renormalized observables would be needed to use the correlation gap theorem on that branch. [Yang–Mills gap paper](../../extra/yang_mills_gap_clay_problem.tex)

Run from the RER root with a fresh output path.

```sh
.venv/bin/python -B code/yang_mills/loop_vacuum_controls.py --output /tmp/ym-loops.json
```
