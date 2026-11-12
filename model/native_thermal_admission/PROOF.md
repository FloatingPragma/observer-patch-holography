# Equilibrium admission for the native population kernel

The source is the unchanged closed three-port/two-face construction in
`model/native_abundance_epochs/model.py`. Its 1,296 states are
`(x0,x1,x2,k)`, with loads zero through five and `k in Z6`. The six original
seam/coin records each have probability `1/6`. Every rejected twist and
no-operation record is retained. No Metropolis acceptance or thermal bath
is substituted for the native law.

## Stationary cost drift

Write `K` for the row kernel and

`E(x,k)=5 sum_i x_i² + 2(C(k)+C(-k))`, with `C=(0,5,8,10,8,5)`.

The native acceptance inequality gives `E(t)<=E(s)` on every positive edge.
Consequently `d(s)=(KE)(s)-E(s)<=0`. Finite push/pull duality gives, for
every stationary probability `p`,

`sum_s p(s)d(s) = <pK,E>-<p,E> = 0`.

If `p` has full support and some `d(s)<0`, the same sum is strictly
negative. This contradiction excludes every full-support stationary law
on that sector. `Lean/Dynamics/StationaryLyapunov.lean` proves the generic
claim by reusing the finite expectation duality in
`Lean/Dynamics/ProtectedCharge.lean`. The positive/unital finite-record
representation in `Lean/Dynamics/PublicMarkov.lean` is the same kernel
convention. The executable supplies exact stochasticity and the drift
premises on the actual complete native records.

For each conserved total load `M=1,...,14`, the certificate gives a literal
strict-descent edge and its nonzero forward probability. Its reverse
probability is zero because the native cost cannot increase. Thus the
full-sector Gibbs family `p_T(s)=exp(-E(s)/T)/Z_M(T)` is neither stationary
nor detailed-balanced for every finite `T>0`. The conclusion is independent
of numerical exponential evaluation and applies to any strictly positive
weights. As a separate exact control, weights `2^(-E(s))` have a strictly
negative mean drift and nonzero stationarity defect in each of these
fourteen sectors. In sectors `M=0,15`, all loads are equal and the kernel is
the identity on its six flux states; every law is stationary and detailed-
balanced, without selecting a temperature.

## Every recurrent class and every stationary law

An exact positive-edge graph calculation gives 1,140 strongly connected
components, of which 46 are closed. The closed states are exactly the
native terminal set. Thirty-six classes are singletons `(q,q,q,k)`, for
`q,k in {0,...,5}`. For each of the ten loads not divisible by three there
is one three-state closed class: all balanced permutations with `k=0`.
These ten classes have transition matrix

`K_C = (1/2)I_3 + (1/6)11^T`.

Its invariant law is uniform and its two nonconstant eigenvalues are
`1/2`. Every singleton is absorbing. The source has 66 recurrent states
and 1,230 transient states.

The existing native absorption theorem gives
`Pr(not terminal at n)<=25^(-floor(n/32))` for every initial law.
If a stationary law assigned positive mass to transient states, applying
this bound to that same law at arbitrarily large `n` would force that mass
to zero, a contradiction. On each closed class the exact displayed kernel
has the unique uniform stationary measure. Therefore every stationary law
is an arbitrary normalized nonnegative mixture of the 46 class-uniform
measures, a simplex of affine dimension 45. At fixed load divisible by
three it has dimension five; at every other load its dimension is zero.
All these measures satisfy detailed balance, because each closed kernel is
symmetric. The kernel supplies no rule selecting their free mixture weights.

The cost is constant within each closed class. Restricting Gibbs weights
to the recurrent set therefore gives a stationary measure for every
positive temperature, with class weight proportional to
`|C| exp(-E_C/T)`. This is a supplied choice among the stationary simplex.
It is distinct from full support on the original sector and is not the
absorption law selected by a general native preparation. In particular a
high-cost frozen singleton is a stationary counterexample to unique
minimum-cost state selection.

## Physical boundary

This settles the specified native equilibrium-admission question. A
finite-temperature detailed-balance model with full support would require
reverse probability for every strict descent, or a changed support. Such
reverse moves and their bath energy exchange would change the declared
source. A physical thermal interpretation additionally needs an energy
conversion, a temperature readout and a clock attachment. The native
Lyapunov cost and attempted-event count supply none of those services.
No early-Universe abundance or expansion law follows from this obstruction.
