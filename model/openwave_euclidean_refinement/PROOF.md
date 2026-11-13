# A common Hamiltonian, Euclidean response and physical clock

Use the full-Haar shared-link source and the complete complement and ground
certificates in [the coupled spectrum](../openwave_coupled_haar_spectrum/PROOF.md).
Let `Omega` be its unique normalized ground, `E0` its energy, `W=W1`, and
`v=(W-<W>)Omega`. The observable is the complete source quantity

`C(t)=<v,exp[-t(H-E0)]v>`, for `t>=0`.

The certified variance proves `v` is nonzero. Its finite energy upper bound
comes from the same native double commutator. The source spectral theorem
therefore gives a positive spectral measure, temporal reflection positivity,
and `0<=C(t)<=C(0) exp(-gap t)`. Positivity of this composite observable
alone is not a proof of nonfree four-dimensional continuum fields.

## Finite computation with the full tail retained

Let `u` be the certified normalized native trial and `P_N` the diagnostic
projection. The represented centered vector uses the exact dyadic rounding
of `P_N W P_N`. Its difference from `(W-<u,Wu>)u` is at most the trial's
boundary norm plus twice the Wilson-rounding error. Both centered vectors
have norm at most one before this rounding, so the quadratic-form error is
bounded by `(2+delta)delta` for this difference `delta`.

The exact full ground projector is at distance at most `epsilon` from the
trial projector. For the positive contraction `B=exp[-t(H-E0)]`, the
uncentered quadratic form changes by at most `2 epsilon` and the squared
mean by at most `4 epsilon`. The approximate trial's centered expression
differs from its uncentered expression minus squared mean by at most
`3 epsilon`: `(B-I)u` has norm at most `epsilon`. Thus the program's
`12 epsilon` ground-comparison budget is conservative.

Write the projected Hamiltonian as its nonnegative diagonal part and the
bounded magnetic hopping part, whose norm is at most `2 beta`. A hopping
insertion changes one loop label by one. Starting with labels at most `m`,
at least `N-m` insertions are required to reach the boundary at `N`.
The Dyson expansion and contractivity of the diagonal semigroups give

`||P_boundary exp(-s P_N H P_N) P_m||
 <= exp(2 beta s) (2 beta s)^(N-m)/(N-m)!`.

Split the represented vector into `P_m` and its explicitly computed tail.
Duhamel's formula, `||Q_N H P_N||<=2 beta`, and contractivity of the full
and finite positive Hamiltonians bound the semigroup quadratic-form error
by

`||v_N|| 2 beta t [||Q_m v_N||
 + exp(2 beta t)(2 beta t)^(N-m)/(N-m)! ||v_N||]`.

Every `m=0,...,N-1` is valid. The code chooses the smallest of these proved
upper bounds. The infinite representation complement is therefore included
before a finite spectral sum is interpreted as a full-source result.

The complete dyadic frame is polar-corrected only inside the existing exact
Gram/residual proof. Its certified operator-norm generator error `delta_H`
and Gram error `eta` bound the finite spectral-sum error by
`||v_N||²(t delta_H+3 eta)`. All stored energy coordinates are positive and
every spectral weight, including the small vacuum weight, is included.
The norm bound follows from Duhamel for two positive generators and the
Gram bound for the frame-to-unitary change. The same frame theorem supplies
the exact `delta_H`; no numerical eigensolver output is trusted as an error.

The full-to-finite comparison uses `exp(t E0)`. The reported finite sum
uses `exp(t mu)`, with `mu` the trial energy coordinate. The ground interval
bounds the normalization error by
`t |E0-mu| exp[t max(E0,mu)] ||v_N||²`. The program propagates the complete
interval and all previous errors through this common amplification.
Taylor bounds and outward dyadic rounding are exact rational inequalities.
The small exponential tail uses the rational certificate `exp(64)>2^92`.

At zero time the result is intersected with the independently certified
variance. At positive time it is intersected with `[0,variance_upper]`.
The three diagnostic resolutions must overlap at every reported time.
The largest resolution must retain a positive lower bound at `t=1/4`.
For the finite spectral witness, reflection forms are sums of nonnegative
squares, with the signed test coefficients `(1,-2,1)` all retained.

The transfer and propagation statements used here are those of the
[source transfer proof](../openwave_su2_transfer/PROOF.md) and
[source propagation proof](../openwave_su2_propagation/PROOF.md): imaginary
time generates Euclidean records. It does not turn postselected cooling
into real-time observer dynamics.

## Source-clock obstruction and continuum target

At fixed coupling the declared Hamiltonian satisfies `H(a)=H(1)/a`.
The ground and Wilson variance are unchanged, while the gap is divided
by `a`. The program evaluates this exact family at six spacings, retaining
the same source clock. Consequently its correlation at any fixed `t>0`
tends to zero as `a` tends to zero, while `C(0)` has a strictly positive
fixed lower bound. This fixed-graph, fixed-coupling sequence cannot have a
nontrivial strongly continuous field limit with this observable and clock.
Changing the clock defines a different physical attachment; it does not
repair this specific sequence.

A viable four-dimensional family instead requires spatial refinements at
fixed physical volume followed by the required volume limit, controlled
coupling flow, common local field limits with regularity and Euclidean
invariance, noncollapse, the pure Yang--Mills short-distance identification,
and transfer of a positive physical gap. A finite nonzero physical gap
requires the corresponding dimensionless gap to scale with `a` in the
declared clock convention. The present result does not prove such a family
exists, nor exclude every changing-coupling family.

The [official Clay target](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)
is nontrivial quantum Yang--Mills on four-dimensional Euclidean space for
every compact simple gauge group, with appropriate axiomatic properties
and a positive mass gap. The existing
[acceptance contract](../particle_yang_mills_gap_contract/PROOF.md) retains
these distinct obligations. An exactly conformal flat vacuum with
vacuum-preserving dilations and a nonzero finite-energy excitation cannot
have a positive spectral gap: dilations scale that excitation below any
positive bound. The elementary energy implication and the inverse-spacing
obstruction are checked in
[`DilationGapObstruction.lean`](../../Lean/ObserverPatchHolography/DilationGapObstruction.lean).
They do not apply to a compact spatial cylinder or a vacuum that breaks
dilation covariance, and do not construct a quantum field theory.
