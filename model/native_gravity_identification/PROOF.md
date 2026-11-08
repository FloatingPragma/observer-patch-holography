# Entropy and modular time on the growing fair source

## Source and observation

Fix the fair connected binary source with its declared population `V`,
conserved total `q`, and native canonical reference

```
tau = I_H / D,        H = span{|x>: |x|=q},        D = binom(V,q).
```

The [owning restriction theorem](../native_capacity_restriction/PROOF.md)
derives the local state from this reference. The source specification may
be read to obtain `V` and `q`. It does not select `q(V)` or prepare this global
reference from an arbitrary input. At federation level `ell`,
`V=240*4^ell`; the integer formula is distinct from an expanded execution at
that level.

Take any fixed set `A` of `v` ports. On the global conserved-number carrier,

```
H = direct_sum_k H_(A,k) tensor H_(Ac,q-k),
d_k = binom(v,k),       e_k = binom(V-v,q-k),
p_k = d_k e_k/D,        w_k = e_k/D.
```

Impossible sectors have `e_k=0` and remain zero outcomes. The actual
regional density and natural number-preserving regional algebra are

```
R_A = direct_sum_(e_k>0) w_k I_(d_k),
M_A = direct_sum_(e_k>0) B(H_(A,k)).                             (1)
```

This algebra permits every matrix within each local number sector; it is
larger than the passive record menu. The theorem therefore applies to all
number-preserving local controls, whether or not every one is implemented.
The multiplicity tensor chart is the literal occupation split, without a
fermionic sign or a physical-spacetime identification.

The code imports the owning hypergeometric restriction instead of
replacing the local state by a chosen two-particle corner. For large `V`,
`w_k=(q)_k(V-q)_(v-k)/(V)_v` uses `v` falling-factorial factors. In particular,

```
w_(k+1)/w_k = (q-k)/(V-v-q+k+1).                               (2)
```

No global density matrix or cosmological binomial coefficient is
materialized. Every local sector is retained, including zero weights.

## A realized finite bulk and edge entropy split

The canonical source itself has the tensor form

```
tau = direct_sum_k p_k (I_(d_k)/d_k tensor I_(e_k)/e_k).
```

Consequently the tracial regional completion fixes this full reference.
There is no replacement of the source state by a different completion, as
can occur for a noncentral closed-history reference. Direct spectral
counting gives

```
S_bulk = H(p) + sum_k p_k log d_k = -sum_k p_k log w_k,
S_edge = sum_k p_k log e_k,
S(tau) = S_bulk + S_edge = log D.                              (3)
```

The central Shannon term belongs to the bulk. Equation (3) instantiates
the convention and hypotheses of the finite direct-sum identity
`OPH.EinsteinBranch.directSumEntropy_eq_bulk_add_edge` in
[Entropy.lean](../../Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean).
The attachment from native occupations to its data is proved here; the
existing Lean theorem does not itself formalize this source attachment.

For any normalized density `rho` on the same global carrier,

```
S(rho) = log D - D(rho || tau).                                (4)
```

Thus the actual canonical reference uniquely maximizes full entropy on that
finite carrier. For the sector-reweighted family
`rho_r = direct_sum_k r_k (I_dk/d_k tensor I_ek/e_k)`, equation (4) reduces to

```
S(rho_r) = log D - D(r || p).                                  (5)
```

These identities follow by inserting `log tau=-(log D)I`; nonnegativity and
the equality case of finite relative entropy give the maximum. They are
finite state statements, not preparation algorithms. The full entropy in
(4), the regional bulk entropy in (3), and a geometric generalized entropy
are different observables.

## Regional first law and its controlled remainder

Let `X=X*` have trace zero and remain on the support of `R_A`. Suppose
`m=lambda_min(R_A)>0` on that support and
`m-|t| ||X||>0`. The usual finite logarithm integral, applied to these
actual regional matrices, yields

```
S(R_A+tX)-S(R_A)
  = t Tr[X K_A] - D(R_A+tX || R_A),       K_A=-log R_A,
0 <= D(R_A+tX || R_A)
  <= t^2 ||X||_HS^2 / (2(m-|t| ||X||)).                        (6)
```

Indeed every matrix on the line segment has eigenvalues at least
`m-|t| ||X||`. The Fréchet derivative
`D log(B)[X]=integral_0^infinity (B+sI)^-1 X (B+sI)^-1 ds`
has nonnegative quadratic form bounded by `||X||_HS^2/lambda_min(B)`.
Integrating this Hessian twice proves (6), including the factor one half.
The matrix identity also covers cross-number noncommuting tangents on the
full regional support; a fixed-global-number source marginal does not prepare
those coherences. Its physical regional tangents are block diagonal and
commute with (1). No scalar Gibbs ansatz is used. The first derivative
depends only on block traces:

```
Tr[X K_A] = -sum_k log(w_k) Tr X_k.                            (7)
```

Every tangent traceless within each sector has zero first derivative and
nonpositive finite entropy change. A tangent transferring probability from
sector `k` to `k+1` has bulk slope `-log(w_(k+1)/w_k)` and edge slope of the
opposite sign. Their cancellation proves stationarity in (5); it does not
assume fixed proper volume or introduce an area coefficient.

## Exact modular flow and the native clock test

Equation (1) puts `K_A` in the center of `M_A`. Hence

```
exp(i s K_A) B exp(-i s K_A) = B       for every B in M_A.       (8)
```

The sign convention for modular time does not affect this identity.
All supported matrix units are covered. On the larger full matrix algebra,
possible inter-number coherences acquire only number-dependent phases;
every occupation observable `n_i` continues to commute with `K_A`.

The native attempt evolution is nontrivial on precisely such occupation
observables. On any internal seam `(a,b)`, one unread fair E attempt has
population map `(I+Swap_ab)/2`. Its adjoint maps the contrast
`n_a-n_b` to zero. In a source basis state with endpoints `10`, that contrast
has expectation one before the attempt and zero afterward. The tests execute
the complete original twenty-five-address native instrument, retain all raw
zero outcomes, and verify the same result. Modular flow leaves the contrast
unchanged for every parameter. Therefore no rescaling of the native attempt
counter identifies this nontrivial evolution with (8) on the retained
occupation algebra.

This excludes the direct state/algebra/clock identification in (1).
It does not exclude a different prepared reference, a different observable
algebra or a proved effective field map. An algebra isomorphism preserving
this faithful state cannot remove its centrality. Omitting every observable
on which native dynamics acts would change the declared observation.
The [Bernoulli source calculation](../native_action/screen_compatibility.py)
has the corresponding charge-only modular flow. The present result retains
the complete finite canonical number law and its capacity dependence.

## What the local preparation error controls

The [fair repair preparation](../native_fair_local_preparation/PROOF.md)
has half-diamond error at most `epsilon=2^-b` after the stated observation
marginalizes mixing-record values. For any global fixed-number input, its
local source marginal target has form

```
Sigma_A = direct_sum_k r_k I_(d_k)/d_k.                         (9)
```

The weights `r_k` are the actual input weights. Tracing the exterior of a
fixed-global-number input removes the exceptional local empty/full
coherence; it need not remove correlations of the retained archive with an
external reference. Formula (9) concerns the source marginal only. Neither
(9) nor the error estimate factors the full archive/reference state. In
particular it does not derive the canonical weights `p_k` of (1).

Let `rho_A` be the actual prepared marginal on the same number support,
`D(rho_A,Sigma_A)<=epsilon`, and `m=min_(r_k>0) r_k/d_k`. If `epsilon<m`,
normalization and Hermiticity give `||rho_A-Sigma_A||<=epsilon`.
The logarithm resolvent identity then gives

```
||log rho_A-log Sigma_A|| <= epsilon/(m-epsilon).
```

Duhamel's formula and (8) imply, for number-preserving `B`,

```
||exp(i s log rho_A) B exp(-i s log rho_A)-B||
 <= min(2, 2 |s| epsilon/(m-epsilon)) ||B||.                    (10)
```

This is a positive finite modular-error estimate with an explicit spectral
premise. A state-distance estimate alone has no uniform logarithm bound
near zero eigenvalues. In the canonical examples with fixed `q=2`, the
smallest retained positive eigenvalue is of order `V^-2`; increasing
capacity can invalidate `2^-64<m`. At fixed density and fixed twelve-port
window the floor is bounded away from zero. The executable evaluates both
cases at levels 0, 6 and 200, without asserting a global preparation at any
of those levels. If the source is exactly canonical, every unread fair
preparation map preserves that marginal and the actual error is zero;
(10) also describes independent preparations with a certified nonzero
error. No geometric modular clock follows from a small value of (10).

The allowed repair word has the original finite resource bounds. The
formal identity for arbitrary `b` is not authorization to exceed them.
Appending measurements must use the same attempt counter and resource
checks. A coherent-service pulse duration requires a separate conversion.

## Area normalization and the Einstein premises

There is a direct finite countercontrol to using a fixed graph cut as the
edge entropy in (3). Take `V=240`, `q=80`, `v=12`, and the two-sided positive
probability tangent `delta p_0=-1`, `delta p_1=1`. For
`|t|<=min(p_0,p_1)/2`, all probabilities stay positive. Equation (2) gives

```
delta S_edge = log(80/149) != 0,
delta S_bulk = -log(80/149),
delta(cut edge count) = 0.                                   (11)
```

No constant multiple of this static cut, of either sign, equals the
required edge variation. This tangent is a mathematical finite state
variation; the result does not assert an implemented controller preparing
it. Fixed combinatorics is not fixed proper volume. A responsive geometry
and its admissible physical tangent family require an independent map.

The following interface distinctions apply to the canonical dictionary.

| Einstein input | Source attachment and boundary |
| --- | --- |
| Lorentzian geometry and count/volume limit | `V` and the native cut are combinatorial data; no proper-volume or responsive Lorentzian metric is supplied. |
| Geometric modular stress | The exact modular action (8) is trivial on the number-preserving regional algebra. It cannot implement nontrivial occupation propagation under the direct clock identification. |
| Independent stress Ward identity | Number conservation alone supplies no conserved symmetric Lorentz stress. |
| Entropy, area and stationarity | Equations (3) and (6) supply the finite entropy split and modular first law. Equations (4) and (5) give finite entropy maxima. Static cut area fails (11); geometric generalized-entropy stationarity is not established. |
| Common small-ball refinement | Formula (10) requires a support floor. It is not a physical-radius or nine-direction stress-response estimate. |
| Universal coupling and positive Newton constant | Neither a metric response nor a calibrated universal area coefficient is supplied. |

These are exact finite results and a failed identification of one specified
state/algebra/area dictionary. They do not refute the conditional Einstein
theorem, the full admitted program class, or OPH as a theory. Every claimed
physical premise must be attached to one common source, state, operations,
clock and refinement domain before its conclusion is invoked.
