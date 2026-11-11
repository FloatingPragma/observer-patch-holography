# Complete-record conjugation and asymmetry bounds

## Native local identity

The source is the integer-load law and class register of
[`defect_integer_lift.model`](../defect_integer_lift/model.py). Loads lie in
`{0,...,5}`. On each addressed seam the repair balances their sum using a
binary coin. For the local state `s=(x_a,x_b,k_+,k_-)`, define

```
C(s) = (5-x_a, 5-x_b, -k_+ mod 6, -k_- mod 6).
```

The complementary coin balances the complementary sum: the repaired loads
are `5-x'_a,5-x'_b` and its current is `-j`. The contrast changes from `d` to
`-d`, so `j(d-j)` is invariant. The cost table `(0,5,8,10,8,5)` is invariant
under class negation. The acceptance decision is consequently identical and
the proposed classes conjugate. This proves `F(C(s),1-c)=C(F(s,c))`, including
zero currents and rejected twists. The executable checks all `6^4 * 2`
local branches against the owning implementation, without dropping records.

The integer link and face lifts change sign because the increment changes
sign. The signed curl commutes with that negation. Any fixed seam schedule
therefore has a paired path with conjugated loads, classes, lifts and complete
records. An independent state-insensitive random schedule also has this
property; an adaptive scheduler must separately satisfy equivariance.

## Ensemble consequence and formal theorem reuse

For unbiased coins and a conjugation-symmetric initial distribution, the
complete kernel satisfies `K(Cx,Cy)=K(x,y)` and the distribution satisfies
`p(Cx)=p(x)`. Reindexing the finite pushforward sum gives
`(pK)(Cy)=(pK)(y)`. Induction gives the same identity after any finite number
of attempts. If `f(Cx)=-f(x)`, reindexing its expectation gives `E f=-E f`,
hence `E f=0`.

[`Dynamics.ChargeConjugation`](../../Lean/Dynamics/ChargeConjugation.lean)
formalizes these statements using `kernelPush` and `chargeExpectation` from
the existing protected-charge library. Positivity and unit row sums place
the kernel under the existing finite public-Markov theorem; symmetry does
not require a thermal equilibrium or a physical clock. The native branch
identity supplies the source-specific equivariance premise, with uniform
loads and unbiased coins supplying the preparation and measure premises.
This is an executable premise check combined with a generic Lean theorem,
not a Lean extraction of the Python source.

The same checker verifies all 7,776 complete branches of the 1,296-state
[`native_abundance_epochs`](../native_abundance_epochs/PROOF.md) population
kernel. Its involution is `(x0,x1,x2,k) -> (5-x0,5-x1,5-x2,-k mod 6)`.
It compares every reflected transition row and verifies the even initial
measure. The generic finite-time theorem therefore applies to that exact
population evolution without a separate Monte Carlo symmetry estimate.

For `k=1,2`, the observable `(n_k-n_(6-k))/N_faces` is odd and bounded by one.
The paired federation run verifies its zero paired mean after each sweep,
while retaining the nonzero counts of individual histories. Classes zero
and three are self-conjugate and require more information to distinguish
conjugate physical interpretations. The pathwise zero sum of the integer
lift follows separately from the closed-surface curl identity; neither that
sum nor the class imbalance is identified with baryon number.

## Sensitivity to broken premises

Let `a_0=TV(p_0,C_*p_0)` and let
`epsilon_t >= sup_x TV(K_t(Cx),C_*K_t(x))`. Couple an initial state and the
conjugated comparison with agreement probability `1-a_0`. On agreement,
maximal coupling preserves agreement at the next step with probability at
least `1-epsilon_t`. Thus

```
TV(p_n,C_*p_n) <= 1-(1-a_0) product_t(1-epsilon_t).
|E_pn f| <= TV(p_n,C_*p_n)             for odd |f| <= 1.
```

The second bound follows by writing `2 E f = E_p f-E_(C_*p) f` and using
the total-variation bound for a function with range at most two. For an
otherwise equivariant step driven by `Bernoulli(p_t)`, pushforward contracts
total variation, so `epsilon_t <= |2p_t-1|`. The exact rational executable
provides a conservative budget; it does not assert that each source attains
the worst case. Fair coins give zero added asymmetry.

Conditioning on a conjugation-even event of probability `q>0` bounds the
conditional odd mean by `min(1,a_n/q)`. An odd selection rule has no such
zero-mean guarantee: a symmetric fair sign conditioned to be positive has
mean one. Complete histories, including rejected events, are therefore
necessary for an asymmetry claim. These bounds concern finite operational
attempts and do not supply seconds, temperature or cosmological expansion.
