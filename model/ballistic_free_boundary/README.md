# A power-law free-continuum class from ballistic propagation

For the same centered source, clock and canonical fields, the sufficient
free-continuum class improves from `U=o((log c)^(3/2))` to

\[
 c=v/a\longrightarrow\infty,\qquad 0\le U(c)=o(c^{2/7}),\qquad v>0.
\]

Every full field limit in this class is the free eight-valley Weyl vacuum
hierarchy. Its particle and hole energies are `v|p|` and its relativistic
rest masses are zero. These are the canonical charged fields on the
**centered** clock `K=H-(U/2)N`; this does not replace the original
uncentered ground or clock. No physical species assignment follows.

The exponent `2/7` is a sufficient bound, not a phase boundary. A coupling
outside this class is not thereby proved interacting. The distinct protected
Heisenberg spin limit remains outside this result. The theorem narrows the
search for a surviving interaction; it does not construct that interaction.

## Retained source inputs

The [centered ground construction](../centered_source_ground/README.md)
gives the exact occupation-defect bound
`Delta<=min(2,11(U/c)^2/8)`. The
[free-limit construction](../continuum_free_limit/README.md) supplies the
canonical composite `d=(2n_opposite-1)c`, the exact centered contact equation,
and its uniform microscopic-time response estimate. Write `u=U/c` and use
unit lattice smearings. If `B(s)` bounds both absolute row and column sums
of the speed-one free propagator, that estimate is

\[
 |F(s)-F_0(s)|
 \le 3u|s|+(300B(s)^3+12B(s))\sqrt{\Delta}.
 \tag{1}
\]

Here `F(s)` is the characteristic function of the full positive signed
addition/removal response measure of `d`, with frequency divided by `c`.
Both charge branches have nonnegative physical excitation energies. The
free composite response satisfies

\[
 \mu_0([-x,x])\le32x^6\quad(x>0).
 \tag{2}
\]

The earlier logarithmic restriction came from the exponential estimate
`B(s)<=exp(3|s|)`. The free source itself has a stronger bound, proved next.
No assumption of exponential clustering, a spectral gap, perturbative
convergence, or a spin-wave approximation is added.

## Polynomial absolute-sum bound for the exact free source

On the normalized Brillouin torus let

\[
 h(k)=-\sum_{j=1}^3\sigma_j\sin k_j,\qquad
 V_s(k)=e^{-is h(k)},\qquad
 K_s(n)=\int e^{ik\cdot n}V_s(k)\,{d^3k\over(2\pi)^3}.
\]

The sign of the time convention will not affect any norm. The first and
second same-coordinate derivatives of `h` have operator norm at most one.
Differentiating the exponential by its Duhamel integral gives

\[
 \|\partial_j V_s\|\le |s|,\qquad
 \|\partial_j^2 V_s\|\le |s|+s^2,\qquad
 \|\Delta_k V_s\|\le3|s|+3s^2.
 \tag{3}
\]

The second derivative has one term differentiating `h` twice and two
ordered first-derivative insertions. Their time-simplex volumes sum to
`s^2`, so (3) does not require commutation of the Pauli matrices.

Put `L=1+|s|` and fix either spin row `alpha`. Parseval, the unitarity of
`V_s`, and (3) give

\[
 \begin{split}
 \sum_{n,\beta}\left(1+{|n|^2\over L^2}\right)^2
       |K_s(n)_{\alpha\beta}|^2
 &\le1+{6s^2\over L^2}
          +{9s^2(1+|s|)^2\over L^4}\\
 &\le16.
 \end{split}\tag{4}
\]

For the fourth spatial moment, Parseval is applied to `Delta_k V_s`,
whose Fourier multiplier is `-|n|^2`; no mixed derivatives are omitted.
The same argument works for a fixed column.

There is an elementary lattice bound

\[
 \sum_{n\in\mathbb Z^3}(1+|n|^2/L^2)^{-2}\le79L^3
 \quad(L\ge1).
 \tag{5}
\]

Indeed let `M=floor L`, so `L/2<=M<=L`. The inner infinity-norm cube
contains at most `(2M+1)^3<=27L^3` sites. The shell `|n|_infinity=r`
has `24r^2+2<=26r^2` sites, and its summand is at most `L^4/r^4`.
The exterior sum is at most
`26L^4 sum_(r>M) r^-2<=26L^4/M<=52L^3`.

Weighted Cauchy--Schwarz in `(n,beta)`, retaining the two spin indices,
now proves

\[
 B(s)\le\sqrt{2\cdot16\cdot79}\,(1+|s|)^{3/2}
       <64(1+|s|)^{3/2}.
 \tag{6}
\]

The derivative bounds also establish the required weighted square
summability before the Cauchy--Schwarz argument. Thus no absolute kernel
summability was assumed in deriving it. This is an infinite-lattice
bound, uniform over the external smearing, not a finite-volume fit.

## Mixed-time Wick comparison without expanding every free field

There is a second improvement to (1). Let `f=e^{ish} delta_(x,sigma)`
be a single freely evolved site mode. The actual centered ground is
translation invariant, so its one-particle covariance is a matrix multiplier
`rho(k)`, whether or not it commutes with `h(k)`. Since
`h(k)P_+(k)=e(k)P_+(k)`,

\[
 \langle P_+f,\rho P_+f\rangle
 =\int e_\sigma^*P_+(k)\rho(k)P_+(k)e_\sigma
       \,{d^3k\over(2\pi)^3}\le\Delta.
\]

The free phases cancel pointwise. The negative-band creator has the same
bound using `I-rho`. Both wrong-band ground-vector norms are therefore
at most `sqrt(Delta)`, uniformly in position, spin and free time.

The ordinary CAR Wick recursion can now be applied directly to any fixed
word of these normalized freely evolved site fields and unevolved site
fields. Moving its first annihilating part through the word gives the same
constants `C_2=2`, `C_4=8` as for site fields at time zero. No position-kernel
expansion of the remaining word is necessary.

In the graded anticommutator of two cubic fields there are at most nine
single CAR contractions, each leaving a four-field word. Only the cross
contraction needs a free propagation kernel. Its external smearing sum is
bounded by `B(s)` using the Schur bound. The spin indices are permutations
of the external spin labels, which preserve that same bound. The mixed
cubic/linear terms have at most three contractions and two remaining fields.
Consequently the state discrepancy is bounded by

\[
 (4\cdot9\cdot8+2\cdot3\cdot2+2\cdot3\cdot2)
 B(s)\sqrt{\Delta}=312B(s)\sqrt{\Delta}.
 \tag{6a}
\]

This is a bound for the actual possibly non-Gaussian ground. It uses neither
Wick factorization in that ground nor a uniform small operator norm for
`rho-P_-`. Finite support proves the smearing calculation first; the exact
CAR norm of the composite extends it to all square-summable smearings.

## Composite spectral escape with a power bound

Since `sqrt(11/8)<6/5`, the microscopic dynamical comparison, (6) and (6a) imply

\[
 |F(s)-F_0(s)|\le 2^{15}u(1+|s|)^{3/2}.
 \tag{7}
\]

The constant follows from the exact rational inequality
`3+(6/5)*312*64<2^15`. For `0<epsilon<=1`, Gaussian smoothing
samples `s=Z/epsilon`, where `Z` is a standard real normal variable. Then

\[
 \mathbb E(1+|Z|/\epsilon)^{3/2}
 \le\epsilon^{-3/2}\mathbb E(1+|Z|)^2
 =\epsilon^{-3/2}(2+2\sqrt{2/\pi})
 <4\epsilon^{-3/2}.
\]

Integrating (2) by parts bounds the free Gaussian response by
`1536 epsilon^6`. Positivity of the full response measure therefore
gives, for either positive-energy GNS charge branch and every unit
lattice smearing `f`,

\[
 \|1_{[0,R]}(L_a)d(f)\Omega_a\|^2\le M_c(R),
\]
\[
 M_c(R)\le\min\left\{1,
 e^{R^2/(2c^2\epsilon^2)}
 [1536\epsilon^6+2^{17}u\epsilon^{-3/2}]\right\}.
 \tag{8}
\]

The corresponding adjoint vector has the same bound. The projection
formula is schematic for the chosen branch, as in the parent theorem;
the signed response frequencies do not redefine physical energy.

Choose, for `c>=1` and eventually `u<=1`,

\[
 R=\sqrt c,\qquad
 \epsilon=\max\{u^{2/15},c^{-1/4}\}.
\]

This also handles `U=0` and arbitrarily fast decreasing positive `U`.
It gives the explicit bound

\[
 U^2M_c(\sqrt c)
 \le e^{1/(2\sqrt c)}\left[
 (1536+2^{17})U^{14/5}c^{-4/5}
       +1536U^2c^{-3/2}\right].
 \tag{9}
\]

For the second term in (8), use `epsilon>=u^(2/15)` when `u>0`;
when `u=0` that term vanishes. For the first use
`epsilon^6<=u^(4/5)+c^(-3/2)`. Thus no division by `U=0` occurs.
Under `U=o(c^(2/7))`, both terms in (9) tend to zero, and
`U/R=U/sqrt(c)->0` as well. At `U=c^(2/7)`, the first expression in
(9) has power zero; this proof does not establish decay there.

## Returning to the full canonical field hierarchy

The exact centered field equation and positive GNS energy give, on an
input kinetic window bounded by `P`, for either field adjoint,

\[
 \|(C_a(f,t)-C_a(e^{it h_a}f,0))\Omega_a\|
 \le {U\over2}\left[|t|\sqrt{M_c(R)}
                    +{2+|t|P\over R}\right]\|f\|.
 \tag{10}
\]

The initial wrong-band annihilator/creator bound is
`(U/2)[sqrt(M_c(R))/b+1/(R+b)]` on a kinetic band with gap `b>0`.
Equations (9)--(10) make all these errors vanish. The
[growing-contact reconstruction](../strong_coupling_refinement/README.md#why-unbounded-u-does-not-silently-reuse-a-bounded-u-compactness-premise)
then applies unchanged: universal bounded-field compactness avoids an
unjustified time-continuity bound containing `U`; equal-time CAR, both
vacuum identities and the onsite-interaction-independent locality bound
identify every full all-word cluster state with the free Weyl vacuum.
It is the same source, centered state, field normalization and clock.
The proof does not replace an interacting state by a guessed free state.

The nonlinear composite retains its exact equal-time CAR norm while its
effect on these finite-energy canonical fields vanishes. This remains
consistent with the separate composite-channel and protected-spin results.

## Exact finite bounds and verification

`bounds.py` checks the rational constant inequalities and implements (8)
using `epsilon=q^2` with rational `q`. It uses
`exp(x)<=1/(1-x)` for `0<=x<1`, or the trivial weight bound otherwise.
For the explicit family

```
c=b^80, U=b^20, R=b^40, epsilon=b^-8, b>=2,
```

the weighted bound before the trivial cap is exactly

\[
 {(1536+2^{17})b^{-8}\over1-1/(2b^{64})}.
\]

Thus `U=c^(1/4)` is a concrete power-growing example in the proved class.
The retained finite table also gives (10) at `|t|=P=1`. These deliberately
conservative scales are mathematical certificates, not proposed physical
parameters or a measurement of nature's particle masses.

Run `python -m pytest -q -W error model/ballistic_free_boundary`, then
`python -m model.ballistic_free_boundary --write` or `--check`. The receipt
binds this proof, its arithmetic and the retained parent derivations. The
analytic argument above, rather than finite examples, establishes the
infinite-lattice and continuum conclusions. No CI workflow is added.
