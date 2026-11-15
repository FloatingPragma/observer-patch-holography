# Exact seven-state control for contact observations

The complete periodic Weyl-contact action remains owned by
[`bound_pair`](../bound_pair/README.md). This package gives a smaller exact
control for one declared observation domain: side `L=4`, total momentum
`K=0`, states in the contact-generated cyclic subspace, and complete
contact/complement instruments between Hamiltonian evolutions. It reuses the
rank-one resolvent theorem; it introduces no new action, physical particle,
reference or continuum claim.

Write `H=T+U|d><d|`, with `T=dGamma(h)`,
`h(k)=-c sigma.sin(k)` and positive rational `U,c` satisfying
`U^2>12c^2`. At `K=0`, `h(-k)=-h(k)`. The singlet contact vector has free
spectral weight equally split between `+2c|sin k|` and `-2c|sin k|`.
On the side-four reciprocal grid, each component has `sin^2 k_j` equal to
zero or one with equal frequency. Therefore `r=sum_j sin^2 k_j` has weights
`binomial(3,r)/8`. The free contact measure has exactly seven atoms:

```
energy: 2c * (-sqrt(3), -sqrt(2), -1, 0, 1, sqrt(2), sqrt(3))
weight:      (1,        3,        3, 2, 3, 3,       1) / 16
```

If `P_lambda` is a free spectral projector, its seven vectors
`P_lambda d / sqrt(weight_lambda)` are orthonormal. Their span is invariant
under `T` and the contact projector, hence under `H` and both contact
instrument outputs. They give the exact reduced matrix
`diag(energy)+U u u*`, where `u=sqrt(weight)` entrywise. This proves the
reduction for arbitrary states in that subspace, arbitrary untouched
references, and any finite sequence of complete contact/complement reads.

The remaining 113 states of the 120-dimensional physical fiber are contact
dark and invariant under this protocol. They are necessary for other
inputs and observables. In particular, this reduction is not a license to
remove scattering states, Weyl pockets, other momentum fibers, or the
many-site occupation pointers used by `band_clock`. All seven reduced
states are retained, including its six scattering states.

The exact free resolvent is the rational function

```
g(E) = 1/(8E) + sum_(r=1)^3 [binomial(3,r)/8] E/(E^2-4c^2 r).
```

Above the free spectrum it decreases strictly. The isolated upper pole
solves `U g(E)=1`. At `E=U` the left side exceeds one: pair the positive
and negative free atoms, for which `E/(E^2-a^2)>1/E` whenever `a>0`.
At `E=U+4c` every positive denominator `E-lambda` exceeds `U`, so the
left side is below one. Rational bisection thus supplies an exact pole
interval with width `4c/2^iterations`; neither floating eigenvalues nor an
irrational comparison enters the certificate. At the canonical couplings
the pole is approximately `1.088594487779588`.

The independent test constructs the seven-dimensional inclusion by spectral
projector polynomials of the complete 120-dimensional physical free
operator, checks intertwining of both the Hamiltonian and contact
instrument, and compares every output of a two-read history with a complex
reference. The dense matrices are numerical diagnostics; the rational
secular signs and the reduction argument supply the exact statements.

This shortcut makes pole isolation and repeated contact protocols cheaper
and independently checkable. It does not simplify the full spatial source
execution or discharge its native reference and physical calibration
premises.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/pair_cyclic_control
```

## A sharper bound on the complete field clock

A second, independent simplification in `spectral_bound.py` applies to the
**full field**, not only the seven-state contact cyclic subspace. Let
`H=T+UP`, with `T=T*`, `P=P*=P²`, `PTP=0`, `U>0`, and `||T||<=b`.
For `E>b`, the convex function `1/(E-x)` lies below its affine endpoint
chord on `[-b,b]`. Functional calculus and compression by P give

```
P(E-T)^-1P <= E/(E²-b²) P.
```

Here the term proportional to T vanishes because `PTP=0`; commutation of
P with T is not assumed. An eigenvector of H with eigenvalue E>b has
nonzero P component satisfying `P(E-T)^-1P v=v/U`. Thus

```
E <= (U+sqrt(U²+4b²))/2.
```

The same spectral upper bound follows in infinite dimension by the resolvent
criterion above that value. The lower bound is `H>=-b I`, since `UP>=0`.
These are bounded-operator statements under their stated hypotheses.

For the canonical complete N=2 field, `U=1` and `b²=3/16`. Rational outward
bounds `E_min>-7/16` and `E_max<7/6` yield

```
||H-(35/96)I|| < 77/96.
```

A common scalar shift changes a unitary by a global phase, so the same
record-derived rate confidence event `|beta_hat-beta|<=3/80` gives the
optional rate-error half-diamond bound `77|s|/2560`, for `|s|<=1`, uniformly
with arbitrary untouched references. This improves the prior conservative
`9|s|/256` term by 14.4%. The splitting and statistical failure terms are separate, and no new preparation, physical clock, or reference claim is
made. The existing `band_clock` receipt retains its valid original bound.

The helper checks the rational edge inequalities exactly. Tests also inspect
all 8,128 physical eigenvalues across every momentum fiber, including
scattering states, and include a noncommuting two-state example saturating
the upper-root expression. Numerical spectral checks are diagnostics; the
operator argument above is the proof.
